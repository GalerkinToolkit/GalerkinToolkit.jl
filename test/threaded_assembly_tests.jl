module ThreadedAssemblyTests

using Test
using LinearAlgebra
using ForwardDiff
import GalerkinToolkit as GT

@testset "Element ranges preserve domain-local IDs" begin
    @test isbitstype(GT.SerialExecution)
    @test GT.SerialExecution().fastmath
    @test GT.CPUThreadedExecution().fastmath
    @test_throws ArgumentError GT.CPUThreadedExecution(workers=0)
    @test_throws ArgumentError GT.CPUThreadedExecution(distribution=:invalid)
    for n in (0,1,7,16), workers in (1,3,20), distribution in (:contiguous,:cyclic)
        count = min(n,workers)
        ranges = [GT.assembly_element_range(GT.AssemblyWorker(i,count,n,Val(distribution))) for i in 1:count]
        @test length(ranges) == min(n,workers)
        @test sort(reduce(vcat,collect.(ranges);init=Int[])) == collect(1:n)
        @test all(!isempty,ranges)
    end
    @test [GT.assembly_element_range(GT.AssemblyWorker(i,3,7,Val(:contiguous))) for i in 1:3] == [1:3,4:5,6:7]
    @test [GT.assembly_element_range(GT.AssemblyWorker(i,3,7,Val(:cyclic))) for i in 1:3] == [1:3:7,2:3:7,3:3:7]
    for distribution in (:contiguous,:cyclic)
        worker = GT.AssemblyWorker(Int32(2),Int32(3),Int32(7),Val(distribution))
        @test isbitstype(typeof(worker))
        @test eltype(GT.assembly_element_range(worker)) === Int32
    end
    worker = GT.AssemblyWorker(Int32(24),Int32(24),typemax(Int32),Val(:contiguous))
    @test last(GT.assembly_element_range(worker)) == typemax(Int32)
end

@testset "Logical worker ID owns scratch and computes its range" begin
    function worker_kernel!(output,workspace,worker)
        workspace[1][1] = worker.id
        yield()
        for element in GT.assembly_element_range(worker)
            output[element] = workspace[1][1]
        end
        nothing
    end
    factory = previous -> previous === nothing ? (zeros(Int,1),) : previous
    prepare = () -> GT.PreparedElementAssembly(worker_kernel!,factory)
    for distribution in (:contiguous,:cyclic),workers in (1,3,20)
        output = zeros(Int,13)
        execution = GT.CPUThreadedExecution(;workers,distribution)
        runner = GT.ThreadedAssemblyRunner(prepare,output,13,execution)
        GT.run_assembly!(runner,())
        @test runner.nworkers == min(workers,13)
        for id in 1:runner.nworkers
            worker = GT.AssemblyWorker(id,runner.nworkers,13,Val(distribution))
            @test all(==(id),output[GT.assembly_element_range(worker)])
        end
        @test all(!iszero,output)
        before = copy(runner.workspaces)
        records = runner.worker_records
        tasks = runner.tasks
        old_tasks = runner.nworkers == 1 ? Task[] : copy(tasks)
        GT.run_assembly!(runner,())
        @test all(i->runner.workspaces[i] === before[i],eachindex(before))
        @test runner.worker_records === records
        @test runner.tasks === tasks
        if runner.nworkers > 1
            @test all(i->records[i].workspace === before[i],eachindex(before))
            @test all(i->tasks[i] !== old_tasks[i] && istaskdone(tasks[i]),eachindex(tasks))
            @test all(task->fetch(task) === nothing,tasks)
        end
    end
end

@testset "Worker records refresh captures and recover after failure" begin
    output = zeros(6)
    completed = Threads.Atomic{Int}(0)
    function prepare(scale)
        factory = previous -> previous === nothing ? (zeros(1),) : previous
        kernel = (output,workspace,worker) -> begin
            workspace[1][1] = scale
            yield()
            for element in GT.assembly_element_range(worker)
                output[element] = workspace[1][1]
            end
            Threads.atomic_add!(completed,1)
            scale < 0 && worker.id == 1 && error("intentional kernel failure")
            nothing
        end
        GT.PreparedElementAssembly(kernel,factory)
    end
    runner = GT.ThreadedAssemblyRunner(prepare,output,6,GT.CPUThreadedExecution(workers=3))
    GT.run_assembly!(runner,(2.0,))
    records,tasks = runner.worker_records,runner.tasks
    old_records = copy(records)
    pool = runner.workspaces
    for scale in (3.0,-1.0,4.0)
        completed[] = 0
        if scale < 0
            @test_throws CompositeException GT.run_assembly!(runner,(scale,))
        else
            GT.run_assembly!(runner,(scale,))
        end
        @test completed[] == 3
        @test output == fill(scale,6)
        @test runner.worker_records === records
        @test all(i->records[i] === old_records[i],eachindex(records))
        @test runner.tasks === tasks
        @test all(istaskdone,tasks)
        @test runner.workspaces === pool
    end
    # A different capture type rebuilds the typed records while retaining the
    # task-handle vector and compatible numerical workspace arrays.
    GT.run_assembly!(runner,(5,))
    @test output == fill(5,6)
    @test runner.worker_records !== records
    @test runner.tasks === tasks
    @test runner.workspaces === pool
    @test all(i->runner.worker_records[i].workspace === records[i].workspace,eachindex(records))
end

@testset "Range lowering selects atomic contribution calls" begin
    expr = :((input,) -> (() -> begin
        alloc -> begin
            for face in 1:length(input)
                GT.contribute_element_block!(alloc,input[face],face,0,2,2)
                GT.contribute_element!(alloc,input[face],1,1,face)
                # TODO: Re-enable DG sides and vector lowering.
                # GT.contribute_element!(alloc,input[face],1,1,face,2,1)
                # GT.contribute!(alloc,input[face],input[face],1)
            end
        end
    end))
    ranged = GT.ast_element_range_assembly(expr,:face)
    calls = Any[]
    function visit(node)
        node isa Expr || return
        node.head === :call && push!(calls,first(node.args))
        foreach(visit,node.args)
    end
    visit(ranged)
    @test :(GT.contribute_element_block_atomic!) in calls
    @test :(GT.contribute_element_atomic!) in calls
    # @test :(GT.contribute_atomic!) in calls
    @test :(GT.assembly_element_range) in calls
    @test !any(c->c in (:(GT.contribute_element_block!),:(GT.contribute_element!),:(GT.contribute!)),calls)
end

@testset "Workspace preparation runs in the assembly task" begin
    caller = current_task()
    factory = previous -> (
        data=previous === nothing ? zeros(Int,1) : previous.data,
        prepared_by=current_task())
    function kernel!(output,workspace,worker)
        workspace.prepared_by === current_task() || error("Workspace prepared in another task")
        worker.count == 1 || current_task() !== caller || error("Preparation ran on caller")
        workspace.data[1] = worker.id
        yield()
        for element in GT.assembly_element_range(worker)
            output[element] = workspace.data[1]
        end
        nothing
    end
    for workers in (1,4),distribution in (:contiguous,:cyclic)
        output = zeros(Int,11)
        prepare = () -> GT.PreparedElementAssembly(kernel!,factory)
        runner = GT.ThreadedAssemblyRunner(prepare,output,11,GT.CPUThreadedExecution(;workers,distribution))
        GT.run_assembly!(runner,())
        old_pool = runner.workspaces
        old_data = [workspace.data for workspace in old_pool]
        for id in 1:runner.nworkers
            worker = GT.AssemblyWorker(id,runner.nworkers,11,Val(distribution))
            @test all(==(id),output[GT.assembly_element_range(worker)])
        end
        GT.run_assembly!(runner,())
        @test runner.workspaces === old_pool
        @test all(id->runner.workspaces[id].data === old_data[id],eachindex(old_data))
    end

    calls = Threads.Atomic{Int}(0)
    failed_factory = previous -> begin
        Threads.atomic_add!(calls,1)
        error("intentional factory failure")
    end
    prepare = () -> GT.PreparedElementAssembly((args...)->nothing,failed_factory)
    runner = GT.ThreadedAssemblyRunner(prepare,nothing,3,GT.CPUThreadedExecution(workers=3))
    @test_throws CompositeException GT.run_assembly!(runner,())
    @test calls[] == 3 # Every task has joined before returning the error.
    @test runner.workspaces === nothing
end

function workspace_tests()
    setup_calls,refresh_calls = Ref(0),Ref(0)
    input = collect(1.0:8.0)
    shared_setup = input -> (setup_calls[] += 1; input)
    refresh = (input,scale) -> (refresh_calls[] += 1; input .* scale)
    expr = :((config,input) -> ((scale,n) -> ((output,worker) -> begin
        shared = config.shared_setup(input)
        updated = config.refresh(shared,scale)
        scratch = alloc_zeros("scratch",Float64,n)
        alias = scratch
        other = Base.zeros(Float64,n)
        for element in GT.assembly_element_range(worker)
            fill!(alias,updated[element])
            fill!(other,2updated[element])
            output[element] = sum(alias)+sum(other)+shared[element]
        end
        nothing
    end)))
    moved = GT.ast_initialization_code_motion(expr;workspaces=true)
    prepare = GT.evaluate(moved,((;shared_setup,refresh),input))
    @test setup_calls[] == 1
    prepared = Base.invokelatest(prepare,2.0,2)
    @test refresh_calls[] == 1
    factory = prepared.workspace_factory
    w1,w2 = Base.invokelatest(factory,nothing),Base.invokelatest(factory,nothing)
    @test length(w1) == length(w2) == 2
    @test !Base.mightalias(w1[1],w2[1])
    @test !Base.mightalias(w1[2],w2[2])
    @test getfield(prepared.kernel,:shared) === input
    @test getfield(prepare,:shared) === input
    @test Base.invokelatest(factory,w1) === w1
    output = zeros(8)
    Base.invokelatest(prepared.kernel,output,w1,GT.AssemblyWorker(1,2,8,Val(:cyclic)))
    Base.invokelatest(prepared.kernel,output,w2,GT.AssemblyWorker(2,2,8,Val(:cyclic)))
    @test output == 13input
    updated = Base.invokelatest(prepare,3.0,2)
    @test Base.invokelatest(updated.workspace_factory,w1) === w1
    @test setup_calls[] == 1
    @test refresh_calls[] == 2
    resized = Base.invokelatest(prepare,3.0,3)
    w3 = Base.invokelatest(resized.workspace_factory,w1)
    @test all(b->length(b)==3,w3)
    @test all(i->w3[i] !== w1[i],eachindex(w1))

    for distribution in (:contiguous,:cyclic)
        runner = GT.ThreadedAssemblyRunner(prepare,output,8,GT.CPUThreadedExecution(;workers=3,distribution))
        GT.run_assembly!(runner,(2.0,2))
        @test output == 13input
        old = copy(runner.workspaces)
        GT.run_assembly!(runner,(3.0,2))
        @test output == 19input
        @test all(i->runner.workspaces[i] === old[i],eachindex(old))
        @test all(i->!Base.mightalias(old[1][1],old[i][1]),2:length(old))
    end
end

@testset "Shared setup and private reusable buffers" begin
    workspace_tests()
end

@testset "Fixed workspace factories are defined in initialization" begin
    expr = :((input,) -> ((scale,) -> ((output,worker) -> begin
        scratch = alloc_zeros("scratch",Float64,2)
        alias = scratch
        for element in GT.assembly_element_range(worker)
            fill!(alias,scale*input[element])
            output[element] = sum(alias)
        end
        nothing
    end)))
    moved = GT.ast_initialization_code_motion(expr;workspaces=true)
    initialization = GT.ast_children(last(GT.ast_children(moved)))
    factory_definitions = filter(initialization) do statement
        GT.ast_is_definition(statement) && GT.ast_head(GT.ast_rhs(statement)) === :(->)
    end
    @test length(factory_definitions) == 1
    factory_name = GT.ast_lhs(only(factory_definitions))
    parameter_body = last(GT.ast_children(last(initialization)))
    prepared_expr = last(GT.ast_children(parameter_body))
    @test last(GT.ast_children(prepared_expr)) === factory_name

    input = [1.0,2.0,3.0]
    prepare = GT.evaluate(moved,(input,))
    prepared = Base.invokelatest(prepare,2.0)
    updated = Base.invokelatest(prepare,3.0)
    workspace = Base.invokelatest(prepared.workspace_factory,nothing)
    other = Base.invokelatest(updated.workspace_factory,nothing)
    @test !Base.mightalias(only(workspace),only(other))
    @test Base.invokelatest(updated.workspace_factory,workspace) === workspace
    output = zeros(3)
    worker = GT.AssemblyWorker(1,1,3,Val(:contiguous))
    Base.invokelatest(prepared.kernel,output,workspace,worker)
    @test output == 4input
    Base.invokelatest(updated.kernel,output,workspace,worker)
    @test output == 6input
    Base.invokelatest(prepared.kernel,output,workspace,worker)
    @test output == 4input
end

@testset "Empty-domain guards and changing workspace types" begin
    expr = :((input,) -> ((n,T) -> ((output,worker) -> begin
        if n == 0
            return
        end
        value = first(input)
        scratch = alloc_zeros("scratch",T,n)
        for element in GT.assembly_element_range(worker)
            fill!(scratch,value)
            output[element] = sum(scratch)
        end
        nothing
    end)))
    prepare = GT.evaluate(GT.ast_initialization_code_motion(expr;workspaces=true),([2.0],))
    output = zeros(3)
    runner = GT.ThreadedAssemblyRunner(prepare,output,3,GT.CPUThreadedExecution(workers=3))
    for (n,T) in ((0,Float64),(2,Float64),(3,Float32),(0,Float32),(1,Float64))
        fill!(output,0)
        GT.run_assembly!(runner,(n,T))
        @test output == fill(2n,3)
        @test all(w->n==0 ? isempty(w) : (length(w[1])==n && eltype(w[1])===T),runner.workspaces)
    end
    empty_prepare = GT.evaluate(GT.ast_initialization_code_motion(expr;workspaces=true),(Float64[],))
    guarded = GT.ThreadedAssemblyRunner(empty_prepare,output,3,GT.CPUThreadedExecution(workers=2))
    GT.run_assembly!(guarded,(0,Float64))
    @test all(isempty,guarded.workspaces)
    empty = GT.ThreadedAssemblyRunner(_->error("must not prepare empty range"),output,0,GT.CPUThreadedExecution(workers=2))
    @test GT.run_assembly!(empty,()) === nothing
    @test empty.workspaces === nothing

    captured_guard = :((input,) -> (() -> ((output,elements) -> begin
        if isempty(input)
            return
        end
        value = first(input)
        scratch = alloc_zeros("scratch",Float64,1)
        fill!(scratch,value)
        output[1] = scratch[1]
    end)))
    guarded_prepare = GT.evaluate(GT.ast_initialization_code_motion(captured_guard;workspaces=true),(Float64[],))
    guarded = GT.ThreadedAssemblyRunner(guarded_prepare,output,1,GT.CPUThreadedExecution(workers=1))
    @test GT.run_assembly!(guarded,()) === nothing
    @test only(guarded.workspaces) === ()
end

const matrix_method = (;matrix=GT.element_csc_assembly())
# TODO: Re-enable vector assembly.
# const vector_method = (;vector=GT.monolithic_assembly(GT.direct_vector_assembly()))

function warmed_update_bytes(update!::F,output,cache,parameters) where F
    for _ in 1:3
        update!(output,cache;parameters)
    end
    @allocated update!(output,cache;parameters)
end

function generated_assembly_tests()
    mesh = GT.cartesian_mesh((0,1,0,1),(3,2);simplexify=false)
    omega = GT.interior(mesh)
    V = GT.lagrange_space(omega,2)
    measure = GT.measure(omega,4)
    scale = Ref(2.0)
    coefficient = GT.analytical_field(x -> scale[],omega)
    # forcing = GT.analytical_field(x -> scale[]*(1+sum(x)),omega)
    gradient = ForwardDiff.gradient
    a = (u,v) -> GT.∫(x->coefficient(x)*(gradient(u,x)⋅gradient(v,x)+u(x)*v(x)),measure)
    # l = v -> GT.∫(x->forcing(x)*v(x),measure)
    reference_A = GT.assemble_matrix(a,Float64,V,V;assembly_method=matrix_method)
    # reference_b = GT.assemble_vector(l,Float64,V;assembly_method=vector_method)
    for (workers,distribution) in ((1,:contiguous),(3,:contiguous),(3,:cyclic))
        A,cache = GT.assemble_matrix(a,Float64,V,V;assembly_method=matrix_method,
            execution=GT.CPUThreadedExecution(;workers,distribution,fastmath=false),parameters=(scale,),reuse=Val(true))
        #= TODO: Re-enable threaded load vectors.
        b,bcache = GT.assemble_vector(l,Float64,V;assembly_method=vector_method,
            execution=GT.CPUThreadedExecution(;workers,distribution),parameters=(scale,),reuse=Val(true))
        =#
        @test A ≈ reference_A
        # @test b ≈ reference_b
        # old_A,old_b = copy(only(cache.loops).workspaces),copy(only(bcache.loops).workspaces)
        old_A = copy(only(cache.loops).workspaces)
        @test !isempty(first(old_A))
        # @test !isempty(first(old_b))
        # for pool in (old_A,old_b),i in 2:length(pool),slot in eachindex(first(pool))
        for pool in (old_A,),i in 2:length(pool),slot in eachindex(first(pool))
            @test !Base.mightalias(first(pool)[slot],pool[i][slot])
        end
        for (update!,output,state,expected,old) in
                # ((GT.update_matrix!,A,cache,reference_A,old_A),(GT.update_vector!,b,bcache,reference_b,old_b))
                ((GT.update_matrix!,A,cache,reference_A,old_A),)
            for factor in (3.0,1.0,2.0)
                scale[] = factor
                @test update!(output,state;parameters=(scale,)) === output
                @test output ≈ (factor/2)*expected
                @test all(i->only(state.loops).workspaces[i] === old[i],eachindex(old))
            end
        end
    end

    # Different domains keep separate domain-local cache columns, including
    # reordered subsets; their kernels contribute to the same numerical output.
    GT.group_faces(mesh,2)["threaded subset"] = Int32[6,2]
    subset = GT.interior(mesh;group_names=["threaded subset"])
    subset_measure = GT.measure(subset,4)
    overlapping = (u,v) -> GT.∫(x->u(x)*v(x),measure)+GT.∫(x->2u(x)*v(x),subset_measure)
    expected = GT.assemble_matrix(overlapping,Float64,V,V;assembly_method=matrix_method)
    A,cache = GT.assemble_matrix(overlapping,Float64,V,V;assembly_method=matrix_method,
        execution=GT.CPUThreadedExecution(workers=4,distribution=:cyclic),reuse=Val(true))
    @test A ≈ expected
    @test length(cache.loops) == 2
    @test all(r->r.argument.allocation.nzval === A.nzval,cache.loops)
    @test GT.update_matrix!(A,cache) ≈ expected

    #= TODO: Re-enable threaded multi-field DG and vector assembly.
    # Distinct local blocks and disabled field pairs must keep independent
    # workspaces while using the existing literal bases and field selection.
    V1 = GT.lagrange_space(omega,1;continuous=Val(false))
    mixed = GT.cartesian_product(V1,V)
    mixed_form = ((u1,u2),(v1,v2)) -> GT.∫(x->(u1(x)+u2(x))*(v1(x)+v2(x)),measure)
    mixed_load = ((v1,v2),) -> GT.∫(x->v1(x)+v2(x),measure)
    options = (;matrix=(;block_mask=Bool[true false; true true]))
    mixed_method = (;matrix=GT.monolithic_assembly(GT.element_csc_assembly()))
    reference = GT.assemble_matrix(mixed_form,Float64,mixed,mixed;
        assembly_method=mixed_method,assembly_options=options)
    A,cache = GT.assemble_matrix(mixed_form,Float64,mixed,mixed;
        assembly_method=mixed_method,assembly_options=options,execution=GT.CPUThreadedExecution(workers=3),reuse=Val(true))
    @test A ≈ reference
    @test GT.update_matrix!(A,cache) ≈ reference
    expected_b = GT.assemble_vector(mixed_load,Float64,mixed;assembly_method=vector_method)
    b,cache = GT.assemble_vector(mixed_load,Float64,mixed;assembly_method=vector_method,
        execution=GT.CPUThreadedExecution(workers=3,distribution=:cyclic),reuse=Val(true))
    @test b ≈ expected_b
    @test GT.update_vector!(b,cache) ≈ expected_b

    =#
    # Zero CSC positions skip entries outside the requested free/Dirichlet selection.
    constrained = GT.lagrange_space(omega,1;dirichlet_boundary=GT.boundary(mesh))
    mass = (u,v) -> GT.∫(x->u(x)*v(x),measure)
    # load = v -> GT.∫(x->v(x),measure)
    for selection in (GT.FREE,GT.DIRICHLET)
        workers = selection == GT.FREE ? 3 : 1
        expected = GT.assemble_matrix(mass,Float64,constrained,constrained;
            assembly_method=matrix_method,free_or_dirichlet=(GT.FREE,selection))
        A,cache = GT.assemble_matrix(mass,Float64,constrained,constrained;
            assembly_method=matrix_method,free_or_dirichlet=(GT.FREE,selection),
            execution=GT.CPUThreadedExecution(;workers),reuse=Val(true))
        @test A ≈ expected
        @test GT.update_matrix!(A,cache) ≈ expected
        if workers == 1
            @test warmed_update_bytes(GT.update_matrix!,A,cache,()) == 0
        end
        #= TODO: Re-enable constrained vector assembly.
        expected_b = GT.assemble_vector(load,Float64,constrained;
            assembly_method=vector_method,free_or_dirichlet=selection)
        b,cache = GT.assemble_vector(load,Float64,constrained;assembly_method=vector_method,
            free_or_dirichlet=selection,execution=GT.CPUThreadedExecution(;workers),reuse=Val(true))
        @test b ≈ expected_b
        @test GT.update_vector!(b,cache) ≈ expected_b
        if workers == 1
            @test warmed_update_bytes(GT.update_vector!,b,cache,()) == 0
        end
        =#
    end

    coo = (;matrix=GT.monolithic_assembly(GT.coo_assembly()))
    @test_throws ArgumentError GT.assemble_matrix(a,Float64,V,V;assembly_method=coo,
        execution=GT.CPUThreadedExecution())
    # The default matrix method (element CSC) supports threaded execution.
    @test GT.assemble_matrix(a,Float64,V,V;execution=GT.CPUThreadedExecution(workers=2),
        parameters=(scale,),reuse=Val(false)) ≈ reference_A
    # @test_throws ArgumentError GT.assemble_vector(l,Float64,V;execution=GT.CPUThreadedExecution())
    @test_throws TypeError GT.assemble_matrix(a,Float64,V,V;execution=:invalid)
end

@testset "Threaded generated single-field CG matrix assembly" begin
    generated_assembly_tests()
end

@testset "Threaded generated complex matrix assembly" begin
    mesh = GT.cartesian_mesh((0,1,0,1,0,1),(2,2,2);simplexify=true)
    omega = GT.interior(mesh)
    gradient = ForwardDiff.gradient
    for k in (1,2)
        V = GT.lagrange_space(omega,k)
        measure = GT.measure(omega,2k)
        a = (u,v)->GT.∫(x->(2-1im)*(gradient(u,x)⋅gradient(v,x))+(0.5+3im)*u(x)*v(x),measure)
        for T in (ComplexF32,ComplexF64)
            expected = GT.assemble_matrix(a,T,V,V;assembly_method=matrix_method)
            @test eltype(expected) === T
            @test !iszero(imag.(expected))
            for (workers,distribution) in ((1,:contiguous),(2,:contiguous),(4,:cyclic))
                A,cache = GT.assemble_matrix(a,T,V,V;assembly_method=matrix_method,
                    execution=GT.CPUThreadedExecution(;workers,distribution,fastmath=false),reuse=Val(true))
                @test eltype(A) === T
                @test A ≈ expected
                @test GT.update_matrix!(A,cache) === A
                @test A ≈ expected
            end
        end
    end
end

end # module

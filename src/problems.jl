
# Measure

function measure(domain::AbstractDomain,degree)
    measure(quadrature,domain,degree)
end

function measure(f,domain::AbstractFaceDomain,degree)
    f(domain,degree)
end

function measure(f,domain::AbstractMeshDomain,degree)
    mesh = GT.mesh(domain)
    d = GT.num_dims(domain)
    dface_to_drefid = GT.face_reference_id(mesh,d)
    drefid_refdface = GT.reference_spaces(mesh,Val(d))
    refid_to_quad = map(drefid_refdface) do refdface
        geo = GT.domain(refdface)
        f(geo,degree)
    end
    mesh_quadrature(;
        domain,
        face_reference_id=dface_to_drefid,
        reference_quadratures=refid_to_quad)
end

# Integrals

function integrate(f,quadrature::AbstractQuadrature)
    coefficient = 1
    contribution = DomainContribution(f,quadrature,coefficient)
    integral(contribution)
end
const ∫ = integrate

struct DomainContribution{A,B,C}
    integrand::A
    quadrature::B
    coefficient::C
end

domain(a::DomainContribution) = domain(a.quadrature)
quadrature(a::DomainContribution) = a.quadrature
coefficient(a::DomainContribution) = a.coefficient

function prototype(a::DomainContribution)
    term = optimize(GT.term(a,index(Val(0))))
    prototype(term)
end

function replace_coefficient(a::DomainContribution,coefficient)
    DomainContribution(a.integrand,a.quadrature,coefficient)
end

function quantity(contribution::DomainContribution)
    (;integrand,quadrature) = contribution
    domain = GT.domain(contribution)
    x = coordinate_quantity(quadrature)
    dV = weight_quantity(quadrature)
    coefficient = GT.coefficient(contribution)
    alpha = uniform_quantity(coefficient)
    alpha*integrand(x)*dV
end

function term(contribution::DomainContribution,index)
    domain = GT.domain(contribution)
    opts = QuantityOptions(domain,index)
    qty = quantity(contribution)
    term(qty,opts)
end

function integral(contributions...)
    Integral(contributions)
end

struct Integral{A} <: AbstractType
    contributions::A
end

contributions(i::Integral) = i.contributions

function contributions(i::Number)
    @assert i == 0
    ()
end

function contribution(i::Integral,domain::AbstractDomain)
    for contribution in contributions(i)
        if domain == GT.domain(contribution)
            return contribution
        end
    end
    error("domain is not found")
end

function Base.:+(int1::Integral,int2::Integral)
    # TODO merge contributions on the same domain?
    contribs = (GT.contributions(int1)...,GT.contributions(int2)...)
    Integral(contribs)
end

function Base.:*(v::Number,int::Integral)
    contribs = map(GT.contributions(int)) do contribution
        coefficient = v*GT.coefficient(contribution)
        replace_coefficient(contribution,coefficient)
    end
    Integral(contribs)
end

function Base.:-(int1::Integral,int2::Integral)
    int1 + (-1)*int2
end

function Base.:-(int1::Real,int2::Integral)
    @assert int1 == 0
    (-1)*int2
end

function Base.:*(int::Integral,v::Number)
    v*int
end

function Base.:/(int::Integral,v::Number)
    (1/v)*int
end

"""
    AbstractAssemblyExecution

Execution policy for assembly, independent of the storage method. The built-in
policies enable `fastmath` by default.
"""
abstract type AbstractAssemblyExecution <: AbstractType end

"""
    SerialExecution(; fastmath=true)

Sequential assembly using the selected storage method. This is the default.
"""
struct SerialExecution <: AbstractAssemblyExecution
    fastmath::Bool
end

SerialExecution(;fastmath::Bool=true) = SerialExecution(fastmath)


# TODO: Support other assembly methods for multithreading.
# Matrix storage currently needs `element_csc_assembly()` and vector storage needs
# `direct_vector_assembly()`; both are the defaults.
"""
    CPUThreadedExecution(; workers=Threads.nthreads(), distribution=:contiguous, fastmath=true)

CPU assembly with one reusable workspace per logical worker and atomic global
updates. `workers` is the number of tasks, capped at the number of elements;
it does not change Julia's thread pool. `distribution` is `:contiguous` or
`:cyclic`, specialized through the policy type.

Hardware properties and GPU launch settings belong to their own backend's
policy/descriptor when needed, rather than fields shared by every backend.
"""
struct CPUThreadedExecution{Distribution} <: AbstractAssemblyExecution
    workers::Int
    fastmath::Bool
    function CPUThreadedExecution{D}(workers::Integer,fastmath::Bool=true) where D
        workers > 0 || throw(ArgumentError("workers must be positive"))
        D in (:contiguous,:cyclic) || throw(ArgumentError("distribution must be :contiguous or :cyclic"))
        new{D}(Int(workers),fastmath)
    end
end

CPUThreadedExecution(;workers=Threads.nthreads(),distribution=:contiguous,fastmath::Bool=true) =
    CPUThreadedExecution{distribution}(workers,fastmath)

assembly_uses_workers(::SerialExecution) = false
assembly_uses_workers(::CPUThreadedExecution) = true

"""
    CUDAExecution(; threads=128, opt_level=3, fastmath=true)

Single-GPU matrix assembly. Load CUDA; the matrix storage method must be
`element_csc_assembly()` (the default). Prepared numeric inputs and matrix storage persist
in the assembly cache. Each GPU thread constructs its own scratch space.

The variational form is compiled automatically for the GPU. Geometry and reference
tabulations are prepared on the host and uploaded once during assembly setup.
`opt_level` selects the GPUCompiler optimization level (0–3) before PTX generation.
CUDA assembly requires a CUDA.jl toolchain whose GPUCompiler supports this option.
"""
struct CUDAExecution <: AbstractAssemblyExecution
    threads::Int
    opt_level::Int
    fastmath::Bool
    function CUDAExecution(;threads::Integer=128,opt_level::Integer=3,fastmath::Bool=true)
        threads > 0 || throw(ArgumentError("threads must be positive"))
        0 <= opt_level <= 3 || throw(ArgumentError("opt_level must be between 0 and 3"))
        new(Int(threads),Int(opt_level),fastmath)
    end
end

assembly_uses_workers(::CUDAExecution) = false

assembly_fastmath(::AbstractAssemblyExecution) = true
assembly_fastmath(execution::Union{SerialExecution,CPUThreadedExecution,CUDAExecution}) = execution.fastmath

function validate_matrix_assembly(execution::AbstractAssemblyExecution)
    throw(ArgumentError("No matrix assembly backend for $(typeof(execution)); load its extension"))
end

validate_matrix_assembly(::Union{SerialExecution,CPUThreadedExecution}) = nothing

prepare_matrix_allocation(::Union{SerialExecution,CPUThreadedExecution},allocation) = allocation

finish_matrix_assembly!(::AbstractAssemblyExecution) = nothing
validate_matrix_parameters(::AbstractAssemblyExecution,parameters) = nothing
validate_matrix_update(execution,parameters,cache) = validate_matrix_parameters(execution,parameters)

"""
    AssemblyWorker(id, count, nelements, Val(distribution))

Compact logical worker coordinates passed to a generated assembly kernel.
`id` is one-based and independent of the CPU thread that runs the task. The
caller supplies valid coordinates. All fields are integers; a device launcher
can use Int32 coordinates and its global work-item ID without CPU APIs.
"""
struct AssemblyWorker{Distribution,I<:Integer}
    id::I
    count::I
    nelements::I
end

AssemblyWorker(id::I,count::I,nelements::I,::Val{D}) where {D,I<:Integer} =
    AssemblyWorker{D,I}(id,count,nelements)

@inline function assembly_element_range(worker::AssemblyWorker{:contiguous,I}) where I
    # Give the first remainder workers one extra element. Products stay within
    # nelements, including when a device uses Int32 coordinates.
    width,remainder = divrem(worker.nelements,worker.count)
    before = worker.id-one(I)
    lo = before*width+min(before,remainder)+one(I)
    hi = worker.id*width+min(worker.id,remainder)
    lo:hi
end

@inline assembly_element_range(worker::AssemblyWorker{:cyclic}) =
    worker.id:worker.count:worker.nelements

# Keep generated closures and immutable allocation arguments inside a mutable
# runner, so only its reference crosses the latest-world boundary.
mutable struct AssemblyRunner{F,A}
    prepare::F
    argument::A
end

function (runner::AssemblyRunner)(parameters::P) where P
    kernel = runner.prepare(parameters...)
    kernel(runner.argument)
    nothing
end

# The generated kernel closes over shared inputs. Only marked allocation sites
# belong to its workspace; the factory reuses those arrays between updates.
struct PreparedElementAssembly{K,F}
    kernel::K
    workspace_factory::F
end

@inline reuse_assembly_buffer(::Nothing,slot,allocator,args...) = allocator(args...)

@inline function reuse_assembly_buffer(previous::Tuple,slot,allocator,args...)
    slot > length(previous) && return allocator(args...)
    reuse_assembly_buffer(previous[slot],allocator,args...)
end

@inline function reuse_assembly_buffer(buffer,::typeof(zeros),::Type{T},dims::Integer...) where T
    eltype(buffer) === T && size(buffer) == dims ? buffer : zeros(T,dims...)
end

@inline reuse_assembly_buffer(buffer,allocator,args...) = allocator(args...)
@inline reuse_assembly_workspace(previous,current) = current
@inline reuse_assembly_workspace(previous::Tuple,current::Tuple) = previous === current ? previous : current

mutable struct AssemblyWorkerBatch{K,F,A}
    kernel::K
    factory::F
    argument::A
end

mutable struct AssemblyWorkerRecord{B,C}
    batch::B
    worker::C
    # Host-side state can change type with the parameters. The worker enters
    # a function barrier before preparing scratch and running the kernel.
    workspace::Any
end

mutable struct ThreadedAssemblyRunner{F,A,E}
    prepare::F
    argument::A
    execution::E
    nelements::Int
    nworkers::Int
    # Parameter updates can change scratch types. A function barrier keeps this
    # host-side state out of the specialized element kernel.
    workspaces::Any
    worker_records::Any
    tasks::Vector{Task}
end

function ThreadedAssemblyRunner(prepare,argument,nelements,execution::CPUThreadedExecution)
    nworkers = min(nelements,execution.workers)
    ThreadedAssemblyRunner(prepare,argument,execution,Int(nelements),nworkers,
        nothing,nothing,Vector{Task}(undef,nworkers))
end

assembly_runner(::SerialExecution,prepare,argument,nelements) = AssemblyRunner(prepare,argument)
assembly_runner(execution::CPUThreadedExecution,prepare,argument,nelements) =
    ThreadedAssemblyRunner(prepare,argument,nelements,execution)

@inline previous_assembly_workspace(::Nothing,id) = nothing
@inline previous_assembly_workspace(workspaces,id) = workspaces[id]

@inline assembly_workspace_pool(previous::Vector{W},workspace::W,n) where W = previous
@inline assembly_workspace_pool(previous,workspace::W,n) where W = Vector{W}(undef,n)

function collect_assembly_workspaces(previous,first_workspace::W,records) where W
    workspaces = assembly_workspace_pool(previous,first_workspace,length(records))
    workspaces[1] = first_workspace
    for id in 2:length(records)
        workspaces[id] = records[id].workspace::W
    end
    workspaces
end

function collect_assembly_workspaces(previous::Vector{W},records) where W
    first_workspace = first(records).workspace
    first_workspace isa W || return collect_assembly_workspaces(previous,first_workspace,records)
    for id in eachindex(records)
        previous[id] = records[id].workspace::W
    end
    previous
end

collect_assembly_workspaces(previous,records) =
    collect_assembly_workspaces(previous,first(records).workspace,records)

@inline function assembly_worker_records(records,kernel::K,factory::F,argument::A,previous,
        ::CPUThreadedExecution{D},nelements,nworkers) where {K,F,A,D}
    B = AssemblyWorkerBatch{K,F,A}
    R = AssemblyWorkerRecord{B,AssemblyWorker{D,Int}}
    if records isa Vector{R}
        batch = first(records).batch
        batch.kernel = kernel
        batch.factory = factory
        batch.argument = argument
        return records
    end
    batch = AssemblyWorkerBatch(kernel,factory,argument)
    [AssemblyWorkerRecord(batch,AssemblyWorker(id,nworkers,nelements,Val(D)),
        previous_assembly_workspace(previous,id)) for id in 1:nworkers]
end

function run_element_worker!(record::AssemblyWorkerRecord)
    run_element_worker!(record,record.workspace)
end

function run_element_worker!(record::AssemblyWorkerRecord,previous)
    batch = record.batch
    workspace = batch.factory(previous)
    batch.kernel(batch.argument,workspace,record.worker)
    workspace === previous || (record.workspace = workspace)
    nothing
end

function run_element_workers!(records,tasks)
    @sync for id in eachindex(records)
        record = records[id]
        tasks[id] = Threads.@spawn run_element_worker!($record)
    end
    nothing
end

function (runner::ThreadedAssemblyRunner)(parameters::P) where P
    runner.nworkers == 0 && return nothing
    runner.workspaces = prepare_element_ranges!(runner,parameters,runner.workspaces)
    nothing
end

function prepare_element_ranges!(runner::ThreadedAssemblyRunner,parameters::P,previous) where P
    # Specialize on the cached record type before refreshing shared closures.
    prepare_element_ranges!(runner,parameters,previous,runner.worker_records)
end

function prepare_element_ranges!(runner::ThreadedAssemblyRunner{F,A,CPUThreadedExecution{D}},
        parameters::P,previous,records) where {F,A,D,P}
    prepared = runner.prepare(parameters...)
    if runner.nworkers == 1
        workspace = prepared.workspace_factory(previous_assembly_workspace(previous,1))
        prepared.kernel(runner.argument,workspace,AssemblyWorker(1,1,runner.nelements,Val(D)))
        workspaces = assembly_workspace_pool(previous,workspace,1)
        workspaces[1] = workspace
        return workspaces
    end
    records = assembly_worker_records(records,prepared.kernel,prepared.workspace_factory,
        runner.argument,previous,runner.execution,runner.nelements,runner.nworkers)
    runner.worker_records = records
    run_element_workers!(records,runner.tasks)
    collect_assembly_workspaces(previous,records)
end

# Store scalar results instead of boxing them on return from invokelatest.
mutable struct ScalarAssemblyRunner{F,T}
    prepare::F
    init::T
    result::T
end

ScalarAssemblyRunner(prepare,init) = ScalarAssemblyRunner(prepare,init,init)

function (runner::ScalarAssemblyRunner)(parameters::P) where P
    kernel = runner.prepare(parameters...)
    runner.result = kernel(runner.init)
    nothing
end

function run_assembly!(runner,parameters)
    Base.invokelatest(runner,parameters)
    nothing
end

function run_assembly!(runner::AssemblyRunner{F,A},argument::A,parameters) where {F,A}
    runner.argument = argument
    run_assembly!(runner,parameters)
end

# Reusing a cache with a different output type needs a matching specialization.
function run_assembly!(runner::AssemblyRunner,argument,parameters)
    run_assembly!(AssemblyRunner(runner.prepare,argument),parameters)
end

# Sample on the quadrature points. This is needed to visualize fields

function sample(f,quadrature::AbstractQuadrature;kwargs...)
    x = coordinate_quantity(quadrature)
    fx = f(x)
    index = GT.index(Val(0))
    domain = GT.domain(quadrature)
    opts = QuantityOptions(domain,index)
    term = GT.optimize(GT.term(fx,opts))
    T = typeof(GT.prototype(term))
    nfaces = num_faces(domain)
    face_npoints = num_points_accessor(quadrature)
    ptrs = zeros(Int32,nfaces+1)
    for face in 1:nfaces
        ptrs[face+1] = face_npoints(face)
    end
    length_to_ptrs!(ptrs)
    ndata = ptrs[end]-1
    data = zeros(T,ndata)
    vals = JaggedArray(data,ptrs)
    sample!(f,vals,quadrature;kwargs...)
end

function sample!(f,vals,quadrature::AbstractQuadrature;
        parameters = (),
        reuse = isempty(parameters) ? Val(false) : Val(true),
    )
    params_loop = generate_sample(f,quadrature;parameters)
    runner = AssemblyRunner(params_loop,vals)
    run_assembly!(runner,parameters)
    if val_parameter(reuse)
        vals, runner
    else
        vals
    end
end

function update_sample!(vals,runner;parameters=())
    run_assembly!(runner,vals,parameters)
    vals
end

# 0-forms
function assemble_scalar(integral::Integral;
        parameters = (),
        reuse = isempty(parameters) ? Val(false) : Val(true),
        execution::SerialExecution = SerialExecution(),
    )
    contributions = GT.contributions(integral)
    loops = map(contributions) do contribution
        init = zero(prototype(contribution))
        params_loop = generate_assemble_scalar(execution,contribution;parameters)
        runner = ScalarAssemblyRunner(params_loop,init)
        run_assembly!(runner,parameters)
        runner
    end
    b = sum(runner -> runner.result,loops)
    if val_parameter(reuse)
        b, loops
    else
        b
    end
end

function update_scalar(loops; parameters=())
    sum(loops) do runner
        run_assembly!(runner,parameters)
        runner.result
    end
end

function Base.sum(int::Integral)
    assemble_scalar(int)
end

function assemble_face_contribution!(face_to_val,contribution::DomainContribution;
        parameters = (),
        reuse = isempty(parameters) ? Val(false) : Val(true),
        execution::SerialExecution = SerialExecution(),
    )
    params_loop = generate_assemble_face_contribution(execution,contribution;parameters)
    runner = AssemblyRunner(params_loop,face_to_val)
    run_assembly!(runner,parameters)
    if val_parameter(reuse)
        face_to_val, runner
    else
        face_to_val
    end
end

function update_face_contribution!(face_to_val,runner;parameters=())
    run_assembly!(runner,face_to_val,parameters)
    face_to_val
end

function face_contribution(int::Integral,domain::AbstractDomain;execution::SerialExecution=SerialExecution())
    contribution = GT.contribution(int,domain)
    nfaces = num_faces(GT.domain(contribution))
    T = typeof(prototype(contribution))
    face_to_val = zeros(T,nfaces)
    assemble_face_contribution!(face_to_val,contribution;execution)
end

function face_contribution!(a,int::Integral,domain::AbstractDomain;execution::SerialExecution=SerialExecution())
    contribution = GT.contribution(int,domain)
    assemble_face_contribution!(a,contribution;execution)
end

# 1-forms

"""
    assemble_vector(f, T, space; execution=SerialExecution(), kwargs...)

Assemble a linear form. `execution` selects hardware execution and scheduling;
`assembly_method.vector` independently selects the output storage/update method.
Use `CPUThreadedExecution(workers=n, distribution=:contiguous)` for threaded
direct-vector assembly. `reuse=Val(true)` returns an update cache.
Set `fastmath=false` on the execution policy to disable Julia's `@fastmath`
annotation on the generated assembly code.
"""
function assemble_vector(f,::Type{T},space;
    parameters = (),
    reuse = isempty(parameters) ? Val(false) : Val(true),
    execution::AbstractAssemblyExecution = SerialExecution(),
    kwargs...
    ) where T

    haskey(kwargs,:fastmath) && throw(ArgumentError("set fastmath on the execution policy"))
    threaded = assembly_uses_workers(execution)
    arg = 1
    dv = form_argument_quantity(space,arg)
    integral = f(dv)
    contributions = GT.contributions(integral)
    domains = map(GT.domain,contributions)
    alloc = allocate_vector(T,space,domains...;kwargs...)
    if threaded && !is_direct_vector_allocation(alloc.allocation)
        throw(ArgumentError("threaded vector assembly requires direct_vector_assembly()"))
    end
    optimize_options = if haskey(kwargs, :optimize_options)
        kwargs[:optimize_options]
    else
        nothing
    end
    loops = map(contributions) do contribution
        params_loop = generate_assemble_vector(execution,contribution,space;
            parameters,optimize_options,allocation=alloc)
        runner = assembly_runner(execution,params_loop,alloc,num_faces(GT.domain(contribution)))
        run_assembly!(runner,parameters)
        runner
    end
    b, vector_cache = compress(alloc;reuse=Val(true))
    cache = (;loops,alloc,vector_cache)
    if val_parameter(reuse)
        b, cache
    else
        b
    end
end

function update_vector!(b,cache;parameters=())
    (;loops,alloc,vector_cache) = cache
    reset!(alloc)
    for runner in loops
        run_assembly!(runner,parameters)
    end
    compress!(alloc,b,vector_cache)
    b
end

# 1-forms (dual)

function assemble_values!(f,uh::DiscreteField,space::AbstractSpace;
    parameters = (),
    reuse = isempty(parameters) ? Val(false) : Val(true),
    free_or_dirichlet = FREE_AND_DIRICHLET, 
    location=1,
    )

    ds = dual_operator_quantity(space)
    qty = ds(f)
    params_loop = generate_assemble_values(space;parameters,free_or_dirichlet,location) do index
        opts = QuantityOptions(domain,index)
        term(qty,opts)
    end
    loop! = params_loop(parameters...)
    loop!(uh)
    if val_parameter(reuse)
        uh, params_loop
    else
        uh
    end
end

function update_values!(uh,params_loop;parameters=())
    loop! = params_loop(parameters...)
    loop!(uh)
    uh
end

# 2-forms

"""
    assemble_matrix(f, T, trial_space, test_space; execution=SerialExecution(), kwargs...)

Assemble a bilinear form. `execution` selects hardware execution and scheduling;
`assembly_method.matrix` independently selects the output storage/update method.
The default is `monolithic_assembly(element_csc_assembly())`, which is also
the method required by `CPUThreadedExecution(workers=n)` and `CUDAExecution()`.
Sequential execution remains the default. `reuse=Val(true)` returns a cache for
`update_matrix!`.

For `CUDAExecution()`, the form is compiled for CUDA automatically. Subsequent updates reuse the
prepared numeric inputs and device matrix storage.

Set `fastmath=false` on the execution policy to disable Julia's `@fastmath`
annotation on the generated assembly code. CUDA also applies the setting to device compilation,
including user functions called by the kernel. An update cache reuses the mode
selected at assembly time.
"""
function assemble_matrix(f,::Type{T},trial_space,test_space;
    parameters = (),
    reuse = isempty(parameters) ? Val(false) : Val(true),
    execution::AbstractAssemblyExecution = SerialExecution(),
    optimize_options = nothing,
    kwargs...
    ) where T
    haskey(kwargs,:fastmath) && throw(ArgumentError("set fastmath on the execution policy"))
    validate_matrix_assembly(execution)
    validate_matrix_parameters(execution,parameters)
    threaded = assembly_uses_workers(execution)
    arg_test = 1
    arg_trial = 2
    du = form_argument_quantity(trial_space,arg_trial)
    dv = form_argument_quantity(test_space,arg_test)
    integral = f(du,dv)
    contributions = GT.contributions(integral)
    domains = map(GT.domain,contributions)
    host_alloc = allocate_matrix(T,test_space,trial_space,domains...;kwargs...)
    alloc = prepare_matrix_allocation(execution,host_alloc)
    if threaded && !is_element_csc_allocation(alloc.allocation)
        throw(ArgumentError("threaded matrix assembly requires element_csc_assembly()"))
    end
    loops = map(contributions,ntuple(identity,length(contributions))) do contribution,domain_id
        domain_alloc = matrix_assembly_contribution_allocation(alloc,domain_id)
        prototype = alloc === host_alloc ? domain_alloc : matrix_assembly_contribution_allocation(host_alloc,domain_id)
        params_loop = generate_assemble_matrix(execution,contribution,trial_space,test_space,
            prototype;parameters,optimize_options)
        runner = assembly_runner(execution,params_loop,domain_alloc,num_faces(GT.domain(contribution)))
        run_assembly!(runner,parameters)
        runner
    end
    b, matrix_cache = compress(alloc;reuse=Val(true))
    finish_matrix_assembly!(execution)
    parameter_count = parameters isa Tuple ? length(parameters) : nothing
    cache = (;loops,alloc,matrix_cache,execution,parameter_count)
    if val_parameter(reuse)
        b, cache
    else
        b
    end
end

function update_matrix!(A,cache;parameters=())
    (;loops,alloc,matrix_cache,execution) = cache
    validate_matrix_update(execution,parameters,cache)
    reset!(alloc)
    for runner in loops
        run_assembly!(runner,parameters)
    end
    compress!(alloc,A,matrix_cache)
    finish_matrix_assembly!(execution)
    A
end

function assemble_matrix_with_free_and_dirichlet_columns(f,::Type{T},trial_space,test_space;
    parameters = (),
    reuse = isempty(parameters) ? Val(false) : Val(true),
    free_or_dirichlet = FREE, # Rows
    kwargs...
    ) where T
    # We are doing two calls, this can be optimized
    A,Acache = assemble_matrix(f,T,trial_space,test_space;
        parameters,reuse=Val(true),free_or_dirichlet=(free_or_dirichlet,FREE),kwargs...)
    Ad,Adcache = assemble_matrix(f,T,trial_space,test_space;
        parameters,reuse=Val(true),free_or_dirichlet=(free_or_dirichlet,DIRICHLET),kwargs...)
    cache = (Acache,Adcache)
    if val_parameter(reuse)
        A,Ad,cache
    else
        A,Ad
    end
end

function update_matrix_with_free_and_dirichlet_columns!(A,Ad,cache;parameters)
    (Acache,Adcache) = cache
    update_matrix!(A,Acache;parameters)
    update_matrix!(Ad,Adcache;parameters)
    A,Ad
end

# 1 and 2 forms.

function assemble_matrix_and_vector(a,l,::Type{T},U,V;
    parameters = (),
    reuse = isempty(parameters) ? Val(false) : Val(true),
    kwargs...
    ) where T
    A,matrix_cache = assemble_matrix(a,T,U,V;parameters,reuse=Val(true),kwargs...)
    b,vector_cache = assemble_vector(l,T,V;parameters,reuse=Val(true),kwargs...)
    cache = (;matrix_cache,vector_cache)
    if val_parameter(reuse)
        A,b,cache
    else
        A,b
    end
end

function update_matrix_and_vector!(A,b,cache;parameters)
    (;matrix_cache,vector_cache) = cache
    update_matrix!(A,matrix_cache;parameters)
    update_vector!(b,vector_cache;parameters)
    A,b
end

function assemble_matrix_and_vector_with_free_and_dirichlet_columns(a,l,::Type{T},U,V;
    parameters = (),
    reuse = isempty(parameters) ? Val(false) : Val(true),
    free_or_dirichlet = FREE, #Rows
    kwargs...
    ) where T
    A,Ad,matrix_cache = assemble_matrix_with_free_and_dirichlet_columns(
        a,T,U,V;parameters,reuse=Val(true),free_or_dirichlet,kwargs...)
    b,vector_cache = assemble_vector(l,T,V;parameters,reuse=Val(true),free_or_dirichlet,kwargs...)
    cache = (;matrix_cache,vector_cache)
    if val_parameter(reuse)
        A,Ad,b,cache
    else
        A,Ad,b
    end
end

function update_matrix_and_vector_with_free_and_dirichlet_columns!(A,Ad,b,cache;parameters)
    (;matrix_cache,vector_cache) = cache
    update_matrix_with_free_and_dirichlet_columns!(A,Ad,matrix_cache;parameters)
    update_vector!(b,vector_cache;parameters)
    A,b
end

# Linear problems

function PartitionedSolvers_linear_problem(uhd::DiscreteField,a,l,V=GT.space(uhd);
        assembly_method = (;),
        kwargs...
    )
    U = GT.space(uhd)
    diri_vals = dirichlet_values(uhd)
    method = vector_assembly_method(;assembly_method...)
    xd = vector_from_values(method,diri_vals)
    T = eltype(xd)
    A,Ad,b = assemble_matrix_and_vector_with_free_and_dirichlet_columns(a,l,T,U,V;assembly_method,kwargs...)
    mul!(b,Ad,xd,-1,1)
    x = similar(b,axes(A,2))
    fill!(x,zero(eltype(x)))
    PS.linear_problem(x,A,b)
end

function PartitionedSolvers_linear_problem(::Type{T},U::AbstractSpace,a,l,V=U;
    kwargs...
    ) where T
    A,b = assemble_matrix_and_vector(a,l,T,U,V;kwargs...)
    x = similar(b,axes(A,2))
    fill!(x,zero(eltype(x)))
    PS.linear_problem(x,A,b)
end

# Nonlinear problems

function PartitionedSolvers_nonlinear_problem(uh::DiscreteField,r,j,V=GT.space(uh);
        assembly_method = (;),
        kwargs...
    )
    U = GT.space(uh)
    parameters = (uh,)
    method = vector_assembly_method(;assembly_method...)
    x = vector_from_values(method,free_values(uh))
    T = eltype(x)
    b,residual_cache = assemble_vector(r(uh),T,V;parameters,assembly_method,kwargs...)
    A,jacobian_cache = assemble_matrix(j(uh),T,U,V;parameters,assembly_method,kwargs...)
    workspace = (;uh,residual_cache,jacobian_cache,assembly_method)
    PS.nonlinear_problem(nonlinear_problem_update,x,b,A,workspace)
end

function nonlinear_problem_update(p)
    (;uh,residual_cache,jacobian_cache,assembly_method) = PS.workspace(p)
    x = PS.solution(p)
    solution_field!(uh,x;assembly_method)
    A = PS.jacobian(p)
    b = PS.residual(p)
    parameters = (uh,)
    if b !== nothing && A !== nothing
        update_vector!(b,residual_cache;parameters)
        update_matrix!(A,jacobian_cache;parameters)
    elseif b !== nothing && A === nothing
        update_vector!(b,residual_cache;parameters)
    elseif b === nothing && A !== nothing
        update_matrix!(A,jacobian_cache;parameters)
    end
    p = PS.update(p,solution=x,residual=b,jacobian=A)
end

# Solution fields

function solution_field(U::AbstractSpace,x::AbstractVector;kwargs...)
    T = eltype(x)
    uhd = zero_dirichlet_field(T,U)
    solution_field(uhd,x;kwargs...)
end

function solution_field(uhd::DiscreteField,x::AbstractVector;
        assembly_method = (;)
    )
    diri_vals = dirichlet_values(uhd)
    T = eltype(diri_vals)
    U = GT.space(uhd)
    method = vector_assembly_method(;assembly_method...)
    free_vals = values_from_vector(method,free_dofs(U),x)
    uh = discrete_field(U,free_vals,diri_vals)
end

function solution_field!(uh::DiscreteField,x::AbstractVector;
        assembly_method = (;)
    )
    U = GT.space(uh)
    free_vals = free_values(uh)
    method = vector_assembly_method(;assembly_method...)
    values_from_vector!(method,free_vals,x)
    uh
end

function solution_field(U::AbstractSpace,p::PS.AbstractProblem)
    # TODO add a workspace in LinearProblem
    #(;assembly_method) = PS.workspace(p)
    assembly_method = (;)
    solution_field(U,PS.solution(p);assembly_method)
end

function solution_field(U::DiscreteField,p::PS.AbstractProblem)
    #(;assembly_method) = PS.workspace(p)
    assembly_method = (;)
    solution_field(U,PS.solution(p);assembly_method)
end

function solution_field(U::AbstractSpace,p::PS.AbstractSolver)
    #(;assembly_method) = PS.workspace(PS.problem(p))
    assembly_method = (;)
    solution_field(U,PS.solution(p);assembly_method)
end

function solution_field(U::DiscreteField,p::PS.AbstractSolver)
    #(;assembly_method) = PS.workspace(PS.problem(p))
    assembly_method = (;)
    solution_field(U,PS.solution(p);assembly_method)
end

function solution_field!(uh::DiscreteField,p::PS.AbstractProblem)
    #(;assembly_method) = PS.workspace(p)
    assembly_method = (;)
    solution_field!(uh,PS.solution(p);assembly_method)
end

function solution_field!(uh::DiscreteField,p::PS.AbstractSolver)
    #(;assembly_method) = PS.workspace(PS.problem(p))
    assembly_method = (;)
    solution_field!(uh,PS.solution(p);assembly_method)
end

# Functions to be extended in the extension modules

function SciMLBase_LinearProblem end
function SciMLBase_NonlinearProblem end
function SciMLBase_ODEProblem end


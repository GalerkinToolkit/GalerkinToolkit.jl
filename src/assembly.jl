
function monolithic_vector_assembly_strategy()
    function init(dofs,::Type{T}) where T
        n_total_dofs = length(dofs)
        offsets = [0]
        (;offsets,n_total_dofs,T)
    end
    function init(block_dofs::BRange,::Type{T}) where T
        n_total_dofs = length(block_dofs)
        offsets = blocklasts(block_dofs) .- map(length,blocks(block_dofs))
        (;offsets,n_total_dofs,T)
    end
    function scalar_type(setup)
        setup.T
    end
    function counter(setup)
        0
    end
    function count(n,setup,i,field_i)
        n+1
    end
    function allocate(n,setup)
        Ti = Int32
        T = setup.T
        I = zeros(Ti,n)
        V = zeros(T,n)
        (;I,V)
    end
    function set!(alloc,n,setup,v,i,field)
        n += 1
        alloc.I[n] = i+setup.offsets[field]
        alloc.V[n] = v
        n
    end
    function compress(alloc,setup)
        I = alloc.I
        V = alloc.V
        n_total_dofs = setup.n_total_dofs
        vec = PartitionedArrays.dense_vector(I,V,n_total_dofs)
        cache = nothing
        (vec, cache)
    end
    function compress!(vec,cache,alloc,setup)
        I = alloc.I
        V = alloc.V
        PartitionedArrays.dense_vector!(vec,I,V)
        vec
    end
    (;init,scalar_type,counter,count,allocate,set!,compress,compress!)
end

function monolithic_matrix_assembly_strategy(;matrix_type=nothing)
    function init(dofs_test,dofs_trial,::Type{T}) where T
        n_total_rows = length(dofs_test)
        n_total_cols = length(dofs_trial)
        offsets_rows = [0]
        offsets_cols = [0]
        (;offsets_rows,offsets_cols,n_total_rows,n_total_cols,T)
    end
    function init(dofs_test::BRange,dofs_trial,::Type{T}) where T
        n_total_rows = length(dofs_test)
        n_total_cols = length(dofs_trial)
        offsets_rows = blocklasts(dofs_test) .- map(length,blocks(dofs_test))
        offsets_cols = blocklasts(dofs_trial) .- map(length,blocks(dofs_trial))
        (;offsets_rows,offsets_cols,n_total_rows,n_total_cols,T)
    end
    function scalar_type(setup)
        setup.T
    end
    function counter(setup)
        0
    end
    function count(n,setup,i,j,field_i,field_j)
        n+1
    end
    function allocate(n,setup)
        Ti = Int32
        T = setup.T
        I = zeros(Ti,n)
        J = zeros(Ti,n)
        V = zeros(T,n)
        (;I,J,V)
    end
    function set!(alloc,n,setup,v,i,j,field_i,field_j)
        n += 1
        alloc.I[n] = i+setup.offsets_rows[field_i]
        alloc.J[n] = j+setup.offsets_cols[field_j]
        alloc.V[n] = v
        n
    end
    function compress(alloc,setup)
        I = alloc.I
        J = alloc.J
        V = alloc.V
        n_total_rows = setup.n_total_rows
        n_total_cols = setup.n_total_cols
        T = setup.T
        M = matrix_type === nothing ? SparseMatrixCSC{T,Int} : matrix_type
        A,cache = sparse_matrix(M,I,J,V,n_total_rows,n_total_cols;reuse=Val(true))
        (A, cache)
    end
    function compress!(A,cache,alloc,setup)
        V = alloc.V
        sparse_matrix!(A,V,cache)
        A
    end
    (;init,scalar_type,counter,count,allocate,set!,compress,compress!)
end

function assemble_vector(f,space,::Type{T};kwargs...) where T
    axis = 1
    dv = GT.form_argument(space,axis)
    integral = f(dv)
    assemble_vector(integral,space,T;kwargs...)
end

function assemble_vector(integral::Number,space,::Type{T};
    reuse = Val(false),
    vector_strategy = monolithic_vector_assembly_strategy(),
    ) where T
    @assert integral == 0
    free_dofs = GT.free_dofs(space)
    setup = vector_strategy.init(free_dofs,T)
    counter = vector_strategy.counter(setup)
    alloc = vector_strategy.allocate(counter,setup)
    vec,cache = vector_strategy.compress(alloc,setup)
    if val_parameter(reuse) == false
        vec
    else
        vec, (cache,alloc,setup,vector_strategy)
    end
end

function assemble_vector(integral::Integral,space,::Type{T};
    reuse = Val(false),
    vector_strategy = monolithic_vector_assembly_strategy(),
    ) where T

    setup = vector_strategy.init(GT.free_dofs(space),T)
    state0 = (;space,vector_strategy,setup)
    state1 = assemble_vector_count(integral,state0)
    state2 = assemble_vector_allocate(state1)
    state3 = assemble_vector_fill!(integral,state2)
    b, bcache = assemble_vector_compress(state3)
    cache = (;bcache,state2)
    if val_parameter(reuse) == false
        b
    else
        b, cache
    end
end

function assemble_vector!(f,space,b,cache)
    dim = 1
    dv = GT.form_argument(space,dim)
    integral = f(dv)
    assemble_vector!(integral,b,cache)
end

function assemble_vector!(integral::Integral,b,cache)
    (;bcache,state2) = cache
    state3 = assemble_vector_fill!(integral,state2)
    assemble_vector_compress!(b,bcache,state3)
    b
end

function assemble_vector!(integral::Integral,V,b,cache)
    (;bcache,state2) = cache
    state3 = assemble_vector_fill!(integral,state2)
    assemble_vector_compress!(b,bcache,state3)
    b
end

function assemble_vector_count(integral,state)
    (;space,vector_strategy,setup) = state
    contributions = GT.contributions(integral)
    nfields = GT.num_fields(space)
    field_to_domain = map(field->domain(space,field),1:nfields)
    field_to_sDface_to_dofs = map(field->face_dofs(space,field),1:nfields)
    field_to_Dface_to_sDface = map(inverse_faces,field_to_domain) # TODO faces or inverse_faces??
    field_to_D = map(num_dims,field_to_domain)
    counter = vector_strategy.counter(setup)
    for measure_and_contribution in contributions
        measure, = measure_and_contribution
        domain = GT.domain(measure)
        sface_to_face = faces(domain)
        topo = topology(mesh(domain))
        d = num_dims(domain)
        nsfaces = length(sface_to_face)
        for sface in 1:nsfaces
            face = sface_to_face[sface]
            for field in 1:nfields
                Dface_to_sDface = field_to_Dface_to_sDface[field]
                D = field_to_D[field]
                face_to_Dfaces = face_incidence(topo,d,D)
                sDface_to_dofs = field_to_sDface_to_dofs[field]
                Dfaces = face_to_Dfaces[face]
                for Dface in Dfaces
                    sDface = Dface_to_sDface[Dface]
                    dofs = sDface_to_dofs[sDface]
                    for dof in dofs
                        counter = vector_strategy.count(counter,setup,dof,field)
                    end
                end
            end
        end
    end
    (;counter,state...)
end

function assemble_vector_allocate(state)
    (;counter,vector_strategy,setup) = state
    alloc = vector_strategy.allocate(counter,setup)
    (;alloc,state...)
end

function max_local_dofs(space,field)
    rid_to_reffe = reference_fes(component(space,field))
    map(num_dofs,rid_to_reffe) |> maximum
end

function max_local_dofs(space)
    nfields = num_fields(space)
    map(field->max_local_dofs(space,field),1:nfields) |> maximum
end

function assemble_vector_fill!(integral,state)
    (;space,vector_strategy,alloc,setup) = state
    contributions = GT.contributions(integral)
    nfields = GT.num_fields(space)
    field_to_domain = map(field->domain(space,field),1:nfields)
    field_to_sDface_to_dofs = map(field->face_dofs(space,field),1:nfields)
    field_to_Dface_to_sDface = map(faces,field_to_domain)
    field_to_D = map(num_dims,field_to_domain)
    counter = vector_strategy.counter(setup)
    T = vector_strategy.scalar_type(setup)
    b = zeros(T,max_local_dofs(space))
    for measure_and_contribution in contributions
        measure,contribution = measure_and_contribution
        domain = GT.domain(measure)
        sface_to_face = faces(domain)
        topo = topology(mesh(domain))
        d = num_dims(domain)
        form_arity = 1
        index = GT.generate_index(domain,form_arity)
        t = GT.term(contribution,index)
        term_npoints = GT.num_points(measure)
        expr_qty = t |> expression |> simplify
        expr_npoints = term_npoints(index) |> expression |> simplify
        face = face_index(index,d)
        point = point_index(index)
        axis = 1
        idof = dof_index(index,axis)
        Dface_around = face_around_index(index,axis)
        field = field_index(index,axis)
        s_qty = GT.topological_sort(expr_qty,(face,field,Dface_around,point,idof))
        s_npoints = GT.topological_sort(expr_npoints,(face,))
        expr = quote
            (counter,args,storage) -> begin
                $(unpack_index_storage(index,:storage))
                $(s_qty[1])
                $(s_npoints[1])
                b = args.b
                nsfaces = length(args.sface_to_face)
                for sface in 1:nsfaces
                    $face = args.sface_to_face[sface]
                    $(s_qty[2])
                    npoints = $(s_npoints[2])
                    for $field in 1:args.nfields
                        $(s_qty[3])
                        Dface_to_sDface = args.field_to_Dface_to_sDface[$field]
                        D = args.field_to_D[$field]
                        face_to_Dfaces = face_incidence(args.topo,args.d,D)
                        sDface_to_dofs = args.field_to_sDface_to_dofs[$field]
                        Dfaces = face_to_Dfaces[$face]
                        for ($Dface_around,Dface) in enumerate(Dfaces)
                            $(s_qty[4])
                            fill!(b,zero(eltype(b)))
                            sDface = Dface_to_sDface[Dface]
                            dofs = sDface_to_dofs[sDface]
                            for $point in 1:npoints
                                $(s_qty[5])
                                for ($idof,dof) in enumerate(dofs)
                                    b[$idof] += $(s_qty[6])
                                end
                            end
                            for ($idof,dof) in enumerate(dofs)
                                counter = args.vector_strategy.set!(args.alloc,counter,args.setup,b[$idof],dof,$field)
                            end
                        end
                    end
                end
                counter
            end
        end
        loop! = eval(expr)
        storage = index_storage(index)
        args = (;d,nfields,alloc,vector_strategy,setup,b,topo,sface_to_face,field_to_D,field_to_Dface_to_sDface,field_to_sDface_to_dofs)
        counter = invokelatest(loop!,counter,args,storage)
    end
    state
end

function assemble_vector_compress(state)
    (;alloc,setup,vector_strategy) = state
    b, bcache = vector_strategy.compress(alloc,setup)
    b, bcache
end

function assemble_vector_compress!(b,bcache,state)
    (;alloc,setup,vector_strategy) = state
    vector_strategy.compress!(b,bcache,alloc,setup)
    b
end

function assemble_matrix(f,trial_space,test_space,::Type{T};kwargs...) where T
    test_dim = 1
    trial_dim = 2
    dv = GT.form_argument(test_space,test_dim)
    du = GT.form_argument(trial_space,trial_dim)
    integral = f(du,dv)
    assemble_matrix(integral,trial_space,test_space,T;kwargs...)
end

function assemble_matrix(integral::Integral,trial_space,test_space,::Type{T};
        reuse=false,
        matrix_strategy = monolithic_matrix_assembly_strategy(),
    ) where T
    setup = matrix_strategy.init(free_dofs(test_space),free_dofs(trial_space),T)
    state0 = (;test_space,trial_space,matrix_strategy,setup)
    state1 = assemble_matrix_count(integral,state0)
    state2 = assemble_matrix_allocate(state1)
    state3 = assemble_matrix_fill!(integral,state2)
    A, Acache = assemble_matrix_compress(state3)
    cache = (;Acache,state2)
    if reuse == false
        A
    else
        A, cache
    end
end

function assemble_matrix!(f,trial_space,test_space,A,cache)
    test_dim = 1
    trial_dim = 2
    dv = GT.form_argument(test_space,test_dim)
    du = GT.form_argument(trial_space,trial_dim)
    integral = f(du,dv)
    assemble_matrix!(integral,A,cache)
end

function assemble_matrix!(integral::Integral,A,cache)
    (;Acache,state2) = cache
    state3 = assemble_matrix_fill!(integral,state2)
    assemble_matrix_compress!(A,Acache,state3)
    A
end

function assemble_matrix!(integral::Integral,U,V,A,cache)
    (;Acache,state2) = cache
    state3 = assemble_matrix_fill!(integral,state2)
    assemble_matrix_compress!(A,Acache,state3)
    A
end

function assemble_matrix_count(integral,state)
    (;test_space,trial_space,matrix_strategy,setup) = state
    contributions = GT.contributions(integral)
    test, trial = 1, 2
    axis_to_space = (test_space,trial_space)
    axis_to_nfields = map(num_fields,axis_to_space)
    axis_to_field_to_domain = map(space->map(field->domain(space,field),1:num_fields(space)),axis_to_space)
    axis_to_field_to_sDface_to_dofs = map(space->map(field->face_dofs(space,field),1:num_fields(space)),axis_to_space)
    axis_to_field_to_Dface_to_sDface = map(field_to_domain->map(faces,field_to_domain),axis_to_field_to_domain)
    axis_to_field_to_D = map(field_to_domain->map(num_dims,field_to_domain),axis_to_field_to_domain)
    counter = matrix_strategy.counter(setup)
    for measure_and_contribution in contributions
        measure, = measure_and_contribution
        domain = GT.domain(measure)
        sface_to_face = faces(domain)
        topo = topology(mesh(domain))
        d = num_dims(domain)
        nsfaces = length(sface_to_face)
        for sface in 1:nsfaces
            face = sface_to_face[sface]
            for field_test in 1:axis_to_nfields[test]
                test_Dface_to_sDface = axis_to_field_to_Dface_to_sDface[field_test][test]
                test_D = axis_to_field_to_D[field_test][test]
                test_face_to_Dfaces = face_incidence(topo,d,test_D)
                test_sDface_to_dofs = axis_to_field_to_sDface_to_dofs[field_test][test]
                test_Dfaces = test_face_to_Dfaces[face]
                for field_trial in 1:axis_to_nfields[trial]
                    trial_Dface_to_sDface = axis_to_field_to_Dface_to_sDface[field_trial][trial]
                    trial_D = axis_to_field_to_D[field_trial][trial]
                    trial_face_to_Dfaces = face_incidence(topo,d,trial_D)
                    trial_sDface_to_dofs = axis_to_field_to_sDface_to_dofs[field_trial][trial]
                    trial_Dfaces = trial_face_to_Dfaces[face]
                    for trial_Dface in trial_Dfaces #TODO Typo this should be test
                        trial_sDface = trial_Dface_to_sDface[trial_Dface]
                        trial_dofs = trial_sDface_to_dofs[trial_sDface]
                        for trial_Dface in trial_Dfaces
                            trial_sDface = trial_Dface_to_sDface[trial_Dface]
                            trial_dofs = trial_sDface_to_dofs[trial_sDface]
                            for test_dof in test_dofs
                                for trial_dof in trial_dofs
                                    counter = matrix_strategy.count(counter,setup,test_dof,trial_dof,test_field,trial_field)
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    (;counter,state...)
end

function assemble_matrix_allocate(state)
    (;matrix_strategy,counter,setup) = state
    alloc = matrix_strategy.allocate(counter,setup)
    (;alloc,state...)
end

function assemble_matrix_fill!(integral,state)
    (;test_space,trial_space,fields_test,fields_trial,alloc,glues_test,glues_trial,matrix_strategy,setup) = state
    contributions = GT.contributions(integral)
    num_fields_test = GT.num_fields(test_space)
    num_fields_trial = GT.num_fields(trial_space)
    counter = matrix_strategy.counter(setup)
    function loop(glue_test,glue_trial,field_per_dim,measure,contribution)
        field_test,field_trial = field_per_dim
        sface_to_faces_test, _,sface_to_faces_around_test = target_face(glue_test)
        sface_to_faces_trial, _,sface_to_faces_around_trial = target_face(glue_trial)
        face_to_dofs_test = face_dofs(test_space,field_test)
        face_to_dofs_trial = face_dofs(trial_space,field_trial)
        term_qty = GT.term(contribution)
        term_npoints = GT.num_points(measure)
        sface = :sface
        point = :point
        idof_test = :idof_test
        idof_trial = :idof_trial
        face_around_test = :face_around_test
        face_around_trial = :face_around_trial
        dof_per_dim = (idof_test,idof_trial)
        face_around_per_dim = (face_around_test,face_around_trial)
        index = GT.index(;face=sface,point,field_per_dim,dof_per_dim,face_around_per_dim)
        expr_qty = term_qty(index) |> simplify
        expr_npoints = term_npoints(index) |> simplify
        # TODO merge statements
        s_qty = GT.topological_sort(expr_qty,(sface,face_around_test,face_around_trial,idof_test,idof_trial,point))
        s_npoints = GT.topological_sort(expr_npoints,(sface,))
        expr = quote
            function loop!(counter,args,storage)
                (;sface_to_faces_test,sface_to_faces_around_test,
                 sface_to_faces_trial,sface_to_faces_around_trial,
                 face_to_dofs_test,face_to_dofs_trial,field_test,field_trial,
                 alloc,matrix_strategy,setup) = args
                $(unpack_storage(index.dict,:storage))
                $(s_qty[1])
                $(s_npoints[1])
                nsfaces = length(sface_to_faces_test)
                z = zero(matrix_strategy.scalar_type(setup))
                for $sface in 1:nsfaces
                    $(s_qty[2])
                    npoints = $(s_npoints[2])
                    faces_test = sface_to_faces_test[sface]
                    faces_trial = sface_to_faces_trial[sface]
                    faces_around_test = sface_to_faces_around_test[sface]
                    faces_around_trial = sface_to_faces_around_trial[sface]
                    for ($face_around_test,face_test) in zip(faces_around_test,faces_test)
                        if face_test == 0
                            continue
                        end
                        $(s_qty[3])
                        dofs_test = face_to_dofs_test[face_test]
                        ndofs_test = length(dofs_test)
                        for ($face_around_trial,face_trial) in zip(faces_around_trial,faces_trial)
                            if face_trial == 0
                                continue
                            end
                            $(s_qty[4])
                            dofs_trial = face_to_dofs_trial[face_trial]
                            ndofs_trial = length(dofs_trial)
                            for $idof_test in 1:ndofs_test
                                $(s_qty[5])
                                dof_test = dofs_test[$idof_test]
                                for $idof_trial in 1:ndofs_trial
                                    $(s_qty[6])
                                    dof_trial = dofs_trial[$idof_trial]
                                    v = z
                                    for $point in 1:npoints
                                        v += $(s_qty[7])
                                    end
                                    counter = matrix_strategy.set!(alloc,counter,setup,v,dof_test,dof_trial,field_test,field_trial)
                                end
                            end
                        end
                    end
                end
                counter
            end
        end
        loop! = eval(expr)
        storage = GT.storage(index)
        args = (;sface_to_faces_test,sface_to_faces_around_test,
                 sface_to_faces_trial,sface_to_faces_around_trial,
                 face_to_dofs_test,face_to_dofs_trial,field_test,field_trial,
                 alloc,matrix_strategy,setup)
        counter = invokelatest(loop!,counter,args,storage)
        nothing
    end
    for (measure_and_contribution,field_to_glue_test,field_to_glue_trial) in zip(contributions,glues_test,glues_trial)
        measure, contribution = measure_and_contribution
        for field_test in fields_test
            glue_test = field_to_glue_test[field_test]
            for field_trial in fields_trial
                glue_trial = field_to_glue_trial[field_trial]
                if glue_test === nothing || glue_trial == nothing
                    continue
                end
                field_per_dim = (field_test,field_trial)
                loop(glue_test,glue_trial,field_per_dim,measure,contribution)
            end
        end
    end
    state
end

function assemble_matrix_compress(state)
    (;alloc,matrix_strategy,setup) = state
    A, Acache = matrix_strategy.compress(alloc,setup)
    A, Acache
end

function assemble_matrix_compress!(A,Acache,state)
    (;alloc,matrix_strategy,setup) = state
    matrix_strategy.compress!(A,Acache,alloc,setup)
    A
end

function assemble_matrix_and_vector(a,l,U,V,::Type{T};
        reuse = Val(false),
        matrix_strategy = monolithic_matrix_assembly_strategy(),
        vector_strategy = monolithic_vector_assembly_strategy(),
    ) where T
    A,Acache = assemble_matrix(a,U,V,T;reuse=Val(true),matrix_strategy)
    b,bcache = assemble_vector(l,V,T;reuse=Val(true),vector_strategy)
    cache = (;Acache,bcache)
    if val_parameter(reuse)
        A,b,cache
    else
        A,b
    end
end

function assemble_matrix_and_vector!(a,l,U,V,A,b,cache)
    (;Acache,bcache) = cache
    assemble_matrix!(a,U,V,A,Acache)
    assemble_vector!(l,V,b,bcache)
    A,b
end

function assemble_matrix_and_vector_with_dirichlet(a,l,U,V,dirichlet_values;kwargs...)
    T = eltype(dirichlet_values)
    free_values = constant_values(zero(T),GT.free_dofs(U))
    ud = discrete_field(U,free_values,dirichlet_values)
    l2(v) = l(v) - a(ud,v)
    assemble_matrix_and_vector(a,l2,U,V,T;kwargs...)
end

function assemble_matrix_and_vector_with_dirichlet!(a,l,U,V,dirichlet_values,A,b,cache)
    T = eltype(diri_vals)
    free_values = constant_values(zero(T),GT.free_dofs(U))
    ud = discrete_field(U,free_values,dirichlet_values)
    l2(v) = l(v) - a(ud,v)
    assemble_matrix_and_vector!(a,l2,U,V,A,b,cache)
end

function linear_problem(uhd::DiscreteField,a,l,V=GT.space(uhd);
        matrix_strategy = monolithic_matrix_assembly_strategy(),
        vector_strategy = monolithic_vector_assembly_strategy(),
    )
    U = GT.space(uhd)
    A,b = assemble_matrix_and_vector_with_dirichlet(a,l,U,V,dirichlet_values(uhd);matrix_strategy,vector_strategy)
    x = similar(b,axes(A,2))
    PS.linear_problem(x,A,b)
    #T = eltype(dirichlet_values(uhd))
    #A = assemble_matrix(a,U,V,T)
    #b = assemble_vector(l,V,T)
    #d = assemble_vector(v->a(uhd,v),V,T)
    #b .= b .- d
    #x = similar(b,axes(A,2))
    #PS.linear_problem(x,A,b)
end

function linear_problem(::Type{T},U::AbstractSpace,a,l,V=U;
        matrix_strategy = monolithic_matrix_assembly_strategy(),
        vector_strategy = monolithic_vector_assembly_strategy(),
    ) where T
    A,b = assemble_matrix_and_vector(a,l,U,V,T;matrix_strategy,vector_strategy)
    x = similar(b,axes(A,2))
    PS.linear_problem(x,A,b)
end

function solution_field(U::AbstractSpace,x::AbstractVector)
    T = eltype(x)
    uhd = zero_dirichlet_field(T,U)
    solution_field(uhd,x)
end

function solution_field(uhd::DiscreteField,x::AbstractVector)
    diri_vals = dirichlet_values(uhd)
    U = GT.space(uhd)
    free_vals = free_values_from_solution(x,free_dofs(U))
    discrete_field(U,free_vals,diri_vals)
end

function solution_field(U::AbstractSpace,p::PS.AbstractProblem)
    solution_field(U,PS.solution(p))
end

function solution_field(U::DiscreteField,p::PS.AbstractProblem)
    solution_field(U,PS.solution(p))
end

function solution_field(U::AbstractSpace,p::PS.AbstractSolver)
    solution_field(U,PS.solution(p))
end

function solution_field(U::DiscreteField,p::PS.AbstractSolver)
    solution_field(U,PS.solution(p))
end

function solution_field(uts,p::PS.AbstractODEProblem)
    x = PS.solution(p)
    t = x[1]
    uas = x[2:end]
    uhs = map((ua,ut)->solution_field(ut(t),ua),uas,uts)
    t, uhs
end

function solution_field(uts,p::PS.AbstractODESolver)
    solution_field(uts,PS.problem(p))
end


function free_values_from_solution(x,dofs)
    x
end

function free_values_from_solution(x,dofs::BRange)
    nfields = blocklength(dofs)
    map(1:nfields) do field
        pend = blocklasts(dofs)[field]
        pini = 1 + pend - length(blocks(dofs)[field])
        view(x,pini:pend)
    end |> BVector
end

function free_values_from_solution(x::BVector,dofs::BRange)
    x
end

function nonlinear_problem(uh0::DiscreteField,r,j,V=GT.space(uh0);
        matrix_strategy = monolithic_matrix_assembly_strategy(),
        vector_strategy = monolithic_vector_assembly_strategy(),
    )
    a0 = j(uh0)
    l0 = r(uh0)
    U = GT.space(uh0)
    x0 = free_values(uh0)
    T = eltype(x0)
    A0,b0,cache = assemble_matrix_and_vector(a0,l0,U,V,T;
       reuse=Val(true),matrix_strategy,vector_strategy)
    PS.nonlinear_problem(x0,b0,A0,cache) do p
        x = PS.solution(p)
        # TODO call to consistent! in parallel code
        uh = solution_field(uh0,x)
        a = j(uh)
        l = r(uh)
        A = PS.jacobian(p)
        b = PS.residual(p)
        ws = PS.workspace(p)
        if b !== nothing && A !== nothing
            assemble_matrix_and_vector!(a,l,U,V,A,b,ws)
        elseif b !== nothing && A === nothing
            assemble_vector!(l,V,b,ws.bcache)
        elseif b === nothing && A !== nothing
            assemble_matrix!(a,U,V,A,ws.Acache)
        end
        p = PS.update(p,solution=x,residual=b,jacobian=A)
    end
end

function semi_discrete_field(T,V::AbstractSpace)
    semi_discrete_field(T,V) do t,uh
        uh
    end
end

function semi_discrete_field(f,T,V::AbstractSpace)
    uh = zero_field(T,V)
    semi_discrete_field(f,uh)
end

function semi_discrete_field(uh::DiscreteField)
    semi_discrete_field(uh) do t,uh
        uh
    end
end

function semi_discrete_field(f,uh::DiscreteField)
    SemiDiscreteField(f,uh)
end

struct SemiDiscreteField{A,B}
    update::A
    discrete_field::B
end

function (u::SemiDiscreteField)(t)
    u.update(t,u.discrete_field)
    u.discrete_field
end

function space(u::SemiDiscreteField)
    space(u.discrete_field)
end

function discrete_field(u::SemiDiscreteField)
    u.discrete_field
end

function nonlinear_ode(
    uts,tspan,res,jacs,V=space(uts[1]);
    matrix_strategy = monolithic_matrix_assembly_strategy(),
    vector_strategy = monolithic_vector_assembly_strategy(),
    )

    U = space(uts[1])
    test=1
    trial=2
    v = form_argument(V,test)
    du = form_argument(U,trial)
    t0 = first(tspan)
    uhs_0 = map(uh->uh(t0),uts)
    coeffs_0 = map(uh_0 -> one(eltype(free_values(uh_0))),uhs_0) 
    r0_int = res(t0,uhs_0)(v)
    js0_int = map((jac,coeff)->coeff*jac(t0,uhs_0)(du,v),jacs,coeffs_0)
    j0_int = sum(js0_int)
    T = eltype(free_values(uhs_0[1]))
    A0,b0,cache = assemble_matrix_and_vector(j0_int,r0_int,U,V,T;
       reuse=Val(true),matrix_strategy,vector_strategy)
    x0 = (t0,map(free_values,uhs_0)...)
    PS.ode_problem(x0,b0,A0,tspan,coeffs_0,cache) do p
        coeffs = PS.coefficients(p)
        x = PS.solution(p)
        A = PS.jacobian(p)
        b = PS.residual(p)
        ws = PS.workspace(p)
        t = x[1]
        uhs_t = map((uh,y)->solution_field(uh(t),y),uts,x[2:end])
        r_int = res(t,uhs_t)(v)
        js_int = map((jac,coeff)->coeff*jac(t,uhs_t)(du,v),jacs,coeffs)
        j_int = sum(js_int)
        if b !== nothing && A !== nothing
            assemble_matrix_and_vector!(j_int,r_int,U,V,A,b,ws)
        elseif b !== nothing && A === nothing
            assemble_vector!(r_int,V,b,ws.bcache)
        elseif b === nothing && A !== nothing
            assemble_matrix!(j_int,U,V,A,ws.Acache)
        end
        p = PS.update(p,solution=x,residual=b,jacobian=A)
    end
end


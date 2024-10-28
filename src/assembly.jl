
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
    dim = 1
    dv = GT.shape_functions(space,dim)
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
    dv = GT.shape_functions(space,dim)
    integral = f(dv)
    assemble_vector!(integral,b,cache)
end

function assemble_vector!(integral::Integral,b,cache)
    (;bcache,state2) = cache
    state3 = assemble_vector_fill!(integral,state2)
    assemble_vector_compress!(b,bcache,state3)
    b
end

function assemble_vector_count(integral,state)
    (;space,vector_strategy,setup) = state
    contributions = GT.contributions(integral)
    num_fields = GT.num_fields(space)
    fields = ntuple(identity,num_fields)
    glues = map(contributions) do measure_and_contribution
        measure, = measure_and_contribution
        domain = GT.domain(measure)
        map(fields) do field
            GT.domain_glue(domain,GT.domain(space,field);strict=false)
        end
    end
    counter = vector_strategy.counter(setup)
    function loop(glue::Nothing,field)
    end
    function loop(glue,field)
        sface_to_faces, = target_face(glue)
        face_to_dofs = face_dofs(space,field)
        nsfaces = length(sface_to_faces)
        for sface in 1:nsfaces
            faces = sface_to_faces[sface]
            for face in faces
                if face == 0
                    continue
                end
                dofs = face_to_dofs[face]
                ndofs = length(dofs)
                for idof in 1:ndofs
                    dof = dofs[idof]
                    counter = vector_strategy.count(counter,setup,dof,field)
                end
            end
        end
    end
    for (domain_and_contribution,field_to_glue) in zip(contributions,glues)
        domain, _ = domain_and_contribution
        map(loop,field_to_glue,fields)
    end
    (;counter,glues,fields,state...)
end

function assemble_vector_allocate(state)
    (;counter,vector_strategy,setup) = state
    alloc = vector_strategy.allocate(counter,setup)
    (;alloc,state...)
end

function assemble_vector_fill!(integral,state)
    (;space,glues,fields,vector_strategy,alloc,setup) = state
    contributions = GT.contributions(integral)
    num_fields = GT.num_fields(space)
    counter = vector_strategy.counter(setup)
    function loop(glue::Nothing,field,measure,contribution)
    end
    function loop(glue,field,measure,contribution)
        sface_to_faces,_,sface_to_faces_around = target_face(glue)
        face_to_dofs = face_dofs(space,field)
        term_qty = GT.term(contribution)
        term_npoints = GT.num_points(measure)
        sface = :sface
        point = :point
        idof = :idof
        face_around = :face_around
        field_per_dim = (field,)
        dof_per_dim = (idof,)
        face_around_per_dim = (face_around,)
        index = GT.index(;face=sface,point,field_per_dim,dof_per_dim,face_around_per_dim)
        expr_qty = term_qty(index) |> simplify
        expr_npoints = term_npoints(index) |> simplify
        # TODO merge statements
        s_qty = GT.topological_sort(expr_qty,(sface,face_around,idof,point))
        s_npoints = GT.topological_sort(expr_npoints,(sface,))
        expr = quote
            function loop!(counter,args,storage)
                (;sface_to_faces,sface_to_faces_around,face_to_dofs,field,vector_strategy,alloc,setup) = args
                $(unpack_storage(index.dict,:storage))
                $(s_qty[1])
                $(s_npoints[1])
                nsfaces = length(sface_to_faces)
                z = zero(vector_strategy.scalar_type(setup))
                for $sface in 1:nsfaces
                    $(s_qty[2])
                    npoints = $(s_npoints[2])
                    faces = sface_to_faces[sface]
                    faces_around = sface_to_faces_around[sface]
                    for ($face_around,face) in zip(faces_around,faces)
                        if face == 0
                            continue
                        end
                        $(s_qty[3])
                        dofs = face_to_dofs[face]
                        ndofs = length(dofs)
                        for $idof in 1:ndofs
                            $(s_qty[4])
                            v = z
                            for $point in 1:npoints
                                v += $(s_qty[5])
                            end
                            dof = dofs[$idof]
                            counter = vector_strategy.set!(alloc,counter,setup,v,dof,field)
                        end
                    end
                end
                counter
            end
        end
        loop! = eval(expr)
        args = (;sface_to_faces,sface_to_faces_around,face_to_dofs,field,vector_strategy,alloc,setup)
        storage = GT.storage(index)
        counter = invokelatest(loop!,counter,args,storage)
        nothing
    end
    for (measure_and_contribution,field_to_glue) in zip(contributions,glues)
        measure, contribution = measure_and_contribution
        map(fields,field_to_glue) do field,glue
            loop(glue,field,measure,contribution)
        end
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
    dv = GT.shape_functions(test_space,test_dim)
    du = GT.shape_functions(trial_space,trial_dim)
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
    dv = GT.shape_functions(test_space,test_dim)
    du = GT.shape_functions(trial_space,trial_dim)
    integral = f(du,dv)
    assemble_matrix!(integral,A,cache)
end

function assemble_matrix!(integral::Integral,A,cache)
    (;Acache,state2) = cache
    state3 = assemble_matrix_fill!(integral,state2)
    assemble_matrix_compress!(A,Acache,state3)
    A
end

function assemble_matrix_count(integral,state)
    (;test_space,trial_space,matrix_strategy,setup) = state
    contributions = GT.contributions(integral)
    num_fields_test = GT.num_fields(test_space)
    num_fields_trial = GT.num_fields(trial_space)
    fields_test = ntuple(identity,num_fields_test)
    fields_trial = ntuple(identity,num_fields_trial)
    glues_test = map(contributions) do measure_and_contribution
        measure, = measure_and_contribution
        domain = GT.domain(measure)
        map(fields_test) do field
            GT.domain_glue(domain,GT.domain(test_space,field),strict=false)
        end
    end
    glues_trial = map(contributions) do measure_and_contribution
        measure, = measure_and_contribution
        domain = GT.domain(measure)
        map(fields_trial) do field
            GT.domain_glue(domain,GT.domain(trial_space,field),strict=false)
        end
    end
    counter = matrix_strategy.counter(setup)
    function loop(glue_test,glue_trial,field_per_dim)
        sface_to_faces_test, = target_face(glue_test)
        sface_to_faces_trial, = target_face(glue_trial)
        face_to_dofs_test = face_dofs(test_space,field_per_dim[1])
        face_to_dofs_trial = face_dofs(trial_space,field_per_dim[2])
        nsfaces = length(sface_to_faces_test)
        field_test,field_trial = field_per_dim
        for sface in 1:nsfaces
            faces_test = sface_to_faces_test[sface]
            faces_trial = sface_to_faces_trial[sface]
            for face_test in faces_test
                if face_test == 0
                    continue
                end
                dofs_test = face_to_dofs_test[face_test]
                ndofs_test = length(dofs_test)
                for face_trial in faces_trial
                    if face_trial == 0
                        continue
                    end
                    dofs_trial = face_to_dofs_trial[face_trial]
                    ndofs_trial = length(dofs_trial)
                    for idof_test in 1:ndofs_test
                        dof_test = dofs_test[idof_test]
                        for idof_trial in 1:ndofs_trial
                            dof_trial = dofs_trial[idof_trial]
                            counter = matrix_strategy.count(counter,setup,dof_test,dof_trial,field_test,field_trial)
                        end
                    end
                end
            end
        end
    end
    for (domain_and_contribution,field_to_glue_test,field_to_glue_trial) in zip(contributions,glues_test,glues_trial)
        domain, _ = domain_and_contribution
        for field_test in fields_test
            glue_test = field_to_glue_test[field_test]
            for field_trial in fields_trial
                glue_trial = field_to_glue_trial[field_trial]
                if glue_test === nothing || glue_trial == nothing
                    continue
                end
                field_per_dim = (field_test,field_trial)
                loop(glue_test,glue_trial,field_per_dim)
            end
        end
    end
    (;counter,glues_test,glues_trial,fields_test,fields_trial,state...)
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
    assemble_vector!(b,V,b,bcache)
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



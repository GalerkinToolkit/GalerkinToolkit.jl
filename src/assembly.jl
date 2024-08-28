
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
    vec,cache = vector_strategy.compress(alloc,counter,setup)
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
    state0 = (;integral,space,vector_strategy,setup)
    state1 = assemble_vector_count(state0)
    state2 = assemble_vector_allocate(state1)
    state3 = assemble_vector_fill(state2)
    result, cache = assemble_vector_compress(state3)
    if val_parameter(reuse) == false
        result
    else
        result, cache
    end
end

function assemble_vector_count(state)
    (;integral,space,vector_strategy,setup) = state
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

function assemble_vector_fill(state)
    (;integral,space,glues,fields,vector_strategy,alloc,setup) = state
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
                (;sface_to_faces,sface_to_faces_around,face_to_dofs,vector_strategy,alloc,setup) = args
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
                            counter = vector_strategy.set!(alloc,counter,setup,v,dof,$field)
                        end
                    end
                end
                counter
            end
        end
        loop! = eval(expr)
        args = (;sface_to_faces,sface_to_faces_around,face_to_dofs,vector_strategy,alloc,setup)
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
    vec, cache = vector_strategy.compress(alloc,setup)
    vec, (cache,alloc,setup,vector_strategy)
end

function assemble_matrix(f,trial_space,test_space;kwargs...)
    test_dim = 1
    trial_dim = 2
    dv = GT.shape_functions(test_space,test_dim)
    du = GT.shape_functions(trial_space,trial_dim)
    integral = f(du,dv)
    assemble_matrix(integral,trial_space,test_space;kwargs...)
end

function assemble_matrix(integral::Integral,trial_space,test_space;reuse=false,Ti=Int32,T=Float64)
    state0 = (;integral,test_space,trial_space,Ti,T)
    state1 = assemble_matrix_count(state0)
    state2 = assemble_matrix_allocate(state1)
    state3 = assemble_matrix_fill(state2)
    result, cache = assemble_matrix_compress(state3)
    if reuse == false
        result
    else
        result, cache
    end
end

function assemble_matrix_count(state)
    (;integral,test_space,trial_space) = state
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
    n = 0
    function loop(glue_test,glue_trial,field_per_dim)
        sface_to_faces_test, = target_face(glue_test)
        sface_to_faces_trial, = target_face(glue_trial)
        face_to_dofs_test = face_dofs(test_space,field_per_dim[1])
        face_to_dofs_trial = face_dofs(trial_space,field_per_dim[2])
        nsfaces = length(sface_to_faces_test)
        for sface in 1:nsfaces
            faces_test = sface_to_faces_test[sface]
            faces_trial = sface_to_faces_trial[sface]
            for face_test in faces_test
                dofs_test = face_to_dofs_test[face_test]
                ndofs_test = length(dofs_test)
                for face_trial in faces_trial
                    dofs_trial = face_to_dofs_trial[face_trial]
                    ndofs_trial = length(dofs_trial)
                    n += ndofs_test*ndofs_trial
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
    (;n,glues_test,glues_trial,fields_test,fields_trial,state...)
end

function assemble_matrix_allocate(state)
    (;n,Ti,T) = state
    I = zeros(Ti,n)
    J = zeros(Ti,n)
    V = zeros(T,n)
    (;I,J,V,state...)
end

function assemble_matrix_fill(state)
    (;integral,test_space,trial_space,fields_test,fields_trial,I,J,V,glues_test,glues_trial) = state
    contributions = GT.contributions(integral)
    num_fields_test = GT.num_fields(test_space)
    num_fields_trial = GT.num_fields(trial_space)
    n = 0
    function loop(glue_test,glue_trial,field_per_dim,measure,contribution)
        sface_to_faces_test, _,sface_to_faces_around_test = target_face(glue_test)
        sface_to_faces_trial, _,sface_to_faces_around_trial = target_face(glue_trial)
        face_to_dofs_test = face_dofs(test_space,field_per_dim[1])
        face_to_dofs_trial = face_dofs(trial_space,field_per_dim[2])
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
            function loop!(n,args,storage)
                (;sface_to_faces_test,sface_to_faces_around_test,
                 sface_to_faces_trial,sface_to_faces_around_trial,
                 face_to_dofs_test,face_to_dofs_trial,I,J,V) = args
                $(unpack_storage(index.dict,:storage))
                $(s_qty[1])
                $(s_npoints[1])
                nsfaces = length(sface_to_faces_test)
                z = zero(eltype(V))
                for $sface in 1:nsfaces
                    $(s_qty[2])
                    npoints = $(s_npoints[2])
                    faces_test = sface_to_faces_test[sface]
                    faces_trial = sface_to_faces_trial[sface]
                    faces_around_test = sface_to_faces_around_test[sface]
                    faces_around_trial = sface_to_faces_around_trial[sface]
                    for ($face_around_test,face_test) in zip(faces_around_test,faces_test)
                        $(s_qty[3])
                        dofs_test = face_to_dofs_test[face_test]
                        ndofs_test = length(dofs_test)
                        for ($face_around_trial,face_trial) in zip(faces_around_trial,faces_trial)
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
                                    n += 1
                                    I[n] = dof_test
                                    J[n] = dof_trial
                                    V[n] = v
                                end
                            end
                        end
                    end
                end
                n
            end
        end
        loop! = eval(expr)
        storage = GT.storage(index)
        args = (;sface_to_faces_test,sface_to_faces_around_test,
                 sface_to_faces_trial,sface_to_faces_around_trial,
                 face_to_dofs_test,face_to_dofs_trial,I,J,V)
        n = invokelatest(loop!,n,args,storage)
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
    (;I,J,V,test_space,trial_space) = state
    free_dofs_test = GT.free_dofs(test_space)
    free_dofs_trial = GT.free_dofs(trial_space)
    m = length(free_dofs_test)
    n = length(free_dofs_trial)
    vec = PartitionedArrays.sparse_matrix(I,J,V,m,n)
    cache = (;I,J,V)
    vec, cache
end


function linear_problem(uh,a,l,V=GT.space(uh))
    U = GT.space(uh)
    x = free_values(uh)
    fill!(x,0)
    T = eltype(x)
    A = assemble_matrix(a,U,V;T)
    b = assemble_vector(l,V;T)
    d = assemble_vector(v->a(uh,v),V)
    b .= b .- d
    x,A,b
end


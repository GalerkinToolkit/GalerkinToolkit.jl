
# TODO monothinic matrix and vector for the moment
# add a strategy to build the matrix and vector

function assemble_vector(f,space;kwargs...)
    dim = 1
    dv = gk.shape_functions(space,dim)
    integral = f(dv)
    assemble_vector(integral,space;kwargs...)
end

function assemble_vector(integral::Number,space;reuse=false,Ti=Int32,T=Float64)
    @assert integral == 0
    free_dofs = gk.free_dofs(space)
    n = length(free_dofs)
    zeros(T,n)
end

function assemble_vector(integral::Integral,space;reuse=false,Ti=Int32,T=Float64)
    state0 = (;integral,space,Ti,T)
    state1 = assemble_vector_count(state0)
    state2 = assemble_vector_allocate(state1)
    state3 = assemble_vector_fill(state2)
    result, cache = assemble_vector_compress(state3)
    if reuse == false
        result
    else
        result, cache
    end
end

function assemble_vector_count(state)
    (;integral,space) = state
    contributions = gk.contributions(integral)
    num_fields = gk.num_fields(space)
    fields = ntuple(identity,num_fields)
    glues = map(contributions) do domain_and_contribution
        domain, = domain_and_contribution
        map(fields) do field
            gk.domain_glue(domain,gk.domain(space,field);strict=false)
        end
    end
    n = 0
    function loop(glue::Nothing,field)
    end
    function loop(glue,field)
        sface_to_faces, = target_face(glue)
        dim = 1
        num_face_dofs = gk.num_face_dofs(space,dim,field)
        field_per_dim = (field,)
        for sface in 1:gk.num_faces(gk.domain(glue))
            faces = sface_to_faces[sface]
            for face in faces
                index = gk.index(;face,field_per_dim)
                ndofs = num_face_dofs(index)
                n += ndofs
            end
        end
    end
    for (domain_and_contribution,field_to_glue) in zip(contributions,glues)
        domain, _ = domain_and_contribution
        map(loop,field_to_glue,fields)
    end
    (;n,glues,fields,state...)
end

function assemble_vector_allocate(state)
    (;n,Ti,T) = state
    I = zeros(Ti,n)
    V = zeros(T,n)
    (;I,V,state...)
end

function assemble_vector_fill(state)
    (;integral,space,glues,I,V,fields) = state
    contributions = gk.contributions(integral)
    num_fields = gk.num_fields(space)
    n = 0
    function loop(glue::Nothing,field,contribution)
    end
    function loop(glue,field,contribution)
        sface_to_faces,_,sface_to_faces_around = target_face(glue)
        dim =1
        dof_map = gk.dof_map(space,dim,field)
        num_face_dofs = gk.num_face_dofs(space,dim,field)
        field_per_dim = (field,)
        for sface in 1:gk.num_faces(gk.domain(glue))
            faces = sface_to_faces[sface]
            faces_around = sface_to_faces_around[sface]
            for (face_around,face) in zip(faces_around,faces)
                face_around_per_dim = (face_around,)
                index = gk.index(;face,field_per_dim)
                ndofs = num_face_dofs(index)
                for dof in 1:ndofs
                    n += 1
                    dof_per_dim = (dof,)
                    index2 = gk.index(;face,dof_per_dim,field_per_dim,face_around_per_dim)
                    I[n] = dof_map(index2)
                    index3 = gk.replace_face(index2,sface)
                    vn = gk.term(contribution)(index3)
                    V[n] = vn
                end
            end
        end
    end
    for (domain_and_contribution,field_to_glue) in zip(contributions,glues)
        domain, contribution = domain_and_contribution
        map(fields,field_to_glue) do field,glue
            loop(glue,field,contribution)
        end
    end
    state
end

function assemble_vector_compress(state)
    (;I,V,space) = state
    free_dofs = gk.free_dofs(space)
    n = length(free_dofs)
    vec = PartitionedArrays.dense_vector(I,V,n)
    cache = (;I,V)
    vec, cache
end

function assemble_matrix(f,trial_space,test_space;kwargs...)
    test_dim = 1
    trial_dim = 2
    dv = gk.shape_functions(test_space,test_dim)
    du = gk.shape_functions(trial_space,trial_dim)
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
    contributions = gk.contributions(integral)
    num_fields_test = gk.num_fields(test_space)
    num_fields_trial = gk.num_fields(trial_space)
    fields_test = ntuple(identity,num_fields_test)
    fields_trial = ntuple(identity,num_fields_trial)
    glues_test = map(contributions) do domain_and_contribution
        domain, = domain_and_contribution
        map(fields_test) do field
            gk.domain_glue(domain,gk.domain(test_space,field),strict=false)
        end
    end
    glues_trial = map(contributions) do domain_and_contribution
        domain, = domain_and_contribution
        map(fields_trial) do field
            gk.domain_glue(domain,gk.domain(trial_space,field),strict=false)
        end
    end
    n = 0
    # TODO some cross terms missing
    # TODO a lot of code duplication
    function loop(glue_test,glue_trial,field_per_dim)
        sface_to_faces_test, = target_face(glue_test)
        sface_to_faces_trial, = target_face(glue_trial)
        test_dim = 1
        trial_dim = 2
        num_face_dofs_test = gk.num_face_dofs(test_space,test_dim,field_per_dim[1])
        num_face_dofs_trial = gk.num_face_dofs(trial_space,trial_dim,field_per_dim[2])
        for sface in 1:gk.num_faces(gk.domain(glue_test))
            faces_test = sface_to_faces_test[sface]
            faces_trial = sface_to_faces_trial[sface]
            for face_test in faces_test
                index_test = gk.index(;face=face_test,field_per_dim)
                ndofs_test = num_face_dofs_test(index_test)
                for face_trial in faces_trial
                    index_trial = gk.index(;face=face_trial,field_per_dim)
                    ndofs_trial = num_face_dofs_trial(index_trial)
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
    contributions = gk.contributions(integral)
    num_fields_test = gk.num_fields(test_space)
    num_fields_trial = gk.num_fields(trial_space)
    n = 0
    function loop(glue_test,glue_trial,field_per_dim,contribution)
        sface_to_faces_test, _,sface_to_faces_around_test = target_face(glue_test)
        sface_to_faces_trial, _,sface_to_faces_around_trial = target_face(glue_trial)
        test_dim = 1
        trial_dim = 2
        num_face_dofs_test = gk.num_face_dofs(test_space,test_dim,field_per_dim[1])
        num_face_dofs_trial = gk.num_face_dofs(trial_space,trial_dim,field_per_dim[2])
        dof_map_test = gk.dof_map(test_space,test_dim,field_per_dim[1])
        dof_map_trial = gk.dof_map(trial_space,trial_dim,field_per_dim[2])
        for sface in 1:gk.num_faces(gk.domain(glue_test))
            faces_test = sface_to_faces_test[sface]
            faces_trial = sface_to_faces_trial[sface]
            faces_around_test = sface_to_faces_around_test[sface]
            faces_around_trial = sface_to_faces_around_trial[sface]
            for (face_around_test,face_test) in zip(faces_around_test,faces_test)
                index_test = gk.index(;face=face_test,field_per_dim)
                ndofs_test = num_face_dofs_test(index_test)
                for (face_around_trial,face_trial) in zip(faces_around_trial,faces_trial)
                    face_around_per_dim = (face_around_test,face_around_trial)
                    index_trial = gk.index(;face=face_trial,field_per_dim)
                    ndofs_trial = num_face_dofs_trial(index_trial)
                    for dof_test in 1:ndofs_test
                        for dof_trial in 1:ndofs_trial
                            n += 1
                            dof_per_dim = (dof_test,dof_trial)
                            index2_test = gk.index(;face=face_test,dof_per_dim,field_per_dim,face_around_per_dim)
                            index2_trial = gk.index(;face=face_trial,dof_per_dim,field_per_dim,face_around_per_dim)
                            I[n] = dof_map_test(index2_test)
                            J[n] = dof_map_trial(index2_trial)
                            index3 = gk.replace_face(index2_test,sface)
                            V[n] = gk.term(contribution)(index3)
                        end
                    end
                end
            end
        end
    end
    for (domain_and_contribution,field_to_glue_test,field_to_glue_trial) in zip(contributions,glues_test,glues_trial)
        domain, contribution = domain_and_contribution
        for field_test in fields_test
            glue_test = field_to_glue_test[field_test]
            for field_trial in fields_trial
                glue_trial = field_to_glue_trial[field_trial]
                if glue_test === nothing || glue_trial == nothing
                    continue
                end
                field_per_dim = (field_test,field_trial)
                loop(glue_test,glue_trial,field_per_dim,contribution)
            end
        end
    end
    state
end

function assemble_matrix_compress(state)
    (;I,J,V,test_space,trial_space) = state
    free_dofs_test = gk.free_dofs(test_space)
    free_dofs_trial = gk.free_dofs(trial_space)
    m = length(free_dofs_test)
    n = length(free_dofs_trial)
    vec = PartitionedArrays.sparse_matrix(I,J,V,m,n)
    cache = (;I,J,V)
    vec, cache
end


function linear_problem(uh,a,l,V=gk.space(uh))
    U = gk.space(uh)
    x = free_values(uh)
    fill!(x,0)
    T = eltype(x)
    A = assemble_matrix(a,U,V;T)
    b = assemble_vector(l,V;T)
    d = assemble_vector(v->a(uh,v),V)
    b .= b .- d
    x,A,b
end


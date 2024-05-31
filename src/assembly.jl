function assemble_vector(f,space;kwargs...)
    num_fields = gk.num_fields(space)
    fields = ntuple(identity,num_fields)
    dim = 1
    dv = gk.shape_functions(space,dim)
    integral = f(dv)
    assemble_vector(integral,space;kwargs...)
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
            gk.domain_glue(domain,gk.domain(space,field))
        end
    end
    n = 0
    function loop(glue,field)
        loop(glue,domain_glue_style(glue),field)
    end
    function loop(glue,style::InteriorGlue,field)
        sface_to_face = target_face(glue)
        dim = 1
        num_face_dofs = gk.num_face_dofs(space,dim)
        field_per_dim = (field,)
        for sface in 1:gk.num_faces(gk.domain(glue))
            face = sface_to_face[sface]
            index = gk.index(;face,field_per_dim)
            ndofs = num_face_dofs(index)
            n += ndofs
        end
    end
    function loop(glue,style::CoboundaryGlue,field)
        sface_to_faces, _ = target_face(glue)
        dim = 1
        num_face_dofs = gk.num_face_dofs(space,dim)
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
    function loop(glue,field,contribution)
        loop(glue,domain_glue_style(glue),field,contribution)
    end
    function loop(glue,style::InteriorGlue,field,contribution)
        sface_to_face = target_face(glue)
        field_per_dim = (field,)
        dim = 1
        num_face_dofs = gk.num_face_dofs(space,dim)
        dof_map = gk.dof_map(space,dim)
        for sface in 1:gk.num_faces(gk.domain(glue))
            face = sface_to_face[sface]
            index = gk.index(;face,field_per_dim)
            ndofs = num_face_dofs(index)
            for dof in 1:ndofs
                n += 1
                dof_per_dim = (dof,)
                index2 = gk.index(;face,dof_per_dim,field_per_dim)
                I[n] = dof_map(index2)
                index3 = gk.replace_face(index2,sface)
                V[n] = gk.term(contribution)(index3)
            end
        end
    end
    function loop(glue,style::CoboundaryGlue,field,contribution)
        sface_to_faces, _ = target_face(glue)
        dim =1
        dof_map = gk.dof_map(space,dim)
        num_face_dofs = gk.num_face_dofs(space,dim)
        field_per_dim = (field,)
        for sface in 1:gk.num_faces(gk.domain(glue))
            faces = sface_to_faces[sface]
            for (face_around,face) in enumerate(faces)
                face_around_per_dim = (face_around,)
                index = gk.index(;face,field_per_dim)
                ndofs = num_face_dofs(index)
                for dof in 1:ndofs
                    n += 1
                    dof_per_dim = (dof,)
                    index2 = gk.index(;face,dof_per_dim,field_per_dim,face_around_per_dim)
                    I[n] = dof_map(index2)
                    index3 = gk.replace_face(index2,sface)
                    V[n] = gk.term(contribution)(index3)
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


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
    function aux(glue,sface,field)
        aux(glue,domain_glue_style(glue),sface,field)
    end
    function aux(glue,style::InteriorGlue,sface,field)
        sface_to_face = target_face(glue)
        face = sface_to_face[sface]
        num_face_dofs = gk.num_face_dofs(space,field)
        index = gk.index(;face)
        ndofs = num_face_dofs(index)
        n += ndofs
    end
    function aux(glue,style::CoboundaryGlue,sface,field)
        sface_to_faces, _ = target_face(glue)
        num_face_dofs = gk.num_face_dofs(space,field)
        faces = sface_to_faces[sface]
        for face in faces
            index = gk.index(;face)
            ndofs = num_face_dofs(index)
            n += ndofs
        end
    end
    for (domain_and_contribution,field_to_glue) in zip(contributions,glues)
        domain, _ = domain_and_contribution
        map(fields,field_to_glue) do field,glue
            for face in 1:gk.num_faces(domain)
                aux(glue,face,field)
            end
        end
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
    function aux(glue,sface,field,contribution)
        aux(glue,domain_glue_style(glue),sface,field,contribution)
    end
    function aux(glue,style::InteriorGlue,sface,field,contribution)
        sface_to_face = target_face(glue)
        face = sface_to_face[sface]
        num_face_dofs = gk.num_face_dofs(space,field)
        index = gk.index(;face)
        ndofs = num_face_dofs(index)
        dof_map = gk.dof_map(space,1)
        field_per_dim = (field,)
        for dof in 1:ndofs
            n += 1
            dof_per_dim = (dof,)
            index2 = gk.index(;face,dof_per_dim,field_per_dim)
            I[n] = dof_map(index2)
            index3 = gk.replace_face(index2,sface)
            V[n] = gk.term(contribution)(index3)
        end
    end
    function aux(glue,style::CoboundaryGlue,sface,field,contribution)
        sface_to_faces, _ = target_face(glue)
        dim =1
        dof_map = gk.dof_map(space,dim)
        num_face_dofs = gk.num_face_dofs(space,field)
        faces = sface_to_faces[sface]
        field_per_dim = (field,)
        for (face_around,face) in enumerate(faces)
            face_around_per_dim = (face_around,)
            index = gk.index(;face)
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
    for (domain_and_contribution,field_to_glue) in zip(contributions,glues)
        domain, contribution = domain_and_contribution
        map(fields,field_to_glue) do field,glue
            for face in 1:gk.num_faces(domain)
                aux(glue,face,field,contribution)
            end
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


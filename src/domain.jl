




function assemble_array(assembler,integral,dim_to_space)

    # Setup
    dim_to_nfields = map(gk.num_fields,dim_to_space)
    field_cis = CartesianIndices(dim_to_nfields)
    D = length(dim_to_space)
    dim_to_dim = ntuple(identity,D)
    contributions = gk.contributions(integral)
    glues = map(contributions) do contribution
        domain, = contribution
        map(dim_to_space) do space
            num_fields = gk.num_fields(space)
            fields = ntuple(identity,num_fields)
            map(fields) do field
                gk.domain_glue(domain,gk.domain(space,field))
            end
        end
    end

    # Count
    counter = gk.counter(assembler)
    for (contribution,dim_to_field_to_glue) in zip(contributions,glues)
        domain, domain_contribution = contribution
        for face in 1:gk.num_faces(domain)
            for local_face in  1:gk.num_local_faces(domain,face)
                source_face = (face,local_face)
                for field_ci in field_cis
                    dim_to_field = Tuple(field_ci)
                    dim_to_glue = map(getindex,dim_to_field_to_glue,dim_to_field)
                    dim_to_num_faces_around = map(dim_to_glue) do glue
                        gk.num_faces_around(glue)(source_face)
                    end
                    faces_around_cis = CartesianIndices(dim_to_num_faces_around)
                    for face_around_ci in faces_around_cis
                        dim_to_face_around = Tuple(face_around_ci)
                        dim_to_target_face = map(dim_to_space,dim_to_glue,dim_to_face_around) do space,glue,face_around
                            gk.target_face(glue)(source_face,face_around)
                        end
                        dim_to_num_face_dofs = map(dim_to_space,dim_to_target_face,dim_to_field) do space, target_face, field
                            gk.num_face_dofs(space)(field,target_face)
                        end
                        face_dof_cis = CartesianIndices(dim_to_num_face_dofs)
                        for face_dof_ci in face_dof_cis
                            dim_to_face_dof = Tuple(face_dof_ci)
                            dim_to_face_DOF = map(dim_to_space,dim_to_target_face,dim_to_field,dim_to_face_dof) do space, target_face, field, face_dof
                                gk.face_dofs(space)(field,target_face,face_dof)
                            end
                            index = gk.index(dim_to_field,dim_to_face_around,dim_to_face_DOF)
                            gk.update_counter!(assembler)(counter,index)
                        end
                    end
                end
            end
        end
    end

    # Fill
    array = gk.array(assembler)(counter)
    counter = gk.counter(assembler)
    z = zero_value(assembler)
    for (contribution,dim_to_field_to_glue) in zip(contributions,glues)
        domain, domain_contribution = contribution
        for face in 1:gk.num_faces(domain)
            face_contribution = domain_contribution(face)
            for local_face in  1:gk.num_local_faces(domain,face)
                local_face_contribution = face_contribution(local_face)
                source_face = (face,local_face)
                for field_ci in field_cis
                    dim_to_field = Tuple(field_ci)
                    dim_to_glue = map(getindex,dim_to_field_to_glue,dim_to_field)
                    dim_to_num_faces_around = map(dim_to_glue) do glue
                        gk.num_faces_around(glue)(source_face)
                    end
                    faces_around_cis = CartesianIndices(dim_to_num_faces_around)
                    for face_around_ci in faces_around_cis
                        dim_to_face_around = Tuple(face_around_ci)
                        dim_to_target_face = map(dim_to_space,dim_to_glue,dim_to_face_around) do space,glue,face_around
                            gk.target_face(glue)(source_face,face_around)
                        end
                        dim_to_constraints = map(dim_to_space,dim_to_target_face,dim_to_field) do space, target_face, field
                            gk.constraints(space)(field,target_face)
                        end
                        dim_to_face_dof_to_ref_dofs = map(dim_to_constraints) do constraints
                            JaggedArray(constraints.rowval,constraints.colptr)
                        end
                        dim_to_face_dof_to_coeffs = map(dim_to_constraints) do constraints
                            JaggedArray(constraints.nzval,constraints.colptr)
                        end
                        dim_to_num_face_dofs = map(dim_to_space,dim_to_target_face,dim_to_field) do space, target_face, field
                            gk.num_face_dofs(space)(field,target_face)
                        end
                        face_dof_cis = CartesianIndices(dim_to_num_face_dofs)
                        for face_dof_ci in face_dof_cis
                            dim_to_face_dof = Tuple(face_dof_ci)
                            dim_to_ref_dofs = map(getindex,dim_to_face_dof_to_ref_dofs,dim_to_face_dof)
                            dim_to_coeffs = map(getindex,dim_to_face_dof_to_coeffs,dim_to_face_dof)
                            dim_to_num_ref_dofs = map(length,dim_to_ref_dofs)
                            ref_dofs_cis = CartesianIndices(dim_to_num_slave_face_dofs)
                            value = z
                            for ref_dof_ci in ref_dofs_cis
                                dim_to_ref_dof = Tuple(ref_dof_ci)
                                coeff = prod(dim_to_dim,dim_to_ref_dof) do dim,ref_dof
                                    dim_to_coeffs[dim][ref_dof]
                                end
                                dim_to_ref_dof = prod(dim_to_dim,dim_to_ref_dof) do dim,ref_dof
                                    dim_to_ref_dofs[dim][ref_dof]
                                end
                                index = gk.index(dim_to_field,dim_to_face_around,dim_to_ref_dof)
                                value += coeff*local_face_contribution(index)
                            end
                            dim_to_face_DOF = map(dim_to_space,dim_to_target_face,dim_to_field,dim_to_face_dof) do space, target_face, field, face_dof
                                gk.face_dofs(space)(field,target_face,face_dof)
                            end
                            index = gk.index(dim_to_field,dim_to_face_around,dim_to_face_DOF)
                            counter = gk.update_counter!(assembler)(counter,index)
                            gk.update_array!(assembler)(value,array,counter,index)
                        end
                    end
                end
            end
        end
    end
    array
end

function basis(space,dim,field)
    face_to_refid = something# depends on field and dim
    refid_to_funs = something# depends on field and dim
    function basis_face(face)
        function basis_local_face(local_face)
            funs = 
            function basis_index(index)
                dof = index.dof[dim]
                funs[dof]
            end
        end
    end
end

function call(f,args...)
    f(args...)
end

struct Quantity{T}
    term::T
end

struct Index{T,N,M}
    face::NTuple{T,N}
    field::NTuple{T,M}
    dof::NTuple{T,M}
    face_around::NTuple{T,M}
    point::T
end

function index(;face=(0,),field=(0,),dof=(0,),face_around=(0,),point=0)
    Index(face,field,dof,face_around,point)
end

function call(g,args::Quantity...)
    fs = map(gk.term,args)
    term = index -> g(map(f->f(index),fs)...)
    Quantity(term)
end

function integrate(f,measure)
    q = gk.coordinates(measure)
    w = gk.weights(measure)
    np = gk.num_points(measure)
    f_q = f(q)
    domain_contribution = gk.quantity() do index
        sum(1:np(index)) do point
            @assert index.point === nothing
            index2 = replace_point(index,point)
            w(index2)*f_q(index2)
        end
    end
    domain = gk.domain(measure)
    contribution = domain => domain_contribution
    contributions = (contribution,)
    gk.integral(contributions)
end

function coordinates(measure)
    function q(face)
        refid = face_to_refid[face]
        function q_face(local_face)
            points = refid_to_local_face_to_points[face][local_face]
            function q_face_local_face(point)
                points[point]
            end
        end
    end
end

function weights(measure)
    function q(face)
        refid = face_to_refid[face]
        function q_face(local_face)
            points = refid_to_local_face_to_weights[face][local_face]
            function q_face_local_face(point)
                points[point]
            end
        end
    end
end

function basis(fe)
end

function interpolate(f,space)
    sigma = gk.dual_basis(space)
    sigma(f)
end















for face in 1:gk.num_faces(domain)
    face_integral = integral(face)
    for local_face in  1:gk.num_local_faces(domain,face)
        local_face_integral = face_integral(local_face)
        for ifield in 1:gk.num_fields(V)
            glue_V = glues[ifield]
            for iface_around in 1:gk.num_faces_around(glue_V,face,local_face)
                face_V = gk.target_face(glue_V,face,local_face,iface_around)
                for idof in 1:gk.num_face_dofs(V,face_V)
                    index = gk.index((ifield,),(iface_around,),(idof,))
                    value = local_face_integral(index)
                    k += 1
                    V[k] = value
                end
            end
        end
    end
end

function linear_problem(a,l,U,V,uhd)
    dim_to_space = (V,U)
    D = length(dim_to_space)
    dim_to_dim = ntuple(identity,D)
    dim_to_basis = map(dim_to_space,dim_to_dim) do space,dim
        gk.basis(space,dim)
    end
    a_integral = a(dim_to_basis...)
    a_contributions = gk.contributions(a_integral)
    a_domains = map(first,a_contributions)
    (Imat,Jmat),Vmat = assemble_array_symbolic(a_domains,dim_to_space)
    assemble_array_numeric!(Vmat,a_contributions,dim_to_space)
    l_integral = l(gk.basis(V,1))
    l_contributions = gk.contributions(l_integral)
    l_domains = map(first,l_contributions)
    (Ivec,),Vvec = assemble_array_symbolic(l_domains,(V,))
    assemble_array_numeric!(Vvec,l_contributions,(V,))
    A = pa.sparse_matrix(Imat,Jmat,Vmat,m,n)
    b = pa.dense_vector(Ivec,Vvec,m)
    # TODO impose dirichlet bcs in b
end

function assemble_array_symbolic(domains,dim_to_space)
    dim_to_nfields = map(gk.num_fields,dim_to_space)
    field_cis = CartesianIndices(dim_to_nfields)
    D = length(dim_to_space)
    dim_to_dim = ntuple(identity,D)
    dim_to_basis = map(dim_to_space,dim_to_dim) do space,dim
        gk.basis(space,dim)
    end
end

function assemble_array(assembler,integral,dim_to_space)

    # Setup
    dim_to_nfields = map(gk.num_fields,dim_to_space)
    field_cis = CartesianIndices(dim_to_nfields)
    D = length(dim_to_space)
    dim_to_dim = ntuple(identity,D)
    contributions = gk.contributions(integral)
    glues = map(contributions) do contribution
        domain, = contribution
        map(dim_to_space) do space
            num_fields = gk.num_fields(space)
            fields = ntuple(identity,num_fields)
            map(fields) do field
                gk.domain_glue(domain,gk.domain(space,field))
            end
        end
    end

    # Count
    counter = gk.counter(assembler)
    for (contribution,dim_to_field_to_glue) in zip(contributions,glues)
        domain, domain_contribution = contribution
        for face in 1:gk.num_faces(domain)
            for local_face in  1:gk.num_local_faces(domain,face)
                source_face = (face,local_face)
                for field_ci in field_cis
                    dim_to_field = Tuple(field_ci)
                    dim_to_glue = map(getindex,dim_to_field_to_glue,dim_to_field)
                    dim_to_num_faces_around = map(dim_to_glue) do glue
                        gk.num_faces_around(glue)(source_face)
                    end
                    faces_around_cis = CartesianIndices(dim_to_num_faces_around)
                    for face_around_ci in faces_around_cis
                        dim_to_face_around = Tuple(face_around_ci)
                        dim_to_target_face = map(dim_to_space,dim_to_glue,dim_to_face_around) do space,glue,face_around
                            gk.target_face(glue)(source_face,face_around)
                        end
                        dim_to_num_face_dofs = map(dim_to_space,dim_to_target_face,dim_to_field) do space, target_face, field
                            gk.num_face_dofs(space)(field,target_face)
                        end
                        face_dof_cis = CartesianIndices(dim_to_num_face_dofs)
                        for face_dof_ci in face_dof_cis
                            dim_to_face_dof = Tuple(face_dof_ci)
                            dim_to_face_DOF = map(dim_to_space,dim_to_target_face,dim_to_field,dim_to_face_dof) do space, target_face, field, face_dof
                                gk.face_dofs(space)(field,target_face,face_dof)
                            end
                            index = gk.index(dim_to_field,dim_to_face_around,dim_to_face_DOF)
                            counter = gk.increment(assembler)(counter,index)
                        end
                    end
                end
            end
        end
    end

    # Fill
    array = gk.array(assembler)(counter)
    counter = gk.counter(assembler)
    z = zero_value(assembler)
    for (contribution,dim_to_field_to_glue) in zip(contributions,glues)
        domain, domain_contribution = contribution
        for face in 1:gk.num_faces(domain)
            face_contribution = domain_contribution(face)
            for local_face in  1:gk.num_local_faces(domain,face)
                local_face_contribution = face_contribution(local_face)
                source_face = (face,local_face)
                for field_ci in field_cis
                    dim_to_field = Tuple(field_ci)
                    dim_to_glue = map(getindex,dim_to_field_to_glue,dim_to_field)
                    dim_to_num_faces_around = map(dim_to_glue) do glue
                        gk.num_faces_around(glue)(source_face)
                    end
                    faces_around_cis = CartesianIndices(dim_to_num_faces_around)
                    for face_around_ci in faces_around_cis
                        dim_to_face_around = Tuple(face_around_ci)
                        dim_to_target_face = map(dim_to_space,dim_to_glue,dim_to_face_around) do space,glue,face_around
                            gk.target_face(glue)(source_face,face_around)
                        end
                        dim_to_constraints = map(dim_to_space,dim_to_target_face,dim_to_field) do space, target_face, field
                            gk.constraints(space)(field,target_face)
                        end
                        dim_to_face_dof_to_ref_dofs = map(dim_to_constraints) do constraints
                            JaggedArray(constraints.rowval,constraints.colptr)
                        end
                        dim_to_face_dof_to_coeffs = map(dim_to_constraints) do constraints
                            JaggedArray(constraints.nzval,constraints.colptr)
                        end
                        dim_to_num_face_dofs = map(dim_to_space,dim_to_target_face,dim_to_field) do space, target_face, field
                            gk.num_face_dofs(space)(field,target_face)
                        end
                        face_dof_cis = CartesianIndices(dim_to_num_face_dofs)
                        for face_dof_ci in face_dof_cis
                            dim_to_face_dof = Tuple(face_dof_ci)
                            dim_to_ref_dofs = map(getindex,dim_to_face_dof_to_ref_dofs,dim_to_face_dof)
                            dim_to_coeffs = map(getindex,dim_to_face_dof_to_coeffs,dim_to_face_dof)
                            dim_to_num_ref_dofs = map(length,dim_to_ref_dofs)
                            ref_dofs_cis = CartesianIndices(dim_to_num_slave_face_dofs)
                            value = z
                            for ref_dof_ci in ref_dofs_cis
                                dim_to_ref_dof = Tuple(ref_dof_ci)
                                coeff = prod(dim_to_dim,dim_to_ref_dof) do dim,ref_dof
                                    dim_to_coeffs[dim][ref_dof]
                                end
                                dim_to_ref_dof = prod(dim_to_dim,dim_to_ref_dof) do dim,ref_dof
                                    dim_to_ref_dofs[dim][ref_dof]
                                end
                                index = gk.index(dim_to_field,dim_to_face_around,dim_to_ref_dof)
                                value += coeff*local_face_contribution(index)
                            end
                            dim_to_face_DOF = map(dim_to_space,dim_to_target_face,dim_to_field,dim_to_face_dof) do space, target_face, field, face_dof
                                gk.face_dofs(space)(field,target_face,face_dof)
                            end
                            index = gk.index(dim_to_field,dim_to_face_around,dim_to_face_DOF)
                            counter = gk.increment(assembler)(counter,index)
                            gk.set(assembler)(value,array,counter,index)
                        end
                    end
                end
            end
        end
    end
    array
end

function assemble_array_numeric!(V,contributions,dim_to_space)
    dim_to_nfields = map(gk.num_fields,dim_to_space)
    field_cis = CartesianIndices(dim_to_nfields)
    D = length(dim_to_space)
    dim_to_dim = ntuple(identity,D)
    dim_to_basis = map(dim_to_space,dim_to_dim) do space,dim
        gk.basis(space,dim)
    end
    integral = a(dim_to_basis...)
    for (domain,contribution) in gk.contributions(integral)
        dim_to_field_to_glue = map(dim_to_space) do space
            nfields = gk.num_fields(space)
            fields = ntuple(identity,nfields)
            map(fields) do field
                gk.domain_glue(domain,gk.domain(space,field))
            end
        end
        for face in 1:gk.num_faces(domain)
            face_contribution = contribution(face)
            for local_face in  1:gk.num_local_faces(domain,face)
                local_face_integral = face_contribution(local_face)
                for field_ci in field_cis
                    dim_to_field = Tuple(field_ci)
                    dim_to_glue = map(getindex,dim_to_field_to_glue,dim_to_field)
                    dim_to_nfaces_around = map(dim_to_glue) do glue
                        gk.num_faces_around(glue,face,local_face)
                    end
                    faces_around_cis = CartesianIndices(dim_to_nfaces_around)
                    for face_around_ci in faces_around_cis
                        dim_to_face_around = Tuple(face_around_ci)
                        dim_to_target_face = map(dim_to_space,dim_to_glue,dim_to_face_around) do space,glue,face_around
                            gk.target_face(glue,face,local_face,face_around)
                        end
                        dim_to_nface_dofs = map(gk.num_face_dofs,dim_to_space,dim_to_target_face)
                        face_dof_cis = CartesianIndices(dim_to_nface_dofs)
                        A = zeros(T,dim_to_slave_face_dofs) # TODO only allocate once
                        for face_dof_ci in face_dof_cis
                            dim_to_face_dof = Tuple(face_dof_ci)
                            index = gk.index(dim_to_field,dim_to_face_around,dim_to_face_dof)
                            A[face_dof_ci] = local_face_integral(index)
                        end
                        dim_to_constraints = map(dim_to_space,dim_to_target_face) do space,target_face
                            gk.constraints(space,target_face)
                        end
                        dim_to_master_to_face_dofs = map(dim_to_constraints) do constraints
                            JaggedArray(constraints.rowval,constraints.colptr)
                        end
                        dim_to_master_to_coeffs = map(dim_to_constraints) do constraints
                            JaggedArray(constraints.nzval,constraints.colptr)
                        end
                        dim_to_num_master_face_dofs = map(i->size(i,2),dim_to_constraints)
                        master_face_dof_cis = CartesianIndices(dim_to_num_master_face_dofs)
                        for master_face_dof_ci in master_face_dof_cis
                            dim_to_master_face_dof = Tuple(master_face_dof_ci)
                            dim_to_slave_to_face_dof = map(getindex,dim_to_master_to_face_dofs,dim_to_master_face_dof)
                            dim_to_slave_to_coeff = map(getindex,dim_to_master_to_coeffs,dim_to_master_face_dof)
                            dim_to_num_slave_face_dofs = map(length,dim_to_slave_to_face_dof)
                            slave_dofs_cis = CartesianIndices(dim_to_num_slave_face_dofs)
                            value = zero(T)
                            for slave_dof_ci in slave_dofs_cis
                                dim_to_slave = Tuple(slave_dof_ci)
                                coeff = prod(dim_to_dim,dim_to_slave) do dim,slave
                                    dim_to_slave_to_coeff[dim][slave]
                                end
                                dim_to_face_dof = map(dim_to_dim,dim_to_slave) do dim,slave
                                    dim_to_slave_to_face_dof[dim][slave]
                                end
                                face_dof_ci = CartesianIndex(dim_to_face_dof)
                                value += coeff*A[face_dof_ci]
                            end
                            k += 1
                            V[k] = value
                        end
                    end
                end
            end
        end
    end
    V
end



function assemble_matrix_coo(a,U,V)
    du = gk.basis(U)
    dv = gk.basis(V)
    integral = a(du,dv)
    for (domain,contribution) in gk.contributions(integral)
        nfields_U = gk.num_fields(U)
        nfields_V = gk.num_fields(V)
        glues_U = map(ifield->domain_glue(domain,gk.domain(U,ifield)),1:nfields_U)
        glues_V = map(ifield->domain_glue(domain,gk.domain(V,ifield)),1:nfields_V)
        faces_U
        for face in 1:gk.num_faces(domain)
            face_integral = integral(face)
            for jfield in 1:nfields_U
                for ifield in 1:nfields_V
                    for jface_around in 1:gk.num_

                    face_U = faces_U[jfield][face]
                    face_V = faces_V[ifield][face]
                    for jdof in 1:gk.num_face_dofs(U,face_U)
                        for idof in 1:gk.num_face_dofs(V,face_V)
                            index = gk.index((ifield,jfield),(idof,jdof))
                            A = face_integral(index)
                        end
                    end
                end
            end
        end
    end
end


abstract type AbstractDomain <: GalerkinToolkitDataType end
domain(a::AbstractDomain) = a

struct Domain{A,B,C,D,E} <: AbstractDomain
    mesh::A
    physical_names::B
    face_dims::C
    reference::D
    permuted::E
end

function domain(mesh::AbstractFEMesh;
    physical_names=gk.physical_names(mesh,gk.num_dims(mesh)),
    face_dims = (gk.num_dims(mesh),),
    reference = (Val(false),),
    permuted = Val(true),
    )

    Domain(
           mesh,
           physical_names,
           face_dims,
           reference,
           permuted,
          )
end

function face_ids(dom::AbstractDomain)
    function barrier(nDfaces,ndfaces,Dface_to_dfaces,physical_dfaces)
        Dface_ndfaces_touched = zeros(Int32,nDfaces)
        dface_touched = fill(false,ndfaces)
        dface_touched[physical_dfaces] .= true
        for Dface in 1:nDfaces
            for dface in Dface_to_dfaces[Dface]
                if dface_touched[dface]
                    Dface_ndfaces_touched[Dface] += 1
                end
            end
        end
        DfaceDom_to_Dface = findall(i->i!=0,Dface_ndfaces_touched)
        nDfacesDom = length(DfaceDom_to_Dface)
        DfaceDom_to_ldfaces_ptrs = zeros(Int32,nDfacesDom+1)
        DfaceDom_to_ldfaces_ptrs[2:end] = Dface_ndfaces_touched
        length_to_ptrs!(DfaceDom_to_ldfaces_ptrs)
        ndata = DfaceDom_to_ldfaces_ptrs[end]-1
        DfaceDom_to_ldfaces_data = zeros(Int32,ndata)
        for DfaceDom in 1:nDfacesDom
            Dface = DfaceDfaceDom_to_Dface[DfaceDom]
            for (ldface,dface) in enumerate(Dface_to_dfaces[Dface])
                if dface_touched[dface]
                    p = DfaceDom_to_ldfaces_ptrs[DfaceDom]
                    DfaceDom_to_ldfaces_data[p] = ldface
                    DfaceDom_to_ldfaces_ptrs[DfaceDom] += 1
                end
            end
        end
        rewind_ptrs!(DfaceDom_to_ldfaces_ptrs)
        DfaceDom_to_ldfaces = JaggedArray(DfaceDom_to_ldfaces_data,DfaceDom_to_ldfaces_ptrs)
        (DfaceDfaceDom_to_Dface, DfaceDom_to_ldfaces)
    end
    mesh = dom.mesh
    tag_to_name = dom |> gk.physical_names
    face_dims = dom.face_dims
    if length(face_dims) == 1
        D = face_dims[1]
        Dface_to_tag = zeros(Int,num_faces(mesh,D))
        gk.classify_mesh_faces!(Dface_to_tag,mesh,D,tag_to_name)
        physical_Dfaces = findall(i->i!=0,Dface_to_tag)
        (physical_Dfaces,)
    elseif length(face_dims) == 2
        D = face_dims[1]
        d = face_dims[2]
        topo = topology(mesh)
        Dface_to_dfaces = face_incidence(topo,D,d)
        nDfaces = num_faces(mesh,D)
        ndfaces = num_faces(mesh,d)
        dface_to_tag = zeros(Int,num_faces(mesh,d))
        gk.classify_mesh_faces!(dface_to_tag,mesh,d,tag_to_name)
        physical_dfaces = findall(i->i!=0,dface_to_tag)
        return barrier(nDfaces,ndfaces,Dface_to_dfaces,physical_dfaces)
    else
        error("Case not implemented")
    end
end

abstract type AbstractField end

struct AnalyticalField{A,B} <: AbstractField
    f::A
    domain::B
end

function analytical_field(f,dom)
    AnalyticalField(f,dom)
end

struct PieceWiseField{A} <: AbstractField
    fields::A
end

# It is up to the user that the provided fields
# are on disjoint domains
# Check?
function piecewise_field(fields)
    PieceWiseField(fields)
end

@enum FreeOrDirichlet FREE=1 DIRICHLET=2

abstract type AbstractSpace end

function zero_field(::Type{T},V::AbstractSpace) where T
    fv = zeros(T,free_dofs(V))
    dv = zeros(T,dirichlet_dofs(V))
    discrete_field(V,fv,dv)
end

struct DiscreteField{A,B,C} <: AbstractField
    space::A
    free_values::B
    dirichlet_values::C
end
free_and_dirichlet_values(a::DiscreteField) = (a.free_values,a.dirichlet_values)
free_values(a::DiscreteField) = a.free_values
dirichlet_values(a::DiscreteField) = a.dirichlet_values

function discrete_field(s,fv,dv)
    DiscreteField(s,fv,dv)
end

function interpolate!(f,uh)
    interpolate_free!(f,uh)
    interpolate_dirichlet!(f,uh)
    uh
end

function interpolate!(f,free_or_diri,uh)
    interpolate!(f,free_or_diri,uh,gk.space(uh))
end

function interpolate_free!(f,uh)
    interpolate!(f,FREE,uh)
end

function interpolate_dirichlet!(f,uh)
    interpolate!(f,DIRICHLET,uh)
end

function iso_parametric_space(dom::AbstractDomain;dirichlet_boundary=nothing)
    IsoParametricSpace(dom,dirichlet_boundary)
end

struct IsoParametricSpace{A,B} <: AbstractSpace
    domain::A
    dirichlet_boundary::B
end

function free_and_dirichlet_nodes(V::IsoParametricSpace)
    free_and_dirichlet_nodes(V,V.dirichlet_boundary)
end

function free_and_dirichlet_nodes(V::IsoParametricSpace,dirichlet_boundary::Nothing)
    n = V |> gk.domain |> gk.mesh |> gk.num_nodes
    gk.partition_from_mask(fill(true,n))
end

function free_and_dirichlet_nodes(V::IsoParametricSpace,dirichlet_boundary::AbstractDomain)
    mesh = V |> gk.domain |> gk.mesh
    n = mesh |> gk.num_nodes
    node_to_tag = fill(Int32(0),n)
    tag_to_name = dirichlet_boundary |> gk.physical_names
    gk.classify_mesh_nodes!(node_to_tag,mesh,tag_to_name)
    gk.partition_from_mask(i->i==0,node_to_tag)
end

function free_dofs(V::IsoParametricSpace)
    n_free = length(first(free_and_dirichlet_nodes))
    1:n_free
end

function dirichlet_dofs(V::IsoParametricSpace)
    n_diri = length(last(free_and_dirichlet_nodes))
    1:n_diri
end

function face_dofs(V::IsoParametricSpace)
    face_dofs(V,V.dirichlet_boundary)
end

function face_dofs(V::IsoParametricSpace,dirichlet_boundary::Nothing)
    domain = gk.domain(V)
    @assert length(gk.face_dims(domain)) == 1 "Not implemented yet"
    @assert length(gk.physical_names(domain)) === :all "Not implemented yet"
    d = gk.face_dims(domain)[1]
    face_nodes(gk.mesh(domain),d)
    face_to_nodes = face_nodes(gk.mesh(domain),d)
    face_to_nodes
end

function face_dofs(V::IsoParametricSpace,dirichlet_boundary::AbstractDomain)
    face_to_nodes = JaggedArray(gk.face_dofs(V,nothing))
    free_and_dirichlet_nodes = gk.free_and_dirichlet_nodes(V)
    node_to_free_node = gk.permutation(free_and_dirichlet_nodes)
    n_free = length(first(free_and_dirichlet_nodes))
    function node_to_dof(node)
        free_node = node_to_free_node[node]
        if free_node <= n_free
            return Int(free_node)
        end
        Int(-(free_node-n_free))
    end
    face_to_dofs_data = node_to_dof.(face_to_nodes.data)
    face_to_dofs = JaggedArray(face_to_dofs_data,face_to_nodes.ptrs)
    face_to_dofs
end

function reference_fes(V::IsoParametricSpace)
    domain = gk.domain(V)
    @assert length(gk.face_dims(domain)) == 1 "Not implemented yet"
    d = gk.face_dims(domain)[1]
    # TODO LagrangianMesh face needs to be a AbstractElement
    reference_faces(gk.mesh(domain),d)
end

function face_reference_id(V::IsoParametricSpace)
    domain = gk.domain(V)
    @assert length(gk.face_dims(domain)) == 1 "Not implemented yet"
    d = gk.face_dims(domain)[1]
    face_reference_id(gk.mesh(domain),d)
end

function transformation(V::IsoParametricSpace)
    nothing
end

function constraints(V::IsoParametricSpace)
    nothing
end

function interpolate!(field::AnalyticalField,free_or_diri::FreeOrDirichlet,uh::DiscreteField,V::IsoParametricSpace)
    # TODO
    # @assert gk.domain(field) ⊂ gk.domain(V)
    mesh = V |> gk.domain |> gk.mesh
    n = mesh |> gk.num_nodes
    node_to_tag = fill(Int32(0),n)
    tag_to_name = field |> gk.domain |> gk.physical_names
    gk.classify_mesh_nodes!(node_to_tag,mesh,tag_to_name)
    free_and_dirichlet_nodes = gk.free_and_dirichlet_nodes(V)
    free_node_to_node = free_and_dirichlet_nodes[free_or_diri]
    free_node_to_tag = node_to_tag[free_node_to_node]
    snode_to_free_node = findall(i->i!=0,free_node_to_tag)
    node_to_coord = mesh |> gk.node_coordinates
    @views snode_to_coord = node_to_coord[free_node_to_node[snode_to_free_node]]
    fv = gk.free_and_dirichlet_values(uh)[free_or_diri]
    fv[snode_to_free_node] .= field.f.(snode_to_coord)
    uh
end

function interpolate!(field::PiecewiseField,free_or_diri::FreeOrDirichlet,uh::DiscreteField,V::IsoParametricSpace)
    for ifield in field.fields
         interpolate!(ifield,free_or_diri,uh,V)
    end
    uh
end

abstract type AbstractPoint end

struct Point{A,B} <: AbstractPoint
    domain::A
    reference_points::B
    face_reference_id::C
end

function return_type(f,x)
    typeof(f(x))
end

function sample(field,point)
    data = allocate_sample(field,point)
    sample!(data,field,point)
end

function sample_length(point::AbstractPoint)
    @assert gk.domain(field) == gk.domain(point)
    domain = gk.domain(field)
    face_dims = gk.face_dims(domain)
    @assert length(face_dims) == 1 "Not implemented yet"
    nfaces = gk.num_faces(domain)
    d, = face_dims
    face_to_refid = gk.face_reference_id(point)
    refid_to_refpoints = gk.reference_points(point)
    i = 0
    for face in faces
        refid = face_to_refid[face]
        refpoints = refid_to_refpoints[refid]
        i += length(refpoints)
    end
    i
end

function allocate_sample(field::AnalyticalField,point::AbstractPoint)
    @assert gk.domain(field) == gk.domain(point)
    domain = gk.domain(field)
    face_dims = gk.face_dims(domain)
    @assert length(face_dims) == 1 "Not implemented yet"
    faces, = gk.faces(domain)
    reference, = gk.reference(domain)
    d, = face_dims
    face_to_refid = gk.face_reference_id(point)
    refid_to_refpoints = gk.reference_points(point)
    if reference
        y = first(first(refid_to_refpoints))
    else
        mesh = domain |> gk.mesh
        node_to_x = gk.node_coordinates(mesh)
        Tx = eltype(node_to_x)
        y = zero(Tx)
    end
    T = gk.return_type(field.f,y)
    i = sample_length(point)
    zeros(T,i)
end

function sample!(output,field::AnalyticalField,point::AbstractPoint)
    @assert gk.domain(field) == gk.domain(point)
    domain = gk.domain(field)
    face_dims = gk.face_dims(domain)
    @assert length(face_dims) == 1 "Not implemented yet"
    nfaces = gk.num_faces(domain)
    faces, = gk.faces(domain)
    reference, = gk.reference(domain)
    d, = face_dims
    face_to_refid = gk.face_reference_id(point)
    refid_to_refpoints = gk.reference_points(point)
    if reference
        refid_to_values = map(refid_to_refpoints) do refpoints
            field.f.(refpoints)
        end
        i = 0
        for face in 1:nfaces
            refid = face_to_refid[face]
            values = refid_to_values[refid]
            for value in values
                i += 1
                output[i] = value
            end
        end
    else
        mesh = domain |> gk.mesh
        face_to_refid2 = gk.face_reference_id(mesh,d)
        refid2_to_refface = gk.reference_faces(mesh,d)
        refid2_refid_to_A = map(refid2_to_refface) do refface
            tabulator = gk.tabulator(refface)
            map(refid_to_refpoints) do refpoints
                tabulator(gk.value,refpoints)
            end
        end
        i = 0
        node_to_x = gk.node_coordinates(mesh)
        face_to_nodes = gk.face_nodes(mesh,d)
        Tx = eltype(node_to_x)
        for face in 1:nfaces
            refid = face_to_refid[face]
            refid2 = face_to_refid2[faces[face]]
            A = refid2_refid_to_A[refid2][refid]
            nodes = face_to_nodes[face]
            for p in 1:size(A,1)
                y = zero(Tx)
                for k in 1:length(nodes)
                    node = nodes[k]
                    x = node_to_x[node]
                    y += A[p,k]*x
                end
                i += 1
                value = field.f(y)
                output[i] = value
            end
        end
    end
    output
end

function allocate_sample(field::DiscreteField,point::AbstractPoint)
    space = field |> gk.space
    @assert gk.domain(field) == gk.domain(point)
    domain = gk.domain(field)
    face_dims = gk.face_dims(domain)
    @assert length(face_dims) == 1 "Not implemented yet"
    nfaces = gk.num_faces(domain)
    reference, = gk.reference(domain)
    d, = face_dims
    refpoints = gk.reference_points(point) |> first
    reffe = gk.reference_fes(space) |> first
    tabulator = gk.tabulator(refface)
    A = tabulator(gk.value,refpoints)
    free_dof_to_x = gk.free_values(field)
    Tx = eltype(free_dof_to_x)
    b = zeros(size(A,2),Tx)
    c = A*b
    T = eltype(c)
    i = sample_length(point)
    zeros(T,i)
end

function sample!(data,field::DiscreteField,point::AbstractPoint)
    constr = field |> gk.space |> gk.constraints
    @assert constr === nothing "Not yet implemented"
    trans = field |> gk.space |> gk.transformation
    sample!(data,field,point,trans)
end

function sample!(data,field::DiscreteField,point::AbstractPoint,transformation::Nothing)
    space = field |> gk.space
    @assert gk.domain(field) == gk.domain(point)
    domain = gk.domain(field)
    face_dims = gk.face_dims(domain)
    @assert length(face_dims) == 1 "Not implemented yet"
    nfaces = gk.num_faces(domain)
    reference, = gk.reference(domain)
    d, = face_dims
    face_to_refid = gk.face_reference_id(point)
    refid_to_refpoints = gk.reference_points(point)
    face_to_refid2 = gk.face_reference_id(space)
    refid2_to_reffe = gk.reference_fes(space)
    refid2_refid_to_A = map(refid2_to_reffe) do refface
        tabulator = gk.tabulator(refface)
        map(refid_to_refpoints) do refpoints
            tabulator(gk.value,refpoints)
        end
    end
    free_dof_to_x, diri_dof_to_x = gk.free_and_dirichlet_values(field)
    face_to_dofs = gk.face_dofs(space)
    T = eltype(data)
    i = 0
    for face in 1:nfaces
        refid = face_to_refid[face]
        refid2 = face_to_refid2[face]
        A = refid2_refid_to_A[refid2][refid]
        dofs = face_to_dofs[face]
        for p in 1:size(A,1)
            y = zero(T)
            for k in 1:length(dofs)
                dof = dofs[k]
                if dof > 0
                    x = free_dof_to_x[dof]
                else
                    x = diri_dof_to_x[-dof]
                end
                y += A[p,k]*x
            end
            i += 1
            output[i] = y
        end
    end
    data
end

function sample!(data,field::DiscreteField,point::AbstractPoint,transformation)
    error("Not implemented")
end

struct VisualizationField{A,B,C}
    f::A
    T::Type{B}
    domain::C
end

function allocate_sample(field::VisualizationField,point::AbstractPoint)
    i = sample_length(point)
    zeros(field.T,i)
end

function sample!(data,field::VisualizationField,point::AbstractPoint)
    field.f(data,point)
end

struct VtkPlot{A,B,C}
    mesh::A
    dim::Int
    vizualization_mesh::B
    vtk::C
end

function vtk_plot(f,filename,mesh::AbstractFEMesh,dim;kwargs...)
    vmesh, vglue = visualization_mesh(mesh,dim;kwargs...)
    vtk_grid(filename,gk.vtk_args(vmesh)...) do vtk
        plt = VtkPlot(mesh,dim,(vmesh,vglue),vtk)
        f(plt)
    end
end

function plot!(plt::VtkPlot,field::AbstractField;label)
    @assert length(gk.face_dims(domain)) == 1 "Not implemented yet"
    domain = field |> gk.domain
    faces = domain |> gk.faces
    @assert faces == 1:gk.num_faces(plt.mesh,plt.dim) "Not implemented yet"
    mesh = domain |> gk.mesh
    vmesh, vglue = plt |> gk.visualization_mesh
    face_to_refid = gk.face_reference_id(mesh)
    refid_to_refpoints = vglue.reference_coordinates
    point = gk.point(domain,refid_to_refpoints,face_to_refid[faces])
    data = sample(field,point)
    plt.vtk[label] = data
    plt
end


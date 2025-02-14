

#domain(space::AbstractSpace,field) = domain(space)
mesh(a::AbstractSpace) = GT.mesh(GT.domain(a))
num_dims(a::AbstractSpace) = num_dims(mesh(a))
num_free_dofs(a::AbstractSpace) = length(free_dofs(a))
num_dirichlet_dofs(a::AbstractSpace) = length(dirichlet_dofs(a))

function max_local_dofs(space::AbstractSpace,field)
    rid_to_reffe = reference_spaces(GT.field(space,field))
    map(num_dofs,rid_to_reffe) |> maximum
end

function max_local_dofs(space::AbstractSpace)
    nfields = num_fields(space)
    map(field->max_local_dofs(space,field),1:nfields) |> maximum
end


#function free_dofs(a::AbstractSpace,field)
#    @assert field == 1
#    free_dofs(a)
#end

function free_dofs(V::AbstractSpace)
    if workspace(V) !== nothing
        return workspace(V).free_dofs
    end
    state = generate_workspace(V)
    state.free_dofs
end

#function dirichlet_dofs(a::AbstractSpace,field)
#    @assert field == 1
#    dirichlet_dofs(a)
#end

function dirichlet_dofs(V::AbstractSpace)
    if workspace(V) !== nothing
        return workspace(V).dirichlet_dofs
    end
    state = generate_workspace(V)
    state.dirichlet_dofs
end

workspace(space::AbstractSpace) = nothing

function generate_workspace(space::AbstractSpace)
    state = generate_dof_ids(space)
    face_dofs = state.Dface_to_dofs
    free_dofs = state.free_dofs
    dirichlet_dofs = state.dirichlet_dofs
    dirichlet_dof_location = state.dirichlet_dof_location
    workspace = (;face_dofs,free_dofs,dirichlet_dofs,dirichlet_dof_location)
end

function setup_space(space::AbstractSpace)
    if GT.workspace(space) !== nothing
        return space
    end
    workspace = generate_workspace(space)
    replace_workspace(space,workspace)
end

function face_dofs(space::AbstractSpace)
    if workspace(space) !== nothing
        return workspace(space).face_dofs
    end
    state = generate_workspace(space)
    state.face_dofs
end

function free_dof_local_indices(space::AbstractSpace)
    if workspace(space) !== nothing
        return workspace(space).free_dof_local_indices
    end
    state = generate_workspace(space)
    state.free_dof_local_indices
end

function dirichlet_dof_local_indices(space::AbstractSpace)
    if workspace(space) !== nothing
        return workspace(space).dirichlet_dof_local_indices
    end
    state = generate_workspace(space)
    state.dirichlet_dof_local_indices
end

function face_dofs(pspace::AbstractSpace{<:AbstractPMesh})
    p_space = partition(pspace)
    values = map(p_space) do space
        face_to_dofs = face_dofs(space)
        face_to_gdofs = JaggedArray(copy(face_to_dofs))
        free_global = local_to_global(free_dof_local_indices(space))
        diri_global = local_to_global(dirichlet_dof_local_indices(space))
        f = dof -> begin
            if dof > 0
                free_global[dof]
            else
                -diri_global[-dof]
            end
        end
        data = face_to_gdofs.data
        data .= f.(data)
        face_to_gdofs
    end
    # TODO this assumes that all faces are active
    # Eventually, this will be correct once we change the meaning of face_dofs
    pmesh = mesh(pspace)
    d = num_dims(domain(pspace))
    ids = face_partition(pmesh,d)
    PVector(values,ids)
end

#function face_dofs(space::AbstractSpace,field)
#    @assert field == 1
#    face_dofs(space)
#end

#function free_and_dirichlet_dofs(V::AbstractSpace)
#    if workspace(V) !== nothing
#        return workspace(V).free_and_dirichlet_dofs
#    end
#    state = generate_workspace(V)
#    state.free_and_dirichlet_dofs
#end

function dirichlet_dof_location(V::AbstractSpace)
    if workspace(V) !== nothing
        return workspace(V).dirichlet_dof_location
    end
    state = generate_workspace(V)
    state.dirichlet_dof_location
end

function num_fields(a::AbstractSpace)
    1
end

function field(a::AbstractSpace,fieldid)
    @assert fieldid == 1
    a
end

function dofs(a::AbstractSpace,free_or_diri::FreeOrDirichlet)
    if free_or_diri == FREE
        free_dofs(a)
    else
        dirichlet_dofs(a)
    end
end

function values(a::AbstractSpace,free_or_diri::FreeOrDirichlet)
    if free_or_diri == FREE
        free_values(a)
    else
        dirichlet_values(a)
    end
end

function push_forward(a::AbstractSpace,qty)
    qty
end

function pull_back(a::AbstractSpace,qty)
    qty
end

function shape_function_quantity(a::AbstractSpace,dof;reference=is_reference_domain(GT.domain(a)))
    face_to_refid = face_reference_id(a)
    refid_to_reffes = reference_spaces(a)
    refid_to_funs = map(GT.shape_functions,refid_to_reffes)
    domain = GT.domain(a)
    qty = shape_function_quantity(refid_to_funs,domain;reference=true,face_reference_id=face_to_refid,dof)
    if reference
        return qty
    end
    D = num_dims(domain)
    ϕinv = GT.inverse_physical_map(mesh(domain),D)
    push_forward(a,qty)∘ϕinv
end

function form_argument_quantity(a::AbstractSpace,axis,field=1)
    face_to_refid = face_reference_id(a)
    refid_to_reffes = reference_spaces(a)
    refid_to_funs = map(GT.shape_functions,refid_to_reffes)
    domain = GT.domain(a)
    qty = form_argument(axis,field,refid_to_funs,domain;reference=true)# TODO,face_reference_id=face_to_refid)
    if is_reference_domain(GT.domain(a))
        return qty
    end
    D = num_dims(domain)
    ϕinv = GT.inverse_physical_map(mesh(domain),D)
    push_forward(a,qty)∘ϕinv
end

function free_or_dirichlet_value(free_vals,diri_vals,dof)
    if dof < 0
        diri_vals[-dof]
    else
        free_vals[dof]
    end
end

function discrete_field_quantity(a::AbstractSpace,free_vals,diri_vals)
    ldof = gensym("fe-dof")
    s = shape_function_quantity(a,ldof;reference=true)
    # TODO ugly
    if is_physical_domain(GT.domain(a))
        s = push_forward(a,s)
    end
    face_to_dofs = GT.face_dofs(a)
    D = num_dims(domain(a))
    qty = quantity() do index
        face_to_dofs_sym = get_symbol!(index,face_to_dofs,"face_to_dofs")
        free_vals_sym = get_symbol!(index,free_vals,"free_vals")
        diri_vals_sym = get_symbol!(index,diri_vals,"diri_vals")
        face = face_index(index,D)
        expr = @term begin
            dof = $face_to_dofs_sym[$face][$ldof]
            GalerkinToolkit.free_or_dirichlet_value($free_vals_sym,$diri_vals_sym,dof)
        end
        p = zero(eltype(free_vals))
        coeffs = expr_term([D],expr,p,index)
        expr = @term begin
            length($face_to_dofs_sym[$face])
        end
        ndofs = expr_term([D],expr,0,index)
        shapes = term(s,index)
        discrete_function_term(coeffs,shapes,ldof,ndofs)
    end
    # TODO ugly
    if is_reference_domain(GT.domain(a))
        return qty
    end
    ϕinv = GT.inverse_physical_map(mesh(domain(a)),D)
    qty∘ϕinv
end

function discrete_field(space::AbstractSpace,free_values,dirichlet_values)
    qty = discrete_field_quantity(space,free_values,dirichlet_values)
    mesh = space |> GT.mesh
    DiscreteField(mesh,space,free_values,dirichlet_values,qty)
end

function dual_basis_quantity(a::AbstractSpace,dof)
    face_to_refid = face_reference_id(a)
    refid_to_reffes = reference_spaces(a)
    refid_to_funs = map(GT.dual_basis,refid_to_reffes)
    domain = GT.domain(a)
    D = num_dims(domain)
    prototype = first(first(refid_to_funs))
    s = GT.quantity() do index
        face_to_sface = get_symbol!(index,inverse_faces(domain),"face_to_sface")
        sface_to_refid_sym = get_symbol!(index,face_to_refid,"sface_to_refid")
        refid_to_funs_sym = get_symbol!(index,refid_to_funs,"refid_to_funs")
        face = face_index(index,D)
        expr = @term begin
            sface = $face_to_sface[$face]
            refid = $sface_to_refid_sym[sface]
            $refid_to_funs_sym[refid][$dof]
        end
        expr_term([D],expr,prototype,index;dof)
    end
    if is_reference_domain(GT.domain(a))
        return s
    end
    ϕ = GT.physical_map(mesh(domain),D)
    s2 = pull_back(a,s)
    call(dual_compose,s2,ϕ)
end

function dual_compose(s,g)
    f -> s(f∘g)
end

function get_symbol!(index,::typeof(dual_compose),name="";prefix=index.data.prefix)
    :(GalerkinToolkit.dual_compose)
end

function generate_dof_ids(space::AbstractSpace)
    state = generate_dof_ids_step_1(space)
    generate_dof_ids_step_2(space,state,space |> GT.dirichlet_boundary)
end

function generate_dof_ids_step_1(space)
    domain = space |> GT.domain
    D = GT.num_dims(domain)
    cell_to_Dface = domain |> GT.faces
    mesh = domain |> GT.mesh
    topology = mesh |> GT.topology
    d_to_ndfaces = map(d->GT.num_faces(topology,d),0:D)
    ctype_to_reference_fe = space |> GT.reference_spaces
    cell_to_ctype = space |> GT.face_reference_id
    d_to_dface_to_dof_offset = map(d->zeros(Int32,GT.num_faces(topology,d)),0:D)
    d_to_ctype_to_ldface_to_own_dofs = map(d->GT.reference_face_own_dofs(space,d),0:D)
    d_to_ctype_to_ldface_to_pindex_to_perm = map(d->GT.reference_face_own_dof_permutations(space,d),0:D)
    d_to_ctype_to_ldface_to_num_own_dofs = map(d->map(ldface_to_own_dofs->length.(ldface_to_own_dofs),d_to_ctype_to_ldface_to_own_dofs[d+1]),0:D)
    d_to_ctype_to_ldface_to_dofs = map(d->map(fe->GT.face_dofs(fe,d),ctype_to_reference_fe),0:D)
    #d_to_ctype_to_ldface_to_pindex_to_perm = map(d->map(fe->GT.face_own_dof_permutations(fe,d),ctype_to_reference_fe),0:D)
    d_to_Dface_to_dfaces = map(d->face_incidence(topology,D,d),0:D)
    d_to_Dface_to_ldface_to_pindex = map(d->face_permutation_ids(topology,D,d),0:D)
    ctype_to_num_dofs = map(GT.num_dofs,ctype_to_reference_fe)
    ncells = length(cell_to_ctype)
    nDfaces = num_faces(topology,D)
    dof_offset = 0
    for d in 0:D
        ctype_to_ldface_to_num_own_dofs = d_to_ctype_to_ldface_to_num_own_dofs[d+1]
        dface_to_dof_offset = d_to_dface_to_dof_offset[d+1]
        Dface_to_dfaces = d_to_Dface_to_dfaces[d+1]
        ndfaces = length(dface_to_dof_offset)
        for cell in 1:ncells
            ctype = cell_to_ctype[cell]
            Dface = cell_to_Dface[cell]
            ldface_to_num_own_dofs = ctype_to_ldface_to_num_own_dofs[ctype]
            ldface_to_dface = Dface_to_dfaces[Dface]
            nldfaces = length(ldface_to_num_own_dofs)
            for ldface in 1:nldfaces
                num_own_dofs = ldface_to_num_own_dofs[ldface]
                dface = ldface_to_dface[ldface]
                dface_to_dof_offset[dface] = num_own_dofs
            end
        end
        for dface in 1:ndfaces
            num_own_dofs = dface_to_dof_offset[dface]
            dface_to_dof_offset[dface] = dof_offset
            dof_offset += num_own_dofs
        end
    end
    ndofs = dof_offset
    Dface_to_ptrs = zeros(Int32,nDfaces+1)
    for cell in 1:ncells
        ctype = cell_to_ctype[cell]
        num_dofs = ctype_to_num_dofs[ctype]
        Dface = cell_to_Dface[cell]
        Dface_to_ptrs[Dface+1] = num_dofs
    end
    length_to_ptrs!(Dface_to_ptrs)
    ndata = Dface_to_ptrs[end]-1
    Dface_to_dofs = JaggedArray(zeros(Int32,ndata),Dface_to_ptrs)
    # TODO we assume non oriented
    # Optimize for the oriented case?
    # TODO add different numbering strategies
    for d in 0:D
        Dface_to_dfaces = d_to_Dface_to_dfaces[d+1]
        ctype_to_ldface_to_own_ldofs = d_to_ctype_to_ldface_to_own_dofs[d+1]
        ctype_to_ldface_to_pindex_to_perm = d_to_ctype_to_ldface_to_pindex_to_perm[d+1]
        dface_to_dof_offset = d_to_dface_to_dof_offset[d+1]
        Dface_to_ldface_to_pindex = d_to_Dface_to_ldface_to_pindex[d+1]
        for cell in 1:ncells
            ctype = cell_to_ctype[cell]
            Dface = cell_to_Dface[cell]
            ldof_to_dof = Dface_to_dofs[Dface]
            ldface_to_dface = Dface_to_dfaces[Dface]
            ldface_to_own_ldofs = ctype_to_ldface_to_own_ldofs[ctype]
            ldface_to_pindex_to_perm = ctype_to_ldface_to_pindex_to_perm[ctype]
            ldface_to_pindex = Dface_to_ldface_to_pindex[Dface]
            nldfaces = length(ldface_to_dface)
            for ldface in 1:nldfaces
                dface = ldface_to_dface[ldface]
                own_ldofs = ldface_to_own_ldofs[ldface]
                dof_offset = dface_to_dof_offset[dface]
                pindex_to_perm = ldface_to_pindex_to_perm[ldface]
                pindex = ldface_to_pindex[ldface]
                perm = pindex_to_perm[pindex]
                n_own_dofs = length(own_ldofs)
                for i in 1:n_own_dofs
                    j = perm[i]
                    dof = j + dof_offset
                    own_dof = own_ldofs[i]
                    ldof_to_dof[own_dof] = dof
                end
            end
        end
    end
    dof_local_indices = PartitionedArrays.block_with_constant_size(1,(1,),(ndofs,))
    (;ndofs,Dface_to_dofs,d_to_Dface_to_dfaces,
     d_to_ctype_to_ldface_to_dofs,d_to_ndfaces,
     cell_to_ctype,cell_to_Dface,dof_local_indices)
end

function generate_dof_ids_step_2(space,state,dirichlet_boundary::Nothing)
    (;ndofs,Dface_to_dofs) = state
    dof_to_tag = zeros(Int32,ndofs)
    dirichlet_dof_location = zeros(Int32,0)
    free_dofs = Base.OneTo(ndofs)
    dirichlet_dofs = Base.OneTo(0)
    (;Dface_to_dofs, free_dofs, dirichlet_dofs, dirichlet_dof_location)
end

function generate_dof_ids_step_2(space,state,dirichlet_boundary::AbstractDomain)
    (;ndofs,Dface_to_dofs,d_to_Dface_to_dfaces,
     d_to_ctype_to_ldface_to_dofs,
     d_to_ndfaces,cell_to_ctype,cell_to_Dface) = state
    dof_to_tag = zeros(Int32,ndofs)
    N = GT.num_dims(dirichlet_boundary)
    #physical_names = dirichlet_boundary |> GT.physical_names
    Nface_to_tag = zeros(Int32,d_to_ndfaces[N+1])
    mesh = dirichlet_boundary |> GT.mesh
    #classify_mesh_faces!(Nface_to_tag,mesh,N,physical_names)
    Nface_to_tag[GT.faces(dirichlet_boundary)] .= 1
    ncells = length(cell_to_Dface)
    let d = N
        Dface_to_dfaces = d_to_Dface_to_dfaces[d+1]
        ctype_to_ldface_to_ldofs = d_to_ctype_to_ldface_to_dofs[d+1]
        for cell in 1:ncells
            ctype = cell_to_ctype[cell]
            Dface = cell_to_Dface[cell]
            ldof_to_dof = Dface_to_dofs[Dface]
            ldface_to_dface = Dface_to_dfaces[Dface]
            ldface_to_ldofs = ctype_to_ldface_to_ldofs[ctype]
            nldfaces = length(ldface_to_dface)
            dofs = Dface_to_dofs[Dface]
            for ldface in 1:nldfaces
                ldofs = ldface_to_ldofs[ldface]
                dface = ldface_to_dface[ldface]
                Nface = dface
                tag = Nface_to_tag[Nface]
                if tag == 0
                    continue
                end
                dof_to_tag[view(dofs,ldofs)] .= tag
            end
        end
    end
    free_and_dirichlet_dofs = GT.partition_from_mask(i->i==0,dof_to_tag)
    dof_permutation = GT.permutation(free_and_dirichlet_dofs)
    n_free_dofs = length(first(free_and_dirichlet_dofs))
    f = dof -> begin
        dof2 = dof_permutation[dof]
        T = typeof(dof2)
        if dof2 > n_free_dofs
            return T(n_free_dofs-dof2)
        end
        dof2
    end
    data = Dface_to_dofs.data
    data .= f.(data)
    ndiri = length(last(free_and_dirichlet_dofs))
    dirichlet_dof_location = ones(Int32,ndiri)
    free_dofs = Base.OneTo(n_free_dofs)
    dirichlet_dofs = Base.OneTo(ndiri)
    (;Dface_to_dofs, free_dofs, dirichlet_dofs,dirichlet_dof_location)
end

function generate_dof_ids_step_2(space,state,q::AbstractField)
    (;ndofs,Dface_to_dofs,d_to_Dface_to_dfaces,
     d_to_ctype_to_ldface_to_dofs,
     d_to_ndfaces,cell_to_ctype,cell_to_Dface) = state
    dirichlet_boundary = domain(q)
    dof_to_tag = zeros(Int32,ndofs)
    N = GT.num_dims(dirichlet_boundary)
    #physical_names = dirichlet_boundary |> GT.physical_names
    Nface_to_tag = zeros(Int32,d_to_ndfaces[N+1])
    mesh = dirichlet_boundary |> GT.mesh
    #classify_mesh_faces!(Nface_to_tag,mesh,N,physical_names)
    Nface_to_tag[GT.faces(dirichlet_boundary)] .= 1
    ncells = length(cell_to_Dface)
    dim = 1
    ldof = :ldof
    sigma = GT.dual_basis_quantity(space,ldof)
    index1 = generate_index(domain(space))
    Dface = face_index(index1,num_dims(domain(space)))
    sigma_expr = term(sigma,index1) |> expression |> simplify
    index2 = GT.generate_index(dirichlet_boundary)
    dface = face_index(index2,num_dims(dirichlet_boundary))
    q_expr = term(q,index2) |> expression |> simplify
    sigma_texpr = topological_sort(sigma_expr,(Dface,ldof))
    # TODO Why empty?
    #q_texpr = topological_sort(q_expr,(dface,))
    expr = quote
        (args,storage1,storage2) -> begin
            $(unpack_index_storage(index1,:storage1))
            $(unpack_index_storage(index2,:storage2))
            d = args.d
            Dface_to_dfaces = args.d_to_Dface_to_dfaces[d+1]
            ctype_to_ldface_to_ldofs = args.d_to_ctype_to_ldface_to_dofs[d+1]
            $(sigma_texpr[1])
            for cell in 1:args.ncells
                ctype = args.cell_to_ctype[cell]
                $Dface = args.cell_to_Dface[cell]
                $(sigma_texpr[2])
                ldface_to_dface = Dface_to_dfaces[$Dface]
                ldface_to_ldofs = ctype_to_ldface_to_ldofs[ctype]
                nldfaces = length(ldface_to_dface)
                Dface = args.cell_to_Dface[cell]
                dofs = args.Dface_to_dofs[Dface]
                for ldface in 1:nldfaces
                    ldofs = ldface_to_ldofs[ldface]
                    $dface = ldface_to_dface[ldface]
                    Nface = $dface
                    tag = args.Nface_to_tag[Nface]
                    if tag == 0
                        continue
                    end
                    qfun = $(q_expr)
                    for ldof in ldofs
                        sfun = $(sigma_texpr[3])
                        mask = sfun(qfun)
                        if !(abs(mask) + 1 ≈ 1)
                            dof = dofs[ldof]
                            args.dof_to_tag[dof] = tag
                        end
                    end
                end
            end
        end
    end
    loop! = eval(expr)
    storage1 = GT.index_storage(index1)
    storage2 = GT.index_storage(index2)
    d = N
    args = (;d,ncells,dof_to_tag,Nface_to_tag,state...)
    invokelatest(loop!,args,storage1,storage2)
    free_and_dirichlet_dofs = GT.partition_from_mask(i->i==0,dof_to_tag)
    dof_permutation = GT.permutation(free_and_dirichlet_dofs)
    n_free_dofs = length(first(free_and_dirichlet_dofs))
    f = dof -> begin
        dof2 = dof_permutation[dof]
        T = typeof(dof2)
        if dof2 > n_free_dofs
            return T(n_free_dofs-dof2)
        end
        dof2
    end
    data = Dface_to_dofs.data
    data .= f.(data)
    ndiri = length(last(free_and_dirichlet_dofs))
    dirichlet_dof_location = ones(Int32,ndiri)
    free_dofs = Base.OneTo(n_free_dofs)
    dirichlet_dofs = Base.OneTo(ndiri)
    (;Dface_to_dofs, free_dofs, dirichlet_dofs, dirichlet_dof_location)
end

function generate_dof_ids_step_2(space,state,dirichlet_boundary_all::PiecewiseDomain)
    (;ndofs,Dface_to_dofs,d_to_Dface_to_dfaces,
     d_to_ctype_to_ldface_to_dofs,
     d_to_ndfaces,cell_to_ctype,cell_to_Dface) = state
    dof_to_location = zeros(Int32,ndofs)
    D = length(d_to_ndfaces)-1
    for d in D:-1:0
        for (location,dirichlet_boundary) in dirichlet_boundary_all.domains |> enumerate
            N = GT.num_dims(dirichlet_boundary)
            if d != N
                continue
            end
            #physical_names = dirichlet_boundary |> GT.physical_names
            Nface_to_tag = zeros(Int32,d_to_ndfaces[N+1])
            mesh = dirichlet_boundary |> GT.mesh
            #classify_mesh_faces!(Nface_to_tag,mesh,N,physical_names)
            Nface_to_tag[GT.faces(dirichlet_boundary)] .= 1
            ncells = length(cell_to_Dface)
            Dface_to_dfaces = d_to_Dface_to_dfaces[d+1]
            ctype_to_ldface_to_ldofs = d_to_ctype_to_ldface_to_dofs[d+1]
            for cell in 1:ncells
                ctype = cell_to_ctype[cell]
                Dface = cell_to_Dface[cell]
                ldof_to_dof = Dface_to_dofs[Dface]
                ldface_to_dface = Dface_to_dfaces[Dface]
                ldface_to_ldofs = ctype_to_ldface_to_ldofs[ctype]
                nldfaces = length(ldface_to_dface)
                dofs = Dface_to_dofs[Dface]
                for ldface in 1:nldfaces
                    ldofs = ldface_to_ldofs[ldface]
                    dface = ldface_to_dface[ldface]
                    Nface = dface
                    tag = Nface_to_tag[Nface]
                    if tag == 0
                        continue
                    end
                    dof_to_location[view(dofs,ldofs)] .= location
                end
            end
        end
    end
    free_and_dirichlet_dofs = GT.partition_from_mask(i->i==0,dof_to_location)
    dof_permutation = GT.permutation(free_and_dirichlet_dofs)
    n_free_dofs = length(first(free_and_dirichlet_dofs))
    f = dof -> begin
        dof2 = dof_permutation[dof]
        T = typeof(dof2)
        if dof2 > n_free_dofs
            return T(n_free_dofs-dof2)
        end
        dof2
    end
    data = Dface_to_dofs.data
    data .= f.(data)
    dirichlet_dof_location = dof_to_location[last(free_and_dirichlet_dofs)]
    free_dofs = Base.OneTo(n_free_dofs)
    ndiri = length(last(free_and_dirichlet_dofs))
    dirichlet_dofs = Base.OneTo(ndiri)
    (;Dface_to_dofs, free_dofs, dirichlet_dofs, dirichlet_dof_location)
end

function generate_dof_ids_step_2(space,state,q_all::PiecewiseField)
    (;ndofs,Dface_to_dofs,d_to_Dface_to_dfaces,
     d_to_ctype_to_ldface_to_dofs,
     d_to_ndfaces,cell_to_ctype,cell_to_Dface) = state
    dof_to_tag = zeros(Int32,ndofs)
    D = length(d_to_ndfaces)-1
    for d in D:-1:0
        for (location,q) in enumerate(q_all.fields)
            dirichlet_boundary = domain(q)
            N = GT.num_dims(dirichlet_boundary)
            if d != N
                continue
            end
            #physical_names = dirichlet_boundary |> GT.physical_names
            Nface_to_tag = zeros(Int32,d_to_ndfaces[N+1])
            mesh = dirichlet_boundary |> GT.mesh
            #classify_mesh_faces!(Nface_to_tag,mesh,N,physical_names)
            Nface_to_tag[GT.faces(dirichlet_boundary)] .= 1
            ncells = length(cell_to_Dface)
            dim = 1
            ldof = :ldof
            sigma = GT.dual_basis_quantity(space,ldof)
            index1 = generate_index(domain(space))
            Dface = face_index(index1,num_dims(domain(space)))
            sigma_expr = term(sigma,index1) |> expression |> simplify
            index2 = GT.generate_index(dirichlet_boundary)
            dface = face_index(index2,num_dims(dirichlet_boundary))
            q_expr = term(q,index2) |> expression |> simplify
            sigma_texpr = topological_sort(sigma_expr,(Dface,ldof))
            # TODO Why empty?
            #q_texpr = topological_sort(q_expr,(dface,))
            expr = quote
                (args,storage1,storage2) -> begin
                    $(unpack_index_storage(index1,:storage1))
                    $(unpack_index_storage(index2,:storage2))
                    d = args.d
                    Dface_to_dfaces = args.d_to_Dface_to_dfaces[d+1]
                    ctype_to_ldface_to_ldofs = args.d_to_ctype_to_ldface_to_dofs[d+1]
                    $(sigma_texpr[1])
                    for cell in 1:args.ncells
                        ctype = args.cell_to_ctype[cell]
                        $Dface = args.cell_to_Dface[cell]
                        $(sigma_texpr[2])
                        ldface_to_dface = Dface_to_dfaces[$Dface]
                        ldface_to_ldofs = ctype_to_ldface_to_ldofs[ctype]
                        nldfaces = length(ldface_to_dface)
                        Dface = args.cell_to_Dface[cell]
                        dofs = args.Dface_to_dofs[Dface]
                        for ldface in 1:nldfaces
                            ldofs = ldface_to_ldofs[ldface]
                            $dface = ldface_to_dface[ldface]
                            Nface = $dface
                            tag = args.Nface_to_tag[Nface]
                            if tag == 0
                                continue
                            end
                            qfun = $(q_expr)
                            for ldof in ldofs
                                sfun = $(sigma_texpr[3])
                                mask = sfun(qfun)
                                if !(abs(mask) + 1 ≈ 1)
                                    dof = dofs[ldof]
                                    args.dof_to_tag[dof] = tag
                                end
                            end
                        end
                    end
                end
            end
            loop! = eval(expr)
            storage1 = GT.index_storage(index1)
            storage2 = GT.index_storage(index2)
            d = N
            args = (;d,ncells,dof_to_tag,Nface_to_tag,state...)
            invokelatest(loop!,args,storage1,storage2)
        end
    end
    free_and_dirichlet_dofs = GT.partition_from_mask(i->i==0,dof_to_tag)
    dof_permutation = GT.permutation(free_and_dirichlet_dofs)
    n_free_dofs = length(first(free_and_dirichlet_dofs))
    f = dof -> begin
        dof2 = dof_permutation[dof]
        T = typeof(dof2)
        if dof2 > n_free_dofs
            return T(n_free_dofs-dof2)
        end
        dof2
    end
    data = Dface_to_dofs.data
    data .= f.(data)
    ndiri = length(last(free_and_dirichlet_dofs))
    dirichlet_dof_location = ones(Int32,ndiri)
    free_dofs = Base.OneTo(n_free_dofs)
    dirichlet_dofs = Base.OneTo(ndiri)
    (;Dface_to_dofs, free_dofs, dirichlet_dofs, dirichlet_dof_location)
end

struct LastDof end

function last_dof()
    LastDof()
end

function generate_dof_ids_step_2(space,state,dirichlet_boundary_all::LastDof)
    (;ndofs,Dface_to_dofs) = state
    dof_to_tag = zeros(Int32,ndofs)
    dof_to_tag[end] = 1
    free_and_dirichlet_dofs = GT.partition_from_mask(i->i==0,dof_to_tag)
    dof_permutation = GT.permutation(free_and_dirichlet_dofs)
    n_free_dofs = length(first(free_and_dirichlet_dofs))
    f = dof -> begin
        dof2 = dof_permutation[dof]
        T = typeof(dof2)
        if dof2 > n_free_dofs
            return T(n_free_dofs-dof2)
        end
        dof2
    end
    data = Dface_to_dofs.data
    data .= f.(data)
    ndiri = length(last(free_and_dirichlet_dofs))
    dirichlet_dof_location = ones(Int32,ndiri)
    free_dofs = Base.OneTo(n_free_dofs)
    dirichlet_dofs = Base.OneTo(ndiri)
    (;Dface_to_dofs, free_dofs, dirichlet_dofs, dirichlet_dof_location)
end

function reference_face_own_dofs(space::AbstractSpace,d)
    ctype_to_reference_fe = reference_spaces(space)
    ctype_to_ldface_to_own_ldofs = map(fe->GT.face_own_dofs(fe,d),ctype_to_reference_fe)
    if GT.conformity(space) === :default
        ctype_to_ldface_to_own_ldofs
    elseif GT.conformity(space) === :L2
        ctype_to_num_dofs = map(GT.num_dofs,ctype_to_reference_fe)
        domain = space |> GT.domain
        D = GT.num_dims(domain)
        map(ctype_to_num_dofs,ctype_to_ldface_to_own_ldofs) do ndofs,ldface_to_own_ldofs
            map(ldface_to_own_ldofs) do own_ldofs
                dofs = if d == D
                    collect(1:ndofs)
                else
                    Int[]
                end
                convert(typeof(own_ldofs),dofs)
            end
        end
    else
        error("This line cannot be reached")
    end
end

function reference_face_own_dof_permutations(space::AbstractSpace,d)
    ctype_to_reference_fe = reference_spaces(space)
    ctype_to_ldface_to_pindex_to_perm = map(fe->GT.face_own_dof_permutations(fe,d),ctype_to_reference_fe)
    if GT.conformity(space) === :default
        ctype_to_ldface_to_pindex_to_perm
    elseif GT.conformity(space) === :L2
        ctype_to_num_dofs = map(GT.num_dofs,ctype_to_reference_fe)
        domain = space |> GT.domain
        D = GT.num_dims(domain)
        map(ctype_to_num_dofs,ctype_to_ldface_to_pindex_to_perm) do ndofs,ldface_to_pindex_to_perm
            map(ldface_to_pindex_to_perm) do pindex_to_perm
                map(pindex_to_perm) do perm
                    dofs = if d == D
                        collect(1:ndofs)
                    else
                        Int[]
                    end
                    convert(typeof(perm),dofs)
                end
            end
        end
    else
        error("This line cannot be reached")
    end
end

partition_from_mask(a) = partition_from_mask(identity,a)

function partition_from_mask(f,node_to_mask)
    T = Vector{Int32}
    free_nodes = convert(T,findall(f,node_to_mask))
    dirichlet_nodes = convert(T,findall(i->!f(i),node_to_mask))
    nfree = length(free_nodes)
    ndiri = length(dirichlet_nodes)
    permutation = T(undef,nfree+ndiri)
    permutation[free_nodes] = 1:nfree
    permutation[dirichlet_nodes] = (1:ndiri) .+ nfree
    TwoWayPartition(free_nodes,dirichlet_nodes,permutation)
end

struct TwoWayPartition{A} <: AbstractVector{A}
    first::A
    last::A
    permutation::A
end

permutation(a::TwoWayPartition) = a.permutation
Base.size(a::TwoWayPartition) = (2,)
Base.IndexStyle(::Type{<:TwoWayPartition}) = IndexLinear()
function Base.getindex(a::TwoWayPartition,i::Int)
    @boundscheck @assert i in (1,2)
    if i == 1
        a.first
    else
        a.last
    end
end

"""
"""
function shape_functions(fe::AbstractSpace)
    primal = primal_basis(fe)
    dual = dual_basis(fe)
    primal_t = permutedims(primal)
    A = value.(dual,primal_t)
    B = A\I
    n = length(primal)
    map(1:n) do i
        Bi = view(B,:,i)
        x->begin
            primal_t_x = map(f->f(x),primal_t)
            (primal_t_x*Bi)[1,1]#TODO
        end
    end
end

"""
"""
function tabulator(fe::AbstractSpace)
    primal = primal_basis(fe)
    dual = dual_basis(fe)
    primal_t = permutedims(primal)
    A = value.(dual,primal_t)
    B = A\I
    (f,x) -> begin
        C = broadcast(f,primal_t,x)
        C*B
    end
end

options(fe::AbstractFaceSpace) = options(domain(fe))
num_dims(fe::AbstractFaceSpace) = num_dims(domain(fe))

"""
"""
function lagrange_space end

const LagrangeFaceDomain = Union{UnitNCube,UnitSimplex}

function lagrange_space(domain::LagrangeFaceDomain, order;
        space_type = default_space_type(domain),
        lib_to_user_nodes = :default,
        major = Val(:component),
        tensor_size = Val(:scalar),
        dirichlet_boundary = nothing,
    )


    D = num_dims(domain)
    order_per_dir = ntuple(d->order,Val(D))
    lagrange_face_space(;
               domain,
               order_per_dir,
               space_type,
               lib_to_user_nodes,
               major,
               tensor_size,
               dirichlet_boundary,
              )
end

function default_space_type(geom::UnitNCube)
    :Q
end

function default_space_type(geom::UnitSimplex)
    :P
end

function lagrange_face_space(;
        domain,
        order_per_dir,
        space_type,
        lib_to_user_nodes,
        major,
        tensor_size,
        dirichlet_boundary,
    )
    contents = (;
        domain,
        order_per_dir,
        space_type,
        lib_to_user_nodes,
        major,
        tensor_size,
        dirichlet_boundary,
       )
    LagrangeFaceSpace(contents)
end

struct LagrangeFaceSpace{A} <: AbstractFaceSpace
    contents::A
end

domain(a::LagrangeFaceSpace) = a.contents.domain
order_per_dir(a::LagrangeFaceSpace) = a.contents.order_per_dir
order(fe::LagrangeFaceSpace) = maximum(order_per_dir(fe);init=0)
space_type(fe::LagrangeFaceSpace) = val_parameter(fe.contents.space_type)
major(fe::LagrangeFaceSpace) = val_parameter(fe.contents.major)
tensor_size(fe::LagrangeFaceSpace) = val_parameter(fe.contents.tensor_size)
dirichlet_boundary(fe::LagrangeFaceSpace) = fe.contents.dirichlet_boundary

function lib_to_user_nodes(fe::LagrangeFaceSpace)
    if val_parameter(fe.contents.lib_to_user_nodes) === :default
        nnodes = num_nodes(fe)
        Ti = int_type(options(fe))
        collect(Ti.(1:nnodes))
    else
        fe.contents.lib_to_user_nodes
    end
end

function reference_spaces(fe::LagrangeFaceSpace)
    (fe,)
end

function face_reference_id(fe::LagrangeFaceSpace)
    [1]
end

function conformity(fe::LagrangeFaceSpace)
    :default
end

function monomial_exponents(a::LagrangeFaceSpace)
    range_per_dir = map(k->0:k,order_per_dir(a))
    exponents_list = map(CartesianIndices(range_per_dir)) do ci
        exponent = Tuple(ci)
        Ti = int_type(options(a))
        D = length(exponent)
        SVector{D,Ti}(exponent)
    end[:]
    if space_type(a) === :Q
        exponents_list = exponents_list
    elseif space_type(a) === :P
        exponents_list = filter(exponents->sum(exponents)<=order(a),exponents_list)
    else
        error("space_type == $(space_type(a)) case not implemented (yet)")
    end
end

num_nodes(fe::LagrangeFaceSpace) = length(monomial_exponents(fe))

function node_coordinates(a::LagrangeFaceSpace)
    if order(a) == 0 && num_dims(a) != 0
        a_linear = lagrange_space(domain(a),1)
        x  = node_coordinates(a_linear)
        return [ sum(x)/length(x) ]
    end
    @assert a |> domain |> is_unitary
    exponents_list = monomial_exponents(a)
    lib_node_to_coords = map(exponents_list) do exponent
        t = map(exponent,order_per_dir(a)) do e,order
            Tv = real_type(options(a))
            if order != 0
               Tv(e/order)
            else
               Tv(e)
            end
        end
        SVector{length(order_per_dir(a)),real_type(options(a))}(t)
    end
    if lib_to_user_nodes(a) === :default
        return lib_node_to_coords
    end
    user_node_to_coords = similar(lib_node_to_coords)
    user_node_to_coords[lib_to_user_nodes(a)] = lib_node_to_coords
    user_node_to_coords
end

function tensor_basis(fe::LagrangeFaceSpace)
    s = tensor_size(fe)
    Tv = fe |> options |> real_type
    if tensor_size(fe) === :scalar
        return Tv(1)
    else
        cis = CartesianIndices(s)
        l = prod(s)
        init_tensor = SArray{Tuple{s...},Tv}
        cis_flat = cis[:]
        return map(cis_flat) do ci
            init_tensor(ntuple(j->cis[j]==ci ? 1 : 0 ,Val(l)))
        end
    end
end

function primal_basis(fe::LagrangeFaceSpace)
    scalar_basis = map(e->(x-> prod(x.^e)),monomial_exponents(fe))
    if tensor_size(fe) === :scalar
        return scalar_basis
    else
        primal_nested = map(scalar_basis) do monomial
            map(tensor_basis(fe))  do e
                x -> monomial(x)*e
            end
        end
        return reduce(vcat,primal_nested)
    end
end

function dual_basis(fe::LagrangeFaceSpace)
    node_coordinates_reffe = node_coordinates(fe)
    scalar_basis = map(x->(f->f(x)),node_coordinates_reffe)
    ts = tensor_size(fe) 
    if ts === :scalar
        return scalar_basis
    else
        if major(fe) === :component
            dual_nested = map(node_coordinates_reffe) do x
                map(tensor_basis(fe)) do e
                    f->contraction(e,f(x))
                end
            end
        elseif major(fe) === :node
            dual_nested = map(tensor_basis(fe)) do e
                map(node_coordinates_reffe) do x
                    f->contraction(e,f(x))
                end
            end
        else
            error("Not Implemented")
        end
        return reduce(vcat,dual_nested)
    end
end

contraction(a,b) = sum(map(*,a,b))

value(f,x) = f(x)

function node_dofs(fe::LagrangeFaceSpace)
    Tv = int_type(options(fe))
    nnodes = num_nodes(fe)
    nodes =  1:nnodes
    s = tensor_size(fe)
    if s === :scalar
        return nodes
    else
        ndofs_per_node = prod(tensor_size(fe))
        init_tensor = SArray{Tuple{tensor_size(fe)...},Tv}
        node_to_dofs = map(nodes) do node
            t = ntuple(Val(ndofs_per_node)) do li
                if major(fe) === :component
                    dof = (node-1)*ndofs_per_node + li
                elseif major(fe) === :node
                    dof = node + (li-1)*nnodes
                else
                    error("Not Implemented")
                end
            end
            init_tensor(t)
        end
        return node_to_dofs
    end
end

function dof_node(fe::LagrangeFaceSpace)
    Tv = int_type(options(fe))
    ndofs = num_dofs(fe)
    if tensor_size(fe) === :scalar
        return collect(Tv,1:ndofs)
    else
        dof_to_node = zeros(Tv,ndofs)
        node_to_dofs = node_dofs(fe)
        for (node,dofs) in enumerate(node_to_dofs)
            for dof in dofs
                dof_to_node[dof] = node
            end
        end
        return dof_to_node
    end
end

function conforming(fe::LagrangeFaceSpace)
    if order(fe) == 0
        return false
    end
    if num_dims(fe) in (0,1)
        return true
    end
    nonconforming = (is_n_cube(domain(fe)) && space_type(fe) === :P)
    ! nonconforming
end

function interior_nodes(fe::LagrangeFaceSpace)
    if !conforming(fe)
        return collect(Int,1:num_nodes(fe))
    end
    nnodes = num_nodes(fe)
    D = num_dims(fe)
    if D == 0
        return collect(1:nnodes)
    else
        mesh = complexify(fe)
        node_is_touched = fill(true,nnodes)
        for d in 0:(D-1)
            face_to_nodes = face_nodes(mesh,d)
            for nodes in face_to_nodes
                node_is_touched[nodes] .= false
            end
        end
        return findall(node_is_touched)
    end
end

function num_interior_nodes(fe::LagrangeFaceSpace)
    length(interior_nodes(fe))
end

function face_nodes(fe::LagrangeFaceSpace,d)
    if conforming(fe)
        face_nodes_from_mesh_face(fe,d)
    else
        D = num_dims(fe)
        if d == D
            [collect(Int32,1:GT.num_nodes(fe))]
        else
            [Int32[] for _ in 1:num_faces(mesh(domain(fe)),d)]
        end
    end
end

function face_nodes_from_mesh_face(fe,d)
    D = num_dims(fe)
    if d == D
        [collect(Int32,1:GT.num_nodes(fe))]
    else
        boundary = GT.complexify(fe)
        GT.face_nodes(boundary,d)
    end
end

function face_interior_nodes(fe::LagrangeFaceSpace,d)
    if conforming(fe)
        face_interior_nodes_from_mesh_face(fe,d)
    else
        face_nodes(fe,d)
    end
end

function face_interior_nodes_from_mesh_face(fe,d)
    D = num_dims(fe)
    if  d == D
        [GT.interior_nodes(fe)]
    else
        boundary = GT.complexify(fe)
        dface_to_lnode_to_node = GT.face_nodes(boundary,d)
        dface_to_ftype = GT.face_reference_id(boundary,d)
        ftype_to_refdface = GT.reference_spaces(boundary,d)
        ftype_to_lnodes = map(GT.interior_nodes,ftype_to_refdface)
        map(dface_to_ftype,dface_to_lnode_to_node) do ftype,lnode_to_node
            lnodes = ftype_to_lnodes[ftype]
            lnode_to_node[lnodes]
        end
    end
end

function face_interior_node_permutations(fe::LagrangeFaceSpace,d)
    if conforming(fe)
        face_interior_node_permutations_from_mesh_face(fe,d)
    else
        D = num_dims(fe)
        if  d == D
            [[ collect(1:num_interior_nodes(fe)) ]]
        else
            [[ Int32[] ] for _ in 1:num_faces(mesh(domain(fe)),d)]
        end
    end
end

function face_interior_node_permutations_from_mesh_face(fe,d)
    D = num_dims(fe)
    if  d == D
        [[ collect(1:num_interior_nodes(fe)) ]]
    else
        boundary = GT.complexify(fe)
        dface_to_ftype = GT.face_reference_id(boundary,d)
        ftype_to_refdface = GT.reference_spaces(boundary,d)
        ftype_to_perms = map(GT.interior_node_permutations,ftype_to_refdface)
        map(dface_to_ftype) do ftype
            perms = ftype_to_perms[ftype]
        end
    end
end

"""
"""
function interior_node_permutations(fe::LagrangeFaceSpace)
    interior_ho_nodes = interior_nodes(fe)
    node_permutations_from_mesh_face(fe,interior_ho_nodes)
end

"""
"""
function node_permutations(fe::LagrangeFaceSpace)
    interior_ho_nodes = 1:num_nodes(fe)
    node_permutations_from_mesh_face(fe,interior_ho_nodes)
end

function node_permutations_from_mesh_face(refface,interior_ho_nodes)
    ho_nodes_coordinates = node_coordinates(refface)
    geo = domain(refface)
    vertex_perms = vertex_permutations(geo)
    if length(interior_ho_nodes) == 0
        return map(i->Int[],vertex_perms)
    end
    if length(vertex_perms) == 1
        return map(i->collect(1:length(interior_ho_nodes)),vertex_perms)
    end
    if order(refface) == 0 # TODO ugly. It assumes the hack above for node coordinates of faces of order 0
        return map(i->collect(1:length(interior_ho_nodes)),vertex_perms)
    end
    geo_mesh = mesh(geo)
    vertex_to_geo_nodes = face_nodes(geo_mesh,0)
    vertex_to_geo_node = map(first,vertex_to_geo_nodes)
    ref_face = lagrange_space(geo,1)
    fun_mesh = complexify(ref_face)
    geo_node_coords = node_coordinates(geo_mesh)
    fun_node_coords = node_coordinates(fun_mesh)
    vertex_coords = geo_node_coords[vertex_to_geo_node]
    q = ho_nodes_coordinates[interior_ho_nodes]
    Tx = eltype(vertex_coords)
    A = zeros(Float64,length(q),length(fun_node_coords))
    A = tabulator(ref_face)(value,q)
    perm_vertex_coords = similar(vertex_coords)
    node_perms = similar(vertex_perms)
    for (iperm,permutation) in enumerate(vertex_perms)
        for (j,cj) in enumerate(permutation)
          perm_vertex_coords[j] = vertex_coords[cj]
        end
        node_to_pnode = fill(INVALID_ID,length(interior_ho_nodes))
        for iq in 1:size(A,1)
            y = zero(Tx)
            for fun_node in 1:size(A,2)
                vertex = fun_node # TODO we are assuming that the vertices and nodes match
                g = A[iq,fun_node]
                x = perm_vertex_coords[vertex]
                y += g*x
            end
            pnode = findfirst(i->(norm(i-y)+1)≈1,q)
            if !isnothing(pnode)
               node_to_pnode[iq] = pnode
            end
        end
        node_perms[iperm] = node_to_pnode
    end
    node_perms
end

function face_dofs(a::LagrangeFaceSpace,d)
    face_to_nodes = face_nodes(a,d)
    face_dofs_from_nodes(a,face_to_nodes)
end

function face_own_dofs(a::LagrangeFaceSpace,d)
    face_to_nodes = face_interior_nodes(a,d)
    face_dofs_from_nodes(a,face_to_nodes)
end

function face_dofs_from_nodes(fe::LagrangeFaceSpace,face_to_nodes)
    if tensor_size(fe) === :scalar
        return face_to_nodes
    else
        node_to_dofs = node_dofs(fe)
        lis = LinearIndices(tensor_size(fe))
        face_to_dofs = map(face_to_nodes) do nodes
            if major(fe) === :component
                nested = map(nodes) do node
                    dofs = node_to_dofs[node]
                    collect(dofs[:])
                end
            elseif major(fe) === :node
                nested = map(lis[:]) do li
                    map(nodes) do node
                        dofs = node_to_dofs[node]
                        dofs[li]
                    end
                end
            else
                error("Not Implemented")
            end
            reduce(vcat,nested;init=Int32[])
        end
    end
end

function face_own_dof_permutations(fe::LagrangeFaceSpace,d)
    face_to_pindex_to_inodes = face_interior_node_permutations(fe,d)
    if tensor_size(fe) === :scalar
        return face_to_pindex_to_inodes
    else
        Tv = fe |> options |> int_type
        node_to_dofs = node_dofs(fe)
        face_to_inodes = face_interior_nodes(fe,d)
        face_to_idofs = face_own_dofs(fe,d)
        nfaces = length(face_to_inodes)
        ndofs = num_dofs(fe)
        dof_to_idof = zeros(Tv,ndofs)
        lis = LinearIndices(tensor_size(fe))
        map(1:nfaces) do face
            pindex_to_inodes = face_to_pindex_to_inodes[face]
            inode_to_node = face_to_inodes[face]
            idof_to_dof = face_to_idofs[face]
            fill!(dof_to_idof,Tv(0))
            dof_to_idof[idof_to_dof] = 1:length(idof_to_dof)
            map(pindex_to_inodes) do inodes
                if major(fe) === :component
                    nested = map(inodes) do inode
                        node = inode_to_node[inode]
                        dofs = node_to_dofs[node]
                        map(dofs) do dof
                            idof = dof_to_idof[dof]
                            @assert idof != 0
                            idof
                        end
                    end
                elseif major(fe) === :node
                    nested = map(lis[:]) do li
                        map(inodes) do inode
                            node = inode_to_node[inode]
                            dofs = node_to_dofs[node]
                            dof = dofs[li]
                            idof = dof_to_idof[dof]
                            @assert idof != 0
                            idof
                        end
                    end
                else
                    error("Not Implemented")
                end
                reduce(vcat,nested;init=Tv[])
            end
        end
    end
end

function num_dofs(a::LagrangeFaceSpace)
    nnodes = num_nodes(a)
    if tensor_size(a) === :scalar
        nnodes
    else
        ndofs_per_node = prod(tensor_size(a))
        nnodes*ndofs_per_node
    end
end

function node_quadrature(fe::LagrangeFaceSpace)
    coordinates = node_coordinates(fe)
    Tv = real_type(options(fe))
    nnodes = length(coordinates)
    weights = fill(Tv(1/nnodes),nnodes)
    domain = GT.domain(fe)
    face_quadrature(;domain,coordinates,weights)
end

function complexify(refface::LagrangeFaceSpace)
    geom = domain(refface)
    mesh = GT.mesh(geom)
    D = num_dims(geom)
    order_inter = order(refface)
    round_coord(coord) = map(xi->round(Int,order_inter*xi),coord)
    node_coordinates = GT.node_coordinates(refface)
    node_coordinates_mapped = map(round_coord,node_coordinates)
    Ti = int_type(options(geom))
    dims = ntuple(d->d-1,Val(D+1))
    face_reference_id = GT.face_reference_id(mesh)
    reference_spaces = map(GT.reference_domains(mesh)) do rid_to_domain
        map(domain->lagrange_space(domain,order_inter),rid_to_domain)
    end
    face_nodes_tuple = map(dims) do d
        domain = GT.domain(mesh,d)
        rid_to_refqua = map(reference_domains(domain)) do refdom
            reffe = lagrange_space(refdom,order_inter)
            node_quadrature(reffe)
        end
        face_to_rid = face_reference_id[d+1]
        quadrature = GT.mesh_quadrature(;
            domain,
            reference_quadratures=rid_to_refqua,
            face_reference_id = face_to_rid
           )
        face_point_x = coordinate_accessor(quadrature)
        face_npoints = num_points_accessor(quadrature)
        nfaces = num_faces(mesh,d)
        dface_to_nodes = Vector{Vector{Ti}}(undef,nfaces)
        for face in 1:nfaces
            npoints = face_npoints(face)
            point_x = face_point_x(face)
            x_mapped = map(1:npoints) do point
                x = point_x(point)
                round_coord(x)
            end
            my_nodes = indexin(x_mapped,node_coordinates_mapped)
            dface_to_nodes[face] = my_nodes
        end
        dface_to_nodes
    end
    face_nodes = collect(face_nodes_tuple)
    outward_normals = GT.outward_normals(mesh)
    GT.mesh(;
        node_coordinates,
        face_nodes,
        face_reference_id,
        reference_spaces,
        is_cell_complex = Val(true),
        outward_normals
       )
end

function simplexify(ref_face::LagrangeFaceSpace)
    mesh_geom  = simplexify(domain(ref_face))
    D = num_dims(mesh_geom)
    node_coordinates_geom = node_coordinates(mesh_geom)
    ref_faces_geom = reference_spaces(mesh_geom,D)
    face_nodes_geom = face_nodes(mesh_geom,D)
    face_ref_id_geom = face_reference_id(mesh_geom,D)
    nfaces = length(face_ref_id_geom)
    # We need the same order in all directions
    # for this to make sense
    my_order = order(ref_face)
    ref_faces_inter = map(r_geom->lagrange_space(domain(r_geom),my_order),ref_faces_geom)
    s_ref = map(ref_faces_geom,ref_faces_inter) do r_geom,r
        m = num_nodes(r)
        n = num_nodes(r_geom)
        x = node_coordinates(r)
        tabulator(r_geom)(value,x)
    end
    node_coordinates_inter = node_coordinates(ref_face)
    node_coordinates_aux = map(xi->map(xii->round(Int,my_order*xii),xi),node_coordinates_inter)
    face_nodes_inter = Vector{Vector{Int}}(undef,nfaces)
    for face in 1:nfaces
        ref_id_geom = face_ref_id_geom[face]
        s = s_ref[ref_id_geom]
        nodes_geom = face_nodes_geom[face]
        nnodes, nnodes_geom = size(s)
        x_mapped = map(1:nnodes) do i
            x = zero(eltype(node_coordinates_inter))
            for k in 1:nnodes_geom
                x += node_coordinates_geom[nodes_geom[k]]*s[i,k]
            end
            map(xi->round(Int,my_order*xi),x)
        end
        my_nodes = indexin(x_mapped,node_coordinates_aux)
        face_nodes_inter[face] = my_nodes
    end
    ref_inter = 
    chain = GT.chain(;
        node_coordinates=node_coordinates_inter,
        face_nodes = face_nodes_inter,
        face_reference_id = face_ref_id_geom,
        reference_spaces = ref_faces_inter,
    )
    mesh = GT.mesh(chain)
    mesh_complex = complexify(mesh)
    pg = physical_faces(mesh_complex)
    pg .= physical_faces(mesh_geom)
    mesh_complex
end

function lagrange_space(domain::AbstractDomain,order;
    conformity = :default,
    dirichlet_boundary=nothing,
    space_type = Val(:default),
    major = Val(:component),
    tensor_size = Val(:scalar),
    workspace = nothing,
    setup = Val(true),
    )

    @assert conformity in (:default,:L2)

    space = lagrange_mesh_space(;
                        domain,
                        order,
                        conformity,
                        dirichlet_boundary,
                        space_type,
                        major,
                        tensor_size,
                        workspace,
                       )
    if val_parameter(setup)
        setup_space(space)
    else
        space
    end
end

function lagrange_mesh_space(;
        domain,
        order,
        conformity,
        dirichlet_boundary,
        space_type,
        major,
        tensor_size,
        workspace,
    )
    contents = (;
        domain,
        order,
        conformity,
        dirichlet_boundary,
        space_type,
        major,
        tensor_size,
        workspace,
       )
    LagrangeMeshSpace(mesh(domain),contents)
end

struct LagrangeMeshSpace{A,B} <: AbstractSpace{A}
    mesh::A
    contents::B
end

function PartitionedArrays.partition(pspace::LagrangeMeshSpace)
    if GT.workspace(pspace) !== nothing
        return GT.workspace(pspace).space_partition
    end
    p_domain = partition(GT.domain(pspace))
    pdirichlet_boundary = GT.dirichlet_boundary(pspace)
    if pdirichlet_boundary isa AbstractDomain
        p_dirichlet_boundary = partition(GT.dirichlet_boundary(pspace))
        map(p_domain,p_dirichlet_boundary) do domain, dirichlet_boundary
            lagrange_space(
                           domain,
                           order(pspace);
                           conformity = conformity(pspace),
                           dirichlet_boundary,
                           space_type = space_type(pspace),
                           major = major(pspace),
                           tensor_size = tensor_size(pspace),
                           setup = Val(false),
                          )
        end

    else
        map(p_domain) do domain
            lagrange_space(
                           domain,
                           order(pspace);
                           conformity = conformity(pspace),
                           dirichlet_boundary = pdirichlet_boundary,
                           space_type = space_type(pspace),
                           major = major(pspace),
                           tensor_size = tensor_size(pspace),
                           setup = Val(false),
                          )
        end
    end
end

conformity(space::LagrangeMeshSpace) = space.contents.conformity
dirichlet_boundary(space::LagrangeMeshSpace) = space.contents.dirichlet_boundary
domain(space::LagrangeMeshSpace) = space.contents.domain
order(space::LagrangeMeshSpace) = space.contents.order
space_type(space::LagrangeMeshSpace) = val_parameter(space.contents.space_type)
major(space::LagrangeMeshSpace) = val_parameter(space.contents.major)
tensor_size(space::LagrangeMeshSpace) = val_parameter(space.contents.tensor_size)
workspace(space::LagrangeMeshSpace) = space.contents.workspace

function replace_workspace(space::LagrangeMeshSpace,workspace)
    contents = (;
        domain = space.contents.domain,
        order = space.contents.order,
        conformity = space.contents.conformity,
        space_type = space.contents.space_type,
        major = space.contents.major,
        tensor_size = space.contents.tensor_size,
        workspace,
       )
    LagrangeMeshSpace(space.mesh,contents)
end

function face_reference_id(space::LagrangeMeshSpace)
    domain = space |> GT.domain
    mesh = domain |> GT.mesh
    cell_to_Dface = domain |> GT.faces
    D = domain |> GT.num_dims
    Dface_to_ctype = GT.face_reference_id(mesh,D)
    cell_to_ctype = Dface_to_ctype[cell_to_Dface]
    cell_to_ctype
end

function reference_spaces(space::LagrangeMeshSpace)
    domain = space |> GT.domain
    mesh = domain |> GT.mesh
    D = domain |> GT.num_dims
    space_type = GT.space_type(space) # TODO Ugly
    ctype_to_refface = GT.reference_spaces(mesh,D)
    ctype_to_geometry = map(GT.domain,ctype_to_refface)
    ctype_to_reffe = map(ctype_to_geometry) do geometry
        space2 = space_type === :default ? default_space_type(geometry) : space_type
        lagrange_space(geometry,order(space);
           space_type=Val(space2),
           major=Val(major(space)),
           tensor_size = Val(tensor_size(space)))
    end
    ctype_to_reffe
end

function face_nodes(a::LagrangeMeshSpace)
    V = lagrange_space(
                       domain(a),
                       order(a);
                       conformity = conformity(a),
                       space_type = space_type(a))

    face_dofs(V)
end

function node_coordinates(a::LagrangeMeshSpace)
    V = lagrange_space(
                       GT.domain(a),
                       GT.order(a);
                       conformity = GT.conformity(a),
                       space_type = GT.space_type(a))
    vrid_to_reffe = reference_spaces(V)
    mface_to_vrid = face_reference_id(V)
    domain = GT.domain(a)
    mesh = GT.mesh(domain)
    d = num_dims(domain)
    mrid_to_refface = reference_spaces(mesh,d)
    mface_to_mrid = face_reference_id(mesh,d)
    vface_to_mface = faces(domain)
    nvfaces = length(vface_to_mface)
    mface_to_nodes = face_dofs(V)
    nnodes = length(free_dofs(V))
    mnode_to_x = node_coordinates(mesh)
    T = eltype(mnode_to_x)
    z = zero(T)
    node_to_x = zeros(T,nnodes)
    mface_to_mnodes = face_nodes(mesh,d)
    vrid_mrid_tabulator = map(vrid_to_reffe) do reffe
        map(mrid_to_refface) do refface
            tabulator(refface)(value,node_coordinates(reffe))
        end
    end
    for vface in 1:nvfaces
        mface = vface_to_mface[vface]
        vrid = mface_to_vrid[mface]
        mrid = mface_to_mrid[mface]
        tab = vrid_mrid_tabulator[vrid][mrid]
        mnodes = mface_to_mnodes[mface]
        nlnodes,nlmnodes = size(tab)
        nodes = mface_to_nodes[mface]
        for lnode in 1:nlnodes
            x = z
            for lmnode in 1:nlmnodes
                x += tab[lnode,lmnode]*mnode_to_x[mnodes[lmnode]]
            end
            node = nodes[lnode]
            node_to_x[node] = x
        end
    end
    node_to_x
end

#TODO these would provably need loop over cells
#function node_dofs(space::LagrangeMeshSpace)
#end
#

function free_dof_node(space::LagrangeMeshSpace)
    free_and_dirichlet_dof_node(space)[1]
end

function dirichlet_dof_node(space::LagrangeMeshSpace)
    free_and_dirichlet_dof_node(space)[2]
end

function free_and_dirichlet_dof_node(space::LagrangeMeshSpace)
    T = Float64
    D = num_ambient_dims(GT.domain(space))
    V = lagrange_space(
                       GT.domain(space),
                       GT.order(space);
                       conformity = GT.conformity(space),
                       space_type = GT.space_type(space))
    face_to_nodes = GT.face_dofs(V)
    face_to_dofs = GT.face_dofs(space)
    nfree = length(free_dofs(space))
    ndiri = length(dirichlet_dofs(space))
    free_dof_to_node = zeros(Int32,nfree)
    diri_dof_to_node = zeros(Int32,ndiri)
    rid_to_lnode_to_ldofs = map(node_dofs,reference_spaces(space))
    face_to_rid = face_reference_id(space)
    nfaces = length(face_to_rid)
    for face in 1:nfaces
        nodes = face_to_nodes[face]
        dofs = face_to_dofs[face]
        rid = face_to_rid[face]
        lnode_to_ldofs = rid_to_lnode_to_ldofs[rid]
        for (lnode,node) in enumerate(nodes)
            ldofs = lnode_to_ldofs[lnode]
            for ldof in ldofs
                dof = dofs[ldof]
                if dof < 0
                    diri_dof_to_node[-dof] = node
                else
                    free_dof_to_node[dof] = node
                end
            end
        end
    end
    free_dof_to_node, diri_dof_to_node
end

"""
"""
function raviart_thomas_space end

function raviart_thomas_space(domain::AbstractFaceDomain,order)
    workspace = rt_setup((;domain,order))
    RaviartThomasFaceSpace(domain,order,workspace)
end

struct RaviartThomasFaceSpace{A,B,C} <: AbstractFaceSpace
    domain::A
    order::B
    workspace::C
end

function rt_setup(fe)
    face_workspace,offset = rt_setup_dual_basis_boundary(fe)
    if is_n_cube(fe.domain)
        cell_workspace,ndofs = rt_setup_dual_basis_interior_n_cube(fe,offset)
        primal_basis = rt_primal_basis_n_cube(fe)
    elseif is_simplex(fe.domain)
        cell_workspace,ndofs = rt_setup_dual_basis_interior_simplex(fe,offset)
        primal_basis = rt_primal_basis_simplex(fe)
    else
        error("case not implemented")
    end
    face_dof_point_moment = map(workspace->workspace.moments,face_workspace)
    push!(face_dof_point_moment,cell_workspace.moments)
    face_point_x = map(workspace->workspace.points,face_workspace)
    push!(face_point_x,cell_workspace.points)
    dual_basis_nested = map(face_dof_point_moment,face_point_x) do dof_point_moment,point_x
        map(dof_point_moment) do point_moment
            npoints = length(point_x)
            v -> sum(1:npoints) do point
                moment = point_moment[point]
                x = point_x[point]
                moment⋅v(x)
            end
        end
    end
    dual_basis = reduce(vcat,dual_basis_nested)
    @assert length(dual_basis) == length(primal_basis)
    (;face_workspace,cell_workspace,primal_basis,dual_basis,ndofs)
end

function rt_primal_basis_n_cube(fe)
    D = num_dims(fe.domain)
    k = fe.order
    nested = map(1:D) do d
        ranges = ntuple(d2-> d==d2 ? (0:k+1) : (0:k) ,Val(D))
        exponents_list = map(Tuple,CartesianIndices(ranges))[:]
        map(exponents_list) do exponents
            x -> ntuple(Val(D)) do d3
                v = prod( x.^exponents )
                d==d3 ? v : zero(v)
            end |> SVector
        end
    end
    reduce(vcat,nested)
end

function rt_primal_basis_simplex(fe)
    D = num_dims(fe.domain)
    k = fe.order
    ranges = ntuple(d2-> (0:k) ,Val(D))
    exponents_list_all = map(Tuple,CartesianIndices(ranges))[:]
    exponents_list = filter(exponents->sum(exponents)<=k,exponents_list_all)
    nested = map(1:D) do d
        map(exponents_list) do exponents
            x -> ntuple(Val(D)) do d3
                v = prod( x.^exponents )
                d==d3 ? v : zero(v)
            end |> SVector
        end
    end
    list1 = reduce(vcat,nested)
    exponents_list = filter(exponents->sum(exponents)==k,exponents_list_all)
    list2 = map(exponents_list) do exponents
        x -> begin
            v = prod( x.^exponents )
            x*v
        end
    end
    # TODO this array does not have concrete eltype
    vcat(list1,list2)
end

function rt_setup_dual_basis_boundary(fe)
    # TODO this function can perhaps be implemented
    # with a higher level API
    k = fe.order # TODO k
    cell = fe.domain
    D = num_dims(cell)
    cell_boundary = mesh(cell)
    face_to_rid = face_reference_id(cell_boundary,D-1)
    rid_to_refface = reference_spaces(cell_boundary,D-1)
    rid_to_fe = map(rid_to_refface) do refface
        lagrange_space(domain(refface),k) # TODO k
    end
    rid_to_quad = map(rid_to_fe) do fe
        quadrature(domain(fe),2*k) # TODO 2*K
    end
    rid_to_xs = map(coordinates,rid_to_quad)
    rid_to_ws = map(weights,rid_to_quad)
    rid_to_tabface = map(rid_to_refface,rid_to_xs) do refface,xs
        tabulator(refface)(value,xs)
    end
    rid_to_tabface_grad = map(rid_to_refface,rid_to_xs) do refface,xs
        tabulator(refface)(ForwardDiff.gradient,xs)
    end
    rid_to_tabfe = map(rid_to_fe,rid_to_xs) do fe,xs
        tabulator(fe)(value,xs)
    end
    rid_to_permutations = map(fe->node_permutations(fe),rid_to_fe)
    face_to_n = GT.outward_normals(cell_boundary)
    face_to_nodes = face_nodes(cell_boundary,D-1)
    node_to_x = node_coordinates(cell_boundary)
    nfaces = length(face_to_n)
    offset = 0
    face_workspace = map(1:nfaces) do face
        rid = face_to_rid[face]
        tabfe = rid_to_tabfe[rid]
        tabface = rid_to_tabface[rid]
        tabface_grad = rid_to_tabface_grad[rid]
        npoints,ndofs = size(tabfe)
        nodes = face_to_nodes[face]
        nnodes = length(nodes)
        n = face_to_n[face]
        point_to_w = rid_to_ws[rid]
        point_to_x = map(1:npoints) do point
            sum(1:nnodes) do node
                x = node_to_x[nodes[node]]
                coeff = tabface[point,node]
                coeff*x
            end
        end
        point_to_dV = map(1:npoints) do point
            w = point_to_w[point]
            J = sum(1:nnodes) do node
                x = node_to_x[nodes[node]]
                coeff = tabface_grad[point,node]
                outer(x,coeff)
            end
            change_of_measure(J)*w
        end
        #scaling = sum(point_to_dV)
        moments = map(1:ndofs) do dof
            map(1:npoints) do point
                x = point_to_x[point]
                dV = point_to_dV[point]
                s = tabfe[point,dof]
                n*(s*dV)#/scaling)
            end
        end
        points = point_to_x
        permutations = rid_to_permutations[rid]
        ids = collect(Int32,1:ndofs) .+ offset
        offset += ndofs
        (;moments,points,ids,permutations)
    end
    face_workspace,offset
end

function rt_setup_dual_basis_interior_n_cube(fe,offset)
    k = fe.order
    cell = fe.domain
    quad = quadrature(cell,2*k) # TODO 2*k
    point_to_x = coordinates(quad)
    point_to_dV = weights(quad)
    npoints = length(point_to_x)
    D = num_dims(cell)
    nested = map(1:D) do d
        ranges = ntuple(d2-> d==d2 ? (0:k-1) : (0:k) ,Val(D))
        exponents_list = map(Tuple,CartesianIndices(ranges))[:]
        map(exponents_list) do exponents
            #scaling = sum(point_to_dV)
            map(1:npoints) do point
                x = point_to_x[point]
                dV = point_to_dV[point]
                s = ntuple(Val(D)) do d3
                    si = prod( x.^exponents )
                    d==d3 ? si : zero(si)
                end |> SVector
                s*dV#/scaling
            end |> Ref
        end
    end
    moments = map(r->r[],reduce(vcat,nested))
    nmoments = length(moments)
    permutations = [collect(Int32,1:nmoments)]
    points = point_to_x
    ids = collect(Int32,1:nmoments).+offset
    offset += nmoments
    (;moments,points,ids,permutations), offset
end

function rt_setup_dual_basis_interior_simplex(fe,offset)
    k = fe.order - 1
    cell = fe.domain
    quad = quadrature(cell,2*k) # TODO 2*k
    point_to_x = coordinates(quad)
    point_to_dV = weights(quad)
    npoints = length(point_to_x)
    D = num_dims(cell)
    nested = map(1:D) do d
        ranges = ntuple(d2->(0:k),Val(D))
        exponents_list = map(Tuple,CartesianIndices(ranges))[:]
        exponents_list = filter(exponents->sum(exponents)<=(k),exponents_list)
        map(exponents_list) do exponents
            #scaling = sum(point_to_dV)
            map(1:npoints) do point
                x = point_to_x[point]
                dV = point_to_dV[point]
                s = ntuple(Val(D)) do d3
                    si = prod( x.^exponents )
                    d==d3 ? si : zero(si)
                end |> SVector
                s*dV#/scaling
            end |> Ref
        end
    end
    moments = map(r->r[],reduce(vcat,nested))
    nmoments = length(moments)
    permutations = [collect(Int32,1:nmoments)]
    points = point_to_x
    ids = collect(Int32,1:nmoments).+offset
    offset += nmoments
    (;moments,points,ids,permutations), offset
end

function num_dofs(fe::RaviartThomasFaceSpace)
    fe.workspace.ndofs
end

function primal_basis(fe::RaviartThomasFaceSpace)
    fe.workspace.primal_basis
end

function dual_basis(fe::RaviartThomasFaceSpace)
    fe.workspace.dual_basis
end

function face_own_dofs(fe::RaviartThomasFaceSpace,d)
    D = num_dims(fe.domain)
    if D == d
        [fe.workspace.cell_workspace.ids]
    elseif D-1 == d
        map(workspace->workspace.ids,fe.workspace.face_workspace)
    else
        [ Int32[] for _ in 1:num_faces(mesh(fe.domain),d)]
    end
end

function face_dofs(fe::RaviartThomasFaceSpace,d)
    face_own_dofs(fe,d)
end

function face_own_dof_permutations(fe::RaviartThomasFaceSpace,d)
    D = num_dims(fe.domain)
    if D == d
        [fe.workspace.cell_workspace.permutations]
    elseif D-1 == d
        map(workspace->workspace.permutations,fe.workspace.face_workspace)
    else
        [ [Int32[]] for _ in 1:num_faces(mesh(fe.domain),d)]
    end
end

function raviart_thomas_space(domain::AbstractMeshDomain,order::Integer;conformity=:default,dirichlet_boundary=nothing)
    mesh = domain |> GT.mesh
    cell_to_Dface = domain |> GT.faces
    D = domain |> GT.num_dims
    Dface_to_ctype = GT.face_reference_id(mesh,D)
    cell_to_ctype = Dface_to_ctype[cell_to_Dface]
    ctype_to_refface = GT.reference_spaces(mesh,D)
    ctype_to_geometry = map(GT.domain,ctype_to_refface)
    ctype_to_reffe = map(ctype_to_geometry) do geometry
        raviart_thomas_space(geometry,order)
    end
    workspace = nothing
    RaviartThomasMeshSpace(
        mesh,
        domain,
        order,
        conformity,
        dirichlet_boundary,
        cell_to_ctype,
        ctype_to_reffe,
        workspace) |> setup_space
end

struct RaviartThomasMeshSpace{M,A,B,C,D,E,F} <: AbstractSpace{M}
    mesh::M
    domain::A
    order::B
    conformity::Symbol
    dirichlet_boundary::C
    face_reference_id::D
    reference_spaces::E
    workspace::F
end

domain(a::RaviartThomasMeshSpace) = a.domain
order(a::RaviartThomasMeshSpace) = a.order
conformity(a::RaviartThomasMeshSpace) = a.conformity
dirichlet_boundary(a::RaviartThomasMeshSpace) = a.dirichlet_boundary
face_reference_id(a::RaviartThomasMeshSpace) = a.face_reference_id
reference_spaces(a::RaviartThomasMeshSpace) = a.reference_spaces
workspace(a::RaviartThomasMeshSpace) = a.workspace

function replace_workspace(space::RaviartThomasMeshSpace,workspace)
    GT.RaviartThomasMeshSpace(
        space.mesh,
        space.domain,
        space.order,
        space.conformity,
        space.dirichlet_boundary,
        space.face_reference_id,
        space.reference_spaces,
        workspace
    )
end

function push_forward(space::RaviartThomasMeshSpace,qty)
    qty_flipped = flip_sign(space,qty)
    mesh = GT.mesh(space)
    D = num_dims(mesh)
    phi = physical_map(mesh,D)
    call(qty_flipped,phi) do q,phi
        x->begin
            J = ForwardDiff.jacobian(phi,x)
            (J/change_of_measure(J))*q(x)
        end
    end
end

function pull_back(space::RaviartThomasMeshSpace,qty)
    qty_flipped = flip_sign(space,qty)
    mesh = GT.mesh(space)
    D = num_dims(mesh)
    phi = physical_map(mesh,D)
    call(qty_flipped,phi) do σ,phi
        u -> begin
            v = x->begin
                J = ForwardDiff.jacobian(phi,x)
                (J/change_of_measure(J))\u(x)
            end
            σ(v)
        end
    end
end

function flip_sign(space,qty)
    D = num_dims(mesh(space))
    quantity() do index
        face_dof_flip = get_symbol!(index,sign_flip_accessor(space),"face_dof_flip")
        face = face_index(index,D)
        t = term(qty,index)
        dof = dof_index(t)
        expr = @term begin
            flip = $face_dof_flip($face)($dof)
            s = $(expression(t))
            x->flip*s(x)
        end
        expr_term(D,expr,x->1*prototype(t)(x),index)
    end
end

function sign_flip_criterion(cell,cells)
    if cell == cells[1]
        1
    else
        -1
    end
end

function sign_flip_accessor(space::RaviartThomasMeshSpace)
    D = num_dims(mesh(space))
    cell_rid = face_reference_id(space)
    rid_fe = reference_spaces(space)
    rid_ldof_lface = map(rid_fe) do fe
        nldofs = num_dofs(fe)
        ldof_lface = zeros(Int32,nldofs)
        lface_to_ldofs = face_own_dofs(fe,D-1)
        for (lface,ldofs) in enumerate(lface_to_ldofs)
            ldof_lface[ldofs] .= lface
        end
        ldof_lface
    end
    topo = topology(mesh(space))
    cell_faces = face_incidence(topo,D,D-1)
    face_cells = face_incidence(topo,D-1,D)
    function cell_dof_flip(cell)
        rid = cell_rid[cell]
        ldof_lface = rid_ldof_lface[rid]
        lface_to_face = cell_faces[cell]
        function dof_flip(ldof)
            lface = ldof_lface[ldof]
            if lface == 0
                return 1
            end
            face = lface_to_face[lface]
            cells = face_cells[face]
            sign_flip_criterion(cell,cells)
        end
    end
end

function cartesian_product(spaces::AbstractSpace...)
    mesh = GT.mesh(first(spaces))
    CartesianProductSpace(mesh,spaces)
end

struct CartesianProductSpace{M,A} <: GT.AbstractSpace{M}
    mesh::M
    spaces::A
end

function mesh(space::CartesianProductSpace)
    GT.mesh(first(space.spaces))
end

function LinearAlgebra.:×(a::AbstractSpace,b::AbstractSpace)
    cartesian_product(a,b)
end

function LinearAlgebra.:×(a::CartesianProductSpace,b::AbstractSpace)
    cartesian_product(a.spaces...,b)
end

function LinearAlgebra.:×(a::AbstractSpace,b::CartesianProductSpace)
    cartesian_product(a,b.spaces...)
end

function fields(a::CartesianProductSpace)
    a.spaces
end

function field(u::CartesianProductSpace,field::Integer)
    u.spaces[field]
end

function num_fields(a::CartesianProductSpace)
    length(a.spaces)
end

function domain(a::CartesianProductSpace)
    domains = map(GT.domain,a.spaces)
    domain = first(domains)
    if all(dom -> dom == domain,domains)
        domain
    else
        # TODO generalize the concept of AbstractQuantity
        # by introducing PiecewiseQuantity and PiecewiseDomain
        # include also domain id in index
        # I am not sure about this. Not needed in practice
        error("Case not yet implemented")
    end
end

#function domain(a::CartesianProductSpace,field)
#    GT.domain(GT.field(a,field))
#end

function face_dofs(a::CartesianProductSpace)
    error("Not implemented, not needed in practice")
end

#function face_dofs(a::CartesianProductSpace,field)
#    face_dofs(a.spaces[field])
#end

function free_dofs(a::CartesianProductSpace)
    nfields = GT.num_fields(a)
    map(1:nfields) do field
        free_dofs(GT.field(a,field))
    end |> BRange
end

#function free_dofs(a::CartesianProductSpace,field)
#    f = GT.field(a,field)
#    free_dofs(f)
#end

function dirichlet_dofs(a::CartesianProductSpace)
    nfields = GT.num_fields(a)
    map(1:nfields) do field
        dirichlet_dofs(GT.field(a,field))
    end |> BRange
end

#function dirichlet_dofs(a::CartesianProductSpace,field)
#    f = GT.field(a,field)
#    dirichlet_dofs(f)
#end

function form_argument_quantity(a::CartesianProductSpace,axis)
    fields = ntuple(identity,GT.num_fields(a))
    map(fields) do field
        GT.form_argument_quantity(GT.field(a,field),axis,field)
    end
end

function discrete_field_quantity(a::CartesianProductSpace,free_vals,diri_vals)
    nothing
end

function shape_function_quantity(a::CartesianProductSpace,field)
    error("Not implemented yet. Not needed in practice.")
end

function dual_basis_quantity(a::CartesianProductSpace)
    error("Not implemented yet. Not needed in practice.")
end

function dual_basis_quantity(a::CartesianProductSpace,field)
    error("Not implemented yet. Not needed in practice.")
end

function shape_function_accessor(f::typeof(value),space::AbstractSpace,measure::Measure)
    mesh = measure.mesh
    dom = measure.domain
    @assert !is_reference_domain(dom)
    d = num_dims(dom)
    @assert num_dims(domain(space)) == d
    # TODO assumes same reference elems for integration and for interpolation
    rid_to_point_to_x = map(coordinates,reference_quadratures(measure))
    rid_to_tab = map(rid_to_point_to_x,reference_spaces(space)) do point_to_x, refface
        tabulator(refface)(f,point_to_x)
    end
    face_to_rid = face_reference_id(mesh,d)
    sface_to_face = faces(dom)
    function face_point_dof_s(sface)
        face = sface_to_face[sface]
        rid = face_to_rid[face]
        tab = rid_to_tab[rid]
        function point_dof_s(point,J=nothing)
            function dof_s(dof)
                tab[point,dof]
            end
        end
    end
end

function shape_function_accessor(f::typeof(ForwardDiff.gradient),space::AbstractSpace,measure::Measure)
    mesh = measure.mesh
    dom = measure.domain
    @assert !is_reference_domain(dom)
    d = num_dims(dom)
    @assert num_dims(domain(space)) == d
    rid_to_point_to_x = map(coordinates,reference_quadratures(measure))
    rid_to_tab = map(rid_to_point_to_x,reference_spaces(space)) do point_to_x, refface
        tabulator(refface)(f,point_to_x)
    end
    rid_to_dof_to_phys = map(rid_to_tab) do tab
        T = eltype(tab)
        ndofs = size(tab,2)
        # This assumes that the jacobian is a square matrix
        zeros(T,ndofs)
    end
    face_to_rid = face_reference_id(mesh,d)
    sface_to_face = faces(dom)
    function face_point_dof_s(sface)
        face = sface_to_face[sface]
        rid = face_to_rid[face]
        tab = rid_to_tab[rid]
        dof_to_phys = rid_to_dof_to_phys[rid]
        ndofs = length(dof_to_phys)
        function point_dof_s(point,J)
            # NB you cannot evaluate this at more than one point at once
            for dof in 1:ndofs
                dof_to_phys[dof] = transpose(J)\tab[point,dof]
            end
            function dof_s(dof)
                dof_to_phys[dof]
            end
        end
    end
end

function dofs_accessor(space::AbstractSpace,dom::AbstractDomain)
    @assert num_fields(space) == 1
    d = num_dims(dom)
    dom2 = domain(space)
    @assert num_dims(dom2) == d
    sface_to_face = faces(dom)
    face_to_dofs_ = face_dofs(space)
    function face_to_dofs(sface)
        face = sface_to_face[sface]
        dofs = face_to_dofs_[face]
        dofs
    end
end

function discrete_field_accessor(f,uh::DiscreteField,measure::Measure)
    dom = domain(measure)
    space = GT.space(uh)
    d = num_dims(dom)
    dom2 = domain(space)
    @assert num_dims(dom2) == d
    sface_to_face = faces(dom)
    face_to_dofs = face_dofs(space)
    face_to_rid = face_reference_id(space)
    free_vals = free_values(uh)
    diri_vals = dirichlet_values(uh)
    face_point_dof_s = shape_function_accessor(f,space,measure)
    function face_point_val(sface)
        face = sface_to_face[sface]
        dofs = face_to_dofs[face]
        rid = face_to_rid[face]
        point_dof_s = face_point_dof_s(sface)
        ndofs = length(dofs)
        function point_val(point,J)
            dof_s = point_dof_s(point,J)
            sum(1:ndofs) do i
                dof = dofs[i]
                s = dof_s(i)
                if dof > 0
                    v = free_vals[dof]
                else
                    v = diri_vals[-dof]
                end
                v*s
            end
        end
    end
end

function dirichlet_accessor(uh::DiscreteField,dom::AbstractDomain)
    space = GT.space(uh)
    d = num_dims(dom)
    dom2 = domain(space)
    @assert num_dims(dom2) == d
    sface_to_face = faces(dom)
    face_to_dofs = face_dofs(space)
    rid_to_u = map(reference_spaces(space)) do fe
        T = eltype(dirichlet_values(uh))
        zeros(T,num_dofs(fe))
    end
    face_to_rid = face_reference_id(space)
    diri_vals = dirichlet_values(uh)
    function face_dirichlet!(sface)
        face = sface_to_face[sface]
        dofs = face_to_dofs[face]
        rid = face_to_rid[face]
        u = rid_to_u[rid]
        fill!(u,zero(eltype(u)))
        for (i,dof) in enumerate(dofs)
            if dof < 0
                u[i] = diri_vals[-dof]
            end
        end
        function dirichlet!(A,b)
            m,n = size(A)
            z = zero(eltype(b))
            for i in 1:m
                bi = z
                for j in 1:n
                    bi += A[i,j]*u[j]
                end
                b[i] -= bi
            end
        end
    end
end

function generate_workspace(space::AbstractSpace{<:PMesh})
    D = num_dims(domain(space))
    mesh = GT.mesh(space)
    spaces = partition(space)
    p_state_1 = map(setup_space_local_step_1,spaces)
    p_d_n_oddofs = map(state->state.d_n_oddofs,p_state_1)
    d_p_n_oddofs = tuple_of_arrays(p_d_n_oddofs)
    d_n_gddofs = map(d->sum(d_p_n_oddofs[d+1]),0:D)
    d_first_gdof = zeros(Int,D+1)
    d_first_gdof[1] = 1
    for d in 1:D
        d_first_gdof[d+1]=d_first_gdof[d]+d_n_gddofs[d]
    end
    d_p_doffset = map(0:D) do d
        p_n_oddofs = d_p_n_oddofs[d+1]
        scan(+,p_n_oddofs,type=:exclusive,init=d_first_gdof[d+1])
    end
    p_d_doffset = array_of_tuples(d_p_doffset)
    p_state_2 = map(setup_space_local_step_2,p_state_1,p_d_doffset)
    p_d_dface_dof_goffset = map(state->state.d_dface_dof_goffset,p_state_2)
    d_p_dface_dof_goffset = tuple_of_arrays(p_d_dface_dof_goffset)
    d_p_dface_ids = map(d->face_partition(mesh,d),0:D)
    for d in 0:D
        p_dface_dof_goffset = d_p_dface_dof_goffset[d+1]
        p_dface_ids = d_p_dface_ids[d+1]
        v = PVector(p_dface_dof_goffset,p_dface_ids)
        wait(consistent!(v))
    end
    ngdofs = d_first_gdof[end]+d_n_gddofs[end]-1
    p_ngdofs = map(s->ngdofs,spaces)
    p_state_3 = map(setup_space_local_step_3,p_state_2,p_ngdofs)
    p_dof_partition = map(state->state.dof_local_indices,p_state_3)
    p_dof_isfree = map(state->state.dof_isfree,p_state_3)
    gdof_isfree = PVector(p_dof_isfree,p_dof_partition)
    wait(consistent!(gdof_isfree))
    gdof_isdiri = .!(gdof_isfree)
    free_gdof_dof, gdof_free_dof = find_local_indices(gdof_isfree)
    diri_gdof_dof, gdof_diri_dof = find_local_indices(gdof_isdiri)
    p_dof_free_dof = partition(gdof_free_dof)
    p_dof_diri_dof = partition(gdof_diri_dof)
    free_dofs = axes(free_gdof_dof,1)
    diri_dofs = axes(diri_gdof_dof,1)
    p_free_dofs_ids = partition(free_dofs)
    p_diri_dofs_ids = partition(diri_dofs)
    p_state_4 = map(setup_space_local_step_4,p_state_3,p_dof_free_dof,p_dof_diri_dof,p_free_dofs_ids,p_diri_dofs_ids)
    space_partition = map(state->state.space_with_setup,p_state_4)
    p_diri_dof_location = map(state->state.dirichlet_dof_location,p_state_4)
    dirichlet_dof_location = PVector(p_diri_dof_location,p_diri_dofs_ids)
    workspace = (;space_partition,free_dofs,dirichlet_dofs=diri_dofs,dirichlet_dof_location)
end

function setup_space_local_step_1(space)
    domain = space |> GT.domain
    D = GT.num_dims(domain)
    cell_Dface = domain |> GT.faces
    mesh = domain |> GT.mesh
    topology = mesh |> GT.topology
    ctype_reference_fe = space |> GT.reference_spaces
    cell_ctype = space |> GT.face_reference_id
    d_dface_dof_goffset = ntuple(t->zeros(Int32,GT.num_faces(topology,t-1)),D+1)
    d_ctype_ldface_own_dofs = map(d->GT.reference_face_own_dofs(space,d),0:D)
    d_ctype_ldface_num_own_dofs = map(d->map(ldface_own_dofs->length.(ldface_own_dofs),d_ctype_ldface_own_dofs[d+1]),0:D)
    d_Dface_dfaces = map(d->face_incidence(topology,D,d),0:D)
    ncells = length(cell_ctype)
    d_n_oddofs = ntuple(D+1) do t
        d = t-1
        dof_offset = 0
        ctype_ldface_num_own_dofs = d_ctype_ldface_num_own_dofs[d+1]
        dface_dof_goffset = d_dface_dof_goffset[d+1]
        Dface_dfaces = d_Dface_dfaces[d+1]
        ndfaces = length(dface_dof_goffset)
        face_ids = face_local_indices(mesh,d)
        part = part_id(face_ids)
        dface_owner = local_to_owner(face_ids)
        for cell in 1:ncells
            ctype = cell_ctype[cell]
            Dface = cell_Dface[cell]
            ldface_num_own_dofs = ctype_ldface_num_own_dofs[ctype]
            ldface_dface = Dface_dfaces[Dface]
            nldfaces = length(ldface_num_own_dofs)
            for ldface in 1:nldfaces
                num_own_dofs = ldface_num_own_dofs[ldface]
                dface = ldface_dface[ldface]
                dface_dof_goffset[dface] = num_own_dofs
            end
        end
        for dface in 1:ndfaces
            owner = dface_owner[dface]
            if owner != part
                continue
            end
            num_own_dofs = dface_dof_goffset[dface]
            dface_dof_goffset[dface] = dof_offset
            dof_offset += num_own_dofs
        end
        dof_offset
    end
    (;space,d_n_oddofs,d_dface_dof_goffset)
end

function setup_space_local_step_2(state,d_doffset)
    (;space,d_n_oddofs,d_dface_dof_goffset) = state
    domain = space |> GT.domain
    D = GT.num_dims(domain)
    mesh = domain |> GT.mesh
    for d in 0:D
        dface_dof_goffset = d_dface_dof_goffset[d+1]
        ndfaces = length(dface_dof_goffset)
        dof_offset = d_doffset[d+1]
        face_ids = face_local_indices(mesh,d)
        part = part_id(face_ids)
        dface_owner = local_to_owner(face_ids)
        for dface in 1:ndfaces
            owner = dface_owner[dface]
            if owner != part
                continue
            end
            dface_dof_goffset[dface] += dof_offset
        end
    end
    state
end

function setup_space_local_step_3(state,ngdofs)
    space = state.space
    dirichlet_boundary = GT.dirichlet_boundary(space)
    state1 = setup_space_local_step_3_1(state,ngdofs)
    state2 = setup_space_local_step_3_2(state1,dirichlet_boundary)
    state2
end

function setup_space_local_step_3_1(state,ngdofs)
    (;space,d_n_oddofs,d_dface_dof_goffset) = state
    domain = space |> GT.domain
    D = GT.num_dims(domain)
    cell_Dface = domain |> GT.faces
    mesh = domain |> GT.mesh
    topology = mesh |> GT.topology
    ctype_reference_fe = space |> GT.reference_spaces
    cell_ctype = space |> GT.face_reference_id
    d_dface_dof_offset = map(d->zeros(Int32,GT.num_faces(topology,d)),0:D)
    d_ctype_ldface_own_dofs = map(d->GT.reference_face_own_dofs(space,d),0:D)
    d_ctype_ldface_pindex_perm = map(d->GT.reference_face_own_dof_permutations(space,d),0:D)
    d_ctype_ldface_num_own_dofs = map(d->map(ldface_own_dofs->length.(ldface_own_dofs),d_ctype_ldface_own_dofs[d+1]),0:D)
    d_ctype_ldface_dofs = map(d->map(fe->GT.face_dofs(fe,d),ctype_reference_fe),0:D)
    d_Dface_dfaces = map(d->face_incidence(topology,D,d),0:D)
    d_Dface_ldface_pindex = map(d->face_permutation_ids(topology,D,d),0:D)
    ctype_num_dofs = map(GT.num_dofs,ctype_reference_fe)
    ncells = length(cell_ctype)
    nDfaces = num_faces(topology,D)
    dof_offset = 0
    for d in 0:D
        ctype_ldface_num_own_dofs = d_ctype_ldface_num_own_dofs[d+1]
        dface_dof_offset = d_dface_dof_offset[d+1]
        Dface_dfaces = d_Dface_dfaces[d+1]
        ndfaces = length(dface_dof_offset)
        for cell in 1:ncells
            ctype = cell_ctype[cell]
            Dface = cell_Dface[cell]
            ldface_num_own_dofs = ctype_ldface_num_own_dofs[ctype]
            ldface_dface = Dface_dfaces[Dface]
            nldfaces = length(ldface_num_own_dofs)
            for ldface in 1:nldfaces
                num_own_dofs = ldface_num_own_dofs[ldface]
                dface = ldface_dface[ldface]
                dface_dof_offset[dface] = num_own_dofs
            end
        end
        for dface in 1:ndfaces
            num_own_dofs = dface_dof_offset[dface]
            dface_dof_offset[dface] = dof_offset
            dof_offset += num_own_dofs
        end
    end
    ndofs = dof_offset
    dof_gdof = zeros(Int,ndofs)
    dof_owner = zeros(Int,ndofs)
    Dface_ptrs = zeros(Int32,nDfaces+1)
    for cell in 1:ncells
        ctype = cell_ctype[cell]
        num_dofs = ctype_num_dofs[ctype]
        Dface = cell_Dface[cell]
        Dface_ptrs[Dface+1] = num_dofs
    end
    length_to_ptrs!(Dface_ptrs)
    ndata = Dface_ptrs[end]-1
    Dface_dofs = JaggedArray(zeros(Int32,ndata),Dface_ptrs)
    for d in 0:D
        Dface_dfaces = d_Dface_dfaces[d+1]
        ctype_ldface_own_ldofs = d_ctype_ldface_own_dofs[d+1]
        ctype_ldface_pindex_perm = d_ctype_ldface_pindex_perm[d+1]
        dface_dof_offset = d_dface_dof_offset[d+1]
        dface_dof_goffset = d_dface_dof_goffset[d+1]
        Dface_ldface_pindex = d_Dface_ldface_pindex[d+1]
        face_ids = face_local_indices(mesh,d)
        part = part_id(face_ids)
        dface_owner = local_to_owner(face_ids)
        for cell in 1:ncells
            ctype = cell_ctype[cell]
            Dface = cell_Dface[cell]
            Dface = cell_Dface[cell]
            ldof_dof = Dface_dofs[Dface]
            ldface_dface = Dface_dfaces[Dface]
            ldface_own_ldofs = ctype_ldface_own_ldofs[ctype]
            ldface_pindex_perm = ctype_ldface_pindex_perm[ctype]
            ldface_pindex = Dface_ldface_pindex[Dface]
            nldfaces = length(ldface_dface)
            for ldface in 1:nldfaces
                dface = ldface_dface[ldface]
                own_ldofs = ldface_own_ldofs[ldface]
                dof_offset = dface_dof_offset[dface]
                dof_goffset = dface_dof_goffset[dface]
                pindex_perm = ldface_pindex_perm[ldface]
                pindex = ldface_pindex[ldface]
                perm = pindex_perm[pindex]
                n_own_dofs = length(own_ldofs)
                owner = dface_owner[dface]
                for i in 1:n_own_dofs
                    j = perm[i]
                    dof = j + dof_offset
                    gdof = j + dof_goffset
                    dof_gdof[dof] = gdof
                    dof_owner[dof] = owner
                    own_dof = own_ldofs[i]
                    ldof_dof[own_dof] = dof
                end
            end
        end
    end
    part = part_id(face_local_indices(mesh,0))
    dof_local_indices = PartitionedArrays.LocalIndices(ngdofs,part,dof_gdof,dof_owner)
    (;ndofs,Dface_dofs,dof_local_indices,state...)
end

function setup_space_local_step_3_2(state,dirichlet_boundary::Nothing)
    (;ndofs,space) = state
    dof_isfree = fill(true,ndofs)
    dof_location = fill(1,ndofs)
    (;dof_isfree,dof_location,state...)
end

function setup_space_local_step_3_2(state,dirichlet_boundary::AbstractDomain)
    (;ndofs,Dface_dofs,space) = state
    domain = space |> GT.domain
    D = GT.num_dims(domain)
    cell_Dface = domain |> GT.faces
    mesh = domain |> GT.mesh
    topology = mesh |> GT.topology
    ctype_reference_fe = space |> GT.reference_spaces
    cell_ctype = space |> GT.face_reference_id
    d_ctype_ldface_dofs = map(d->map(fe->GT.face_dofs(fe,d),ctype_reference_fe),0:D)
    d_Dface_dfaces = map(d->face_incidence(topology,D,d),0:D)
    dof_isfree = fill(true,ndofs)
    dof_location = fill(1,ndofs)
    N = GT.num_dims(dirichlet_boundary)
    nNfaces = num_faces(topology,N)
    Nface_tag = zeros(Int32,nNfaces)
    mesh = dirichlet_boundary |> GT.mesh
    Nface_tag[GT.faces(dirichlet_boundary)] .= 1
    ncells = length(cell_Dface)
    let d = N
        Dface_dfaces = d_Dface_dfaces[d+1]
        ctype_ldface_ldofs = d_ctype_ldface_dofs[d+1]
        for cell in 1:ncells
            ctype = cell_ctype[cell]
            Dface = cell_Dface[cell]
            ldof_dof = Dface_dofs[Dface]
            ldface_dface = Dface_dfaces[Dface]
            ldface_ldofs = ctype_ldface_ldofs[ctype]
            nldfaces = length(ldface_dface)
            dofs = Dface_dofs[Dface]
            for ldface in 1:nldfaces
                ldofs = ldface_ldofs[ldface]
                dface = ldface_dface[ldface]
                Nface = dface
                tag = Nface_tag[Nface]
                if tag == 0
                    continue
                end
                dof_isfree[view(dofs,ldofs)] .= false
            end
        end
    end
    (;dof_isfree,dof_location,state...)
end

function setup_space_local_step_4(state,dof_free_dof,dof_diri_dof,free_dofs_ids,diri_dofs_ids)
    (;Dface_dofs,dof_isfree,ndofs,dof_location,space) = state
    f = dof -> begin
        isfree = dof_isfree[dof]
        if isfree
            dof2 = dof_free_dof[dof]
        else
            dof2 = -dof_diri_dof[dof]
        end
    end
    data = Dface_dofs.data
    data .= f.(data)
    nfree = local_length(free_dofs_ids)
    ndiri = local_length(diri_dofs_ids)
    dirichlet_dof_location = zeros(Int32,ndiri)
    for dof in 1:ndofs
        if !dof_isfree[dof]
            diri_dof = dof_diri_dof[dof]
            location = dof_location[diri_dof]
            dirichlet_dof_location[diri_dof] = location
        end
    end
    free_dofs = Base.OneTo(nfree)
    dirichlet_dofs = Base.OneTo(ndiri)
    face_dofs = Dface_dofs
    free_dof_local_indices = free_dofs_ids
    dirichlet_dof_local_indices = diri_dofs_ids
    workspace = (;face_dofs,free_dofs,dirichlet_dofs,dirichlet_dof_location,free_dof_local_indices,dirichlet_dof_local_indices)
    space_with_setup = replace_workspace(space,workspace)
    (;space_with_setup,dirichlet_dof_location,state...)
end


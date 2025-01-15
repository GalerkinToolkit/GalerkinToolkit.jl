

domain(space::AbstractSpace,field) = domain(space)
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


function free_dofs(a::AbstractSpace,field)
    @assert field == 1
    free_dofs(a)
end

function free_dofs(a::AbstractSpace)
    nfree = length(first(free_and_dirichlet_dofs(a)))
    Base.OneTo(nfree)
end

function dirichlet_dofs(a::AbstractSpace,field)
    @assert field == 1
    dirichlet_dofs(a)
end

function dirichlet_dofs(a::AbstractSpace)
    ndiri = length(last(free_and_dirichlet_dofs(a)))
    Base.OneTo(ndiri)
end

function setup_space(space::AbstractSpace)
    if GT.workspace(space) !== nothing
        return space
    end
    state = generate_dof_ids(space)
    face_dofs = state.cell_to_dofs
    free_and_dirichlet_dofs = state.free_and_dirichlet_dofs
    dirichlet_dof_location = state.dirichlet_dof_location
    workspace = (;face_dofs,free_and_dirichlet_dofs,dirichlet_dof_location)
    replace_workspace(space,workspace)
end

function face_dofs(space::AbstractSpace)
    if workspace(space) !== nothing
        return workspace(space).face_dofs
    end
    state = generate_dof_ids(space)
    state.cell_to_dofs # TODO rename face_dofs ?
end

function face_dofs(space::AbstractSpace,field)
    @assert field == 1
    face_dofs(space)
end

function free_and_dirichlet_dofs(V::AbstractSpace)
    if workspace(V) !== nothing
        return workspace(V).free_and_dirichlet_dofs
    end
    state = generate_dof_ids(V)
    state.free_and_dirichlet_dofs
end

function dirichlet_dof_location(V::AbstractSpace)
    if workspace(V) !== nothing
        return workspace(V).dirichlet_dof_location
    end
    state = generate_dof_ids(V)
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
    cell_to_ptrs = zeros(Int32,ncells+1)
    for cell in 1:ncells
        ctype = cell_to_ctype[cell]
        num_dofs = ctype_to_num_dofs[ctype]
        cell_to_ptrs[cell+1] = num_dofs
    end
    length_to_ptrs!(cell_to_ptrs)
    ndata = cell_to_ptrs[end]-1
    cell_to_dofs = JaggedArray(zeros(Int32,ndata),cell_to_ptrs)
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
            ldof_to_dof = cell_to_dofs[cell]
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
    (;ndofs,cell_to_dofs,d_to_Dface_to_dfaces,
     d_to_ctype_to_ldface_to_dofs,d_to_ndfaces,
     cell_to_ctype,cell_to_Dface)
end

function generate_dof_ids_step_2(space,state,dirichlet_boundary::Nothing)
    (;ndofs,cell_to_dofs) = state
    dof_to_tag = zeros(Int32,ndofs)
    free_and_dirichlet_dofs = GT.partition_from_mask(i->i==0,dof_to_tag)
    dirichlet_dof_location = zeros(Int32,0)
    (;cell_to_dofs, free_and_dirichlet_dofs,dirichlet_dof_location)
end

function generate_dof_ids_step_2(space,state,dirichlet_boundary::AbstractDomain)
    (;ndofs,cell_to_dofs,d_to_Dface_to_dfaces,
     d_to_ctype_to_ldface_to_dofs,
     d_to_ndfaces,cell_to_ctype,cell_to_Dface) = state
    dof_to_tag = zeros(Int32,ndofs)
    N = GT.num_dims(dirichlet_boundary)
    #physical_names = dirichlet_boundary |> GT.physical_names
    Nface_to_tag = zeros(Int32,d_to_ndfaces[N+1])
    mesh = dirichlet_boundary |> GT.mesh
    #classify_mesh_faces!(Nface_to_tag,mesh,N,physical_names)
    Nface_to_tag[GT.faces(dirichlet_boundary)] .= 1
    ncells = length(cell_to_dofs)
    let d = N
        Dface_to_dfaces = d_to_Dface_to_dfaces[d+1]
        ctype_to_ldface_to_ldofs = d_to_ctype_to_ldface_to_dofs[d+1]
        for cell in 1:ncells
            ctype = cell_to_ctype[cell]
            Dface = cell_to_Dface[cell]
            ldof_to_dof = cell_to_dofs[cell]
            ldface_to_dface = Dface_to_dfaces[Dface]
            ldface_to_ldofs = ctype_to_ldface_to_ldofs[ctype]
            nldfaces = length(ldface_to_dface)
            dofs = cell_to_dofs[cell]
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
    data = cell_to_dofs.data
    data .= f.(data)
    ndiri = length(last(free_and_dirichlet_dofs))
    dirichlet_dof_location = ones(Int32,ndiri)
    (;cell_to_dofs, free_and_dirichlet_dofs,dirichlet_dof_location)
end

function generate_dof_ids_step_2(space,state,q::AbstractField)
    (;ndofs,cell_to_dofs,d_to_Dface_to_dfaces,
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
    ncells = length(cell_to_dofs)
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
                dofs = args.cell_to_dofs[cell]
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
    data = cell_to_dofs.data
    data .= f.(data)
    ndiri = length(last(free_and_dirichlet_dofs))
    dirichlet_dof_location = ones(Int32,ndiri)
    (;cell_to_dofs, free_and_dirichlet_dofs,dirichlet_dof_location)
end

function generate_dof_ids_step_2(space,state,dirichlet_boundary_all::PiecewiseDomain)
    (;ndofs,cell_to_dofs,d_to_Dface_to_dfaces,
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
            ncells = length(cell_to_dofs)
            Dface_to_dfaces = d_to_Dface_to_dfaces[d+1]
            ctype_to_ldface_to_ldofs = d_to_ctype_to_ldface_to_dofs[d+1]
            for cell in 1:ncells
                ctype = cell_to_ctype[cell]
                Dface = cell_to_Dface[cell]
                ldof_to_dof = cell_to_dofs[cell]
                ldface_to_dface = Dface_to_dfaces[Dface]
                ldface_to_ldofs = ctype_to_ldface_to_ldofs[ctype]
                nldfaces = length(ldface_to_dface)
                dofs = cell_to_dofs[cell]
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
    data = cell_to_dofs.data
    data .= f.(data)
    dirichlet_dof_location = dof_to_location[last(free_and_dirichlet_dofs)]
    (;cell_to_dofs, free_and_dirichlet_dofs, dirichlet_dof_location)
end

function generate_dof_ids_step_2(space,state,q_all::PiecewiseField)
    (;ndofs,cell_to_dofs,d_to_Dface_to_dfaces,
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
            ncells = length(cell_to_dofs)
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
                        dofs = args.cell_to_dofs[cell]
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
    data = cell_to_dofs.data
    data .= f.(data)
    ndiri = length(last(free_and_dirichlet_dofs))
    dirichlet_dof_location = ones(Int32,ndiri)
    (;cell_to_dofs, free_and_dirichlet_dofs,dirichlet_dof_location)
end

struct LastDof end

function last_dof()
    LastDof()
end

function generate_dof_ids_step_2(space,state,dirichlet_boundary_all::LastDof)
    (;ndofs,cell_to_dofs) = state
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
    data = cell_to_dofs.data
    data .= f.(data)
    ndiri = length(last(free_and_dirichlet_dofs))
    dirichlet_dof_location = ones(Int32,ndiri)
    (;cell_to_dofs, free_and_dirichlet_dofs,dirichlet_dof_location)
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

const LagrangeFaceDomain = Union{UnitNCube,UnitSimplex}

function lagrange_space(domain::LagrangeFaceDomain, order;
        space_type = default_space_type(domain),
        lib_to_user_nodes = :default,
        major = Val(:component),
        tensor_size = Val(:scalar),
    )


    D = num_dims(domain)
    order_per_dir = ntuple(d->order,Val(D))
    lagrange_face_space(;
               domain,
               order_per_dir,
               space_type,
               lib_to_user_nodes,
               major,
               tensor_size)
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
    )
    contents = (;
        domain,
        order_per_dir,
        space_type,
        lib_to_user_nodes,
        major,
        tensor_size)
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

function lib_to_user_nodes(fe::LagrangeFaceSpace)
    if val_parameter(fe.contents.lib_to_user_nodes) === :default
        nnodes = num_nodes(fe)
        Ti = int_type(options(fe))
        collect(Ti.(1:nnodes))
    else
        fe.contents.lib_to_user_nodes
    end
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
    )

    @assert conformity in (:default,:L2)

    lagrange_mesh_space(;
                        domain,
                        order,
                        conformity,
                        dirichlet_boundary,
                        space_type,
                        major,
                        tensor_size,
                        workspace,
                       ) |> setup_space

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
    LagrangeMeshSpace(contents)
end

struct LagrangeMeshSpace{A} <: AbstractSpace
    contents::A
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
    LagrangeMeshSpace(contents)
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
    vface_to_vrid = face_reference_id(V)
    domain = GT.domain(a)
    mesh = GT.mesh(domain)
    d = num_dims(domain)
    mrid_to_refface = reference_spaces(mesh,d)
    mface_to_mrid = face_reference_id(mesh,d)
    vface_to_mface = faces(domain)
    nvfaces = length(vface_to_vrid)
    vface_to_nodes = face_dofs(V)
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
        vrid = vface_to_vrid[vface]
        mface = vface_to_mface[vface]
        mrid = mface_to_mrid[mface]
        tab = vrid_mrid_tabulator[vrid][mrid]
        mnodes = mface_to_mnodes[mface]
        nlnodes,nlmnodes = size(tab)
        nodes = vface_to_nodes[vface]
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

function raviart_thomas_space(domain::AbstractFaceDomain,order)
    workspace = rt_setup((;domain,order))
    RaviartThomasFaceSpace(domain,order,workspace)
end

struct RaviartThomasFaceSpace{A,B,C} <: AbstractSpace
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
        domain,
        order,
        conformity,
        dirichlet_boundary,
        cell_to_ctype,
        ctype_to_reffe,
        workspace) |> setup_space
end

struct RaviartThomasMeshSpace{A,B,C,D,E,F} <: AbstractSpace
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
    CartesianProductSpace(spaces)
end

struct CartesianProductSpace{A} <: GT.AbstractSpace
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

function domain(a::CartesianProductSpace,field)
    GT.domain(GT.field(a,field))
end

function face_dofs(a::CartesianProductSpace)
    error("Not implemented, not needed in practice")
end

function face_dofs(a::CartesianProductSpace,field)
    face_dofs(a.spaces[field])
end

function free_dofs(a::CartesianProductSpace)
    nfields = GT.num_fields(a)
    map(1:nfields) do field
        free_dofs(a,field)
    end |> BRange
end

function free_dofs(a::CartesianProductSpace,field)
    f = GT.field(a,field)
    free_dofs(f)
end

function dirichlet_dofs(a::CartesianProductSpace)
    nfields = GT.num_fields(a)
    map(1:nfields) do field
        dirichlet_dofs(a,field)
    end |> BRange
end

function dirichlet_dofs(a::CartesianProductSpace,field)
    f = GT.field(a,field)
    dirichlet_dofs(f)
end

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

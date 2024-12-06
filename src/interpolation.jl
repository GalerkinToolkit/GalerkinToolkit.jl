
abstract type AbstractLagrangeFE <: AbstractMeshFace end

struct GenericLagrangeFE{A,B,C,D} <: AbstractLagrangeFE
    geometry::A
    order_per_dir::B
    space::Symbol
    lib_to_user_nodes::C
    major::Symbol
    shape::D
end

# TODO remove SCALAR_SHAPE and simply use Number?
# Scalar case should be a separate reffe?
struct ScalarShape end
const SCALAR_SHAPE = ScalarShape()

"""
"""
lagrangian_fe(args...) = GenericLagrangeFE(args...)

function lagrangian_fe(geometry,order;
        space = default_space(geometry),
        lib_to_user_nodes = int_type(geometry)[],
        major = :component,
        shape = SCALAR_SHAPE)
    D = num_dims(geometry)
    order_per_dir = repeat_per_dir(geometry,order)
    lagrangian_fe(
               geometry,
               order_per_dir,
               space,
               lib_to_user_nodes,
               major,
               shape)
end

order(fe::AbstractLagrangeFE) = maximum(order_per_dir(fe),init=0)

function lib_to_user_nodes(fe::AbstractLagrangeFE)
    if length(fe.lib_to_user_nodes) == 0
        nnodes = num_nodes(fe)
        Ti = int_type(geometry(fe))
        collect(Ti.(1:nnodes))
    else
        fe.lib_to_user_nodes
    end
end

function monomial_exponents(fe::AbstractLagrangeFE)
    monomial_exponents_from_space(fe.space,fe.order_per_dir,fe.geometry |> int_type)
end

function node_coordinates(fe::AbstractLagrangeFE)
    @assert fe |> geometry |> is_unitary
    mexps = monomial_exponents(fe)
    node_coordinates_from_monomials_exponents(mexps,fe.order_per_dir,fe.geometry |> real_type)
end

function tensor_basis(fe::AbstractLagrangeFE)
    Tv = fe |> geometry |> real_type
    if fe.shape == SCALAR_SHAPE
        return Tv(1)
    else
        cis = CartesianIndices(val_parameter(fe.shape))
        l = prod(val_parameter(fe.shape))
        init_tensor = SArray{Tuple{val_parameter(fe.shape)...},Tv}
        cis_flat = cis[:]
        return map(cis_flat) do ci
            init_tensor(ntuple(j->cis[j]==ci ? 1 : 0 ,Val(l)))
        end
    end
end

function primal_basis(fe::AbstractLagrangeFE)
    scalar_basis = map(e->(x-> prod(x.^e)),fe|>monomial_exponents)
    if fe.shape == SCALAR_SHAPE
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

function dual_basis(fe::AbstractLagrangeFE)
    node_coordinates_reffe = fe|>node_coordinates
    scalar_basis = map(x->(f->f(x)),node_coordinates_reffe)
    if fe.shape == SCALAR_SHAPE
        return scalar_basis
    else
        if fe.major === :component
            dual_nested = map(node_coordinates_reffe) do x
                map(tensor_basis(fe)) do e
                    f->inner(e,f(x))
                end
            end
        elseif fe.major === :node
            dual_nested = map(tensor_basis(fe)) do e
                map(node_coordinates_reffe) do x
                    f->inner(e,f(x))
                end
            end
        else
            error("Not Implemented")
        end
        return reduce(vcat,dual_nested)
    end
end

function node_dofs(fe::AbstractLagrangeFE)
    Tv = fe |> geometry |> int_type
    node_dofs(Tv,fe)
end

function node_dofs(::Type{Tv},fe) where Tv
    nnodes = num_nodes(fe)
    nodes =  1:nnodes
    if fe.shape == SCALAR_SHAPE
        return nodes
    else
        ndofs_per_node = prod(val_parameter(fe.shape))
        init_tensor = SArray{Tuple{val_parameter(fe.shape)...},Tv}
        node_to_dofs = map(nodes) do node
            t = ntuple(Val(ndofs_per_node)) do li
                if fe.major === :component
                    dof = (node-1)*ndofs_per_node + li
                elseif fe.major === :node
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

function dof_node(fe::AbstractLagrangeFE)
    Tv = fe |> geometry |> int_type
    dof_node(Tv,fe)
end

function dof_node(::Type{Tv},fe::AbstractLagrangeFE) where Tv
    ndofs = num_dofs(fe)
    if a.shape == SCALAR_SHAPE
        return 1:ndofs
    else
        dof_to_node = zeros(Tv,ndofs)
        node_to_dofs = node_dofs(fe)
        for (node,dofs) in enumerate(node_to_dofs)
            for dof in dofs
                dof_to_node[dof] = node
            end
        end
        return node_to_dofs
    end
end

function face_dofs_from_nodes(fe::AbstractLagrangeFE,face_to_nodes)
    if fe.shape == SCALAR_SHAPE
        return face_to_nodes
    else
        node_to_dofs = node_dofs(fe)
        lis = LinearIndices(val_parameter(fe.shape))
        face_to_dofs = map(face_to_nodes) do nodes
            if fe.major === :component
                nested = map(nodes) do node
                    dofs = node_to_dofs[node]
                    collect(dofs[:])
                end
            elseif fe.major === :node
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

function conforming(fe::AbstractLagrangeFE)
    nonconforming = (order(fe) == 0) || (is_n_cube(fe.geometry) && fe.space === :P)
    ! nonconforming
end

function interior_nodes(fe::AbstractLagrangeFE)
    if conforming(fe)
        interior_nodes_from_mesh_face(fe)
    else
        collect(Int,1:num_nodes(fe))
    end
end

function face_nodes(fe::AbstractLagrangeFE,d)
    if conforming(fe)
        face_nodes_from_mesh_face(fe,d)
    else
        D = num_dims(fe)
        if d == D
            [collect(Int32,1:GT.num_nodes(fe))]
        else
            [Int32[] for _ in 1:num_faces(boundary(fe.geometry),d)]
        end
    end
end

function face_interior_nodes(fe::AbstractLagrangeFE,d)
    if conforming(fe)
        face_interior_nodes_from_mesh_face(fe,d)
    else
        face_nodes(fe,d)
    end
end

function face_interior_node_permutations(fe::AbstractLagrangeFE,d)
    if conforming(fe)
        face_interior_node_permutations_from_mesh_face(fe,d)
    else
        D = num_dims(fe)
        if  d == D
            [[ collect(1:num_interior_nodes(fe)) ]]
        else
            [[ Int32[] ] for _ in 1:num_faces(boundary(fe.geometry),d)]
        end
    end
end

function face_dofs(a::AbstractLagrangeFE,d)
    face_to_nodes = face_nodes(a,d)
    face_dofs_from_nodes(a,face_to_nodes)
end

function face_own_dofs(a::AbstractLagrangeFE,d)
    face_to_nodes = face_interior_nodes(a,d)
    face_dofs_from_nodes(a,face_to_nodes)
end

function face_own_dof_permutations(fe::AbstractLagrangeFE,d)
    face_to_pindex_to_inodes = face_interior_node_permutations(fe,d)
    if fe.shape == SCALAR_SHAPE
        return face_to_pindex_to_inodes
    else
        Tv = fe |> geometry |> int_type
        node_to_dofs = node_dofs(fe)
        face_to_inodes = face_interior_nodes(fe,d)
        face_to_idofs = face_own_dofs(fe,d)
        nfaces = length(face_to_inodes)
        ndofs = num_dofs(fe)
        dof_to_idof = zeros(Tv,ndofs)
        lis = LinearIndices(val_parameter(fe.shape))
        map(1:nfaces) do face
            pindex_to_inodes = face_to_pindex_to_inodes[face]
            inode_to_node = face_to_inodes[face]
            idof_to_dof = face_to_idofs[face]
            fill!(dof_to_idof,Tv(0))
            dof_to_idof[idof_to_dof] = 1:length(idof_to_dof)
            map(pindex_to_inodes) do inodes
                if fe.major === :component
                    nested = map(inodes) do inode
                        node = inode_to_node[inode]
                        dofs = node_to_dofs[node]
                        map(dofs) do dof
                            idof = dof_to_idof[dof]
                            @assert idof != 0
                            idof
                        end
                    end
                elseif fe.major === :node
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

function num_dofs(a::AbstractLagrangeFE)
    nnodes = num_nodes(a)
    if a.shape == SCALAR_SHAPE
        nnodes
    else
        ndofs_per_node = prod(val_parameter(a.shape))
        nnodes*ndofs_per_node
    end
end

abstract type AbstractSpace <: GT.AbstractType end

function setup_space(space::AbstractSpace)
    if space.cache !== nothing
        return space
    end
    state = generate_dof_ids(space)
    face_dofs = state.cell_to_dofs
    free_and_dirichlet_dofs = state.free_and_dirichlet_dofs
    dirichlet_dof_location = state.dirichlet_dof_location
    cache = space_cache(;face_dofs,free_and_dirichlet_dofs,dirichlet_dof_location)
    replace_cache(space,cache)
end

function space_cache(;kwargs...)
    (;kwargs...)
end

Base.iterate(m::AbstractSpace) = iterate(components(m))
Base.iterate(m::AbstractSpace,state) = iterate(components(m),state)
Base.getindex(m::AbstractSpace,field::Integer) = component(m,field)
Base.length(m::AbstractSpace) = num_fields(m)

mesh(a::AbstractSpace) = GT.mesh(GT.domain(a))

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

# Single field by default

function num_fields(a::AbstractSpace)
    1
end

function components(a::AbstractSpace)
    (a,)
end

function component(a::AbstractSpace,field)
    @assert field == 1
    a
end

function domain(space::AbstractSpace,field)
    @assert field == 1
    space |> GT.domain
end

@enum FreeOrDirichlet FREE=1 DIRICHLET=2

function dofs(a,free_or_diri::FreeOrDirichlet)
    if free_or_diri == FREE
        free_dofs(a)
    else
        dirichlet_dofs(a)
    end
end

function values(a,free_or_diri::FreeOrDirichlet)
    if free_or_diri == FREE
        free_values(a)
    else
        dirichlet_values(a)
    end
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
    ctype_to_reference_fe = space |> GT.reference_fes
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
    physical_names = dirichlet_boundary |> GT.physical_names
    Nface_to_tag = zeros(Int32,d_to_ndfaces[N+1])
    mesh = dirichlet_boundary |> GT.mesh
    classify_mesh_faces!(Nface_to_tag,mesh,N,physical_names)
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
    physical_names = dirichlet_boundary |> GT.physical_names
    Nface_to_tag = zeros(Int32,d_to_ndfaces[N+1])
    mesh = dirichlet_boundary |> GT.mesh
    classify_mesh_faces!(Nface_to_tag,mesh,N,physical_names)
    ncells = length(cell_to_dofs)
    dim = 1
    ldof = :ldof
    sigma = GT.dual_basis(space,ldof)
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
            physical_names = dirichlet_boundary |> GT.physical_names
            Nface_to_tag = zeros(Int32,d_to_ndfaces[N+1])
            mesh = dirichlet_boundary |> GT.mesh
            classify_mesh_faces!(Nface_to_tag,mesh,N,physical_names)
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
            physical_names = dirichlet_boundary |> GT.physical_names
            Nface_to_tag = zeros(Int32,d_to_ndfaces[N+1])
            mesh = dirichlet_boundary |> GT.mesh
            classify_mesh_faces!(Nface_to_tag,mesh,N,physical_names)
            ncells = length(cell_to_dofs)
            dim = 1
            ldof = :ldof
            sigma = GT.dual_basis(space,ldof)
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
    ctype_to_reference_fe = reference_fes(space)
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
    ctype_to_reference_fe = reference_fes(space)
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

function face_dofs(space::AbstractSpace)
    if space.cache !== nothing
        return space.cache.face_dofs
    end
    state = generate_dof_ids(space)
    state.cell_to_dofs # TODO rename face_dofs ?
end

function face_dofs(space::AbstractSpace,field)
    @assert field == 1
    face_dofs(space)
end

function free_and_dirichlet_dofs(V::AbstractSpace)
    if V.cache !== nothing
        return V.cache.free_and_dirichlet_dofs
    end
    state = generate_dof_ids(V)
    state.free_and_dirichlet_dofs
end

function dirichlet_dof_location(V::AbstractSpace)
    if V.cache !== nothing
        return V.cache.dirichlet_dof_location
    end
    state = generate_dof_ids(V)
    state.dirichlet_dof_location
end

#function num_face_dofs(a::AbstractSpace,dim)
#    field = 1
#    num_face_dofs(a,dim,field)
#end
#
#function num_face_dofs(a::AbstractSpace,dim,field)
#    face_to_dofs = face_dofs(a)
#    index -> begin
#        dict = index.dict
#        field_sym = get!(dict,field,gensym("field"))
#        face_to_dofs_sym = get!(dict,face_to_dofs,gensym("face_to_dofs"))
#        dim_sym = get!(dict,dim,gensym("dim"))
#        face = index.face
#        field_per_dim = index.field_per_dim
#        @assert face !== nothing
#        @assert field_per_dim !== nothing
#        @term begin
#            f = length($face_to_dofs_sym[$face])
#            bool = ($field_per_dim[$dim_sym] != $field_sym) | ($face == 0)
#            ifelse(bool,0,f)
#        end
#    end
#end
#
#function dof_map(a::AbstractSpace,dim)
#    field = 1
#    dof_map(a,dim,field)
#end
#
#function dof_map(a::AbstractSpace,dim,field)
#    face_to_dofs = face_dofs(a)
#    T = eltype(eltype(face_to_dofs))
#    z = T(-1)
#    index -> begin
#        dict = index.dict
#        field_sym = get!(dict,field,gensym("field"))
#        z_sym = get!(dict,z,gensym("z"))
#        face_to_dofs_sym = get!(dict,face_to_dofs,gensym("face_to_dofs"))
#        dim_sym = get!(dict,dim,gensym("dim"))
#        face = index.face
#        dof = index.dof
#        dof_per_dim = index.dof_per_dim
#        field_per_dim = index.field_per_dim
#        @assert face !== nothing
#        @assert dof !== nothing
#        @assert dof_per_dim !== nothing
#        @assert field_per_dim !== nothing
#        @term begin
#            dof = $dof_per_dim[$dim_sym]
#            f = $face_to_dofs_sym[$face][dof]
#            bool = ($field_per_dim[$dim_sym] != $field_sym) | ($face == 0)
#            ifelse(bool,$z_sym,f)
#        end
#    end
#end

## TODO I guess this is not needed anymore
#function shape_function_mask(f,face_around_per_dim,face_around,dim)
#    @assert face_around !== nothing
#    x -> face_around_per_dim[dim] == face_around ? f(x) : zero(f(x))
#end
#
#function shape_function_mask(f,face_around_per_dim::Nothing,face_around,dim)
#    f
#end
#
#function shape_function_mask(f,face_around_per_dim::Nothing,face_around::Nothing,dim)
#    f
#end
#
#function shape_function_mask(f,face_around_per_dim,face_around::Nothing,dim)
#    f
#end

function primal_map(a::AbstractSpace)
    nothing
end

function dual_map(a::AbstractSpace)
    nothing
end

#function mask_function(f,bool)
#    x->bool*f(x)
#end
#
## TODO a better name
#function face_shape_function(rid_to_fs,face_to_rid,face,dof)
#    if dof == 0 || face == 0
#        return first(first(rid_to_fs))
#    end
#    fs = reference_value(rid_to_fs,face_to_rid,face)
#    fs[dof]
#end

function shape_functions(a::AbstractSpace,dof;reference=is_reference_domain(GT.domain(a)))
    @assert primal_map(a) === nothing
    face_to_refid = face_reference_id(a)
    refid_to_reffes = reference_fes(a)
    refid_to_funs = map(GT.shape_functions,refid_to_reffes)
    domain = GT.domain(a)
    qty = shape_function_quantity(refid_to_funs,domain;reference=true,face_reference_id=face_to_refid,dof)
    if reference
        return qty
    end
    D = num_dims(domain)
    ϕinv = GT.inverse_physical_map(mesh(domain),D)
    qty∘ϕinv
end

function form_argument(a::AbstractSpace,axis,field=1)
    @assert primal_map(a) === nothing
    face_to_refid = face_reference_id(a)
    refid_to_reffes = reference_fes(a)
    refid_to_funs = map(GT.shape_functions,refid_to_reffes)
    domain = GT.domain(a)
    qty = form_argument(axis,field,refid_to_funs,domain;reference=true,face_reference_id=face_to_refid)
    if is_reference_domain(GT.domain(a))
        return qty
    end
    D = num_dims(domain)
    ϕinv = GT.inverse_physical_map(mesh(domain),D)
    qty∘ϕinv
end

## TODO a better name
#function face_function_free_and_dirichlet(
#    rid_to_fs,
#    face_to_rid,
#    face_to_dofs,
#    free_dofs_to_value,
#    diri_dofs_to_value,
#    face)
#    fs = reference_value(rid_to_fs,face_to_rid,face)
#    dofs = face_to_dofs[face]
#    x -> begin
#        ndofs = length(dofs)
#        sum(1:ndofs) do idof
#            dof = dofs[idof]
#            coeff =  dof > 0 ? free_dofs_to_value[dof] : diri_dofs_to_value[-dof]
#            f = fs[idof]
#            f(x)*coeff
#        end
#    end
#end

function free_or_dirichlet_value(free_vals,diri_vals,dof)
    if dof < 0
        diri_vals[-dof]
    else
        free_vals[dof]
    end
end

function discrete_field_qty(a::AbstractSpace,free_vals,diri_vals)
    ldof = gensym("fe-dof")
    s = shape_functions(a,ldof;reference=true)
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
    #face_to_refid = face_reference_id(a)
    #refid_to_reffes = reference_fes(a)
    #refid_to_funs = map(GT.shape_functions,refid_to_reffes)
    #domain = GT.domain(a)
    #f_prototype = first(first(refid_to_funs))
    #z = zero(eltype(free_vals))
    #prototype = x-> z*f_prototype(x) + z*f_prototype(x)
    #face_to_dofs = GT.face_dofs(a)
    #qty = GT.quantity(prototype,domain) do index
    #    dict = index.dict
    #    face_to_refid_sym = get!(dict,face_to_refid,gensym("face_to_refid"))
    #    refid_to_funs_sym = get!(dict,refid_to_funs,gensym("refid_to_funs"))
    #    face_to_dofs_sym = get!(dict,face_to_dofs,gensym("face_to_dofs"))
    #    free_vals_sym = get!(dict,free_vals,gensym("free_vals"))
    #    diri_vals_sym = get!(dict,diri_vals,gensym("diri_vals"))
    #    face = index.face
    #    @assert face !== nothing
    #    @term begin
    #        face_function_free_and_dirichlet(
    #                                         $refid_to_funs_sym,
    #                                         $face_to_refid_sym,
    #                                         $face_to_dofs_sym,
    #                                         $free_vals_sym,
    #                                         $diri_vals_sym,
    #                                         $face)
    #    end
    #end
    if is_reference_domain(GT.domain(a))
        return qty
    end
    ϕinv = GT.inverse_physical_map(mesh(domain(a)),D)
    qty∘ϕinv
end

function dual_basis(a::AbstractSpace,dof)
    @assert dual_map(a) === nothing
    face_to_refid = face_reference_id(a)
    refid_to_reffes = reference_fes(a)
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
        expr_term([D],expr,prototype,index)
    end
    if is_reference_domain(GT.domain(a))
        return s
    end
    ϕ = GT.physical_map(mesh(domain),D)
    call(dual_compose,s,ϕ)
    ##call(dual_compose,s,ϕ)
    #ϕ_proto = GT.prototype(ϕ)
    #s_proto = GT.prototype(s)
    #proto = dual_compose(s_proto,ϕ_proto)
    #quantity(proto,Ω) do index
    #    ϕ_expr = term(ϕ)(index)
    #    s_expr = term(s)(index)
    #    @term begin
    #        dual_compose($(s_expr),$(ϕ_expr))
    #    end
    #end
    ####u -> qty(compose(u,ϕ))
end

function dual_compose(s,g)
    f -> s(f∘g)
end

function get_symbol!(index,::typeof(dual_compose),name="";prefix=index.data.prefix)
    :(GalerkinToolkit.dual_compose)
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

function components(a::CartesianProductSpace)
    a.spaces
end

function component(u::CartesianProductSpace,field::Integer)
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
    GT.domain(GT.component(a,field))
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
    f = component(a,field)
    free_dofs(f)
end

function dirichlet_dofs(a::CartesianProductSpace)
    nfields = GT.num_fields(a)
    map(1:nfields) do field
        dirichlet_dofs(a,field)
    end |> BRange
end

function dirichlet_dofs(a::CartesianProductSpace,field)
    f = component(a,field)
    dirichlet_dofs(f)
end

function form_argument(a::CartesianProductSpace,axis)
    fields = ntuple(identity,GT.num_fields(a))
    map(fields) do field
        GT.form_argument(a,dim,field)
    end
end

#function shape_functions(a::CartesianProductSpace,dim)
#    fields = ntuple(identity,GT.num_fields(a))
#    map(fields) do field
#        GT.shape_functions(a,dim,field)
#    end
#end
#
#function shape_functions(a::CartesianProductSpace,dim,field)
#    f = component(a,field)
#    qty = shape_functions(f,dim)
#    t = GT.term(qty)
#    GT.quantity(GT.prototype(qty),GT.domain(qty)) do index
#        # TODO field_per_dim[dim] has possibly already a numeric value
#        field_per_dim = index.field_per_dim
#        dof_per_dim = index.dof_per_dim
#        dof = @term begin
#            coeff = $field == $(field_per_dim[dim])
#            coeff*$(dof_per_dim[dim])
#        end
#        dof_per_dim2 = Base.setindex(dof_per_dim,dof,dim)
#        index2 = replace_dof_per_dim(index,dof_per_dim2)
#        @assert field_per_dim !== nothing
#        @term begin
#            bool = $field == $(field_per_dim[dim])
#            mask_function($(t(index2)),bool)
#        end
#    end
#end

function discrete_field_qty(a::CartesianProductSpace,free_vals,diri_vals)
    nothing
end

function dual_basis(a::CartesianProductSpace)
    error("Not implemented yet. Not needed in practice.")
end

function dual_basis(a::CartesianProductSpace,field)
    error("Not implemented yet. Not needed in practice.")
end

function discrete_field(space,free_values,dirichlet_values)
    qty = discrete_field_qty(space,free_values,dirichlet_values)
    mesh = space |> GT.mesh
    DiscreteField(mesh,space,free_values,dirichlet_values,qty)
end

struct DiscreteField{A,B,C,D,E} <: GT.AbstractField
    mesh::A
    space::B
    free_values::C
    dirichlet_values::D
    quantity::E
end

Base.iterate(m::DiscreteField) = iterate(components(m))
Base.iterate(m::DiscreteField,state) = iterate(components(m),state)
Base.getindex(m::DiscreteField,field::Integer) = component(m,field)
Base.length(m::DiscreteField) = num_fields(m)

function allocate_values(::Type{T},dofs) where T
    n = length(dofs)
    Vector{T}(undef,n)
end

function allocate_values(::Type{T},dofs::BRange) where T
    map(dofs.blocks) do dofs_i
        allocate_values(T,dofs_i)
    end |> BVector
end

function rand_values(::Type{T},dofs) where T
    n = length(dofs)
    rand(T,n)
end

function rand_values(::Type{T},dofs::BRange) where T
    map(dofs.blocks) do dofs_i
        rand_values(T,dofs_i)
    end |> BVector
end

function constant_values(f,dofs)
    n = length(dofs)
    Fill(f,n)
end

function constant_values(f,dofs::BRange)
    map(dofs.blocks) do dofs_i
        constant_values(f,dofs_i)
    end |> BVector
end

function undef_field(::Type{T},space::AbstractSpace) where T
    free_values = allocate_values(T,GT.free_dofs(space))
    dirichlet_values = allocate_values(T,GT.dirichlet_dofs(space))
    discrete_field(space,free_values,dirichlet_values)
end

function rand_field(::Type{T},space::AbstractSpace) where T
    free_values = rand_values(T,GT.free_dofs(space))
    dirichlet_values = zeros(T,GT.dirichlet_dofs(space))
    discrete_field(space,free_values,dirichlet_values)
end

function zero_field(::Type{T},space::AbstractSpace) where T
    fill_field(zero(T),space)
end

function fill_field!(u::DiscreteField,z)
    fill!(free_values(u),z)
    fill!(dirichlet_values(u),z)
    u
end

function fill_field(z,space::AbstractSpace)
    u = undef_field(typeof(z),space)
    fill_field!(u,z)
end

function undef_dirichlet_field(::Type{T},space::AbstractSpace) where T
    free_values = constant_values(zero(T),GT.free_dofs(space))
    dirichlet_values = allocate_values(T,GT.dirichlet_dofs(space))
    discrete_field(space,free_values,dirichlet_values)
end

function zero_dirichlet_field(::Type{T},space::AbstractSpace) where T
    uh = undef_dirichlet_field(T,space)
    fill!(dirichlet_values(uh),zero(T))
    uh
end

function dirichlet_field(::Type{T},space::AbstractSpace) where T
    zero_dirichlet_field(T,space)
end

function num_fields(u::DiscreteField)
    num_fields(GT.space(u))
end

function components(u)
    ntuple(GT.num_fields(u)) do field
        component(u,field)
    end
end

function component(u::DiscreteField,field::Integer)
    space = GT.space(u)
    space_field = component(space,field)
    free_values_field = view(free_values(u),Block(field))
    dirichlet_values_field = view(dirichlet_values(u),Block(field))
    discrete_field(space_field,free_values_field,dirichlet_values_field)
end

function space(x::DiscreteField)
    x.space
end

function free_values(x::DiscreteField)
    x.free_values
end

function dirichlet_values(x::DiscreteField)
    x.dirichlet_values
end

function domain(u::DiscreteField)
    GT.domain(GT.space(u))
end

#function prototype(u::DiscreteField)
#    GT.prototype(u.quantity)
#    #space = GT.space(u)
#    #dim = 2
#    #basis = GT.shape_functions(space,dim)
#    #f = GT.prototype(basis)
#    #T = eltype(GT.free_values(u))
#    #z = zero(T)
#    ## TODO use also dirichlet values type?
#    #x -> f(x)*z + f(x)*z
#end

#function term(u::DiscreteField,i...)
#    GT.term(u.quantity,i...)
#    #space = GT.space(u)
#    #free_vals = GT.free_values(u)
#    #diri_vals = GT.dirichlet_values(u)
#    #space = GT.space(u)
#    #free_vals = GT.free_values(u)
#    #diri_vals = GT.dirichlet_values(u)
#    #dim = 2
#    #dof_map = GT.dof_map(space,dim)
#    #num_face_dofs = GT.num_face_dofs(space,dim)
#    #basis = GT.shape_functions(space,dim)
#    #basis_term = GT.term(basis)
#    #index -> begin
#    #    x -> begin
#    #        index2 = GT.index(
#    #                          face=index.face,
#    #                          local_face=index.local_face,
#    #                          face_around=index.face_around)
#    #        ndofs = num_face_dofs(index2)
#    #        sum(1:ndofs) do dof
#    #            dof_per_dim = (nothing,dof)
#    #            index3 = replace_dof_per_dim(index2,dof_per_dim)
#    #            dof = dof_map(index3)
#    #            f = basis_term(index3)
#    #            coeff =  dof > 0 ? free_vals[dof] : diri_vals[-dof]
#    #            f(x)*coeff
#    #        end
#    #    end
#    #end
#end

# TODO
# With an extra param we can allow interpolation
# in another domain
# what about coBoundary?
# Or should we detect the domain of f?
#GT.interpolate!(uref,v;domain=GT.domain(uref))
#GT.interpolate!(uref,v;domain_map=)
#GT.interpolate!(uref,v;domain_glue=)
#GT.interpolate!(u,v;domain=Ω)
# We could automatic discover the domain of f
# but we can ask for an addition kwarg domain
# to confirm that we want the domain of f
# and we can also provide another one witht he glue?
# Solution: Two options: u is defined either on the domain of the space
# or on the Dirichlet boundary

function interpolate!(f,u::DiscreteField)
    interpolate!(f,u,nothing)
end

function interpolate!(f,u::DiscreteField,field)
    interpolate!(f,u,nothing,field)
end

function interpolate_free!(f,u::DiscreteField)
    interpolate!(f,u,FREE)
end

function interpolate_free!(f,u::DiscreteField,field)
    interpolate!(f,u,FREE,field)
end

function interpolate_dirichlet!(f,u::DiscreteField)
    # TODO for dirichlet we perhaps want to allow integrating on boundaries?
    interpolate!(f,u,DIRICHLET)
end

function interpolate_dirichlet!(f,u::DiscreteField,field)
    # TODO for dirichlet we perhaps want to allow integrating on boundaries?
    interpolate!(f,u,DIRICHLET,field)
end

function interpolate!(f,u::DiscreteField,free_or_diri::Union{Nothing,FreeOrDirichlet})
    interpolate_impl!(f,u,free_or_diri)
end

function interpolate!(f,u::DiscreteField,free_or_diri::Union{Nothing,FreeOrDirichlet},field)
    ui = component(u,field)
    interpolate_impl!(f,ui,free_or_diri)
    u
end

function interpolate_impl!(f,u,free_or_diri;location=1)
    free_vals = GT.free_values(u)
    diri_vals = GT.dirichlet_values(u)
    space = GT.space(u)
    dof = gensym("fe-dof")
    sigma = GT.dual_basis(space,dof)
    face_to_dofs = GT.face_dofs(space)
    vals = sigma(f)
    domain = GT.domain(space)
    dirichlet_dof_location = GT.dirichlet_dof_location(space)
    index = generate_index(domain)
    sface_to_face = get_symbol!(index,faces(domain),"sface_to_face")
    t = term(vals,index)
    expr_qty = t |> expression |> simplify
    D = num_dims(domain)
    face = face_index(index,D)
    s_qty = GT.topological_sort(expr_qty,(face,dof))
    expr = quote
        (args,storage) -> begin
            (;face_to_dofs,free_vals,diri_vals,nfaces,dirichlet_dof_location,location,free_or_diri) = args
            $(unpack_index_storage(index,:storage))
            $(s_qty[1])
            for sface in 1:nfaces
                $face = $sface_to_face[sface]
                $(s_qty[2])
                dofs = face_to_dofs[$face]
                ndofs = length(dofs)
                for $dof in 1:ndofs
                    v = $(s_qty[3])
                    gdof = dofs[$dof]
                    if gdof > 0
                        if free_or_diri != DIRICHLET
                            free_vals[gdof] = v
                        end
                    else
                        diri_dof = -gdof
                        if free_or_diri != FREE && dirichlet_dof_location[diri_dof] == location
                            diri_vals[diri_dof] = v
                        end
                    end
                end
            end
        end
    end
    loop! = eval(expr)
    storage = GT.index_storage(index)
    nfaces = GT.num_faces(domain)
    args = (;face_to_dofs,free_vals,diri_vals,nfaces,dirichlet_dof_location,location,free_or_diri)
    invokelatest(loop!,args,storage)
    u
end

function interpolate_impl!(f::PiecewiseField,u,free_or_diri)
    for (location,field) in f.fields |> enumerate
        interpolate_impl!(field,u,free_or_diri;location)
    end
    u
end

# This space is not suitable for assembling problems in general, as there might be gaps in the dofs

function iso_parametric_space(dom::ReferenceDomain;
        dirichlet_boundary=nothing,
        cache=nothing)
    IsoParametricSpace(dom,dirichlet_boundary,cache) |> setup_space
end

struct IsoParametricSpace{A,B,C} <: AbstractSpace
    domain::A
    dirichlet_boundary::B
    cache::C
end

#function setup_space(V::IsoParametricSpace)
#    if V.cache !== nothing
#        return V
#    end
#    face_dofs = GT.face_dofs(V)
#    cache = space_cache(;face_dofs)
#    replace_cache(V,cache)
#end

function replace_cache(V::IsoParametricSpace,cache)
    GT.iso_parametric_space(
        V.domain;
        V.dirichlet_boundary,
        cache
    )
end

function free_and_dirichlet_nodes(V::IsoParametricSpace)
    free_and_dirichlet_nodes(V,V.dirichlet_boundary)
end

function free_and_dirichlet_nodes(V::IsoParametricSpace,dirichlet_boundary::Nothing)
    n = V |> GT.domain |> GT.mesh |> GT.num_nodes
    GT.partition_from_mask(fill(true,n))
end

function free_and_dirichlet_nodes(V::IsoParametricSpace,dirichlet_boundary::AbstractDomain)
    mesh = V |> GT.domain |> GT.mesh
    n = mesh |> GT.num_nodes
    node_to_tag = fill(Int32(0),n)
    tag_to_name = dirichlet_boundary |> GT.physical_names
    GT.classify_mesh_nodes!(node_to_tag,mesh,tag_to_name)
    GT.partition_from_mask(i->i==0,node_to_tag)
end

function free_and_dirichlet_dofs(V::IsoParametricSpace)
    free_and_dirichlet_nodes(V)
end

function dirichlet_dof_location(V::IsoParametricSpace)
    ndiri = length(last(free_and_dirichlet_nodes(V)))
    ones(Int32,ndiri)
end

function face_dofs(V::IsoParametricSpace)
    if V.cache !== nothing
        return V.cache.face_dofs
    end
    face_dofs(V,V.dirichlet_boundary)
end

function face_dofs(V::IsoParametricSpace,dirichlet_boundary::Nothing)
    domain = GT.domain(V)
    d = GT.face_dim(domain)
    faces = GT.faces(domain)
    face_nodes(GT.mesh(domain),d)
    face_to_nodes = view(face_nodes(GT.mesh(domain),d),faces)
    face_to_nodes
end

function face_dofs(V::IsoParametricSpace,dirichlet_boundary::AbstractDomain)
    face_to_nodes = JaggedArray(GT.face_dofs(V,nothing))
    free_and_dirichlet_nodes = GT.free_and_dirichlet_nodes(V)
    node_to_free_node = GT.permutation(free_and_dirichlet_nodes)
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
    domain = GT.domain(V)
    d = GT.face_dim(domain)
    # TODO LagrangianMesh face needs to be a AbstractElement
    reference_faces(GT.mesh(domain),d)
end

function face_reference_id(V::IsoParametricSpace)
    domain = GT.domain(V)
    faces = GT.faces(domain)
    d = GT.face_dim(domain)
    face_refid = face_reference_id(GT.mesh(domain),d)
    view(face_refid,faces)
end

function domain(a::IsoParametricSpace)
    a.domain
end

# TODO rename kwarg space?
function lagrange_space(domain,order;
    conformity = :default,
    dirichlet_boundary=nothing,
    space=nothing,
    major=:component,
    shape=SCALAR_SHAPE,
    cache=nothing,
    )

    @assert conformity in (:default,:L2)

    LagrangeSpace(
                  domain,
                  order,
                  conformity,
                  dirichlet_boundary,
                  space,
                  major,
                  shape,
                  cache) |> setup_space

end

struct LagrangeSpace{A,B,C,F,G,H,I} <: AbstractSpace
    domain::A
    order::B
    conformity::Symbol
    dirichlet_boundary::C
    space::F # TODO rename this one?
    major::G
    shape::H
    cache::I
end

#function setup_space(space::LagrangeSpace)
#    if space.cache !== nothing
#        return space
#    end
#    face_dofs = GT.face_dofs(space)
#    cache = space_cache(;face_dofs)
#    replace_cache(space,cache)
#end

function replace_cache(space::LagrangeSpace,cache)
    GT.lagrange_space(
        space.domain,
        space.order;
        space.conformity,
        space.dirichlet_boundary,
        space.space,
        space.major,
        space.shape,
        cache
    )
end

function face_nodes(a::LagrangeSpace)
    major = :component
    shape = SCALAR_SHAPE
    dirichlet_boundary = nothing
    cache = nothing
    V = LagrangeSpace(
                      a.domain,
                      a.order,
                      a.conformity,
                      dirichlet_boundary,
                      a.space,
                      major,
                      shape,
                      cache,
                     )
    face_dofs(V)
end

function node_coordinates(a::LagrangeSpace)
    major = :component
    shape = SCALAR_SHAPE
    dirichlet_boundary = nothing
    cache = nothing
    V = LagrangeSpace(
                      a.domain,
                      a.order,
                      a.conformity,
                      dirichlet_boundary,
                      a.space,
                      major,
                      shape,
                      cache,
                     )
    vrid_to_reffe = reference_fes(V)
    vface_to_vrid = face_reference_id(V)
    domain = a.domain
    mesh = GT.mesh(domain)
    d = num_dims(domain)
    mrid_to_refface = reference_faces(mesh,d)
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
#function node_dofs(space::LagrangeSpace)
#end
#

function free_dof_node(space::LagrangeSpace)
    free_and_dirichlet_dof_node(space)[1]
end

function dirichlet_dof_node(space::LagrangeSpace)
    free_and_dirichlet_dof_node(space)[2]
end

function free_and_dirichlet_dof_node(space::LagrangeSpace)
    T = Float64
    D = num_ambient_dims(space.domain)
    major = :component
    shape = SCALAR_SHAPE
    dirichlet_boundary = nothing
    cache = nothing
    V = LagrangeSpace(
                      space.domain,
                      space.order,
                      space.conformity,
                      dirichlet_boundary,
                      space.space,
                      major,
                      shape,
                      cache,
                     )
    face_to_nodes = GT.face_dofs(V)
    face_to_dofs = GT.face_dofs(space)
    nfree = length(free_dofs(space))
    ndiri = length(dirichlet_dofs(space))
    free_dof_to_node = zeros(Int32,nfree)
    diri_dof_to_node = zeros(Int32,ndiri)
    rid_to_lnode_to_ldofs = map(node_dofs,reference_fes(space))
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

conformity(space::LagrangeSpace) = space.conformity

function domain(space::LagrangeSpace)
    space.domain
end

function dirichlet_boundary(space::LagrangeSpace)
    space.dirichlet_boundary
end

function face_reference_id(space::LagrangeSpace)
    domain = space |> GT.domain
    mesh = domain |> GT.mesh
    cell_to_Dface = domain |> GT.faces
    D = domain |> GT.num_dims
    Dface_to_ctype = GT.face_reference_id(mesh,D)
    cell_to_ctype = Dface_to_ctype[cell_to_Dface]
    cell_to_ctype
end

function reference_fes(space::LagrangeSpace)
    domain = space |> GT.domain
    mesh = domain |> GT.mesh
    D = domain |> GT.num_dims
    order = space.order
    major = space.major
    shape = space.shape
    space3 = space.space # TODO Ugly
    ctype_to_refface = GT.reference_faces(mesh,D)
    ctype_to_geometry = map(GT.geometry,ctype_to_refface)
    ctype_to_reffe = map(ctype_to_geometry) do geometry
        space2 = space3 === nothing ? default_space(geometry) : space3
        lagrangian_fe(geometry,order;space=space2,major,shape)
    end
    ctype_to_reffe
end


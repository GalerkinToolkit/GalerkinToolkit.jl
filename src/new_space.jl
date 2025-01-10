

"""
    abstract type AbstractSpace <: AbstractType end

# Basic queries

[`domain`](@ref)
[`num_dofs`](@ref)
[`num_nodes`](@ref)
[`face_dofs`](@ref)
[`face_nodes`](@ref)
[`face_reference_id`](@ref)
[`reference_spaces`](@ref)
[`interior_nodes`](@ref)
[`interior_nodes_permutations`](@ref)
[`geometry_own_dofs`](@ref)
[`geometry_own_dofs_permutations`](@ref)
[`geometry_interior_nodes`](@ref)
[`geometry_interior_nodes_permutations`](@ref)
[`geometry_nodes`](@ref)
[`geometry_nodes_permutations`](@ref)

# Basic constructors

[`lagrange_space`](@ref)
[`raviart_thomas_space`](@ref)

"""
abstract type AbstractSpace <: AbstractType end

abstract type AbstractFaceSpace <: AbstractSpace end

"""
"""
function shape_functions(fe::AbstractFaceSpace)
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
function tabulator(fe::AbstractFaceSpace)
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

function lagrange_space(domain::LagrangeFaceDomain;
        order = 1,
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

function default_space_type(geom::LagrangeFaceDomain)
    if is_simplex(geom)
        :P
    elseif is_n_cube(geom)
        :Q
    else
        error("Not implemented")
    end
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
space_type(fe::LagrangeFaceSpace) = fe.contents.space_type
lib_to_user_nodes(fe::LagrangeFaceSpace) = fe.contents.lib_to_user_nodes
major(fe::LagrangeFaceSpace) = val_parameter(fe.contents.major)
tensor_size(fe::LagrangeFaceSpace) = val_parameter(fe.contents.tensor_size)

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
        error("Case not implemented (yet)")
    end
end

num_nodes(fe::LagrangeFaceSpace) = length(monomial_exponents(fe))

function node_coordinates(a::LagrangeFaceSpace)
    if order(a) == 0 && num_dims(a) != 0
        a_linear = lagrange_space(domain(a))
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

num_dofs(fe::LagrangeFaceSpace) = length(primal_basis(fe))

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

function interior_nodes(fe::LagrangeFaceSpace)
    nnodes = num_nodes(fe)
    D = num_dims(fe)
    if D == 0
        return collect(1:nnodes)
    else
        mesh = mesh_from_space(fe)
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

function node_quadrature(fe::LagrangeFaceSpace)
    coordinates = node_coordinates(fe)
    Tv = real_type(options(fe))
    nnodes = length(coordinates)
    weights = fill(Tv(1/nnodes),nnodes)
    domain = GT.domain(fe)
    face_quadrature(;domain,coordinates,weights)
end

function mesh_from_space(refface::LagrangeFaceSpace)
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
        map(domain->lagrange_space(domain;order=order_inter),rid_to_domain)
    end
    face_nodes_tuple = map(dims) do d
        domain = GT.domain(mesh,d)
        rid_to_refqua = map(reference_domains(domain)) do refdom
            reffe = lagrange_space(refdom;order=order_inter)
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



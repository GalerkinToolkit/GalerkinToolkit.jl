

"""
    abstract type AbstractType end

Parent of all types defined in GalerkinToolkit.
"""
abstract type AbstractType end

function Base.show(io::IO,data::GT.AbstractType)
    print(io,"GalerkinToolkit.$(nameof(typeof(data)))(…)")
end

"""
    push(a,ai)

Like `push!`, but creates a new object to store the result. 
This function is used to push to immutable collections such as tuples.
"""
function push end

function push(a::AbstractVector,x)
    b = copy(a)
    push!(b,x)
    b
end

function push(a::Tuple,x)
    (a...,x)
end

"""
    val_parameter(a)

For `a::Val{A}` it returns `A`. Otherwise, it returns `a`.
"""
val_parameter(a) = a
val_parameter(::Val{a}) where a = a

"""
    options(;kwargs...) -> Options

Create an object representing the default options for the current simulation.
This object can be used as an optional argument in several object constructors in GalerkinToolkit,
such as the mesh constructors `cartesian_mesh` and `mesh_from_gmsh`.
In this case, the computations using the generated mesh, will use the given options by default.
"""
function options(;
    reference_int_type=Int16,
    int_type=Int32,
    global_int_type=Int,
    real_type=Float64,
    )
    contents = (;reference_int_type,int_type,global_int_type,real_type)
    Options(contents)
end

options(object::AbstractType) = object.options

"""
    struct Options{...} <: AbstractType

Type of the objects returned by function `options`.
All properties and type parameters are private.

# Basic queries

- [`reference_int_type`](@ref)
- [`int_type`](@ref)
- [`global_int_type`](@ref)
- [`real_type`](@ref)
"""
struct Options{A} <: AbstractType
    contents::A
end

"""
    reference_int_type(options::Options)

Return the type of the integers used to enumerate reference quantities.
"""
reference_int_type(options::Options) = options.contents.reference_int_type

"""
    int_type(options::Options)

Return the default integer type used in the computation except for reference and global quantities.
"""
int_type(options::Options) = options.contents.int_type

"""
    global_int_type(options::Options)

Return the type of the integers used to enumerate global quantities.
"""
global_int_type(options::Options) = options.contents.global_int_type

"""
    real_type(options::Options)

Return the default real type used in the computation.
"""
real_type(options::Options) = options.contents.real_type

abstract type AbstractDomain <: AbstractType end

domain(a::AbstractDomain) = a

"""
    abstract type AbstractFaceDomain <: AbstractType end

Abstract type representing the geometry of a single mesh face, typically one of the reference faces.

# Basic queries

- [`num_dims`](@ref)
- [`is_axis_aligned`](@ref)
- [`is_simplex`](@ref)
- [`is_n_cube`](@ref)
- [`is_unit_n_cube`](@ref)
- [`is_unit_simplex`](@ref)
- [`is_unitary`](@ref)
- [`bounding_box`](@ref)
- [`boundary`](@ref)
- [`vertex_permutations`](@ref)

# Basic constructors

- [`unit_simplex`](@ref)
- [`unit_n_cube`](@ref)

"""
abstract type AbstractFaceDomain <: AbstractDomain end

function is_unit_n_cube(geo::AbstractFaceDomain)
    is_n_cube(geo) && is_unitary(geo)
end

"""
"""
function is_unit_simplex(geo::AbstractFaceDomain)
    is_simplex(geo) && is_unitary(geo)
end

"""
    unit_n_cube(d)
    unit_n_cube(Val(d))

Return an object representing a unit `d`-cube.
"""
function unit_n_cube(d;options=GT.options())
    D = val_parameter(d)
    num_dims = Val{D}()
    UnitNCube(num_dims,options)
end

struct UnitNCube{D,A} <: AbstractFaceDomain
    num_dims::Val{D}
    options::A
end

num_dims(geo::UnitNCube{D}) where D = D
is_n_cube(geo::UnitNCube) = true
is_simplex(geo::UnitNCube) = false
is_simplex(geo::UnitNCube{0}) = true
is_simplex(geo::UnitNCube{1}) = true
is_axis_aligned(geo::UnitNCube) = true
is_unitary(geo::UnitNCube) = true

function bounding_box(geo::UnitNCube)
    D = num_dims(geo)
    Tv = real_type(options(geo))
    p0 = ntuple(i->zero(real_type(options(geo))),Val(D)) |> SVector{D,Tv}
    p1 = ntuple(i->one(real_type(options(geo))),Val(D)) |> SVector{D,Tv}
    (p0,p1)
end

"""
    unit_simplex(d)
    unit_simplex(Val(d))

Return an object representing a unit simplex of dimension `d`.
"""
function unit_simplex(d;options=GT.options())
    D = val_parameter(d)
    num_dims = Val{D}()
    UnitSimplex(num_dims,options)
end

struct UnitSimplex{D,A} <: AbstractFaceDomain
    num_dims::Val{D}
    options::A
end

num_dims(geo::UnitSimplex{D}) where D = D
is_simplex(geo::UnitSimplex) = true
is_n_cube(geo::UnitSimplex) = false
is_n_cube(geo::UnitSimplex{0}) = true
is_n_cube(geo::UnitSimplex{1}) = true
is_axis_aligned(geo::UnitSimplex) = true
is_unitary(geo::UnitSimplex) = true
bounding_box(geo::UnitSimplex) = bounding_box(UnitNCube(geo.num_dims,geo.options))

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

function node_quadrature(fe::AbstractFaceSpace)
    coordinates = node_coordinates(fe)
    Tv = real_type(options(fe))
    nnodes = length(coordinates)
    weights = fill(Tv(1/nnodes),nnodes)
end

function lagrange_space(domain::AbstractFaceDomain;
        order = 1,
        space_type = default_space_type(domain),
        lib_to_user_nodes = :default,
        major = :component,
        tensor_size = :scalar,
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

function default_space_type(geom::AbstractFaceDomain)
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
major(fe::LagrangeFaceSpace) = fe.contents.major
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
    if tensor_size(fe) === :scalar
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

abstract type AbstractQuadrature <: AbstractType end

"""
    abstract type AbstractFaceQuadrature

# Basic queries

- [`domain`](@ref)
- [`coordinates`](@ref)
- [`weights`](@ref)

# Basic constructors

- [`quadrature`](@ref)
- [`duffy_quadrature`](@ref)
- [`tensor_product_quadrature`](@ref)

# Supertype hierarchy

    AbstractQuadrature <: GT.AbstractType
"""
abstract type AbstractFaceQuadrature <: AbstractQuadrature end

function face_quadrature(;domain,coordinates,weights)
    contents = (;domain,coordinates,weights)
    FaceQuadrature(contents)
end

struct FaceQuadrature{A} <: AbstractFaceQuadrature
    contents::A
end

domain(a::FaceQuadrature) = a.contents.domain
coordinates(a::FaceQuadrature) = a.contents.coordinates
weights(a::FaceQuadrature) = a.contents.weights

"""
"""
function quadrature(geo::AbstractFaceDomain,degree)
    if is_n_cube(geo) && is_axis_aligned(geo)
        D = num_dims(geo)
        tensor_product_quadrature(geo,degree)
    elseif is_unit_simplex(geo)
        duffy_quadrature(geo,degree)
    else
        error("Not implemented")
    end
end

function tensor_product_quadrature(geo::AbstractFaceDomain,degree::Integer)
    D = num_dims(geo)
    degree_per_dir = ntuple(d->degree,Val(D))
    tensor_product_quadrature(geo,degree_per_dir)
end

"""
"""
function tensor_product_quadrature(geo::AbstractFaceDomain,degree_per_dir)
    @assert is_n_cube(geo)
    @assert is_axis_aligned(geo)
    my_bounding_box = bounding_box(geo)
    D = num_dims(geo)
    limits_per_dir = ntuple(i->(my_bounding_box[1][i],my_bounding_box[2][i]),Val(D))
    n_per_dir = map(d->ceil(Int,(d+1)/2),degree_per_dir)
    function quadrature_1d(n,limits)
        x,w = FastGaussQuadrature.gausslegendre(n)
        a,b = limits
        x .= (0.5*(b-a)) .*x .+ (0.5*(b+a))
        w .*= 0.5*(b-a)
        x,w
    end
    quad_per_dir = map(quadrature_1d,n_per_dir,limits_per_dir)
    coords_per_dir = map(first,quad_per_dir)
    weights_per_dir = map(last,quad_per_dir)
    m = prod(map(length,weights_per_dir))
    Tv = real_type(options(geo))
    w = zeros(Tv,m)
    x = zeros(SVector{D,Tv},m)
    tensor_product!(identity,x,coords_per_dir)
    tensor_product!(prod,w,weights_per_dir)
    face_quadrature(;domain=geo,coordinates=x,weights=w)
end

function tensor_product!(f,result,values_per_dir)
    shape = Tuple(map(length,values_per_dir))
    cis = CartesianIndices(shape)
    lis = LinearIndices(cis)
    for ci in cis
        li = lis[ci]
        result[li] = f(map((q,i)->q[i],values_per_dir,Tuple(ci)))
    end
    result
end

"""
"""
function duffy_quadrature(geo,degree)
    @assert is_unit_simplex(geo)
    Tv = real_type(options(geo))
    D = num_dims(geo)
    if D == 0
        x = zeros(SVector{0,Tv},1)
        w = ones(Tv,1)
        return face_quadrature(;domain=geo,coordinates=x,weights=w)
    end
    function map_to(a,b,(points,weights))
      points_ab = similar(points)
      weights_ab = similar(weights)
      points_ab .= 0.5*(b-a)*points .+ 0.5*(a+b)
      weights_ab .= 0.5*(b-a)*weights
      (points_ab, weights_ab)
    end
    function duffy_map(q)
        D = length(q)
        a = 1.0
        m = ntuple(Val(D)) do i
            if i == 1
                q[i]
            else
                a *= (1-q[i-1])
                a*q[i]
            end
        end
        typeof(q)(m)
    end
    Ti = int_type(options(geo))
    n = ceil(Ti, (degree + 1.0) / 2.0 )
    beta = 0
    dim_to_quad_1d = map(1:(D-1)) do d
        alpha = (D-1)-(d-1)
        map_to(0,1,gaussjacobi(n,alpha,beta))
    end
    quad_1d = map_to(0,1,gausslegendre(n))
    push!(dim_to_quad_1d,quad_1d)
    coords_per_dir = map(first,dim_to_quad_1d)
    weights_per_dir =  map(last,dim_to_quad_1d)
    a = 0.5
    for d in (D-1):-1:1
        ws_1d = weights_per_dir[d]
        ws_1d[:] *= a
        a *= 0.5
    end
    m = prod(map(length,weights_per_dir))
    w = zeros(Tv,m)
    x = zeros(SVector{D,Tv},m)
    tensor_product!(identity,x,coords_per_dir)
    tensor_product!(prod,w,weights_per_dir)
    x .= duffy_map.(x)
    face_quadrature(;domain=geo,coordinates=x,weights=w)
end

"""
    abstract type AbstractMesh

# Basic queries

- [`node_coordinates`](@ref)
- [`face_nodes`](@ref)
- [`face_reference_id`](@ref)
- [`reference_spaces`](@ref)
- [`periodic_nodes`](@ref)
- [`physical_faces`](@ref)
- [`physical_nodes`](@ref)
- [`outward_normals`](@ref)

# Basic constructors

- [`mesh_from_arrays`](@ref)
- [`mesh_from_gmsh`](@ref)
- [`cartesian_mesh`](@ref)
- [`mesh_from_chain`](@ref)

"""
abstract type AbstractMesh <: AbstractType end

num_dims(m::AbstractMesh) = length(reference_spaces(m))
num_ambient_dims(m::AbstractMesh) = length(eltype(node_coordinates(m)))
options(m::AbstractMesh) = options(first(last(reference_spaces(m))))

function mesh_from_arrays(;
        node_coordinates,
        face_nodes,
        face_reference_id,
        reference_spaces,
        periodic_nodes = :default,
        physical_faces = :default,
        outward_normals = :default,
    )
    contents = (;
                node_coordinates,
                face_nodes,
                face_reference_id,
                reference_spaces,
                periodic_nodes,
                physical_faces,
                outward_normals,
               )
    Mesh(contents)
end

struct Mesh{A} <: AbstractMesh
    contents::A
end

node_coordinates(m::Mesh) = m.contents.node_coordinates
face_nodes(m::Mesh) = m.contents.face_nodes
face_nodes(m::Mesh,d) = m.contents.face_nodes[d+1]
face_reference_id(m::Mesh) = m.contents.face_reference_id
face_reference_id(m::Mesh,d) = m.contents.face_reference_id[d+1]
reference_spaces(m::Mesh) = m.contents.reference_spaces
reference_spaces(m::Mesh,d) = m.contents.reference_spaces[d+1]

function reference_domains(a::AbstractMesh,d)
    map(domain,reference_spaces(a,d))
end

function reference_domains(a::AbstractMesh)
    ntuple(d->reference_domains(a,d-1),Val{D+1})
end

function periodic_nodes(m::Mesh)
    if m.contents.periodic_nodes === :default
        Ti = int_type(options(mesh))
        Ti[] => Ti[]
    else
        m.contents.periodic_nodes
    end
end

function physical_faces(m::Mesh,d)
    if m.contents.physical_faces === :default
        Ti = int_type(options(mesh))
        Dict{String,Vector{Ti}}("$d-faces"=>collect(Ti,1:num_faces(m,d)))
    else
        m.contents.physical_faces[d+1]
    end
end

function physical_faces(m::Mesh)
    if m.contents.physical_faces === :default
        ntuple(d->physical_faces(m,d-1),Val{D+1})
    else
        m.contents.physical_faces
    end
end

function outward_normals(m::Mesh)
    @assert m.contents.outward_normals !== :default
    m.contents.outward_normals
end

"""
"""
function physical_names(mesh,d)
    groups = physical_faces(mesh,d)
    Set(keys(groups))
end

function physical_names(mesh;merge_dims=Val(false))
    D = num_dims(mesh)
    d_to_names = [ physical_names(mesh,d) for d in 0:D]
    if val_parameter(merge_dims) == false
        return d_to_names
    end
    reduce(union,d_to_names)
end

function complexify(geom::Union{UnitNCube{0},UnitSimplex{0}})
    Tv = real_type(options(geom))
    Ti = int_type(options(geom))
    Tr = reference_int_type(options(geom))
    node_coordinates = [SVector{0,Tv}()]
    face_nodes = [Ti[1]]
    face_reference_id = [Tr[1]]
    space = lagrange_space(geom)
    reference_spaces = ((space,),)
    mesh_from_arrays(;
                node_coordinates,
                face_nodes,
                face_reference_id,
                reference_spaces)
end

function complexify(geom::Union{UnitNCube{1},UnitSimplex{1}})
    Tv = real_type(options(geom))
    Ti = int_type(options(geom))
    Tr = reference_int_type(options(geom))
    geom0 = if is_unit_n_cube(geom)
        unit_n_cube(Val(0);options=options(geom))
    else
        unit_simplex(Val(0);options=options(geom))
    end
    geom1 = geom
    space0 = lagrange_space(geom0)
    space1 = lagrange_space(geom1)
    node_coordinates = SVector{1,Tv}[(0,),(1,)]
    face_nodes = Vector{Vector{Ti}}[[[1],[2]],[[1,2]]]
    face_reference_id = Vector{Tr}[[1,1],[1]]
    reference_spaces = ((space0,),(space1,))
    outwards_normals = SVector{1,Tv}[(-1,),(1,)]
    mesh_from_arrays(;
                node_coordinates,
                face_nodes,
                face_reference_id,
                reference_spaces,
                outward_normals
               )
end

function complexify(geom::UnitSimplex{2})
    Tv = real_type(options(geom))
    Ti = int_type(options(geom))
    Tr = reference_int_type(options(geom))
    geom0 = unit_simplex(Val(0);options=options(geom))
    geom1 = unit_simplex(Val(1);options=options(geom))
    geom2 = geom
    space0 = lagrange_space(geom0)
    space1 = lagrange_space(geom1)
    space2 = lagrange_space(geom2)
    node_coordinates = SVector{2,Tv}[(0,0),(1,0),(0,1)]
    face_nodes = Vector{Vector{Ti}}[[[1],[2],[3]],[[1,2],[1,3],[2,3]],[[1,2,3]]]
    face_reference_id = Vector{Tr}[[1,1,1],[1,1,1],[1]]
    reference_spaces = ((space0,),(space1,),(space2,))
    n1 = sqrt(2)/2
    outwards_normals = SVector{2,Tv}[(0,-1),(-1,0),(n1,n1)]
    mesh_from_arrays(;
                node_coordinates,
                face_nodes,
                face_reference_id,
                reference_spaces,
                outward_normals
               )
end

function complexify(geom::UnitNCube{2})
    Tv = real_type(options(geom))
    Ti = int_type(options(geom))
    Tr = reference_int_type(options(geom))
    geom0 = unit_n_cube(Val(0);options=options(geom))
    geom1 = unit_n_cube(Val(1);options=options(geom))
    geom2 = geom
    space0 = lagrange_space(geom0)
    space1 = lagrange_space(geom1)
    space2 = lagrange_space(geom2)
    node_coordinates = SVector{2,Tv}[(0,0),(1,0),(0,1),(1,1)]
    face_nodes = Vector{Vector{Ti}}[[[1],[2],[3],[4]],[[1,2],[3,4],[1,3],[2,4]],[[1,2,3,4]]]
    face_reference_id = Vector{Tr}[[1,1,1,1],[1,1,1,1],[1]]
    outwards_normals = SVector{2,Tv}[(0,-1),(0,1),(-1,0),(1,0)]
    reference_spaces = ((space0,),(space1,),(space2,))
    mesh_from_arrays(;
                node_coordinates,
                face_nodes,
                face_reference_id,
                reference_spaces,
                outward_normals
               )
end

function complexify(geom::UnitSimplex{3})
    Tv = real_type(options(geom))
    Ti = int_type(options(geom))
    Tr = reference_int_type(options(geom))
    geom0 = unit_simplex(Val(0);options=options(geom))
    geom1 = unit_simplex(Val(1);options=options(geom))
    geom2 = unit_simplex(Val(2);options=options(geom))
    geom3 = geom
    space0 = lagrange_space(geom0)
    space1 = lagrange_space(geom1)
    space2 = lagrange_space(geom2)
    space3 = lagrange_space(geom3)
    node_coordinates = SVector{3,Tv}[(0,0,0),(1,0,0),(0,1,0),(0,0,1)]
    face_nodes = [
                  Vector{Ti}[[1],[2],[3],[4]],
                  Vector{Ti}[[1,2],[1,3],[2,3],[1,4],[2,4],[3,4]],
                  Vector{Ti}[[1,2,3],[1,2,4],[1,3,4],[2,3,4]],
                  Vector{Ti}[[1,2,3,4]]
                 ]
    face_reference_id = [ones(Tr,4),ones(Tr,6),ones(Tr,4),ones(Tr,1)]
    n1 = sqrt(3)/3
    outwards_normals = SVector{3,Tv}[(0,0,-1),(0,-1,0),(-1,0,0),(n1,n1,n1)]
    reference_spaces = ((space0,),(space1,),(space2,),(space3,))
    mesh_from_arrays(;
                node_coordinates,
                face_nodes,
                face_reference_id,
                reference_spaces,
                outward_normals
               )
end

function complexify(geom::UnitNCube{3})
    Tv = real_type(options(geom))
    Ti = int_type(options(geom))
    Tr = reference_int_type(options(geom))
    geom0 = unit_n_cube(Val(0);options=options(geom))
    geom1 = unit_n_cube(Val(1);options=options(geom))
    geom2 = unit_n_cube(Val(2);options=options(geom))
    geom3 = geom
    space0 = lagrange_space(geom0)
    space1 = lagrange_space(geom1)
    space2 = lagrange_space(geom2)
    space3 = lagrange_space(geom3)
    node_coordinates = SVector{3,Tv}[(0,0,0),(1,0,0),(0,1,0),(1,1,0),(0,0,1),(1,0,1),(0,1,1),(1,1,1)]
    face_nodes = [
                  Vector{Ti}[[1],[2],[3],[4],[5],[6],[7],[8]],
                  Vector{Ti}[[1,2],[3,4],[1,3],[2,4],[5,6],[7,8],[5,7],[6,8],[1,5],[3,7],[2,6],[4,8]],
                  Vector{Ti}[[1,2,3,4],[5,6,7,8],[1,2,5,6],[3,4,7,8],[1,3,5,7],[2,4,6,8]],
                  Vector{Ti}[[1,2,3,4,5,6,7,8]],
                 ]
    face_reference_id = [ones(Tr,8),ones(Tr,12),ones(Tr,6),ones(Tr,1)]
    outwards_normals = SVector{3,Tv}[(0,0,-1),(0,0,1),(0,-1,0),(0,1,0),(-1,0,0),(1,0,0)]
    reference_spaces = ((space0,),(space1,),(space2,),(space3,))
    mesh_from_arrays(;
                node_coordinates,
                face_nodes,
                face_reference_id,
                reference_spaces,
                outward_normals
               )
end

function opposite_faces(geom::Union{UnitNCube{0},UnitSimplex{0}})
    Ti = int_type(options(geom))
    Vector{Ti}[[1]]
end

function opposite_faces(geom::Union{UnitNCube{1},UnitSimplex{1}})
    Ti = int_type(options(geom))
    Vector{Ti}[[2,1],[1]]
end

function opposite_faces(geom::UnitNCube{2})
    Ti = int_type(options(geom))
    Vector{Ti}[[4,3,2,1],[2,1,4,3],[1]]
end

function opposite_faces(geom::UnitNCube{3})
    Ti = int_type(options(geom))
    Vector{Ti}[[8,7,6,5,4,3,2,1],[6,5,8,7,2,1,4,3,12,11,10,9],[2,1,4,3,6,5],[1]]
end

abstract type AbstractMeshDomain <: AbstractDomain end

mesh(a::AbstractMeshDomain) = a.mesh
mesh_id(a::AbstractMeshDomain) = a.mesh_id
physical_names(a::AbstractMeshDomain) = a.physical_names
num_dims(a::AbstractMeshDomain) = GT.val_parameter(a.num_dims)
num_ambient_dims(a::AbstractMeshDomain) = num_ambient_dims(mesh(a))
is_physical_domain(a::AbstractMeshDomain) = ! is_reference_domain(a)
face_around(a::AbstractMeshDomain) = a.face_around

function is_boundary(dom::AbstractMeshDomain)
    dom.face_around !== nothing && (num_dims(dom) + 1) == num_dims(mesh(dom))
end

function inverse_faces(domain::AbstractMeshDomain)
    if domain.cache !== nothing
        return domain.cache.inverse_faces
    end
    Ti = int_type(options(domain))
    d = num_dims(domain)
    ndfaces = num_faces(mesh(domain),d)
    dface_to_face = zeros(Ti,ndfaces)
    face_to_dface = faces(domain)
    dface_to_face[face_to_dface] = 1:length(face_to_dface)
    dface_to_face
end

function faces(domain::AbstractMeshDomain)
    if domain.cache !== nothing
        return domain.cache.faces
    end
    Ti = int_type(options(domain))
    mesh = domain |> GT.mesh
    D = GT.num_dims(domain)
    Dface_to_tag = zeros(Ti,GT.num_faces(mesh,D))
    tag_to_name = GT.physical_names(domain)
    GT.classify_mesh_faces!(Dface_to_tag,mesh,D,tag_to_name)
    physical_Dfaces = findall(i->i!=0,Dface_to_tag)
    Ti.(physical_Dfaces)
end

"""
"""
function classify_mesh_faces!(dface_to_tag,mesh::AbstractMesh,d,tag_to_name)
    fill!(dface_to_tag,zero(eltype(dface_to_tag)))
    face_groups = physical_faces(mesh,d)
    for (tag,name) in enumerate(tag_to_name)
        for (name2,faces) in face_groups
            if name != name2
                continue
            end
            dface_to_tag[faces] .= tag
        end
    end
    dface_to_tag
end

function face_reference_id(a::AbstractMeshDomain)
    d = num_dims(a)
    face_reference_id(mesh(a),d)
end

function reference_domains(a::AbstractMeshDomain)
    d = num_dims(a)
    reference_domains(mesh(a),d)
end

function Base.:(==)(a::AbstractMeshDomain,b::AbstractMeshDomain)
    flag = true
    flag = flag && (GT.mesh_id(a) == GT.mesh_id(b))
    flag = flag && (GT.physical_names(a) == GT.physical_names(b))
    flag = flag && (GT.num_dims(a) == GT.num_dims(b))
    flag = flag && (GT.is_reference_domain(a) == GT.is_reference_domain(b))
    flag
end

function domain(mesh::AbstractMesh,d;
    mesh_id = objectid(mesh),
    physical_names=GT.physical_names(mesh,d),
    is_reference_domain = Val(false))
    mesh_domain(;mesh,num_dims=Val(val_parameter(d)),face_around=nothing,mesh_id,physical_names,is_reference_domain)
end

function interior(mesh::AbstractMesh;
    mesh_id = objectid(mesh),
    physical_names=GT.physical_names(mesh,num_dims(mesh)),
    is_reference_domain = Val(false))
    D = num_dims(mesh)
    mesh_domain(;mesh,num_dims=D,face_around=1,mesh_id,physical_names,is_reference_domain)
end

function skeleton(mesh::AbstractMesh;
    mesh_id = objectid(mesh),
    physical_names=GT.physical_names(mesh,num_dims(mesh)-1),
    is_reference_domain = Val(false))
    D = num_dims(mesh)
    mesh_domain(;mesh,num_dims=D-1,face_around=nothing,mesh_id,physical_names,is_reference_domain)
end

function boundary(mesh::AbstractMesh;
    face_around=1,
    mesh_id = objectid(mesh),
    physical_names=GT.physical_names(mesh,num_dims(mesh)-1),
    is_reference_domain = Val(false))
    D = num_dims(mesh)
    mesh_domain(;mesh,num_dims=D-1,face_around,mesh_id,physical_names,is_reference_domain)
end

function mesh_domain(;
    mesh,
    mesh_id = objectid(mesh),
    num_dims = Val(GT.num_dims(mesh)),
    physical_names=GT.physical_names(mesh,num_dims),
    is_reference_domain = Val(false),
    face_around=nothing,
    cache=nothing,
    )

    if val_parameter(is_reference_domain)
        ReferenceDomain(
                        mesh,
                        mesh_id,
                        physical_names,
                        Val(val_parameter(num_dims)),
                        face_around,
                        cache,
                       )
    else
        PhysicalDomain(
                        mesh,
                        mesh_id,
                        physical_names,
                        Val(val_parameter(num_dims)),
                        face_around,
                        cache,
                       )
    end |> setup_domain
end

function setup_domain(domain::AbstractMeshDomain)
    if domain.cache !== nothing
        return domain
    end
    faces = GT.faces(domain)
    inverse_faces = GT.inverse_faces(domain)
    cache = (;faces,inverse_faces)
    replace_cache(domain,cache)
end

function replace_cache(domain::AbstractMeshDomain,cache)
    mesh = GT.mesh(domain)
    num_dims = GT.num_dims(domain)
    mesh_id = GT.mesh_id(domain)
    physical_names = GT.physical_names(domain)
    is_reference_domain = GT.is_reference_domain(domain)
    face_around = GT.face_around(domain)
    GT.mesh_domain(;mesh,num_dims,mesh_id,physical_names,is_reference_domain,face_around,cache)
end

struct PhysicalDomain{A,B,C,D,E,F} <: AbstractMeshDomain{A}
    mesh::A
    mesh_id::B
    physical_names::C
    num_dims::Val{D}
    face_around::E
    cache::F
end
is_reference_domain(a::PhysicalDomain) = false

struct ReferenceDomain{A,B,C,D,E,F} <: AbstractMeshDomain{A}
    mesh::A
    mesh_id::B
    physical_names::C
    num_dims::Val{D}
    face_around::E
    cache::F
end
is_reference_domain(a::ReferenceDomain) = true

function reference_domain(domain::PhysicalDomain)
    mesh = GT.mesh(domain)
    num_dims = GT.num_dims(domain)
    mesh_id = GT.mesh_id(domain)
    physical_names = GT.physical_names(domain)
    is_reference_domain = Val(true)
    face_around = GT.face_around(domain)
    cache = domain.cache
    GT.domain(;mesh,num_dims,mesh_id,physical_names,is_reference_domain,face_around,cache)
end

function reference_domain(domain::ReferenceDomain)
    domain
end

function physical_domain(domain::ReferenceDomain)
    mesh = GT.mesh(domain)
    num_dims = GT.num_dims(domain)
    mesh_id = GT.mesh_id(domain)
    physical_names = GT.physical_names(domain)
    face_around = GT.face_around(domain)
    is_reference_domain = Val(false)
    cache = domain.cache
    GT.domain(;mesh,num_dims,mesh_id,physical_names,is_reference_domain,face_around,cache)
end

function physical_domain(domain::PhysicalDomain)
    domain
end

abstract type AbstractMeshQuadrature <: AbstractQuadrature end

function quadrature(domain::AbstractMeshDomain,degree)
    rid_to_dom = reference_domains(domain)
    reference_quadratures = map(dom->quadrature(dom,degree),rid_to_dom)
    face_reference_id = GT.face_reference_id(domain)
    mesh_quadrature(;domain,face_reference_id,reference_quadratures)
end

function node_quadrature(domain::AbstractMeshDomain)
    d = num_dims(domain)
    rid_to_dom = reference_spaces(mesh(domain),d)
    reference_quadratures = map(node_quadrature,rid_to_dom)
    face_reference_id = GT.face_reference_id(domain)
    mesh_quadrature(;domain,face_reference_id,reference_quadratures)
end

function mesh_quadrature(;domain,face_reference_id,reference_quadratures)
    contents = (;domain,face_reference_id,reference_quadratures)
    MeshQuadrature(contents)
end

struct MeshQuadrature{A} <: AbstractMeshQuadrature
    contents::A
end

domain(q::MeshQuadrature) = q.contents.domain
face_reference_id(q::MeshQuadrature) = q.contents.face_reference_id
reference_quadratures(q::MeshQuadrature) = q.contents.reference_quadratures

function num_points_accessor(measure::MeshQuadrature)
    mesh = measure.mesh
    dom = measure.domain
    d = num_dims(dom)
    rid_to_point_to_w = map(weights,reference_quadratures(measure))
    face_to_rid = face_reference_id(measure)
    sface_to_face = faces(dom)
    function face_npoints(sface)
        face = sface_to_face[sface]
        rid = face_to_rid[face]
        point_to_w = rid_to_point_to_w[rid]
        length(point_to_w)
    end
end

function coordinate_accessor(measure::MeshQuadrature)
    mesh = measure.mesh
    dom = measure.domain
    @assert is_physical_domain(dom)
    d = num_dims(dom)
    rid_to_point_to_x = map(coordinates,reference_quadratures(measure))
    rid_to_tab = map(rid_to_point_to_x,reference_faces(mesh,d)) do point_to_x, refface
        tabulator(refface)(value,point_to_x)
    end
    face_to_rid = face_reference_id(measure)
    face_to_nodes = face_nodes(mesh,d)
    node_to_x = node_coordinates(mesh)
    sface_to_face = faces(dom)
    function face_point_x(sface)
        face = sface_to_face[sface]
        rid = face_to_rid[face]
        tab = rid_to_tab[rid]
        nodes = face_to_nodes[face]
        function point_x(point)
            nnodes = length(nodes)
            sum(1:nnodes) do i
                node = nodes[i]
                x = node_to_x[node]
                x*tab[point,i]
            end
        end
    end
end

outer(a,b) = a*transpose(b)

function jacobian_accessor(measure::MeshQuadrature)
    mesh = measure.mesh
    dom = measure.domain
    @assert is_physical_domain(dom)
    d = num_dims(dom)
    rid_to_point_to_x = map(coordinates,reference_quadratures(measure))
    rid_to_tab = map(rid_to_point_to_x,reference_faces(mesh,d)) do point_to_x, refface
        tabulator(refface)(ForwardDiff.gradient,point_to_x)
    end
    face_to_rid = face_reference_id(measure)
    face_to_nodes = face_nodes(mesh,d)
    node_to_x = node_coordinates(mesh)
    sface_to_face = faces(dom)
    function face_point_x(sface)
        face = sface_to_face[sface]
        rid = face_to_rid[face]
        tab = rid_to_tab[rid]
        nodes = face_to_nodes[face]
        function point_x(point)
            nnodes = length(nodes)
            sum(1:nnodes) do i
                node = nodes[i]
                x = node_to_x[node]
                outer(x,tab[point,i])
            end
        end
    end
end

function weight_accessor(measure::MeshQuadrature)
    mesh = measure.mesh
    dom = measure.domain
    @assert is_physical_domain(dom)
    d = num_dims(dom)
    rid_to_point_to_w = map(weights,reference_quadratures(measure))
    face_to_rid = face_reference_id(measure)
    sface_to_face = faces(dom)
    function face_point_x(sface)
        face = sface_to_face[sface]
        rid = face_to_rid[face]
        point_to_w = rid_to_point_to_w[rid]
        function point_x(point,J)
            w = point_to_w[point]
            change_of_measure(J)*w
        end
    end
end




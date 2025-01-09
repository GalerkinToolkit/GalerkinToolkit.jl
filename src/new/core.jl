

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
    contents = (;
                reference_int_type=Val(reference_int_type),
                int_type=Val(int_type),
                global_int_type=Val(global_int_type),
                real_type=Val(real_type),
               )
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
reference_int_type(options::Options) = val_parameter(options.contents.reference_int_type)

"""
    int_type(options::Options)

Return the default integer type used in the computation except for reference and global quantities.
"""
int_type(options::Options) = val_parameter(options.contents.int_type)

"""
    global_int_type(options::Options)

Return the type of the integers used to enumerate global quantities.
"""
global_int_type(options::Options) = val_parameter(options.contents.global_int_type)

"""
    real_type(options::Options)

Return the default real type used in the computation.
"""
real_type(options::Options) = val_parameter(options.contents.real_type)

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
    domain = GT.domain(fe)
    face_quadrature(;domain,coordinates,weights)
end

function lagrange_space(domain::AbstractFaceDomain;
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
- [`node_quadrature`](@ref)

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
        if degree in strang_quadrature_degrees(geo)
            strang_quadrature(geo,degree)
        else
            duffy_quadrature(geo,degree)
        end
    else
        error("Not implemented")
    end
end

abstract type AbstractDiscretization end

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

- [`mesh`](@ref)
- [`mesh_from_gmsh`](@ref)
- [`cartesian_mesh`](@ref)

"""
abstract type AbstractMesh <: AbstractDiscretization end

num_dims(m::AbstractMesh) = length(reference_spaces(m))-1
num_ambient_dims(m::AbstractMesh) = length(eltype(node_coordinates(m)))
options(m::AbstractMesh) = options(first(last(reference_spaces(m))))
num_faces(m::AbstractMesh,d) = length(face_reference_id(m,d))

struct Mesh{A} <: AbstractMesh
    contents::A
end

function mesh(;
        node_coordinates,
        face_nodes,
        face_reference_id,
        reference_spaces,
        periodic_nodes = default_periodic_nodes(reference_spaces),
        physical_faces = default_physical_faces(reference_spaces),
        outward_normals = nothing,
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

node_coordinates(m::Mesh) = m.contents.node_coordinates
face_nodes(m::Mesh) = m.contents.face_nodes
face_nodes(m::Mesh,d) = m.contents.face_nodes[d+1]
face_reference_id(m::Mesh) = m.contents.face_reference_id
face_reference_id(m::Mesh,d) = m.contents.face_reference_id[d+1]
reference_spaces(m::Mesh) = m.contents.reference_spaces
reference_spaces(m::Mesh,d) = m.contents.reference_spaces[d+1]
physical_faces(m::Mesh) = m.contents.physical_faces
physical_faces(m::Mesh,d) = m.contents.physical_faces[d+1]
periodic_nodes(m::Mesh) = m.contents.periodic_nodes

function default_physical_faces(reference_spaces)
    [ Dict{String,Vector{int_type(options(first(last(reference_spaces))))}}() for _ in 1:length(reference_spaces) ]
end

function default_periodic_nodes(reference_spaces)
    Ti = int_type(options(first(last(reference_spaces))))
    Ti[] => Ti[]
end

function reference_domains(a::AbstractMesh,d)
    map(domain,reference_spaces(a,d))
end

function reference_domains(a::AbstractMesh)
    D = num_dims(a)
    ntuple(d->reference_domains(a,d-1),Val(D+1))
end

function outward_normals(m::Mesh)
    @assert m.contents.outward_normals !== nothing
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

abstract type AbstractChain <: AbstractDiscretization end

num_dims(m::AbstractChain) = num_dims(domain(first(reference_spaces(m))))
num_ambient_dims(m::AbstractChain) = length(eltype(node_coordinates(m)))
options(m::AbstractChain) = options(first(reference_spaces(m)))
num_faces(m::AbstractChain) = length(face_reference_id(m))

"""
"""
function mesh(chain::AbstractChain)
    D = num_dims(chain)
    cell_nodes = face_nodes(chain)
    cell_reference_id = face_reference_id(chain)
    reference_cells = reference_spaces(chain)
    node_coords = node_coordinates(chain)
    face_to_nodes = Vector{typeof(cell_nodes)}(undef,D+1)
    face_to_refid = Vector{typeof(cell_reference_id)}(undef,D+1)
    for d in 0:D-1
        face_to_nodes[d+1] = Vector{Int}[]
        face_to_refid[d+1] = Int[]
    end
    face_to_nodes[end] = cell_nodes
    face_to_refid[end] = cell_reference_id
    ref_cell = first(reference_cells)
    ref_faces = reference_spaces(complexify(ref_cell))
    refid_to_refface = push(ref_faces[1:end-1],reference_cells)
    cell_groups = physical_faces(chain)
    groups = [ typeof(cell_groups)() for d in 0:D]
    groups[end] = cell_groups
    pnodes = periodic_nodes(chain)
    onormals = outward_normals(chain)
    GT.mesh(;
      node_coordinates = node_coords,
      face_nodes = face_to_nodes,
      face_reference_id = face_to_refid,
      reference_spaces = refid_to_refface,
      periodic_nodes = pnodes,
      physical_faces = groups,
      outward_normals = onormals)
end

struct Chain{A} <: AbstractChain
    contents::A
end

function chain(;
        node_coordinates,
        face_nodes,
        face_reference_id,
        reference_spaces,
        periodic_nodes = default_periodic_nodes((reference_spaces,)),
        physical_faces = default_physical_faces((reference_spaces,))[end],
        outward_normals = nothing,
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
    Chain(contents)
end

node_coordinates(m::Chain) = m.contents.node_coordinates
face_nodes(m::Chain) = m.contents.face_nodes
face_reference_id(m::Chain) = m.contents.face_reference_id
reference_spaces(m::Chain) = m.contents.reference_spaces
physical_faces(m::Chain) = m.contents.physical_faces
periodic_nodes(m::Chain) = m.contents.periodic_nodes
outward_normals(m::Chain) = m.contents.outward_normals

abstract type AbstractMeshDomain{A} <: AbstractDomain end

num_ambient_dims(a::AbstractMeshDomain) = num_ambient_dims(mesh(a))
is_physical_domain(a::AbstractMeshDomain) = ! is_reference_domain(a)
options(a::AbstractMeshDomain) = options(mesh(a))

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

function is_boundary(dom::AbstractMeshDomain)
    face_around(dom) !== nothing && (num_dims(dom) + 1) == num_dims(mesh(dom))
end

function label_faces_in_dim!(m::AbstractMesh,d;physical_name="__$d-FACES__")
    groups = physical_faces(m,d)
    if haskey(groups,physical_name)
        return m
    end
    Ti = int_type(options(m))
    faces = collect(Ti,1:num_faces(m,d))
    groups[physical_name] = faces
    m
end

function domain(mesh::AbstractMesh,d;is_reference_domain=Val(false),physical_name="__$d-FACES__")
    label_faces_in_dim!(mesh,d;physical_name)
    mesh_domain(;
        mesh,
        num_dims=Val(val_parameter(d)),
        physical_names=[physical_name],
        is_reference_domain)
end

struct PhysicalDomain{A,B} <: AbstractMeshDomain{A}
    mesh::A
    contents::B
end
is_reference_domain(a::PhysicalDomain) = false

struct ReferenceDomain{A,B} <: AbstractMeshDomain{A}
    mesh::A
    contents::B
end
is_reference_domain(a::ReferenceDomain) = true

const MeshDomain = Union{PhysicalDomain{A,B},ReferenceDomain{A,B}} where {A,B} 

function mesh_domain(;
    mesh,
    mesh_id = objectid(mesh),
    num_dims = Val(GT.num_dims(mesh)),
    physical_names=GT.physical_names(mesh,num_dims),
    is_reference_domain = Val(false),
    face_around=nothing,
    workspace=nothing,
    )

    if val_parameter(is_reference_domain)
        constructor = ReferenceDomain
    else
        constructor = PhysicalDomain
    end
    contents = (;
                mesh_id,
                physical_names,
                num_dims = Val(val_parameter(num_dims)),
                face_around,
                workspace,
               )
    constructor(mesh,contents) |> setup_domain
end

function reference_domain(domain::PhysicalDomain)
    mesh = GT.mesh(domain)
    num_dims = GT.num_dims(domain)
    mesh_id = GT.mesh_id(domain)
    physical_names = GT.physical_names(domain)
    is_reference_domain = Val(true)
    face_around = GT.face_around(domain)
    workspace = domain.workspace
    GT.mesh_domain(;mesh,num_dims,mesh_id,physical_names,is_reference_domain,face_around,workspace)
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
    workspace = domain.workspace
    GT.mesh_domain(;mesh,num_dims,mesh_id,physical_names,is_reference_domain,face_around,workspace)
end

function physical_domain(domain::PhysicalDomain)
    domain
end

mesh(a::MeshDomain) = a.mesh
mesh_id(a::MeshDomain) = a.contents.mesh_id
physical_names(a::MeshDomain) = a.contents.physical_names
face_around(a::MeshDomain) = a.contents.face_around
num_dims(a::MeshDomain) = GT.val_parameter(a.contents.num_dims)
workspace(a::MeshDomain) = a.contents.workspace

function setup_domain(domain::MeshDomain)
    if GT.workspace(domain) !== nothing
        return domain
    end
    faces = GT.faces(domain)
    inverse_faces = GT.inverse_faces(domain)
    workspace = (;faces,inverse_faces)
    replace_workspace(domain,workspace)
end

function replace_workspace(domain::MeshDomain,workspace)
    mesh = GT.mesh(domain)
    num_dims = GT.num_dims(domain)
    mesh_id = GT.mesh_id(domain)
    physical_names = GT.physical_names(domain)
    is_reference_domain = GT.is_reference_domain(domain)
    face_around = GT.face_around(domain)
    GT.mesh_domain(;mesh,num_dims,mesh_id,physical_names,is_reference_domain,face_around,workspace)
end

function faces(domain::MeshDomain)
    if workspace(domain) !== nothing
        return workspace(domain).faces
    end
    Ti = int_type(options(domain))
    mesh = domain |> GT.mesh
    D = GT.num_dims(domain)
    Dface_to_tag = zeros(Ti,GT.num_faces(mesh,D))
    tag_to_name = GT.physical_names(domain)
    fill!(Dface_to_tag,zero(eltype(Dface_to_tag)))
    face_groups = physical_faces(mesh,D)
    for (tag,name) in enumerate(tag_to_name)
        for (name2,faces) in face_groups
            if name != name2
                continue
            end
            Dface_to_tag[faces] .= tag
        end
    end
    physical_Dfaces = findall(i->i!=0,Dface_to_tag)
    Ti.(physical_Dfaces)
end

function inverse_faces(domain::MeshDomain)
    if workspace(domain) !== nothing
        return workspace(domain).inverse_faces
    end
    Ti = int_type(options(domain))
    d = num_dims(domain)
    ndfaces = num_faces(mesh(domain),d)
    dface_to_face = zeros(Ti,ndfaces)
    face_to_dface = faces(domain)
    dface_to_face[face_to_dface] = 1:length(face_to_dface)
    dface_to_face
end

abstract type AbstractMeshQuadrature{A} <: AbstractQuadrature end

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

struct MeshQuadrature{A,B} <: AbstractMeshQuadrature{A}
    mesh::A
    contents::B
end

function mesh_quadrature(;domain,face_reference_id,reference_quadratures)
    contents = (;domain,face_reference_id,reference_quadratures)
    MeshQuadrature(mesh(domain),contents)
end

mesh(q::MeshQuadrature) = q.mesh
domain(q::MeshQuadrature) = q.contents.domain
face_reference_id(q::MeshQuadrature) = q.contents.face_reference_id
reference_quadratures(q::MeshQuadrature) = q.contents.reference_quadratures

function num_points_accessor(measure::MeshQuadrature)
    mesh = GT.mesh(measure)
    dom = GT.domain(measure)
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
    ## TODO the following ones assume that the same reference ids are for the mesh
    # as for the quadrature. This needs to be generalized for adaptive quadrature schemes.
    mesh = GT.mesh(measure)
    dom = GT.domain(measure)
    @assert is_physical_domain(dom)
    d = num_dims(dom)
    rid_to_point_to_x = map(coordinates,reference_quadratures(measure))
    rid_to_tab = map(rid_to_point_to_x,reference_spaces(mesh,d)) do point_to_x, refface
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
    mesh = GT.mesh(measure)
    dom = GT.domain(measure)
    @assert is_physical_domain(dom)
    d = num_dims(dom)
    rid_to_point_to_x = map(coordinates,reference_quadratures(measure))
    rid_to_tab = map(rid_to_point_to_x,reference_spaces(mesh,d)) do point_to_x, refface
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
    mesh = GT.mesh(measure)
    dom = GT.domain(measure)
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

function change_of_measure(J)
    sqrt(det(transpose(J)*J))
end



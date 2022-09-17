"""
   const INVALID_ID::Int

Stores an invalid linear index for 1-based indexed objects.
`INVALID_ID` is some number less than `1`. The user should not rely on the
specific value.
"""
const INVALID_ID = -1

"""
    num_dims(geo)
Return the number of space dimensions of the geometrical object `geo`.
"""
function num_dims end

"""
    num_ambient_dims(geo)
Return the number of dimensions of the ambient space,
where the geometrical object `geo` is embedded.
"""
function num_ambient_dims end
num_ambient_dims(geo) = length(eltype(node_coordinates(geo)))

"""
    num_faces(geo,dim)
Number of faces of dimension `dim` in `geo`.
"""
function num_faces end
num_faces(geo,dim) = length(face_ref_id(geo,dim))

"""
    face_nodes(geo,dim)
Return the face node connectivity for faces of dimension `dim` in `geo`.
The result is a vector of vectors of outer length `num_faces(geo,dim)`
that can be indexed with a face id, returning the node ids of the
nodes on the *closure* of this face.
"""
function face_nodes end

"""
    face_nodes!(geo,nodes,dim)
Sets the face node connectivity for faces of dimension `dim` in `geo`.
"""
function face_nodes! end

"""
    face_own_nodes(geo,dim)
Like `face_nodes(geo,dim)`, but it returns node ids in the *interior*
of the faces in `geo`.
"""
function face_own_nodes end

"""
    face_faces(geo,dim_from,dim_to)
Return the face-to-face connectivity of faces of dimension `dim_from`
against faces of dimension `dim_to` in `geo`.
This is represented with a vector of vectors
of outer length `num_faces(geo,dim_from) `. For `dim_from > dim_to`,
the inner vectors contain the face ids of faces of dimension `dim_to` on
the *boundary* of the corresponding face of dimension `dim_from`.
For `dim_from > dim_to`, the results is for the *co-boundary* instead
of the boundary.
For `dim_from == dim_to`, the result is the identity.
"""
function face_faces end

"""
    face_faces!(geo,connec,dim_from,dim_to)
Set the face-to-face connectivity of faces of dimension `dim_from`
against faces of dimension `dim_to` in `geo`.
"""
function face_faces! end

"""
    face_vertices(geo,dim)
Return the face vertex connectivity for faces of dimension
`dim` in `geo`. The result is a vector of vectors
of outer length `num_faces(geo,dim)`
that can be indexed with a face id, returning the node ids of the nodes
on the *closure* of this face.
"""
function face_vertices end
face_vertices(geo) = face_faces(geo,num_dims(geo),0)

"""
    num_vertices(geo)
Number of vertices in `geo`.
"""
function num_vertices end
num_vertices(geo) = num_faces(geo,0)

"""
    num_edges(geo)
Number of edges in `geo`.
"""
function num_edges end
num_edges(geo) = num_faces(geo,1)

"""
    num_nodes(geo)
Number of nodes in `geo`.
"""
function num_nodes end
num_nodes(geo) = length(node_coordinates(geo))

"""
    node_coordinates(geo)
Return a vector of length `num_nodes(geo)` with the nodal coordinates
of the nodes in `geo`.
"""
function node_coordinates end

"""
    ref_faces(geo,dim)

Return a vector with the reference faces of dimension `dim`
associated with the geometrical object `geo`.
"""
function ref_faces end

"""
    ref_faces!(geo,reffaces,dim)

Sets the vector with the reference faces of dimension `dim`
associated with the geometrical object `geo`.
"""
function ref_faces! end

"""
    face_ref_id(geo,dim)
Return a vector of integer ids, indicating which is the reference face
for each face of dimension `dim` in the geometrical object `geo`. I.e.,

    id_to_ref_face = ref_faces(geo,dim)
    face_to_id = face_ref_id(geo,dim)
    ref_face = id_to_ref_face[face_to_id[face]]

`ref_face` is the reference face associated with face number `face`
of dimension `dim` in `geo`.
"""
function face_ref_id end

"""
    face_ref_id!(geo,refids,dim)
Set the vector of face reference ids.
"""
function face_ref_id! end

"""
    ref_face_faces(geo,dim,m,n)
Collect the `m`-face to `n`-face connectivity for all reference faces
of dimension `dim` in `geo`.
"""
function ref_face_faces(mesh,dim,m,n)
  refid_to_refface = ref_faces(mesh,dim)
  map(rf->face_faces(rf,m,n),refid_to_refface)
end

"""
    ref_face_nodes(geo,dim,n)
Collect the `n`-face node connectivity for all reference faces
of dimension `dim` in `geo`.
"""
function ref_face_nodes(mesh,dim,m)
  refid_to_refface = ref_faces(mesh,dim)
  map(rf->face_nodes(rf,m),refid_to_refface)
end

"""
    ref_face_own_nodes(geo,dim,n)
Collect the `n`-face own node connectivity for all reference faces
of dimension `dim` in `geo`.
"""
function ref_face_own_nodes(mesh,dim,m)
  refid_to_refface = ref_faces(mesh,dim)
  map(rf->face_own_nodes(rf,m),refid_to_refface)
end

"""
    nodes = periodic_nodes(geo)
Returns an array of instances of `PeriodicNode` for `geo`.
The returned object is stored in memory as an struct of arrays
(instead of an array of structs) and it is possible to get a given field
for all entries simultaneously.  For instance, for a given index  `i`
nodes[i].free == nodes.free[i]
"""
function periodic_nodes end

"""
    periodic_nodes!(geo,info)
Sets the information about periodic nodes in `geo`.
"""
function periodic_nodes! end

"""
    struct PeriodicNode{Ti,T}

Description of a periodic node.

# Properties
    dependent::Ti
    free::Ti
    coeff::T

# Supertype hierarchy

    PeriodicNode <: Any

"""
struct PeriodicNode{Ti,T}
    dependent::Ti
    free::Ti
    coeff::T
end

struct PeriodicNodeVector{Ti,T} <: AbstractVector{PeriodicNode{Ti,T}}
    dependent::Vector{Ti}
    free::Vector{Ti}
    coeff::Vector{T}
end

Base.size(a::PeriodicNodeVector) = (length(a.dependent),)
function Base.getindex(a::PeriodicNodeVector,i::Integer)
    PeriodicNode(a.dependent[i],a.free[i],a.coeff[i])
end

"""
    num_periodic_nodes(geo)
Number of periodic nodes in `geo`.
"""
function num_periodic_nodes end
num_periodic_nodes(geo) = length(periodic_nodes(geo))

"""
    has_periodic_nodes(geo)
`true` if `num_periodic_nodes(geo)>0`, `false` otherwise.
"""
function has_periodic_nodes end
has_periodic_nodes(geo) = num_periodic_nodes(geo)>0

"""
    nodes = hanging_nodes(geo)
Returns an array of instances of `HangingNode` for `geo`.
The returned object is stored in memory as an struct of arrays
(instead of an array of structs) and it is possible to get a given field
for all nodes simultaneously.  For instance, for a given index  `i`
nodes[i].free == nodes.free[i]
"""
function hanging_nodes end

"""
    struct HangingNode{Ti,T}

Description of a hanging node.

# Properties
    dependent::Ti
    free::Ti
    coeff::T

# Supertype hierarchy

    HangingNode <: Any

"""
struct HangingNode{Ti,T}
  dependent::Ti
  free::SubArray{Ti,1,Vector{Ti},Tuple{UnitRange{Int32}},true}
  coeff::SubArray{T,1,Vector{T},Tuple{UnitRange{Int32}},true}
end

struct HangingNodeVector{Ti,T} <: AbstractVector{HangingNode{Ti,T}}
  dependent::Vector{Ti}
  free::JaggedArray{Ti,Int32}
  coeff::JaggedArray{T,Int32}
end

Base.size(a::HangingNodeVector) = (length(a.dependent),)
function Base.getindex(a::HangingNodeVector,i::Integer)
    HangingNode(a.dependent[i],a.free[i],a.coeff[i])
end

"""
    hanging_nodes!(geo,info)
Sets the information about hanging nodes in `geo`.
"""
function hanging_nodes! end

"""
    num_hanging_nodes(geo)
Number of hanging nodes in `geo`.
"""
function num_hanging_nodes end
num_hanging_nodes(geo) = length(hanging_nodes(geo))

"""
    has_hanging_nodes(geo)
`true` if `num_hanging_nodes(geo)>0`, `false` otherwise.
"""
function has_hanging_nodes end
has_hanging_nodes(geo) = num_hanging_nodes(geo)>0

"""
    is_simplex(geo)
`true` if `geo` represents a simplex, `false` otherwise.
"""
function is_simplex end
is_simplex(geo) = false

"""
    is_hypercube(geo)
`true` if `geo` represents a hypercube, `false` otherwise.
"""
function is_hypercube end
is_hypercube(geo) = false

"""
    groups = physical_groups(geo)
Get an array of all physical groups in `geo`.
Each item in the array is an instance of `PhysicalGroup`.
"""
function physical_groups end

"""
    struct PhysicalGroup

Data defining a physical group.

# Properties

- `faces::Vector{Int32}` Ids of the faces defining the physical group.
- `dim::Int` Dimension of the physical group.
- `name::String` Name of the physical group.
- `id::Int32` Id (i.e., the numeric name) of the physical group.

# Supertype hierarchy

    PhysicalGroup <: Any

"""
struct PhysicalGroup
    faces::Vector{Int32}
    dim::Int
    name::String
    id::Int32
end

"""
    groups = physical_groups(geo,dim)
Get an array of all physical groups in `geo` of dimension `dim`.
"""
function physical_groups(geo,dim)
    groups = physical_groups(geo)
    filter(group->group.dim == dim,groups)
end

"""
    physical_groups!(geo,groups)
Sets the array of physical groups in `geo`.
"""
function physical_groups! end

"""
    vertex_node(geo)
Vertex to node mapping for `geo`.
The returned vector is of length `num_vertices(geo)`.
"""
function vertex_node end

"""
    node_vertex(geo)
Node to vertex mapping for `geo`.
The returned vector is of length `num_nodes(geo)`.
"""
function node_vertex end

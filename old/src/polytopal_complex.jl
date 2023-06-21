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
    ref_faces(geo,Val(dim))
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

function PeriodicNodeVector{Ti,T}()
    PeriodicNodeVector(
        Ti[],
        Ti[],
        T[])
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

function HangingNodeVector{Ti,T}()
    HangingNodeVector(
        Ti[],
        JaggedArray{Ti,Int32}(),
        JaggedArray{T,Int32}())
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
Get the array of all physical groups in `geo`.
Each item in the array is an instance of `PhysicalGroup`.
It returns an empty array by default.
The resulting object `groups` takes ownership of the groups stored internally in `geo`.
"""
function physical_groups end
#@memoize physical_groups(a) = PhysicalGroup[]

"""
    struct PhysicalGroup{T}

Data defining a physical group.

# Properties

- `faces::Vector{T}` Ids of the faces defining the physical group.
- `dim::Int` Dimension of the physical group.
- `name::String` Name of the physical group.

# Supertype hierarchy

    PhysicalGroup{T} <: Any

"""
struct PhysicalGroup{T}
    faces::Vector{T}
    dim::Int
    name::String
end

"""
    groups = physical_groups(geo,dim)
Get an array of all physical groups in `geo` of dimension `dim`.
"""
function physical_groups(geo,dim::Integer)
    groups = physical_groups(geo)
    filter(group->group.dim == dim,groups)
end

"""
    physical_groups!(geo,groups)
Sets the array of physical groups in `geo`.
"""
function physical_groups! end

"""
    num_physical_groups(geo)
Return number of physical groups in `geo`.
"""
function num_physical_groups end
num_physical_groups(geo) = length(physical_groups(geo))

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

"""
    polytopal_complex(geo)
Return the polytopal complex associated with `geo`.
"""
function polytopal_complex end

struct MeshMetaData{Ti,T}
    physical_groups::Vector{PhysicalGroup{Ti}}
end

struct NodeInfo{D,T,Ti}
    node_coordinates::Vector{SVector{D,T}}
    hanging_nodes::HangingNodeVector{Ti,T}
    periodic_nodes::PeriodicNodeVector{Ti,T}
end

struct FaceInfo{Ti,A}
    face_nodes::Vector{JaggedArray{Ti,Int32}}
    face_faces::Matrix{JaggedArray{Ti,Int32}}
    face_ref_id::Vector{Vector{Ti}}
    ref_faces::A
end

struct SimpleMesh{D,T,Ti,A}
    node_coordinates::Vector{SVector{D,T}}
    cell_nodes::JaggedArray{Ti,Int32}
    cell_ref_id::{Vector{Ti}}
    ref_cell::A
    metadata::MeshMetaData{Ti,T}
    polytopal_complex::Ref{T}
end

struct FEMesh{D,T,Ti,A}
    nodes::NodeInfo
    faces::
end

struct PolytopalComplex
    node_coordinates::Vector{SVector{D,T}}
    face_nodes::Vector{JaggedArray{Ti,Int32}}
    face_faces::Matrix{JaggedArray{Ti,Int32}}
    face_ref_id::Vector{Vector{Ti}}
    ref_faces::A
    metadata::MeshMetaData{Ti,T}
    face_faces::Matrix{JaggedArray{Ti,Int32}}
end

femesh

femesh = generate_faces(femesh,1)



struct PolytopalComplex{D,T,Ti,A}
    node_coordinates::Vector{SVector{D,T}}
    face_nodes::Vector{JaggedArray{Ti,Int32}}
    face_faces::Matrix{JaggedArray{Ti,Int32}}
    face_ref_id::Vector{Vector{Ti}}
    ref_faces::A
    hanging_nodes::HangingNodeVector{Ti,T}
    periodic_nodes::PeriodicNodeVector{Ti,T}
    physical_groups::Vector{PhysicalGroup{Ti}}
end

function PolytopalComplex{D,T,Ti}(ref_faces) where {D,T,Ti}
    node_coordinates::Vector{SVector{D,T}}
    face_nodes::Vector{JaggedArray{Ti,Int32}}
    face_faces::Matrix{JaggedArray{Ti,Int32}}
    face_ref_id::Vector{Vector{Ti}}
    ref_faces::A
    hanging_nodes::HangingNodeVector{Ti,T}
    periodic_nodes::PeriodicNodeVector{Ti,T}
    physical_groups::Vector{PhysicalGroup{Ti}}

end


polytopal_complex(a::PolytopalComplex) = a
num_dims(a::PolytopalComplex) = length(a.face_to_id)-1
face_nodes(a::PolytopalComplex,dim) = a.face_nodes[dim+1]
face_faces(a::PolytopalComplex,m,n) = a.face_faces[m+1,n+1]
node_coordinates(a::PolytopalComplex) = a.node_coordinates
ref_faces(a::PolytopalComplex,::Val{dim}) = a.ref_faces[dim+1]
ref_faces(a::PolytopalComplex,dim) = a.ref_faces[dim+1]
face_ref_id(a::PolytopalComplex,dim) = a.face_ref_id[dim+1]
periodic_nodes(a::PolytopalComplex) = a.periodic_nodes
hanging_nodes(a::PolytopalComplex) = a.hanging_nodes
is_simplex(a::PolytopalComplex) = num_dims(a)+1 == num_vertices(a)
is_hypercube(a::PolytopalComplex) = 2^num_dims(a) == num_vertices(a)
physical_groups(a::PolytopalComplex) = a.physical_groups

function vtk_mesh_cell(a::PolytopalComplex)
    t = if num_dims(a) == 0
       WriteVTK.VTKCellTypes.VTK_VERTEX
    elseif num_dims(a) == 1
        if num_nodes(a) == 2
            WriteVTK.VTKCellTypes.VTK_LINE
        else
            error("Case not implemented")
        end
    elseif num_dims(a) == 2
        if num_nodes(a) == 3
            WriteVTK.VTKCellTypes.VTK_TRIANGLE
        else
            error("Case not implemented")
        end
    elseif num_dims(a) == 3
        if num_nodes(a) == 4
            WriteVTK.VTKCellTypes.VTK_TETRA
        else
            error("Case not implemented")
        end
    else
        error("Case not implemented")
    end
    nodes -> WriteVTK.MeshCell(t,nodes)
end


struct Polytope{B}
    boundary::B
end

polytopal_complex(a::Polytope) = a
num_dims(a::Polytope) = num_dims(a.boundary)+1
function face_nodes(a::Polytope,dim)
    if dim == num_dims(a)
        
    else
    end
end
face_faces(a::Polytope,m,n) = a.face_faces[m+1,n+1]
node_coordinates(a::Polytope) = a.node_coordinates
ref_faces(a::Polytope,::Val{dim}) = a.ref_faces[dim+1]
ref_faces(a::Polytope,dim) = a.ref_faces[dim+1]
face_ref_id(a::Polytope,dim) = a.face_ref_id[dim+1]
periodic_nodes(a::Polytope) = a.periodic_nodes
hanging_nodes(a::Polytope) = a.hanging_nodes
is_simplex(a::Polytope) = num_dims(a)+1 == num_vertices(a)
is_hypercube(a::Polytope) = 2^num_dims(a) == num_vertices(a)
physical_groups(a::Polytope) = a.physical_groups

#ref_point = Point()
#
#ref_line_boundary = fe_mesh(
#   node_coordinates=[(0),(1)],
#   cell_nodes=[[1],[2]]
#   ref_cell=ref_point)
#
#ref_line = Polytope(ref_line_boundary)
#
#ref_triangle_boundary = fe_mesh(
#    node_coordinates=[(0,0),(1,0),(0,1)],
#    face_nodes=[ [[1],[2],[3]], [[1,2],[2,3],[3,1]] ]
#    face_ref_id=[[1,1,1],[1,1,1]]
#    ref_faces=([ref_point],[ref_line]))
#
#ref_triangle = Polytope(ref_triangle_boundary)
#
#mesh = fe_mesh(
#  node_coordinates=[(0,0),(1,0),(0,1),(1,1)],
#  face_nodes=[ Vector{Int}[], [[1,2],[2,4]], [[1,2,3],[2,3,4]] ]
#  face_ref_id=[ Int[], [1,1],[1,1]],
#  ref_faces=( [],[ref_line],[ref_triangle]),
#  generate_faces = true)
#
#mesh, newids = generate_faces(mesh,1)
#mesh, newids = generate_faces(mesh,0)


#grid = CartesianGrid(3,3)
#mesh = fe_mesh(grid)
#mesh = fe_mesh(grid,isperiodic=(0,1))
#mesh = fe_mesh(grid,generate_faces=false)





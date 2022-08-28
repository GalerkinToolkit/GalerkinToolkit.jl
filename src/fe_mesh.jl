
"""
   const INVALID_ID::Int

Stores an invalid linear index for 1-based indexed objects.
`INVALID_ID` is some number less than `1`. The user should not rely on the
specific value.
"""
const INVALID_ID = -1

"""
    domain_dim(geo)
Return the number of space dimensions of the geometrical object `geo`.
"""
function domain_dim end

"""
    ambient_dim(geo)
Return the number of dimensions of the ambient space,
where the geometrical object `geo` is embedded.
"""
function ambient_dim end
ambient_dim(geo) = length(eltype(node_coordinates(geo)))

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
face_vertices(geo) = face_faces(geo,domain_dim(geo),0)

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
    info = periodic_nodes(geo)
Returns information about periodic nodes in `geo`.
- `info.periodic` is a vector with the node ids of periodic nodes.
- `info.master` is a vector containing the corresponding master ids.
- `info.coeff` is a vector containing the scaling a master node over
its periodic node (usually a vector of `1.0`).
"""
function periodic_nodes end

"""
    periodic_nodes!(geo,info)
Sets the information about periodic nodes in `geo`.
"""
function periodic_nodes! end

"""
    num_periodic_nodes(geo)
Number of periodic nodes in `geo`.
"""
function num_periodic_nodes end
num_periodic_nodes(geo) = length(periodic_nodes(geo).periodic)

"""
    has_periodic_nodes(geo)
`true` if `num_periodic_nodes(geo)>0`, `false` otherwise.
"""
function has_periodic_nodes end
has_periodic_nodes(geo) = num_periodic_nodes(geo)>0

"""
    info = hanging_nodes(geo)
Returns information about hanging nodes in `geo`.
- `info.hanging` is a vector with the node ids of hanging nodes.
- `info.masters` is a vector containing the corresponding master ids.
- `info.coeffs` is a vector containing the coefficients to compute
the value at a hanging node from its corresponding masters.

These vectors are of length `num_hanging_nodes(geo)`.
"""
function hanging_nodes end

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
num_hanging_nodes(geo) = length(hanging_nodes(geo).hanging)

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
Get a handle to the physical groups in `geo`.
"""
function physical_groups end

"""
    physical_groups!(geo,groups)
Sets the handle to the physical groups in `geo`.
"""
function physical_groups! end

"""
    physical_group_faces(groups,dim,label)

Return the face ids for faces of dimension `dim` with label `label`
in the physical groups `groups`.
The label can be either an integer group id or a string with the group name.
"""
function physical_group_faces end

"""
    physical_group_faces!(groups,faces,dim,label)

Set the face ids `faces` for faces of dimension `dim` with label `label`
in the physical groups `groups`.
The label can be either an integer group id or a string with the group name.
"""
function physical_group_faces! end

"""
    has_physical_group(groups,dim,label)
`true` if there is a group of dimension `dim` with label `label` in `groups`.
"""
function has_physical_group end

"""
    physical_group_name(groups,dim,label)
Return the name of the physical group of dimension `dim` and label `label`.
The label can be either an integer group id or a string with the group name.
In the latter case, this function will raise an exception if the given group
name is not found.
"""
function physical_group_name end

"""
    physical_group_id(groups,dim,label)
Return the id of the physical group of dimension `dim` and label `label`.
The label can be either an integer group id or a string with the group name.
In the former case, this function will raise an exception if the given group
id is not found.
"""
function physical_group_id end

"""
    physical_group_names(groups,dim)
Return a vector with the names of physical groups in `groups` of dimension `dim`.
"""
function physical_group_names end

"""
    physical_group_ids(groups,dim)
Return a vector with the ids of physical groups in `groups` of dimension `dim`.
"""
function physical_group_ids end

"""
    physical_group!(groups,dim,name[,id])
Add a new empty physical group to `groups` of dimension
`dim` name `name`  and id `id`. If `id` is not provided, an id will be chosen
automatically using `new_physical_group_id(groups)`.
"""
function physical_group! end

"""
    new_physical_group_id(groups)
Find a new available id in `groups`.
"""
function new_physical_group_id end

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

#"""
#    physical_group_nodes(
#        mesh,dim,label;
#        boundary=true,groups=physical_groups(mesh))
#Return the node ids in the physical group of dimension `dim` and
#label `label` in the mesh `mesh`. If `boundary==true`, include also
#nodes on the boundary (default) otherwise take only nodes on the interior
#of the physical group. If `groups`  is omitted, one takes the physical
#groups stored in the given mesh (default).
#"""
#function physical_group_nodes(mesh,dim,label;boundary=true,groups=physical_groups(mesh))
#  id = physical_group_id(groups,dim,label)
#  if boundary
#    faces_in_physical_group = physical_group_faces(groups,dim,id)
#    face_to_nodes = face_nodes(mesh,dim)
#    nnodes = num_nodes(mesh)
#    _physical_group_ids_with_boundary(nnodes,face_to_nodes,faces_in_physical_group)
#  else
#    error("Not implemented.")
#    # 1. Compute nodes with boundary
#    # 2. Identify d-1 faces on the boundary of the group
#    # 3. Collect nodes from faces computed in 2.
#    # 4. Subtract from nodes computed in 1.
#    # Setp 2. would require to compute the polytopal complex
#    #polycom = polytopal_complex(mesh)
#  end
#end
#
#function _physical_group_ids_with_boundary(nnodes,face_to_nodes,faces_in_physical_group)
#  node_to_mask = fill(false,nnodes)
#  for face in faces_in_physical_group
#    nodes = face_to_nodes[face]
#    for node in nodes
#      node_to_mask[node] = true
#    end
#  end
#  collect(Int32,findall(node_to_mask))
#end
#
#"""
#    physical_group_low_dim_faces(
#        mesh,dim,n,label;
#        boundary=true,groups=physical_groups(mesh))
#Return the low dimensional `n`-face ids (`n<dim`) in the physical group of dimension `dim` and
#label `label` in the mesh `mesh`. If `boundary==true`, include also
#faces on the boundary (default) otherwise take only nodes on the interior
#of the physical group. If `groups`  is omitted, one takes the physical
#groups stored in the given mesh (default).
#"""
#function physical_group_low_dim_faces(mesh,dim,n,label;boundary=true,groups=physical_groups(mesh))
#  id = physical_group_id(groups,dim,label)
#  if boundary
#    faces_in_physical_group = physical_group_faces(groups,dim,id)
#    face_to_nfaces = face_faces(mesh,dim,n)
#    nnfaces = num_faces(mesh,n)
#    _physical_group_ids_with_boundary(nnfaces,face_to_nfaces,faces_in_physical_group)
#  else
#    error("Not implemented.")
#    # 1. Compute nodes with boundary
#    # 2. Identify d-1 faces on the boundary of the group
#    # 3. Collect nodes from faces computed in 2.
#    # 4. Subtract from nodes computed in 1.
#    # Setp 2. would require to compute the polytopal complex
#    #polycom = polytopal_complex(mesh)
#  end
#end
#
#function classify_nodes(mesh,ids;boundary=fill(true,length(ids)))
#  D = domain_dim(mesh)
#  groups = physical_groups(mesh)
#  i_to_id = ids
#  ni = length(i_to_id)
#  nnodes = num_nodes(mesh)
#  node_to_i = fill(Int32(INVALID_ID),nnodes)
#  for d = D:-1:0
#    for i in 1:ni
#      id = i_to_id[i]
#      if !has_physical_group(groups,d,id)
#        continue
#      end
#      nodes = physical_group_nodes(mesh,d,id;boundary=boundary[i])
#      node_to_i[nodes] .= i
#    end
#  end
#  node_to_i
#end

struct PeriodicNodeCollection{Ti,T}
  periodic::Vector{Ti}
  master::Vector{Ti}
  coeff::Vector{T}
end

function PeriodicNodeCollection{Ti,T}() where {Ti,T}
  PeriodicNodeCollection(Ti[],Ti[],T[])
end

struct HangingNodeCollection{Ti,T}
  hanging::Vector{Ti}
  masters::JaggedArray{Ti,Int32}
  coeffs::JaggedArray{T,Int32}
end

function HangingNodeCollection{Ti,T}() where {Ti,T}
  HangingNodeCollection(
    Ti[],
    JaggedArray([Ti[]]),
    JaggedArray([T[]]))
end

struct PhysicalGroupCollection{Ti}
  physical_group_name::Vector{Dict{Int,String}}
  physical_group_id::Vector{Dict{String,Int}}
  physical_group_faces::Vector{Dict{Int,Vector{Ti}}}
end

function PhysicalGroupCollection(dim::Integer)
  PhysicalGroupCollection{Int32}(dim)
end

function PhysicalGroupCollection{Ti}(dim::Integer) where Ti
  PhysicalGroupCollection(
    [Dict{Int,String}() for d in 0:dim],
    [Dict{String,Int}() for d in 0:dim],
    [Dict{Int,Vector{Ti}}() for d in 0:dim])
end

physical_groups(a::PhysicalGroupCollection) = a
domain_dim(a::PhysicalGroupCollection) = length(a.physical_group_id)-1

function physical_group_faces(g::PhysicalGroupCollection,dim,label)
  g.physical_group_faces[dim+1][physical_group_id(g,dim,label)]
end

function physical_group_faces!(g::PhysicalGroupCollection,faces,dim,label)
  (g.physical_group_faces[dim+1][physical_group_id(g,dim,label)] = faces)
end

function has_physical_group(g::PhysicalGroupCollection,dim,name::AbstractString)
  haskey(g.physical_group_id[dim+1],name)
end

function has_physical_group(g::PhysicalGroupCollection,dim,id::Integer)
  haskey(g.physical_group_name[dim+1],id)
end

function physical_group_name(g::PhysicalGroupCollection,dim,id::Integer)
  g.physical_group_name[dim+1][Int(id)]
end

function physical_group_name(g::PhysicalGroupCollection,dim,name::AbstractString)
  @assert has_physical_group(g,dim,name)
  name
end

function physical_group_id(g::PhysicalGroupCollection,dim,name::AbstractString)
  g.physical_group_id[dim+1][String(name)]
end

function physical_group_id(g::PhysicalGroupCollection,dim,id::Integer)
  @assert has_physical_group(g,dim,id)
  id
end

function physical_group_names(g::PhysicalGroupCollection,dim)
  ids = physical_group_ids(g,dim)
  [ g.physical_group_name[dim+1][i] for i in ids]
end

function physical_group_ids(g::PhysicalGroupCollection,dim)
  sort(collect(keys(g.physical_group_name[dim+1])))
end

function physical_group!(g::PhysicalGroupCollection,dim,name,id=new_physical_group_id(g))
  haskey(g.physical_group_name[dim+1],id) && error(
    "id $id already present in PhysicalGroupCollection for dimension $dim")
  haskey(g.physical_group_id[dim+1],name) && error(
    "Name $name already present in PhysicalGroupCollection for dimension $dim")
  g.physical_group_name[dim+1][id] = name
  g.physical_group_id[dim+1][name] = id
  g.physical_group_faces[dim+1][id] = Int32[]
end

function new_physical_group_id(g::PhysicalGroupCollection)
  function maxid(dict)
    if length(dict) == 0
      id = 0
    else
      id = maximum(values(dict))
    end
  end
  id = maximum(map(maxid,g.physical_group_id))+1
  id
end

"""
    fe_mesh(geo)
Return the FE mesh associated with geo.
"""
function fe_mesh end

# We can add this one if we need generic
# fields
#mutable struct GenericFEMesh{}
#end

mutable struct FEMesh{N,Ti,Tf}
  node_coordinates::Vector{SVector{N,Tf}}
  face_nodes::Vector{JaggedArray{Ti,Int32}}
  face_ref_id::Vector{Vector{Int8}}
  ref_faces::Vector{Vector{Any}}
  periodic_nodes::PeriodicNodeCollection{Ti,Tf}
  hanging_nodes::HangingNodeCollection{Ti,Tf}
  physical_groups::PhysicalGroupCollection{Ti}
  buffer::Dict{Symbol,Any}
end
fe_mesh(a::FEMesh) = a

function FEMesh{N}(dim) where N
  FEMesh{N,Int32,Float64}(dim)
end

function FEMesh{N,Ti,Tf}(dim) where {N,Ti,Tf}
  FEMesh(
    zeros(SVector{N,Tf},0),
    [ JaggedArray([Ti[]]) for i in 1:(dim+1)],
    [ Int8[] for i in 1:(dim+1)],
    [ Any[] for i in 1:(dim+1)],
    PeriodicNodeCollection{Ti,Tf}(),
    HangingNodeCollection{Ti,Tf}(),
    PhysicalGroupCollection{Ti}(dim),
    Dict{Symbol,Any}()
   )
end

domain_dim(m::FEMesh) = length(m.face_ref_id)-1
node_coordinates(m::FEMesh) = m.node_coordinates
face_nodes(m::FEMesh,dim) = m.face_nodes[dim+1]
face_ref_id(m::FEMesh,dim) = m.face_ref_id[dim+1]
ref_faces(m::FEMesh,dim) = m.ref_faces[dim+1]
periodic_nodes(m::FEMesh) = m.periodic_nodes
hanging_nodes(m::FEMesh) = m.hanging_nodes
physical_groups(m::FEMesh) = m.physical_groups

node_coordinates!(m::FEMesh,v) = (m.node_coordinates = v)
face_nodes!(m::FEMesh,v,dim) = (m.face_nodes[dim+1] = v)
face_ref_id!(m::FEMesh,v,dim) = (m.face_ref_id[dim+1] = v)
ref_faces!(m::FEMesh,v,dim) = (m.ref_faces[dim+1] = v)
periodic_nodes!(m::FEMesh,v) = (m.periodic_nodes = v)
hanging_nodes!(m::FEMesh,v) = (m.hanging_nodes = v)
physical_groups!(m::FEMesh,v) = (m.physical_groups = v)

function polytopal_complex(m::FEMesh)
  if !haskey(m.buffer,:polytopal_complex)
    polycomplex = PolytopalComplex(m)
    m.buffer[:polytopal_complex] = polycomplex
  end
  m.buffer[:polytopal_complex]
end

function node_vertex(mesh::FEMesh)
  if !haskey(mesh.buffer,:node_vertex)
    _setup_vertices!(mesh)
  end
  mesh.buffer[:node_vertex]
end

function vertex_node(mesh::FEMesh)
  if !haskey(mesh.buffer,:vertex_node)
    _setup_vertices!(mesh)
  end
  mesh.buffer[:vertex_node]
end

function _setup_vertices!(mesh)
  D = domain_dim(mesh)
  nnodes = num_nodes(mesh)
  d_refid_to_refface = [ref_faces(mesh,d) for d in 0:D]
  d_refid_to_lvertex_to_lnodes = [ref_face_nodes(mesh,d,0) for d in 0:D]
  d_dface_to_refid = [ face_ref_id(mesh,d) for d in 0:D]
  d_dface_to_nodes = [ face_nodes(mesh,d) for d in 0:D]
  node_to_vertex, vertex_to_node = _setup_vertices(
    D,
    nnodes,
    d_refid_to_refface,
    d_refid_to_lvertex_to_lnodes,
    d_dface_to_refid,
    d_dface_to_nodes)
  mesh.buffer[:node_vertex] = node_to_vertex
  mesh.buffer[:vertex_node] = vertex_to_node
end

function _setup_vertices(
  D,
  nnodes,
  d_refid_to_refface,
  d_refid_to_lvertex_to_lnodes,
  d_dface_to_refid,
  d_dface_to_nodes)

  node_to_touched = fill(false,nnodes)
  for d in D:-1:0
    refid_to_refface = d_refid_to_refface[d+1]
    refid_to_lvertex_to_lnodes = d_refid_to_lvertex_to_lnodes[d+1]
    dface_to_refid = d_dface_to_refid[d+1]
    dface_to_nodes = d_dface_to_nodes[d+1]
    for dface in 1:length(dface_to_refid)
      refid = dface_to_refid[dface]
      nodes = dface_to_nodes[dface]
      lvertex_to_lnodes = refid_to_lvertex_to_lnodes[refid]
      for lnodes in lvertex_to_lnodes
        lnode = first(lnodes)
        node = nodes[lnode]
        node_to_touched[node] = true
      end
    end
  end
  vertex_to_node = collect(Int32,findall(node_to_touched))
  nvertices = length(vertex_to_node)
  node_to_vertex = fill(Int32(INVALID_ID),nnodes)
  node_to_vertex[vertex_to_node] .= Int32(1):Int32(nvertices)
  node_to_vertex, vertex_to_node
end

function face_faces(mesh::FEMesh,dim_from,dim_to)
  @assert dim_to == 0 "Case not implemented"
  if !haskey(mesh.buffer,:face_vertices)
    J = typeof(GenericJaggedArray(Vector{Int32}[]))
    mesh.buffer[:face_vertices] = Vector{J}(undef,domain_dim(mesh)+1)
  end
  if !isassigned(mesh.buffer[:face_vertices],dim_from+1)
    _setup_mesh_face_vertices!(mesh,dim_from)
  end
  mesh.buffer[:face_vertices][dim_from+1]
end

function _setup_mesh_face_vertices!(poly,d)
  mesh = poly.mesh
  node_to_vertex = node_vertex(poly)
  refid_to_lvertex_to_lnodes = ref_face_nodes(mesh,d,0)
  dface_to_refid = face_ref_id(mesh,d)
  dface_to_nodes = face_nodes(mesh,d)
  dface_to_vertices = _setup_mesh_face_vertices(
    node_to_vertex,
    refid_to_lvertex_to_lnodes,
    dface_to_refid,
    dface_to_nodes)
  poly.buffer[:mesh_face_vertices][d+1] = dface_to_vertices
end

function _setup_mesh_face_vertices(
  node_to_vertex,
  refid_to_lvertex_to_lnodes,
  dface_to_refid,
  dface_to_nodes)

  ptrs = zeros(Int32,length(dface_to_refid)+1)
  for dface in 1:length(dface_to_refid)
    refid = dface_to_refid[dface]
    lvertex_to_lnodes = refid_to_lvertex_to_lnodes[refid]
    ptrs[dface+1] += length(lvertex_to_lnodes)
  end
  prefix!(ptrs)
  ndata = ptrs[end]-1
  data = zeros(Int32,ndata)
  for dface in 1:length(dface_to_refid)
    refid = dface_to_refid[dface]
    nodes = dface_to_nodes[dface]
    lvertex_to_lnodes = refid_to_lvertex_to_lnodes[refid]
    p = ptrs[dface]-Int32(1)
    for (lvertex,lnodes) in enumerate(lvertex_to_lnodes)
      lnode = first(lnodes)
      node = nodes[lnode]
      vertex = node_to_vertex[node]
      data[p+lvertex] = vertex
    end
  end
  dface_to_vertices = GenericJaggedArray(data,ptrs)
  dface_to_vertices
end


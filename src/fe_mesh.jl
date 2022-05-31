
const INVALID = -1

function domain_dim end
function ambient_dim end
function face_ref_id end
function ref_faces end
function face_nodes end
function face_own_nodes end
function node_coordinates end
function periodic_nodes end
function hanging_nodes end
function is_hanging end
function is_periodic end
function num_faces end
function num_nodes end
function is_simplex end
function is_hypercube end
function group_faces end
function group_name end
function group_dim end
function group_dims end
function group_names end
function group_id end
function group_ids end
function add_group! end
function group_faces! end
function face_incidence end
function face_vertices end

function ref_face_incidence(mesh,rank,m,n)
  refid_to_refface = ref_faces(mesh,rank)
  map(rf->face_incidence(rf,m,n),refid_to_refface)
end
function ref_face_nodes(mesh,rank,m)
  refid_to_refface = ref_faces(mesh,rank)
  map(rf->face_nodes(rf,m),refid_to_refface)
end
function ref_face_own_nodes(mesh,rank,m)
  refid_to_refface = ref_faces(mesh,rank)
  map(rf->face_own_nodes(rf,m),refid_to_refface)
end

default_is_simplex(geo) = false
default_is_hypercube(geo) = false
default_is_hanging(geo) =  length(first(hanging_nodes(geo)))>0
default_hanging_nodes(geo) = (Int32[],JaggedArray(Vector{Int32}[]))
default_is_periodic(geo) = length(first(periodic_nodes(geo)))>0
default_periodic_nodes(geo) = (Int32[],Int32[])
default_num_faces(geo,rank) = length(face_ref_id(geo,rank))
default_num_nodes(geo) = length(node_coordinates(geo))
default_ambient_dim(geo) = length(eltype(node_coordinates(geo)))
default_physical_groups(geo) = GroupCollection(VOID,domain_dim(geo))

is_simplex(geo) = default_is_simplex(geo)
is_hypercube(geo) = default_is_hypercube(geo)
is_hanging(geo) = default_is_hanging(geo)
hanging_nodes(geo) = default_hanging_nodes(geo)
is_periodic(geo) = default_is_periodic(geo)
periodic_nodes(geo) = default_periodic_nodes(geo)
num_faces(geo,rank) = default_num_faces(geo,rank)
num_nodes(geo) = default_num_nodes(geo)
ambient_dim(geo) = default_ambient_dim(geo)

struct EmptyInitializer end

const VOID = EmptyInitializer()

# Allow same ids in different dims? (not possible in gmsh)
# TO-think about the union of groups in different dims (not possible in gmsh)
# Allow same name for groups in different dims? (possible in gmsh!)
struct GroupCollection
  domain_dim::Int
  group_dim::Dict{Int,Int}
  group_names::Dict{Int,String}
  group_ids::Dict{String,Int}
  group_faces::Dict{Int,Vector{Int32}}
end
physical_groups(a::GroupCollection) = a

function GroupCollection(::EmptyInitializer,dim::Integer)
  GroupCollection(
    dim,
    Dict{Int,Int}(),
    Dict{Int,String}(),
    Dict{String,Int}(),
    Dict{Int,Vector{Int32}}())
end

function add_group!(g::GroupCollection,dim,name,id)
  haskey(g.group_names,id) && error("id $id already present in GroupCollection")
  haskey(g.group_ids,name) && error("Name $name alreade present in GroupCollection")
  g.group_dim[id] = dim
  g.group_names[id] = name
  g.group_ids[name] = id
  g.group_faces[id] = Int32[]
end

function add_group!(g::GroupCollection,dim,name)
  if length(g.group_ids) == 0
    id = 1
  else
    id = maximum(values(g.group_ids)) + 1
  end
  add_group!(g,dim,name,id)
end

group_faces(g::GroupCollection,name) = g.group_faces[group_id(g,name)]
group_faces!(g::GroupCollection,faces,name) = (g.group_faces[group_id(g,name)] = faces)
domain_dim(g::GroupCollection) = g.domain_dim
group_dim(g::GroupCollection,id) = g.group_dim[group_id(g,id)]
group_id(g::GroupCollection,name::String) = g.group_ids[name]
group_id(g::GroupCollection,id::Int) = id
group_name(g::GroupCollection,name::String) = name
group_name(g::GroupCollection,id::Int) = g.group_names[id]
function group_names(g::GroupCollection)
  ids = group_ids(g)
  [ g.group_names[i] for i in ids]
end
function group_ids(g::GroupCollection)
  pairs = sort(collect(g.group_names))
  map(first,pairs)
end
function group_dims(g::GroupCollection)
  ids = group_ids(g)
  [ g.group_dim[i] for i in ids]
end
function group_names(g::GroupCollection,d)
  ids = group_ids(g,d)
  [ g.group_names[i] for i in ids]
end
function group_ids(g::GroupCollection,d)
  pairs = sort(collect(g.group_names))
  ids = map(first,pairs)
  filter(i->g.group_dim[i]==d,ids)
end
function group_dims(g::GroupCollection,d)
  ids = group_ids(g,d)
  [ g.group_dim[i] for i in ids]
end

function classify_nodes(mesh,ids;boundary=fill(true,length(ids)))
  D = ambient_dim(mesh)
  groups = physical_groups(mesh)
  i_to_id = ids
  ni = length(i_to_id)
  nnodes = num_nodes(mesh)
  node_to_i = fill(Int32(INVALID),nnodes)
  for d = D:-1:0
    for i in 1:ni
      id = i_to_id[i]
      if group_dim(groups,id) != d
        continue
      end
      nodes = group_nodes(mesh,id;boundary=boundary[i])
      node_to_i[nodes] .= i
    end
  end
  node_to_i
end

function group_nodes(mesh,id;boundary=true)
  if boundary
    group_nodes_with_boundary(mesh,id)
  else
    group_nodes_without_boundary(mesh,id)
  end
end

function group_nodes_with_boundary(mesh,id)
  function _barrier(nnodes,face_to_nodes,faces_in_group)
    node_to_mask = fill(false,nnodes)
    for face in faces_in_group
      nodes = face_to_nodes[face]
      for node in nodes
        node_to_mask[node] = true
      end
    end
    collect(Int32,findall(node_to_mask))
  end
  groups = physical_groups(mesh)
  faces_in_group = group_faces(groups,id)
  d = group_dim(groups,id)
  face_to_nodes = face_nodes(mesh,d)
  nnodes = num_nodes(mesh)
  _barrier(nnodes,face_to_nodes,faces_in_group)
end

function group_nodes_without_boundary(mesh,id)
  error("Not implemented.")
  # 1. Compute nodes with boundary
  # 2. Identify d-1 faces on the boundary of the group
  # 3. Collect nodes from faces computed in 2.
  # 4. Subtract from nodes computed in 1.
  # Setp 2. would require to compute the polytopal complex
  #polycom = polytopal_complex(mesh)
end

function fe_mesh end

mutable struct SimpleFEMesh{T}
  domain_dim::Int
  node_coordinates::Vector{T}
  face_nodes::Vector{JArray{Int32}}
  face_ref_id::Vector{Vector{Int8}}
  ref_faces::Vector{Vector{Any}}
  periodic_nodes::Tuple{Vector{Int32},Vector{Int32}}
  hanging_nodes::Tuple{Vector{Int32},JArray{Int32}}
  physical_groups::GroupCollection
  buffer::Dict{Symbol,Any}
end
fe_mesh(a::SimpleFEMesh) = a

function SimpleFEMesh{T}(::EmptyInitializer,rank) where T
  SimpleFEMesh(
    rank,
    zeros(T,0),
    [ JaggedArray(Vector{Int32}[]) for i in 1:(rank+1)],
    [ Int8[] for i in 1:(rank+1)],
    [ Any[] for i in 1:(rank+1)],
    (Int32[],Int32[]),
    (Int32[],JaggedArray(Vector{Int32}[])),
    GroupCollection(VOID,rank),
    Dict{Symbol,Any}()
   )
end

domain_dim(m::SimpleFEMesh) = m.domain_dim
node_coordinates(m::SimpleFEMesh) = m.node_coordinates
face_nodes(m::SimpleFEMesh,rank) = m.face_nodes[rank+1]
face_ref_id(m::SimpleFEMesh,rank) = m.face_ref_id[rank+1]
ref_faces(m::SimpleFEMesh,rank) = m.ref_faces[rank+1]
periodic_nodes(m::SimpleFEMesh) = m.periodic_nodes
hanging_nodes(m::SimpleFEMesh) = m.periodic_nodes
physical_groups(m::SimpleFEMesh) = m.physical_groups

node_coordinates!(m::SimpleFEMesh,v) = (m.node_coordinates = v)
face_nodes!(m::SimpleFEMesh,rank,v) = (m.face_nodes[rank+1] = v)
face_ref_id!(m::SimpleFEMesh,rank,v) = (m.face_ref_id[rank+1] = v)
ref_faces!(m::SimpleFEMesh,rank,v) = (m.ref_faces[rank+1] = v)
periodic_nodes!(m::SimpleFEMesh,v) = (m.periodic_nodes = v)
hanging_nodes!(m::SimpleFEMesh,v) = (m.periodic_nodes = v)
physical_groups!(m::SimpleFEMesh,v) = (m.physical_groups = v)


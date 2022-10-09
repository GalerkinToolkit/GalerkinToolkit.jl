
















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
#  D = num_dims(mesh)
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

struct PeriodicNodeVector{Ti,T}
  periodic::Vector{Ti}
  master::Vector{Ti}
  coeff::Vector{T}
end

function PeriodicNodeVector{Ti,T}() where {Ti,T}
  PeriodicNodeVector(Ti[],Ti[],T[])
end

struct HangingNodeVector{Ti,T}
  hanging::Vector{Ti}
  masters::JaggedArray{Ti,Int32}
  coeffs::JaggedArray{T,Int32}
end

function HangingNodeVector{Ti,T}() where {Ti,T}
  HangingNodeVector(
    Ti[],
    JaggedArray([Ti[]]),
    JaggedArray([T[]]))
end

#struct PhysicalGroupCollection{Ti}
#  physical_group_name::Vector{Dict{Int,String}}
#  physical_group_id::Vector{Dict{String,Int}}
#  physical_group_faces::Vector{Dict{Int,Vector{Ti}}}
#end
#
#function PhysicalGroupCollection(dim::Integer)
#  PhysicalGroupCollection{Int32}(dim)
#end
#
#function PhysicalGroupCollection{Ti}(dim::Integer) where Ti
#  PhysicalGroupCollection(
#    [Dict{Int,String}() for d in 0:dim],
#    [Dict{String,Int}() for d in 0:dim],
#    [Dict{Int,Vector{Ti}}() for d in 0:dim])
#end
#
#physical_groups(a::PhysicalGroupCollection) = a
#num_dims(a::PhysicalGroupCollection) = length(a.physical_group_id)-1
#
#function physical_group_faces(g::PhysicalGroupCollection,dim,label)
#  g.physical_group_faces[dim+1][physical_group_id(g,dim,label)]
#end
#
#function physical_group_faces!(g::PhysicalGroupCollection,faces,dim,label)
#  (g.physical_group_faces[dim+1][physical_group_id(g,dim,label)] = faces)
#end
#
#function has_physical_group(g::PhysicalGroupCollection,dim,name::AbstractString)
#  haskey(g.physical_group_id[dim+1],name)
#end
#
#function has_physical_group(g::PhysicalGroupCollection,dim,id::Integer)
#  haskey(g.physical_group_name[dim+1],id)
#end
#
#function physical_group_name(g::PhysicalGroupCollection,dim,id::Integer)
#  g.physical_group_name[dim+1][Int(id)]
#end
#
#function physical_group_name(g::PhysicalGroupCollection,dim,name::AbstractString)
#  @assert has_physical_group(g,dim,name)
#  name
#end
#
#function physical_group_id(g::PhysicalGroupCollection,dim,name::AbstractString)
#  g.physical_group_id[dim+1][String(name)]
#end
#
#function physical_group_id(g::PhysicalGroupCollection,dim,id::Integer)
#  @assert has_physical_group(g,dim,id)
#  id
#end
#
#function physical_group_names(g::PhysicalGroupCollection,dim)
#  ids = physical_group_ids(g,dim)
#  [ g.physical_group_name[dim+1][i] for i in ids]
#end
#
#function physical_group_ids(g::PhysicalGroupCollection,dim)
#  sort(collect(keys(g.physical_group_name[dim+1])))
#end
#
#function physical_group!(g::PhysicalGroupCollection,dim,name,id=next_physical_group_id(g))
#  haskey(g.physical_group_name[dim+1],id) && error(
#    "id $id already present in PhysicalGroupCollection for dimension $dim")
#  haskey(g.physical_group_id[dim+1],name) && error(
#    "Name $name already present in PhysicalGroupCollection for dimension $dim")
#  g.physical_group_name[dim+1][id] = name
#  g.physical_group_id[dim+1][name] = id
#  g.physical_group_faces[dim+1][id] = Int32[]
#end
#
#function next_physical_group_id(g::PhysicalGroupCollection)
#  function maxid(dict)
#    if length(dict) == 0
#      id = 0
#    else
#      id = maximum(values(dict))
#    end
#  end
#  id = maximum(map(maxid,g.physical_group_id))+1
#  id
#end
#
#"""
#    fe_mesh(geo)
#Return the FE mesh associated with geo.
#"""
#function fe_mesh end
#
## We can add this one if we need generic
## fields
##mutable struct GenericFEMesh{}
##end
#
#mutable struct FEMesh{N,Ti,Tf}
#  node_coordinates::Vector{SVector{N,Tf}}
#  face_nodes::Vector{JaggedArray{Ti,Int32}}
#  face_ref_id::Vector{Vector{Int8}}
#  ref_faces::Vector{Vector{Any}}
#  periodic_nodes::PeriodicNodeVector{Ti,Tf}
#  hanging_nodes::HangingNodeVector{Ti,Tf}
#  physical_groups::PhysicalGroupCollection{Ti}
#  buffer::Dict{Symbol,Any}
#end
#fe_mesh(a::FEMesh) = a
#
#function FEMesh{N}(dim) where N
#  FEMesh{N,Int32,Float64}(dim)
#end
#
#function FEMesh{N,Ti,Tf}(dim) where {N,Ti,Tf}
#  FEMesh(
#    zeros(SVector{N,Tf},0),
#    [ JaggedArray([Ti[]]) for i in 1:(dim+1)],
#    [ Int8[] for i in 1:(dim+1)],
#    [ Any[] for i in 1:(dim+1)],
#    PeriodicNodeVector{Ti,Tf}(),
#    HangingNodeVector{Ti,Tf}(),
#    PhysicalGroupCollection{Ti}(dim),
#    Dict{Symbol,Any}()
#   )
#end
#
#num_dims(m::FEMesh) = length(m.face_ref_id)-1
#node_coordinates(m::FEMesh) = m.node_coordinates
#face_nodes(m::FEMesh,dim) = m.face_nodes[dim+1]
#face_ref_id(m::FEMesh,dim) = m.face_ref_id[dim+1]
#ref_faces(m::FEMesh,dim) = m.ref_faces[dim+1]
#periodic_nodes(m::FEMesh) = m.periodic_nodes
#hanging_nodes(m::FEMesh) = m.hanging_nodes
#physical_groups(m::FEMesh) = m.physical_groups
#
#node_coordinates!(m::FEMesh,v) = (m.node_coordinates = v)
#face_nodes!(m::FEMesh,v,dim) = (m.face_nodes[dim+1] = v)
#face_ref_id!(m::FEMesh,v,dim) = (m.face_ref_id[dim+1] = v)
#ref_faces!(m::FEMesh,v,dim) = (m.ref_faces[dim+1] = v)
#periodic_nodes!(m::FEMesh,v) = (m.periodic_nodes = v)
#hanging_nodes!(m::FEMesh,v) = (m.hanging_nodes = v)
#physical_groups!(m::FEMesh,v) = (m.physical_groups = v)
#
#function polytopal_complex(m::FEMesh)
#  if !haskey(m.buffer,:polytopal_complex)
#    polycomplex = PolytopalComplex(m)
#    m.buffer[:polytopal_complex] = polycomplex
#  end
#  m.buffer[:polytopal_complex]
#end
#
#function node_vertex(mesh::FEMesh)
#  if !haskey(mesh.buffer,:node_vertex)
#    _setup_vertices!(mesh)
#  end
#  mesh.buffer[:node_vertex]
#end
#
#function vertex_node(mesh::FEMesh)
#  if !haskey(mesh.buffer,:vertex_node)
#    _setup_vertices!(mesh)
#  end
#  mesh.buffer[:vertex_node]
#end
#
#function _setup_vertices!(mesh)
#  D = num_dims(mesh)
#  nnodes = num_nodes(mesh)
#  d_refid_to_refface = [ref_faces(mesh,d) for d in 0:D]
#  d_refid_to_lvertex_to_lnodes = [ref_face_nodes(mesh,d,0) for d in 0:D]
#  d_dface_to_refid = [ face_ref_id(mesh,d) for d in 0:D]
#  d_dface_to_nodes = [ face_nodes(mesh,d) for d in 0:D]
#  node_to_vertex, vertex_to_node = _setup_vertices(
#    D,
#    nnodes,
#    d_refid_to_refface,
#    d_refid_to_lvertex_to_lnodes,
#    d_dface_to_refid,
#    d_dface_to_nodes)
#  mesh.buffer[:node_vertex] = node_to_vertex
#  mesh.buffer[:vertex_node] = vertex_to_node
#end
#
#function _setup_vertices(
#  D,
#  nnodes,
#  d_refid_to_refface,
#  d_refid_to_lvertex_to_lnodes,
#  d_dface_to_refid,
#  d_dface_to_nodes)
#
#  node_to_touched = fill(false,nnodes)
#  for d in D:-1:0
#    refid_to_refface = d_refid_to_refface[d+1]
#    refid_to_lvertex_to_lnodes = d_refid_to_lvertex_to_lnodes[d+1]
#    dface_to_refid = d_dface_to_refid[d+1]
#    dface_to_nodes = d_dface_to_nodes[d+1]
#    for dface in 1:length(dface_to_refid)
#      refid = dface_to_refid[dface]
#      nodes = dface_to_nodes[dface]
#      lvertex_to_lnodes = refid_to_lvertex_to_lnodes[refid]
#      for lnodes in lvertex_to_lnodes
#        lnode = first(lnodes)
#        node = nodes[lnode]
#        node_to_touched[node] = true
#      end
#    end
#  end
#  vertex_to_node = collect(Int32,findall(node_to_touched))
#  nvertices = length(vertex_to_node)
#  node_to_vertex = fill(Int32(INVALID_ID),nnodes)
#  node_to_vertex[vertex_to_node] .= Int32(1):Int32(nvertices)
#  node_to_vertex, vertex_to_node
#end
#
#function face_faces(mesh::FEMesh,dim_from,dim_to)
#  @assert dim_to == 0 "Case not implemented"
#  if !haskey(mesh.buffer,:face_vertices)
#    J = typeof(GenericJaggedArray(Vector{Int32}[]))
#    mesh.buffer[:face_vertices] = Vector{J}(undef,num_dims(mesh)+1)
#  end
#  if !isassigned(mesh.buffer[:face_vertices],dim_from+1)
#    _setup_mesh_face_vertices!(mesh,dim_from)
#  end
#  mesh.buffer[:face_vertices][dim_from+1]
#end
#
#function _setup_mesh_face_vertices!(poly,d)
#  mesh = poly.mesh
#  node_to_vertex = node_vertex(poly)
#  refid_to_lvertex_to_lnodes = ref_face_nodes(mesh,d,0)
#  dface_to_refid = face_ref_id(mesh,d)
#  dface_to_nodes = face_nodes(mesh,d)
#  dface_to_vertices = _setup_mesh_face_vertices(
#    node_to_vertex,
#    refid_to_lvertex_to_lnodes,
#    dface_to_refid,
#    dface_to_nodes)
#  poly.buffer[:mesh_face_vertices][d+1] = dface_to_vertices
#end
#
#function _setup_mesh_face_vertices(
#  node_to_vertex,
#  refid_to_lvertex_to_lnodes,
#  dface_to_refid,
#  dface_to_nodes)
#
#  ptrs = zeros(Int32,length(dface_to_refid)+1)
#  for dface in 1:length(dface_to_refid)
#    refid = dface_to_refid[dface]
#    lvertex_to_lnodes = refid_to_lvertex_to_lnodes[refid]
#    ptrs[dface+1] += length(lvertex_to_lnodes)
#  end
#  prefix!(ptrs)
#  ndata = ptrs[end]-1
#  data = zeros(Int32,ndata)
#  for dface in 1:length(dface_to_refid)
#    refid = dface_to_refid[dface]
#    nodes = dface_to_nodes[dface]
#    lvertex_to_lnodes = refid_to_lvertex_to_lnodes[refid]
#    p = ptrs[dface]-Int32(1)
#    for (lvertex,lnodes) in enumerate(lvertex_to_lnodes)
#      lnode = first(lnodes)
#      node = nodes[lnode]
#      vertex = node_to_vertex[node]
#      data[p+lvertex] = vertex
#    end
#  end
#  dface_to_vertices = GenericJaggedArray(data,ptrs)
#  dface_to_vertices
#end


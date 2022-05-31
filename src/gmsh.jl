
function default_gmsh_options()
  [
    "General.Terminal"=>1,
    "Mesh.SaveAll"=>1,
    "Mesh.MedImportGroupsOfNodes"=>1
  ]
end

function with_gmsh(f;options=default_gmsh_options(),renumber=true)
  gmsh.initialize()
  for (k,v) in options
    gmsh.option.setNumber(k,v)
  end
  renumber && gmsh.model.mesh.renumberNodes()
  renumber && gmsh.model.mesh.renumberElements()
  try
    return f()
  finally
    gmsh.finalize()
  end
end

struct MshFile{T}
  path::String
  kwargs::T
end

msh_file(args...;kwargs...) = MshFile(joinpath(args...),kwargs)

function fe_mesh(f::MshFile)
  with_gmsh(;f.kwargs...) do
    gmsh.open(f.path)
    fe_mesh_from_gmsh()
  end
end

function fe_mesh_from_gmsh()
  entities = gmsh.model.getEntities()
  nodeTags, coord, parametricCoord = gmsh.model.mesh.getNodes()

  # find domain_dim
  ddim = -1
  for e in entities
    ddim = max(ddim,e[1])
  end
  if ddim == -1
    error("No entities in the msh file.")
  end
  D = ddim

  # find ambient_dim
  dtouched = [false,false,false]
  for node in nodeTags
    if !(coord[(node-1)*3+1] + 1 ≈ 1)
      dtouched[1] = true
    end
    if !(coord[(node-1)*3+2] + 1 ≈ 1)
      dtouched[2] = true
    end
    if !(coord[(node-1)*3+3] + 1 ≈ 1)
      dtouched[3] = true
    end
  end
  if dtouched[3]
    adim = 3
  elseif dtouched[2]
    adim = 2
  elseif dtouched[1]
    adim = 1
  else
    adim = 0
  end

  # Allocate an empty mesh
  mesh = SimpleFEMesh{SVector{adim,Float64}}(VOID,D)

  # Setup node coords
  nmin = minimum(nodeTags)
  nmax = maximum(nodeTags)
  nnodes = length(nodeTags)
  if !(nmax == nnodes && nmin == 1)
    error("Only consecutive node tags allowed.")
  end
  node_to_coords = zeros(SVector{adim,Float64},nnodes)
  m = zero(MVector{adim,Float64})
  for node in nodeTags
    for j in 1:adim
      k = (node-1)*3 + j
      xj = coord[k]
      m[j] = xj
    end
    node_to_coords[node] = m
  end
  node_coordinates!(mesh,node_to_coords)

  # Setup face nodes
  offsets = zeros(Int32,D+2)
  for d in 0:D
    elemTypes, elemTags, nodeTags = gmsh.model.mesh.getElements(d)
    ndfaces = 0
    for t in 1:length(elemTypes)
      ndfaces += length(elemTags[t])
    end
    offsets[d+2] = ndfaces
    ptrs = zeros(Int32,ndfaces+1)
    dface = 0
    for t in 1:length(elemTypes)
      elementName, dim, order, numNodes, nodeCoord =
        gmsh.model.mesh.getElementProperties(elemTypes[t])
      for e in 1:length(elemTags[t])
        dface += 1
        ptrs[dface+1] = numNodes
      end
    end
    prefix!(ptrs)
    ndata = ptrs[end]-1
    data = zeros(Int32,ndata)
    dface = 1
    for t in 1:length(elemTypes)
      p = ptrs[dface]-Int32(1)
      for (i,node) in enumerate(nodeTags[t])
        data[p+i] = node
      end
      dface += length(elemTags[t])
    end
    face_nodes!(mesh,d,JaggedArray(data,ptrs))
  end
  prefix!(offsets)

  # Setup face_ref_id
  for d in 0:D
    elemTypes, elemTags, nodeTags = gmsh.model.mesh.getElements(d)
    ndfaces = length(face_nodes(mesh,d))
    dface_to_refid = zeros(Int8,ndfaces)
    refid = 0
    dface = 0
    for t in 1:length(elemTypes)
      refid += 1
      for e in 1:length(elemTags[t])
        dface += 1
        dface_to_refid[dface] = refid
      end
    end
    face_ref_id!(mesh,d,dface_to_refid)
  end

  # Setup reference faces
  for d in 0:D
    elemTypes, elemTags, nodeTags = gmsh.model.mesh.getElements(d)
    refdfaces = []
    for t in 1:length(elemTypes)
      refface = ref_face_from_gmsh_eltype(elemTypes[t])
      push!(refdfaces,refface)
    end
    ref_faces!(mesh,d,refdfaces)
  end

  # Setup periodic nodes
  node_to_main_node = fill(Int32(INVALID),nnodes)
  for (dim,tag) in entities
    tagMaster, nodeTags, nodeTagsMaster, = gmsh.model.mesh.getPeriodicNodes(dim,tag)
    for i in 1:length(nodeTags)
      node = nodeTags[i]
      main_node = nodeTagsMaster[i]
      node_to_main_node[node] = main_node
    end
  end
  periodic_dep = collect(Int32,findall(i->i!=INVALID,node_to_main_node))
  periodic_indep = node_to_main_node[periodic_dep]
  periodic_nodes!(mesh,(periodic_dep,periodic_indep,ones(length(periodic_indep))))

  # Setup physical groups
  groups = physical_groups(mesh)
  for d in 0:D
    offset = Int32(offsets[d+1]-1)
    dimTags = gmsh.model.getPhysicalGroups(d)
    for (dim,tag) in dimTags
      @boundscheck @assert dim == d
      g_entities = gmsh.model.getEntitiesForPhysicalGroup(dim,tag)
      ndfaces_in_group = 0
      for entity in g_entities
        elemTypes, elemTags, nodeTags = gmsh.model.mesh.getElements(dim,entity)
        for t in 1:length(elemTypes)
          ndfaces_in_group += length(elemTags[t])
        end
      end
      dfaces_in_group = zeros(Int32,ndfaces_in_group)
      ndfaces_in_group = 0
      for entity in g_entities
        elemTypes, elemTags, nodeTags = gmsh.model.mesh.getElements(dim,entity)
        for t in 1:length(elemTypes)
          for etag in elemTags[t]
            ndfaces_in_group += 1
            dfaces_in_group[ndfaces_in_group] = Int32(etag)-offset
          end
        end
      end
      groupname = gmsh.model.getPhysicalName(dim,tag)
      add_group!(groups,d,groupname,tag)
      group_faces!(groups,dfaces_in_group,d,Int(tag))
    end
  end

  mesh
end

function ref_face_from_gmsh_eltype(eltype)
  if eltype == 1
    Meshes.Segment(Meshes.Point(0),Meshes.Point(1))
  elseif eltype == 2
    Meshes.Triangle(Meshes.Point.([(0,0),(1,0),(0,1)]))
  elseif eltype == 3
    Meshes.Quadrangle(Meshes.Point.([(0,0),(1,0),(1,1),(0,1)]))
  elseif eltype == 4
    Meshes.Tetrahedron(Meshes.Point.([(0,0,0),(1,0,0),(0,1,0),(0,0,1)]))
  elseif eltype == 5
    Meshes.Hexahedron(Meshes.Point.([(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)]))
  elseif eltype == 15
    Meshes.Point(SVector{0,Float64}())
  else
    error("Unsupported element type. elemType: $etype")
  end
end

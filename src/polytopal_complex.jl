

function linear_polytope end

default_linear_polytope(a) = a
linear_polytope(a) = default_linear_polytope(a)

function default_polytope_boundary(p)
  x = node_coordinates(p)
  T = eltype(x)
  D = domain_dim(p)-1
  mesh = GenericFEMesh{T}(VOID,D)
  node_coordinates!(mesh,x)
  for d in 0:D
    face_nodes!(mesh,face_nodes(p,d),d)
    face_ref_id!(mesh,face_ref_id(p,d),d)
    ref_faces!(mesh,ref_faces(p,d),d)
  end
  polytopal_complex(mesh)
end

polytope_boundary(p) = default_polytope_boundary(p)

function polytope_face_incedence(a,d1,d2)
  d = domain_dim(a)
  if d==d1
    nd2faces = num_faces(a,d2)
    GenericJaggedArray([[ Int32(d2face) for d2face in 1:nd2faces]])
  elseif d==d2
    nd1faces = num_faces(a,d1)
    GenericJaggedArray([[ Int32(1)] for d1face in 1:nd1faces])
  else
    poly = polytope_boundary(a)
    face_incidence(poly,d1,d2)
  end
end

function vertex_node end
function node_vertex end
function mesh_faces end
function mesh_face_vertices end
function polytopal_complex end

polytopal_complex(mesh) = PolyComplexFromFEMesh(mesh)

struct PolyComplexFromFEMesh
  mesh::Any
  buffer::Dict{Symbol,Any}
end

function PolyComplexFromFEMesh(mesh)
  buffer = Dict{Symbol,Any}()
  PolyComplexFromFEMesh(mesh,buffer)
end

domain_dim(m::PolyComplexFromFEMesh) = domain_dim(m.mesh)
ambient_dim(m::PolyComplexFromFEMesh) = ambient_dim(m.mesh)
node_coordinates(m::PolyComplexFromFEMesh) = node_coordinates(m.mesh)
periodic_nodes(m::PolyComplexFromFEMesh) = periodic_nodes(m.mesh)
hanging_nodes(m::PolyComplexFromFEMesh) = hanging_nodes(m.mesh)
is_hanging(m::PolyComplexFromFEMesh) = is_hanging(m.mesh)
is_periodic(m::PolyComplexFromFEMesh) = is_periodic(m.mesh)
is_simplex(m::PolyComplexFromFEMesh) = is_simplex(m.mesh)
is_hypercube(m::PolyComplexFromFEMesh) = is_hypercube(m.mesh)
face_vertices(m::PolyComplexFromFEMesh,rank) = face_incidence(m,rank,0)

function mesh_faces!(m::PolyComplexFromFEMesh,f,d)
  if ! haskey(m.buffer,:mesh_faces)
    D = domain_dim(m)
    m.buffer[:mesh_faces] = Vector{Vector{Int32}}(undef,D+1)
  end
  m.buffer[:mesh_faces][d+1] = f
end


function node_vertex(mesh::PolyComplexFromFEMesh)
  if !haskey(mesh.buffer,:node_vertex)
    _setup_vertices!(mesh)
  end
  mesh.buffer[:node_vertex]
end

function vertex_node(mesh::PolyComplexFromFEMesh)
  if !haskey(mesh.buffer,:vertex_node)
    _setup_vertices!(mesh)
  end
  mesh.buffer[:vertex_node]
end

function _setup_vertices!(poly)
  mesh = poly.mesh
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
  poly.buffer[:node_vertex] = node_to_vertex
  poly.buffer[:vertex_node] = vertex_to_node
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
  node_to_vertex = fill(Int32(INVALID),nnodes)
  node_to_vertex[vertex_to_node] .= Int32(1):Int32(nvertices)
  node_to_vertex, vertex_to_node
end

function mesh_face_vertices(mesh::PolyComplexFromFEMesh,d)
  if !haskey(mesh.buffer,:mesh_face_vertices)
    J = typeof(GenericJaggedArray(Vector{Int32}[]))
    mesh.buffer[:mesh_face_vertices] = Vector{J}(undef,domain_dim(mesh)+1)
  end
  if !isassigned(mesh.buffer[:mesh_face_vertices],d+1)
    _setup_mesh_face_vertices!(mesh,d)
  end
  mesh.buffer[:mesh_face_vertices][d+1]
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

function mesh_faces(mesh::PolyComplexFromFEMesh,d)
  if !haskey(mesh.buffer,:mesh_faces)
    mesh.buffer[:mesh_faces] = Vector{Vector{Int32}}(undef,domain_dim(mesh)+1)
  end
  if !isassigned(mesh.buffer[:mesh_faces],d+1)
    _setup_mesh_faces!(mesh,d)
  end
  mesh.buffer[:mesh_faces][d+1]
end

function  _setup_mesh_faces!(poly,d)
  mesh = poly.mesh
  if d != 0
    poly.buffer[:mesh_faces][d+1] = collect(Int32,1:num_faces(mesh,d))
  else
    poly.buffer[:mesh_faces][d+1] = map(first,mesh_face_vertices(poly,0))
  end
end

function group_collection(m::PolyComplexFromFEMesh)
  if !haskey(m.buffer,:group_collection)
    mesh_groups = group_collection(m.mesh)
    groups = GroupCollection(VOID,domain_dim(m))
    D = domain_dim(m)
    for d in 0:D
      for id in group_ids(mesh_groups,d)
        name = group_name(mesh_groups,d,id)
        add_group!(groups,d,name,id)
        mesh_face_to_face = mesh_faces(m,d)
        mesh_face_in_group = group_faces(mesh_groups,d,id)
        face_in_group = mesh_face_to_face[mesh_face_in_group]
        group_faces!(groups,face_in_group,d,id)
      end
    end
    m.buffer[:group_collection] = groups
  end
  m.buffer[:group_collection]
end
function group_collection!(m::PolyComplexFromFEMesh,groups)
  m.buffer[:group_collection] = groups
end

function num_faces(m::PolyComplexFromFEMesh,d)
  @boundscheck @assert d <= domain_dim(m)
  if ! haskey(m.buffer,:num_faces)
    m.buffer[:num_faces] = fill(INVALID,domain_dim(m)+1)
  end
  if m.buffer[:num_faces][d+1] == INVALID
    if d == domain_dim(m)
      m.buffer[:num_faces][d+1] = num_faces(m.mesh,d)
    elseif d == 0
      m.buffer[:num_faces][d+1] = length(vertex_node(m))
    else
      face_incidence(m,d+1,d)
      @boundscheck @assert m.buffer[:num_faces][d+1] != INVALID
    end
  end
  m.buffer[:num_faces][d+1]
end

function face_ref_id(a::PolyComplexFromFEMesh,m)
  d = domain_dim(a)
  if !haskey(a.buffer,:face_ref_id)
    a.buffer[:face_ref_id] = Vector{Vector{Int8}}(undef,d+1)
    a.buffer[:ref_faces] = Vector{Vector{Any}}(undef,d+1)
  end
  if !isassigned(a.buffer[:face_ref_id],m+1)
    if m == d
       a.buffer[:face_ref_id][m+1] = face_ref_id(a.mesh,m)
    else
      _ref_faces!(a,m)
    end
  end
  a.buffer[:face_ref_id][m+1]
end

function ref_faces(a::PolyComplexFromFEMesh,m)
  d = domain_dim(a)
  if !haskey(a.buffer,:ref_faces)
    a.buffer[:face_ref_id] = Vector{Vector{Int8}}(undef,d+1)
    a.buffer[:ref_faces] = Vector{Vector{Any}}(undef,d+1)
  end
  if !isassigned(a.buffer[:ref_faces],m+1)
    if m == d
       a.buffer[:ref_faces][m+1] = ref_faces(a.mesh,m)
    else
      _ref_faces!(a,m)
    end
  end
  a.buffer[:ref_faces][m+1]
end

function _ref_faces!(a,m)
  @assert m<domain_dim(a)
  d = m+1
  nmfaces = num_faces(a,m)
  drefid_refdface = ref_faces(a,d)
  drefid_lmface_refmface = map(i->ref_faces(i,m)[face_ref_id(i,m)],drefid_refdface)
  dface_to_lmface_to_mface = face_incidence(a,d,m)
  dface_to_drefid = face_ref_id(a,d)
  mface_to_u, u_to_refmface = _ref_faces(
    nmfaces,
    drefid_refdface,
    drefid_lmface_refmface,
    dface_to_lmface_to_mface,
    dface_to_drefid)
  a.buffer[:ref_faces][m+1] = u_to_refmface
  a.buffer[:face_ref_id][m+1] = mface_to_u
end

function _ref_faces(
  nmfaces,
  drefid_refdface,
  drefid_lmface_refmface,
  dface_to_lmface_to_mface,
  dface_to_drefid)

  i_to_drefid = Int[]
  i_to_lmface = Int[]
  i_to_refmface = Any[]
  for (drefid,lmface_refmface) in enumerate(drefid_lmface_refmface)
    for (lmface,refmface) in enumerate(lmface_refmface)
      push!(i_to_drefid,drefid)
      push!(i_to_lmface,lmface)
      push!(i_to_refmface,refmface)
    end
  end
  u_to_refmface = unique(Base.isequal,i_to_refmface)
  i_to_u = indexin(i_to_refmface,u_to_refmface)
  drefid_lmface_u = Vector{Int8}[]
  for (drefid,lmface_refmface) in enumerate(drefid_lmface_refmface)
    push!(drefid_lmface_u,zeros(Int8,length(lmface_refmface)))
  end
  for i in 1:length(i_to_drefid)
    drefid = i_to_drefid[i]
    lmface = i_to_lmface[i]
    drefid_lmface_u[drefid][lmface] = i_to_u[i]
  end
  mface_to_u = zeros(Int8,nmfaces)
  for (dface,drefid) in enumerate(dface_to_drefid)
   lmface_u = drefid_lmface_u[drefid]
   lmface_to_mface = dface_to_lmface_to_mface[dface]
   for (lmface,u) in enumerate(lmface_u)
     mface = lmface_to_mface[lmface]
     mface_to_u[mface] = u
   end
  end
  mface_to_u, u_to_refmface
end

function face_incidence(a::PolyComplexFromFEMesh,m,n)
  d = domain_dim(a)
  if !haskey(a.buffer,:face_incidence)
    J = typeof(GenericJaggedArray(Vector{Int32}[]))
    a.buffer[:face_incidence] = Matrix{J}(undef,d+1,d+1)
  end

  if !isassigned(a.buffer[:face_incidence],m+1,n+1)
    if m==d && n==0
      a.buffer[:face_incidence][m+1,n+1] = mesh_face_vertices(a,d)
    elseif m==d && n==d
      a.buffer[:face_incidence][m+1,n+1] = _face_interior(num_faces(a.mesh,d))
    elseif n==d && m==0
       cell_to_vertices = mesh_face_vertices(a,d)
       nvertices = length(vertex_node(a))
       a.buffer[:face_incidence][m+1,n+1] = _face_coboundary(cell_to_vertices,nvertices)
    elseif n==0 && m==0
       nvertices = length(vertex_node(a))
       a.buffer[:face_incidence][m+1,n+1] = _face_interior(nvertices)
    elseif m==d && n==(d-1)
      _face_boundary!(a,mesh_face_vertices(a,n),m,n)
    elseif n==0
      _face_vertices!(a,m)
    elseif m==n
      a.buffer[:face_incidence][m+1,n+1] = _face_interior(num_faces(a,n))
    elseif n+1==m
      _face_boundary!(a,mesh_face_vertices(a,n),m,n)
    elseif m>n
      _face_boundary!(a,face_vertices(a,n),m,n)
    else
       nface_to_mfaces = face_incidence(a,n,m)
       nmfaces = num_faces(a,m)
       a.buffer[:face_incidence][m+1,n+1] = _face_coboundary(nface_to_mfaces,nmfaces)
    end
  end
  a.buffer[:face_incidence][m+1,n+1]
end

function face_nodes(a::PolyComplexFromFEMesh,m)
  d = domain_dim(a)
  if !haskey(a.buffer,:face_nodes)
    J = typeof(GenericJaggedArray(Vector{Int32}[]))
    a.buffer[:face_nodes] = Vector{J}(undef,d+1)
  end
  if !isassigned(a.buffer[:face_nodes],m+1)
    if m == d
      a.buffer[:face_nodes][m+1] = face_nodes(a.mesh,d)
    else
      _face_nodes!(a,m)
    end
  end
  a.buffer[:face_nodes][m+1]
end

function _face_interior(nmfaces)
  ptrs = zeros(Int32,nmfaces+1)
  for mface in 1:nmfaces
    ptrs[mface+1] = 1
  end
  prefix!(ptrs)
  data = collect(Int32,1:nmfaces)
  mface_to_mfaces = GenericJaggedArray(data,ptrs)
  mface_to_mfaces
end

function _face_coboundary(nface_to_mfaces,nmfaces)
  ptrs = zeros(Int32,nmfaces+1)
  nnfaces = length(nface_to_mfaces)
  for nface in 1:nnfaces
    mfaces = nface_to_mfaces[nface]
    for mface in mfaces
      ptrs[mface+1] += Int32(1)
    end
  end
  prefix!(ptrs)
  ndata = ptrs[end]-1
  data = zeros(Int32,ndata)
  for nface in 1:nnfaces
    mfaces = nface_to_mfaces[nface]
    for mface in mfaces
      p = ptrs[mface]
      data[p] = nface
      ptrs[mface] += Int32(1)
    end
  end
  rewind!(ptrs)
  mface_to_nfaces = GenericJaggedArray(data,ptrs)
  mface_to_nfaces
end

function _face_boundary!(femesh,dface_to_vertices,D,d)
  Dface_to_vertices = face_incidence(femesh,D,0)
  vertex_to_Dfaces = face_incidence(femesh,0,D)
  nvertices = length(vertex_to_Dfaces)
  vertex_to_dfaces = _face_coboundary(dface_to_vertices,nvertices)
  Dface_to_refid = face_ref_id(femesh,D)
  Drefid_to_ldface_to_lvertices = ref_face_incidence(femesh,D,d,0)
  Dface_to_dfaces, nnewdface = _face_boundary(
    Dface_to_vertices,
    vertex_to_Dfaces,
    dface_to_vertices,
    vertex_to_dfaces,
    Dface_to_refid,
    Drefid_to_ldface_to_lvertices)
  if !haskey(femesh.buffer,:num_faces)
    femesh.buffer[:num_faces] = fill(INVALID,domain_dim(femesh)+1)
  end
  femesh.buffer[:num_faces][d+1] = nnewdface
  femesh.buffer[:face_incidence][D+1,d+1] = Dface_to_dfaces
end

function _face_boundary(
  Dface_to_vertices,
  vertex_to_Dfaces,
  dface_to_vertices,
  vertex_to_dfaces,
  Dface_to_refid,
  Drefid_to_ldface_to_lvertices)

  # Count
  ndfaces = length(dface_to_vertices)
  nDfaces = length(Dface_to_vertices)
  nvertices = length(vertex_to_Dfaces)
  maxldfaces = 0
  for ldface_to_lvertices in Drefid_to_ldface_to_lvertices
    maxldfaces = max(maxldfaces,length(ldface_to_lvertices))
  end
  maxDfaces = 0
  for vertex in 1:length(vertex_to_Dfaces)
    Dfaces = vertex_to_Dfaces[vertex]
    maxDfaces = max(maxDfaces,length(Dfaces))
  end
  # Allocate output
  ptrs = zeros(Int32,nDfaces+1)
  for Dface in 1:nDfaces
    Drefid = Dface_to_refid[Dface]
    ldface_to_lvertices = Drefid_to_ldface_to_lvertices[Drefid]
    ptrs[Dface+1] = length(ldface_to_lvertices)
  end
  prefix!(ptrs)
  ndata = ptrs[end]-1
  data = fill(Int32(INVALID),ndata)
  Dface_to_dfaces = GenericJaggedArray(data,ptrs)
  # Main loop
  Dfaces1 = fill(Int32(INVALID),maxDfaces)
  Dfaces2 = fill(Int32(INVALID),maxDfaces)
  ldfaces1 = fill(Int32(INVALID),maxDfaces)
  nDfaces1 = 0
  nDfaces2 = 0
  newdface = Int32(ndfaces)
  for Dface in 1:nDfaces
    Drefid = Dface_to_refid[Dface]
    ldface_to_lvertices = Drefid_to_ldface_to_lvertices[Drefid]
    lvertex_to_vertex = Dface_to_vertices[Dface]
    ldface_to_dface = Dface_to_dfaces[Dface]
    for (ldface,lvertices) in enumerate(ldface_to_lvertices)
      # Do nothing if this local face has already been processed by
      # a neighbor
      if ldface_to_dface[ldface] != Int32(INVALID)
        continue
      end
      # Find if there is already a global d-face for this local d-face
      # if yes, then use the global id of this d-face
      # if not, create a new one
      dface2 = Int32(INVALID)
      fill!(Dfaces1,Int32(INVALID))
      fill!(Dfaces2,Int32(INVALID))
      vertices = view(lvertex_to_vertex,lvertices)
      for (i,lvertex) in enumerate(lvertices)
        vertex = lvertex_to_vertex[lvertex]
        dfaces = vertex_to_dfaces[vertex]
        for dface1 in dfaces
          vertices1 = dface_to_vertices[dface1]
          if _contain_same_ids(vertices,vertices1)
            dface2 = dface1
            break
          end
        end
        if dface2 != Int32(INVALID)
          break
        end
      end
      if dface2 == Int32(INVALID)
        newdface += Int32(1)
        dface2 = newdface
      end
      # Find all D-faces around this local d-face
      for (i,lvertex) in enumerate(lvertices)
        vertex = lvertex_to_vertex[lvertex]
        Dfaces = vertex_to_Dfaces[vertex]
        if i == 1
          copyto!(Dfaces1,Dfaces)
          nDfaces1 = length(Dfaces)
        else
          copyto!(Dfaces2,Dfaces)
          nDfaces2 = length(Dfaces)
          _intersect_ids!(Dfaces1,nDfaces1,Dfaces2,nDfaces2)
        end
      end
      # Find their correspondent local d-face and set the d-face
      for Dface1 in Dfaces1
        if Dface1 != INVALID
          Drefid1 = Dface_to_refid[Dface1]
          lvertex1_to_vertex1 = Dface_to_vertices[Dface1]
          ldface1_to_lvertices1 = Drefid_to_ldface_to_lvertices[Drefid1]
          ldface2 = Int32(INVALID)
          for (ldface1,lvertices1) in enumerate(ldface1_to_lvertices1)
            vertices1 = view(lvertex1_to_vertex1,lvertices1)
            if _contain_same_ids(vertices,vertices1)
              ldface2 = ldface1
              break
            end
          end
          Dface_to_dfaces[Dface1][ldface2] = dface2
        end
      end
    end # (ldface,lvertices)
  end # Dface
  Dface_to_dfaces, newdface
end

function _intersect_ids!(a,na,b,nb)
  function findeq!(i,a,b,nb)
    for j in 1:nb
      if a[i] == b[j]
        return
      end
    end
    a[i] = INVALID
    return
  end
  for i in 1:na
    if a[i] == INVALID
      continue
    end
    findeq!(i,a,b,nb)
  end
end

function _contain_same_ids(a,b)
  function is_subset(a,b)
    for i in 1:length(a)
      v = a[i]
      if v == INVALID
        continue
      end
      c = find_eq(v,b)
      if c == false; return false; end
    end
    return true
  end
  function find_eq(v,b)
    for vs in b
      if v == vs
        return true
      end
    end
    return false
  end
  c = is_subset(a,b)
  if c == false; return false; end
  c = is_subset(b,a)
  if c == false; return false; end
  return true
end

function _face_vertices!(femesh,d)
  D = d+1
  dface_to_vertices = mesh_face_vertices(femesh,d)
  Dface_to_dfaces = face_incidence(femesh,D,d)
  Dface_to_vertices = face_incidence(femesh,D,0)
  Dface_to_refid = face_ref_id(femesh,D)
  Drefid_to_ldface_to_lvertices = ref_face_incidence(femesh,D,d,0)
  nnewdfaces = num_faces(femesh,d)
  femesh.buffer[:face_incidence][d+1,0+1] = _face_vertices(
    nnewdfaces,
    dface_to_vertices,
    Dface_to_dfaces,
    Dface_to_vertices,
    Dface_to_refid,
    Drefid_to_ldface_to_lvertices)
end

function _face_nodes!(femesh,d)
  D = d+1
  dface_to_vertices = face_nodes(femesh.mesh,d)
  Dface_to_dfaces = face_incidence(femesh,D,d)
  Dface_to_vertices = face_nodes(femesh,D)
  Dface_to_refid = face_ref_id(femesh,D)
  Drefid_to_ldface_to_lnodes = ref_face_nodes(femesh,D,d)
  nnewdfaces = num_faces(femesh,d)
  femesh.buffer[:face_nodes][d+1] = _face_vertices(
    nnewdfaces,
    dface_to_vertices,
    Dface_to_dfaces,
    Dface_to_vertices,
    Dface_to_refid,
    Drefid_to_ldface_to_lnodes)
end

function _face_vertices(
    nnewdfaces,
    dface_to_vertices,
    Dface_to_dfaces,
    Dface_to_vertices,
    Dface_to_refid,
    Drefid_to_ldface_to_lvertices)

  dface_to_touched = fill(false,nnewdfaces)
  ptrs = zeros(Int32,nnewdfaces+1)
  for (dface,vertices) in enumerate(dface_to_vertices)
    ptrs[dface+1] = length(vertices)
    dface_to_touched[dface] = true
  end
  for (Dface,ldface_to_dface) in enumerate(Dface_to_dfaces)
    Drefid = Dface_to_refid[Dface]
    ldface_to_lvertices = Drefid_to_ldface_to_lvertices[Drefid]
    for (ldface,dface) in enumerate(ldface_to_dface)
      if !dface_to_touched[dface]
        lvertices = ldface_to_lvertices[ldface]
        ptrs[dface+1] = length(lvertices)
        dface_to_touched[dface] = true
      end
    end
  end
  prefix!(ptrs)
  ndata = ptrs[end]-1
  data = zeros(Int32,ndata)
  fill!(dface_to_touched,false)
  for (dface,vertices) in enumerate(dface_to_vertices)
    p = ptrs[dface] - Int32(1)
    for (lvertex,vertex) in enumerate(vertices)
      data[p+lvertex] = vertex
    end
    dface_to_touched[dface] = true
  end
  for (Dface,ldface_to_dface) in enumerate(Dface_to_dfaces)
    Drefid = Dface_to_refid[Dface]
    lvertex_to_vertex = Dface_to_vertices[Dface]
    ldface_to_lvertices = Drefid_to_ldface_to_lvertices[Drefid]
    for (ldface,dface) in enumerate(ldface_to_dface)
      if !dface_to_touched[dface]
        lvertices = ldface_to_lvertices[ldface]
        p = ptrs[dface] - Int32(1)
        for (i,lvertex) in enumerate(lvertices)
          data[p+i] = lvertex_to_vertex[lvertex]
        end
        dface_to_touched[dface] = true
      end
    end
  end
  GenericJaggedArray(data,ptrs)
end


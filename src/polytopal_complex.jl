
function polytopal_complex end

polytopal_complex(mesh) = GenericPolyComplex(mesh)
function polytopal_complex(m::SimpleFEMesh)
  if !haskey(m.buffer,:polytopal_complex)
    polycomplex = GenericPolyComplex(m)
    m.buffer[:polytopal_complex] = polycomplex
  end
  m.buffer[:polytopal_complex]
end

struct GenericPolyComplex
  mesh::Any
  buffer::Dict{Symbol,Any}
end

function GenericPolyComplex(mesh)
  buffer = Dict{Symbol,Any}()
  GenericPolyComplex(mesh,buffer)
end

domain_dim(m::GenericPolyComplex) = domain_dim(m.mesh)
ambient_dim(m::GenericPolyComplex) = ambient_dim(m.mesh)
node_coordinates(m::GenericPolyComplex) = node_coordinates(m.mesh)
periodic_nodes(m::GenericPolyComplex) = periodic_nodes(m.mesh)
hanging_nodes(m::GenericPolyComplex) = hanging_nodes(m.mesh)
is_hanging(m::GenericPolyComplex) = is_hanging(m.mesh)
is_periodic(m::GenericPolyComplex) = is_periodic(m.mesh)
is_simplex(m::GenericPolyComplex) = is_simplex(m.mesh)
is_hypercube(m::GenericPolyComplex) = is_hypercube(m.mesh)
node_vertex(m::GenericPolyComplex) = node_vertex(m.mesh)
vertex_node(m::GenericPolyComplex) = vertex_node(m.mesh)
face_vertices(m::GenericPolyComplex,rank) = face_incidence(m,rank,0)

function physical_groups(m::GenericPolyComplex)
  if !haskey(m.buffer,:physical_groups)
    m.buffer[:physical_groups] = deepcopy(physical_groups(m.mesh))
  end
  m.buffer[:physical_groups]
end
function physical_groups!(m::GenericPolyComplex,groups)
  m.buffer[:physical_groups] = groups
end

function num_faces(m::GenericPolyComplex,d)
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

function face_ref_id(a::GenericPolyComplex,m)
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

function ref_faces(a::GenericPolyComplex,m)
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
  drefid_lmface_refmface = map(i->ref_faces(i,m),drefid_refdface)
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
      push!(i_to_drefid,drefif)
      push!(i_to_lmface,drefif)
      push!(i_to_refmface,refmface)
    end
  end
  u_to_refmface = unique(i_to_refmface)
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

function face_incidence(a::GenericPolyComplex,m,n)
  d = domain_dim(a)
  if !haskey(a.buffer,:face_incidence)
    J = typeof(JaggedArray(Vector{Int32}[]))
    a.buffer[:face_incidence] = Matrix{J}(undef,d+1,d+1)
  end

  if !isassigned(a.buffer[:face_incidence],m+1,n+1)
    if m==d && n==0
      a.buffer[:face_incidence][m+1,n+1] = face_vertices(a.mesh,d)
    elseif m==d && n==d
      a.buffer[:face_incidence][m+1,n+1] = _face_interior(num_faces(a.mesh,d))
    elseif n==d && m==0
       cell_to_vertices = face_vertices(a.mesh,d)
       nvertices = length(vertex_node(a.mesh))
       a.buffer[:face_incidence][m+1,n+1] = _face_coboundary(cell_to_vertices,nvertices)
    elseif n==0 && m==0
       nvertices = length(vertex_node(a.mesh))
       a.buffer[:face_incidence][m+1,n+1] = _face_interior(nvertices)
    elseif m==d && n==(d-1)
      _face_boundary!(a,a.mesh,m,n)
    elseif n==0
      _face_vertices!(a,m)
    elseif m==n
      a.buffer[:face_incidence][m+1,n+1] = _face_interior(num_faces(a,n))
    elseif m>n
      _face_boundary!(a,a,m,n)
    else
       nface_to_mfaces = face_incidence(a,n,m)
       nmfaces = num_faces(a,m)
       a.buffer[:face_incidence][m+1,n+1] = _face_coboundary(nface_to_mfaces,nmfaces)
    end
  end
  a.buffer[:face_incidence][m+1,n+1]
end

function face_nodes(a::GenericPolyComplex,m)
  d = domain_dim(a)
  if !haskey(a.buffer,:face_nodes)
    J = typeof(JaggedArray(Vector{Int32}[]))
    a.buffer[:face_nodes] = Vector{J}(undef,d+1)
  end
  if !isassigned(a.buffer[:face_nodes],m+1)
    _face_nodes!(a,m)
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
  mface_to_mfaces = JaggedArray(data,ptrs)
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
  mface_to_nfaces = JaggedArray(data,ptrs)
  mface_to_nfaces
end

function _face_boundary!(femesh,mesh,D,d)
  Dface_to_vertices = face_incidence(femesh,D,0)
  vertex_to_Dfaces = face_incidence(femesh,0,D)
  dface_to_vertices = face_vertices(mesh,d)
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
  Dface_to_dfaces = JaggedArray(data,ptrs)
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
  dface_to_vertices = face_vertices(femesh.mesh,d)
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

function _face_nodes!(a,d)
  D = d+1
  dface_to_vertices = face_nodes(a.mesh,d)
  Dface_to_dfaces = face_incidence(a,D,d)
  Dface_to_vertices = face_nodes(a,D)
  Dface_to_refid = face_ref_id(femesh,D)
  Drefid_to_ldface_to_lnodes = ref_face_nodes(femesh,D,d)
  nnewdfaces = num_faces(femesh,d)
  a.buffer[:face_nodes][d+1] = _face_vertices(
    nnewdfaces,
    dface_to_vertices,
    Dface_to_dfaces,
    Dface_to_vertices,
    Dface_to_refid,
    Drefid_to_ldface_to_lvertices)
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
  JaggedArray(data,ptrs)
end


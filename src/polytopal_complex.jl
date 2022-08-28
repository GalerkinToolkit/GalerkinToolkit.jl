

"""
   linear_polytope(poly)
Return the linear version of polytope `poly`.
"""
function linear_polytope end

"""
    polytope_boundary(poly)
Return a polytopal complex representing the boundary
of polytope `poly`.
"""
function polytope_boundary end
polytope_boundary(p) = default_polytope_boundary(p)
function default_polytope_boundary(p)
  x = node_coordinates(p)
  T = eltype(x)
  D = domain_dim(p)-1
  mesh = FEMesh{T}(D)
  node_coordinates!(mesh,x)
  for d in 0:D
    face_nodes!(mesh,face_nodes(p,d),d)
    face_ref_id!(mesh,face_ref_id(p,d),d)
    ref_faces!(mesh,ref_faces(p,d),d)
  end
  polytopal_complex(mesh)
end
function default_polytope_face_faces(poly,dim_from,dim_to)
  d = domain_dim(poly)
  if d==dim_from
    nd2faces = num_faces(poly,dim_to)
    JaggedArray([[ Int32(d2face) for d2face in 1:nd2faces]])
  elseif d==dim_to
    nd1faces = num_faces(poly,d1dim_from)
    JaggedArray([[ Int32(1)] for d1face in 1:nd1faces])
  else
    poly = polytope_boundary(poly)
    face_faces(poly,dim_from,dim_to)
  end
end

"""
    mesh_faces(complex,dim)
Return the face mapping between the faces in the mesh `fe_mesh(complex)`
and the faces in the polytopal complex `complex`.
"""
function mesh_faces end

"""
    mesh_faces!(complex,meshfaces,dim)
Sets the face mapping between the faces in the mesh `fe_mesh(complex)`
and the faces in the polytopal complex `complex`.
"""
function mesh_faces! end

"""
    polytopal_complex(geo)
Return the polytopal complex associated with `geo`.
"""
function polytopal_complex end

struct PolytopalComplex{M}
  fe_mesh::M
  buffer::Dict{Symbol,Any}
end

fe_mesh(a::PolytopalComplex) = a.fe_mesh

function PolytopalComplex(femesh)
  buffer = Dict{Symbol,Any}()
  PolytopalComplex(femesh,buffer)
end

domain_dim(m::PolytopalComplex) = domain_dim(m.fe_mesh)
ambient_dim(m::PolytopalComplex) = ambient_dim(m.fe_mesh)
node_coordinates(m::PolytopalComplex) = node_coordinates(m.fe_mesh)
periodic_nodes(m::PolytopalComplex) = periodic_nodes(m.fe_mesh)
hanging_nodes(m::PolytopalComplex) = hanging_nodes(m.fe_mesh)
is_simplex(m::PolytopalComplex) = is_simplex(m.fe_mesh)
is_hypercube(m::PolytopalComplex) = is_hypercube(m.fe_mesh)
node_vertex(m::PolytopalComplex) = node_vertex(m.fe_mesh)
vertex_node(m::PolytopalComplex) = vertex_node(m.fe_mesh)

function mesh_faces!(m::PolytopalComplex,f,d)
  if ! haskey(m.buffer,:mesh_faces)
    D = domain_dim(m)
    m.buffer[:mesh_faces] = Vector{Vector{Int32}}(undef,D+1)
  end
  m.buffer[:mesh_faces][d+1] = f
end

function mesh_faces(mesh::PolytopalComplex,d)
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
    poly.buffer[:mesh_faces][d+1] = map(first,face_vertices(mesh,0))
  end
end

function physical_groups(m::PolytopalComplex)
  if !haskey(m.buffer,:physical_groups)
    mesh_physical_groups = physical_groups(m.mesh)
    groups = PhysicalGroupCollection(VOID,domain_dim(m))
    D = domain_dim(m)
    for d in 0:D
      for id in physical_group_ids(mesh_physical_groups,d)
        name = physical_group_name(mesh_physical_groups,d,id)
        new_physical_group!(groups,d,name,id)
        mesh_face_to_face = mesh_faces(m,d)
        mesh_face_in_physical_group = physical_group_faces(mesh_physical_groups,d,id)
        face_in_physical_group = mesh_face_to_face[mesh_face_in_physical_group]
        physical_group_faces!(groups,face_in_physical_group,d,id)
      end
    end
    m.buffer[:physical_groups] = groups
  end
  m.buffer[:physical_groups]
end
function physical_groups!(m::PolytopalComplex,groups)
  m.buffer[:physical_groups] = groups
end

function num_faces(m::PolytopalComplex,d)
  @boundscheck @assert d <= domain_dim(m)
  if ! haskey(m.buffer,:num_faces)
    m.buffer[:num_faces] = fill(INVALID_ID,domain_dim(m)+1)
  end
  if m.buffer[:num_faces][d+1] == INVALID_ID
    if d == domain_dim(m)
      m.buffer[:num_faces][d+1] = num_faces(m.mesh,d)
    elseif d == 0
      m.buffer[:num_faces][d+1] = length(vertex_node(m))
    else
      face_faces(m,d+1,d)
      @boundscheck @assert m.buffer[:num_faces][d+1] != INVALID_ID
    end
  end
  m.buffer[:num_faces][d+1]
end

function face_ref_id(a::PolytopalComplex,m)
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

function ref_faces(a::PolytopalComplex,m)
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
  dface_to_lmface_to_mface = face_faces(a,d,m)
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

function face_faces(a::PolytopalComplex,m,n)
  d = domain_dim(a)
  if !haskey(a.buffer,:face_faces)
    J = typeof(GenericJaggedArray(Vector{Int32}[]))
    a.buffer[:face_faces] = Matrix{J}(undef,d+1,d+1)
  end

  if !isassigned(a.buffer[:face_faces],m+1,n+1)
    if m==d && n==0
      a.buffer[:face_faces][m+1,n+1] = face_vertices(a.fe_mesh,d)
    elseif m==d && n==d
      a.buffer[:face_faces][m+1,n+1] = _face_interior(num_faces(a.mesh,d))
    elseif n==d && m==0
       cell_to_vertices = face_vertices(a.fe_mesh,d)
       nvertices = length(vertex_node(a))
       a.buffer[:face_faces][m+1,n+1] = _face_coboundary(cell_to_vertices,nvertices)
    elseif n==0 && m==0
       nvertices = length(vertex_node(a))
       a.buffer[:face_faces][m+1,n+1] = _face_interior(nvertices)
    elseif m==d && n==(d-1)
      _face_boundary!(a,face_vertices(a.fe_mesh,n),m,n)
    elseif n==0
      _face_vertices!(a,m)
    elseif m==n
      a.buffer[:face_faces][m+1,n+1] = _face_interior(num_faces(a,n))
    elseif n+1==m
      _face_boundary!(a,face_vertices(a.fe_mesh,n),m,n)
    elseif m>n
      _face_boundary!(a,face_vertices(a,n),m,n)
    else
       nface_to_mfaces = face_faces(a,n,m)
       nmfaces = num_faces(a,m)
       a.buffer[:face_faces][m+1,n+1] = _face_coboundary(nface_to_mfaces,nmfaces)
    end
  end
  a.buffer[:face_faces][m+1,n+1]
end

function face_nodes(a::PolytopalComplex,m)
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
  Dface_to_vertices = face_faces(femesh,D,0)
  vertex_to_Dfaces = face_faces(femesh,0,D)
  nvertices = length(vertex_to_Dfaces)
  vertex_to_dfaces = _face_coboundary(dface_to_vertices,nvertices)
  Dface_to_refid = face_ref_id(femesh,D)
  Drefid_to_ldface_to_lvertices = ref_face_faces(femesh,D,d,0)
  Dface_to_dfaces, nnewdface = _face_boundary(
    Dface_to_vertices,
    vertex_to_Dfaces,
    dface_to_vertices,
    vertex_to_dfaces,
    Dface_to_refid,
    Drefid_to_ldface_to_lvertices)
  if !haskey(femesh.buffer,:num_faces)
    femesh.buffer[:num_faces] = fill(INVALID_ID,domain_dim(femesh)+1)
  end
  femesh.buffer[:num_faces][d+1] = nnewdface
  femesh.buffer[:face_faces][D+1,d+1] = Dface_to_dfaces
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
  data = fill(Int32(INVALID_ID),ndata)
  Dface_to_dfaces = GenericJaggedArray(data,ptrs)
  # Main loop
  Dfaces1 = fill(Int32(INVALID_ID),maxDfaces)
  Dfaces2 = fill(Int32(INVALID_ID),maxDfaces)
  ldfaces1 = fill(Int32(INVALID_ID),maxDfaces)
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
      if ldface_to_dface[ldface] != Int32(INVALID_ID)
        continue
      end
      # Find if there is already a global d-face for this local d-face
      # if yes, then use the global id of this d-face
      # if not, create a new one
      dface2 = Int32(INVALID_ID)
      fill!(Dfaces1,Int32(INVALID_ID))
      fill!(Dfaces2,Int32(INVALID_ID))
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
        if dface2 != Int32(INVALID_ID)
          break
        end
      end
      if dface2 == Int32(INVALID_ID)
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
        if Dface1 != INVALID_ID
          Drefid1 = Dface_to_refid[Dface1]
          lvertex1_to_vertex1 = Dface_to_vertices[Dface1]
          ldface1_to_lvertices1 = Drefid_to_ldface_to_lvertices[Drefid1]
          ldface2 = Int32(INVALID_ID)
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
    a[i] = INVALID_ID
    return
  end
  for i in 1:na
    if a[i] == INVALID_ID
      continue
    end
    findeq!(i,a,b,nb)
  end
end

function _contain_same_ids(a,b)
  function is_subset(a,b)
    for i in 1:length(a)
      v = a[i]
      if v == INVALID_ID
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
  dface_to_vertices =face_vertices(femesh.fe_mesh,d)
  Dface_to_dfaces = face_faces(femesh,D,d)
  Dface_to_vertices = face_faces(femesh,D,0)
  Dface_to_refid = face_ref_id(femesh,D)
  Drefid_to_ldface_to_lvertices = ref_face_faces(femesh,D,d,0)
  nnewdfaces = num_faces(femesh,d)
  femesh.buffer[:face_faces][d+1,0+1] = _face_vertices(
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
  Dface_to_dfaces = face_faces(femesh,D,d)
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


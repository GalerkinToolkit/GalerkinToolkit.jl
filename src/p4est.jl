

p4est_quadrants(p4est_ptr) = QuadrantIterator(p4est_ptr)

struct QuadrantIterator{T}
  p4est_ptr::Ptr{T}
end

function Base.iterate(iter::QuadrantIterator)
  p4est_ptr = iter.p4est_ptr
  p4est = unsafe_wrap(p4est_ptr)
  trees_ptr = p4est.trees
  trees = unsafe_wrap(trees_ptr)
  ntrees = Int(trees.elem_count)
  ntrees == 0 && return nothing
  itree = 1
  iquadrant = 0
  iterate(iter,(itree,iquadrant,trees))
end

function Base.iterate(iter::QuadrantIterator,state)
  itree, iquadrant, trees = state
  ntrees = Int(trees.elem_count)
  @assert itree <= ntrees
  tree = unsafe_load(Ptr{P4est.p4est_tree_t}(trees.array),itree)
  quadrants =  tree.quadrants
  nquadrants = Int(quadrants.elem_count)
  iquadrant += 1
  if nquadrants < iquadrant
    itree += 1
    ntrees < itree && return nothing
    iquadrant = 0
    return iterate(iter,(itree,iquadrant,trees))
  end
  quadrant = unsafe_load(Ptr{P4est.p4est_quadrant_t}(quadrants.array),iquadrant)
  quadrant, (itree,iquadrant,trees)
end

P4EST_VERTEX_PERMUTATION = [1,2,4,3]
P4EST_EDGE_PERMUTATION = [4,2,1,3]
p4est_vertex_permutation(q::Meshes.Quadrangle) = P4EST_VERTEX_PERMUTATION
function p4est_face_permutation(q::Meshes.Quadrangle,d)
  d==0 && return P4EST_VERTEX_PERMUTATION
  d==1 && return P4EST_EDGE_PERMUTATION
  d==2 && return [1]
  throw(DomainError(d))
end

struct P4estQuad{T}
  node_coordinates::Vector{SVector{2,T}}
end
function P4estQuad{T}() where T
  P4estQuad(SVector{2,T}[(0,0),(1,0),(0,1),(1,1)])
end
domain_dim(a::P4estQuad) = 2
is_hypercube(a::P4estQuad) = true
node_coordinates(a::P4estQuad) = a.node_coordinates
function face_ref_id(a::P4estQuad,d)
  d==0 && return fill(Int8(1),4)
  d==1 && return fill(Int8(1),4)
  d==2 && return fill(Int8(1),1)
  throw(DomainError(d))
end
function ref_faces(::Type{<:P4estQuad},d)
  d==0 && return [Meshes.Point(SVector{0,Float64}())]
  d==1 && return [Meshes.Segment(Meshes.Point(0),Meshes.Point(1))]
  d==2 && return [P4estQuad()]
  throw(DomainError(d))
end
ref_faces(a::P4estQuad,d) = ref_faces(typeof(a),d)
function face_nodes(a::P4estQuad,d)
  d==0 && return JaggedArray(Vector{Int32}[[1],[2],[3],[4]])
  d==1 && return JaggedArray(Vector{Int32}[[1,3],[2,4],[1,2],[3,4]])
  d==2 && return JaggedArray(Vector{Int32}[[1,2,3,4]])
  throw(DomainError(d))
end
const _P4EST_BUFFER = Dict{Symbol,Any}()
_P4EST_BUFFER[:P4estQuad] = Dict{Symbol,Any}()
function polytope_boundary(a::P4estQuad)
  if !haskey(_P4EST_BUFFER[:P4estQuad],:polytope_boundary)
    _P4EST_BUFFER[:P4estQuad][:polytope_boundary] = default_polytope_boundary(a)
  end
  _P4EST_BUFFER[:P4estQuad][:polytope_boundary]
end
face_incidence(a::P4estQuad,d1,d2) = polytope_face_incedence(a,d1,d2)
function vtk_mesh_cell(a::P4estQuad)
  nodes -> WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_QUAD,nodes[P4EST_VERTEX_PERMUTATION])
end

mutable struct P4estMeshRefiner{A,B}
  p4est_ptr::Ptr{A}
  coarse_mesh::B
  state::Symbol
end

domain_dim(a::P4estMeshRefiner) = domain_dim(a.coarse_mesh)
ambient_dim(a::P4estMeshRefiner) = ambient_dim(a.coarse_mesh)

function p4est_corner_neighbor(o,c)
  T = Cint
  level = o.level
  h = 2^(P4est.P4EST_MAXLEVEL-o.level)
  x = o.x + (2(T(c-1)&T(1))-1)*h
  y = o.y + (1(T(c-1)&T(2))-1)*h
  (;level,x,y)
end

function p4est_face_neighbor(o,f)
  T = Cint
  level = o.level
  h = 2^(P4est.P4EST_MAXLEVEL-o.level)
  x = o.x + ((f-1)==0 ? -1h : ((f-1)==1 ? 1h : 0h))
  y = o.y + ((f-1)==2 ? -1h : ((f-1)==3 ? 1h : 0h))
  (;level,x,y)
end

function p4est_mesh_refiner(geo;initial_level=0)
  if !MPI.Initialized()
    MPI.Init()
  end
  mesh = fe_mesh(geo)
  D = domain_dim(mesh)
  @assert num_faces(mesh,D) == 1 "For the moment connectivity creation only for a single tree"
  refface = first(ref_faces(mesh,D))
  perm = p4est_vertex_permutation(refface)
  ntrees = num_faces(mesh,D)
  nvertices = num_nodes(mesh)
  vertices = zeros(Cdouble,3,nvertices)
  for (node,x) in enumerate(node_coordinates(mesh))
    for d in 1:D
      vertices[d,node] = x[d]
    end
  end
  coff = [P4est.p4est_topidx_t(0)]
  nodes = face_nodes(mesh,D)[1]
  if D == 2
    tree_to_vertex = P4est.p4est_topidx_t[0,1,2,3][nodes][perm]
    tree_to_tree = P4est.p4est_topidx_t[0,0,0,0]
    tree_to_face = P4est.LibP4est.int8_t[0,1,2,3]
    conn_ptr = P4est.p4est_connectivity_new_copy(
      nvertices, ntrees, 0,
      vertices, tree_to_vertex,
      tree_to_tree, tree_to_face,
      C_NULL, coff, C_NULL, C_NULL);
    min_quadrants = 0
    fill_uniform = 1
    p4est_ptr = P4est.p4est_new_ext(
      MPI.COMM_WORLD,
      conn_ptr,
      min_quadrants,
      initial_level,
      fill_uniform,
      sizeof(Int),
      C_NULL,
      C_NULL)
  elseif D == 3
    error("Not yet implemented")
  else
    error("P4est only for 2d and 3d")
  end
  amr = P4estMeshRefiner(p4est_ptr,mesh,:usable)
  finalizer(destroy!,amr)
end

function destroy!(a::P4estMeshRefiner)
  if a.state === :usable
    p4est = unsafe_wrap(a.p4est_ptr)
    conn_ptr = p4est.connectivity
    if domain_dim(a) == 2
      P4est.p4est_connectivity_destroy(conn_ptr)
      P4est.p4est_destroy(a.p4est_ptr)
    elseif domain_dim(a) == 3
      P4est.p8est_connectivity_destroy(conn_ptr)
      P4est.p8est_destroy(a.p4est_ptr)
    else
      error("P4est only for 2d and 3d.")
    end
    a.state = :destroyed
  end
end

function _refine_fn(p4est_ptr,which_tree,quadrant_ptr)
  p4est = unsafe_wrap(p4est_ptr)
  quadrant = unsafe_wrap(quadrant_ptr)
  flags_ptr = Base.unsafe_convert(Ptr{Bool},p4est.user_pointer)
  elem_ptr = Base.unsafe_convert(Ptr{Int},quadrant.p.user_data)
  elem = unsafe_load(elem_ptr,1)
  elem == INVALID && return Int32(1)
  flag = unsafe_load(flags_ptr,elem)
  i = flag ? 1 : 0
  Int32(i)
end

function _coarsen_fn(p4est_ptr,which_tree,quadrants_ptr)
  p4est = unsafe_wrap(p4est_ptr)
  flags_ptr = Base.unsafe_convert(Ptr{Bool},p4est.user_pointer)
  quadrants = unsafe_wrap(Array,quadrants_ptr,4)
  r = true
  for i in 1:4
    quadrant = quadrants[i]
    elem_ptr = Base.unsafe_convert(Ptr{Int},quadrant.p.user_data)
    elem = unsafe_load(elem_ptr,1)
    if elem == INVALID
      flag = false
    else
      flag = unsafe_load(flags_ptr,elem)
    end
    r = r && (!flag)
  end
  r ? Int32(1) : Int32(0)
end

function _init_fn(p4est_ptr,which_tree,quadrant_ptr)
  quadrant = unsafe_wrap(quadrant_ptr)
  elem_ptr = Base.unsafe_convert(Ptr{Int},quadrant.p.user_data)
  unsafe_store!(elem_ptr,INVALID)
  Cvoid()
end

function refine!(amr::P4estMeshRefiner,_flags;num_levels=1)
  flags = convert(Vector{Bool},_flags)
  D = domain_dim(amr)
  if D == 3
    error("Not implemented yet")
  end
  p4est_ptr = amr.p4est_ptr
  p4est = unsafe_wrap(p4est_ptr)
  p4est.user_pointer = Base.unsafe_convert(Ptr{Nothing},flags)
  trees_ptr = p4est.trees
  trees = unsafe_wrap(trees_ptr)
  ntrees = Int(trees.elem_count)
  elem = 0
  for itree in 1:ntrees
    tree = unsafe_load(Ptr{P4est.p4est_tree_t}(trees.array),itree)
    quadrants =  tree.quadrants
    nquadrants = Int(quadrants.elem_count)
    for iquadrant in 1:nquadrants
      elem += 1
      quadrant = unsafe_load(Ptr{P4est.p4est_quadrant_t}(quadrants.array),iquadrant)
      elem_ptr = Base.unsafe_convert(Ptr{Int},quadrant.p.user_data)
      unsafe_store!(elem_ptr,elem)
    end
  end
  refine_fn_ptr = @cfunction(_refine_fn,Cint,
    (Ptr{P4est.p4est_t},P4est.p4est_topidx_t,Ptr{P4est.p4est_quadrant_t}))
  init_fn_ptr = @cfunction(_init_fn,Cvoid,
    (Ptr{P4est.p4est_t},P4est.p4est_topidx_t,Ptr{P4est.p4est_quadrant_t}))
  for i in 1:num_levels
    P4est.p4est_refine(p4est_ptr,0,refine_fn_ptr,init_fn_ptr)
  end
  amr
end

function balance!(amr)
  @assert domain_dim(amr) == 2 "Not yet implemented for 3d"
  p4est_ptr = amr.p4est_ptr
  init_fn_ptr = @cfunction(_init_fn,Cvoid,
    (Ptr{P4est.p4est_t},P4est.p4est_topidx_t,Ptr{P4est.p4est_quadrant_t}))
  P4est.p4est_balance(p4est_ptr,P4est.P4EST_CONNECT_CORNER,init_fn_ptr)
  amr
end

function coarsen!(amr::P4estMeshRefiner,_flags;num_levels=1)
  flags = convert(Vector{Bool},_flags)
  D = domain_dim(amr)
  if D == 3
    error("Not implemented yet")
  end
  p4est_ptr = amr.p4est_ptr
  p4est = unsafe_wrap(p4est_ptr)
  p4est.user_pointer = Base.unsafe_convert(Ptr{Nothing},flags)
  trees_ptr = p4est.trees
  trees = unsafe_wrap(trees_ptr)
  ntrees = Int(trees.elem_count)
  elem = 0
  for itree in 1:ntrees
    tree = unsafe_load(Ptr{P4est.p4est_tree_t}(trees.array),itree)
    quadrants =  tree.quadrants
    nquadrants = Int(quadrants.elem_count)
    for iquadrant in 1:nquadrants
      elem += 1
      quadrant = unsafe_load(Ptr{P4est.p4est_quadrant_t}(quadrants.array),iquadrant)
      elem_ptr = Base.unsafe_convert(Ptr{Int},quadrant.p.user_data)
      unsafe_store!(elem_ptr,elem)
    end
  end
  coarsen_fn_ptr = @cfunction(_coarsen_fn,Cint,
    (Ptr{P4est.p4est_t},P4est.p4est_topidx_t,Ptr{Ptr{P4est.p4est_quadrant_t}}))
  init_fn_ptr = @cfunction(_init_fn,Cvoid,
    (Ptr{P4est.p4est_t},P4est.p4est_topidx_t,Ptr{P4est.p4est_quadrant_t}))
  for i in 1:num_levels
    P4est.p4est_coarsen(p4est_ptr,0,coarsen_fn_ptr,init_fn_ptr)
  end
  amr
end

function p4est_lnodes_decode!(hanging_corner,face_code::P4est.p4est_lnodes_code_t)
  # This is for 2D
  # Copied from p4est_step4.c without understanding the logic
  fill!(hanging_corner,-1)
  _ones = Cint(P4est.P4EST_CHILDREN - 1)
  if face_code != Int8(0)
    c = Cint( face_code & _ones)
    work = Cint(face_code >> P4est.P4EST_DIM)
    hanging_corner[c+1] = -1
    hanging_corner[xor(c,_ones)+1] = -1
    for i in Cint(1):Cint(P4est.P4EST_DIM)
      h = xor(c,Cint(1) << Cint(i-1))
      hanging_corner[xor(h,_ones)+1] = (work & Cint(1))!=0 ? c+1 : -1
      work = work >> Cint(1)
    end
    true
  else
    false
  end
end

# This is just a proof of concept implementation
# Do not expect it to be efficient in its current state
# Main missing points
# - function barriers
# - extension to 3d
# - split long function
# - creation of a polytopal complex directly without using an intermediate fe mesh
#   by trying to take more advantadge of p4est queries.
function fe_mesh(amr::P4estMeshRefiner)
  balance!(amr)
  p4est_ptr = amr.p4est_ptr
  D = domain_dim(amr)
  @assert D == 2 "Only implemented for 2d"
  Da = ambient_dim(amr)
  order = 1
  nchildren = Int(P4est.P4EST_CHILDREN)
  ghost_ptr = P4est.p4est_ghost_new(p4est_ptr,P4est.P4EST_CONNECT_FULL)
  lnodes_ptr = P4est.p4est_lnodes_new(p4est_ptr,ghost_ptr,order)
  lnodes = unsafe_wrap(lnodes_ptr)
  nelems = Int(lnodes.num_local_elements)
  elem_to_nodes_mat = unsafe_wrap(
    Array,lnodes.element_nodes,(nchildren,nelems))
  ptrs = zeros(Int32,nelems+1)
  for elem in 1:nelems
    ptrs[elem+1] = nchildren
  end
  prefix!(ptrs)
  ndata = ptrs[end]-1
  data = zeros(Int32,ndata)
  for i in 1:length(elem_to_nodes_mat)
    data[i] = Int32(elem_to_nodes_mat[i])+Int32(1)
  end
  elem_to_nodes = JaggedArray(data,ptrs)
  nnodes = lnodes.num_local_nodes
  p4est = unsafe_wrap(p4est_ptr)
  trees_ptr = p4est.trees
  trees = unsafe_wrap(trees_ptr)
  ntrees = Int(trees.elem_count)
  hanging_corner = zeros(Int,nchildren)
  face_code_ptr = lnodes.face_code
  face_code = unsafe_wrap(Array,face_code_ptr,(nelems,))
  ihang = 0
  elem = 0
  elem_to_tree = zeros(Int32,nelems)
  for itree in 1:ntrees
    tree = unsafe_load(Ptr{P4est.p4est_tree_t}(trees.array),itree)
    quadrants =  tree.quadrants
    nquadrants = Int(quadrants.elem_count)
    for iquadrant in 1:nquadrants
      elem += 1
      elem_to_tree[elem] = itree
      quadrant = unsafe_load(Ptr{P4est.p4est_quadrant_t}(quadrants.array),iquadrant)
      level = Int(quadrant.level)
      q0 = SVector(quadrant.x,quadrant.y)
      jvertex = 0
      anyhang = p4est_lnodes_decode!(hanging_corner,face_code[elem])
      for icorner = 1:nchildren
        if anyhang && hanging_corner[icorner] != -1
          ihang += 1
        end
      end
    end
  end
  maxindeps = 2
  hang_to_indeps = fill(Int32(INVALID),maxindeps,ihang)
  hang_to_corner_and_elem = zeros(Int32,2,ihang)
  ihang = 0
  elem = 0
  for itree in 1:ntrees
    tree = unsafe_load(Ptr{P4est.p4est_tree_t}(trees.array),itree)
    quadrants =  tree.quadrants
    nquadrants = Int(quadrants.elem_count)
    for iquadrant in 1:nquadrants
      elem += 1
      nodes = elem_to_nodes[elem]
      quadrant = unsafe_load(Ptr{P4est.p4est_quadrant_t}(quadrants.array),iquadrant)
      level = Int(quadrant.level)
      q0 = SVector(quadrant.x,quadrant.y)
      jvertex = 0
      anyhang = p4est_lnodes_decode!(hanging_corner,face_code[elem])
      for icorner = 1:nchildren
        if anyhang && hanging_corner[icorner] != -1
          indep1 = nodes[icorner]
          indep2 = nodes[hanging_corner[icorner]]
          ihang += 1
          hang_to_indeps[1,ihang] = indep1
          hang_to_indeps[2,ihang] = indep2
          hang_to_corner_and_elem[1,ihang] = icorner
          hang_to_corner_and_elem[2,ihang] = elem
        end
      end
    end
  end
  ptrs = zeros(Int32,nnodes+1)
  for indeps in hang_to_indeps
    for indep in indeps
      ptrs[indep+1] += Int32(1)
    end
  end
  maxhang = maximum(view(ptrs,2:length(ptrs)))
  hangs1 = fill(Int32(INVALID),maxhang)
  hangs2 = fill(Int32(INVALID),maxhang)
  prefix!(ptrs)
  ndata = ptrs[end]-1
  data = zeros(Int32,ndata)
  for hang in 1:size(hang_to_indeps,2)
    for k in 1:size(hang_to_indeps,1)
      indep = hang_to_indeps[k,hang]
      p = ptrs[indep]
      data[p] = hang
      ptrs[indep] += Int32(1)
    end
  end
  rewind!(ptrs)
  indep_to_hang = JaggedArray(data,ptrs)
  jhang = nnodes
  elem = 0
  elem_touch = fill(false,(nchildren,nelems))
  hang_to_finalhang = fill(Int32(INVALID),ihang)
  for itree in 1:ntrees
    tree = unsafe_load(Ptr{P4est.p4est_tree_t}(trees.array),itree)
    quadrants =  tree.quadrants
    nquadrants = Int(quadrants.elem_count)
    for iquadrant in 1:nquadrants
      elem += 1
      nodes = elem_to_nodes[elem]
      quadrant = unsafe_load(Ptr{P4est.p4est_quadrant_t}(quadrants.array),iquadrant)
      level = Int(quadrant.level)
      q0 = SVector(quadrant.x,quadrant.y)
      jvertex = 0
      anyhang = p4est_lnodes_decode!(hanging_corner,face_code[elem])
      for icorner = 1:nchildren
        if anyhang && hanging_corner[icorner] != -1
          if elem_touch[icorner,elem]
            continue
          end
          indep1 = nodes[icorner]
          indep2 = nodes[hanging_corner[icorner]]
          jhang += 1
          fill!(hangs1,Int32(INVALID))
          fill!(hangs2,Int32(INVALID))
          copyto!(hangs1,indep_to_hang[indep1])
          nhangs1 = length(indep_to_hang[indep1])
          copyto!(hangs2,indep_to_hang[indep2])
          nhangs2 = length(indep_to_hang[indep2])
          _intersect_ids!(hangs1,nhangs1,hangs2,nhangs2)
          for ihang in hangs1
            if ihang != Int32(INVALID)
              hang_to_finalhang[ihang] = jhang
              break
            end
          end
          for ihang in hangs1
            if ihang != Int32(INVALID)
              jcorner = hang_to_corner_and_elem[1,ihang]
              jelem = hang_to_corner_and_elem[2,ihang]
              elem_touch[jcorner,jelem] = true
              elem_to_nodes[jelem][jcorner] = jhang
            end
          end
        end
      end
    end
  end
  ids = findall(i->i!=Int32(INVALID),hang_to_finalhang)
  hanging_ids = hang_to_finalhang[ids]
  ptrs = zeros(Int32,length(hanging_ids)+1)
  for (i,id) in enumerate(hanging_ids)
    for j in 1:maxindeps
      indep = hang_to_indeps[j,ids[i]]
      if indep != Int32(INVALID)
        ptrs[i+1] += Int32(1)
      end
    end
  end
  prefix!(ptrs)
  ndata = ptrs[end]-1
  data = zeros(Int32,ndata)
  for (i,id) in enumerate(hanging_ids)
    p = ptrs[i]-Int32(1)
    for j in 1:maxindeps
      indep = hang_to_indeps[j,ids[i]]
      if indep != Int32(INVALID)
        data[p+j] = indep
      end
    end
  end
  hanging_indeps = JaggedArray(data,ptrs)
  data2 = fill(0.5,ndata)
  hanging_coeffs = JaggedArray(data2,ptrs)
  conn_ptr = p4est.connectivity
  conn = unsafe_wrap(conn_ptr)
  nvertices = Int(conn.num_vertices)
  v_ptr = conn.vertices
  v = unsafe_wrap(Array,v_ptr,(3,nvertices,))
  T = eltype(v)
  vertex_to_coords = [ SVector(v[1,i],v[2,i]) for i in 1:nvertices ]
  tree_to_vertex_ptr = conn.tree_to_vertex
  tree_to_vertices = unsafe_wrap(Array,tree_to_vertex_ptr,(nchildren,ntrees))
  tree_to_vertices = copy(tree_to_vertices)
  tree_to_vertices .+= 1
  node_to_coordinates = zeros(SVector{Da,Float64},jhang)
  s = zeros(SVector{2,Float64},2)
  elem = 0
  for itree in 1:ntrees
    tree = unsafe_load(Ptr{P4est.p4est_tree_t}(trees.array),itree)
    quadrants =  tree.quadrants
    nquadrants = Int(quadrants.elem_count)
    vertices = tree_to_vertices[:,itree]
    coords = vertex_to_coords[vertices]
    for iquadrant in 1:nquadrants
      elem += 1
      nodes = elem_to_nodes[elem]
      quadrant = unsafe_load(Ptr{P4est.p4est_quadrant_t}(quadrants.array),iquadrant)
      level = Int(quadrant.level)
      q0 = SVector(quadrant.x,quadrant.y)
      jvertex = 0
      for jy in 0:1
        for jx in 0:1
          jvertex += 1
          node = nodes[jvertex]
          q = q0 + 2^(P4est.P4EST_MAXLEVEL-level)*SVector(jx,jy)
          s[2] = q ./ 2^P4est.P4EST_MAXLEVEL
          s[1] = 1. .- s[2]
          ivertex = 0
          x = zero(SVector{2,Float64})
          for iy in 1:2
            yfactor = s[iy][2]
            for ix in 1:2
              ivertex += 1
              xfactor = yfactor*s[ix][1]
              x += xfactor*coords[ivertex]
            end
          end
          node_to_coordinates[node] = x
        end
      end
    end
  end
  P4est.p4est_lnodes_destroy(lnodes_ptr)
  elem_to_refid = fill(Int8(1),length(elem_to_nodes))
  refid_to_refelem = [P4estQuad{Float64}()]
  mesh = GenericFEMesh{SVector{Da,Float64}}(VOID,D)
  node_coordinates!(mesh,node_to_coordinates)
  face_nodes!(mesh,elem_to_nodes,D)
  face_ref_id!(mesh,elem_to_refid,D)
  ref_faces!(mesh,refid_to_refelem,D)
  hanging_nodes!(mesh,(hanging_ids,hanging_indeps,hanging_coeffs))
  cpoly = polytopal_complex(amr.coarse_mesh)
  fpoly = polytopal_complex(mesh)
  groups_cpoly = physical_groups(cpoly)
  groups_fpoly = physical_groups(fpoly)
  groups_mesh = physical_groups(mesh)
  fcell_to_ccell = elem_to_tree
  max_coord = 2^P4est.P4EST_MAXLEVEL
  min_coord = 0*max_coord
  refcell = first(ref_faces(amr.coarse_mesh,D))
  for d in 0:D
    ncfaces = num_faces(cpoly,d)
    nffaces = num_faces(fpoly,d)
    fcell_to_ffaces = face_incidence(fpoly,D,d)
    ccell_to_cfaces = face_incidence(cpoly,D,d)
    lfface_to_lcface = p4est_face_permutation(refcell,d)
    @boundscheck @assert D == 2 "Not implemented for 3d yet"
    if d==0
      fface_to_cface = fill(Int32(INVALID),nffaces)
      for (fcell,quadrant) in enumerate(p4est_quadrants(p4est_ptr))
        for lfface in 1:4
          neigh = p4est_corner_neighbor(quadrant,lfface)
          lcface = lfface_to_lcface[lfface]
          if (neigh.x < min_coord && neigh.y < min_coord) ||
            (neigh.x < min_coord && neigh.y >= max_coord) ||
            (neigh.x >= max_coord && neigh.y < min_coord) ||
            (neigh.x >= max_coord && neigh.y >= max_coord)
            ccell = fcell_to_ccell[fcell]
            cface = ccell_to_cfaces[ccell][lcface]
            fface = fcell_to_ffaces[fcell][lfface]
            fface_to_cface[fface] = cface
          end
        end
      end
    elseif d==1
      fface_to_cface = fill(Int32(INVALID),nffaces)
      for (fcell,quadrant) in enumerate(p4est_quadrants(p4est_ptr))
        for lfface in 1:4
          neigh = p4est_face_neighbor(quadrant,lfface)
          lcface = lfface_to_lcface[lfface]
          if neigh.x < min_coord || neigh.x>=max_coord || neigh.y < min_coord || neigh.y>=max_coord
            ccell = fcell_to_ccell[fcell]
            cface = ccell_to_cfaces[ccell][lcface]
            fface = fcell_to_ffaces[fcell][lfface]
            fface_to_cface[fface] = cface
          end
        end
      end
    else
      fface_to_cface = fcell_to_ccell
    end
    cface_to_mask = fill(false,ncfaces)
    fface_to_mask = fill(false,nffaces)
    for id in group_ids(groups_cpoly,d)
      fill!(cface_to_mask,false)
      cfaces_in_group = group_faces(groups_cpoly,d,id)
      for cface in cfaces_in_group
        cface_to_mask[cface] = true
      end
      fill!(fface_to_mask,false)
      for fface in 1:nffaces
        cface = fface_to_cface[fface]
        if cface != Int32(INVALID)
          fface_to_mask[fface] = cface_to_mask[cface]
        end
      end
      ffaces_in_group = collect(Int32,findall(fface_to_mask))
      name = group_name(groups_cpoly,d,id)
      add_group!(groups_fpoly,d,name,id)
      group_faces!(groups_fpoly,ffaces_in_group,d,id)
    end
    if d!=D
      fill!(fface_to_mask,false)
      for id in group_ids(groups_fpoly,d)
        ffaces_in_group = group_faces(groups_fpoly,d,id)
        for fface in ffaces_in_group
          fface_to_mask[fface] = true
        end
      end
      fface_to_nodes = face_nodes(fpoly,d)
      mesh_fface_to_fface = collect(Int32,findall(fface_to_mask))
      fface_to_mesh_fface = fill(Int32(INVALID),nffaces)
      n_mesh_ffaces = length(mesh_fface_to_fface)
      ptrs = zeros(Int32,n_mesh_ffaces+1)
      for (mesh_fface,fface) in enumerate(mesh_fface_to_fface)
        ptrs[mesh_fface+1] = length(fface_to_nodes[fface])
        fface_to_mesh_fface[fface] = mesh_fface
      end
      prefix!(ptrs)
      ndata = ptrs[end]-1
      data = zeros(Int32,ndata)
      for (mesh_fface,fface) in enumerate(mesh_fface_to_fface)
        nodes = fface_to_nodes[fface]
        p = ptrs[mesh_fface]-1
        for (lnode,node) in enumerate(nodes)
          data[p+lnode] = node
        end
      end
      mesh_fface_to_nodes = JaggedArray(data,ptrs)
      face_nodes!(mesh,mesh_fface_to_nodes,d)
      face_ref_id!(mesh,fill(Int8(1),length(mesh_fface_to_nodes)),d)
      ref_faces!(mesh,ref_faces(P4estQuad{Float64}(),d),d)
      mesh_faces!(fpoly,mesh_fface_to_fface,d)
      for id in group_ids(groups_fpoly,d)
        ffaces_in_group = group_faces(groups_fpoly,d,id)
        mesh_ffaces_in_group = fface_to_mesh_fface[ffaces_in_group]
        name = group_name(groups_fpoly,d,id)
        add_group!(groups_mesh,d,name,id)
        group_faces!(groups_mesh,mesh_ffaces_in_group,d,id)
      end
    else
      mesh_fface_to_fface = collect(Int32,1:nffaces)
      mesh_faces!(fpoly,mesh_fface_to_fface,d)
      for id in group_ids(groups_fpoly,d)
        ffaces_in_group = group_faces(groups_fpoly,d,id)
        mesh_ffaces_in_group = ffaces_in_group
        name = group_name(groups_fpoly,d,id)
        add_group!(groups_mesh,d,name,id)
        group_faces!(groups_mesh,mesh_ffaces_in_group,d,id)
      end
    end
  end
  mesh
end



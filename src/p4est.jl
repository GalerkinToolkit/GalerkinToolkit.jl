
_p4est_vertex_perm = [1,2,4,3]
p4est_vertex_perm(q::Meshes.Quadrangle) = _p4est_vertex_perm

struct P4estQuad end
function vtk_mesh_cell(a::P4estQuad)
  nodes -> WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_QUAD,nodes[_p4est_vertex_perm])
end

mutable struct P4estAMR{A,B}
  p4est_ptr::Ptr{A}
  coarse_mesh::B
  state::Symbol
end

domain_dim(a::P4estAMR) = domain_dim(a.coarse_mesh)
ambient_dim(a::P4estAMR) = ambient_dim(a.coarse_mesh)

function p4est_amr(geo;initial_level=0)
  if !MPI.Initialized()
    MPI.Init()
  end
  mesh = fe_mesh(geo)
  D = domain_dim(mesh)
  @assert num_faces(mesh,D) == 1 "For the moment connectivity creation only for a single tree"
  refface = first(ref_faces(mesh,D))
  perm = p4est_vertex_perm(refface)
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
    conn_ptr = P4est.p4est_connectivity_new_unitsquare()
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
  amr = P4estAMR(p4est_ptr,mesh,:usable)
  finalizer(destroy!,amr)
end

function destroy!(a::P4estAMR)
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

function _set_all_quadrants_data!(p4est_ptr)
  p4est = unsafe_wrap(p4est_ptr)
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
end

function _fine_to_coarse_after_refinement(p4est_ptr)
  p4est = unsafe_wrap(p4est_ptr)
  trees_ptr = p4est.trees
  trees = unsafe_wrap(trees_ptr)
  ntrees = Int(trees.elem_count)
  elem = 0
  for itree in 1:ntrees
    tree = unsafe_load(Ptr{P4est.p4est_tree_t}(trees.array),itree)
    quadrants =  tree.quadrants
    nquadrants = Int(quadrants.elem_count)
    elem += nquadrants
  end
  new_to_old = zeros(Int32,elem)
  new = 0
  old = 0
  ichild = 0
  for itree in 1:ntrees
    tree = unsafe_load(Ptr{P4est.p4est_tree_t}(trees.array),itree)
    quadrants =  tree.quadrants
    nquadrants = Int(quadrants.elem_count)
    for iquadrant in 1:nquadrants
      new += 1
      quadrant = unsafe_load(Ptr{P4est.p4est_quadrant_t}(quadrants.array),iquadrant)
      elem_ptr = Base.unsafe_convert(Ptr{Int},quadrant.p.user_data)
      elem = unsafe_load(elem_ptr)
      if elem == INVALID
        ichild +=1
        if ichild == 1
          old += 1
        end
      else
        ichild = 0
        old +=1
      end
      new_to_old[new] = old
    end
  end
  new_to_old
end

function _init_fn(p4est_ptr,which_tree,quadrant_ptr)
  quadrant = unsafe_wrap(quadrant_ptr)
  elem_ptr = Base.unsafe_convert(Ptr{Int},quadrant.p.user_data)
  unsafe_store!(elem_ptr,INVALID)
  Cvoid()
end
const init_fn_ptr = @cfunction(_init_fn,Cvoid,
  (Ptr{P4est.p4est_t},P4est.p4est_topidx_t,Ptr{P4est.p4est_quadrant_t}))

function _refine_fn(p4est_ptr,which_tree,quadrant_ptr)
  p4est = unsafe_wrap(p4est_ptr)
  quadrant = unsafe_wrap(quadrant_ptr)
  flags_ptr = Base.unsafe_convert(Ptr{Bool},p4est.user_pointer)
  elem_ptr = Base.unsafe_convert(Ptr{Int},quadrant.p.user_data)
  elem = unsafe_load(elem_ptr,1)
  flag = unsafe_load(flags_ptr,elem)
  i = flag ? 1 : 0
  Int32(i)
end
const refine_fn_ptr = @cfunction(_refine_fn,Cint,
  (Ptr{P4est.p4est_t},P4est.p4est_topidx_t,Ptr{P4est.p4est_quadrant_t}))

function refine!(amr::P4estAMR,_flags)
  flags = convert(Vector{Bool},_flags)
  D = domain_dim(amr)
  if D == 3
    error("Not implemented yet")
  end
  p4est_ptr = amr.p4est_ptr
  p4est = unsafe_wrap(p4est_ptr)
  p4est.user_pointer = Base.unsafe_convert(Ptr{Nothing},flags)
  _set_all_quadrants_data!(p4est_ptr)
  P4est.p4est_refine(p4est_ptr,0,refine_fn_ptr,init_fn_ptr)
  _fine_to_coarse_after_refinement(p4est_ptr)
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
    flag = unsafe_load(flags_ptr,elem)
    r = r && (!flag)
  end
  r ? Int32(1) : Int32(0)
end
const coarsen_fn_ptr = @cfunction(_coarsen_fn,Cint,
    (Ptr{P4est.p4est_t},P4est.p4est_topidx_t,Ptr{Ptr{P4est.p4est_quadrant_t}}))

function coarsen!(amr::P4estAMR,_flags)
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
  P4est.p4est_coarsen(p4est_ptr,0,coarsen_fn_ptr,C_NULL)
  amr
end

function balance!(amr)
  @assert domain_dim(amr) == 2 "Not yet implemented for 3d"
  p4est_ptr = amr.p4est_ptr
  _set_all_quadrants_data!(p4est_ptr)
  P4est.p4est_balance(p4est_ptr,P4est.P4EST_CONNECT_CORNER,init_fn_ptr)
  _fine_to_coarse_after_refinement(p4est_ptr)
end

function _lnodes_decode!(hanging_corner,face_code::P4est.p4est_lnodes_code_t)
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

function fe_mesh(amr::P4estAMR)
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
  for itree in 1:ntrees
    tree = unsafe_load(Ptr{P4est.p4est_tree_t}(trees.array),itree)
    quadrants =  tree.quadrants
    nquadrants = Int(quadrants.elem_count)
    for iquadrant in 1:nquadrants
      elem += 1
      quadrant = unsafe_load(Ptr{P4est.p4est_quadrant_t}(quadrants.array),iquadrant)
      level = Int(quadrant.level)
      q0 = SVector(quadrant.x,quadrant.y)
      jvertex = 0
      anyhang = _lnodes_decode!(hanging_corner,face_code[elem])
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
      anyhang = _lnodes_decode!(hanging_corner,face_code[elem])
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
      anyhang = _lnodes_decode!(hanging_corner,face_code[elem])
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
  refid_to_refelem = [P4estQuad()]
  mesh = SimpleFEMesh{SVector{Da,Float64}}(VOID,D)
  node_coordinates!(mesh,node_to_coordinates)
  face_nodes!(mesh,D,elem_to_nodes)
  face_ref_id!(mesh,D,elem_to_refid)
  ref_faces!(mesh,D,refid_to_refelem)
  hanging_nodes!(mesh,(hanging_ids,hanging_indeps,hanging_coeffs))
  mesh
end



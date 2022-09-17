
function fe_mesh(mesh::Meshes.Mesh;kwargs...)
  _fe_mesh_from_mesh(mesh;kwargs...)
end

function _fe_mesh_from_mesh(
  mesh::Meshes.Mesh;
  physical_groups=nothing,
  hanging_nodes=nothing,
  periodic_nodes=nothing)

  D = num_dims(mesh)
  fe_mesh = FEMesh{num_ambient_dims(mesh)}(D)
  node_coordinates!(fe_mesh,node_coordinates(mesh))
  for d in 0:D
    face_nodes!(fe_mesh,face_nodes(mesh,d),d)
    face_to_refid, refid_to_refface = _meshes_setup_ref_faces(mesh,d)
    face_ref_id!(fe_mesh,face_to_refid,d)
    ref_faces!(fe_mesh,refid_to_refface,d)
  end
  if physical_groups !== nothing
    physical_groups!(fe_mesh,physical_groups)
  end
  if hanging_nodes !== nothing
    hanging_nodes!(fe_mesh,hanging_nodes)
  end
  if periodic_nodes !== nothing
    periodic_nodes!(fe_mesh,periodic_nodes)
  end
  fe_mesh
end

function _meshes_setup_ref_faces(mesh,d)
  if isconcretetype(eltype(mesh))
    face_ref_id(mesh,d), ref_faces(mesh,d)
  else
    if num_dims(mesh) == d
      # This is inefficient
      i_to_rf = []
      ct_to_i = Dict{WriteVTK.VTKCellTypes.VTKCellType,Int}()
      iface_to_i = zeros(Int8,length(mesh))
      i = 0
      for (iface,face) in enumerate(mesh)
        @boundscheck @assert num_dims(face)==d "This case is not implemented"
        rf = first(ref_faces(face,num_dims(face)))
        ct = vtk_cell_type(rf)
        if !haskey(ct_to_i,ct)
          i += 1
          ct_to_i[ct] = i
          push!(i_to_rf,rf)
        end
        iface_to_i[iface] = ct_to_i[ct]
      end
      iface_to_i, i_to_rf
    else
      fill(Int8(1),0), []
    end
  end
end

num_dims(::Type{<:Meshes.Ngon}) = 2
num_dims(a::Meshes.Mesh) = num_dims(eltype(a))
num_ambient_dims(a::Meshes.Mesh) = Meshes.embeddim(a)
function face_ref_id(a::Meshes.Mesh,d)
  @assert isconcretetype(eltype(a))
  if d == num_dims(a)
    fill(Int8(1),length(a))
  else
    fill(Int8(1),0)
  end
end
function ref_faces(a::Meshes.Mesh,d)
  @assert isconcretetype(eltype(a))
  if d == num_dims(a)
    ref_faces(eltype(a),d)
  else
    []
  end
end
function face_nodes(a::Meshes.Mesh,d)
  if d == num_dims(a)
    vs = convert(Meshes.FullTopology,Meshes.topology(a)).connec
    j = map(vs) do v
      convert(SVector{length(v.indices),Int32}, v.indices)
    end
    JaggedArray(j)
  else
    JaggedArray(Vector{Int32}[])
  end
end
function node_coordinates(a::Meshes.Mesh)
  collect(Meshes.coordinates.(Meshes.vertices(a)))
end

const _MESHES_BUFFER = Dict{Symbol,Any}()

# Point
num_dims(a::Meshes.Point) = 0
num_dims(a::Type{<:Meshes.Point}) = 0
num_ambient_dims(a::Meshes.Point) = Meshes.embeddim(a)
is_simplex(a::Meshes.Point) = true
is_hypercube(a::Meshes.Point) = true
function face_ref_id(a::Meshes.Point,d)
  d!=0 && throw(DomainError(d))
  [Int8(1)]
end
function ref_faces(::Type{<:Meshes.Point},d)
  d!=0 && throw(DomainError(d))
  [Meshes.Point(SVector{0,Float64}())]
end
ref_faces(a::Meshes.Point,d) = ref_faces(typeof(a),d)
function face_faces(a::Meshes.Point,d1,d2)
  d1!=0 && throw(DomainError(d1))
  d2!=0 && throw(DomainError(d2))
  JaggedArray([[Int32(1)]])
end
face_nodes(a::Meshes.Point,d) = face_faces(a,d,0)
#face_own_nodes(a::Meshes.Point,d) = face_nodes(a,d)
node_coordinates(a::Meshes.Point) = [Meshes.coordinates(a)]
function vtk_mesh_cell(a::Meshes.Point)
  nodes -> WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_VERTEX,nodes)
end
vtk_cell_type(a::Meshes.Point) = WriteVTK.VTKCellTypes.VTK_VERTEX

# Segment
num_dims(a::Meshes.Segment) = 1
num_dims(a::Type{<:Meshes.Segment}) = 1
num_ambient_dims(a::Meshes.Segment) = Meshes.embeddim(a)
is_simplex(a::Meshes.Segment) = true
is_hypercube(a::Meshes.Segment) = true
function face_ref_id(a::Meshes.Segment,d)
  d==0 && return [Int8(1),Int8(1)]
  d==1 && return [Int8(1)]
  throw(DomainError(d))
end
function ref_faces(::Type{<:Meshes.Segment},d)
  d==0 && return [Meshes.Point(SVector{0,Float64}())]
  d==1 && return [Meshes.Segment(Meshes.Point(0),Meshes.Point(1))]
  throw(DomainError(d))
end
ref_faces(a::Meshes.Segment,d) = ref_faces(typeof(a),d)
function face_faces(a::Meshes.Segment,d1,d2)
  (d1==0 && d2==0) && return JaggedArray([[Int32(1)],[Int32(2)]])
  (d1==1 && d2==0) && return JaggedArray([[Int32(1)],[Int32(2)]])
  (d1==0 && d2==1) && return JaggedArray([[Int32(1)],[Int32(1)]])
  (d1==1 && d2==1) && return JaggedArray([[Int32(1)],[Int32(2)]])
  throw(DomainError((d1,d2)))
end
face_nodes(a::Meshes.Segment,d) = face_faces(a,d,0)
#face_own_nodes(a::Meshes.Segment,d) = face_nodes(a,d)
node_coordinates(a::Meshes.Segment) = collect(Meshes.coordinates.(Meshes.vertices(a)))
function vtk_mesh_cell(a::Meshes.Segment)
  nodes -> WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_LINE,nodes)
end
vtk_cell_type(a::Meshes.Segment) = WriteVTK.VTKCellTypes.VTK_LINE

# Quadrangle
num_dims(a::Meshes.Quadrangle) = 2
num_dims(a::Type{<:Meshes.Quadrangle}) = 2
num_ambient_dims(a::Meshes.Quadrangle) = Meshes.embeddim(a)
is_hypercube(a::Meshes.Quadrangle) = true
function face_ref_id(a::Meshes.Quadrangle,d)
  d==0 && return fill(Int8(1),4)
  d==1 && return fill(Int8(1),4)
  d==2 && return fill(Int8(1),1)
  throw(DomainError(d))
end
function ref_faces(::Type{<:Meshes.Quadrangle},d)
  d==0 && return [Meshes.Point(SVector{0,Float64}())]
  d==1 && return [Meshes.Segment(Meshes.Point(0),Meshes.Point(1))]
  d==2 && return [Meshes.Quadrangle(Meshes.Point.([(0,0),(1,0),(1,1),(0,1)]))]
  throw(DomainError(d))
end
ref_faces(a::Meshes.Quadrangle,d) = ref_faces(typeof(a),d)
function face_faces(a::Meshes.Quadrangle,d1,d2)
  (d1==0 && d2==0) && return JaggedArray(Vector{Int32}[[1],[2],[3],[4]])
  (d1==1 && d2==0) && return JaggedArray(Vector{Int32}[[1,2],[2,3],[3,4],[4,1]])
  (d1==2 && d2==0) && return JaggedArray(Vector{Int32}[[1,2,3,4]])
  (d1==0 && d2==1) && return JaggedArray(Vector{Int32}[[1,4],[1,2],[2,3],[3,4]])
  (d1==1 && d2==1) && return JaggedArray(Vector{Int32}[[1],[2],[3],[4]])
  (d1==2 && d2==1) && return JaggedArray(Vector{Int32}[[1,2,3,4]])
  (d1==0 && d2==2) && return JaggedArray(Vector{Int32}[[1],[1],[1],[1]])
  (d1==1 && d2==2) && return JaggedArray(Vector{Int32}[[1],[1],[1],[1]])
  (d1==2 && d2==2) && return JaggedArray(Vector{Int32}[[1]])
  throw(DomainError((d1,d2)))
end
face_nodes(a::Meshes.Quadrangle,d) = face_faces(a,d,0)
#face_own_nodes(a::Meshes.Quadrangle,d) = face_nodes(a,d)
node_coordinates(a::Meshes.Quadrangle) = collect(Meshes.coordinates.(Meshes.vertices(a)))
function vtk_mesh_cell(a::Meshes.Quadrangle)
  nodes -> WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_QUAD,nodes)
end
vtk_cell_type(a::Meshes.Quadrangle) = WriteVTK.VTKCellTypes.VTK_QUAD

# Triangle
num_dims(a::Meshes.Triangle) = 2
num_dims(a::Type{<:Meshes.Triangle}) = 2
num_ambient_dims(a::Meshes.Triangle) = Meshes.embeddim(a)
is_simplex(a::Meshes.Triangle) = true
function face_ref_id(a::Meshes.Triangle,d)
  d==0 && return fill(Int8(1),3)
  d==1 && return fill(Int8(1),3)
  d==2 && return fill(Int8(1),1)
  throw(DomainError(d))
end
function ref_faces(::Type{<:Meshes.Triangle},d)
  d==0 && return [Meshes.Point(SVector{0,Float64}())]
  d==1 && return [Meshes.Segment(Meshes.Point(0),Meshes.Point(1))]
  d==2 && return [Meshes.Triangle(Meshes.Point.([(0,0),(1,0),(0,1)]))]
  throw(DomainError(d))
end
ref_faces(a::Meshes.Triangle,d) = ref_faces(typeof(a),d)
function face_faces(a::Meshes.Triangle,d1,d2)
  (d1==0 && d2==0) && return JaggedArray(Vector{Int32}[[1],[2],[3]])
  (d1==1 && d2==0) && return JaggedArray(Vector{Int32}[[1,2],[2,3],[3,1]])
  (d1==2 && d2==0) && return JaggedArray(Vector{Int32}[[1,2,3]])
  (d1==0 && d2==1) && return JaggedArray(Vector{Int32}[[1,3],[1,2],[2,3]])
  (d1==1 && d2==1) && return JaggedArray(Vector{Int32}[[1],[2],[3]])
  (d1==2 && d2==1) && return JaggedArray(Vector{Int32}[[1,2,3]])
  (d1==0 && d2==2) && return JaggedArray(Vector{Int32}[[1],[1],[1]])
  (d1==1 && d2==2) && return JaggedArray(Vector{Int32}[[1],[1],[1]])
  (d1==2 && d2==2) && return JaggedArray(Vector{Int32}[[1]])
  throw(DomainError((d1,d2)))
end
face_nodes(a::Meshes.Triangle,d) = face_faces(a,d,0)
#face_own_nodes(a::Meshes.Triangle,d) = face_nodes(a,d)
node_coordinates(a::Meshes.Triangle) = collect(Meshes.coordinates.(Meshes.vertices(a)))
function vtk_mesh_cell(a::Meshes.Triangle)
  nodes -> WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_TRIANGLE,nodes)
end
vtk_cell_type(a::Meshes.Triangle) = WriteVTK.VTKCellTypes.VTK_TRIANGLE

# Tetrahedron
num_dims(a::Meshes.Tetrahedron) = 3
num_dims(a::Type{<:Meshes.Tetrahedron}) = 3
num_ambient_dims(a::Meshes.Tetrahedron) = Meshes.embeddim(a)
is_simplex(a::Meshes.Tetrahedron) = true
node_coordinates(a::Meshes.Tetrahedron) = collect(Meshes.coordinates.(Meshes.vertices(a)))
function face_ref_id(a::Meshes.Tetrahedron,d)
  d==0 && return fill(Int8(1),4)
  d==1 && return fill(Int8(1),6)
  d==2 && return fill(Int8(1),4)
  d==3 && return fill(Int8(1),1)
  throw(DomainError(d))
end
function ref_faces(::Type{<:Meshes.Tetrahedron},d)
  d==0 && return [Meshes.Point(SVector{0,Float64}())]
  d==1 && return [Meshes.Segment(Meshes.Point(0),Meshes.Point(1))]
  d==2 && return [Meshes.Triangle(Meshes.Point.([(0,0),(1,0),(0,1)]))]
  d==3 && return [Meshes.Tetrahedron(Meshes.Point.([(0,0,0),(1,0,0),(0,1,0),(0,0,1)]))]
  throw(DomainError(d))
end
ref_faces(a::Meshes.Tetrahedron,d) = ref_faces(typeof(a),d)
function face_nodes(a::Meshes.Tetrahedron,d)
  d==0 && return JaggedArray(Vector{Int32}[[1],[2],[3],[4]])
  d==1 && return JaggedArray(Vector{Int32}[[1,2],[2,3],[3,1],[1,4],[2,4],[3,4]])
  d==2 && return JaggedArray(Vector{Int32}[[1,3,2],[1,2,4],[2,3,4],[3,1,4]])
  d==3 && return JaggedArray(Vector{Int32}[[1,2,3,4]])
  throw(DomainError(d))
end
_MESHES_BUFFER[:Tetrahedron] = Dict{Symbol,Any}()
function polytope_boundary(a::Meshes.Tetrahedron)
  if !haskey(_MESHES_BUFFER[:Tetrahedron],:polytope_boundary)
    _MESHES_BUFFER[:Tetrahedron][:polytope_boundary] = default_polytope_boundary(a)
  end
  _MESHES_BUFFER[:Tetrahedron][:polytope_boundary]
end
face_faces(a::Meshes.Tetrahedron,d1,d2) = default_polytope_face_faces(a,d1,d2)
function vtk_mesh_cell(a::Meshes.Tetrahedron)
  nodes -> WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_TETRA,nodes)
end
vtk_cell_type(a::Meshes.Tetrahedron) = WriteVTK.VTKCellTypes.VTK_TERA

# Hexahedron
num_dims(a::Meshes.Hexahedron) = 3
num_dims(a::Type{<:Meshes.Hexahedron}) = 3
num_ambient_dims(a::Meshes.Hexahedron) = Meshes.embeddim(a)
is_hypercube(a::Meshes.Hexahedron) = true
node_coordinates(a::Meshes.Hexahedron) = collect(Meshes.coordinates.(Meshes.vertices(a)))
function face_ref_id(a::Meshes.Hexahedron,d)
  d==0 && return fill(Int8(1),8)
  d==1 && return fill(Int8(1),12)
  d==2 && return fill(Int8(1),6)
  d==3 && return fill(Int8(1),1)
  throw(DomainError(d))
end
function ref_faces(::Type{<:Meshes.Hexahedron},d)
  d==0 && return [Meshes.Point(SVector{0,Float64}())]
  d==1 && return [Meshes.Segment(Meshes.Point(0),Meshes.Point(1))]
  d==2 && return [Meshes.Quadrangle(Meshes.Point.([(0,0),(1,0),(1,1),(0,1)]))]
  d==3 && return [Meshes.Hexahedron(Meshes.Point.([(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)]))]
  throw(DomainError(d))
end
ref_faces(a::Meshes.Hexahedron,d) = ref_faces(typeof(a),d)
function face_nodes(a::Meshes.Hexahedron,d)
  d==0 && return JaggedArray(Vector{Int32}[[1],[2],[3],[4],[5],[6],[7],[8]])
  d==1 && return JaggedArray(Vector{Int32}[[1,2],[2,3],[3,4],[4,1],[5,6],[6,7],[7,8],[8,5],[1,5],[2,6],[3,7],[4,8]])
  d==2 && return JaggedArray(Vector{Int32}[[1,4,3,2],[1,2,6,5],[2,3,7,6],[3,4,8,7],[4,1,5,8],[5,6,7,8]])
  d==3 && return JaggedArray(Vector{Int32}[[1,2,3,4,5,6,7,8]])
  throw(DomainError(d))
end
_MESHES_BUFFER[:Hexahedron] = Dict{Symbol,Any}()
function polytope_boundary(a::Meshes.Hexahedron)
  if !haskey(_MESHES_BUFFER[:Hexahedron],:polytope_boundary)
    _MESHES_BUFFER[:Hexahedron][:polytope_boundary] = default_polytope_boundary(a)
  end
  _MESHES_BUFFER[:Hexahedron][:polytope_boundary]
end
face_faces(a::Meshes.Hexahedron,d1,d2) = default_polytope_face_faces(a,d1,d2)
function vtk_mesh_cell(a::Meshes.Hexahedron)
  nodes -> WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_HEXAHEDRON,nodes)
end
vtk_cell_type(a::Meshes.Hexahedron) = WriteVTK.VTKCellTypes.VTK_HEXAHEDRON

function fe_mesh(
  mesh::Meshes.CartesianGrid;
  is_periodic=map(i->false,mesh.dims))

  fe_mesh = _fe_mesh_from_mesh(mesh)
  groups, faces = _default_physical_groups_cartesian_grid(fe_mesh)
  face_to_nodes, face_to_refid, refid_to_refface = faces
  physical_groups!(fe_mesh,groups)
  D = num_dims(fe_mesh)
  for d in 0:(D-1)
    face_nodes!(fe_mesh,face_to_nodes[d+1],d)
    face_ref_id!(fe_mesh,face_to_refid[d+1],d)
    ref_faces!(fe_mesh,refid_to_refface[d+1],d)
  end
  periodicnodes = _periodic_nodes_cartesian_grid(mesh.dims,is_periodic)
  periodic_nodes!(fe_mesh,periodicnodes)
  fe_mesh
end

function _default_physical_groups_cartesian_grid(fe_mesh)
  D = num_dims(fe_mesh)
  cell_to_nodes = face_nodes(fe_mesh,D)
  refcell = first(ref_faces(fe_mesh,D))
  nnodes = num_nodes(fe_mesh)
  d_to_ldface_to_lnodes = [ face_nodes(refcell,d) for d in 0:(D-1)]
  groups, face_to_nodes = _default_physical_groups_cartesian_grid(
    D,
    cell_to_nodes,
    nnodes,
    d_to_ldface_to_lnodes)
  face_to_refid = [ ones(Int8,length(face_to_nodes[d+1]))  for d in 0:(D-1)]
  refid_to_faces = [ ref_faces(refcell,d)  for d in 0:(D-1)]
  groups, (face_to_nodes,face_to_refid,refid_to_faces)
end

function _default_physical_groups_cartesian_grid(
  D,
  cell_to_nodes,
  nnodes,
  d_to_ldface_to_lnodes)

  node_to_n = zeros(Int32,nnodes)
  for nodes in cell_to_nodes
    for node in nodes
      node_to_n[node] += Int32(1)
    end
  end
  J = typeof(JaggedArray(Vector{Int32}[]))
  face_to_nodes = Vector{J}(undef,D)
  groups = PhysicalGroupCollection(D)
  ngroups = 0
  for d in 0:(D-1)
    nmax = 2^d
    ldface_to_lnodes = d_to_ldface_to_lnodes[d+1]
    ndfaces = 0
    for nodes in cell_to_nodes
      for (ldface,lnodes) in enumerate(ldface_to_lnodes)
        isboundary = true
        for lnode in lnodes
          node = nodes[lnode]
          if node_to_n[node] > nmax
            isboundary = false
            break
          end
        end
        if isboundary
          ndfaces += 1
        end
      end
    end
    ptrs = zeros(Int32,ndfaces+1)
    for dface in 1:ndfaces
      ptrs[dface+1] += Int32(nmax)
    end
    prefix!(ptrs)
    ndata = ptrs[end]-1
    data = zeros(Int32,ndata)
    dface_to_physical_group = zeros(Int32,ndfaces)
    ndfaces = 0
    for nodes in cell_to_nodes
      for (ldface,lnodes) in enumerate(ldface_to_lnodes)
        isboundary = true
        for lnode in lnodes
          node = nodes[lnode]
          if node_to_n[node] > nmax
            isboundary = false
            break
          end
        end
        if isboundary
          ndfaces += 1
          group = ngroups + ldface
          dface_to_physical_group[ndfaces] = group
          p = ptrs[ndfaces]-Int32(1)
          for (i,lnode) in enumerate(lnodes)
            node = nodes[lnode]
            data[p+i] = node
          end
        end
      end
    end
    nldfaces = length(ldface_to_lnodes)
    face_to_nodes[d+1] = JaggedArray(data,ptrs)
    for ldface in 1:nldfaces
      group = ngroups + ldface
      physical_group!(groups,d,"$(d)-face-$ldface",group)
      faces_in_physical_group = findall(g->g==group,dface_to_physical_group)
      physical_group_faces!(groups,faces_in_physical_group,d,group)
    end
    ngroups += nldfaces
  end # d
  ngroups += 1
  physical_group!(groups,D,"$(D)-face-1",ngroups)
  ncells = length(cell_to_nodes)
  physical_group_faces!(groups,collect(Int32,1:ncells),D,ngroups)
  groups, face_to_nodes
end

function _periodic_nodes_cartesian_grid(dim_ncells,is_periodic)
  dim_nnodes = Tuple(dim_ncells) .+ 1
  cis = CartesianIndices(dim_nnodes)
  lis = LinearIndices(dim_nnodes)
  D = length(dim_nnodes)
  nperiodic = 0
  node_touched = fill(false,axes(cis))
  for d in 1:D
    if is_periodic[d]
      for ci in cis
        if node_touched[ci]
          continue
        end
        t = Tuple(ci)
        td = t[d]
        if td == 1
          s = ntuple(i->(i==d ? dim_nnodes[i] : t[i]),Val{D}())
          ci2 = CartesianIndex(s)
          nperiodic += 1
          node_touched[ci] = true
          node_touched[ci2] = true
        end
      end
    end
  end
  periodic_indep = zeros(Int32,nperiodic)
  periodic_dep = zeros(Int32,nperiodic)
  nperiodic = 0
  fill!(node_touched,false)
  for d in 1:D
    if is_periodic[d]
      for ci in cis
        if node_touched[ci]
          continue
        end
        t = Tuple(ci)
        td = t[d]
        if td == 1
          s = ntuple(i->(i==d ? dim_nnodes[i] : t[i]),Val{D}())
          ci2 = CartesianIndex(s)
          nperiodic += 1
          node_touched[ci] = true
          node_touched[ci2] = true
          li = lis[ci]
          li2 = lis[ci2]
          periodic_indep[nperiodic] = li
          periodic_dep[nperiodic] = li2
        end
      end
    end
  end
  PeriodicNodeVector(periodic_dep, periodic_indep, ones(length(periodic_dep)))
end

function polytopal_complex(mesh::Meshes.Mesh)
  polytopal_complex(fe_mesh(mesh))
end




fe_mesh(args...;kwargs...) = default_fe_mesh(args...;kwargs...)

function default_fe_mesh(
  mesh::Meshes.Mesh;
  physical_groups=nothing,
  hanging_nodes=nothing,
  periodic_nodes=nothing)

  D = domain_dim(mesh)
  T = SVector{ambient_dim(mesh),Float64}
  fe_mesh = SimpleFEMesh{T}(void,D)
  node_coordinates!(fe_mesh,node_coordinates(mesh))
  for d in 0:D
    face_nodes!(fe_mesh,d,face_nodes(mesh,d))
    face_ref_id!(fe_mesh,d,face_ref_id(mesh,d))
    ref_faces!(fe_mesh,d,ref_faces(mesh,d))
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

domain_dim(a::Meshes.Mesh) = domain_dim(eltype(a))
ambient_dim(a::Meshes.Mesh) = Meshes.embeddim(a)
function face_ref_id(a::Meshes.Mesh,d)
  @assert isconcretetype(eltype(a))
  if d == domain_dim(a)
    fill(Int8(1),length(a))
  else
    fill(Int8(1),0)
  end
end
function ref_faces(a::Meshes.Mesh,d)
  @assert isconcretetype(eltype(a))
  if d == domain_dim(a)
    ref_faces(eltype(a),d)
  else
    []
  end
end
function face_nodes(a::Meshes.Mesh,d)
  if d == domain_dim(a)
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

# Point
domain_dim(a::Meshes.Point) = 0
domain_dim(a::Type{<:Meshes.Point}) = 0
ambient_dim(a::Meshes.Point) = Meshes.embeddim(a)
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
function face_incidence(a::Meshes.Point,d1,d2)
  d1!=0 && throw(DomainError(d1))
  d2!=0 && throw(DomainError(d2))
  JaggedArray([[Int32(1)]])
end
face_nodes(a::Meshes.Point,d) = face_incidence(a,d,0)
face_own_nodes(a::Meshes.Point,d) = face_nodes(a,d)
node_coordinates(a::Meshes.Point) = [Meshes.coordinates(a)]

# Segment
domain_dim(a::Meshes.Segment) = 1
domain_dim(a::Type{<:Meshes.Segment}) = 1
ambient_dim(a::Meshes.Segment) = Meshes.embeddim(a)
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
function face_incidence(a::Meshes.Segment,d1,d2)
  (d1==0 && d2==0) && return JaggedArray([[Int32(1)],[Int32(2)]])
  (d1==1 && d2==0) && return JaggedArray([[Int32(1)],[Int32(2)]])
  (d1==0 && d2==1) && return JaggedArray([[Int32(1)],[Int32(1)]])
  (d1==1 && d2==1) && return JaggedArray([[Int32(1)],[Int32(2)]])
  throw(DomainError((d1,d2)))
end
face_nodes(a::Meshes.Segment,d) = face_incidence(a,d,0)
face_own_nodes(a::Meshes.Segment,d) = face_nodes(a,d)
node_coordinates(a::Meshes.Segment) = collect(Meshes.coordinates.(Meshes.vertices(a)))

# Quadrangle
domain_dim(a::Meshes.Quadrangle) = 2
domain_dim(a::Type{<:Meshes.Quadrangle}) = 2
ambient_dim(a::Meshes.Quadrangle) = Meshes.embeddim(a)
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
function face_incidence(a::Meshes.Quadrangle,d1,d2)
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
face_nodes(a::Meshes.Quadrangle,d) = face_incidence(a,d,0)
face_own_nodes(a::Meshes.Quadrangle,d) = face_nodes(a,d)
node_coordinates(a::Meshes.Quadrangle) = collect(Meshes.coordinates.(Meshes.vertices(a)))

# Triangle
domain_dim(a::Meshes.Triangle) = 2
domain_dim(a::Type{<:Meshes.Triangle}) = 2
ambient_dim(a::Meshes.Triangle) = Meshes.embeddim(a)
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
function face_incidence(a::Meshes.Triangle,d1,d2)
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
face_nodes(a::Meshes.Triangle,d) = face_incidence(a,d,0)
face_own_nodes(a::Meshes.Triangle,d) = face_nodes(a,d)
node_coordinates(a::Meshes.Triangle) = collect(Meshes.coordinates.(Meshes.vertices(a)))

function fe_mesh(
  mesh::Meshes.CartesianGrid;
  is_periodic=fill(false,domain_dim(mesh)))

  any(is_periodic) && error("Periodic BCs not yer implemented")
  fe_mesh = default_fe_mesh(mesh,periodic_nodes=default_periodic_nodes(mesh))
  groups, faces = default_groups_cartesian_grid(fe_mesh)
  face_to_nodes, face_to_refid, refid_to_refface = faces
  physical_groups!(fe_mesh,groups)
  D = domain_dim(fe_mesh)
  for d in 0:(D-1)
    face_nodes!(fe_mesh,d,face_to_nodes[d+1])
    face_ref_id!(fe_mesh,d,face_to_refid[d+1])
    ref_faces!(fe_mesh,d,refid_to_refface[d+1])
  end
  fe_mesh
end

function default_groups_cartesian_grid(fe_mesh)
  D = domain_dim(fe_mesh)
  cell_to_nodes = face_nodes(fe_mesh,D)
  refcell = first(ref_faces(fe_mesh,D))
  nnodes = num_nodes(fe_mesh)
  d_to_ldface_to_lnodes = [ face_nodes(refcell,d) for d in 0:(D-1)]

  node_to_n = zeros(Int32,nnodes)
  for nodes in cell_to_nodes
    for node in nodes
      node_to_n[node] += Int32(1)
    end
  end
  J = typeof(JaggedArray(Vector{Int32}[]))
  face_to_nodes = Vector{J}(undef,D)
  groups = GroupCollection(VOID,D)
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
    dface_to_group = zeros(Int32,ndfaces)
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
          dface_to_group[ndfaces] = group
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
      add_group!(groups,d,"$(d)-face-$ldface",group)
      faces_in_group = findall(g->g==group,dface_to_group)
      group_faces!(groups,faces_in_group,group)
    end
    ngroups += nldfaces
  end # d
  ngroups += 1
  add_group!(groups,D,"$(D)-face-1",ngroups)
  ncells = length(cell_to_nodes)
  group_faces!(groups,collect(Int32,1:ncells),ngroups)


  face_to_refid = [ ones(Int8,length(face_to_nodes[d+1]))  for d in 0:(D-1)]
  refid_to_faces = [ ref_faces(refcell,d)  for d in 0:(D-1)]
  groups, (face_to_nodes,face_to_refid,refid_to_faces)
end






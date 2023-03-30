#TODO
# Physical groups as collection of Pairs
# Visualize all dims in vtk_args if dim not provided
# Remove Generic e.g GenericMesh -> Mesh and in the future: const Mesh = Union{MeshGeneric,MeshNative}
# Groups as vector of pairs or dict of pairs? First more lightweight second more general
# Goups as Dict{String,Vector{Int32}} ?
# Groupname as Symbol or String?
# Replace name by tag? in groups?
# SimpleMesh
# remove new_mesh and use GenericMesh
# rename gmsh_mesh to mesh_from_gmsh
# Val(d) also in other functions? YES
# better way to represent hanging node constraints? Solution: global constraints
# follow the same approach for periodic
# groups for faces and nodes: no. Node groups only make sense for for Lagrangian spaces of the same order of the mesh
# rename physical groups by  physical_faces: NO.
# num_faces(a) by number_of_faces(a) ?
# allow to use custom integer and Float precision in mesh_from_forest
# HangingNodeConstraints to GenericHangingNodeconstraints
# Create new structs to avoid burden of type parameters
# node coordinates in p4est: use the original node numeration and the transpose of the leaf constraints. then apply the global constraints to get the final coordinates
# Index groups by symbol and remove GenericPhysicalGroup
# float_type
# integer_type
# no hard-coded int nor float types
# TODO implement efficient mul!

# Nobody should overwrite these
has_periodic_nodes(a) = has_periodic_nodes(typeof(a))
has_hanging_nodes(a) = has_hanging_nodes(typeof(a))
has_physical_groups(a) = has_physical_groups(typeof(a))

# Defaults
has_periodic_nodes(::Type) = false
has_hanging_nodes(::Type) = false
has_physical_groups(::Type) = false
is_simplex(a) = false
is_hypercube(a) = false
num_faces(a,d) = length(face_reference_id(a,d))
num_faces(a) = map(d->num_faces(a,d),0:dimension(a))
function face_offsets(a)
    D = dimension(a)
    offsets = zeros(Int,D+1)
    for d in 1:D
        offsets[d+1] = offsets[d] + num_faces(a,d-1)
    end
    offsets
end
num_nodes(a) = length(node_coordinates(a))
embedded_dimension(a) = length(eltype(node_coordinates(a)))
reference_faces(a,d) = reference_faces(a,Val(d))
reference_faces(a,::Val{d}) where d = error("reference_faces($(typeof(a)),::Val{$d}) not implemented.")

abstract type AbstractMeshWithData{A} end
has_periodic_nodes(::Type{<:AbstractMeshWithData{A}}) where A = has_periodic_nodes(A)
has_hanging_nodes(::Type{<:AbstractMeshWithData{A}}) where A = has_hanging_nodes(A)
has_physical_groups(::Type{<:AbstractMeshWithData{A}}) where A = has_physical_groups(A)
is_simplex(a::AbstractMeshWithData) = is_simplex(a.mesh)
is_hypercube(a::AbstractMeshWithData) = is_hypercube(a.mesh)
num_faces(a::AbstractMeshWithData,d) = num_faces(a.mesh,d)
dimension(a::AbstractMeshWithData) = dimension(a.mesh)
embedded_dimension(a::AbstractMeshWithData) = embedded_dimension(a.mesh)
node_coordinates(a::AbstractMeshWithData) = node_coordinates(a.mesh)
face_nodes(a::AbstractMeshWithData,d) = face_nodes(a.mesh,d)
face_reference_id(a::AbstractMeshWithData,d) = face_reference_id(a.mesh,d)
reference_faces(a::AbstractMeshWithData,::Val{d}) where d = reference_faces(a.mesh,Val{d}())
num_faces(a::AbstractMeshWithData) = num_faces(a.mesh)
periodic_node_constraints(a::AbstractMeshWithData) = periodic_node_constraints(a.mesh)
hanging_node_constraints(a::AbstractMeshWithData) = hanging_node_constraints(a.mesh)
physical_groups(a::AbstractMeshWithData,d) = physical_groups(a.mesh,d)

struct GenericMesh{A,B,C,D}
    node_coordinates::A
    face_nodes::Vector{B}
    face_reference_id::Vector{C}
    reference_faces::D
end
dimension(a::GenericMesh) = length(a.reference_faces)-1
node_coordinates(a::GenericMesh) = a.node_coordinates
face_nodes(a::GenericMesh,d) = a.face_nodes[d+1]
face_reference_id(a::GenericMesh,d) = a.face_reference_id[d+1]
reference_faces(a::GenericMesh,::Val{d}) where d = a.reference_faces[d+1]

set_phyisical_groups(mesh,groups) = MeshWithPhysicalGroups(mesh,groups)
struct MeshWithPhysicalGroups{A,B} <: AbstractMeshWithData{A}
    mesh::A
    physical_groups::B
end
has_physical_groups(::Type{<:MeshWithPhysicalGroups}) = true
physical_groups(a::MeshWithPhysicalGroups) = a.physical_groups
physical_groups(a::MeshWithPhysicalGroups,d) = a.physical_groups[d+1]

function add_trivial_physical_groups(a)
    D = dimension(a)
    groups = map(0:D) do d
        ["$(d)face$(f)"=>[Int32(f)] for f in 1:num_faces(a,d)]
    end
    set_phyisical_groups(a,groups)
end

partition_from_mask(a) = partition_from_mask(identity,a)

function partition_from_mask(f,node_to_mask)
    T = Vector{Int32}
    free_nodes = convert(T,findall(f,node_to_mask))
    dirichlet_nodes = convert(T,findall(i->!f(i),node_to_mask))
    nfree = length(free_nodes)
    ndiri = length(dirichlet_nodes)
    permutation = T(undef,nfree+ndiri)
    permutation[free_nodes] = 1:nfree
    permutation[dirichlet_nodes] = (1:ndiri) .+ nfree
    TwoPartPartition(free_nodes,dirichlet_nodes,permutation)
end

struct TwoPartPartition{A} <: AbstractVector{A}
    first::A
    last::A
    permutation::A
end

permutation(a::TwoPartPartition) = a.permutation
Base.size(a::TwoPartPartition) = (2,)
Base.IndexStyle(::Type{<:TwoPartPartition}) = IndexLinear()
function Base.getindex(a::TwoPartPartition,i::Int)
    @boundscheck @assert i in (1,2)
    if i == 1
        a.first
    else
        a.last
    end
end
struct GenericPeriodicNodeConstraints{A,B,C,T} <: AbstractMatrix{T}
    nrows::Int
    ncols::Int
    master_nodes::A
    master_coeffs::B
    free_and_periodic_nodes::TwoPartPartition{C}
    function GenericPeriodicNodeConstraints(
            nrows,
            ncols,
            master_nodes,
            master_coeffs,
            free_and_periodic_nodes::TwoPartPartition{C}) where C
        A = typeof(master_nodes)
        B = typeof(master_coeffs)
        T = eltype(master_coeffs)
        @assert nrows == length(permutation(free_and_periodic_nodes))
        new{A,B,C,T}(nrows,ncols,master_nodes,master_coeffs,free_and_periodic_nodes)
    end
end

has_periodic_nodes(::Type{<:GenericPeriodicNodeConstraints}) = true
free_and_periodic_nodes(a::GenericPeriodicNodeConstraints) = a.free_and_periodic_nodes
periodic_nodes(a::GenericPeriodicNodeConstraints) = last(free_and_periodic_nodes(a))
free_nodes(a::GenericPeriodicNodeConstraints) = first(free_and_periodic_nodes(a))
node_permutation(a::GenericPeriodicNodeConstraints) = permutation(free_and_periodic_nodes(a))
master_nodes(a::GenericPeriodicNodeConstraints) = a.master_nodes
master_coeffs(a::GenericPeriodicNodeConstraints) = a.master_coeffs

Base.size(a::GenericPeriodicNodeConstraints) = (a.nrows,a.ncols)
Base.IndexStyle(::Type{<:GenericPeriodicNodeConstraints}) = IndexCartesian()

function Base.getindex(a::GenericPeriodicNodeConstraints,i::Int,j::Int)
    T = eltype(a)
    m,n = size(a)
    @boundscheck @assert i in 1:m
    @boundscheck @assert j in 1:n

    free_nodes, master_nodes = free_and_periodic_nodes(a)
    node_permutation = permutation(free_and_periodic_nodes(a))
    p = node_permutation[i]
    nfree = length(free_nodes)
    if p > nfree
        h = p-nfree
        master_nodes[h] == j ? a.master_coeffs[h] : zero(T)
    else
        free_nodes[p] == j ? one(T) : zero(T)
    end
end

#TODO implement mul!

set_periodic_node_constraints(mesh,constraints) = MeshWithPeriodicNodeConstraints(mesh,constraints)
struct MeshWithPeriodicNodeConstraints{A,B} <: AbstractMeshWithData{A}
    mesh::A
    constraints::B
end
has_periodic_nodes(::Type{<:MeshWithPeriodicNodeConstraints{A,B}}) where {A,B} = has_periodic_nodes(B)
periodic_node_constraints(a::MeshWithPeriodicNodeConstraints) = a.constraints

set_haning_node_constraints(mesh,constraints) = MeshWithHangingNodeConstraints(mesh,constraints)
struct MeshWithHangingNodeConstraints{A,B} <: AbstractMeshWithData{A}
    mesh::A
    constraints::B
end
has_hanging_nodes(::Type{<:MeshWithHangingNodeConstraints{A,B}}) where {A,B} = has_hanging_nodes(B)
hanging_node_constraints(a::MeshWithHangingNodeConstraints) = a.constraints

struct HangingNodeConstraints{A,B,C,T} <: AbstractMatrix{T}
    nrows::Int
    ncols::Int
    master_nodes::A
    master_coeffs::B
    free_and_hanging_nodes::TwoPartPartition{C}
    function HangingNodeConstraints(
            nrows,
            ncols,
            master_nodes,
            master_coeffs,
            free_and_hanging_nodes::TwoPartPartition{C}) where C
        A = typeof(master_nodes)
        B = typeof(master_coeffs)
        T = eltype(eltype(master_coeffs))
        @assert nrows == length(permutation(free_and_hanging_nodes))
        new{A,B,C,T}(nrows,ncols,master_nodes,master_coeffs,free_and_hanging_nodes)
    end
end

free_and_hanging_nodes(a::HangingNodeConstraints) = a.free_and_hanging_nodes
hanging_nodes(a::HangingNodeConstraints) = last(free_and_hanging_nodes(a))
free_nodes(a::HangingNodeConstraints) = first(free_and_periodic_nodes(a))
node_permutation(a::HangingNodeConstraints) = permutation(free_and_hanging_nodes(a))
master_nodes(a::HangingNodeConstraints) = a.master_nodes
master_coeffs(a::HangingNodeConstraints) = a.master_coeffs

Base.size(a::HangingNodeConstraints) = (a.nrows,a.ncols)
Base.IndexStyle(::Type{<:HangingNodeConstraints}) = IndexCartesian()

function Base.getindex(a::HangingNodeConstraints,i::Int,j::Int)
    T = eltype(a)
    m,n = size(a)
    @boundscheck @assert i in 1:m
    @boundscheck @assert j in 1:n
    free_nodes, hanging_nodes = free_and_hanging_nodes(a)
    node_permutation = permutation(free_and_hanging_nodes(a))
    p = node_permutation[i]
    nfree = length(free_nodes)
    if p > nfree
        h = p-nfree
        masters = a.master_nodes[h]
        l = findfirst(k->k==j,masters)
        if l === nothing
            zero(T)
        else
            coeffs = a.master_coeffs[h]
            coeffs[l]
        end
    else
        free_nodes[p]==j ? one(T) : zero(T)
    end
end

## TODO implement mul and also for its transpose

# Node
# Face aka mesh from Boundary
# Mesh

# NodeTopology
# FaceTopology
# MeshTopology

struct NodeTopology end
dimension(a::NodeTopology) = 0
function face_incidence(a::NodeTopology,d1,d2)
    @boundscheck @assert d1 == d2 == 0
    [[Int32(1)]]
end
function face_reference_id(a::NodeTopology,d1)
    @boundscheck @assert d1 == 0
    [Int32(1)]
end
function reference_faces(a::NodeTopology,::Val{0})
    (a,)
end

struct FaceTopology{A}
    boundary::A
end
dimension(a::FaceTopology) = dimension(a.boundary) + 1
function face_incidence(a::FaceTopology,d1,d2)
    d = dimension(a)
    @boundscheck @assert d1 <= d && d2 <= d
    if d1 == d && d2 == d
        [[Int32(1)]] 
    elseif d1 == d && d2 < d
        [collect(Int32,1:num_faces(a.boundary,d2))]
    elseif d1 < d && d2 == d
        [ [Int32(i)] for i in 1:num_faces(a.boundary,d1)]
    else
        face_incidence(a.boundary,d1,d2)
    end
end
function face_reference_id(a::FaceTopology,d1)
    d = dimension(a)
    @boundscheck @assert d1 <= d 
    if d1 == d
        [Int32(1)]
    else
        a.face_reference_id[d1+1]
    end
end
function reference_faces(a::FaceTopology,::Val{d1}) where d1
    d = dimension(a)
    @assert  d >= d1
    if d1 == d
        (a,)
    else
        reference_faces(a.boundary,Val(d1))
    end
end

struct MeshTopology{A,B,C}
    face_incidence::A
    face_reference_id::B
    reference_faces::C
end
dimension(a::MeshTopology) = length(a.face_reference_id)
face_incidence(a::MeshTopology,d1,d2) = a.face_incidence[d1+1,d2+1]
face_reference_id(a::MeshTopology,d1) = a.face_reference_id[d1+1]
function reference_faces(a::MeshTopology,::Val{d}) where d
    @assert dimension(a) >= d
    a.reference_faces[d+1]
end

function face_boundary(a)
    D = dimension(a)
    my_coords = node_coordinates(a)
    my_face_nodes = [face_nodes(a,d) for d in 0:(D-1)]
    my_face_reference_id = [face_reference_id(a,d) for d in 0:(D-1)]
    my_reference_faces = Tuple([reference_faces(a,d) for d in 0:(D-1)])
    GenericMesh(my_coords,my_face_nodes,my_face_reference_id,my_reference_faces)
end

function face_topology(a)
    if dimension(a) == 0
        NodeTopology()
    else
        mesh = face_boundary(a)
        topo = mesh_topology(mesh)
        FaceTopology(topo)
    end
end

function mesh_topology(a)
    # Assumes that the input is a cell complex
    T = JaggedArray{Int32,Int32}
    D = dimension(a)
    my_face_incidence = Matrix{T}(undef,D+1,D+1)
    my_face_reference_id  = [ face_reference_id(a,d) for d in 0:D ]
    my_reference_faces = Tuple([ map(face_topology,reference_faces(a,d)) for d in 0:D ])
    topo = MeshTopology(my_face_incidence,my_face_reference_id,my_reference_faces)
    for d in 0:D
        fill_face_interior!(topo,a,d)
    end
    for d in 1:D
        fill_face_vertices!(topo,a,d)
        fill_face_coboundary!(topo,a,d,0)
    end
    for d in 1:(D-1)
        for n in (D-d):-1:1
            m = n+d
            @show (m,n)
            fill_face_boundary!(topo,a,m,n)
            display(face_incidence(topo,m,n))
            fill_face_coboundary!(topo,a,m,n)
        end
    end
    topo
end

function fill_face_interior!(topo::MeshTopology,mesh,d)
    n = num_faces(mesh,d)
    ptrs = collect(Int32,1:(n+1))
    data = collect(Int32,1:n)
    topo.face_incidence[d+1,d+1] = JaggedArray(data,ptrs)
end

function fill_face_coboundary!(topo::MeshTopology,mesh,n,m)
    function barrier(nface_to_mfaces,nmfaces)
        ptrs = zeros(Int32,nmfaces+1)
        nnfaces = length(nface_to_mfaces)
        for nface in 1:nnfaces
            mfaces = nface_to_mfaces[nface]
            for mface in mfaces
                ptrs[mface+1] += Int32(1)
            end
        end
        length_to_ptrs!(ptrs)
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
        rewind_ptrs!(ptrs)
        mface_to_nfaces = GenericJaggedArray(data,ptrs)
        mface_to_nfaces
    end
    nmfaces = num_faces(mesh,m)
    nface_to_mfaces = face_incidence(topo,n,m)
    topo.face_incidence[m+1,n+1] = barrier(nface_to_mfaces,nmfaces)
end

function fill_face_vertices!(topo::MeshTopology,mesh,d)
    function barrier(nnodes,vertex_to_nodes,dface_to_nodes,dface_to_refid,refid_to_lvertex_to_lnodes)
        vertex_to_node = JaggedArray(vertex_to_nodes).data
        #if vertex_to_node == 1:nnodes
        #    return dface_to_nodes
        #end
        node_to_vertex = zeros(Int32,nnodes)
        nvertices = length(vertex_to_nodes)
        node_to_vertex[vertex_to_node] = 1:nvertices
        ndfaces = length(dface_to_nodes)
        dface_to_vertices_ptrs = zeros(Int32,ndfaces+1)
        for dface in 1:ndfaces
            refid = dface_to_refid[dface]
            nlvertices = length(refid_to_lvertex_to_lnodes[refid])
            dface_to_vertices_ptrs[dface+1] = nlvertices
        end
        length_to_ptrs!(dface_to_vertices_ptrs)
        ndata = dface_to_vertices_ptrs[end]-1
        dface_to_vertices_data = zeros(Int32,ndata)
        for dface in 1:ndfaces
            refid = dface_to_refid[dface]
            lvertex_to_lnodes = refid_to_lvertex_to_lnodes[refid]
            nlvertices = length(lvertex_to_lnodes)
            lnode_to_node = dface_to_nodes[dface]
            offset = dface_to_vertices_ptrs[dface]-1
            for lvertex in 1:nlvertices
                lnode = first(lvertex_to_lnodes[lvertex])
                node = lnode_to_node[lnode]
                vertex = node_to_vertex[node]
                dface_to_vertices_data[offset+lvertex] = vertex
            end
        end
        dface_to_vertices = JaggedArray(dface_to_vertices_data,dface_to_vertices_ptrs)
    end
    nnodes = num_nodes(mesh)
    vertex_to_nodes = face_nodes(mesh,0)
    dface_to_nodes = face_nodes(mesh,d)
    dface_to_refid = face_reference_id(mesh,d)
    refid_refface = reference_faces(mesh,d)
    refid_to_lvertex_to_lnodes = map(refface->face_nodes(refface,0),refid_refface)
    face_to_vertices = barrier(nnodes,vertex_to_nodes,dface_to_nodes,dface_to_refid,refid_to_lvertex_to_lnodes)
    topo.face_incidence[d+1,0+1] = face_to_vertices
end

function fill_face_boundary!(topo::MeshTopology,mesh,D,d)
    function barrier(
            Dface_to_vertices,
            vertex_to_Dfaces,
            dface_to_vertices,
            vertex_to_dfaces,
            Dface_to_refid,
            Drefid_to_ldface_to_lvertices)

        # Count
        ndfaces = length(dface_to_vertices)
        nDfaces = length(Dface_to_vertices)
        # Allocate output
        ptrs = zeros(Int32,nDfaces+1)
        for Dface in 1:nDfaces
            Drefid = Dface_to_refid[Dface]
            ldface_to_lvertices = Drefid_to_ldface_to_lvertices[Drefid]
            ptrs[Dface+1] = length(ldface_to_lvertices)
        end
        length_to_ptrs!(ptrs)
        ndata = ptrs[end]-1
        data = fill(Int32(INVALID_ID),ndata)
        Dface_to_dfaces = JaggedArray(data,ptrs)
        for Dface in 1:nDfaces
            Drefid = Dface_to_refid[Dface]
            ldface_to_lvertices = Drefid_to_ldface_to_lvertices[Drefid]
            lvertex_to_vertex = Dface_to_vertices[Dface]
            ldface_to_dface = Dface_to_dfaces[Dface]
            for (ldface,lvertices) in enumerate(ldface_to_lvertices)
                # Find the global d-face for this local d-face
                dface2 = Int32(INVALID_ID)
                vertices = view(lvertex_to_vertex,lvertices)
                for (i,lvertex) in enumerate(lvertices)
                    vertex = lvertex_to_vertex[lvertex]
                    dfaces = vertex_to_dfaces[vertex]
                    for dface1 in dfaces
                        vertices1 = dface_to_vertices[dface1]
                        if same_valid_ids(vertices,vertices1)
                            dface2 = dface1
                            break
                        end
                    end
                    if dface2 != Int32(INVALID_ID)
                        break
                    end
                end
                @boundscheck @assert dface2 != Int32(INVALID_ID)
                ldface_to_dface[ldface] = dface2
            end # (ldface,lvertices)
        end # Dface
        Dface_to_dfaces
    end
    Dface_to_vertices = face_incidence(topo,D,0)
    vertex_to_Dfaces = face_incidence(topo,0,D)
    dface_to_vertices = face_incidence(topo,d,0)
    vertex_to_dfaces = face_incidence(topo,0,d)
    Dface_to_refid = face_reference_id(topo,D)
    refid_refface = reference_faces(topo,D)
    Drefid_to_ldface_to_lvertices = map(refface->face_incidence(refface,d,0),refid_refface)
    Dface_to_dfaces = barrier(
            Dface_to_vertices,
            vertex_to_Dfaces,
            dface_to_vertices,
            vertex_to_dfaces,
            Dface_to_refid,
            Drefid_to_ldface_to_lvertices)
    topo.face_incidence[D+1,d+1] = Dface_to_dfaces
end

const INVALID_ID = 0

function intersection!(a,b,na,nb)
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

function same_valid_ids(a,b)
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

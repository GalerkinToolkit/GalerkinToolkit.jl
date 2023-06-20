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
            fill_face_boundary!(topo,a,m,n)
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

function generate_face_coboundary(nface_to_mfaces,nmfaces)
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
    mface_to_nfaces = JaggedArray(data,ptrs)
    mface_to_nfaces
end

function fill_face_coboundary!(topo::MeshTopology,mesh,n,m)
    nmfaces = num_faces(mesh,m)
    nface_to_mfaces = face_incidence(topo,n,m)
    topo.face_incidence[m+1,n+1] = generate_face_coboundary(nface_to_mfaces,nmfaces)
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
    topo.face_incidence[d+1,0+1] = barrier(nnodes,vertex_to_nodes,dface_to_nodes,dface_to_refid,refid_to_lvertex_to_lnodes)
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

function face_complex(mesh)
    Ti = Int32
    T = JaggedArray{Ti,Ti}
    D = dimension(mesh)
    oldface_to_newvertices = Vector{T}(undef,D+1)
    newvertex_to_oldfaces = Vector{T}(undef,D+1)
    newface_incidence = Matrix{T}(undef,D+1,D+1)
    nnewfaces = zeros(Int,D+1)
    newface_refid = Vector{Vector{Ti}}(undef,D+1)
    newreffaces = Vector{Any}(undef,D+1)
    newface_nodes = Vector{T}(undef,D+1)
    old_to_new = Vector{Vector{Ti}}(undef,D+1)
    node_to_newvertex, n_new_vertices = find_node_to_vertex(mesh) # Optimizable for linear meshes
    for d in 0:D
        oldface_to_newvertices[d+1] = fill_face_vertices(mesh,d,node_to_newvertex) # Optimizable for linear meshes
        newvertex_to_oldfaces[d+1] = generate_face_coboundary(oldface_to_newvertices[d+1],n_new_vertices) # Optimizable for linear meshes
    end
    newface_incidence[D+1,0+1] = oldface_to_newvertices[D+1]
    newface_incidence[0+1,D+1] = newvertex_to_oldfaces[D+1]
    nnewfaces[D+1] = length(oldface_to_newvertices[D+1])
    newface_refid[D+1] = face_reference_id(mesh,D)
    newreffaces[D+1] = reference_faces(mesh,D)
    newface_nodes[D+1] = face_nodes(mesh,D)
    old_to_new[D+1] = collect(Ti,1:length(newface_nodes[D+1]))
    # TODO optimize for d==0
    for d in (D-1):-1:0
        n = d+1
        new_nface_to_new_vertices = newface_incidence[n+1,0+1]
        new_vertex_to_new_nfaces = newface_incidence[0+1,n+1]
        old_dface_to_new_vertices = oldface_to_newvertices[d+1]
        new_vertex_to_old_dfaces = newvertex_to_oldfaces[d+1]
        new_nface_to_nrefid = newface_refid[n+1]
        old_dface_to_drefid = face_reference_id(mesh,d)
        drefid_to_ref_dface = reference_faces(mesh,d)
        old_dface_to_nodes = face_nodes(mesh,d)
        new_nface_to_nodes = newface_nodes[n+1]
        nrefid_to_ldface_to_lvertices = map(a->face_incidence(face_topology(a),d,0),newreffaces[n+1])
        nrefid_to_ldface_to_lnodes = map(a->face_nodes(a,d),newreffaces[n+1])
        nrefid_to_ldface_to_drefrefid = map(a->face_reference_id(a,d),newreffaces[n+1])
        nrefid_to_drefrefid_to_ref_dface = map(a->reference_faces(a,d),newreffaces[n+1])
        new_nface_to_new_dfaces, n_new_dfaces, old_dface_to_new_dface = generate_face_boundary(
            new_nface_to_new_vertices,
            new_vertex_to_new_nfaces,
            old_dface_to_new_vertices,
            new_vertex_to_old_dfaces,
            new_nface_to_nrefid,
            nrefid_to_ldface_to_lvertices)
        new_dface_to_new_vertices = generate_face_vertices(
            n_new_dfaces,
            old_dface_to_new_dface,
            old_dface_to_new_vertices,
            new_nface_to_new_vertices,
            new_nface_to_new_dfaces,
            new_nface_to_nrefid,
            nrefid_to_ldface_to_lvertices)
        new_vertex_to_new_dfaces = generate_face_coboundary(new_dface_to_new_vertices,n_new_vertices)
        new_dface_to_nodes = generate_face_vertices(
            n_new_dfaces,
            old_dface_to_new_dface,
            old_dface_to_nodes,
            new_nface_to_nodes,
            new_nface_to_new_dfaces,
            new_nface_to_nrefid,
            nrefid_to_ldface_to_lnodes)
        new_dface_to_new_drefid, new_refid_to_ref_dface = generate_reference_faces(
            n_new_dfaces,
            old_dface_to_new_dface,
            old_dface_to_drefid,
            drefid_to_ref_dface,
            new_nface_to_new_dfaces,
            new_nface_to_nrefid,
            nrefid_to_ldface_to_drefrefid,
            nrefid_to_drefrefid_to_ref_dface)
        newface_incidence[n+1,d+1] = new_nface_to_new_dfaces
        newface_incidence[d+1,0+1] = new_dface_to_new_vertices
        newface_incidence[0+1,d+1] = new_vertex_to_new_dfaces
        newface_refid[d+1] = new_dface_to_new_drefid
        newreffaces[d+1] = new_refid_to_ref_dface
        nnewfaces[d+1] = n_new_dfaces
        newface_nodes[d+1] = new_dface_to_nodes
        old_to_new[d+1] = old_dface_to_new_dface
    end
    node_to_coords = node_coordinates(mesh)
    mesh = GenericMesh(node_to_coords,newface_nodes,newface_refid,Tuple(newreffaces))
    mesh, old_to_new
end

function generate_face_vertices(
    n_new_dfaces,
    old_dface_to_new_dface,
    old_dface_to_new_vertices,
    new_nface_to_new_vertices,
    new_nface_to_new_dfaces,
    new_nface_to_nrefid,
    nrefid_to_ldface_to_lvertices
    )

    Ti = eltype(eltype(old_dface_to_new_vertices))
    new_dface_to_touched = fill(false,n_new_dfaces)
    new_dface_to_new_vertices_ptrs = zeros(Ti,n_new_dfaces+1)
    n_old_dfaces = length(old_dface_to_new_dface)
    for old_dface in 1:n_old_dfaces
        new_dface = old_dface_to_new_dface[old_dface]
        new_vertices = old_dface_to_new_vertices[old_dface]
        new_dface_to_new_vertices_ptrs[new_dface+1] = length(new_vertices)
        new_dface_to_touched[new_dface] = true
    end
    n_new_nfaces = length(new_nface_to_new_vertices)
    for new_nface in 1:n_new_nfaces
        nrefid = new_nface_to_nrefid[new_nface]
        ldface_to_lvertices = nrefid_to_ldface_to_lvertices[nrefid]
        ldface_to_new_dface = new_nface_to_new_dfaces[new_nface]
        n_ldfaces = length(ldface_to_new_dface)
        for ldface in 1:n_ldfaces
            new_dface = ldface_to_new_dface[ldface]
            if new_dface_to_touched[new_dface]
                continue
            end
            lvertices = ldface_to_lvertices[ldface]
            new_dface_to_new_vertices_ptrs[new_dface+1] = length(lvertices)
            new_dface_to_touched[new_dface] = true
        end

    end
    length_to_ptrs!(new_dface_to_new_vertices_ptrs)
    ndata = new_dface_to_new_vertices_ptrs[end]-1
    new_dface_to_new_vertices_data = zeros(Ti,ndata)
    new_dface_to_new_vertices = JaggedArray(new_dface_to_new_vertices_data,new_dface_to_new_vertices_ptrs)
    fill!(new_dface_to_touched,false)
    for old_dface in 1:n_old_dfaces
        new_dface = old_dface_to_new_dface[old_dface]
        new_vertices_in = old_dface_to_new_vertices[old_dface]
        new_vertices_out = new_dface_to_new_vertices[new_dface]
        for i in 1:length(new_vertices_in)
            new_vertices_out[i] = new_vertices_in[i]
        end
        new_dface_to_touched[new_dface] = true
    end
    for new_nface in 1:n_new_nfaces
        nrefid = new_nface_to_nrefid[new_nface]
        ldface_to_lvertices = nrefid_to_ldface_to_lvertices[nrefid]
        ldface_to_new_dface = new_nface_to_new_dfaces[new_nface]
        n_ldfaces = length(ldface_to_new_dface)
        new_vertices_in = new_nface_to_new_vertices[new_nface]
        for ldface in 1:n_ldfaces
            new_dface = ldface_to_new_dface[ldface]
            if new_dface_to_touched[new_dface]
                continue
            end
            new_vertices_out = new_dface_to_new_vertices[new_dface]
            lvertices = ldface_to_lvertices[ldface]
            for i in 1:length(lvertices)
                new_vertices_out[i] = new_vertices_in[lvertices[i]]
            end
            new_dface_to_touched[new_dface] = true
        end

    end
    new_dface_to_new_vertices
end

function generate_reference_faces(
        n_new_dfaces,
        old_dface_to_new_dface,
        old_dface_to_drefid,
        drefid_to_ref_dface,
        new_nface_to_new_dfaces,
        new_nface_to_nrefid,
        nrefid_to_ldface_to_drefrefid,
        nrefid_to_drefrefid_to_ref_dface)

    i_to_ref_dface = collect(Any,drefid_to_ref_dface)
    drefid_to_i = collect(1:length(drefid_to_ref_dface))
    i = length(i_to_ref_dface)
    Ti = Int32
    nrefid_to_drefrefid_to_i = map(a->zeros(Ti,length(a)),nrefid_to_drefrefid_to_ref_dface)
    for (nrefid,drefrefid_to_ref_dface) in enumerate(nrefid_to_drefrefid_to_ref_dface)
        for (drefrefid, ref_dface) in enumerate(drefrefid_to_ref_dface)
            push!(i_to_ref_dface,ref_dface)
            i += 1
            nrefid_to_drefrefid_to_i[nrefid][drefrefid] = i
        end
    end
    u_to_ref_dface = unique(i_to_ref_dface)
    i_to_u = indexin(i_to_ref_dface,u_to_ref_dface)
    new_dface_to_u = zeros(Ti,n_new_dfaces)
    new_dface_to_touched = fill(false,n_new_dfaces)
    for (old_dface,new_dface) in enumerate(old_dface_to_new_dface)
        drefid = old_dface_to_drefid[old_dface]
        i = drefid_to_i[drefid]
        u = i_to_u[i]
        new_dface_to_u[new_dface] = u
        new_dface_to_touched[new_dface] = true
    end
    for (new_nface,nrefid) in enumerate(new_nface_to_nrefid)
        ldface_to_drefrefid = nrefid_to_ldface_to_drefrefid[nrefid]
        ldface_to_new_dface = new_nface_to_new_dfaces[new_nface]
        drefrefid_to_i = nrefid_to_drefrefid_to_i[nrefid]
        for (ldface,new_dface) in enumerate(ldface_to_new_dface)
            new_dface = ldface_to_new_dface[ldface]
            if new_dface_to_touched[new_dface]
                continue
            end
            drefrefid = ldface_to_drefrefid[ldface]
            i = drefrefid_to_i[drefrefid]
            u = i_to_u[i]
            new_dface_to_u[new_dface] = u
            new_dface_to_touched[new_dface] = true
        end
    end
    new_dface_to_u, Tuple(u_to_ref_dface)
end

function generate_face_boundary(
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
  length_to_ptrs!(ptrs)
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
  old_to_new = collect(Int32,1:ndfaces)
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
          if same_valid_ids(vertices,vertices1)
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
          intersection!(Dfaces1,Dfaces2,nDfaces1,nDfaces2)
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
            if same_valid_ids(vertices,vertices1)
              ldface2 = ldface1
              break
            end
          end
          @boundscheck @assert ldface2 != INVALID_ID
          Dface_to_dfaces[Dface1][ldface2] = dface2
        end
      end
    end # (ldface,lvertices)
  end # Dface
  Dface_to_dfaces, newdface, old_to_new
end

function fill_face_vertices(mesh,d,node_to_vertex)
    function barrier(node_to_vertex,face_to_nodes,face_to_refid,refid_to_lvertex_to_lnodes)
        Ti = eltype(node_to_vertex)
        nfaces = length(face_to_nodes)
        face_to_vertices_ptrs = zeros(Ti,nfaces+1)
        for face in 1:nfaces
            nodes = face_to_nodes[face]
            refid = face_to_refid[face]
            lvertex_to_lnodes = refid_to_lvertex_to_lnodes[refid]
            nlvertices = length(lvertex_to_lnodes)
            face_to_vertices_ptrs[face+1] = nlvertices
        end
        length_to_ptrs!(face_to_vertices_ptrs)
        ndata = face_to_vertices_ptrs[end]-1
        face_to_vertices_data = zeros(Ti,ndata)
        face_to_vertices = JaggedArray(face_to_vertices_data,face_to_vertices_ptrs)
        for face in 1:nfaces
            vertices = face_to_vertices[face]
            nodes = face_to_nodes[face]
            refid = face_to_refid[face]
            lvertex_to_lnodes = refid_to_lvertex_to_lnodes[refid]
            for (lvertex,lnodes) in enumerate(lvertex_to_lnodes)
                lnode = first(lnodes)
                vertex = node_to_vertex[nodes[lnode]]
                @boundscheck @assert vertex != INVALID_ID
                vertices[lvertex] = vertex
            end
        end
        face_to_vertices
    end
    face_to_nodes = face_nodes(mesh,d)
    face_to_refid = face_reference_id(mesh,d)
    refid_to_lvertex_to_lnodes = map(a->face_nodes(a,0),reference_faces(mesh,d))
    barrier(node_to_vertex,face_to_nodes,face_to_refid,refid_to_lvertex_to_lnodes)
end

function find_node_to_vertex(mesh)
    function barrier!(node_to_vertex,face_to_nodes,face_to_refid,refid_to_lvertex_to_lnodes)
        valid_id = one(eltype(node_to_vertex))
        for face in eachindex(face_to_nodes)
            nodes = face_to_nodes[face]
            refid = face_to_refid[face]
            lvertex_to_lnodes = refid_to_lvertex_to_lnodes[refid]
            for lnodes in lvertex_to_lnodes
                lnode = first(lnodes)
                node = nodes[lnode]
                node_to_vertex[node] = valid_id
            end
        end
    end
    Ti = Int32
    nnodes = num_nodes(mesh)
    node_to_vertex = zeros(Ti,nnodes)
    fill!(node_to_vertex,Ti(INVALID_ID))
    D = dimension(mesh)
    for d in 0:D
        face_to_nodes = face_nodes(mesh,d)
        face_to_refid = face_reference_id(mesh,d)
        refid_to_lvertex_to_lnodes = map(a->face_nodes(a,0),reference_faces(mesh,d))
        barrier!(node_to_vertex,face_to_nodes,face_to_refid,refid_to_lvertex_to_lnodes)
    end
    vertex = Ti(0)
    for node in eachindex(node_to_vertex)
        if node_to_vertex[node] != Ti(INVALID_ID)
            vertex += Ti(1)
            node_to_vertex[node] = vertex
        end
    end
    node_to_vertex, vertex
end








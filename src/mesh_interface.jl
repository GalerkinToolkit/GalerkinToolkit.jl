#TODO
# Physical groups as collection of Pairs
# Visualize all dims in vtk_args if dim not provided
# Remove Generic e.g GenericMesh -> Mesh and in the future: const Mesh = Union{MeshGeneric,MeshNative}
# Groups as vector of pairs or dict of pairs? First more lightweight second more general
# Goups as Dict{String,Vector{Int32}} ?
# Groupname as Symbol or String?
# SimpleMesh
# remove new_mesh and use GenericMesh
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
physical_groups(a::MeshWithPhysicalGroups,d) = a.physical_groups[d+1]

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


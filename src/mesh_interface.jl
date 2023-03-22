#TODO
# SimpleMesh
# remove new_mesh and use GenericMesh
# Val(d) also in other functions? YES
# better way to represent hanging node constraints? Solution: global constraints
# follow the same approach for periodic
# groups for faces and nodes
# rename physical groups by groups
# num_faces(a) by number_of_faces(a) ?
# allow to use custom integer and Float precision in mesh_from_forest
# HangingNodeConstraints to GenericHangingNodeconstraints
# Create new structs to avoid burden of type parameters

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
periodic_nodes(a::AbstractMeshWithData) = periodic_nodes(a.mesh)
periodic_to_master(a::AbstractMeshWithData) = periodic_nodes(a.mesh)
periodic_to_coeff(a::AbstractMeshWithData) = periodic_to_coeff(a.mesh)
hanging_nodes(a::AbstractMeshWithData) = hanging_nodes(a.mesh)
hanging_to_masters(a::AbstractMeshWithData) = hanging_to_masters(a.mesh)
hanging_to_coeffs(a::AbstractMeshWithData) = hanging_to_coeffs(a.mesh)
physical_groups(a::AbstractMeshWithData,d) = physical_groups(a.mesh,d)

new_mesh(args...) = GenericMesh(args...)
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

physical_group(args...) = GenericPhysicalGroup(args...)
struct GenericPhysicalGroup{A}
    faces_in_group::A
    group_name::String
end

faces_in_group(a) = a.faces_in_group
group_name(a) = a.group_name

set_phyisical_groups(mesh,groups) = MeshWithPhysicalGroups(mesh,groups)
struct MeshWithPhysicalGroups{A,B} <: AbstractMeshWithData{A}
    mesh::A
    physical_groups::Vector{Dict{Int,B}}
end
has_physical_groups(::Type{<:MeshWithPhysicalGroups}) = true
physical_groups(a::MeshWithPhysicalGroups,d) = a.physical_groups[d+1]

periodic_node_constraints(args...) = GenericPeriodicNodeConstraints(args...)
struct GenericPeriodicNodeConstraints{A,B,C}
    periodic_nodes::A
    periodic_to_master::B
    periodic_to_coeff::C
end
has_periodic_nodes(::Type{<:GenericPeriodicNodeConstraints}) = true

set_periodic_node_constraints(mesh,constraints) = MeshWithPeriodicNodeConstraints(mesh,constraints)
struct MeshWithPeriodicNodeConstraints{A,B} <: AbstractMeshWithData{A}
    mesh::A
    constraints::B
end
has_periodic_nodes(::Type{<:MeshWithPeriodicNodeConstraints{A,B}}) where {A,B} = has_periodic_nodes(B)
periodic_nodes(a::MeshWithPeriodicNodeConstraints) = a.constraints.periodic_nodes
periodic_to_master(a::MeshWithPeriodicNodeConstraints) = a.constraints.periodic_to_master
periodic_to_coeff(a::MeshWithPeriodicNodeConstraints) = a.constraints.periodic_to_coeff

set_haning_node_constraints(mesh,constraints) = MeshWithHangingNodeConstraints(mesh,constraints)
struct MeshWithHangingNodeConstraints{A,B} <: AbstractMeshWithData{A}
    mesh::A
    constraints::B
end
has_hanging_nodes(::Type{<:MeshWithHangingNodeConstraints{A,B}}) where {A,B} = has_hanging_nodes(B)
hanging_node_constraints(a::MeshWithHangingNodeConstraints) = a.constraints

struct HangingNodeConstraints{A,B,C,D,T} <: AbstractMatrix{T}
    maxid::Int
    free_nodes::A
    hanging_nodes::A
    master_nodes::B
    master_coeffs::C
    permutation::D
    function HangingNodeConstraints(maxid,free_nodes,hanging_nodes,master_nodes,master_coeffs,permutation)
        A = typeof(free_nodes)
        B = typeof(master_nodes)
        C = typeof(master_coeffs)
        D = typeof(permutation)
        T = eltype(eltype(master_coeffs))
        new{A,B,C,D,T}(maxid,free_nodes,hanging_nodes,master_nodes,master_coeffs,permutation)
    end
end

free_nodes(a) = a.free_nodes
hanging_nodes(a) = a.hanging_nodes
master_nodes(a) = a.master_nodes
master_coeffs(a) = a.master_coeffs
permutation(a) = a.permutation

function Base.size(a::HangingNodeConstraints)
    m = length(a.free_nodes) + length(a.hanging_nodes)
    n = a.maxid
    (m,n)
end

Base.IndexStyle(::Type{<:HangingNodeConstraints}) = IndexCartesian()

function Base.getindex(a::HangingNodeConstraints,i::Int,j::Int)
    T = eltype(a)
    m,n = size(a)
    @boundscheck @assert i in 1:m
    @boundscheck @assert j in 1:n
    p = a.permutation[i]
    nfree = length(a.free_nodes)
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
        a.free_nodes[p]==j ? one(T) : zero(T)
    end
end

## TODO implement mul and also for its transpose


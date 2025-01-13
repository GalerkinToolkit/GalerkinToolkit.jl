"""
    abstract type AbstractDomain <: AbstractType end

Abstract type representing the geometry of a single mesh face, typically one of the reference faces.

# Basic queries

- [`num_dims`](@ref)
- [`is_axis_aligned`](@ref)
- [`is_simplex`](@ref)
- [`is_n_cube`](@ref)
- [`is_unit_n_cube`](@ref)
- [`is_unit_simplex`](@ref)
- [`is_unitary`](@ref)
- [`bounding_box`](@ref)
- [`vertex_permutations`](@ref)
- [`mesh`](@ref)
- [`faces`](@ref)
- [`inverse_faces`](@ref)
- [`geometries`](@ref)
- [`topology`](@ref)
- [`geometries`](@ref)
- [`options`](@ref)
- [`is_boundary`](@ref)
- [`face_around`](@ref)

# Basic constructors

- [`unit_simplex`](@ref)
- [`unit_n_cube`](@ref)
- [`domain`](@ref)

"""
abstract type AbstractDomain <: AbstractType end

abstract type AbstractFaceDomain <: AbstractDomain end

"""
    abstract type AbstractMesh

# Basic queries

- [`node_coordinates`](@ref)
- [`face_nodes`](@ref)
- [`face_reference_id`](@ref)
- [`reference_spaces`](@ref)
- [`periodic_nodes`](@ref)
- [`physical_faces`](@ref)
- [`outward_normals`](@ref)

# Basic constructors

- [`mesh`](@ref)
- [`mesh_from_gmsh`](@ref)
- [`cartesian_mesh`](@ref)

"""
abstract type AbstractMesh end

abstract type AbstractChain <: AbstractType end

"""
    abstract type AbstractTopology

# Basic queries
- [`face_incidence`](@ref)
- [`face_reference_id`](@ref)
- [`face_permutation_ids`](@ref)
- [`reference_topologies`](@ref)
- [`vertex_permutations`](@ref)

# Basic constructors

- [`topology`](@ref)

"""
abstract type AbstractTopology <: AbstractType end

abstract type AbstractFaceTopology <: AbstractTopology end

abstract type AbstractMeshTopology <:AbstractTopology end

"""
    abstract type AbstractSpace <: AbstractType end

# Basic queries

[`domain`](@ref)
[`num_dofs`](@ref)
[`num_nodes`](@ref)
[`face_dofs`](@ref)
[`face_nodes`](@ref)
[`face_reference_id`](@ref)
[`reference_spaces`](@ref)
[`interior_nodes`](@ref)
[`interior_nodes_permutations`](@ref)
[`geometry_own_dofs`](@ref)
[`geometry_own_dofs_permutations`](@ref)
[`geometry_interior_nodes`](@ref)
[`geometry_interior_nodes_permutations`](@ref)
[`geometry_nodes`](@ref)
[`geometry_nodes_permutations`](@ref)

# Basic constructors

[`lagrange_space`](@ref)
[`raviart_thomas_space`](@ref)

"""
abstract type AbstractSpace <: AbstractType end

abstract type AbstractFaceSpace <: AbstractSpace end

"""
    abstract type AbstractQuadrature

# Basic queries

- [`domain`](@ref)
- [`coordinates`](@ref)
- [`weights`](@ref)
- [`num_points`]
- [`face_reference_id`](@ref)
- [`reference_quadratures`](@ref)

# Basic constructors

- [`quadrature`](@ref)
- [`duffy_quadrature`](@ref)
- [`tensor_product_quadrature`](@ref)
- [`node_quadrature`](@ref)

# Supertype hierarchy

    AbstractQuadrature <: GT.AbstractType
"""
abstract type AbstractQuadrature <: AbstractType end

abstract type AbstractFaceQuadrature <: AbstractQuadrature end

abstract type AbstractMeshQuadrature{A} <: AbstractQuadrature end

abstract type AbstractQuantity <: GT.AbstractType end

abstract type AbstractTerm <: GT.AbstractType end

# A quantity representing a field on a domain
abstract type AbstractField <: AbstractQuantity end


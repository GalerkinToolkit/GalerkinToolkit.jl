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

"""
"""
function num_dims end

"""
"""
function is_axis_aligned end

"""
"""
function is_simplex end

"""
"""
function is_n_cube end

"""
"""
function is_unit_n_cube end

"""
"""
function is_unitary end

"""
"""
function bounding_box end

"""
"""
function vertex_permutations end

"""
"""
function faces end

"""
"""
function inverse_faces end

"""
"""
function geometries end

"""
"""
function is_boundary end

"""
"""
function face_around end

"""
"""
function domain end

"""
"""
function node_coordinates end

"""
"""
function face_nodes end

"""
"""
function face_reference_id end

"""
"""
function reference_spaces end

"""
"""
function periodic_nodes end

"""
"""
function physical_faces end

"""
"""
function outward_normals end

"""
"""
function coordinates end

"""
"""
function weights end

"""
"""
function reference_quadratures end

"""
"""
function node_quadrature end

"""
"""
function num_dofs end

"""
"""
function interior_nodes end

"""
"""
function num_nodes end

"""
"""
function face_dofs end

"""
"""
function interior_nodes_permutations end

"""
"""
function geometry_own_dofs end

"""
"""
function geometry_own_dofs_permutations end

"""
"""
function geometry_interior_nodes end

"""
"""
function geometry_interior_nodes_permutations end

"""
"""
function geometry_nodes end

"""
"""
function geometry_nodes_permutations end

"""
"""
function face_incidence end

"""
"""
function face_permutation_ids end

"""
"""
function reference_topologies end


@doc raw"""
    abstract type AbstractDomain <: AbstractType end

Abstract type representing a subset of $\mathbb{R}^d$, typically $d\in\{0,1,2,3\}$.
Domains are defined using an underlying computational mesh.

See also [`AbstractMesh`](@ref).

# Level

Beginner

# Basic constructors

- [`unit_simplex`](@ref)
- [`unit_n_cube`](@ref)
- [`domain`](@ref)
- [`interior`](@ref)
- [`boundary`](@ref)
- [`skeleton`](@ref)

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

"""
abstract type AbstractDomain <: AbstractType end

"""
    abstract type AbstractFaceDomain <: AbstractDomain end

A domain defined on a single mesh face. Typically used as helper to identify cases that only make sense for a single mesh face.

# Level

Advanced

# Basic constructors

- [`unit_simplex`](@ref)
- [`unit_n_cube`](@ref)

"""
abstract type AbstractFaceDomain <: AbstractDomain end

@doc raw"""
    abstract type AbstractMesh <: AbstractType end

Abstract type representing a triangulation of a subset of $\mathbb{R}^d$, typically $d\in\{0,1,2,3\}$,
plus metadata useful in finite element computations, such as physical groups for imposing boundary conditions.

# Notation

Each of of the elements of the triangulation is referred to as a `face`. A mesh can contain faces of different dimensions.
A mesh might or might not represent a cell complex (all possible low dimensional faces are present in the mesh), but
it is often assumed that it represents a cell complex.

# Level

Beginner

# Basic constructors

- [`mesh`](@ref)
- [`chain`](@ref)
- [`mesh_from_gmsh`](@ref)
- [`mesh_from_space`](@ref)
- [`cartesian_mesh`](@ref)
- [`complexify`](@ref)
- [`simplexify`](@ref)

# Basic queries

- [`num_dims`](@ref)
- [`num_faces`](@ref)
- [`num_nodes`](@ref)
- [`node_coordinates`](@ref)
- [`face_nodes`](@ref)
- [`face_reference_id`](@ref)
- [`reference_spaces`](@ref)
- [`periodic_nodes`](@ref)
- [`physical_faces`](@ref)
- [`outward_normals`](@ref)

"""
abstract type AbstractMesh <: AbstractType end

"""
    abstract type AbstractTopology

Abstract type representing the incidence relations in a cell complex.

See also [`AbstractFaceTopology`](@ref).

# Level

Intermediate

# Basic constructors

- [`topology`](@ref)

# Basic queries

- [`face_incidence`](@ref)
- [`face_reference_id`](@ref)
- [`face_permutation_ids`](@ref)
- [`reference_topologies`](@ref)
- [`vertex_permutations`](@ref)

"""
abstract type AbstractTopology <: AbstractType end

"""
abstract type AbstractFaceTopology <: AbstractTopology end

Like [`AbstractTopology`](@ref), but for a single mesh face. Typically used as helper to identify cases that only make sense for a single mesh face.

# Level

Advanced
"""
abstract type AbstractFaceTopology <: AbstractTopology end

"""
    abstract type AbstractSpace <: AbstractType end

Abstract type representing a finite element space.

# Level

Basic

# Basic constructors

[`lagrange_space`](@ref)
[`raviart_thomas_space`](@ref)

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
"""
abstract type AbstractSpace <: AbstractType end

"""
abstract type AbstractFaceSpace <: AbstractSpace end

Like [`AbstractSpace`](@ref), but for a single mesh face. Typically used as helper to identify cases that only make sense for a single mesh face.

# Level

Advanced
"""
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
    num_dims(x)

Return the parametric dimension of `x`.

See also [`num_ambient_dims`](@ref), [`num_codims`](@ref).

# Level

Beginner
"""
function num_dims end

"""
    num_ambient_dims(x)

Return the ambient dimension where object `x` lives.

See also [`num_codims`](@ref), [`num_dims`](@ref).

# Level

Beginner
"""
function num_ambient_dims end

"""
    num_codims(x)

Return `num_ambient_dims(x)-num_dims(x)`.

See also [`num_ambient_dims`](@ref), [`num_dims`](@ref).

# Level

Beginner
"""
function num_codims end

num_codims(x) = num_ambient_dims(x) - num_dims(x)

"""
    is_axis_aligned(x)

True if `x` is a unit simplex or a unit cube.

# Level

Advanced
"""
function is_axis_aligned end

"""
    is_simplex(x)

True if `x` is a simplex.

See also [`is_n_cube`](@ref), [`is_unit_simplex`](@ref), [`is_unit_n_cube`](@ref).

# Level

Intermediate
"""
function is_simplex end

"""
    is_n_cube(x)

True if `x` is a n-cube (hypercube).

See also [`is_simplex`](@ref), [`is_unit_simplex`](@ref), [`is_unit_n_cube`](@ref).

# Level

Intermediate
"""
function is_n_cube end

"""
    is_unit_n_cube(x)

True if `x` is a unit n-cube.

See also [`is_n_cube`](@ref), [`is_unit_simplex`](@ref), [`is_simplex`](@ref).

# Level

Intermediate
"""
function is_unit_n_cube end

"""
    is_unit_simplex(x)

True if `x` is a unit simplex.

See also [`is_n_cube`](@ref), [`is_unit_n_cube`](@ref), [`is_simplex`](@ref).

# Level

Intermediate
"""
function is_unit_simplex end

"""
    is_unitary(x)

True `bounding_box(x)` coincides with a unit n-cube.

# Level

Advanced
"""
function is_unitary end

"""
    p0,p1 = bounding_box(x)

Return a tuple of two vectors, where the vectors `p0` and `p1`
define the span of the bounding box of `x`.

# Level

Intermediate
"""
function bounding_box end

"""
    geometries(x,d)
    geometries(x,Val(d))

Return a vector of domains representing the geometrical entities of `x` of dimension `d`. The returned domains
and `x` are defined in the same mesh. That is, `faces(geometries(x,1)[1])` are the face ids in `mesh(x)` representing the
first edge of `x`.

# Notation

`geometries(x,Val(0))` returns the vertices of `x`.
`geometries(x,Val(1))` returns the edges of `x`.
`geometries(x,Val(d))` returns the `d`-faces of `x`.

# Level

Advanced
"""
function geometries end

"""
    vertex_permutations(x)

Return a list of permutations representing the admissible re-labelings of the vertices
of `x`.

# Level

Advanced
"""
function vertex_permutations end

"""
    faces(x)

Return the subset of face ids in `mesh(x)` of dimension `num_dims(x)` defining the domain `x`.
This is effectively the map from domain face id to mesh face id.

See also [`inverse_faces`](@ref).

# Level

Intermediate
"""
function faces end

"""
    inverse_faces(x)

Return the inverse integer of `faces(x)`. This is effectively the map from mesh face id to domain face id.

See also [`faces`](@ref).

# Level

Intermediate
"""
function inverse_faces end

"""
    is_boundary(x)

True if `x` represent an (internal) boundary. Faces in an internal boundary "point" to only of the two faces around of a dimension higher.

See also [`face_around`](@ref).

# Level

Intermediate
"""
function is_boundary end

"""
    face_around(x)

Return an integer that allows to break ties when faces in `x` need to point to faces around of one dimension higher.
Return nothing otherwise.

Note: This function will eventually return a vector of integers.

# Level

Intermediate
"""
function face_around end

"""
    mesh_from_space(space)

Return the mesh induced by `space`. For instance, a (high order) Lagrange space can be interpreted as a mesh
using this function.
"""
function mesh_from_space end

"""
"""
function mesh_from_chain end

"""
"""
function complexify end

"""
"""
function simplexify end

"""
"""
function domain end

"""
"""
function boundary end

"""
"""
function interior end

"""
"""
function skeleton end

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


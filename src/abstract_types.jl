abstract type NewAbstractAccessor <: AbstractType end

abstract type AbstractAccessor <: Function end

function Base.show(io::IO,data::GT.AbstractAccessor)
    print(io,"GalerkinToolkit.$(nameof(typeof(data)))(…)")
end

@doc raw"""
    abstract type AbstractMesh <: AbstractType end

Abstract type representing a computational mesh.

# Level

Beginner

# Constructors

- [`create_mesh`](@ref)
- [`mesh_from_msh`](@ref)
- [`mesh_from_gmsh`](@ref)
- [`cartesian_mesh`](@ref)
- [`moebius_strip`](@ref)
- [`complexify`](@ref)
- [`simplexify`](@ref)

# Iteration

- [`each_face`](@ref)

# Data retrival

- [`num_dims`](@ref)
- [`num_ambient_dims`](@ref)
- [`num_codims`](@ref)
- [`num_faces`](@ref)
- [`num_nodes`](@ref)
- [`node_coordinates`](@ref)
- [`face_nodes`](@ref)
- [`face_reference_id`](@ref)
- [`reference_spaces`](@ref)
- [`group_faces`](@ref)
- [`group_names`](@ref)
- [`is_cell_complex`](@ref)
- [`normals`](@ref)

"""
abstract type AbstractMesh <: AbstractType end

#"""
#"""
#abstract type AbstractPMesh <: AbstractMesh end

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
- [`num_ambient_dims`](@ref)
- [`num_codims`](@ref)
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

"""
    abstract type AbstractTopology

Abstract type representing the incidence relations in a cell complex.

# Level

Intermediate

# Constructors

- [`topology`](@ref)

# Queries

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

- [`domain`](@ref)
- [`num_dofs`](@ref)
- [`face_dofs`](@ref)
- [`face_nodes`](@ref)
- [`face_reference_id`](@ref)
- [`reference_spaces`](@ref)
- [`geometry_own_dofs`](@ref)
- [`geometry_own_dofs_permutations`](@ref)

# Additional queries

For spaces, used as reference spaces in [`AbstractMesh`](@ref) specializations.

- [`num_nodes`](@ref)
- [`interior_nodes`](@ref)
- [`interior_nodes_permutations`](@ref)
- [`geometry_interior_nodes`](@ref)
- [`geometry_interior_nodes_permutations`](@ref)
- [`geometry_nodes`](@ref)
- [`geometry_nodes_permutations`](@ref)
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
- [`num_points`](@ref)
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

abstract type AbstractMeshQuadrature <: AbstractQuadrature end

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
    num_faces(x)
    num_faces(x,d)

Return the number of faces  of dimension `d` in `mesh(x)`.
If `d` is omitted, return a vector with the number of faces in each dimension,
starting from dimension 0 up to `num_dims(x)`.

# Level

Beginner
"""
function num_faces end

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
and `x` are defined on the same mesh. That is, `faces(geometries(x,1)[2])` are the face ids in `mesh(x)` representing the
second edge of `x`.

# Notation

`geometries(x,Val(0))` are referred to as the vertices of `x`.
`geometries(x,Val(1))` are referred to as the edges of `x`.
`geometries(x,Val(d))` are referred to as the `d`-faces of `x`.

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

Return the inverse integer mas of `faces(x)`. This is effectively the map from mesh face id to domain face id.
Mesh faces not present in the domain, receive an invalid index id.

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
Return nothing if `x` does not break such ties.

Note: This function will eventually return a vector of integers.

# Level

Intermediate
"""
function face_around end

"""
    mesh_from_space(space)

Return the mesh induced by `space`. For instance, a (high order) Lagrange space can be interpreted as a mesh
using this function.

# Level

Advanced
"""
function mesh_from_space end

"""
    complexify(x)

Convert `x` into a mesh representing a cell complex.

See also [`is_cell_complex`](@ref).

# Level

Intermediate
"""
function complexify end

"""
    simplexify(x)

Convert `x` into a mesh made of simplex cells.

# Level

Intermediate
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
    node_coordinates(x)

Return the vector of node coordinates associated with `x.
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
    group_faces(mesh)
    group_faces(mesh,d)

Return the dictionary containing the faces in each group in dimension `d`.
If `d` is omitted, it returns the dictionaries for all dimensions in a vector. I.e.,
calling `group_faces(mesh,d)` is equivalent to `group_faces(mesh)[d+1]`.

The faces of dimension `d` in group `group` are `group_faces(mesh,d)[group_name]`, where `group_name` is a string with the group name.
One can create new groups by adding new keys to these dictionaries as long as the key is not already present.
Calling `group_faces(mesh,d)[new_group_name] = faces_in_newgroup` will add a new group to dimension `d` with name equal to the string `new_group_name` with faces in vector `faces_in_newgroup`.

See also [`group_names`](@ref).

# Level

Beginner

"""
function group_faces end

"""
    group_names(mesh)
    group_names(mesh,d)

Return the names of the groups in dimension `d`. Calling `group_names(mesh,d)` is equivalent to `group_names(mesh)[d+1]`.

See also [`group_faces`](@ref).

# Level

Beginner


"""
function group_names end

"""
"""
function is_cell_complex end

"""
"""
function normals end

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
    num_dofs(x)

Return the number of degrees of freedom of `x`.

See also [`num_nodes`](@ref).

# Level

Beginner
"""
function num_dofs end

"""
"""
function interior_nodes end

"""
    num_nodes(x)

Return the number of nodes of `x`.

See also [`num_dofs`](@ref).

# Level

Beginner
"""
function num_nodes end

"""
    num_points(x)

Return the number of integration points in `x`.

See also [`num_nodes`](@ref).

# Level

Beginner
"""
function num_points end

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
    face_incidence(x,d)
"""
function face_incidence end

"""
"""
function face_permutation_ids end

"""
    reference_topologies(topo)
    reference_topologies(topo,d)
    reference_topologies(topo,Val(d))

Return the list (a vector or a tuple) of reference topologies in `topo` of dimension `d`.
If the second argument is omitted, then the function returns a collection such that `reference_topologies(topo)[d+1]` is equivalent to `reference_topologies(topo,Val(d))`.

The face reference topology of face `f` of dimension `d`, is accessed as `reference_topologies(topo,d)[r]` with `r=face_reference_id(topo,d)[f]`.

See also [`face_reference_id`](@ref).

# Level

Intermediate
"""
function reference_topologies end


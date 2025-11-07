# ```@meta
# CurrentModule=GalerkinToolkit
# ```
# # Mesh geometry
#
# !!! note
#     TODOs
#     - `AbstractMeshFace` is currently called `MeshAccessor`.
#     - `AbstractMeshFace` should be `AbstractDomain` or `AbstractAccessor`  ?
#     - API for the referene face. Now a reference face `F` isa `AbstractDomain`, and not the equivalent of `AbstractMeshFace`.
#     - Reference face vs reference domain?
#     - face complex is currently called cell complex.
#
# A mesh object in GalerkinToolkit contains all geometrical information needed in a finite element (FE) computation.
# This includes the discretization of computational domains as well as data to impose different types of boundary conditions.
# It is worth noting that GalerkinToolkit is not a mesh generation library.
# The mesh API is designed to provide the rich geometrical information needed in FE methods, rather than mesh generation.
# Meshes are often generated with external tools and then transformed into GalerkinToolkit objects with helper functions such as [`mesh_from_gmsh`](@ref).
# In this page, we discuss the geometry representation of computational meshes in GalerkinToolkit and in particular:
#
# - the definition of mesh and its code representation,
# - reference and physical faces,
# - the physical map,
# - face, node, and reference ids, and
# - face groups.
#
# In this page, we assume that
# we have already created a mesh object.
# The construction of mesh object is not detailed here, but in Section Mesh Generation.
#
# ## Definition
#
# A *mesh* $M$ in GalerkinToolkit is defined as set of *physical faces*
# embedded in the Euclidean space $\mathbb{R}^D$, with $D$ often being $D=1,2,3$.
# A physical face is defined in terms of a reference face 
# and a reference space as detailed later below.
# In the API, a mesh $M$ is represented with a mesh object `M`,
# whose type is a subtype of [`AbstractMesh`](@ref). Even though our math notation
# defines a mesh $M$ as a set, a mesh object `M::AbstractMesh` has not the API
# of a set, but rich API providing the information encoding the set of faces $M$.
# There is a one-to-one relation between the mathematical mesh $M$ and the API
# mesh `M` and we often refer to them simply as "a mesh". The same is true for other
# math definitions and their corresponding API.
#
# A face $F\in M$ in a mesh is represented in the code with an object `F::AbstractMeshFace`.
# Given `M::AbstractMesh`, function `GT.each_face(M,d)` creates an iterator
# used to traverse all faces of dimensions `d` in mesh `M`. 
# Using this iterator, faces can be accessed with the 
# Julia loop syntax.
# ```julia
# for F in GT.each_face(M,d)
#     # F isa AbstractMeshFace
# end
# ```
#
#
# ## Dimensions
#
# Faces $F\in M$ in a mesh are *open* $d$-dimensional manifolds embedded in the Euclidean space $\mathbb{R}^D$. We call 
# $D$ the number of *ambient* dimensions of the mesh $M$  and of faces $F\in M$, where as
# $d$ is the number of dimensions of face $F$, which might be $d=0,\ldots,D$. We define the number of
# *co-dimensions* as the number of ambient dimensions minus the number of dimensions. 
# A mesh might contain faces of different dimensions and we define the number of dimensions of a mesh
# as the maximum number of dimensions
# of their faces. In the API, the number of dimensions, ambient dimensions, and co-dimensions are
# obtained with `num_dims(X)`, `num_ambient_dims(X)`, and `num_codims(X)` respectively for an instance
# `X::AbstractMeshFace` or `X::AbstractMesh`.
#
# This is some extra notation that we
# often use in the library. We call a $d$-face to a face
# of $d$ dimensions. We call *vertices*, *edges*, *surfaces*, and *volumes* to faces
# of 0, 1, 2, and 3 dimensions respectively. We call *chain* to a mesh, whose faces are all of the same dimension.
# The next figure shows a two-dimensional mesh embedded in a three-dimensional space.
# This mesh contains vertices, edges, and surfaces, but not volumes. This mesh is *not* a chain, but a *face complex*
# (see section [Mesh topology](@ref)).
#
# ![](fig_meshes_defs_1.png)
#
# **Figure:** Visualization of the faces in the mesh of a Möbius strip.
# We shrink
# the mesh faces to illustrate that faces are open sets.
# Otherwise faces of hider dimensions would hide faces of lower dimensions
# in the figure.
#
module MeshGeo1 # hide
#Code used to generate the figure
import GalerkinToolkit as GT
import GLMakie as Makie
import FileIO # hide
cells = (4,40)
mesh = GT.moebius_strip(cells;width=0.6)
fig = Makie.Figure()
elevation = 0.24π
azimuth = -0.55π
aspect = :data
ax = Makie.Axis3(fig[1,1];aspect,elevation,azimuth)
shrink = 0.8
shading = Makie.NoShading
GT.makie_surfaces!(mesh;shrink,shading,dim=2)
GT.makie_edges!(mesh;shrink,dim=1)
GT.makie_vertices!(mesh;shrink,dim=0)
FileIO.save(joinpath(@__DIR__,"fig_meshes_defs_1.png"),Makie.current_figure()) # hide
end # hide
nothing # hide

#
# ## Face ids
#
# For each $d$, we assign a unique integer in to each $d$-face in $M$, called the *face id*.
# Face ids are assigned per dimension (two faces of different dimension might have the same id).
# Thus, a face is uniquely identified by its face id *and* its dimension $d$.
# The face ids are arbitrary as long as they are consecutive integers starting by one.
# In the API, `num_faces(M,d)` returns the number
# of `d`-faces and `faces(M,d)` is the range `1:num_faces(M,d)` containing all face ids
# in dimension `d`. There is a one-to-one relation between face
# objects and face ids.
# For a face object `F::AbstractMeshFace`,
# `face=id(F)` is the face id of `F` and `F=each_face(M,d)[face]` is the face `F`
# for a face id `face::Integer`.
# We often refer to face objects `F` and face ids `face` simply
# as "a face" since they are equivalent.
# In the API, the same operation can be done often in two different
# ways, one using face ids (integers) and another using face objects.
#
# ## Node ids
#
# Like in many other FE codes, the node coordinates of a face are
# encoded here using the vector of node coordinates of the mesh and
# the *node ids* of the face. 
# The vector of node coordinates of a mesh is accessed
# with `node_coordinates(M)` for `M::AbstractMesh`. The length of this vector is
# `num_nodes(M)` and `nodes(M)` is the range `1:num_nodes(M)` containing all possible
# node ids.
#
# For a given `d`, we collect the node ids of all `d`-faces in a mesh `M` in the vector returned
# by `face_nodes(M,d)`. 
# The node ids of a face with id `face::Integer` and dimension `d` are accessed with
# `face_nodes(M,d)[face]`. 
# The vector `face_nodes(M,d)` is often called  the *local-to-global* (index) map or the face *connectivity*.
# We call them the face node ids.
#  The value `g = face_nodes(M,d)[face][l]` is called the *global* node id
#  associated with the *face-local* (or simply local) node id `l` in face `face`.
#  The coordinates of the local node `l` are then computed by indexing the mesh coordinates
#  with the global node id `g`, namely
#  `node_coordinates(M)[g]`.
#  The vector `face_nodes(M,d)` is a long vector of small vectors of integers with possibly different lengths.
#  It is often represented using a `PartitionedArrays.JaggedArray` object that uses continuous linear memory for performance.
#
#  It is also possible to access
#  node ids and node coordinates from a face object `F::AbstractMeshFace` as `nodes(F)`
#  and `node_coordinates(F)` respectively. Note that 
#  `nodes(F)` and `face_nodes(M,num_dims(F))[id(F)]` are equivalent.
#
#  ## Reference ids
#
# Each physical face in a mesh is defined by means of a reference FE space. Several faces often share the
# same reference space and, in the limit, all faces of the same dimensions share the same reference
# space. For each dimension `d`, we collect the unique reference spaces of faces of dimension `d`
# in a tuple. This tuple is returned by `reference_spaces(M,d)` for a mesh object `M::AbstractMesh`.
#
# The reference space assigned to a face is then obtained as
# `reference_spaces(M,d)[r]`, where `r` is called the *reference* id of the face.
# The reference id is obtained from a face id `face::Integer`, by indexing the vector
# `face_reference_id(M,d)`, namely `r = face_reference_id(M,d)[face]`.
# It can also be obtained from a face object `F::AbstractMeshFace`
# with `r=reference_id(F)`, and the corresponding reference space with `reference_space(F)`.
# Note that `reference_id(F)`
# and `face_reference_id(M,num_dims(F))[id(F)]` are equivalent.
#
#  The notion of reference id is introduced since different face
#  typologies such as simplices and hyper-cubes might be in the same mesh.
#  If all $d$-faces are topologically equivalent (which is often the case),
#  there is only one reference  space for all $d$-faces and their reference id is one.
#
#  ## Reference spaces
#
# An item `Vref` in `reference_spaces(M,d)` is of a type that specializes the [`AbstractSpace`](@ref)
# interface.
# The `AbstractSpace` interface 
# is detailed in section [Interpolation](@ref). In this page,
# we only need to consider that a reference space has a basis of scalar shape functions,
# which is accessed as `shape_functions(Vref)`. The result is a vector of functions, where
# each function is evaluated at a point `x::AbstractVector (often x::SVector)`,
# and returns a scalar value `s::Real`.
# 
# A reference $\hat V$ space is a scalar-valued (possibly high-order) Lagrange FE space 
# defined on a reference face $\hat F$. Reference spaces can be get from a mesh object
# as shown above or created from scratch with function [lagrange_space](@ref).
# E.g., `Vref = lagrange_space(Fref,order)` creates a reference space
# of order `order` on the reference face `Fref`.
#
# ## Reference faces
#
# A reference $d$-face $\hat F$ is a $d$-dimensional [polytope](https://en.wikipedia.org/wiki/Polytope)
# embedded in the Euclidian space $\mathbb{R}^d$. In particular, $\hat F$ is a segment, polygon, and a polyhedron for $d=1,2,3$ respectively.
# For $d=0$, we define a reference face $\hat F:=\{v\}$ as a set containing the only point  $v\in\mathbb{R}^0$. For $d>0$.
# We define the boundary $\partial\hat F$ of a reference $d$-face $\hat F$
# as the union $\partial\hat F := U_{f\in \hat M_{d-1}} \bar f$,
# where $\hat M_{d-1}$ is a chain of faces of dimension $d-1$
# and $\bar f$ is the [closure](https://en.wikipedia.org/wiki/Closure_(topology)) of a face $f$.
# E.g., the boundary of a segment is the union of two vertices. The boundary of a square is the union
# of four segments, and the boundary of a cube is the union of four squares.
# Assuming that $\partial\hat F$ is closed, we define the reference face $\hat F$ as
# the open bounded subset of $\mathbb{R}^d$ with boundary $\partial\hat F$. 
# i.e., $\hat F$ is the space "inside" $\partial\hat F$.
#
# Starting from the chain $\hat M_{d-1}$, we create a mesh $\hat M$ that contains $n$-faces of dimensions
# $n=0,\ldots,d$ as follows.
# Given two faces $A$ and $B$, let $\Gamma_{A,B}$ be the *interface* of $A$ and $B$. We
# define the interface
# $\Gamma_{AB}$ as the open set such that $\bar \Gamma_{A,B}=\bar A \cap \bar B$. With this notation
# we define the mesh $\hat M$ as the set containing
# 1. the reference face $\hat F$,
# 2. the faces in the chain $\hat M_{d-1}$,
# 3. and all possible interfaces $\Gamma_{A,B}$ for $A$ and $B$ in $\hat M_{d-1}$.
#
# The set $\hat M$ contains all faces in the boundary of the polytope $\hat F$. E.g.,
# for a reference edge, it contains the edge and the two end vertices.
# For a reference square, it contains the square, four edges, and the four vertices
# at the intersection of the edges. For a reference cube, it contains the cube, six surfaces, twelve
# edges, and eight vertices (see the next figure).  We call the faces in $\hat M$ the
# *face-local* faces (or simply the local faces) of $\hat F$.
#
# In the API, reference faces are often built as `Fref = unit_n_cube(Val(d))` or
# `Fref = unit_simplex(Val(d))` that create a unit hypercube or a unit simplex
# of `d` dimensions respectively. One can also access the reference face associated with
# a physical face in a mesh, by getting the reference space, and then the domain of this space, e.g.,
# `Fref = domain(reference_space(F))`.
#
# The type of the returned object `Fref`
# is a subtype of [`AbstractDomain`](@ref). One can access the mesh of local faces
# from the reference face `Fref` with `Mref = mesh(Fref)`. The mesh
# `Mref::AbstractMesh` is like any other mesh used in the code and, e.g., it can be visualized as
# any other mesh shown in next Figure. See section [Domains](@ref) for information
# about the `AbstractDomain` interface.
# 
# ![](fig_meshes_defs_5.png)
#
# **Figure**: Visualization of the local faces of a reference cube, square, and segment.
# Faces are shrunk in the visualization to illustrate that they are open sets.
#
module Meshes05 # hide
#Code used to generate the figure
import GalerkinToolkit as GT
import GLMakie as Makie
import FileIO # hide
domain = (0,2,0,1,0,2)
cells = (4,2,4)
F3 = GT.unit_n_cube(Val(3))
F2 = GT.unit_n_cube(Val(2))
F1 = GT.unit_n_cube(Val(1))
M3 = GT.mesh(F3)
M2 = GT.mesh(F2)
M1 = GT.mesh(F1)
fig = Makie.Figure()
aspect = :data
ax = Makie.Axis3(fig[1,1];aspect)
Makie.hidespines!(ax)
Makie.hidedecorations!(ax)
shrink = 0.8
GT.makie_surfaces!(ax,M3;dim=2,shrink)
GT.makie_surfaces!(ax,M3;dim=3,shrink)
GT.makie_edges!(ax,M3;dim=1,shrink)
GT.makie_vertices!(ax,M3;dim=0,shrink)
aspect = Makie.DataAspect()
axis = (;aspect)
ax, = GT.makie_surfaces(fig[1,2],M2;dim=2,shrink,axis)
Makie.hidespines!(ax)
Makie.hidedecorations!(ax)
GT.makie_edges!(M2;dim=1,shrink)
GT.makie_vertices!(M2;dim=0,shrink)
ax = Makie.Axis(fig[1,3])
Makie.hidespines!(ax)
Makie.hidedecorations!(ax)
GT.makie_edges!(ax,M1;dim=1,shrink)
GT.makie_vertices!(ax,M1;dim=0,shrink)
FileIO.save(joinpath(@__DIR__,"fig_meshes_defs_5.png"),fig) # hide
end # hide
nothing # hide
#
#
#  ## Physical faces
#
# A physical face $F\in M$ is defined as the image $\phi(\hat F)$ of its reference face $\hat F$ via a map $\phi: \mathbb{R}^d \rightarrow \mathbb{R}^D$, where $d$ is the dimension of $F$ and $D$ is the ambient dimension of the mesh, see next Figure. This map is called the *physical map* and it is defined using
# the reference space $\hat V$ 
# and the node coordinates of the face.
#
# The physical map is defined as follows:
# 
# $\phi(\hat x) := \sum_{n=1}^{N^F} x^F_n s^{\hat V}_n(\hat x),$
#
# where $N^F$, $x^F_n$, and $s^{\hat V}_n$ are the number of nodes in face `F`, the coordinate vector
# of the local node $n$ if face $F$ and the shape function number $n$ in the reference space $\hat V$.
# These quantities can be programmatically  obtained, e.g., as
# `num_nodes(F)`, `node_coordinates(F)[n]` and `shape_functions(reference_space(F))[n]`
# for a mesh object `F::AbstractMeshFace` an integer `n`.
#
# ![](fig_meshes_defs_3.png)
#
# **Figure:** Effect of mapping a reference cube with a third order physical map.
# Orange dots illustrate the nodes before and after the map.

module MeshGeo3 # hide
#Code used to generate the figure
import GalerkinToolkit as GT
import GLMakie as Makie
import StaticArrays
import FileIO # hide
order = 3
Fref = GT.unit_n_cube(Val(3))
Mref = GT.mesh(Fref)
Vref = GT.lagrange_space(Fref,order)
xref = GT.node_coordinates(Vref)
T = StaticArrays.SVector{3,Float64}
x = map(y->y + 0.1*(rand(T) .- 0.5) ,xref)
node_coordinates = x
N = GT.num_nodes(Vref)
face_nodes = [collect(1:N)]
reference_spaces = (Vref,)
M = GT.create_chain(;
    node_coordinates,
    face_nodes,
    reference_spaces,
   )
F = GT.domain(M,3)
fig = Makie.Figure()
aspect = :data
ax = Makie.Axis3(fig[1,1];aspect)
GT.makie_surfaces!(ax,Mref)
Makie.scatter!(xref;color=Makie.Cycled(2))
ax = Makie.Axis3(fig[1,2];aspect)
refinement = 30
GT.makie_surfaces!(ax,F;refinement)
Makie.scatter!(x;color=Makie.Cycled(2))
FileIO.save(joinpath(@__DIR__,"fig_meshes_defs_3.png"),Makie.current_figure()) # hide
end # hide
nothing # hide
#
#
#
#
# ## Face groups
#
# For a given $d$, we call a *face group* to a subset 
# $G\subset M$ of the $d$-faces in a mesh $M$.
# A mesh is typically endowed with several
# of these groups to identify particular faces of the mesh for modeling purposes,
# e.g., to impose boundary conditions, or define position-dependent material properties.
# Each group is given a *group name*, which identifies the group.
#
# In the API, `group_faces(M,d)` provides access to the face groups for faces of dimension `d`.
# It is a Julia `Dict`. The keys are `String` objects for the group names, and the value
# `group_faces(M,d)[group]` is a vector of integers containing the face ids for faces inside the group
# with name `name`. The keys of this dictionary can also be accessed as
# `group_names(M,d)`. 
# Face groups are defined per dimension
# and it is accepted to have the same group name in two or more dimensions.
# It is also possible to add new groups by adding new key-value pairs to this
# dictionary.
#
# ### Example
#
# We illustrate how a new face group is added to an existing mesh.
# We create a new group with all 2-faces whose center is
# inside the ball centered at the origin and radius 1.
# The example uses part of the API described above to find
# the faces to be added in the group. We color code faces inside
# the group with value 1, and outside with value 0.

module MeshGeo4 # hide
import GalerkinToolkit as GT
import GLMakie as Makie
using LinearAlgebra
import FileIO # hide

#Create a Cartesian mesh mesh
domain = (0,1,-1,1,0,1)
cells = (10,20,10)
mesh = GT.cartesian_mesh(domain,cells)

#Find faces in new group
d = 3
mesh_faces = GT.each_face(mesh,d)
new_group_faces = findall(mesh_faces) do F
    xs = GT.node_coordinates(F)
    x = sum(xs) / length(xs)
    norm(x) < 1
end

#Add new group to mesh
GT.group_faces(mesh,d)["foo"] = new_group_faces

#Visualize
color = GT.FaceColor("foo")
fig = Makie.Figure()
aspect = :data
ax = Makie.Axis3(fig[1,1];aspect)
surfs = GT.makie_surfaces!(ax,mesh;dim=d,color)
GT.makie_edges!(ax,mesh;color=:black)
Makie.Colorbar(fig[1,2],surfs)
FileIO.save(joinpath(@__DIR__,"fig_meshes_defs_4.png"),Makie.current_figure()) # hide
end # hide
nothing # hide

# ![](fig_meshes_defs_4.png)
#
# ## Summary
#
# We discussed how computational meshes are defined in GalerkinToolkit
# and the core API of the `AbstractMesh` interface.
# See the docsting of [`AbstractMesh`](@ref) for a summary list of the API functions.
# Some of them will be introduced in other sections.
#
# ```@docs; canonical=false
# AbstractMesh
# ```
#

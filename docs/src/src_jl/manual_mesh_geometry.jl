# ```@meta
# CurrentModule=GalerkinToolkit
# ```
# # Mesh geometry
#
# A mesh object in GalerkinToolkit contains all geometrical information needed in a finite element (FE) computation.
# This includes the discretization of computational domains as well as data to impose different types of boundary conditions.
# It is worth noting that GalerkinToolkit is not a mesh generation library.
# The mesh API is designed to provide the rich geometrical information needed in FE methods, rather than mesh generation.
# Meshes are often generated with external tools and then transformed into GalerkinToolkit objects with helper functions such as [`mesh_from_gmsh`](@ref).
# In this page, we discuss the geometry representation of computational meshes in GalerkinToolkit.
#
# ## Definition
#
# A *mesh* $M$ in GalerkinToolkit is defined as set of *physical faces* embedded in the Euclidean space $\mathbb{R}^D$, with $D$ often being $D=1,2,3$. A physical face is defined in terms of a reference face 
# and a reference space as detailed later in this page.
# In the API, A mesh $M$ is represented with a mesh object `M`,
# whose type is a subtype of [`AbstractMesh`](@ref). Even though our math notation
# defines a mesh $M$ as a set, a mesh object `M::AbstractMesh` has not the API
# of a set, but an API providing the information encoding the set of faces $M$.
# Thus, there is a one-to-one relation between the mathematical mesh $M$ and the API
# mesh `M` and we often refer to them as the same thing. The same is true for other
# math definitions and their corresponding API.
#
# A face $F\in M$ in a mesh is represented in the code as an object `F::AbstractMeshFace` (TODO now it is called `MeshAccessor`).
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
# The construction of mesh object is detailed in Section Mesh Generation. In this page, we assume that
# we have already created a mesh object.
#
# ## Dimensions
#
# Faces $F\in M$ in a mesh are $d$-dimensional manifolds embedded in the Euclidean space $\mathbb{R}^D$. The value
# $D$ is called the number of *ambient* dimensions of the mesh $M$.
# We call $d$ the number of dimensions of face $F$, which might be $d=0,\ldots,D$.
# *co-dimensions* as the number of ambient dimensions minus the number of dimensions. For continence,
# we define the number of dimensions of a mesh as the maximum number of dimensions
# of their faces. In the API, the number of dimensions, ambient dimensions, and co-dimensions are
# obtained as `num_dims(X)`, `num_ambient_dims(X)`, and `num_codims(X)` respectively for an instance
# `X::AbstractMeshFace` or `X::AbstractMesh`.
#
#  We call $d$-face a face
# of $d$ dimensions. We call *vertices*, *edges*, *surfaces*, and *volumes* to faces of 0, 1, 2, and 3 dimensions respectively. We call *chain* to a mesh $M$, whose faces are all of the same dimension.
# The next figure shows a two-dimensional mesh embedded in a three-dimensional space.
# This mesh contains vertices, edges, and surfaces. This mesh is *not* a chain.
#
# ![](fig_meshes_defs_1.png)
#
# **Figure:** Visualization of the faces in the mesh of a Moebius strip.
# We shrink
# the mesh faces for visualization purposes.
# Otherwise faces of hider dimensions would hide faces of lower dimensions.
#
module MeshGeo1 # hide
#Code to generate the figure
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
# Given $d$, We assign a unique integer in to each $d$-face in $M$, called the *face id*.
# Note that face ids are assigned per dimension (two faces of different dimension might have the same id).
# A face is uniquely identified by its face id *and* its dimension $d$.
# The face ids are arbitrary as long as they are consecutive integers starting by one. Thus, the largest
# face id for $d$-faces is the number of $d$-faces. In the API, `num_faces(M,d)` returns the number
# of `d`-faces and `faces(M,d)` is the range `1:num_faces(M,d)`. For a face object `F::AbstractMeshFace`,
# `id(F)` is the face id of `F`.
#
# ## Reference faces
#
# A reference $d$-face $\hat F$ is a $d$-dimensional [polytope](https://en.wikipedia.org/wiki/Polytope)
# embedded in the Euclidian space $\mathbb{R}^d$. In particular, $\hat F$ is a segment, polygon, and a polyhedron for $d=1,2,3$ respectively.
# For $d=0$, we define a reference face $\hat F:=\{v\}$ as a set containing the only point  $v\in\mathbb{R}^0$. For $d>0$.
# We define the boundary $\partial\hat F$ of a reference $d$-face $\hat F$ as the union $\partial\hat F := U_{f\in \hat M_{d-1}} \bar f$,
# where $\hat M_{d-1}$ is a set of faces of dimension $(d-1)$ and $\bar f$ is the [closure](https://en.wikipedia.org/wiki/Closure_(topology)) of a face $f\in \hat M_{d-1}$.
# E.g., the boundary of a segment is the union of two vertices. The boundary of a square is the union of four segments and the boundary of a cube is the union of four squares (see next Figure).
# Assuming that $\partial\hat F$ is closed, we define the reference face $\hat F$ as
# the open bounded subset of $\mathbb{R}^d$ with boundary $\partial\hat F$.
#
# ![](fig_meshes_defs_2.png)
#
# **Figure:** Visualization of reference hypercubes of dimension 1,2,3,
# along the co-dimension 1 faces on the boundary.

module MeshGeo2 # hide
#Code to generate the figure
import GalerkinToolkit as GT
import GLMakie as Makie
import FileIO # hide
F1 = GT.unit_n_cube(Val(1))
F2 = GT.unit_n_cube(Val(2))
F3 = GT.unit_n_cube(Val(3))
M1 = GT.mesh(F1)
M2 = GT.mesh(F2)
M3 = GT.mesh(F3)
fig = Makie.Figure()
aspect = :data
shrink = 0.8
ax = Makie.Axis3(fig[1,1];aspect)
Makie.hidespines!(ax)
Makie.hidedecorations!(ax)
GT.makie_surfaces!(ax,M3;dim=3,shrink)
GT.makie_surfaces!(ax,M3;dim=2,shrink)
aspect = Makie.DataAspect()
axis = (;aspect)
ax, = GT.makie_surfaces(fig[1,2],M2;dim=2,shrink,axis)
Makie.hidespines!(ax)
Makie.hidedecorations!(ax)
GT.makie_edges!(M2;dim=1,shrink)
ax = Makie.Axis(fig[1,3])
Makie.hidespines!(ax)
Makie.hidedecorations!(ax)
GT.makie_edges!(ax,M1;dim=1,shrink)
GT.makie_vertices!(ax,M1;dim=0,shrink)
FileIO.save(joinpath(@__DIR__,"fig_meshes_defs_2.png"),fig) # hide
end # hide
nothing # hide
#
#
# ## Reference spaces
#
# Each face $F\in M$ is associated with a reference face $\hat F$ and a reference space $\hat V$.
# The reference $\hat V$ space is a scalar-valued (possibly high-order) Lagrange FE space 
# defined on a reference face $\hat F$. The full details
# about FE spaces are given in another section. Here, we only need
# to consider the shape functions of space $\hat V$. In the API, a FE
# space is represented with an instance `V::AbstractSpace`. From this object, we obtain
# the shape functions as `shape_functions(V).
#
# The reference space assigned a face `F::AbstractMeshFace` is obtained
# with `reference_space(F)`. The map between a face and its reference space is encoded using a
# compressed format, where several faces share the same space.  
# For a given dimension $d$, the unique reference spaces for the $d$-faces
# are collected in a tuple. In the API, this tuple is returned by `reference_spaces(M,d)` for
# a mesh `M::AbstractMesh`. The reference space assigned a face with id `face` is then obtained as
# `reference_spaces(M,d)[r]` with `r=face_reference_id(M,d)[face]`. The vector `face_reference_id(M,d)`
# The value `r` is called
# the *reference* id of face id `face`.
#  The notion of reference id is introduced since different face
#  typologies such as simplixes and hyper-cubes might be in the same mesh.
#  If all $d$-faces are topologically equivalent (which is often the case),
#  there is only one reference  space for all $d$-faces and their reference id is one.
#
#  ## Physical faces
#
# A physical face $F\in M$ is defined as the image $\phi(\hat F)$ of its reference face $\hat F$ via a map $\phi: \mathbb{R}^d \rightarrow \mathbb{R}^D$, where $d$ is the dimension of $F$ and $D$ is the ambient dimension of the mesh, see next Figure. This map is called the *physical map* and it is defined using
# the reference space $\hat V$ 
# and the *physical* nodes for this face. The coordinates for the physical nodes are in vector `node_coordinates(F)`
# for an instance `F::AbstractMeshFace`.
#
# The physical map is defined as follows:
# 
# $\phi(\hat x) := \sum_{n=1}^{N^F} x^F_n s^{\hat V}_n(\hat x),$
#
# where $N^F$, $x^F_n$, and $s^{\hat V}_n$ are programmatically  obtained as
# `num_nodes(F)`, `node_coordinates(F)[n]` and `shape_functions(reference_space(F))[n]`
# for an integer `n`.
#
# ![](fig_meshes_defs_3.png)
#
# **Figure:** Effect of mapping a reference cube with a third order physical map.
# Orange dots illustrate the nodes before and after the map.

module MeshGeo3 # hide
#Code to generate the figure
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
#

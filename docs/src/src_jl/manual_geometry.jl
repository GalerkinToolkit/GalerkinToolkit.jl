#
# ```@meta
# CurrentModule=GalerkinToolkit
# ```
#
# # Geometry
#
# We use the following dependencies in code snippets in this page.

import GalerkinToolkit as GT
import GLMakie
import Makie
using StaticArrays
import FileIO # hide

#
# ## Meshes
#
#  ### Definitions
#
# A mesh object in GalerkinToolkit contains all geometrical information needed in a FEM computation.#
# Formally:
#
# - A *mesh* $M$ is a set of faces plus a collection of *face labels* used to identify particular faces in the mesh, e.g., to impose boundary conditions.
# - A *face* $F \in M$ is a manifold, subset of $\mathbb{R}^D$, where $D$ is the number of *ambient* dimensions.
# - Each face $F \in M$ is called a *physical* face.
# - A physical face is defined as a transformation $\varphi(\hat F)$ of a *reference* face $\hat F$ via a map mapping $\varphi: \hat F \rightarrow \mathbb{R}^D$.
# - A reference face $\hat F$ is a subset of $\mathbb{R}^d$ with non-zero measure. In practice, it is  either a unit simplex of $d$ dimensions, or a unit hypercube of $d$-dimensions.
# - The integer $d$ is the number of (parametric) dimensions of the reference face $\hat F \subset \mathbb{R}^d$ and the corresponding physical face $F$.
# - A face of $d$ dimensions is called a $d$-face.
# - A *vertex* is a $0$-face, an *edge* is a $1$-face.
# - The number of dimensions of a mesh is the maximum number of dimensions of its faces.
# - The mapping $\varphi(\hat x) = \sum_i \hat s_i(\hat x) x_{(F,i)}$ is defined as a linear compinations of scalar-valued basis functions $\hat s_i: \hat F \rightarrow \mathbb{R}$ and coordinates $x_{(F,i)}\in\mathbb{R}^D$.
# - The space spanned by the reference basis functions $\hat s_i$ is called the *reference space*.
# - The vector $x_{(F,i)}$ is a *node* coordinate of the physical face $F$.
# - The face node coordinates are stored using two arrays as in many other FEM packages: a vector of node coordinates and an array of *face node ids* (also known as face connectivities), or simply *face nodes*.
#
# Important features of a mesh include:
#
# - The number of ambient dimensions is the same for all faces in a mesh.
#
# - A mesh can potentially contain faces with different number of dimensions. I.e., the same mesh object can include vertices, edges, $2$-faces, and $3$-faces.
#
# - The number of ambient dimensions might be larger or equal than the maximum dimension of the mesh faces. I.e., it is possible for a mesh embedded in a 3D space to contain faces only up to dimension $2$.
#
# - A mesh might or might not represent a cell complex (all possible low dimensional faces are present in the mesh). However, many algorithms require to work with a cell complex.
#
# - Multiple reference spaces and reference faces are allowed for faces of a given dimension. This is useful to support meshes with mixed face topologies and/or mixed interpolation order.
#
#
#  ### Interface
#
# All types implementing the meshes are subtypes of [`AbstractMesh`](@ref). All of them implement the *mesh interface*:
#
# ```@docs; canonical=false
# AbstractMesh
# ```
#
# ### Building a mesh from scratch.
#
# Function [`create_mesh`](@ref) builds a mesh object from its constituent data: node coordinates, face nodes, reference spaces, etc.
#
# ```@docs; canonical=false
# create_mesh
# ```
#
# In the following example, we generate and visualize a mesh of three first order triangles. Only 2-faces are present in this example. The arrays for vertices and faces are empty.
#

#Node coordinates
node_coordinates = SVector{2,Float64}[(0,0),(1,0),(0,1),(1,1),(2,0)]

#Face nodes
face_nodes_0 = Vector{Int}[]
face_nodes_1 = Vector{Int}[]
face_nodes_2 = [[1,2,3],[2,3,4],[2,4,5]]
face_nodes = [
    face_nodes_0,
    face_nodes_1,
    face_nodes_2]

#Reference spaces
reference_spaces_0 = ()
reference_spaces_1 = ()
order = 1
triangle = GT.unit_simplex(Val(2))
triangle3 = GT.lagrange_space(triangle,order)
reference_spaces_2 = (triangle3,)
reference_spaces = (
    reference_spaces_0,
    reference_spaces_1,
    reference_spaces_2)

#Create mesh
mesh = GT.create_mesh(;
    node_coordinates,
    face_nodes,
    reference_spaces)

#Visualize
axis = (;aspect=Makie.DataAspect())
GT.makie_surface(mesh;axis)
GT.makie_lines!(mesh;color=:black)
FileIO.save(joinpath(@__DIR__,"fig_meshes_1.png"),Makie.current_figure()) # hide
nothing # hide

# ![](fig_meshes_1.png)

# In this other slightly more complex example, we define a mesh including
# faces of different dimensions: linear triangles, edges and vertices. To be able to see
# all faces in the visualization, we need to "shrink" them. Otherwise, the triangles would hide the edges and vertices.

#Face nodes
face_nodes_0 = [[1],[3]]
face_nodes_1 = [[1,2],[2,5],[5,4]]
face_nodes_2 = [[1,2,3],[2,3,4],[2,4,5]]
face_nodes = [
    face_nodes_0,
    face_nodes_1,
    face_nodes_2]

#Reference spaces
vertex = GT.unit_simplex(Val(0))
vertex1 = GT.lagrange_space(vertex,order)
segment = GT.unit_simplex(Val(1))
segment2 = GT.lagrange_space(segment,order)
reference_spaces_0 = (vertex1,)
reference_spaces_1 = (segment2,)
reference_spaces_2 = (triangle3,)
reference_spaces = (
    reference_spaces_0,
    reference_spaces_1,
    reference_spaces_2)

#Create mesh
mesh = GT.create_mesh(;
    node_coordinates,
    face_nodes,
    reference_spaces)

#Visualize
axis = (;aspect=Makie.DataAspect())
shrink = 0.8
GT.makie_surface(mesh;axis,shrink)
GT.makie_lines!(mesh;dim=1,shrink)
GT.makie_points!(mesh;dim=0)
FileIO.save(joinpath(@__DIR__,"fig_meshes_2.png"),Makie.current_figure()) # hide
nothing # hide


# ![](fig_meshes_2.png)


#
#
#
#
# ## Domains
#
# A domain in GalerkinToolkit is a subset of $\mathbb{R}^d$, where $d$ is the number of space dimensions. A domain has enough information perform two key operations: numerical integration, and definition of interpolation spaces. All types that implement the *domain interface* are sub-types of [`AbstractDomain`](@ref).
#
# ### Single-face domains
#
# The simplest domain types are defined as a single primitive geometrical shape. One can build such domains with the functions [`unit_n_cube`](@ref) and [`unit_simplex`](@ref). Both function take `N` or `Val(N)`, being `N` an integer equal to the number of spatial dimensions of the object we want to create.

# For example:
cube = GT.unit_n_cube(3)
nothing # hide

# and

triangle = GT.unit_simplex(Val(2))
nothing # hide

# create a unit cube and an unit triangle respectively.

#
# !!! note
#     Many functions in GalerkinToolkit accept `N` or `Val(N)`, when expecting `N` be the number of spatial dimensions. The latter case is needed for type-stability since `N` will be stored as a type parameter in many cases. The former is provided for convenience.
#
#
# Functions [`unit_n_cube`](@ref) and [`unit_simplex`](@ref) build unit hypercubes and unit simplices of any arbitrary dimension. These are the two primitive domains needed in GalerkinToolkit. The rest of the domains are obtained by transforming these primitive ones in multiple ways.

# We can visualize the domains generated above using the Makie recipes provided by GalerkinToolkit. For the cube:

axis = (;aspect=:data)
GT.makie_surface(cube;axis)
GT.makie_lines!(cube;color=:black)
FileIO.save(joinpath(@__DIR__,"fig_geom_1.png"),Makie.current_figure()) # hide
nothing # hide
#
# ![](fig_geom_1.png)

# and for the triangle:
#
axis = (;aspect=Makie.DataAspect())
GT.makie_surface(triangle;axis)
GT.makie_lines!(triangle;color=:black)
FileIO.save(joinpath(@__DIR__,"fig_geom_2.png"),Makie.current_figure()) # hide
nothing # hide
#
# ![](fig_geom_2.png)
#
# The full details about visualization are given in the Visualization section.
#
#
# ### Multi-face domains
# 
# More complex domains are defined as the union of (a subset of) the faces in a computational mesh. Each one of the faces in a mesh is in turn defined as the mapping of a single-face domain (a unit hypercupe or unit simplex in practice). We explain here how to create multi-face domains from a given mesh (assuming that we already have a mesh object).  The details on how to build mesh objects are given in the [Meshes](@ref) section.
#
# ### Selecting all faces of a given dimension
# 
#
# To build a domain from a mesh one needs to define which faces from within the mesh to use. This is done by specifying the dimension of the faces (0, 1, 2, or 3) plus a vector of strings containing the names of the face groups (called *physical faces*) to use. These face groups are included in the mesh object and can be defined as explained in section [Meshes](@ref). This action is performed with function [`domain`](@ref). 
#
# For instance, the following build a computational domain taking all faces of dimension 3 in the computational mesh.

assets_dir = normpath(joinpath(@__DIR__,"..","..","..","assets"))
msh_file = joinpath(assets_dir,"model.msh")
mesh = GT.mesh_from_msh(msh_file)
Ω = GT.domain(mesh,Val(3))

# The domain can be visualized as

axis = (;aspect=:data)
GT.makie_surface(Ω;axis)
GT.makie_lines!(Ω;color=:black)
FileIO.save(joinpath(@__DIR__,"fig_geom_3.png"),Makie.current_figure()) # hide
nothing # hide

# ![](fig_geom_3.png)

# In this other example, we take the faces of dimension 2 including in the face groups named `"circle"`, `"triangle"`, and `"square"`.
physical_names = ["circle", "triangle", "square"]
Γ = GT.domain(mesh,Val(2);physical_names)

#
# We also visualize this domain:

GT.makie_surface(Γ;axis)
GT.makie_lines!(Γ;color=:black)
FileIO.save(joinpath(@__DIR__,"fig_geom_4.png"),Makie.current_figure()) # hide
nothing # hide

# ![](fig_geom_4.png)

#
# ### Using computational domains
#
#

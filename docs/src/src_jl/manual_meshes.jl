#
# ```@meta
# CurrentModule=GalerkinToolkit
# ```
#
# # Meshes
#
#
#
# We use the following dependencies in code snippets in this page.

import GalerkinToolkit as GT
import GLMakie
import Makie
using StaticArrays
import FileIO # hide

#
# ## Mesh specification
#
# A mesh object in GalerkinToolkit contains all geometrical information needed in a finite element (FE) computation.
# A mesh is a set of polygons (or polytopes in general) which we refer to as *faces*. Meshes also
# include additional metadata, including *face groups* used to identify particular faces in the mesh, e.g., to impose boundary conditions.
# All types implementing the meshes are subtypes of [`AbstractMesh`](@ref).
#
# Important features of a mesh include:
#
# - A mesh can potentially contain faces with different number of dimensions. I.e., the same mesh object can include vertices, edges, surfaces, and volumes. The number of dimensions of a mesh is the maximum number of dimension of its faces.
# - The number of dimensions of a mesh can be smaller or equal to the number of *ambient* dimension. The latter is the number of components in a node coordinate vector.
# - A mesh might or might not represent a cell complex (all possible low dimensional faces are present in the mesh). However, many algorithms require to work with a cell complex.
# - A *physical* face $F=\varphi(\hat F)$ in the mesh is defined by transforming a *reference* face $\hat F$ with a mapping $\varphi: \hat F \rightarrow \mathbb{R}^D$, where D is the number of ambient dimensions.  The mapping is defined as $\varphi(\hat x) = \sum_i \hat s_i(\hat x) x_{(F,i)}$. Function $\hat s_i: \hat F \rightarrow \mathbb{R}$ is the scalar basis function number $i$ in the reference (interpolation) space of $F$. The vector $x_{(F,i)}$ contains the coordinates of the local node $i$ in face $F$.
#
#
# Arbitrary mesh objects are defined from low-level quantities with function [`create_mesh`](@ref).
#
# ```@docs; canonical=false
# create_mesh
# ```
#
# The reference interpolation spaces are defined with functions like [`lagrange_space`](@ref). The shape 
# function $\hat s_i$ defining the mapping $\varphi$ can be recovered programmatically from the keyword arguments above using function [`shape_functions`](@ref) as `GT.shape_functions(ref_space)[i]` with `ref_space=reference_spaces[d+1][r]` and `r=face_reference_id[d+1][f]` for the face number `f` in dimension `d`. The vector $x_{(F,i)}$ with the coordinates of the local node $i$ in face $F$ is accessed in the code as `node_coordinates[n]` with `n=face_nodes[d+1][f][i]` being `f` the face number of $F$.
#
# ### Example
#
# In the following example, we generate and visualize a mesh of three first order triangles. Only faces of dimension 2 are present in this example. The arrays for vertices and edges are empty.
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
#
# ### Example
#
# In this other slightly more complex example, we define a mesh including
# faces of different dimensions: surfaces, edges and vertices. To be able to see
# all faces in the visualization, we need to "shrink" them. Otherwise, the surfaces would hide the edges and vertices.

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
# ## Cell complexes
#
# The mesh object described above is general enough to describe cell or [polyhedral complexes](https://en.wikipedia.org/wiki/Polyhedral_complex), but it is not guaranteed that the mesh is indeed a cell complex. For instance, none of the meshes in the two last examples is a cell complex. The first one has no vertices nor edges. The second one has only few vertices and edges, but many are missing.
# One can complete a given mesh with all low-dimensional faces needed to be a cell complex with function [`complexify`](@ref). Calling [`is_cell_complex`](@ref) on the returned mesh, will give `true`.
#
# For function `complexify` to work, neighboring faces should share node ids. I.e., duplicated nodes are not allowed in the input mesh. Otherwise, duplicated faces might be generated, or the functions might not work at all.
#
# ### Example
#
# Let us complete the mesh we generated in the last example into a cell complex.

#Convert
mesh2 = GT.complexify(mesh)
@assert GT.is_cell_complex(mesh2)

#Visualize
GT.makie_surface(mesh2;axis,shrink)
GT.makie_lines!(mesh2;dim=1,shrink)
GT.makie_points!(mesh2;dim=0)
FileIO.save(joinpath(@__DIR__,"fig_meshes_3.png"),Makie.current_figure()) # hide
nothing # hide

# ![](fig_meshes_3.png)
#
# Note that the mesh contains now all low-dimensional faces.


# ## Mesh topology
#
#
# When a mesh is a cell complex, there is a well defined partial order of their faces, or face incidence.
# All face incidence relations are stored in an object called mesh *topology*. A mesh topology is obtained with
# function [`topology`](@ref) called on a given mesh object: `topo = GT.topology(mesh)`. The mesh needs to be a cell complex for this to work.  Then, one uses function [`face_incidence`](@ref) on the topology object to get the incidence relations. 
# `GT.face_incidence(topo,D,d)` is a long vector of small vectors of integers, often implemented with a `JaggedArray`.
# - For `d<D`, `GT.face_incidence(topo,D,d)[F]` is a vector of integers containing the ids of the faces of dimension `d` on the boundary of face number `F` of dimension `D`.
# - For `d<D`, `GT.face_incidence(topo,d,D)[f]` is a vector of integers containing the ids of the faces of dimension `D` around the face number `f` of dimension `d`.
# - For `d==D`, `GT.face_incidence(topo,D,D)[F] == [F]`.
#
# ### Example
#
# Let us get some of the incidence relations for the cell complex we generated above.

topo = GT.topology(mesh2)
surface_to_edges = GT.face_incidence(topo,2,1)

# This cutout is read as follows. Surface 1 has edges 1, 4, and 5 on its boundary; surface 2 has edges 5, 6, 7 on its boundary; etc.

edge_to_surfaces = GT.face_incidence(topo,1,2)

# According to this output, edge 1 touches surface 1, edge 2 touches surface 3, etc. We can also see that there are two interior edges touching two surfaces. Edge 5 touches surfaces 1 and 2, and edge 6 touches surfaces 2 and 3.
#
# ## Face groups
#
# Face groups allow us to select specific faces in a mesh for different modeling purposes: impose boundary conditions, define different equations in different parts of the mesh etc. A face group is a vector of integers containing the ids of the faces in this group plus a string containing a name for this group. This groups are stored using a dictionary that maps strings to vectors (group names to group definitions) in a per dimension basis (one dictionary per dimension). The vector contains faces of the same dimension, but it is possible define groups containing faces of different dimensions by splitting them in a vector per dimension.
# Face groups can overlap and can be added after the mesh object is created.
#
# Face groups are accessed and added using function [`group_faces`](@ref).
#
# ```@docs; canonical=false
# group_faces
# ```
#
# ### Common face groups
#
# GalerkinToolkit provides a number of functions that generate commonly used face groups.
# - [`label_boundary_faces!`](@ref)
# - [`label_interior_faces!`](@ref)
# - [`label_faces_in_dim!`](@ref)
#
# These function do what the name suggests (see the docstrings for further details). The first one is often use to impose boundary
# conditions and the second one in discontinuous Galerkin methods to define interior penalty terms. They are often called under the hood when calling functions like [`boundary`](@ref) and [`skeleton`](@ref).
#
#
#
#
# ### Example
#
# Let us add some face groups to the last mesh we created.



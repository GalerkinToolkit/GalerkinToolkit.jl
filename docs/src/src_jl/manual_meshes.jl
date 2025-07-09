#
# ```@meta
# CurrentModule=GalerkinToolkit
# ```
#
# # Meshes
#
#
# A mesh object in GalerkinToolkit contains all geometrical information needed in a FEM computation.
# A mesh is a set of polygons (or polytopes in general) which we refer to as *faces*. Meshes also
# include additional meta data, including *face labels* used to identify particular faces in the mesh, e.g., to impose boundary conditions.
# All types implementing the meshes are subtypes of [`AbstractMesh`](@ref).
#
# Important features of a mesh include:
#
# - A mesh can potentially contain faces with different number of dimensions. I.e., the same mesh object can include vertices, edges, surfaces, and volumes.
#
# - If the mesh contains faces of the same dimension, it is called a *chain*.
#
# - A mesh might or might not represent a cell complex (all possible low dimensional faces are present in the mesh). However, many algorithms require to work with a cell complex.
#
#
# We use the following dependencies in code snippets in this page.

import GalerkinToolkit as GT
import GLMakie
import Makie
using StaticArrays
import FileIO # hide

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



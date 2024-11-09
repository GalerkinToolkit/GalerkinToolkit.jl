# # Mesh generation
#


import GalerkinToolkit as GT
import PartitionedArrays as PA
import GLMakie as Makie
import Gmsh
import Metis
import FileIO # hide

# ## Cartesian meshes
#
# Generate a Cartesian mesh of the 3D domain $(0,1)\times(2,3)\times(-1,1)$ using 5 cells per direction and visualize it.

domain = (0,1,1,3,-1,1)
cells = (5,5,5)
mesh = GT.cartesian_mesh(domain,cells)
Makie.plot(mesh,color=:pink,strokecolor=:blue)
FileIO.save(joinpath(@__DIR__,"fig_mg_cm_3d.png"),Makie.current_figure()) # hide

# ![](fig_mg_cm_3d.png)


# Generate a Cartesian mesh of the 2D domain $(0,1)\times(2,3)$ and visualize it.

domain = (0,1,2,3)
cells = (5,5)
mesh = GT.cartesian_mesh(domain,cells)
Makie.plot(mesh,color=:pink,strokecolor=:blue)
FileIO.save(joinpath(@__DIR__,"fig_mg_cm_2d.png"),Makie.current_figure()) # hide

# ![](fig_mg_cm_2d.png)


# Now visualize all objects (vertices, edges, faces) in the mesh.

Makie.plot(mesh,color=:pink,strokecolor=:blue,shrink=0.6,dim=(0:2))
FileIO.save(joinpath(@__DIR__,"fig_mg_cm_2d_a.png"),Makie.current_figure()) # hide

# ![](fig_mg_cm_2d_a.png)


# Now, do not generate low-dimensional objects on the interior of the mesh

domain = (0,1,2,3)
cells = (5,5)
mesh = GT.cartesian_mesh(domain,cells;complexify=false)
Makie.plot(mesh,color=:pink,strokecolor=:blue,shrink=0.6,dim=(0:2))
FileIO.save(joinpath(@__DIR__,"fig_mg_cm_2d_b.png"),Makie.current_figure()) # hide

# ![](fig_mg_cm_2d_b.png)


# !!! note
#     Most algorithms require working with a polytopal complex (i.e., a mesh containing all low dimensional objects). Thus using the option `complexify=false` is not recommended, unless you know what you are doing.
#

# Now, use triangles instead of squares.

domain = (0,1,2,3)
cells = (5,5)
mesh = GT.cartesian_mesh(domain,cells;simplexify=true)
Makie.plot(mesh,color=:pink,strokecolor=:blue)
FileIO.save(joinpath(@__DIR__,"fig_mg_cm_2d_c.png"),Makie.current_figure()) # hide

# ![](fig_mg_cm_2d_c.png)

# ## Gmsh meshes

# Read a mesh from a ".msh" file and visualize it.

repodir = joinpath(@__DIR__,"..","..","..")
fn = joinpath(repodir,"assets","mesh1.msh")
mesh = GT.mesh_from_gmsh(fn)
Makie.plot(mesh,color=:pink,strokecolor=:blue)
FileIO.save(joinpath(@__DIR__,"fig_mg_gmsh.png"),Makie.current_figure()) # hide

# ![](fig_mg_gmsh.png)

# Now, define the mesh using GMSH DSL

# !!! warning
#     TODO
#     
#

# Now, define the mesh using GMSH julia API

# !!! warning
#     TODO
#     

# ## Parallel meshes
#

# Generate a Cartesian mesh of the 2D domain $(0,1)\times(2,3)$ using 10 cells in each direction.
# Partitioned it into 2 parts per direction and visualize with face color according to part owner.

domain = (0,1,2,3)
cells_per_dir = (10,10)
parts_per_dir = (2,2)
parts = LinearIndices((prod(parts_per_dir),))
pmesh = GT.cartesian_mesh(domain,cells_per_dir;parts_per_dir,parts)
Makie.plot(pmesh;color=GT.FaceData("__OWNER__"),strokecolor=:blue)
FileIO.save(joinpath(@__DIR__,"fig_mg_pmesh.png"),Makie.current_figure()) # hide

# ![](fig_mg_pmesh.png)


# !!! warning
#     TODO Maybe a different API for parts_per_dir and parts?
#     TODO Maybe a different API for GT.FaceData("__OWNER__")
#     


# Generate a mesh on a single machine (using Gmsh in this case), partition it using Metis into 4 parts, and visualize it.

np = 4
parts = LinearIndices((np,))
pmesh = PA.map_main(parts) do parts
    fn = joinpath(repodir,"assets","mesh1.msh")
    mesh = GT.mesh_from_gmsh(fn)
    graph = GT.mesh_graph(mesh)
    graph_partition = Metis.partition(graph,np)
    GT.partition_mesh(mesh,np;graph,graph_partition)
end |> GT.scatter_mesh
Makie.plot(pmesh,color=GT.FaceData("__OWNER__"),strokecolor=:blue)
FileIO.save(joinpath(@__DIR__,"fig_mg_pmesh_a.png"),Makie.current_figure()) # hide

# ![](fig_mg_pmesh_a.png)

# !!! warning
#     * TODO `color=GT.FaceData("__PART__")` not working
#     * TODO better syntax for `color=GT.FaceData("__PART__")` ?
#     

# ## Meshes from arrays
#

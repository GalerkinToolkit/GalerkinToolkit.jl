# # Mesh generation
#

repodir = joinpath(@__DIR__,"..","..","..")

import GalerkinToolkit as GT
import PartitionedArrays as PA
import GLMakie
import Gmsh
import Metis

# ## Cartesian meshes
#
# Generate a Cartesian mesh of the 3D domain $(0,1)\times(2,3)\times(-1,1)$ using 10 cells per direction and visualize it.

function ex()
    domain = (0,1,1,3,-1,1)
    cells = (10,10,10)
    mesh = GT.cartesian_mesh(domain,cells)
    GLMakie.plot(mesh,color=:pink,strokecolor=:blue)
end

 #ex()

# !!! warning
#     TODO: Visualization only working for simplices at the moment

# Generate a Cartesian mesh of the 2D domain $(0,1)\times(2,3)$ and visualize it.


function ex()
    domain = (0,1,2,3)
    cells = (10,10)
    mesh = GT.cartesian_mesh(domain,cells)
    GLMakie.plot(mesh,color=:pink,strokecolor=:blue)
end

 #ex()


# Now visualize all objects (vertices, edges, faces) in the mesh.

function ex()
    domain = (0,1,2,3)
    cells = (10,10)
    mesh = GT.cartesian_mesh(domain,cells)
    GLMakie.plot(mesh,color=:pink,shrink=0.6,dim=(0:2))
end

 #ex()


# Now, do not generate low-dimensional objects on the interior of the mesh

function ex()
    domain = (0,1,2,3)
    cells = (10,10)
    mesh = GT.cartesian_mesh(domain,cells;complexify=false)
    GLMakie.plot(mesh,color=:pink,shrink=0.6,dim=(0:2))
end

 #ex()

# !!! warning
#     TODO: Error if dim=0:2
#

# Now, use triangles instead of squares.

function ex()
    domain = (0,1,2,3)
    cells = (10,10)
    mesh = GT.cartesian_mesh(domain,cells;simplexify=true)
    GLMakie.plot(mesh,color=:pink,strokecolor=:blue)
end

ex()


# ## Gmsh meshes

# Read a mesh from a ".msh" file and visualize it.


function ex()
    fn = joinpath(repodir,"assets","mesh1.msh")
    mesh = GT.mesh_from_gmsh(fn)
    GLMakie.plot(mesh,color=:pink,strokecolor=:blue)
end

ex()

# Now, define the mesh using GMSH DSL

function ex()
    code = """
    SetFactory(\"OpenCASCADE\");
    Box(1) = {0, 0, 0, 1, 1, 1};
    Cylinder(2) = {0, 0.5, 0.5, 1, 0, 0, 0.4, 2*Pi};
    Cylinder(3) = {0.5, 0.0, 0.5, 0, 1, 0, 0.4, 2*Pi};
    BooleanUnion{ Volume{3}; Delete; }{ Volume{2}; Delete; }
    BooleanDifference{ Volume{1}; Delete; }{ Volume{2}; Delete; }
    Physical Surface(\"top\") = {14};
    Physical Surface(\"bottom\") = {15};
    Physical Volume(\"volume\") = {1};
    """
    #TODO using GMSH.jl to generate the mesh
    mesh = GT.mesh_from_gmsh_module()
    GLMakie.plot(mesh,color=:pink,strokecolor=:blue)
end

 # ex()

# !!! warning
#     TODO using GMSH.jl to generate the mesh
#     
#

# Now, define the mesh using GMSH julia API

function ex()
    #TODO use julia API to generate the mesh above
    #TODO using GMSH.jl to generate the mesh
    mesh = GT.mesh_from_gmsh_module()
    GLMakie.plot(mesh,color=:pink,strokecolor=:blue)
end

 # ex()

# !!! warning
#     TODO using GMSH.jl to generate the mesh
#     

# ## Parallel meshes
#

# Generate a Cartesian mesh of the 2D domain $(0,1)\times(2,3)$ using 10 cells in each direction.
# Partitioned it into 2 parts per direction and visualize with face color according to part owner.

function ex()
    domain = (0,1,2,3)
    cells_per_dir = (10,10)
    parts_per_dir = (2,2)
    parts = LinearIndices((prod(parts_per_dir),))
    mesh = GT.cartesian_mesh(domain,cells_per_dir;parts_per_dir,parts)
    GLMakie.plot(mesh,color=GT.FaceData("__OWNER__"),strokecolor=:blue)
end

 #ex()


# !!! warning
#     TODO Maybe a different API for parts_per_dir and parts?
#     TODO Maybe a different API for GT.FaceData("__OWNER__")
#     


# Now, by defining a partition strategy.

function ex()
    domain = (0,1,2,3)
    cells_per_dir = (10,10)
    parts_per_dir = (2,2)
    parts = LinearIndices((prod(parts_per_dir),))
    partition_strategy = GT.partition_strategy(
        graph_nodes=:cells,
        graph_edges=:nodes,
        ghost_layers=1)
    mesh = GT.cartesian_mesh(domain,cells_per_dir;
       parts_per_dir,
       parts,
       partition_strategy)
    GLMakie.plot(mesh,color=GT.FaceData("__OWNER__"),strokecolor=:blue)
end

 #ex()

# Generate a mesh on a single machine (using Gmsh in this case), partition it using Metis into 4 parts, and visualize it.

function ex()
    np = 4
    parts = LinearIndices((np,))
    pmesh = PA.map_main(parts) do parts
        fn = joinpath(repodir,"assets","mesh1.msh")
        mesh = GT.mesh_from_gmsh(fn)
        graph = GT.mesh_graph(mesh)
        graph_partition = Metis.partition(graph,np)
        GT.partition_mesh(mesh,np;graph,graph_partition)
    end |> GT.scatter_mesh
    #GLMakie.plot(pmesh,color=GT.FaceData("__OWNER__"),strokecolor=:blue)
    plt = GT.plot(pmesh)
    GLMakie.plot(plt,color=GT.FaceData("__OWNER__"),strokecolor=:blue)
end

ex()

# !!! warning
#     * TODO plot recipie for PMesh missing
#     * TODO `color=GT.FaceData("__PART__")` not working
#     * TODO better syntax for `color=GT.FaceData("__PART__")` ?
#     


# Now, by defining a partition strategy.

function ex()
    np = 4
    parts = LinearIndices((np,))
    pmesh = PA.map_main(parts) do parts
        fn = joinpath(repodir,"assets","mesh1.msh")
        mesh = GT.mesh_from_gmsh(fn)
        partition_strategy = GT.partition_strategy(
             graph_nodes=:cells,
             graph_edges=:nodes,
             ghost_layers=0)
        graph = GT.mesh_graph(mesh;partition_strategy)
        graph_partition = Metis.partition(graph,np)
        GT.partition_mesh(mesh,np;graph,graph_partition,partition_strategy)
    end |> GT.scatter_mesh
    #GLMakie.plot(pmesh,color=GT.FaceData("__OWNER__"),strokecolor=:blue)
    plt = GT.plot(pmesh)
    GLMakie.plot(plt,color=GT.FaceData("__OWNER__"),strokecolor=:blue)
end

ex()

# !!! warning
#     TODO Why colors are not right?
#     


# ## Mesh transformations
#
# !!! warning
#     TODO
#     * Simplexify
#     * Complexify
#     * Refine uniformly
#     

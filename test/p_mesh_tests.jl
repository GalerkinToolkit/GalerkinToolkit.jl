module PMeshTests

import GalerkinToolkit as GT
using Test
using PartitionedArrays
using Metis

np = 2
parts = DebugArray(LinearIndices((np,)))
mesh = GT.cartesian_mesh((0,1,0,1),(2,2))
cell_to_color = [1,1,2,2]
pmesh = GT.partition_mesh(mesh,np;partition_strategy=GT.partition_strategy(graph_nodes=:cells,graph_edges=:nodes,ghost_layers=0),parts,graph_partition=cell_to_color)
pmesh = GT.partition_mesh(mesh,np;partition_strategy=GT.partition_strategy(graph_nodes=:cells,graph_edges=:nodes,ghost_layers=1),parts,graph_partition=cell_to_color)

node_to_color = [1,1,1,1,2,2,2,2,2]
pmesh = GT.partition_mesh(mesh,np;partition_strategy=GT.partition_strategy(graph_nodes=:nodes,graph_edges=:cells,ghost_layers=0),parts,graph_partition=node_to_color)
pmesh = GT.partition_mesh(mesh,np;partition_strategy=GT.partition_strategy(graph_nodes=:nodes,graph_edges=:cells,ghost_layers=1),parts,graph_partition=node_to_color)

## In this case, we do the hard work only on the main proc
msh =  joinpath(@__DIR__,"..","assets","demo.msh")
np = 4
parts = DebugArray(LinearIndices((np,)))
pmesh = map_main(parts) do parts
    mesh = GT.mesh_from_msh(msh)
    partition_strategy = GT.partition_strategy(graph_nodes=:nodes,graph_edges=:cells)
    graph = GT.mesh_graph(mesh;partition_strategy)
    graph_partition = Metis.partition(graph,np)
    GT.partition_mesh(mesh,np;partition_strategy,graph,graph_partition)
end |> GT.scatter_mesh

# In this one, we do only the graph partition on the main
# but we load the mesh everywhere
mesh = GT.mesh_from_msh(msh)
partition_strategy = GT.partition_strategy(graph_nodes=:nodes,graph_edges=:cells,ghost_layers=1)
graph = GT.mesh_graph(mesh;partition_strategy)
graph_partition = map_main(parts) do parts
     Metis.partition(graph,np)
end |> multicast |> PartitionedArrays.getany
pmesh = GT.partition_mesh(mesh,np;partition_strategy,parts,graph,graph_partition)

graph_nodes = :cells
partition_strategy = GT.partition_strategy(graph_nodes=:cells,graph_edges=:nodes)
graph = GT.mesh_graph(mesh;partition_strategy)
graph_partition = map_main(parts) do parts
     Metis.partition(graph,np)
end |> multicast |> PartitionedArrays.getany
pmesh = GT.partition_mesh(mesh,np;partition_strategy,parts,graph,graph_partition)
pmesh = GT.partition_mesh(mesh,np;partition_strategy,parts,graph,graph_partition)

domain = (0,1,0,1)
cells_per_dir = (4,4)
parts_per_dir = (2,2)
pmesh = GT.cartesian_pmesh(domain,cells_per_dir,parts,parts_per_dir;partition_strategy=GT.partition_strategy(graph_nodes=:cells,graph_edges=:nodes,ghost_layers=1))
pmesh = GT.cartesian_pmesh(domain,cells_per_dir,parts,parts_per_dir;partition_strategy=GT.partition_strategy(graph_nodes=:cells,graph_edges=:nodes,ghost_layers=0))
np = prod(parts_per_dir)
parts = DebugArray(LinearIndices((np,)))
pmesh = GT.cartesian_pmesh(domain,cells_per_dir,parts,parts_per_dir;partition_strategy=GT.partition_strategy(graph_nodes=:cells,graph_edges=:nodes,ghost_layers=1))
pmesh = GT.cartesian_pmesh(domain,cells_per_dir,parts,parts_per_dir;partition_strategy=GT.partition_strategy(graph_nodes=:cells,graph_edges=:nodes,ghost_layers=0))
@test GT.partition_strategy(pmesh).ghost_layers == 0
@test GT.partition_strategy(pmesh).graph_nodes === :cells

end # module

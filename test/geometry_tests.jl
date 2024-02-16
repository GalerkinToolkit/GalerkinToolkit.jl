module MeshInterfaceTests

using Test
import GalerkinToolkit as gk
using WriteVTK
using PartitionedArrays
using Metis

spx0 = gk.unit_simplex(0)
spx1 = gk.unit_simplex(1)
spx2 = gk.unit_simplex(2)
spx3 = gk.unit_simplex(3)
display(spx3)

cube0 = gk.unit_n_cube(0)
cube1 = gk.unit_n_cube(1)
cube2 = gk.unit_n_cube(2)
cube3 = gk.unit_n_cube(3)
display(cube3)

@show typeof(spx0)
@show typeof(cube0)
@show typeof(spx1)
@show typeof(cube1)

degree = 4
quad = gk.default_quadrature(spx0,degree)
quad = gk.default_quadrature(spx1,degree)
quad = gk.default_quadrature(spx2,degree)
quad = gk.default_quadrature(spx3,degree)

quad = gk.default_quadrature(cube0,degree)
quad = gk.default_quadrature(cube1,degree)
quad = gk.default_quadrature(cube2,degree)
quad = gk.default_quadrature(cube3,degree)

order = 1
fe = gk.lagrange_mesh_face(spx0,order)
fe = gk.lagrange_mesh_face(spx1,order)
fe = gk.lagrange_mesh_face(spx2,order)
fe = gk.lagrange_mesh_face(spx3,order)
display(fe)

fe = gk.lagrange_mesh_face(cube0,order)
fe = gk.lagrange_mesh_face(cube1,order)
fe = gk.lagrange_mesh_face(cube2,order)
fe = gk.lagrange_mesh_face(cube3,order)
display(fe)

fe = gk.lagrange_mesh_face(cube0,order)
@show gk.monomial_exponents(fe)
@show gk.node_coordinates(fe)
fe = gk.lagrange_mesh_face(cube2,order)
@show gk.node_coordinates(fe)

spx2 = gk.unit_simplex(2)
quad = gk.default_quadrature(spx2,degree)
fe = gk.lagrange_mesh_face(spx2,order)
funs = gk.shape_functions(fe)
x = gk.coordinates(quad)
B = broadcast(gk.value,permutedims(funs),x)
display(B)
tabulator = gk.tabulator(fe)
A = tabulator(gk.value,x)
@test A≈B
x = gk.node_coordinates(fe)
A = tabulator(gk.value,x)

outdir = mkpath(joinpath(@__DIR__,"..","output"))

msh =  joinpath(@__DIR__,"..","assets","demo.msh")
mesh = gk.mesh_from_gmsh(msh;complexify=false)
vtk_grid(joinpath(outdir,"demo"),gk.vtk_args(mesh)...) do vtk
    gk.vtk_physical_faces!(vtk,mesh)
    gk.vtk_physical_nodes!(vtk,mesh)
end
for d in 0:gk.num_dims(mesh)
    vtk_grid(joinpath(outdir,"demo_$d"),gk.vtk_args(mesh,d)...) do vtk
        gk.vtk_physical_faces!(vtk,mesh,d)
        gk.vtk_physical_nodes!(vtk,mesh,d)
    end
end

@show gk.unit_simplex(0) |> gk.boundary
@show gk.unit_simplex(1) |> gk.boundary
@show gk.unit_simplex(2) |> gk.boundary
@show gk.unit_simplex(3) |> gk.boundary

@show gk.unit_n_cube(0) |> gk.boundary
@show gk.unit_n_cube(1) |> gk.boundary
@show gk.unit_n_cube(2) |> gk.boundary
@show gk.unit_n_cube(3) |> gk.boundary

gk.unit_n_cube(0) |> gk.opposite_faces
gk.unit_n_cube(1) |> gk.opposite_faces
gk.unit_n_cube(2) |> gk.opposite_faces
gk.unit_n_cube(3) |> gk.opposite_faces

gk.opposite_faces(cube3,0)
gk.opposite_faces(cube3,1)
gk.opposite_faces(cube3,2)
gk.opposite_faces(cube3,3)

order = 2
gk.lagrange_mesh_face(spx1,order) |> gk.boundary |> gk.topology

mesh = gk.mesh_from_gmsh(msh;complexify=false)

new_mesh, old_to_new = gk.complexify(mesh)

mesh = gk.mesh_from_gmsh(msh)

domain = (0,1,0,1)
cells = (2,2)
mesh = gk.cartesian_mesh(domain,cells)
mesh = gk.cartesian_mesh(domain,cells,boundary=false)
mesh = gk.cartesian_mesh(domain,cells,simplexify=true)
mesh = gk.cartesian_mesh(domain,cells,boundary=false,simplexify=true)

mesh = gk.cartesian_mesh(domain,cells)
vtk_grid(joinpath(outdir,"cartesian"),gk.vtk_args(mesh)...) do vtk
    gk.vtk_physical_faces!(vtk,mesh)
    gk.vtk_physical_nodes!(vtk,mesh)
end
for d in 0:gk.num_dims(mesh)
    vtk_grid(joinpath(outdir,"cartesian_$d"),gk.vtk_args(mesh,d)...) do vtk
        gk.vtk_physical_faces!(vtk,mesh,d)
        gk.vtk_physical_nodes!(vtk,mesh,d)
    end
end

mesh = gk.mesh_from_gmsh(msh)
face_groups = gk.physical_faces(mesh)
group_names = gk.physical_names(mesh,2)
group_names = gk.physical_names(mesh)
group_names = gk.physical_names(mesh;merge_dims=true)
node_groups = gk.physical_nodes(mesh;merge_dims=true)
node_groups = gk.physical_nodes(mesh;merge_dims=true,disjoint=true)

vmesh, vglue = gk.visualization_mesh(mesh)

mesh = gk.cartesian_mesh(domain,cells)
gk.mesh_graph(mesh;graph_nodes=:nodes)
gk.mesh_graph(mesh;graph_nodes=:cells)
gk.mesh_graph(mesh;graph_nodes=:nodes,graph_edges=:cells)
gk.mesh_graph(mesh;graph_nodes=:nodes,graph_edges=:faces)
gk.mesh_graph(mesh;graph_nodes=:cells,graph_edges=:nodes)
gk.mesh_graph(mesh;graph_nodes=:cells,graph_edges=:faces)

np = 2
parts = DebugArray(LinearIndices((np,)))
mesh = gk.cartesian_mesh((0,1,0,1),(2,2))
cell_to_color = [1,1,2,2]
pmesh = gk.partition_mesh(mesh,np)
pmesh = gk.partition_mesh(mesh,np;parts)
pmesh = gk.partition_mesh(mesh,np;parts,graph_partition=cell_to_color,graph_nodes=:cells,ghost_layers=0)
pmesh = gk.partition_mesh(mesh,np;parts,graph_partition=cell_to_color,graph_nodes=:cells,ghost_layers=1)

node_to_color = [1,1,1,1,2,2,2,2,2]
pmesh = gk.partition_mesh(mesh,np;parts,graph_partition=node_to_color,graph_nodes=:nodes,ghost_layers=0)
pmesh = gk.partition_mesh(mesh,np;parts,graph_partition=node_to_color,graph_nodes=:nodes,ghost_layers=1)

function setup(mesh,ids,rank)
    face_to_owner = zeros(Int,sum(gk.num_faces(mesh)))
    D = gk.num_dims(mesh)
    for d in 0:D
        face_to_owner[gk.face_range(mesh,d)] = local_to_owner(gk.face_indices(ids,d))
    end
    pvtk_grid(joinpath(outdir,"pmesh"),gk.vtk_args(mesh)...;part=rank,nparts=np) do vtk
        gk.vtk_physical_faces!(vtk,mesh)
        gk.vtk_physical_nodes!(vtk,mesh)
        vtk["piece"] = fill(rank,sum(gk.num_faces(mesh)))
        vtk["owner"] = local_to_owner(gk.node_indices(ids))
        vtk["owner"] = face_to_owner
    end
end
map(setup,partition(pmesh),gk.index_partition(pmesh),parts)

## In this case, we do the hard work only on the main proc
np = 4
parts = DebugArray(LinearIndices((np,)))
pmesh = map_main(parts) do parts
    mesh = gk.mesh_from_gmsh(msh)
    graph_nodes = :nodes
    graph = gk.mesh_graph(mesh;graph_nodes=:nodes)
    graph_partition = Metis.partition(graph,np)
    gk.partition_mesh(mesh,np;graph,graph_nodes,graph_partition)
end |> gk.scatter_mesh

function setup(mesh,ids,rank)
    face_to_owner = zeros(Int,sum(gk.num_faces(mesh)))
    D = gk.num_dims(mesh)
    for d in 0:D
        face_to_owner[gk.face_range(mesh,d)] = local_to_owner(gk.face_indices(ids,d))
    end
    pvtk_grid(joinpath(outdir,"pmesh"),gk.vtk_args(mesh)...;part=rank,nparts=np) do vtk
        gk.vtk_physical_faces!(vtk,mesh)
        gk.vtk_physical_nodes!(vtk,mesh)
        vtk["piece"] = fill(rank,sum(gk.num_faces(mesh)))
        vtk["owner"] = local_to_owner(gk.node_indices(ids))
        vtk["owner"] = face_to_owner
    end
end
map(setup,partition(pmesh),gk.index_partition(pmesh),parts)

# In this one, we do only the graph partition on the main
# but we load the mesh everywhere
mesh = gk.mesh_from_gmsh(msh)
graph_nodes = :nodes
graph = gk.mesh_graph(mesh;graph_nodes)
graph_partition = map_main(parts) do parts
     Metis.partition(graph,np)
end |> multicast |> PartitionedArrays.getany
pmesh = gk.partition_mesh(mesh,np;parts,graph,graph_nodes,graph_partition,ghost_layers=0)
pmesh = gk.partition_mesh(mesh,np;parts,graph,graph_nodes,graph_partition,ghost_layers=1)

graph_nodes = :cells
graph = gk.mesh_graph(mesh;graph_nodes)
graph_partition = map_main(parts) do parts
     Metis.partition(graph,np)
end |> multicast |> PartitionedArrays.getany
pmesh = gk.partition_mesh(mesh,np;parts,graph,graph_nodes,graph_partition,ghost_layers=0)
pmesh = gk.partition_mesh(mesh,np;parts,graph,graph_nodes,graph_partition,ghost_layers=1)

# coarse pmesh testing and visualization
domain = (0,1,0,1)
cells_per_dir = (4,4)
parts_per_dir = (2,2)
pmesh = gk.cartesian_mesh(domain,cells_per_dir,parts_per_dir)
pmesh = gk.cartesian_mesh(domain,cells_per_dir,parts_per_dir;ghost_layers=0)
np = prod(parts_per_dir)
parts = DebugArray(LinearIndices((np,)))
pmesh = gk.cartesian_mesh(domain,cells_per_dir,parts_per_dir;parts,ghost_layers=0)

function setup(mesh,ids,rank)
    face_to_owner = zeros(Int,sum(gk.num_faces(mesh)))
    D = gk.num_dims(mesh)
    for d in 0:D
        face_to_owner[gk.face_range(mesh,d)] = local_to_owner(gk.face_indices(ids,d))
    end
    pvtk_grid(joinpath(outdir, "pmesh-cartesian"), gk.vtk_args(mesh)...; part=rank, nparts=np) do vtk
        gk.vtk_physical_faces!(vtk,mesh)
        gk.vtk_physical_nodes!(vtk,mesh)
        vtk["piece"] = fill(rank,sum(gk.num_faces(mesh)))
        vtk["owner"] = local_to_owner(gk.node_indices(ids))
        vtk["owner"] = face_to_owner
    end
end
map(setup,partition(pmesh),gk.index_partition(pmesh),parts)

domain = (0,1,0,1)
cells = (10,10)
fine_mesh = gk.cartesian_mesh(domain,cells)

domain = (0,30,0,10)
cells = (2,2)
coarse_mesh = gk.cartesian_mesh(domain,cells)

final_mesh, glue = gk.two_level_mesh(coarse_mesh,fine_mesh)

vtk_grid(joinpath(outdir,"two-level-mesh"),gk.vtk_args(final_mesh)...) do vtk
    gk.vtk_physical_faces!(vtk,final_mesh)
    gk.vtk_physical_nodes!(vtk,final_mesh)
    vtk["node_ids"] = 1:gk.num_nodes(final_mesh)
end

# Simple two level parallel mesh test and visualization
fine_mesh = gk.cartesian_mesh((0, 1, 0, 1), (2, 2)) # very simple fine mesh
domain = (0,30,0,10)
cells = (4,4)
parts_per_dir = (2,2)
np = prod(parts_per_dir)
parts = DebugArray(LinearIndices((np,)))
coarse_mesh = gk.cartesian_mesh(domain,cells,parts_per_dir; parts, ghost_layers=0)
final_pmesh, final_pglue = gk.two_level_mesh(coarse_mesh,fine_mesh)

function final_pmesh_setup(mesh, ids, rank)
    face_to_owner = zeros(Int, sum(gk.num_faces(mesh)))
    D = gk.num_dims(mesh)
    for d in 0:D
        face_to_owner[gk.face_range(mesh, D)] = local_to_owner(gk.face_indices(ids, D))
    end
    pvtk_grid(
        joinpath(outdir, "final-pmesh-cartesian"), 
        gk.vtk_args(mesh)...; 
        part=rank, nparts=np) do vtk

        gk.vtk_physical_faces!(vtk, mesh)
        gk.vtk_physical_nodes!(vtk, mesh)
        vtk["piece"] = fill(rank, sum(gk.num_faces(mesh)))
        vtk["owner"] = local_to_owner(gk.node_indices(ids))
        vtk["owner"] = face_to_owner
    end
end

map(
    final_pmesh_setup, 
    partition(final_pmesh), 
    gk.index_partition(final_pmesh), 
    parts)

end # module

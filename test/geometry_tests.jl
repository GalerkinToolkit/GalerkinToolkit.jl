module MeshInterfaceTests

using Test
import GalerkinToolkit as GT
using WriteVTK
using PartitionedArrays
using Metis

###################################################################
# Helper functions for visualization
###################################################################
function visualize_mesh(mesh, outpath, glue = nothing, d = nothing)
    node_ids = collect(1:GT.num_nodes(mesh))

    # NOTE: without this check, attempting to visualize a 3D periodic mesh fails because
    # the `periodic_nodes` property is a Vector and not a Pair
    valid_periodic_mesh = GT.has_periodic_nodes(mesh) && GT.periodic_nodes(mesh) isa Pair
    if valid_periodic_mesh
        # Get periodic node info aboutcell
        periodic_nodes = GT.periodic_nodes(mesh)
        fine_pnode_to_fine_node = periodic_nodes.first
        fine_pnode_to_master_fine_node = periodic_nodes.second

        # Labeling periodic nodes
        fine_node_to_master_fine_node = copy(node_ids)
        fine_node_to_master_fine_node[fine_pnode_to_fine_node] =
            fine_pnode_to_master_fine_node

        # Handles periodic vertices whose master node is the master of another periodic
        # vertex assuming one level of indirection for periodicity
        # (e.g., vertex -> master -> master)
        n_fine_nodes = length(node_ids)
        for fine_node = 1:n_fine_nodes
            master_fine_node = fine_node_to_master_fine_node[fine_node]
            if fine_node == master_fine_node
                continue
            end

            master_master_fine_node = fine_node_to_master_fine_node[master_fine_node]
            fine_node_to_master_fine_node[fine_node] = master_master_fine_node
        end
    end

    # Visualize in paraview
    if isnothing(d)
        vtk_grid(outpath, GT.vtk_args(mesh)...) do vtk
            GT.vtk_physical_faces!(vtk, mesh)
            GT.vtk_physical_nodes!(vtk, mesh)
            valid_periodic_mesh &&
                (vtk["periodic_master_id"] = fine_node_to_master_fine_node)
            vtk["node_id"] = node_ids
            if !isnothing(glue)
                vtk["coarse_cell_id"] = glue.final_cell_to_coarse_cell
            end
        end
    else
        vtk_grid(outpath, GT.vtk_args(mesh, d)...) do vtk
            GT.vtk_physical_faces!(vtk, mesh, d)
            GT.vtk_physical_nodes!(vtk, mesh, d)
            valid_periodic_mesh &&
                (vtk["periodic_master_id"] = fine_node_to_master_fine_node)
            vtk["node_id"] = node_ids
            if !isnothing(glue)
                vtk["coarse_cell_id"] = glue.final_cell_to_coarse_cell
            end
        end
    end
end

"""
    visualize_pmesh(pmesh::GT.PMesh, parts, nparts, pmesh_vtk_fpath::String)

Writes a parallel mesh on `nparts` partitions to `pmesh_vtk_fpath`.
"""
function visualize_pmesh(pmesh::GT.PMesh, parts, nparts, pmesh_vtk_fpath::String)
    map(
        visualize_pmesh_setup(nparts, pmesh_vtk_fpath),
        partition(pmesh),
        GT.index_partition(pmesh),
        parts,
    )
end

"""
    visualize_pmesh_setup(nparts, outpath)

Return function that writes vtk for each part in a parallel mesh.
"""
function visualize_pmesh_setup(nparts, outpath)
    @assert isabspath(outpath) "abspath with pvtk_grid ensures function" # PR this?
    function setup(mesh, ids, rank)
        face_to_owner = zeros(Int, sum(GT.num_faces(mesh)))
        D = GT.num_dims(mesh)
        face_to_owner[GT.face_range(mesh, D)] = local_to_owner(GT.face_indices(ids, D))
        pvtk_grid(
            outpath,
            GT.vtk_args(mesh)...;
            part = rank,
            nparts = nparts,
            append = false,
            ascii = true,
        ) do vtk
            GT.vtk_physical_faces!(vtk, mesh)
            GT.vtk_physical_nodes!(vtk, mesh)
            vtk["piece"] = fill(rank, sum(GT.num_faces(mesh)))
            vtk["owner"] = local_to_owner(GT.node_indices(ids))
            vtk["owner"] = face_to_owner
            vtk["node"] = local_to_global(GT.node_indices(ids))
        end
    end
    setup
end


"""
    node_coordinates(mesh, face_id, d)

Return node coordinates corresponding to the `d`-dimensional face with `face_id`

Variables matching the pattern `mesh_node*` correspond to the granularity of the supplied
`mesh`. For example, if `mesh` is a `final_mesh`, then `mesh_node_to_coordinates`
is understood as `final_mesh_node_to_coordinates`.
"""
function node_coordinates(mesh, face, d)
    n_dfaces = num_faces(mesh, d)
    dface_to_local_node_to_mesh_node = face_nodes(mesh, d)
    mesh_node_to_coordinates = node_coordinates(mesh)
    @assert face <= n_dfaces "face id is in 1:n_dfaces, got $(face) ∉ 1:$(n_dfaces)"
    local_node_to_mesh_node = dface_to_local_node_to_mesh_node[face]
    coordinates = mesh_node_to_coordinates[local_node_to_mesh_node]
    coordinates
end


function print_pmesh_coordinates(
    pmesh::GT.PMesh,
    parts,
    face_id::Int,
    face_dim::Int,
    part_id::Int = 0,
)
    @show face_id face_dim
    map(partition(pmesh), parts) do mesh, part
        # print coordinates for all parts if part_id == 0
        if part_id == 0 || part == part_id
            @show part
            display(node_coordinates(mesh, face_id, face_dim))
        end
    end
end

"""
    show_num_nodes(pmesh::GT.PMesh, parts)

Show the number of nodes on each partition of a `pmesh`.
"""
function show_num_nodes(pmesh::GT.PMesh, parts)
    @show GT.num_nodes(pmesh)
    map(partition(pmesh), parts) do mesh, part
        @show GT.num_nodes(mesh) part
    end
end

###################################################################
# Tests
###################################################################
spx0 = GT.unit_simplex(0)
spx1 = GT.unit_simplex(1)
spx2 = GT.unit_simplex(2)
spx3 = GT.unit_simplex(3)
display(spx3)

cube0 = GT.unit_n_cube(0)
cube1 = GT.unit_n_cube(1)
cube2 = GT.unit_n_cube(2)
cube3 = GT.unit_n_cube(3)
display(cube3)

@show typeof(spx0)
@show typeof(cube0)
@show typeof(spx1)
@show typeof(cube1)

degree = 4
quad = GT.default_quadrature(spx0, degree)
quad = GT.default_quadrature(spx1, degree)
quad = GT.default_quadrature(spx2, degree)
quad = GT.default_quadrature(spx3, degree)

quad = GT.default_quadrature(cube0, degree)
quad = GT.default_quadrature(cube1, degree)
quad = GT.default_quadrature(cube2, degree)
quad = GT.default_quadrature(cube3, degree)

order = 1
fe = GT.lagrange_mesh_face(spx0, order)
fe = GT.lagrange_mesh_face(spx1, order)
fe = GT.lagrange_mesh_face(spx2, order)
fe = GT.lagrange_mesh_face(spx3, order)
display(fe)

fe = GT.lagrange_mesh_face(cube0, order)
fe = GT.lagrange_mesh_face(cube1, order)
fe = GT.lagrange_mesh_face(cube2, order)
fe = GT.lagrange_mesh_face(cube3, order)
display(fe)

fe = GT.lagrange_mesh_face(cube0, order)
@show GT.monomial_exponents(fe)
@show GT.node_coordinates(fe)
fe = GT.lagrange_mesh_face(cube2, order)
@show GT.node_coordinates(fe)

spx2 = GT.unit_simplex(2)
quad = GT.default_quadrature(spx2, degree)
fe = GT.lagrange_mesh_face(spx2, order)
funs = GT.shape_functions(fe)
x = GT.coordinates(quad)
B = broadcast(GT.value, permutedims(funs), x)
display(B)
tabulator = GT.tabulator(fe)
A = tabulator(GT.value, x)
@test A ≈ B
x = GT.node_coordinates(fe)
A = tabulator(GT.value, x)

outdir = mkpath(joinpath(@__DIR__, "..", "output"))
msh = joinpath(@__DIR__, "..", "assets", "quad.msh")
mesh = GT.mesh_from_gmsh(msh)
vtk_grid(joinpath(outdir, "quad"), GT.vtk_args(mesh)...) do vtk
    GT.vtk_physical_faces!(vtk, mesh)
    GT.vtk_physical_nodes!(vtk, mesh)
end
for d = 0:GT.num_dims(mesh)
    vtk_grid(joinpath(outdir, "quad_$d"), GT.vtk_args(mesh, d)...) do vtk
        GT.vtk_physical_faces!(vtk, mesh, d)
        GT.vtk_physical_nodes!(vtk, mesh, d)
    end
end

msh = joinpath(@__DIR__, "..", "assets", "demo.msh")
mesh = GT.mesh_from_gmsh(msh; complexify = false)
vtk_grid(joinpath(outdir, "demo"), GT.vtk_args(mesh)...) do vtk
    GT.vtk_physical_faces!(vtk, mesh)
    GT.vtk_physical_nodes!(vtk, mesh)
end
for d = 0:GT.num_dims(mesh)
    vtk_grid(joinpath(outdir, "demo_$d"), GT.vtk_args(mesh, d)...) do vtk
        GT.vtk_physical_faces!(vtk, mesh, d)
        GT.vtk_physical_nodes!(vtk, mesh, d)
    end
end

@show GT.unit_simplex(0) |> GT.boundary
@show GT.unit_simplex(1) |> GT.boundary
@show GT.unit_simplex(2) |> GT.boundary
@show GT.unit_simplex(3) |> GT.boundary

@show GT.unit_n_cube(0) |> GT.boundary
@show GT.unit_n_cube(1) |> GT.boundary
@show GT.unit_n_cube(2) |> GT.boundary
@show GT.unit_n_cube(3) |> GT.boundary

GT.unit_n_cube(0) |> GT.opposite_faces
GT.unit_n_cube(1) |> GT.opposite_faces
GT.unit_n_cube(2) |> GT.opposite_faces
GT.unit_n_cube(3) |> GT.opposite_faces

GT.opposite_faces(cube3, 0)
GT.opposite_faces(cube3, 1)
GT.opposite_faces(cube3, 2)
GT.opposite_faces(cube3, 3)

order = 2
GT.lagrange_mesh_face(spx1, order) |> GT.boundary |> GT.topology

mesh = GT.mesh_from_gmsh(msh; complexify = false)

new_mesh, old_to_new = GT.complexify(mesh)

mesh = GT.mesh_from_gmsh(msh)

domain = (0, 1, 0, 1)
cells = (2, 2)
mesh = GT.cartesian_mesh(domain, cells)
mesh = GT.cartesian_mesh(domain, cells, boundary = false)
mesh = GT.cartesian_mesh(domain, cells, simplexify = true)
mesh = GT.cartesian_mesh(domain, cells, boundary = false, simplexify = true)

mesh = GT.cartesian_mesh(domain, cells)
vtk_grid(joinpath(outdir, "cartesian"), GT.vtk_args(mesh)...) do vtk
    GT.vtk_physical_faces!(vtk, mesh)
    GT.vtk_physical_nodes!(vtk, mesh)
end
for d = 0:GT.num_dims(mesh)
    vtk_grid(joinpath(outdir, "cartesian_$d"), GT.vtk_args(mesh, d)...) do vtk
        GT.vtk_physical_faces!(vtk, mesh, d)
        GT.vtk_physical_nodes!(vtk, mesh, d)
    end
end

mesh = GT.mesh_from_gmsh(msh)
face_groups = GT.physical_faces(mesh)
group_names = GT.physical_names(mesh, 2)
group_names = GT.physical_names(mesh)
group_names = GT.physical_names(mesh; merge_dims = true)
node_groups = GT.physical_nodes(mesh; merge_dims = true)
node_groups = GT.physical_nodes(mesh; merge_dims = true, disjoint = true)

vmesh, vglue = GT.visualization_mesh(mesh, 2)

mesh = GT.cartesian_mesh(domain, cells)
GT.mesh_graph(
    mesh;
    partition_strategy = GT.partition_strategy(graph_nodes = :nodes, graph_edges = :cells),
)
GT.mesh_graph(
    mesh;
    partition_strategy = GT.partition_strategy(
        graph_nodes = :nodes,
        graph_edges = :faces,
        graph_edges_dim = :all,
    ),
)
GT.mesh_graph(
    mesh;
    partition_strategy = GT.partition_strategy(
        graph_nodes = :nodes,
        graph_edges = :faces,
        graph_edges_dim = 1,
    ),
)
GT.mesh_graph(
    mesh;
    partition_strategy = GT.partition_strategy(graph_nodes = :cells, graph_edges = :nodes),
)
GT.mesh_graph(
    mesh;
    partition_strategy = GT.partition_strategy(
        graph_nodes = :cells,
        graph_edges = :faces,
        graph_edges_dim = 1,
    ),
)

np = 2
parts = DebugArray(LinearIndices((np,)))
mesh = GT.cartesian_mesh((0, 1, 0, 1), (2, 2))
cell_to_color = [1, 1, 2, 2]
pmesh = GT.partition_mesh(
    mesh,
    np;
    partition_strategy = GT.partition_strategy(
        graph_nodes = :cells,
        graph_edges = :nodes,
        ghost_layers = 0,
    ),
    parts,
    graph_partition = cell_to_color,
)
pmesh = GT.partition_mesh(
    mesh,
    np;
    partition_strategy = GT.partition_strategy(
        graph_nodes = :cells,
        graph_edges = :nodes,
        ghost_layers = 1,
    ),
    parts,
    graph_partition = cell_to_color,
)

node_to_color = [1, 1, 1, 1, 2, 2, 2, 2, 2]
pmesh = GT.partition_mesh(
    mesh,
    np;
    partition_strategy = GT.partition_strategy(
        graph_nodes = :nodes,
        graph_edges = :cells,
        ghost_layers = 0,
    ),
    parts,
    graph_partition = node_to_color,
)
pmesh = GT.partition_mesh(
    mesh,
    np;
    partition_strategy = GT.partition_strategy(
        graph_nodes = :nodes,
        graph_edges = :cells,
        ghost_layers = 1,
    ),
    parts,
    graph_partition = node_to_color,
)

function setup(mesh, ids, rank)
    face_to_owner = zeros(Int, sum(GT.num_faces(mesh)))
    D = GT.num_dims(mesh)
    for d = 0:D
        face_to_owner[GT.face_range(mesh, d)] = local_to_owner(GT.face_indices(ids, d))
    end
    pvtk_grid(
        joinpath(outdir, "pmesh"),
        GT.vtk_args(mesh)...;
        part = rank,
        nparts = np,
    ) do vtk
        GT.vtk_physical_faces!(vtk, mesh)
        GT.vtk_physical_nodes!(vtk, mesh)
        vtk["piece"] = fill(rank, sum(GT.num_faces(mesh)))
        vtk["owner"] = local_to_owner(GT.node_indices(ids))
        vtk["owner"] = face_to_owner
    end
end
map(setup, partition(pmesh), GT.index_partition(pmesh), parts)

## In this case, we do the hard work only on the main proc
np = 4
parts = DebugArray(LinearIndices((np,)))
pmesh =
    map_main(parts) do parts
        mesh = GT.mesh_from_gmsh(msh)
        partition_strategy =
            GT.partition_strategy(graph_nodes = :nodes, graph_edges = :cells)
        graph = GT.mesh_graph(mesh; partition_strategy)
        graph_partition = Metis.partition(graph, np)
        GT.partition_mesh(mesh, np; partition_strategy, graph, graph_partition)
    end |> GT.scatter_mesh

function setup(mesh, ids, rank)
    face_to_owner = zeros(Int, sum(GT.num_faces(mesh)))
    D = GT.num_dims(mesh)
    for d = 0:D
        face_to_owner[GT.face_range(mesh, d)] = local_to_owner(GT.face_indices(ids, d))
    end
    pvtk_grid(
        joinpath(outdir, "pmesh"),
        GT.vtk_args(mesh)...;
        part = rank,
        nparts = np,
    ) do vtk
        GT.vtk_physical_faces!(vtk, mesh)
        GT.vtk_physical_nodes!(vtk, mesh)
        vtk["piece"] = fill(rank, sum(GT.num_faces(mesh)))
        vtk["owner"] = local_to_owner(GT.node_indices(ids))
        vtk["owner"] = face_to_owner
    end
end
map(setup, partition(pmesh), GT.index_partition(pmesh), parts)

# In this one, we do only the graph partition on the main
# but we load the mesh everywhere
mesh = GT.mesh_from_gmsh(msh)
partition_strategy =
    GT.partition_strategy(graph_nodes = :nodes, graph_edges = :cells, ghost_layers = 1)
graph = GT.mesh_graph(mesh; partition_strategy)
graph_partition = map_main(parts) do parts
    Metis.partition(graph, np)
end |> multicast |> PartitionedArrays.getany
pmesh = GT.partition_mesh(mesh, np; partition_strategy, parts, graph, graph_partition)

graph_nodes = :cells
partition_strategy = GT.partition_strategy(graph_nodes = :cells, graph_edges = :nodes)
graph = GT.mesh_graph(mesh; partition_strategy)
graph_partition = map_main(parts) do parts
    Metis.partition(graph, np)
end |> multicast |> PartitionedArrays.getany
pmesh = GT.partition_mesh(mesh, np; partition_strategy, parts, graph, graph_partition)
pmesh = GT.partition_mesh(mesh, np; partition_strategy, parts, graph, graph_partition)

# coarse pmesh testing and visualization
domain = (0, 1, 0, 1)
cells_per_dir = (4, 4)
parts_per_dir = (2, 2)
pmesh = GT.cartesian_mesh(
    domain,
    cells_per_dir;
    parts_per_dir,
    partition_strategy = GT.partition_strategy(
        graph_nodes = :cells,
        graph_edges = :nodes,
        ghost_layers = 1,
    ),
)
pmesh = GT.cartesian_mesh(
    domain,
    cells_per_dir;
    parts_per_dir,
    partition_strategy = GT.partition_strategy(
        graph_nodes = :cells,
        graph_edges = :nodes,
        ghost_layers = 0,
    ),
)
np = prod(parts_per_dir)
parts = DebugArray(LinearIndices((np,)))
pmesh = GT.cartesian_mesh(
    domain,
    cells_per_dir;
    parts_per_dir,
    parts,
    partition_strategy = GT.partition_strategy(
        graph_nodes = :cells,
        graph_edges = :nodes,
        ghost_layers = 1,
    ),
)
pmesh = GT.cartesian_mesh(
    domain,
    cells_per_dir;
    parts_per_dir,
    parts,
    partition_strategy = GT.partition_strategy(
        graph_nodes = :cells,
        graph_edges = :nodes,
        ghost_layers = 0,
    ),
)
@test GT.partition_strategy(pmesh).ghost_layers == 0
@test GT.partition_strategy(pmesh).graph_nodes === :cells

function setup(mesh, ids, rank)
    face_to_owner = zeros(Int, sum(GT.num_faces(mesh)))
    D = GT.num_dims(mesh)
    for d = 0:D
        face_to_owner[GT.face_range(mesh, d)] = local_to_owner(GT.face_indices(ids, d))
    end
    pvtk_grid(
        joinpath(outdir, "pmesh-cartesian"),
        GT.vtk_args(mesh)...;
        part = rank,
        nparts = np,
    ) do vtk
        GT.vtk_physical_faces!(vtk, mesh)
        GT.vtk_physical_nodes!(vtk, mesh)
        vtk["piece"] = fill(rank, sum(GT.num_faces(mesh)))
        vtk["owner"] = local_to_owner(GT.node_indices(ids))
        vtk["owner"] = face_to_owner
    end
end
map(setup, partition(pmesh), GT.index_partition(pmesh), parts)

### 3D periodic puzzle piece visualizations
#=
# The file naming convention for geo files
unit_cell_<dimension>_<periodic OR nonperiodic> \
    _<geometry, e.g., puzzlepiece, square>_geometry.geo

# The file naming convention for msh files:
<unit cell geo file name>_<refcell type, e.g., quad, triangular>\
    _[, mesh dims]_refcell.msh

# The file naming convention for unit cell/coarse cell vtk files:
<unit_cell_mesh OR coarse_cell_mesh>_<dimension>_<periodic OR nonperiodic>_<gmsh OR glk>\
    _<geometry, e.g., puzzlepiece, square>_geometry_<refcell type e.g., quad, triangular>\
    _[mesh dims]_refcell.vtu

# The file naming convention for final meshes vtu (i.e., two level meshes)
[,coarse_cell_id]_final_mesh_<unit cell vtk file name>_<coarse cell vtk fname>
=#
assetsdir = joinpath(@__DIR__, "..", "assets")

## Sequential
# Load periodic fine (unit cell) mesh with triangular refcells
unit_cell_mesh_fpath =
    joinpath(assetsdir, "unit_cell_3D_periodic_puzzlepiece_geometry_triangular_refcell.msh")
unit_cell_mesh = GT.mesh_from_gmsh(unit_cell_mesh_fpath)

# visualize the periodic gmsh unit cell with triangular refcells
unit_cell_vtk_fname = "unit_cell_mesh_3D_periodic_gmsh_puzzlepiece_geometry_triangular_refcell"
visualize_mesh(unit_cell_mesh, joinpath(outdir, unit_cell_vtk_fname))

# test 4x4x4 coarse mesh
coarse_domain = (0, 10, 0, 10, 0, 10)
coarse_mesh_dims = (4, 4, 4)
coarse_mesh_4x4x4 = GT.cartesian_mesh(coarse_domain, coarse_mesh_dims)
coarse_cell_vtk_fname_4x4x4 = "coarse_cell_mesh_3D_nonperiodic_glk_box_geometry_quad_4x4x4_refcell"
visualize_mesh(coarse_mesh_4x4x4, joinpath(outdir, coarse_cell_vtk_fname_4x4x4))

# visualize final mesh with 4x4x4 coarse mesh and unit cell
periodic_final_mesh, glue = GT.two_level_mesh(coarse_mesh_4x4x4, unit_cell_mesh)
final_mesh_vtk_fname = "final_mesh_$(unit_cell_vtk_fname)_$(coarse_cell_vtk_fname_4x4x4)"
visualize_mesh(periodic_final_mesh, joinpath(outdir, final_mesh_vtk_fname))

## Parallel
# 2x2x2 part per dir
domain = (0, 10, 0, 10, 0, 10)
cells = (4, 4, 4)
parts_per_dir = (2, 2, 2)
nparts = prod(parts_per_dir)
parts = DebugArray(LinearIndices((nparts,)))
coarse_pmesh = GT.cartesian_mesh(
    domain,
    cells;
    parts_per_dir,
    parts,
    partition_strategy = GT.partition_strategy(; ghost_layers = 0),
)
coarse_pmesh_vtk_fname = "coarse_cell_pmesh_3D_nonperiodic_glk_square_geometry_quad_4x4x4_refcell_2x2x2_parts_per_direction"

final_pmesh, _ = GT.two_level_mesh(coarse_pmesh, unit_cell_mesh)

# visualize the parallel mesh
pmesh_vtk_fpath = joinpath(
    @__DIR__,
    outdir,
    "final_pmesh_$(unit_cell_vtk_fname)_$(coarse_pmesh_vtk_fname)",
)
visualize_pmesh(final_pmesh, parts, nparts, pmesh_vtk_fpath)

end # module

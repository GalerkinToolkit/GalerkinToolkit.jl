"""
To easily write sequential meshes to vtk, simply call the `visualize_mesh` function.
Similarly, call `visualize_pmesh` for parallel meshes.

# Notes for Visualizations 
For visualizations, use threshold and or translations of coordinates to aid in viewing
only desired parts of the coarse cell ids. 

In GMSH visibility you can change to elementary entities to debug visualization 

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
"""

import GalerkinToolkit as gk
using WriteVTK
using Test
using PartitionedArrays

function test_two_level_mesh_with_periodic_3D_puzzlepiece_unit_cell(outdir::String, assetsdir::String)
    # corresponds to 2D cell in glk mesh
    cell_dim = 3
    cell_id_to_inspect = 42

    ## Sequential 
    # test 1x1x1 coarse mesh
    coarse_domain = (0, 10, 0, 10, 0, 10)
    coarse_mesh_dims = (1, 1, 1)
    coarse_mesh_1x1x1 = gk.cartesian_mesh(coarse_domain, coarse_mesh_dims)
    coarse_cell_vtk_fname_1x1x1 = "coarse_cell_mesh_3D_nonperiodic_glk_box_geometry_quad_1x1x1_refcell"

    # Load periodic fine (unit cell) mesh with triangular refcells 
    unit_cell_mesh_fpath = joinpath(
      assetsdir,
        "unit_cell_3D_periodic_puzzlepiece_geometry_triangular_refcell.msh")
    unit_cell_mesh = gk.mesh_from_gmsh(unit_cell_mesh_fpath)

    # visualize the periodic gmsh unit cell with triangular refcells 
    unit_cell_vtk_fname = "unit_cell_mesh_3D_periodic_gmsh_puzzlepiece_geometry_triangular_refcell"
    visualize_mesh(unit_cell_mesh, joinpath(outdir, unit_cell_vtk_fname))

    # visualize final mesh with 1x1x1 coarse mesh and periodic unit cell 
    periodic_final_mesh, _ = gk.two_level_mesh(coarse_mesh_1x1x1, unit_cell_mesh)
    final_mesh_vtk_fname = "final_mesh_$(unit_cell_vtk_fname)_$(coarse_cell_vtk_fname_1x1x1)"
    visualize_mesh(periodic_final_mesh, joinpath(outdir, final_mesh_vtk_fname))

    # Coordinate check for unit cell in a 1x1x1 coarse mesh 
    expected_coordinates = [
        [8.31698729810778, 8.316987298107781, 10.0],
        [10.0, 8.453395893624393, 8.436185836985588],
        [10.0, 6.76776695296637, 8.232233047033633],
        [6.601456609551659, 6.630776452002998, 6.5285120063889135]
    ]
    example_coordinates = gk.node_coordinates(
        periodic_final_mesh, cell_id_to_inspect, cell_dim)
    @test all(example_coordinates .≈ expected_coordinates)

    # test 4x4x4 coarse mesh 
    coarse_domain = (0, 10, 0, 10, 0, 10)
    coarse_mesh_dims = (4, 4, 4)
    coarse_mesh_4x4x4 = gk.cartesian_mesh(coarse_domain, coarse_mesh_dims)
    coarse_cell_vtk_fname_4x4x4 = "coarse_cell_mesh_3D_nonperiodic_glk_box_geometry_quad_4x4x4_refcell"
    visualize_mesh(coarse_mesh_4x4x4, joinpath(outdir, coarse_cell_vtk_fname_4x4x4))

    # visualize final mesh with 4x4x4 coarse mesh and unit cell 
    periodic_final_mesh, glue = gk.two_level_mesh(coarse_mesh_4x4x4, unit_cell_mesh)
    final_mesh_vtk_fname = "final_mesh_$(unit_cell_vtk_fname)_$(coarse_cell_vtk_fname_4x4x4)"
    visualize_mesh(periodic_final_mesh, joinpath(outdir, final_mesh_vtk_fname))

    # coordinate check 
    expected_coordinates = [
        [2.079246824526945, 2.0792468245269453, 2.5],
        [2.5, 2.113348973406098, 2.109046459246397],
        [2.5, 1.6919417382415924, 2.0580582617584082],
        [1.6503641523879147, 1.6576941130007494, 1.6321280015972284],
    ]
    example_coordinates = gk.node_coordinates(
        periodic_final_mesh, cell_id_to_inspect, cell_dim)
    @test all(example_coordinates .≈ expected_coordinates)

    ## Parallel 
    # 1 part per dir (i.e., no parallelism)
    domain = (0, 10, 0, 10, 0, 10)
    cells = (4, 4, 4)
    parts_per_dir = (1, 1, 1)
    nparts = prod(parts_per_dir)
    parts = DebugArray(LinearIndices((nparts,)))
    coarse_pmesh = gk.cartesian_mesh(
        domain, cells;
        parts_per_dir, parts,
        partition_strategy=gk.partition_strategy(; ghost_layers=0))
    coarse_pmesh_vtk_fname = "coarse_cell_pmesh_3D_nonperiodic_glk_square_geometry_quad_4x4x4_refcell_1x1x1_parts_per_direction"

    final_pmesh, _ = gk.two_level_mesh(coarse_pmesh, unit_cell_mesh)

    # visualize the parallel mesh
    pmesh_vtk_fpath = joinpath(
        @__DIR__,
        outdir,
        "final_pmesh_$(unit_cell_vtk_fname)_$(coarse_pmesh_vtk_fname)")
    visualize_pmesh(final_pmesh, parts, nparts, pmesh_vtk_fpath)

    # coordinate check 
    expected_coordinates = DebugArray([
        [
            [2.079246824526945, 2.0792468245269453, 2.5],
            [2.5, 2.113348973406098, 2.109046459246397],
            [2.5, 1.6919417382415924, 2.0580582617584082],
            [1.6503641523879147, 1.6576941130007494, 1.6321280015972284],
        ]
    ])
    test_pmesh_coordinates(expected_coordinates, final_pmesh, cell_id_to_inspect, cell_dim)

    # 2x2x2 part per dir 
    domain = (0, 10, 0, 10, 0, 10)
    cells = (4, 4, 4)
    parts_per_dir = (2, 2, 2)
    nparts = prod(parts_per_dir)
    parts = DebugArray(LinearIndices((nparts,)))
    coarse_pmesh = gk.cartesian_mesh(
        domain, cells;
        parts_per_dir, parts,
        partition_strategy=gk.partition_strategy(; ghost_layers=0))
    coarse_pmesh_vtk_fname = "coarse_cell_pmesh_3D_nonperiodic_glk_square_geometry_quad_4x4x4_refcell_2x2x2_parts_per_direction"

    final_pmesh, _ = gk.two_level_mesh(coarse_pmesh, unit_cell_mesh)

    # visualize the parallel mesh
    pmesh_vtk_fpath = joinpath(
        @__DIR__,
        outdir,
        "final_pmesh_$(unit_cell_vtk_fname)_$(coarse_pmesh_vtk_fname)")
    visualize_pmesh(final_pmesh, parts, nparts, pmesh_vtk_fpath)

    # coordinate check 
    expected_coordinates = DebugArray([
        # p1
        [
            [2.079246824526945, 2.0792468245269453, 2.5],                                                                                                                                                
            [2.5, 2.113348973406098, 2.109046459246397],                                                                                                                                                 
            [2.5, 1.6919417382415924, 2.0580582617584082],                                                                                                                                               
            [1.6503641523879147, 1.6576941130007494, 1.6321280015972284]     
        ],
        # p2 
        [
            [7.079246824526946, 2.0792468245269453, 2.5],                                                                                                                                                
            [7.5, 2.113348973406098, 2.109046459246397],                                                                                                                                                 
            [7.5, 1.6919417382415924, 2.0580582617584082],                                                                                                                                               
            [6.650364152387915, 1.6576941130007494, 1.6321280015972284]   
        ],
        # p3 
        [
            [2.079246824526945, 7.079246824526946, 2.5],                                                                                                                                                 
            [2.5, 7.113348973406098, 2.109046459246397],                                                                                                                                                 
            [2.5, 6.691941738241592, 2.0580582617584082],                                                                                                                                                
            [1.6503641523879147, 6.657694113000749, 1.6321280015972284]
        ],
        # p4 
        [
            [7.079246824526946, 7.079246824526946, 2.5],                                                                                                                                                 
            [7.5, 7.113348973406098, 2.109046459246397],                                                                                                                                                 
            [7.5, 6.691941738241592, 2.0580582617584082],                                                                                                                                               
            [6.650364152387915, 6.657694113000749, 1.6321280015972284]
        ],
        # p5  
        [
            [2.079246824526945, 2.0792468245269453, 7.5],
            [2.5, 2.113348973406098, 7.109046459246397],
            [2.5, 1.6919417382415924, 7.058058261758408],
            [1.6503641523879147, 1.6576941130007494, 6.632128001597229]
        ],
        # p6 
        [
            [7.079246824526946, 2.0792468245269453, 7.5],
            [7.5, 2.113348973406098, 7.109046459246397],
            [7.5, 1.6919417382415924, 7.058058261758408],
            [6.650364152387915, 1.6576941130007494, 6.632128001597229]
        ],
        # p7 
        [ 
            [2.079246824526945, 7.079246824526946, 7.5],
            [2.5, 7.113348973406098, 7.109046459246397],
            [2.5, 6.691941738241592, 7.058058261758408],
            [1.6503641523879147, 6.657694113000749, 6.632128001597229]
        ],
        # p8
        [
            [7.079246824526946, 7.079246824526946, 7.5],
            [7.5, 7.113348973406098, 7.109046459246397],
            [7.5, 6.691941738241592, 7.058058261758408],
            [6.650364152387915, 6.657694113000749, 6.632128001597229]
        ]
    ])
    test_pmesh_coordinates(expected_coordinates, final_pmesh, cell_id_to_inspect, cell_dim)
end

function visualize_mesh(mesh, outpath, glue = nothing, d = nothing)
    node_ids = collect(1:gk.num_nodes(mesh))

    # NOTE: without this check, attempting to visualize a 3D periodic mesh fails because
    # the `periodic_nodes` property is a Vector and not a Pair 
    valid_periodic_mesh = gk.has_periodic_nodes(mesh) && gk.periodic_nodes(mesh) isa Pair
    if valid_periodic_mesh
        # Get periodic node info aboutcell
        periodic_nodes = gk.periodic_nodes(mesh)
        fine_pnode_to_fine_node = periodic_nodes.first
        fine_pnode_to_master_fine_node = periodic_nodes.second

        # Labeling periodic nodes 
        fine_node_to_master_fine_node = copy(node_ids)
        fine_node_to_master_fine_node[
            fine_pnode_to_fine_node] = fine_pnode_to_master_fine_node

        # Handles periodic vertices whose master node is the master of another periodic 
        # vertex assuming one level of indirection for periodicity 
        # (e.g., vertex -> master -> master)
        n_fine_nodes = length(node_ids)
        for fine_node in 1:n_fine_nodes
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
        vtk_grid(
            outpath,
            gk.vtk_args(mesh)...) do vtk
            gk.vtk_physical_faces!(vtk, mesh)
            gk.vtk_physical_nodes!(vtk, mesh)
            valid_periodic_mesh && (
                vtk["periodic_master_id"] = fine_node_to_master_fine_node)
            vtk["node_id"] = node_ids
            if !isnothing(glue)
                vtk["coarse_cell_id"] = glue.final_cell_to_coarse_cell
            end 
        end
    else
        vtk_grid(
            outpath,
            gk.vtk_args(mesh, d)...) do vtk
            gk.vtk_physical_faces!(vtk, mesh, d)
            gk.vtk_physical_nodes!(vtk, mesh, d)
            valid_periodic_mesh && (
                vtk["periodic_master_id"] = fine_node_to_master_fine_node)
            vtk["node_id"] = node_ids
            if !isnothing(glue)
                vtk["coarse_cell_id"] = glue.final_cell_to_coarse_cell
            end 
        end
    end 
end

"""
    visualize_pmesh(pmesh::gk.PMesh, parts, nparts, pmesh_vtk_fpath::String)

Writes a parallel mesh on `nparts` partitions to `pmesh_vtk_fpath`.
"""
function visualize_pmesh(pmesh::gk.PMesh, parts, nparts, pmesh_vtk_fpath::String)
    map(
        visualize_pmesh_setup(nparts, pmesh_vtk_fpath),
        partition(pmesh),
        gk.index_partition(pmesh),
        parts)
end

"""
    visualize_pmesh_setup(nparts, outpath) 

Return function that writes vtk for each part in a parallel mesh.
"""
function visualize_pmesh_setup(nparts, outpath)
    @assert isabspath(outpath) "abspath with pvtk_grid ensures function" # PR this?
    function setup(mesh, ids, rank)
        face_to_owner = zeros(Int, sum(gk.num_faces(mesh)))
        D = gk.num_dims(mesh)
        face_to_owner[gk.face_range(mesh, D)] = local_to_owner(gk.face_indices(ids, D))
        pvtk_grid(
            outpath, gk.vtk_args(mesh)...;
            part=rank, nparts=nparts, append=false, ascii=true) do vtk
            gk.vtk_physical_faces!(vtk, mesh)
            gk.vtk_physical_nodes!(vtk, mesh)
            vtk["piece"] = fill(rank, sum(gk.num_faces(mesh)))
            vtk["owner"] = local_to_owner(gk.node_indices(ids))
            vtk["owner"] = face_to_owner
            vtk["node"] = local_to_global(gk.node_indices(ids))
        end
    end
    setup
end

""" 
    test_pmesh_coordinates(
        expected_pcoordinates::DebugArray, pmesh::gk.PMesh, face_id::Int, face_dim::Int)

Use `@test` to check whether `expected_pcoordinates` matches the coordinates of a particular 
`face_id` of dimension `face_dim` in the mesh of `pmesh`.
"""
function test_pmesh_coordinates(
    expected_pcoordinates::DebugArray, pmesh::gk.PMesh, face_id::Int, face_dim::Int)
    map(partition(pmesh), expected_pcoordinates) do mesh, expected_coordinates
        example_coordinates = gk.node_coordinates(mesh, face_id, face_dim)
        @test all(example_coordinates .≈ expected_coordinates)
    end
end

function print_pmesh_coordinates(
    pmesh::gk.PMesh, parts, face_id::Int, face_dim::Int, part_id::Int = 0)
    @show face_id face_dim 
    map(partition(pmesh), parts) do mesh, part 
        # print coordinates for all parts if part_id == 0
        if part_id == 0 || part == part_id 
            @show part 
            display(gk.node_coordinates(mesh, face_id, face_dim))
        end 
    end
end

"""
    show_num_nodes(pmesh::gk.PMesh, parts)

Show the number of nodes on each partition of a `pmesh`.
"""
function show_num_nodes(pmesh::gk.PMesh, parts)
    @show gk.num_nodes(pmesh)
    map(partition(pmesh), parts) do mesh, part 
        @show gk.num_nodes(mesh) part 
    end 
end 
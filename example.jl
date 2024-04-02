"""
TODOs:
puzzle piece 2D test
parallel puzzle piece 2D test
puzzle piece 3D test
parallel puzzle piece 3D test

using DrWatson

partitionedarryasbenchmarks
check the helpers and the `runjobs`
analysis.ipynb
setup ssh tunnel with jupyter notebook and das5
jupyter-notebook --no-browser
history | grep ssh # on client
careful on job submission to DAS

TODO:
3D cartesian and naming and gmsh

geometry.jl
geometry_test.jl # clean up example.jl
descriptive names for assets and geo files
2D square_<square, triangl>_4x4
3D  

# example naming convention, only include refcell info in the mesh fname, not the geo fname
[unit_cell_mesh, two_level_mesh]_<dimension>_[periodic,nonperiodic]_[gmsh,glk]\
    _<shape, e.g., puzzlepiece, square>_shape_[quad, triangle]\
    _[mesh dims]_refcell.[geo, msh, vtu]

for testing correctness of outputs:
   given a julia object, store the integer for correct object contents (find this func, not ash)   
    take some example node ids and corresponding coordinates of interface nodes against
        the hardcoded physical coordinates that would be expected for those points

can use PartitionedArraysBenchmarks/benchmarks.jl as reference

mesh generation should be of order of multigrid solver example02

loading of sequential total mesh 

parallel mesh support in gmsh?

Go through assets/ and remove unneeded geos/meshes 

Test using puzzle piece mesh??

Make 3D periodic cube
"""
module TMP

import GalerkinToolkit as gk
using WriteVTK
using Test
using PartitionedArrays

function test_two_level_mesh_with_nonperiodic_square_unit_cell()
    # corresponds to 2D cell in glk mesh
    cell_dim = 2

    # initialize 4x4 glk unit cell mesh  
    unit_cell_domain = (0, 1, 0, 1)
    unit_cell_dims = (4, 4)
    unit_cell_mesh = gk.cartesian_mesh(unit_cell_domain, unit_cell_dims)

    # visualize the glk unit cell 
    unit_cell_vtk_fname = "unit_cell_mesh_2D_nonperiodic_glk_square_geometry_quad_4x4_refcell"
    visualize_mesh(unit_cell_mesh, joinpath("output", unit_cell_vtk_fname))

    # initialize 1x1 coarse mesh
    coarse_domain = (0, 10, 0, 10)
    coarse_mesh_dims = (1, 1)
    coarse_mesh = gk.cartesian_mesh(coarse_domain, coarse_mesh_dims)
    coarse_cell_vtk_fname_1x1 = "coarse_cell_mesh_2D_nonperiodic_glk_square_geometry_quad_1x1_refcell"

    # construct the final mesh with 1x1 coarse mesh 
    final_mesh, _ = gk.two_level_mesh(coarse_mesh, unit_cell_mesh)
    final_mesh_vtk_fname = "two_level_mesh_$(unit_cell_vtk_fname)_$(coarse_cell_vtk_fname_1x1)"
    visualize_mesh(final_mesh, joinpath("output", final_mesh_vtk_fname))

    # Coordinate check for unit cell in 1x1 coarse mesh
    final_cell_to_inspect = 12 # arbitrary 
    final_cell_to_inspect_coordinates = [
        [7.5, 5.0],
        [10.0, 5.0],
        [7.5, 7.5],
        [10.0, 7.5]
    ]
    example_coordinates = gk.node_coordinates(final_mesh, final_cell_to_inspect, cell_dim)
    @test example_coordinates == final_cell_to_inspect_coordinates

    # initialize 4x4 coarse mesh 
    coarse_domain = (0, 10, 0, 10)
    coarse_mesh_dims = (4, 4)
    coarse_mesh = gk.cartesian_mesh(coarse_domain, coarse_mesh_dims)
    coarse_cell_vtk_fname_4x4 = "coarse_cell_mesh_2D_nonperiodic_glk_square_geometry_quad_4x4_refcell"

    # construct the final mesh with a 4x4 coarse mesh     
    final_mesh, _ = gk.two_level_mesh(coarse_mesh, unit_cell_mesh)
    final_mesh_vtk_fname = "two_level_mesh_$(unit_cell_vtk_fname)_$(coarse_cell_vtk_fname_4x4)"
    visualize_mesh(final_mesh, joinpath("output", final_mesh_vtk_fname))

    # Coordinate check for unit cell in 4x4 coarse mesh 
    final_cell_to_inspect = 114 # arbitrary 
    final_cell_to_inspect_coordinates = [
        [8.125, 2.5],
        [8.75, 2.5],
        [8.125, 3.125],
        [8.75, 3.125]
    ]
    example_coordinates = gk.node_coordinates(final_mesh, final_cell_to_inspect, cell_dim)
    @test example_coordinates == final_cell_to_inspect_coordinates
end

function test_two_level_mesh_with_nonperiodic_box_unit_cell()
    # corresponds to 3D cell in glk mesh
    cell_dim = 3

    # initialize 4x4x4 glk unit cell mesh  
    unit_cell_domain = (0, 1, 0, 1, 0, 1)
    unit_cell_dims = (4, 4, 4)
    unit_cell_mesh = gk.cartesian_mesh(unit_cell_domain, unit_cell_dims)

    # visualize the glk unit cell 
    unit_cell_vtk_fname = "unit_cell_mesh_3D_nonperiodic_glk_box_geometry_quad_4x4x4_refcell"
    visualize_mesh(unit_cell_mesh, joinpath("output", unit_cell_vtk_fname))

    # Coarse 1x1x1 mesh 
    coarse_domain = (0, 10, 0, 10, 0, 10)
    coarse_mesh_dims = (1, 1, 1)
    coarse_mesh_1x1x1 = gk.cartesian_mesh(coarse_domain, coarse_mesh_dims)
    coarse_cell_vtk_fname_1x1x1 = "coarse_cell_mesh_3D_nonperiodic_glk_box_geometry_quad_1x1x1_refcell"

    # visualize final mesh with 1x1x1 coarse mesh and periodic unit cell 
    final_mesh, _ = gk.two_level_mesh(coarse_mesh_1x1x1, unit_cell_mesh)
    final_mesh_vtk_fname = "two_level_mesh_$(unit_cell_vtk_fname)_$(coarse_cell_vtk_fname_1x1x1)"
    visualize_mesh(final_mesh, joinpath("output", final_mesh_vtk_fname))

    # Coordinate check for unit cell in 1x1x1 coarse grid 
    final_cell_to_inspect = 32
    final_cell_to_inspect_coordinates = [
        [7.5, 7.5, 2.5],
        [10.0, 7.5, 2.5],
        [7.5, 10.0, 2.5],
        [10.0, 10.0, 2.5],
        [7.5, 7.5, 5.0],
        [10.0, 7.5, 5.0],
        [7.5, 10.0, 5.0],
        [10.0, 10.0, 5.0]
    ]
    example_coordinates = gk.node_coordinates(
        final_mesh, final_cell_to_inspect, cell_dim)
    @test example_coordinates == final_cell_to_inspect_coordinates

    # Coarse 4x4x4 mesh 
    coarse_domain = (0, 10, 0, 10, 0, 10)
    coarse_mesh_dims = (4, 4, 4)
    coarse_mesh_4x4x4 = gk.cartesian_mesh(coarse_domain, coarse_mesh_dims)
    coarse_cell_vtk_fname_4x4x4 = "coarse_cell_mesh_3D_nonperiodic_glk_box_geometry_quad_4x4x4_refcell"

    # visualize final mesh with 4x4x4 coarse mesh and unit cell 
    final_mesh, _ = gk.two_level_mesh(coarse_mesh_4x4x4, unit_cell_mesh)
    final_mesh_vtk_fname = "two_level_mesh_$(unit_cell_vtk_fname)_$(coarse_cell_vtk_fname_4x4x4)"
    visualize_mesh(final_mesh, joinpath("output", final_mesh_vtk_fname))

    # coordinate check for 4x4x4 coarse mesh 
    final_cell_to_inspect = 2048
    final_cell_to_inspect_coordinates = [
        [9.375, 9.375, 4.375],
        [10.0, 9.375, 4.375],
        [9.375, 10.0, 4.375],
        [10.0, 10.0, 4.375],
        [9.375, 9.375, 5.0],
        [10.0, 9.375, 5.0],
        [9.375, 10.0, 5.0],
        [10.0, 10.0, 5.0]
    ]
    example_coordinates = gk.node_coordinates(
        final_mesh, final_cell_to_inspect, cell_dim)
    @test example_coordinates == final_cell_to_inspect_coordinates
end

function test_two_level_mesh_with_periodic_square_unit_cell()
    # corresponds to 2D cell in glk mesh
    cell_dim = 2

    # Load periodic fine (unit cell) mesh with triangular refcells 
    unit_cell_mesh_fpath = joinpath(
        @__DIR__,
        "assets",
        "unit_cell_2D_periodic_square_geometry_triangular_refcell.msh")
    unit_cell_mesh = gk.mesh_from_gmsh(unit_cell_mesh_fpath)

    # visualize the periodic gmsh unit cell with triangular refcells 
    unit_cell_vtk_fname = "unit_cell_mesh_2D_periodic_gmsh_square_geometry_triangular_refcell"
    visualize_mesh(unit_cell_mesh, joinpath("output", unit_cell_vtk_fname))

    # Initialize coarse 1x1 mesh 
    coarse_domain = (0, 10, 0, 10)
    coarse_mesh_dims = (1, 1)
    coarse_mesh_1x1 = gk.cartesian_mesh(coarse_domain, coarse_mesh_dims)
    coarse_cell_vtk_fname_1x1 = "coarse_cell_mesh_2D_nonperiodic_glk_square_geometry_quad_1x1_refcell"

    # visualize final mesh with 1x1 coarse mesh and periodic unit cell 
    periodic_final_mesh, _ = gk.two_level_mesh(coarse_mesh_1x1, unit_cell_mesh)
    final_mesh_vtk_fname = "two_level_mesh_$(unit_cell_vtk_fname)_$(coarse_cell_vtk_fname_1x1)"
    visualize_mesh(periodic_final_mesh, joinpath("output", final_mesh_vtk_fname))

    # Coordinate check for unit cell in 1x1 coarse mesh 
    final_cell_to_inspect = 23
    final_cell_to_inspect_coordinates = [
        [2.4999999999999982, 0.0],
        [1.830127018922196, 1.8301270189221992],
        [0.0, 0.0]
    ]
    example_coordinates = gk.node_coordinates(
        periodic_final_mesh, final_cell_to_inspect, cell_dim)
    @test example_coordinates == final_cell_to_inspect_coordinates

    # Initialize coarse 4x4 mesh 
    coarse_domain = (0, 10, 0, 10)
    coarse_mesh_dims = (4, 4)
    coarse_mesh_4x4 = gk.cartesian_mesh(coarse_domain, coarse_mesh_dims)
    coarse_cell_vtk_fname_4x4 = "coarse_cell_mesh_2D_nonperiodic_glk_square_geometry_quad_4x4_refcell"

    # visualize final mesh with 4x4 coarse mesh and periodic unit cell 
    periodic_final_mesh, _ = gk.two_level_mesh(coarse_mesh_4x4, unit_cell_mesh)
    final_mesh_vtk_fname = "two_level_mesh_$(unit_cell_vtk_fname)_$(coarse_cell_vtk_fname_4x4)"
    visualize_mesh(periodic_final_mesh, joinpath("output", final_mesh_vtk_fname))

    # Coordinate check for unit cell in 4x4 coarse mesh
    final_cell_to_inspect = 272
    final_cell_to_inspect_coordinates = [
        [6.874999999999997, 2.5],
        [6.562499999999998, 3.0412658773652734],
        [6.2499999999999964, 2.5]
    ]
    example_coordinates = gk.node_coordinates(
        periodic_final_mesh, final_cell_to_inspect, cell_dim)
    @test example_coordinates == final_cell_to_inspect_coordinates
end

function test_two_level_mesh_with_periodic_box_unit_cell()
    # corresponds to 3D cell in glk mesh
    cell_dim = 3

    # Coarse 1x1x1 mesh 
    coarse_domain = (0, 10, 0, 10, 0, 10)
    coarse_mesh_dims = (1, 1, 1)
    coarse_mesh_1x1x1 = gk.cartesian_mesh(coarse_domain, coarse_mesh_dims)
    coarse_cell_vtk_fname_1x1x1 = "coarse_cell_mesh_3D_periodic_glk_box_geometry_quad_1x1x1_refcell"

    # Load periodic fine (unit cell) mesh with triangular refcells 
    unit_cell_mesh_fpath = joinpath(
        @__DIR__,
        "assets",
        "unit_cell_3D_periodic_box_geometry_triangular_refcell.msh")
    unit_cell_mesh = gk.mesh_from_gmsh(unit_cell_mesh_fpath)

    # visualize the periodic gmsh unit cell with triangular refcells 
    unit_cell_vtk_fname = "unit_cell_mesh_3D_periodic_gmsh_box_geometry_triangular_refcell"
    visualize_mesh(unit_cell_mesh, joinpath("output", unit_cell_vtk_fname))

    # visualize final mesh with 1x1x1 coarse mesh and periodic unit cell 
    periodic_final_mesh, _ = gk.two_level_mesh(coarse_mesh_1x1x1, unit_cell_mesh)
    final_mesh_vtk_fname = "two_level_mesh_$(unit_cell_vtk_fname)_$(coarse_cell_vtk_fname_1x1x1)"
    visualize_mesh(periodic_final_mesh, joinpath("output", final_mesh_vtk_fname))

    # Coordinate check for unit cell in a 1x1x1 coarse mesh 
    final_cell_to_inspect = 112
    final_cell_to_inspect_coordinates = [
        [10.0, 5.007190394081739, 2.4743991828389467],
        [10.0, 3.333333333333324, 0.0],
        [10.0, 2.423197548516857, 2.423197548516857],
        [7.576802451483146, 2.423197548516848, 0.0]
    ]
    example_coordinates = gk.node_coordinates(
        periodic_final_mesh, final_cell_to_inspect, cell_dim)
    @test example_coordinates == final_cell_to_inspect_coordinates

    # Coarse 4x4x4 mesh 
    coarse_domain = (0, 10, 0, 10, 0, 10)
    coarse_mesh_dims = (4, 4, 4)
    coarse_mesh_4x4x4 = gk.cartesian_mesh(coarse_domain, coarse_mesh_dims)
    coarse_cell_vtk_fname_4x4x4 = "coarse_cell_mesh_3D_periodic_glk_box_geometry_quad_4x4x4_refcell"

    # visualize final mesh with 4x4x4 coarse mesh and unit cell 
    periodic_final_mesh, _ = gk.two_level_mesh(coarse_mesh_4x4x4, unit_cell_mesh)
    final_mesh_vtk_fname = "two_level_mesh_$(unit_cell_vtk_fname)_$(coarse_cell_vtk_fname_4x4x4)"
    visualize_mesh(periodic_final_mesh, joinpath("output", final_mesh_vtk_fname))

    # Coordinate check for unit cell in 4x4x4 coarse grid 
    final_cell_to_inspect = 1096
    final_cell_to_inspect_coordinates = [
        [3.751797598520435, 3.118599795709737, 0.0],
        [3.1057993871292147, 3.1057993871292147, 0.0],
        [3.7501979036595494, 2.5, 0.6185997957097367],
        [3.4252014608865786, 3.2937309998746542, 0.7928884928858897]
    ]
    example_coordinates = gk.node_coordinates(
        periodic_final_mesh, final_cell_to_inspect, cell_dim)
    @test example_coordinates == final_cell_to_inspect_coordinates
end

function test_two_level_mesh_with_periodic_2D_puzzlepiece_unit_cell()
    # corresponds to 2D cell in glk mesh
    cell_dim = 2
    cell_id_to_inspect = 42

    ## Sequential 
    # Load periodic fine (unit cell) mesh
    unit_cell_mesh_fpath = joinpath(
        @__DIR__,
        "assets",
        "unit_cell_2D_periodic_puzzlepiece_geometry_triangular_refcell.msh")
    unit_cell_mesh = gk.mesh_from_gmsh(unit_cell_mesh_fpath)

    # visualize the periodic gmsh unit cell 
    unit_cell_vtk_fname = "unit_cell_mesh_2D_periodic_gmsh_puzzlepiece_geometry_triangular_refcell"
    visualize_mesh(unit_cell_mesh, joinpath("output", unit_cell_vtk_fname))

    # coarse 1x1 mesh
    coarse_domain = (0, 10, 0, 10)
    coarse_mesh_dims = (1, 1)
    coarse_mesh_1x1 = gk.cartesian_mesh(coarse_domain, coarse_mesh_dims)
    coarse_cell_vtk_fname_1x1 = "coarse_cell_mesh_2D_nonperiodic_glk_square_geometry_quad_1x1_refcell"

    # visualize final mesh with 1x1 coarse mesh and puzzle piece unit cell
    periodic_final_mesh, _ = gk.two_level_mesh(coarse_mesh_1x1, unit_cell_mesh)
    final_mesh_vtk_fname = "two_level_mesh_$(unit_cell_vtk_fname)_$(coarse_cell_vtk_fname_1x1)"
    visualize_mesh(periodic_final_mesh, joinpath("output", final_mesh_vtk_fname))

    # Check coordinates of example refcell in final mesh
    final_cell_to_inspect = 64
    final_cell_to_inspect_coordinates = [
        [7.1650635094611, 8.750000000000007],
        [8.148393354515463, 9.262620317047551],
        [7.5, 10.0]
    ]
    example_coordinates = gk.node_coordinates(
        periodic_final_mesh, final_cell_to_inspect, cell_dim)
    @test example_coordinates == final_cell_to_inspect_coordinates

    # coarse 4x4 mesh 
    coarse_domain = (0, 10, 0, 10)
    coarse_mesh_dims = (4, 4)
    coarse_mesh_4x4 = gk.cartesian_mesh(coarse_domain, coarse_mesh_dims)
    coarse_cell_vtk_fname_4x4 = "coarse_cell_mesh_2D_nonperiodic_glk_square_geometry_quad_4x4_refcell"

    # visualize final mesh with 4x4 coarse mesh and puzzle piece unit cell
    periodic_final_mesh, _ = gk.two_level_mesh(coarse_mesh_4x4, unit_cell_mesh)
    final_mesh_vtk_fname = "two_level_mesh_$(unit_cell_vtk_fname)_$(coarse_cell_vtk_fname_4x4)"
    visualize_mesh(periodic_final_mesh, joinpath("output", final_mesh_vtk_fname))

    # Check coordinates  
    final_cell_to_inspect = 256
    final_cell_to_inspect_coordinates = [
        [5.0, 2.187500000000001],
        [5.229888333973762, 2.270164229732711],
        [5.0, 2.5]
    ]
    example_coordinates = gk.node_coordinates(
        periodic_final_mesh, final_cell_to_inspect, cell_dim)
    @test example_coordinates == final_cell_to_inspect_coordinates

    ## Parallel
    # 1 part per dir (i.e., no parallelism)
    domain = (0, 10, 0, 10)
    cells = (4, 4)
    parts_per_dir = (1, 1)
    nparts = prod(parts_per_dir)
    parts = DebugArray(LinearIndices((nparts,)))
    coarse_pmesh = gk.cartesian_mesh(
        domain, cells;
        parts_per_dir, parts,
        partition_strategy=gk.partition_strategy(; ghost_layers=0))
    coarse_pmesh_vtk_fname = "coarse_cell_pmesh_2D_nonperiodic_glk_square_geometry_quad_4x4_refcell_1x1_parts_per_direction"

    # generate parallel mesh 
    periodic_final_pmesh, _ = gk.two_level_mesh(coarse_pmesh, unit_cell_mesh)

    # visualize the parallel mesh
    pmesh_vtk_fpath = joinpath(
        @__DIR__,
        "output",
        "final_pmesh_$(unit_cell_vtk_fname)_$(coarse_pmesh_vtk_fname)")

    visualize_pmesh(periodic_final_pmesh, parts, nparts, pmesh_vtk_fpath)

    # Coordinate check for parallel single mesh 
    expected_coordinates = DebugArray([
        [
            [2.187499999999999, 2.5],
            [2.2701477152063223, 2.270121200692494],
            [2.5, 2.5]
        ]
    ])

    test_pmesh_coordinates(
        expected_coordinates, periodic_final_pmesh, cell_id_to_inspect, cell_dim)

    # parallel coarse mesh 2x2 parts per dir 
    domain = (0, 10, 0, 10)
    cells = (4, 4)
    parts_per_dir = (2, 2)
    nparts = prod(parts_per_dir)
    parts = DebugArray(LinearIndices((nparts,)))
    coarse_pmesh = gk.cartesian_mesh(
        domain, cells;
        parts_per_dir, parts,
        partition_strategy=gk.partition_strategy(; ghost_layers=0))
    coarse_pmesh_vtk_fname = "coarse_cell_pmesh_2D_nonperiodic_glk_square_geometry_quad_4x4_refcell_2x2_parts_per_direction"

    # generate parallel mesh 
    periodic_final_pmesh, _ = gk.two_level_mesh(coarse_pmesh, unit_cell_mesh)

    # visualize the parallel mesh
    pmesh_vtk_fpath = joinpath(
        @__DIR__,
        "output",
        "final_pmesh_$(unit_cell_vtk_fname)_$(coarse_pmesh_vtk_fname)")

    visualize_pmesh(periodic_final_pmesh, parts, nparts, pmesh_vtk_fpath)

    # Check whether coordinates on different partitions match the expected coordinates 
    expected_coordinates = DebugArray([
        # mesh 1
        [
            [2.187499999999999, 2.5],
            [2.2701477152063223, 2.270121200692494],
            [2.5, 2.5]
        ],
        # mesh 2
        [
            [7.187499999999999, 2.5],
            [7.270147715206322, 2.270121200692494],
            [7.5, 2.5]
        ],
        # mesh 3
        [
            [2.187499999999999, 7.5],
            [2.2701477152063223, 7.270121200692493],
            [2.5, 7.5]
        ],
        # mesh 4
        [
            [7.187499999999999, 7.5],
            [7.270147715206322, 7.270121200692493],
            [7.5, 7.5]
        ],
    ])

    test_pmesh_coordinates(
        expected_coordinates, periodic_final_pmesh, cell_id_to_inspect, cell_dim)

    # TODO: Assert ownership of nodes on boundaries 
    test_pmesh_node_ownership(periodic_final_pmesh, parts)

    # test_pmesh_cell_ownership(cell_ownership, periodic_final_pmesh)
end

function test_two_level_mesh_with_periodic_3D_puzzlepiece_unit_cell()
    # corresponds to 2D cell in glk mesh
    cell_dim = 3
    cell_id_to_inspect = 42

    ## Sequential 
    # test 1x1x1 coarse mesh
    coarse_domain = (0, 10, 0, 10, 0, 10)
    coarse_mesh_dims = (1, 1, 1)
    coarse_mesh_1x1x1 = gk.cartesian_mesh(coarse_domain, coarse_mesh_dims)
    coarse_cell_vtk_fname_1x1x1 = "coarse_cell_mesh_3D_periodic_glk_box_geometry_quad_1x1x1_refcell"

    # Load periodic fine (unit cell) mesh with triangular refcells 
    unit_cell_mesh_fpath = joinpath(
        @__DIR__,
        "assets",
        "unit_cell_3D_periodic_puzzlepiece_geometry_triangular_refcell.msh")
    unit_cell_mesh = gk.mesh_from_gmsh(unit_cell_mesh_fpath)

    # visualize the periodic gmsh unit cell with triangular refcells 
    unit_cell_vtk_fname = "unit_cell_mesh_3D_periodic_gmsh_puzzlepiece_geometry_triangular_refcell"
    visualize_mesh(unit_cell_mesh, joinpath("output", unit_cell_vtk_fname))

    # visualize final mesh with 1x1x1 coarse mesh and periodic unit cell 
    periodic_final_mesh, _ = gk.two_level_mesh(coarse_mesh_1x1x1, unit_cell_mesh)
    final_mesh_vtk_fname = "two_level_mesh_$(unit_cell_vtk_fname)_$(coarse_cell_vtk_fname_1x1x1)"
    visualize_mesh(periodic_final_mesh, joinpath("output", final_mesh_vtk_fname))

    # TODO Coordinate check for unit cell in a 1x1x1 coarse mesh 
    expected_coordinates = [
        [2.462298894969888, 4.330483324970414, 4.387785104676822],
        [2.6085458181788184, 5.1674805956749905, 3.037271336319951],
        [3.4823714583516914, 6.455029775322512, 3.3900969204290567],
        [4.588614680940672, 5.402023185913057, 3.459281419337512]
    ]
    example_coordinates = gk.node_coordinates(
        periodic_final_mesh, cell_id_to_inspect, cell_dim)
    @test example_coordinates == expected_coordinates

    # test 2x2x2 coarse mesh
    coarse_domain = (0, 10, 0, 10, 0, 10)
    coarse_mesh_dims = (2, 2, 2)
    coarse_mesh_2x2x2 = gk.cartesian_mesh(coarse_domain, coarse_mesh_dims)
    coarse_cell_vtk_fname_2x2x2 = "coarse_cell_mesh_3D_periodic_glk_box_geometry_quad_2x2x2_refcell"

    # visualize final mesh with 2x2x2 coarse mesh and unit cell 
    periodic_final_mesh, _ = gk.two_level_mesh(coarse_mesh_2x2x2, unit_cell_mesh)
    final_mesh_vtk_fname = "two_level_mesh_$(unit_cell_vtk_fname)_$(coarse_cell_vtk_fname_2x2x2)"
    visualize_mesh(periodic_final_mesh, joinpath("output", final_mesh_vtk_fname))

    throw("problem with 3D periodic puzzlepiece sequential mesh visualization")

    # coarse_domain = (0, 10, 0, 10, 0, 10)
    # coarse_mesh_dims = (4, 4, 4)
    # coarse_mesh_4x4x4 = gk.cartesian_mesh(coarse_domain, coarse_mesh_dims)
    # coarse_cell_vtk_fname_4x4x4 = "coarse_cell_mesh_3D_periodic_glk_box_geometry_quad_4x4x4_refcell"


    ## Parallel 
    # test one mesh one process

    # test 2x2x2 parts per dir 
end

function visualize_mesh(mesh, outpath)
    node_ids = collect(1:gk.num_nodes(mesh))

    # TODO: without this check, attempting to visualize a 3D periodic mesh fails because
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
    vtk_grid(
        outpath,
        gk.vtk_args(mesh)...) do vtk
        gk.vtk_physical_faces!(vtk, mesh)
        gk.vtk_physical_nodes!(vtk, mesh)
        valid_periodic_mesh && (
            vtk["periodic_master_id"] = fine_node_to_master_fine_node)
        vtk["node_id"] = node_ids
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
    @assert isabspath(outpath) "abspath with pvtk_grid ensures function" # TODO: PR this?
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
        @test example_coordinates == expected_coordinates
    end
end

function print_pmesh_coordinates(
    pmesh::gk.PMesh, parts, part_id::Int, face_id::Int, face_dim::Int)
    map(partition(pmesh), parts) do mesh, part 
        if part == part_id 
            @show face_id face_dim 
            display(gk.node_coordinates(mesh, face_id, face_dim))
        end 
    end
end

function show_num_nodes(pmesh::gk.PMesh, parts)
    @show gk.num_nodes(pmesh)
    map(partition(pmesh), parts) do mesh, part 
        @show gk.num_nodes(mesh) part 
    end 
end 

"""

In the current scheme, the partition ids in the 2D case are

```
[3 4
 1 2]
```
and ownership of a node should belong to the partition id with the highest partition id 
at the interfaces. So do something like 

for each mesh in pmesh 
    myid = mesh.id 
    function intersection_owned_by_mesh_4(interface)
        @test count(owner -> owner == 4, owners(interface)) == 1 
        # test that the intersection node is at the center of the partitions ???
    end

    # could get physical group 2 (top) and 4 (right) for mesh 1
    # top owned by mesh 3 and right owned by mesh 2 
    if myid == 1
        top = get_top(mesh)
        right = get_right(mesh)
        intersection_owned_by_mesh_4(top)
        intersection_owned_by_mesh_4(right)
        @test all(map(owner -> owner == 3 || owner == 4, unique(owners(top))))
        @test all(map(owner -> owner == 2 || owner == 4, unique(owners(right))))

    # could get physical group 3 (left) and 2 (top) for mesh 2 
    # left owned by mesh 2 and top owned by mesh 4
    elseif myid == 2
        left = get_left(mesh)
        top = get_top(mesh)
        intersection_owned_by_mesh_4(left)
        @test all(map(owner -> owner == 2 || owner == 4, unique(owners(left))))
        @test all(unique(owners(top)) .== 4)

    # could get physical group 1 (bottom) and 4 (right) for mesh 3
    # bottom owned by mesh 3 and right owned by mesh 4 
    elseif myid == 3
        bottom = get_bottom(mesh)
        right = get_right(mesh)
        intersection_owned_by_mesh_4(bottom)
        @test all(map(owner -> owner == 3 || owner == 4, unique(owners(bottom))))
        @test all(unique(owners(right)) .== 4)

    # could get physical group 3 (left) and (1) bottom for mesh 4
    # all owned by mesh 4
    elseif myid == 4
        left = get_left(mesh)
        bottom = get_bottom(mesh)
        @test all(unique(owners(left)) .== 4)
        @test all(unique(owners(bottom)) .== 4)

    else 
        throw("notimplemented")
    end


point at (5, 5, 0) --- given current domain -- should have 4 has owner for all meshes 
"""
function test_pmesh_node_ownership(pmesh, parts)
    map(partition(pmesh), gk.index_partition(pmesh), parts) do mesh, ids, part
        node_indices = gk.node_indices(ids)
        node_index_to_owner = local_to_owner(node_indices)

    end
end

function test_pmesh_cell_ownership(cell_ownership, pmesh)
    # this is guaranteed currently by uniform mesh partition, but could iterate 
    # all cells on partition and verify ownership
end

#TMP.test_two_level_mesh_with_nonperiodic_square_unit_cell()
#TMP.test_two_level_mesh_with_periodic_square_unit_cell()

#TMP.test_two_level_mesh_with_nonperiodic_box_unit_cell()
#TMP.test_two_level_mesh_with_periodic_box_unit_cell()

# TMP.test_two_level_mesh_with_periodic_2D_puzzlepiece_unit_cell()
TMP.test_two_level_mesh_with_periodic_3D_puzzlepiece_unit_cell()

end # module TMP

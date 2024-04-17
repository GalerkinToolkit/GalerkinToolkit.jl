module TMP

import GalerkinToolkit as gk
using GalerkinToolkitExamples: Example002, Example001
using Test
using PartitionedArrays
using PetscCall
using Metis
using WriteVTK

function test_solver_periodic_2D_puzzle_piece_mesh()
    tol = 1.0e-8
    
    # Puzzle piece 2D
    # The puzzle piece has gaps for which there is no labeling of phhysical groups 
    # and therefore the solution on these gaps essentially does not occur, resulting
    # in failure of thetests 
    unit_cell_mesh_fpath = joinpath(
        @__DIR__,
        "..",
        "..",
        "assets",
        "unit_cell_2D_periodic_puzzlepiece_geometry_triangular_refcell.msh")
    fine_mesh = gk.mesh_from_gmsh(unit_cell_mesh_fpath)
    domain = (0,10,0,10)
    cells = (2,2)
    coarse_mesh = gk.cartesian_mesh(domain,cells)
    final_mesh, = gk.two_level_mesh(coarse_mesh,fine_mesh)
    gk.label_boundary_faces!(final_mesh; physical_name="boundary")

    params = Dict{Symbol,Any}()
    params[:mesh] = final_mesh
    params[:dirichlet_tags] = ["boundary"] 
    params[:example_path] = joinpath(@__DIR__, "output", "puzzle_piece_2D_serial_test")
    params[:export_vtu] = true 
    results = Example001.main(params)
    @test results[:eh1] < tol
    @test results[:el2] < tol

end 

function test_solver_periodic_2D_puzzle_piece_pmesh()
    tol = 1.0e-8

    # Coarse 2D pmesh
    domain = (0, 10, 0, 10)
    cells = (4, 4)
    parts_per_dir = (2, 2)
    nparts = prod(parts_per_dir)
    parts = DebugArray(LinearIndices((nparts,)))
    coarse_pmesh = gk.cartesian_mesh(
        domain, cells;
        parts_per_dir, parts,
        partition_strategy=gk.partition_strategy(; ghost_layers=0))

    # Periodic 2D puzzle tests 
    params = Dict{Symbol,Any}()

    # load unit cell 
    unit_cell_mesh_fpath = joinpath(
        @__DIR__,
        "..",
        "..",
        "assets",
        "unit_cell_2D_periodic_puzzlepiece_geometry_triangular_refcell.msh")
    fine_mesh = gk.mesh_from_gmsh(unit_cell_mesh_fpath)

    final_pmesh, _ = gk.two_level_mesh(coarse_pmesh, fine_mesh)
    gk.label_boundary_faces!(final_pmesh; physical_name="boundary")

    params[:mesh] = final_pmesh 
    params[:dirichlet_tags] = ["boundary"]
    params[:example_path] = joinpath(@__DIR__, "output", "puzzle_piece_2D_np_4_test")
    params[:export_vtu] = true 
    results = Example002.main(params)
    @test results[:eh1] < tol 
    @test results[:el2] < tol
end

function test_solver_periodic_3D_puzzle_piece_mesh()
    tol = 1.0e-8

    # Load periodic fine (unit cell) mesh with triangular refcells 
    unit_cell_mesh_fpath = joinpath(
        @__DIR__,
        "..",
        "..",
        "assets",
        "unit_cell_3D_periodic_puzzlepiece_geometry_triangular_refcell.msh")
    fine_mesh = gk.mesh_from_gmsh(unit_cell_mesh_fpath)

    # test 4x4x4 coarse mesh 
    coarse_domain = (0, 10, 0, 10, 0, 10)
    coarse_mesh_dims = (4, 4, 4)
    coarse_mesh = gk.cartesian_mesh(coarse_domain, coarse_mesh_dims)

    # Construct mesh
    final_mesh, = gk.two_level_mesh(coarse_mesh,fine_mesh)
    gk.label_boundary_faces!(final_mesh; physical_name="boundary")

    # Call solver
    params = Dict{Symbol,Any}()
    params[:mesh] = final_mesh
    params[:dirichlet_tags] = ["boundary"] 
    params[:example_path] = joinpath(@__DIR__, "output", "puzzle_piece_3D_serial_test")
    params[:export_vtu] = true 
    results = Example001.main(params)
    @test results[:eh1] < tol
    @test results[:el2] < tol

end 

function test_solver_periodic_3D_puzzle_piece_pmesh()

    tol = 1.0e-8

    # Load periodic fine (unit cell) mesh with triangular refcells 
    unit_cell_mesh_fpath = joinpath(
        @__DIR__,
        "..",
        "..",
        "assets",
        "unit_cell_3D_periodic_puzzlepiece_geometry_triangular_refcell.msh")
    fine_mesh = gk.mesh_from_gmsh(unit_cell_mesh_fpath)

    # Coarse 4x4x4 domain with 2x2x2 parts 
    domain = (0, 10, 0, 10, 0, 10)
    cells = (4, 4, 4)
    parts_per_dir = (2, 2, 2)
    nparts = prod(parts_per_dir)
    parts = DebugArray(LinearIndices((nparts,)))
    coarse_pmesh = gk.cartesian_mesh(
        domain, cells;
        parts_per_dir, parts,
        partition_strategy=gk.partition_strategy(; ghost_layers=0))

    # Construct final mesh
    final_pmesh, _ = gk.two_level_mesh(coarse_pmesh, fine_mesh) 
    gk.label_boundary_faces!(final_pmesh; physical_name="boundary")

    # Call solver
    params = Dict{Symbol,Any}()
    params[:mesh] = final_pmesh 
    params[:dirichlet_tags] = ["boundary"]
    params[:example_path] = joinpath(@__DIR__, "output", "puzzle_piece_3D_np_4_test")
    params[:export_vtu] = true 
    results = Example002.main(params)
    @test results[:eh1] < tol 
    @test results[:el2] < tol
end 

function test_solver_periodic_3D_box_mesh()
    tol = 1.0e-8
    
    # load unit cell 
    unit_cell_mesh_fpath = joinpath(
        @__DIR__,
        "..",
        "..",
        "assets",
        "unit_cell_3D_periodic_box_geometry_triangular_refcell.msh")
    fine_mesh = gk.mesh_from_gmsh(unit_cell_mesh_fpath)

    # Coarse 4x4x4 mesh 
    coarse_domain = (0, 10, 0, 10, 0, 10)
    coarse_mesh_dims = (4, 4, 4)
    coarse_mesh = gk.cartesian_mesh(coarse_domain, coarse_mesh_dims)

    # Construct final mesh 
    final_mesh, _  = gk.two_level_mesh(coarse_mesh,fine_mesh)

    # Call solver
    params = Dict{Symbol,Any}()
    params[:mesh] = final_mesh
    params[:dirichlet_tags] = ["boundary"] 
    params[:example_path] = joinpath(@__DIR__, "output", "periodic_box_3D_serial_test")
    params[:export_vtu] = true 
    results = Example001.main(params)
    @test results[:eh1] < tol
    @test results[:el2] < tol

end 

function test_solver_periodic_3D_box_pmesh()
    tol = 1.0e-8    

    # load unit cell 
    unit_cell_mesh_fpath = joinpath(
        @__DIR__,
        "..",
        "..",
        "assets",
        "unit_cell_3D_periodic_box_geometry_triangular_refcell.msh")
    fine_mesh = gk.mesh_from_gmsh(unit_cell_mesh_fpath)

    # Initialize coarse 4x4x4 mesh with 2x2x2 parts 
    domain = (0, 10, 0, 10, 0, 10)
    cells = (4, 4, 4)
    parts_per_dir = (2, 2, 2)
    nparts = prod(parts_per_dir)
    parts = DebugArray(LinearIndices((nparts,)))
    coarse_pmesh = gk.cartesian_mesh(
        domain, cells;
        parts_per_dir, parts,
        partition_strategy=gk.partition_strategy(; ghost_layers=0))

    # Construct final pmesh
    final_pmesh, _ = gk.two_level_mesh(coarse_pmesh, fine_mesh)

    # Call solver
    params = Dict{Symbol,Any}()
    params[:mesh] = final_pmesh 
    params[:dirichlet_tags] = ["boundary"]
    params[:example_path] = joinpath(@__DIR__, "output", "periodic_box_3D_np_4_test")
    params[:export_vtu] = true 
    results = Example002.main(params)
    @test results[:eh1] < tol 
    @test results[:el2] < tol
end 

TMP.test_solver_periodic_2D_puzzle_piece_mesh()
TMP.test_solver_periodic_2D_puzzle_piece_pmesh()

TMP.test_solver_periodic_3D_puzzle_piece_mesh()
TMP.test_solver_periodic_3D_puzzle_piece_pmesh()

TMP.test_solver_periodic_3D_box_mesh()
TMP.test_solver_periodic_3D_box_pmesh()

end # module TMP

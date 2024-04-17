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
    mesh, = gk.two_level_mesh(coarse_mesh,fine_mesh)
    gk.label_boundary_faces!(mesh; physical_name="boundary")

    params = Dict{Symbol,Any}()
    params[:mesh] = mesh
    params[:dirichlet_tags] = ["boundary"] 
    params[:example_path] = joinpath(@__DIR__, "output", "puzzle_piece_2D_serial_test")
    params[:export_vtu] = true 
    results = Example001.main(params)
    @test results[:eh1] < tol
    @test results[:el2] < tol

end 

function test_solver_periodic_2D_puzzle_piece_pmesh()
    tol = 1.0e-8

    # ## Coarse 2D pmesh
    domain = (0, 10, 0, 10)
    cells = (4, 4)
    parts_per_dir = (2, 2)
    nparts = prod(parts_per_dir)
    parts = DebugArray(LinearIndices((nparts,)))
    coarse_pmesh = gk.cartesian_mesh(
        domain, cells;
        parts_per_dir, parts,
        partition_strategy=gk.partition_strategy(; ghost_layers=0))


    ## Periodic 2D puzzle tests 
    ## TODO: This test fails 
    # using default solver 
    params = Dict{Symbol,Any}()

    # load unit cell 
    unit_cell_mesh_fpath = joinpath(
        @__DIR__,
        "..",
        "..",
        "assets",
        "unit_cell_2D_periodic_puzzlepiece_geometry_triangular_refcell.msh")
    unit_cell_mesh = gk.mesh_from_gmsh(unit_cell_mesh_fpath)

    final_pmesh, _ = gk.two_level_mesh(coarse_pmesh, unit_cell_mesh)
    gk.label_boundary_faces!(final_pmesh; physical_name="boundary")

    params[:mesh] = final_pmesh 
    params[:dirichlet_tags] = ["boundary"]
    params[:example_path] = joinpath(@__DIR__, "output", "puzzle_piece_2D_np_4_test")
    params[:export_vtu] = true 
    results = Example002.main(params)
    @test results[:eh1] < tol 
    @test results[:el2] < tol
end 

end # module TMP

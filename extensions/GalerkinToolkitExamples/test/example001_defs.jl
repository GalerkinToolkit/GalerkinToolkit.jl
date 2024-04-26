import GalerkinToolkit as gk
using GalerkinToolkitExamples: Example001
using Test
using PetscCall

function test_solver_periodic_3D_puzzle_piece_mesh(outdir::String, assetsdir::String, tol)
    # Load periodic fine (unit cell) mesh with triangular refcells 
    unit_cell_mesh_fpath = joinpath(
       assetsdir,
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
    params[:example_path] = joinpath(outdir, "puzzle_piece_3D_serial_test")
    params[:export_vtu] = true 
    results = Example001.main(params)
    @test results[:eh1] < tol
    @test results[:el2] < tol
end 

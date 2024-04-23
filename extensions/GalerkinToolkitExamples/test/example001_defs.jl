import GalerkinToolkit as gk
using GalerkinToolkitExamples: Example001
using Test
using PetscCall

function test_solver_periodic_2D_square_mesh(outdir::String, assetsdir::String, tol) 
    # load unit cell 
    unit_cell_mesh_fpath = joinpath(
       assetsdir,
        "unit_cell_2D_periodic_square_geometry_triangular_refcell.msh")
    fine_mesh = gk.mesh_from_gmsh(unit_cell_mesh_fpath)

    # coarse mesh 
    domain = (0,10,0,10)
    cells = (2,2)
    coarse_mesh = gk.cartesian_mesh(domain,cells)

    # final mesh 
    final_mesh, = gk.two_level_mesh(coarse_mesh,fine_mesh)

    # call solver 
    params = Dict{Symbol,Any}()
    params[:mesh] = final_mesh
    params[:dirichlet_tags] = ["boundary"] 
    params[:example_path] = joinpath(outdir, "periodic_square_2D_serial_test")
    params[:export_vtu] = true 
    results = Example001.main(params)
    @test results[:eh1] < tol
    @test results[:el2] < tol
end

function test_solver_periodic_2D_puzzle_piece_mesh(outdir::String, assetsdir::String, tol)
    # Puzzle piece 2D
    # The puzzle piece has gaps for which there is no labeling of phhysical groups 
    # and therefore the solution on these gaps essentially does not occur, resulting
    # in failure of thetests 
    unit_cell_mesh_fpath = joinpath(
       assetsdir,
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
    params[:example_path] = joinpath(outdir, "puzzle_piece_2D_serial_test")
    params[:export_vtu] = true 
    results = Example001.main(params)
    @test results[:eh1] < tol
    @test results[:el2] < tol

end 

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

function test_solver_periodic_3D_box_mesh(outdir::String, assetsdir::String, tol)
    # load unit cell 
    unit_cell_mesh_fpath = joinpath(
       assetsdir,
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
    params[:example_path] = joinpath(outdir, "periodic_box_3D_serial_test")
    params[:export_vtu] = true 
    results = Example001.main(params)
    @test results[:eh1] < tol
    @test results[:el2] < tol

end 
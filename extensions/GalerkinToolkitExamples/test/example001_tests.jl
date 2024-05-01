module Example001Tests

import GalerkinToolkit as gk
using GalerkinToolkitExamples: Example001
using Test

tol = 1.0e-10

params = Dict{Symbol,Any}()
params[:mesh] = gk.cartesian_mesh((0,1,0,1),(3,3))
results = Example001.main(params)
@test results[:eh1] < tol
@test results[:el2] < tol

params = Dict{Symbol,Any}()
params[:mesh] = gk.cartesian_mesh((0,3,0,2),(5,10),complexify=false)
params[:dirichlet_tags] = ["1-face-1","1-face-3","1-face-4"]
params[:neumann_tags] = ["1-face-2"]
params[:u] = (x) -> sum(x)
params[:f] = (x) -> 0.0
params[:g] = (x) -> 1.0
results = Example001.main(params)
@test results[:eh1] < tol
@test results[:el2] < tol
@test results[:ncells] == 10*5

params = Dict{Symbol,Any}()
params[:mesh] = gk.cartesian_mesh((0,3,0,2,0,1),(5,5,5))
results = Example001.main(params)
@test results[:eh1] < tol
@test results[:el2] < tol
@test results[:ncells] == 5*5*5

params = Dict{Symbol,Any}()
params[:mesh] = gk.cartesian_mesh((0,3,0,2),(5,5),simplexify=true)
results = Example001.main(params)
@test results[:eh1] < tol
@test results[:el2] < tol

params = Dict{Symbol,Any}()
params[:mesh] = gk.cartesian_mesh((0,3,0,2,0,1),(5,5,5),simplexify=true)
results = Example001.main(params)
@test results[:eh1] < tol
@test results[:el2] < tol

params = Dict{Symbol,Any}()
params[:hi] = 1
results = Example001.main(params)
@test results[:eh1] < tol
@test results[:el2] < tol

domain = (0,1,0,1)
cells = (10,10)
fine_mesh = gk.cartesian_mesh(domain,cells)
domain = (0,30,0,10)
cells = (2,2)
coarse_mesh = gk.cartesian_mesh(domain,cells)
mesh, = gk.two_level_mesh(coarse_mesh,fine_mesh)

params = Dict{Symbol,Any}()
params[:mesh] = mesh
results = Example001.main(params)
@test results[:eh1] < tol
@test results[:el2] < tol

params = Dict{Symbol,Any}()
params[:mesh] = gk.cartesian_mesh((0,3,0,2,0,1),(80,80,80))
params[:solver] = Example001.cg_amg_solver(;verbose=true)
params[:export_vtu] = false
Example001.main(params)
Example001.main(params)

# Serial tests of periodic geometries 
outdir = joinpath(@__DIR__, "..", "output")
assetsdir = joinpath(@__DIR__, "..", "..", "..", "assets")

# Load periodic fine (unit cell) puzzle piece mesh with triangular refcells 
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

end # module

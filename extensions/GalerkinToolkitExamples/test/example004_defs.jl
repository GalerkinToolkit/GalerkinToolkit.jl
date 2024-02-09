
import GalerkinToolkit as gk
using GalerkinToolkitExamples: Example004
using Test
using PartitionedArrays
using PetscCall

function example004_tests_np_4(distribute)
    tol = 1.0e-8
    params = Dict{Symbol,Any}()
    domain = (0,10,0,10)
    cells_per_dir = (20,20)
    parts_per_dir = (2,2)
    np = prod(parts_per_dir)
    parts = distribute(LinearIndices((np,)))
    ghost_layers = 0
    mesh = gk.cartesian_mesh(domain,cells_per_dir,parts_per_dir;parts,ghost_layers)
    params[:mesh] = mesh
    results = Example004.main(params)
    @test results[:eh1] < tol
    @test results[:el2] < tol

    params = Dict{Symbol,Any}()
    domain = (0,10,0,10,0,10)
    cells_per_dir = (20,20,20)
    parts_per_dir = (2,2,2)
    np = prod(parts_per_dir)
    parts = distribute(LinearIndices((np,)))
    ghost_layers = 0
    mesh = gk.cartesian_mesh(domain,cells_per_dir,parts_per_dir;parts,ghost_layers)
    params[:mesh] = mesh
    results = Example004.main(params)
    @test results[:eh1] < tol
    @test results[:el2] < tol

end

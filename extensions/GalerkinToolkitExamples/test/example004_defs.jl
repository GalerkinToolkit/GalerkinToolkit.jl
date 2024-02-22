
import GalerkinToolkit as gk
using GalerkinToolkitExamples: Example004, Example002
using Test
using PartitionedArrays
using PetscCall
using TimerOutputs

function example004_tests_np_4(distribute)
    tol = 1.0e-8
    params = Dict{Symbol,Any}()
    domain = (0,10,0,10)
    cells_per_dir = (20,20)
    parts_per_dir = (2,2)
    np = prod(parts_per_dir)
    parts = distribute(LinearIndices((np,)))
    partition_strategy = gk.partition_strategy(;graph_nodes=:cells,graph_edges=:nodes,ghost_layers=0)
    mesh = gk.cartesian_mesh(domain,cells_per_dir;parts_per_dir,parts,partition_strategy)
    params[:mesh] = mesh
    results = Example004.main(params)
    @test results[:eh1] < tol
    @test results[:el2] < tol

    params = Dict{Symbol,Any}()
    domain = (0,10,0,10,0,10)
    cells_per_dir = (10,10,10)
    parts_per_dir = (2,2,1)
    np = prod(parts_per_dir)
    parts = distribute(LinearIndices((np,)))
    mesh = gk.cartesian_mesh(domain,cells_per_dir;parts_per_dir,parts,partition_strategy)
    params[:mesh] = mesh
    results = Example004.main(params)
    @test results[:eh1] < tol
    @test results[:el2] < tol


    params = Dict{Symbol,Any}()
    domain = (0,10,0,10,0,10)
    cells_per_dir = (50,50,50)
    parts_per_dir = (2,2,1)
    np = prod(parts_per_dir)
    parts = distribute(LinearIndices((np,)))
    mesh = gk.cartesian_mesh(domain,cells_per_dir;parts_per_dir,parts,partition_strategy)
    params[:mesh] = mesh
    params[:export_vtu] = false
    options = "-pc_type gamg -ksp_type cg -ksp_error_if_not_converged true -ksp_converged_reason -ksp_rtol 1.0e-6"
    PetscCall.init(args=split(options))
    # TODO a more comfortable way of passing a timer?
    linear_solver = Example002.ksp_solver()
    timer = TimerOutput()
    params[:solver] = Example004.nlsolve_solver(;timer,linear_solver,method=:newton,show_trace=true)
    params[:timer] = timer
    results = Example004.main(params)
    @test results[:eh1] < tol
    @test results[:el2] < tol

end

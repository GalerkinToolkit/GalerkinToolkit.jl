
import GalerkinToolkit as gk
using GalerkinToolkitExamples: Example002
using Test
using PartitionedArrays
using PetscCall
using Metis
using WriteVTK

function example002_tests_np_4(distribute)
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
    results = Example002.main(params)
    results = Example002.main(params)
    @test results[:eh1] < tol
    @test results[:el2] < tol

    # This is how to repartition a general unstructured grid
    pmesh = map_main(parts) do parts
        mesh = gk.cartesian_mesh(domain,cells_per_dir)
        graph = gk.mesh_graph(mesh;partition_strategy)
        graph_partition = Metis.partition(graph,np)
        gk.partition_mesh(mesh,np;partition_strategy,graph,graph_partition)
    end |> gk.scatter_mesh
    params[:mesh] = pmesh
    results = Example002.main(params)
    @test results[:eh1] < tol
    @test results[:el2] < tol

    partition_strategy = gk.partition_strategy(;graph_nodes=:nodes,graph_edges=:cells,ghost_layers=1)
    pmesh = map_main(parts) do parts
        mesh = gk.cartesian_mesh(domain,(3,3))
        graph = gk.mesh_graph(mesh;partition_strategy)
        graph_partition = Metis.partition(graph,np)
        gk.partition_mesh(mesh,np;partition_strategy,graph,graph_partition)
    end |> gk.scatter_mesh
    params[:mesh] = pmesh
    results = Example002.main(params)
    @test results[:eh1] < tol
    @test results[:el2] < tol

    options = "-pc_type gamg -ksp_type cg -ksp_error_if_not_converged true -ksp_converged_reason -ksp_rtol 1.0e-6"
    PetscCall.init(args=split(options))
    params = Dict{Symbol,Any}()
    domain = (0,10,0,10,0,10)
    cells_per_dir = (40,40,40)
    parts_per_dir = (2,2,1)
    np = prod(parts_per_dir)
    parts = distribute(LinearIndices((np,)))
    partition_strategy = gk.partition_strategy(;graph_nodes=:cells,graph_edges=:nodes,ghost_layers=0)
    mesh = gk.cartesian_mesh(domain,cells_per_dir;parts_per_dir,parts,partition_strategy)
    params[:mesh] = mesh
    params[:export_vtu] = false
    params[:solver] = Example002.ksp_solver()
    results = Example002.main(params)
    results = Example002.main(params)
end

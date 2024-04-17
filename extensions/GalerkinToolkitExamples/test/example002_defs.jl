
import GalerkinToolkit as gk
using GalerkinToolkitExamples: Example002
using Test
using PartitionedArrays
using PetscCall
using Metis
using WriteVTK

# NOTE: Global initialization prevents errors associated with initialization/finalization in 
# function bodies 
options = "-pc_type gamg -ksp_type cg -ksp_error_if_not_converged true -ksp_converged_reason -ksp_rtol 1.0e-6"
PetscCall.init(args=split(options))

function example002_tests_np_4(distribute)
    tol = 1.0e-8
    # options = "-pc_type gamg -ksp_type cg -ksp_error_if_not_converged true -ksp_converged_reason -ksp_rtol 1.0e-6"
    # PetscCall.init(args=split(options))

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
    results = Example002.main(params) # TODO: why call twice?
    results = Example002.main(params) # no need check convergence due to petsc args

    # PetscCall.finalize()
end

function example002_advanced_tests_np_4(distribute)
    tol = 1.0e-8
    # options = "-pc_type gamg -ksp_type cg -ksp_error_if_not_converged true -ksp_converged_reason -ksp_rtol 1.0e-6"
    # PetscCall.init(args=split(options))

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

    ## Periodic 2D square tests 
    # using default solver 
    params = Dict{Symbol,Any}()

    # load unit cell 
    unit_cell_mesh_fpath = joinpath(
        @__DIR__,
        "..", # projectdir?
        "..",
        "..",
        "assets",
        "unit_cell_2D_periodic_square_geometry_triangular_refcell.msh")
    unit_cell_mesh = gk.mesh_from_gmsh(unit_cell_mesh_fpath)

    final_pmesh, _ = gk.two_level_mesh(coarse_pmesh, unit_cell_mesh)

    params[:mesh] = final_pmesh 
    results = Example002.main(params)
    @test results[:eh1] < tol
    @test results[:el2] < tol

    # using kspsolver
    params[:export_vtu] = false
    params[:solver] = Example002.ksp_solver()
    results = Example002.main(params)
    results = Example002.main(params)

    ## Periodic 2D puzzle tests 
    ## TODO: This test fails 
    # using default solver 
    params = Dict{Symbol,Any}()

    # load unit cell 
    unit_cell_mesh_fpath = joinpath(
        @__DIR__,
        "..", # projectdir?
        "..",
        "..",
        "assets",
        "unit_cell_2D_periodic_puzzlepiece_geometry_triangular_refcell.msh")
    unit_cell_mesh = gk.mesh_from_gmsh(unit_cell_mesh_fpath)

    final_pmesh, _ = gk.two_level_mesh(coarse_pmesh, unit_cell_mesh)
    final_pmesh = gk.label_boundary_faces!(final_pmesh; physical_name="boundary")

    params[:mesh] = final_pmesh 
    params[:dirichlet_tags] = ["boundary"]
    params[:example_path] = joinpath("..", "output")
    params[:export_vtu] = true 
    results = Example002.main(params)
    @test results[:eh1] < tol 
    @test results[:el2] < tol

    # # using kspsolver
    # params[:export_vtu] = false
    # params[:solver] = Example002.ksp_solver()
    # results = Example002.main(params)
    # results = Example002.main(params)

    # ## Coarse 3D pmesh
    # domain = (0, 10, 0, 10, 0, 10)
    # cells = (4, 4, 4)
    # parts_per_dir = (2, 2, 2)
    # nparts = prod(parts_per_dir)
    # parts = DebugArray(LinearIndices((nparts,)))
    # coarse_pmesh = gk.cartesian_mesh(
    #     domain, cells;
    #     parts_per_dir, parts,
    #     partition_strategy=gk.partition_strategy(; ghost_layers=0))

    # ## Periodic 3D box 
    # params = Dict{Symbol,Any}()

    # unit_cell_mesh_fpath = joinpath(
    #     @__DIR__,
    #     "..", 
    #     "..",
    #     "..",
    #     "assets",
    #     "unit_cell_3D_periodic_box_geometry_triangular_refcell.msh")
    # unit_cell_mesh = gk.mesh_from_gmsh(unit_cell_mesh_fpath)

    # final_pmesh, _ = gk.two_level_mesh(coarse_pmesh, unit_cell_mesh)

    # params[:mesh] = final_pmesh 
    # results = Example002.main(params)
    # @test results[:eh1] < tol
    # @test results[:el2] < tol

    ## Periodic 3D puzzlepiece

    # PetscCall.finalize()
end 

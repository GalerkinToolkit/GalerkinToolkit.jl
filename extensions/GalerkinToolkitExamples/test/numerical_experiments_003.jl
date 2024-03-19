module numerical_experiments_003

import GalerkinToolkit as gk
using GalerkinToolkitExamples: Example003, Example001
using Test
using TimerOutputs
using ProfileView
using DataFrames
using Serialization

######## OPTIONS THAT CAN BE CONFIGURED ########
# 1. Parallelized loops: options on parallel are cell, quad, elem_j, elem_ij and full
# 2. Memory layout: for coalesced access options are cell_major and dof_major
# 3. Jacobian function: original, cpu_extension and gpu_extension.
# 4. autodiff to use: hand, flux or energy.
# 5. Precision of the jacobian: double, single or half.

# EXAMPLES
# params[:parallelization_level] = :cell
# params[:mem_layout] = :cell_major
# params[:jacobian] = :cpu_extension
# params[:autodiff] = :hand
# params[:precision] = Dict(:Float => Float16, :Int => Int16)

# Each row is one experiment configuration.
double = Dict(:Float => Float64, :Int => Int64)
single = Dict(:Float => Float32, :Int => Int32)
half = Dict(:Float => Float16, :Int => Int16)

# Define the parameter settings for each Jacobian function and get the combinations. 
autodiff = [:hand,:flux,:energy]
mesh = [gk.cartesian_mesh((0,3,0,2,0,1),(15,15,15))]
jacobian_orig = [:original]
prl_level_orig = [:cell]
mem_layout_orig = [:cell_major]

jacobian_cpu = [:cpu_extension]
prl_level_cpu = [:cell]
mem_layout_cpu = [:cell_major]

jacobian_gpu = [:gpu_extension]
prl_level_gpu = [:cell, :coalesce]
mem_layout_gpu = [:cell_major, :dof_major]

# Get all combinations
combinations_orig = collect(Iterators.product(autodiff, mesh, jacobian_orig, prl_level_orig, mem_layout_orig))
combinations_cpu = collect(Iterators.product(autodiff, mesh, jacobian_cpu, prl_level_cpu, mem_layout_cpu))
#combinations_gpu = collect(Iterators.product(autodiff, mesh, jacobian_gpu, prl_level_gpu, mem_layout_gpu))
#combinations_gpu = filter(comb -> !(comb[4] == :coalesce && comb[5] == :cell_major), combinations_gpu)

#param_combinations = vcat(DataFrame(combinations_orig), DataFrame(combinations_cpu), DataFrame(combinations_gpu))
param_combinations = vcat(DataFrame(combinations_orig), DataFrame(combinations_cpu))
rename!(param_combinations, :1 => :Autodiff, :2 => :Mesh, :3 => :Jacobian, :4 => :Prl_level, :5 => :Mem_layout)
param_combinations.Precision .= :double

# Set the different experiment parameters and variables.
params = Dict{Symbol,Any}()
meshes = [gk.cartesian_mesh((0,3,0,2,0,1),(15,15,15))]
params[:export_vtu] = false

tol = 1.0e-10
timer = TimerOutput()

simulations = 3
ns_to_s = 1e+9

# Set up the data structures to collect results and timer outputs.
result_df = DataFrame(Experiment = Int[], SimID = Int[], NewtonIterations = Int[], 
        TimeJacobian = Float64[], TimeSolve = Float64[], TimeSetup = Float64[], nCalls = Int[],
        nCells = Int[], Precision = DataType[], AutoDiff = [], Jacobian = [])

# Go for the loop of the experiments through the parameter sets.
for (index, experiment) in enumerate(eachrow(param_combinations))
    # Set the parameters
    params[:mesh] = experiment.Mesh
    params[:autodiff] = experiment.Autodiff
    params[:jacobian] = experiment.Jacobian
    params[:parallelization_level] = experiment.Prl_level
    params[:mem_layout] = experiment.Mem_layout
    params[:precision] = double

    for sim in 1:simulations

        result, time = Example003.main(params)

        result[:time_jacobian] = TimerOutputs.todict(time)["inner_timers"]["solve_problem"]["inner_timers"]["solver.solve!"]["inner_timers"]["jacobian_cells!"]["time_ns"]
        result[:time_solve] = TimerOutputs.todict(time)["inner_timers"]["solve_problem"]["time_ns"]
        result[:time_setup] = TimerOutputs.todict(time)["inner_timers"]["setup"]["time_ns"]
        result[:ncalls] = TimerOutputs.todict(time)["inner_timers"]["solve_problem"]["inner_timers"]["solver.solve!"]["inner_timers"]["jacobian_cells!"]["n_calls"]

        sim_result = (index, sim, result[:iterations], result[:time_jacobian]/ns_to_s, result[:time_solve]/ns_to_s, 
            result[:time_setup]/ns_to_s, result[:ncalls], result[:ncells], params[:precision][:Float], params[:autodiff], 
            params[:jacobian])
        push!(result_df, sim_result)

        reset_timer!(timer)
    end

end
println(result_df)

# # Save the results to the folder for further analysis. 
# filename = "experimental_result.dat"

# # Serialize the DataFrame
# open(joinpath(@__DIR__, filename), "w") do io
#     serialize(io, result_df)
# end


end
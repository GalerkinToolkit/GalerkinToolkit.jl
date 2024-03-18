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

# Each column is one experiment configuration.
double = Dict(:Float => Float64, :Int => Int64)
single = Dict(:Float => Float32, :Int => Int32)
half = Dict(:Float => Float16, :Int => Int16)

configurations_dict = Dict(
    :jacobian => [:original, :original, :original, :cpu_extension, :cpu_extension, :cpu_extension],
    :autodiff => [:hand, :flux, :energy, :hand, :flux, :energy],
    :mem_layout => [:cell_major, :cell_major, :cell_major, :dof_major, :dof_major, :dof_major],
    :parallelization_level => [:cell, :cell, :cell, :coalesce, :coalesce, :coalesce],
    :precision => [double, double, double, single, single, single]
)

# Set the different experiment parameters and variables.
params = Dict{Symbol,Any}()
meshes = [gk.cartesian_mesh((0,3,0,2,0,1),(15,15,15))]
params[:export_vtu] = false

tol = 1.0e-10
timer = TimerOutput()
experiments = length(configurations_dict[:jacobian])
simulations = 3
ns_to_s = 1e+9

# Set up the data structures to collect results and timer outputs.
result_df = DataFrame(Experiment = Int[], SimID = Int[], NewtonIterations = Int[], 
        TimeJacobian = Float64[], TimeSolve = Float64[], TimeSetup = Float64[], nCalls = Int[],
        nCells = Int[], Precision = DataType[], AutoDiff = [], Jacobian = [])

# Go for the loop of the experiments through the parameter sets.
for mesh in meshes
    params[:mesh] = mesh
    for e in 1:experiments
        # Set the parameters
        params[:autodiff] = configurations_dict[:autodiff][e]
        params[:jacobian] = configurations_dict[:jacobian][e]
        params[:parallelization_level] = configurations_dict[:parallelization_level][e]
        params[:mem_layout] = configurations_dict[:mem_layout][e]
        params[:precision] = configurations_dict[:precision][e]

        for sim in 1:simulations

            result, time = Example003.main(params)

            result[:time_jacobian] = TimerOutputs.todict(time)["inner_timers"]["solve_problem"]["inner_timers"]["solver.solve!"]["inner_timers"]["jacobian_cells!"]["time_ns"]
            result[:time_solve] = TimerOutputs.todict(time)["inner_timers"]["solve_problem"]["time_ns"]
            result[:time_setup] = TimerOutputs.todict(time)["inner_timers"]["setup"]["time_ns"]
            result[:ncalls] = TimerOutputs.todict(time)["inner_timers"]["solve_problem"]["inner_timers"]["solver.solve!"]["inner_timers"]["jacobian_cells!"]["n_calls"]

            sim_result = (e, sim, result[:iterations], result[:time_jacobian]/ns_to_s, result[:time_solve]/ns_to_s, 
                result[:time_setup]/ns_to_s, result[:ncalls], result[:ncells], params[:precision][:Float], params[:autodiff], 
                params[:jacobian])
            push!(result_df, sim_result)

            reset_timer!(timer)
        end
    end
end
println(result_df)

# Save the results to the folder for further analysis. 
filename = "experimental_result.dat"

# Serialize the DataFrame
open(joinpath(@__DIR__, filename), "w") do io
    serialize(io, result_df)
end


end
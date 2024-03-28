module test_functions

import GalerkinToolkit as gk
using GalerkinToolkitExamples: Example003, Example001
using Test
using TimerOutputs
#using ProfileView


tol = 1.0e-10
timer = TimerOutput()
params = Dict{Symbol,Any}()

#params[:mesh] = gk.cartesian_mesh((0,3,0,2,0,1),(5,5,5)) # 64000 cells. Out of memory Error.
#params[:mesh] = gk.cartesian_mesh((0,3,0,2,0,1),(5,5,5)) # Around 2.5 seconds per Jacobian. 27000 cells
#params[:mesh] = gk.cartesian_mesh((0,10,0,10),(2,2))
params[:mesh] = gk.cartesian_mesh((0,3,0,2,0,1),(20,20,20))
# (80,80,80) is about 18 seconds for jacobian with cpu_v1 500k cells | flux
# (60,60,60) is about 7 seconds for jacobian with cpu_v1 | flux

# First test original correct solution to compare the rest to
params[:autodiff] = :flux
params[:jacobian_implementation] = :original
params[:export_vtu] = false
ns_to_s = 1e+9
params[:float_type] = Dict(:Float => Float32, :Int => Int32)
#println("Jacobian original ")

# You want to loop this one to get several measurements. 
# results1, time  = Example003.main(params)
# print_timer(time)
# #From here you pull out the time measurement and the number of Newton iterations for convergence and store it in a df.

# # Test the kernel generic of the cpu extension
params[:jacobian_implementation] = :cpu_v1

println("Jacobian extension ", params[:jacobian_implementation])
reset_timer!(timer)
results2, time  = Example003.main(params)
print_timer(time)
# iters_extension = results2[:iterations]
jacobian = TimerOutputs.todict(time)["inner_timers"]["solve_problem"]["inner_timers"]["solver.solve!"]["inner_timers"]["jacobian_cells!"]["time_ns"]/ns_to_s
println("Time for jacobian ", params[:jacobian_implementation], " is ",jacobian)
# println("Number of cells: ", results2[:ncells])



# To test coalesce
params[:jacobian_implementation] = :cpu_v3
println("Jacobian extension ", params[:jacobian_implementation])

reset_timer!(timer)
results3, time = Example003.main(params)
print_timer(time)
jacobian = TimerOutputs.todict(time)["inner_timers"]["solve_problem"]["inner_timers"]["solver.solve!"]["inner_timers"]["jacobian_cells!"]["time_ns"]/ns_to_s
println("Time for jacobian ", params[:jacobian_implementation], " is ",jacobian)


###### To profile the extension #########
# println("The profiling part")
# reset_timer!(timer)
# ProfileView.@profview Example003.main(params)  # run once to trigger compilation (ignore this one)
# ProfileView.closeall()
# reset_timer!(timer)
# ProfileView.@profview Example003.main(params)

# gpu testing
# params[:jacobian_implementation] = :gpu_v1
# println("Jacobian extension ", params[:jacobian_implementation])

# reset_timer!(timer)
# results4, time = Example003.main(params)
# jacobian = TimerOutputs.todict(time)["inner_timers"]["solve_problem"]["inner_timers"]["solver.solve!"]["inner_timers"]["jacobian_cells!"]["time_ns"]/ns_to_s
# println("Time for jacobian ", params[:jacobian_implementation], " is ",jacobian)
# println("ncells is ", results4[:ncells])
# print_timer(time)


# params[:jacobian_implementation] = :gpu_v4
# println("Jacobian extension ", params[:jacobian_implementation])
# reset_timer!(timer)
# results4, time = Example003.main(params)
# print_timer(time)
# jacobian = TimerOutputs.todict(time)["inner_timers"]["solve_problem"]["inner_timers"]["solver.solve!"]["inner_timers"]["jacobian_cells!"]["time_ns"]/ns_to_s
# println("Time for jacobian ", params[:jacobian_implementation], " is ",jacobian)





end # module

module test_functions

import GalerkinToolkit as gk
using GalerkinToolkitExamples: Example003, Example001
using Test


tol = 1.0e-10
params = Dict{Symbol,Any}()

#params[:mesh] = gk.cartesian_mesh((0,3,0,2,0,1),(5,5,5)) # 64000 cells. Out of memory Error.
#params[:mesh] = gk.cartesian_mesh((0,3,0,2,0,1),(5,5,5)) # Around 2.5 seconds per Jacobian. 27000 cells
#params[:mesh] = gk.cartesian_mesh((0,10,0,10),(2,2))
params[:mesh] = gk.cartesian_mesh((0,3,0,2,0,1),(60,60,60))
# (80,80,80) is about 18 seconds for jacobian with cpu_v1 500k cells | flux
# (60,60,60) is about 7 seconds for jacobian with cpu_v1 | flux

# First test original correct solution to compare the rest to
params[:autodiff] = :flux
params[:export_vtu] = false
ns_to_s = 1e+9
params[:float_type] = Dict(:Float => Float64, :Int => Int32)
#println("Jacobian original ")

# You want to loop this one to get several measurements. 
# results, jacobian_time, gpu_transfer_time, gpu_setup_time = Example003.main(params)
# #From here you pull out the time measurement and the number of Newton iterations for convergence and store it in a df.

# # Test the kernel generic of the cpu extension
params[:jacobian_implementation] = :cpu_v1

println("Jacobian extension ", params[:jacobian_implementation])
t = @elapsed results, jacobian_time, gpu_transfer_time, gpu_setup_time  = Example003.main(params)
println(t)
println(jacobian_time, " ",gpu_transfer_time, " ",gpu_setup_time)
# println("Number of cells: ", results2[:ncells])



# # To test coalesce
# params[:jacobian_implementation] = :gpu_v1
# println("Jacobian extension ", params[:jacobian_implementation])

# results, jacobian_time, gpu_transfer_time, gpu_setup_time = Example003.main(params)


###### To profile the extension #########
# println("The profiling part")
# ProfileView.@profview Example003.main(params)  # run once to trigger compilation (ignore this one)
# ProfileView.closeall()
# ProfileView.@profview Example003.main(params)

# gpu testing
# params[:jacobian_implementation] = :gpu_v2
# println("Jacobian extension ", params[:jacobian_implementation])

# reset_timer!(timer)
# results, jacobian_time, gpu_transfer_time, gpu_setup_time = Example003.main(params)

# println("ncells is ", results4[:ncells])


# params[:jacobian_implementation] = :gpu_v3
# println("Jacobian extension ", params[:jacobian_implementation])
# results, jacobian_time, gpu_transfer_time, gpu_setup_time = Example003.main(params)




end # module

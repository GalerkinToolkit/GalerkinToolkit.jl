module test_functions

import GalerkinToolkit as gk
using GalerkinToolkitExamples: Example003, Example001
using Test
using TimerOutputs
using ProfileView

######## OPTIONS THAT CAN BE CONFIGURED ########
# 1. Jacobian function: original, cpu_extension and gpu_extension.
# 2. autodiff to use: hand, flux or energy.
# 3. Precision -> double, single or half

# EXAMPLES
# params[:jacobian] = :cpu_v1
# params[:autodiff] = :hand
# params[:precision] = :double

tol = 1.0e-10
timer = TimerOutput()

# Initialize the mesh to be used
params = Dict{Symbol,Any}()
# 
#params[:mesh] = gk.cartesian_mesh((0,3,0,2,0,1),(40,40,40)) # 64000 cells. Out of memory Error.
#params[:mesh] = gk.cartesian_mesh((0,3,0,2,0,1),(30,30,30)) # Around 2.5 seconds per Jacobian. 27000 cells
#params[:mesh] = gk.cartesian_mesh((0,10,0,10),(2,2))
params[:mesh] = gk.cartesian_mesh((0,3,0,2,0,1),(10,10,10))

# First test original correct solution to compare the rest to
params[:autodiff] = :flux
params[:jacobian_implementation] = :original
params[:export_vtu] = false
params[:precision] = Dict(:Float => Float64, :Int => Int32)
println("Jacobian original ")

# You want to loop this one to get several measurements. 
results1, time  = Example003.main(params)
print_timer(time,allocations=false)
# From here you pull out the time measurement and the number of Newton iterations for convergence and store it in a df.

iters_original = results1[:iterations]

@test results1[:eh1] < tol
@test results1[:el2] < tol

# Test the kernel generic of the cpu extension
params[:jacobian_implementation] = :cpu_v1

println("Jacobian extension ", params[:jacobian_implementation])

reset_timer!(timer)
results2, time  = Example003.main(params)
print_timer(time,allocations=false)
iters_extension = results2[:iterations]

println("Number of cells: ", results2[:ncells])
@test results2[:eh1] < tol
@test results2[:el2] < tol
@test iters_original == iters_extension


# To test coalesce
params[:jacobian_implementation] = :cpu_v2
println("Jacobian extension ", params[:jacobian_implementation])

reset_timer!(timer)
results3, time = Example003.main(params)
print_timer(time,allocations=false)
iters_coalesce = results2[:iterations]

@test results3[:eh1] < tol
@test results3[:el2] < tol
@test iters_original == iters_coalesce

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
# print_timer(time,allocations=false)
# iters_gpu_v1 = results3[:iterations]

# @test results4[:eh1] < tol
# @test results4[:el2] < tol
# @test iters_original == iters_gpu_v1

# params[:jacobian_implementation] = :gpu_v2
# println("Jacobian extension ", params[:jacobian_implementation])

# reset_timer!(timer)
# results4, time = Example003.main(params)
# print_timer(time,allocations=false)
# iters_coalesce = results4[:iterations]

# @test results4[:eh1] < tol
# @test results4[:el2] < tol
# @test iters_original == iters_gpu_v2

end # module

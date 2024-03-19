module test_functions

import GalerkinToolkit as gk
using GalerkinToolkitExamples: Example003, Example001
using Test
using TimerOutputs
using ProfileView

######## OPTIONS THAT CAN BE CONFIGURED ########
# 1. Parallelized loops: options on parallel are cell, quad, elem_j, elem_ij and full
# 2. Then you have the coalesced access kernel function coalesce
# 3. Memory layout: for coalesced access options are cell_major and dof_major
# 4. Jacobian function: original, cpu_extension and gpu_extension.
# 5. autodiff to use: hand, flux or energy.

# EXAMPLES
# params[:parallelization_level] = :cell
# params[:mem_layout] = :cell_major
# params[:jacobian] = :cpu_extension
# params[:autodiff] = :hand

tol = 1.0e-10
timer = TimerOutput()

# Initialize the mesh to be used
params = Dict{Symbol,Any}()
# 
#params[:mesh] = gk.cartesian_mesh((0,3,0,2,0,1),(40,40,40)) # 64000 cells. Out of memory Error.
#params[:mesh] = gk.cartesian_mesh((0,3,0,2,0,1),(30,30,30)) # Around 2.5 seconds per Jacobian. 27000 cells
#params[:mesh] = gk.cartesian_mesh((0,10,0,10),(2,2))
params[:mesh] = gk.cartesian_mesh((0,3,0,2,0,1),(5,5,5))

# First test original correct solution to compare the rest to
params[:autodiff] = :flux
params[:jacobian] = :original
params[:export_vtu] = false
params[:precision] = Dict(:Float => Float64, :Int => Int64)
println("Jacobian original ")

# You want to loop this one to get several measurements. 
results1, time  = Example003.main(params)
print_timer(time,allocations=false)
# From here you pull out the time measurement and the number of Newton iterations for convergence and store it in a df.

iters_original = results1[:iterations]

@test results1[:eh1] < tol
@test results1[:el2] < tol

# Test the kernel generic of the cpu extension
params[:jacobian] = :cpu_extension
params[:parallelization_level] = :split_kernel_nq
params[:mem_layout] = :cell_major

println("Jacobian extension ", params[:parallelization_level])

reset_timer!(timer)
results2, time  = Example003.main(params)
print_timer(time,allocations=false)
iters_extension = results2[:iterations]

println("Number of cells: ", results2[:ncells])
@test results2[:eh1] < tol
@test results2[:el2] < tol
@test iters_original == iters_extension


# # To test coalesce
# params[:parallelization_level] = :coalesce
# params[:mem_layout] = :dof_major
# println("Jacobian extension ", params[:parallelization_level])

# reset_timer!(timer)
# results2 = Example003.main(params)
# iters_coalesce = results2[:iterations]

# @test results2[:eh1] < tol
# @test results2[:el2] < tol
# @test iters_original == iters_coalesce

# ###### To profile the extension #########
# # println("The profiling part")
# # reset_timer!(timer)
# # ProfileView.@profview Example003.main(params)  # run once to trigger compilation (ignore this one)
# # ProfileView.closeall()
# # reset_timer!(timer)
# # ProfileView.@profview Example003.main(params)

# # Test the different parallelization levels for correctness.

# params[:parallelization_level] = :elem_j
# params[:mem_layout] = :cell_major
# println("Jacobian extension ", params[:parallelization_level])

# reset_timer!(timer)
# results3 = Example003.main(params)
# iters_element_j = results3[:iterations]

# @test results3[:eh1] < tol
# @test results3[:el2] < tol
# @test iters_original == iters_element_j


# # Fr the ij parallelisation
# params[:parallelization_level] = :elem_ij
# println("Jacobian extension ", params[:parallelization_level])

# reset_timer!(timer)
# results4 = Example003.main(params)
# iters_element_ij = results4[:iterations]

# @test results4[:eh1] < tol
# @test results4[:el2] < tol
# @test iters_original == iters_element_ij


# # For the quad 
# params[:parallelization_level] = :quad
# println("Jacobian extension ", params[:parallelization_level])

# reset_timer!(timer)
# results5 = Example003.main(params)
# iters_element_quad = results5[:iterations]

# @test results5[:eh1] < tol
# @test results5[:el2] < tol
# @test iters_original == iters_element_quad


# # For full parallelisation
# params[:parallelization_level] = :full
# println("Jacobian extension ", params[:parallelization_level])

# reset_timer!(timer)
# results6 = Example003.main(params)
# iters_element_full = results6[:iterations]

# @test results6[:eh1] < tol
# @test results6[:el2] < tol
# @test iters_original == iters_element_full



end # module

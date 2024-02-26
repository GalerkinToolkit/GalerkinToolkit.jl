module test_functions

import GalerkinToolkit as gk
using GalerkinToolkitExamples: Example003, Example001
using Test
using TimerOutputs
using ProfileView

# alternative 1: outer loop cells, inner loop dofs(quadrature points) ..extension!() + kernel generic
# alternative 2: parallelize the different loops.

tol = 1.0e-10

params = Dict{Symbol,Any}()
params[:mesh] = gk.cartesian_mesh((0,3,0,2,0,1),(30,30,30))
params[:autodiff] = :flux

# # Options are original, cpu_extension and gpu_extension
# params[:jacobian] = :original
# println("Jacobian original ")
# results = Example003.main(params)
# iters_original = results[:iterations]

# @test results[:eh1] < tol
# @test results[:el2] < tol

# ####### Second test #############
# params[:jacobian] = :cpu_extension
# params[:parallelization_level] = :cell
# println("Jacobian extension ", params[:parallelization_level])
# results = Example003.main(params)
# iters_extension = results[:iterations]

# @test results[:eh1] < tol
# @test results[:el2] < tol

# @test iters_original == iters_extension

# ###### To profile the extension #########
# ProfileView.@profview Example003.main(params)  # run once to trigger compilation (ignore this one)
# ProfileView.closeall()
# ProfileView.@profview Example003.main(params)

# ####### Loops parallelized OR cell/dof reordering ###########
# # options on parallel are cell, elem_j, elem_ij
# # params[:parallelization_level] = :elem_ij


params[:parallelization_level] = :reverse 
println("Jacobian extension ", params[:parallelization_level])
results = Example003.main(params)





end # module

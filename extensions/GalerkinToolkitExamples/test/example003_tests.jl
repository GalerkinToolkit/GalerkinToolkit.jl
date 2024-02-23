module Example003Tests

import GalerkinToolkit as gk
using GalerkinToolkitExamples: Example003, Example001
using Test
using TimerOutputs

tol = 1.0e-10

params = Dict{Symbol,Any}()
params[:mesh] = gk.cartesian_mesh((0,10,0,10),(2,2))
results = Example003.main(params)
@test results[:eh1] < tol
@test results[:el2] < tol

params[:autodiff] = :hand
results = Example003.main(params)
@test results[:eh1] < tol
@test results[:el2] < tol

params[:autodiff] = :flux
results = Example003.main(params)
@test results[:eh1] < tol
@test results[:el2] < tol

params[:autodiff] = :energy
results = Example003.main(params)
@test results[:eh1] < tol
@test results[:el2] < tol

params = Dict{Symbol,Any}()
params[:mesh] = gk.cartesian_mesh((0,3,0,2),(5,10),complexify=false)
params[:dirichlet_tags] = ["1-face-1","1-face-3","1-face-4"]
params[:neumann_tags] = ["1-face-2"]
params[:u] = (x) -> sum(x)
params[:f] = (x) -> 0.0
params[:g] = (x) -> 1.0
results = Example003.main(params)
@test results[:eh1] < tol
@test results[:el2] < tol
@test results[:ncells] == 10*5

mesh = gk.cartesian_mesh((0,3,0,2,0,1),(5,5,5))
timer = TimerOutput()
linear_solver = Example001.cg_amg_solver(;verbose=true)
solver = Example003.nlsolve_solver(;timer,linear_solver,show_trace=true,method=:newton)

params = Dict{Symbol,Any}()
params[:timer] = timer
params[:solver] = solver
params[:mesh] = mesh
params[:export_vtu] = false
params[:autodiff] = :hand
Example003.main(params)

reset_timer!(timer)
params[:autodiff] = :flux
Example003.main(params)

reset_timer!(timer)
params[:autodiff] = :energy
Example003.main(params)

mesh = gk.cartesian_mesh((0,3,0,2,0,1),(70,70,70))
reset_timer!(timer)
params[:mesh] = mesh

params[:autodiff] = :hand
Example003.main(params)

reset_timer!(timer)
params[:autodiff] = :flux
Example003.main(params)

reset_timer!(timer)
params[:autodiff] = :energy
Example003.main(params)


end # module

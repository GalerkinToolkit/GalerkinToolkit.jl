module Example003Tests

import GalerkinToolkit as gk
using GalerkinToolkitExamples: Example003
using Test

tol = 1.0e-10

params = Dict{Symbol,Any}()
params[:mesh] = gk.cartesian_mesh((0,1,0,1),(3,3))
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

end # module
module Example001Tests

import GalerkinToolkit as gk
using GalerkinToolkitExamples: Example001
using Test

tol = 1.0e-10

params = Dict{Symbol,Any}()
n = 10
params[:mesh] = gk.cartesian_mesh((0,1,0,1),(n,n))
results = Example001.main(params)
@test results[:error_h1_norm] < tol
@test results[:error_l2_norm] < tol

end # module


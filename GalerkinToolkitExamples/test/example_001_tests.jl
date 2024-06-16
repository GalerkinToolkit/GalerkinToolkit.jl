module Example001Tests

import GalerkinToolkit as gk
using GalerkinToolkitExamples: Example001
using Test

tol = 1.0e-10

n = 2
mesh = gk.cartesian_mesh((0,2,0,2),(n,n))
for dirichlet_method in (:strong,)#:nitsche)
    params = Dict{Symbol,Any}()
    params[:mesh] = mesh
    params[:dirichlet_tags] = ["1-face-1","1-face-3","1-face-4"]
    params[:neumann_tags] = ["1-face-2"]
    params[:dirichlet_method] = dirichlet_method
    results = Example001.main(params)
    @test results[:error_h1_norm] < tol
    @test results[:error_l2_norm] < tol
end

end # module


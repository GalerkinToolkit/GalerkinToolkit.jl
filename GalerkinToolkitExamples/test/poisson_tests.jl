module PoissonTests

import GalerkinToolkit as GT
using GalerkinToolkitExamples: Poisson
using Test

tol = 1.0e-10

# TODO do not use a large value of n here.
# The code now is VERY slow.
n = 2
mesh = GT.cartesian_mesh((0,2,0,2),(n,n))

for discretization_method in (:continuous_galerkin,)
    for dirichlet_method in (:strong,)
        params = Dict{Symbol,Any}()
        params[:implementation] = :hand_written
        params[:mesh] = mesh
        params[:dirichlet_tags] = ["1-face-1","1-face-2","1-face-3","1-face-4"]
        params[:discretization_method] = discretization_method
        params[:dirichlet_method] = dirichlet_method
        params[:interpolation_degree] = 2
        results = Poisson.main(params)
        @test results[:error_l2_norm] < tol
        @test results[:error_h1_norm] < tol
    end
end

for discretization_method in (:interior_penalty,:continuous_galerkin)
    for dirichlet_method in (:multipliers,:nitsche,:strong)
        params = Dict{Symbol,Any}()
        params[:implementation] = :automatic
        params[:mesh] = mesh
        params[:dirichlet_tags] = ["1-face-1","1-face-3","1-face-4"]
        params[:neumann_tags] = ["1-face-2"]
        params[:discretization_method] = discretization_method
        params[:dirichlet_method] = dirichlet_method
        params[:interpolation_degree] = 2
        results = Poisson.main(params)
        @test results[:error_l2_norm] < tol
        @test results[:error_h1_norm] < tol
    end
end

end # module


module PoissonTest

using Test
import GalerkinToolkit as gt
# TODO normal vectors (reference normals plus push)
# Triangular meshes

include("../examples/poisson.jl")

tol = 1.0e-10

params = Dict{Symbol,Any}()
params[:mesh] = gt.cartesian_mesh((0,3,0,2),(10,5),complexify=false)
params[:dirichlet_tags] = ["1-face-1","1-face-3","1-face-4"]
params[:neumann_tags] = ["1-face-2"]
params[:u] = (x) -> sum(x)
params[:f] = (x) -> 0.0
params[:g] = (x) -> 1.0
@time results = poisson(params)
@test results[:eh1] < tol
@test results[:el2] < tol
@test results[:ncells] == 10*5

params = Dict{Symbol,Any}()
params[:mesh] = gt.cartesian_mesh((0,3,0,2,0,1),(5,5,5))
@time results = poisson(params)
@test results[:eh1] < tol
@test results[:el2] < tol
@test results[:ncells] == 5*5*5

params = Dict{Symbol,Any}()
params[:mesh] = gt.cartesian_mesh((0,3,0,2),(5,5),simplexify=true)
@time results = poisson(params)
@test results[:eh1] < tol
@test results[:el2] < tol

params = Dict{Symbol,Any}()
params[:mesh] = gt.cartesian_mesh((0,3,0,2,0,1),(5,5,5),simplexify=true)
@time results = poisson(params)
@test results[:eh1] < tol
@test results[:el2] < tol

params = Dict{Symbol,Any}()
params[:hi] = 1
@time results = poisson(params)
@test results[:eh1] < tol
@test results[:el2] < tol


end # module

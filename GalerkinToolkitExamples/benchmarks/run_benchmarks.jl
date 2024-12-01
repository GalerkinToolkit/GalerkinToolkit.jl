module GalerkinToolkitBenchmarkTests

using BenchmarkTools

import GalerkinToolkit as GT
using GalerkinToolkitExamples: Poisson


function handwritten_poisson(n)
	mesh = GT.cartesian_mesh((0,2,0,2,0,2), (n,n,n))

	params = Dict{Symbol,Any}()
	params[:implementation] = :hand_written
	params[:mesh] = mesh
	params[:dirichlet_tags] = ["1-face-1","1-face-2","1-face-3","1-face-4"]
	params[:discretization_method] = :continuous_galerkin
	params[:dirichlet_method] = :strong
	params[:integration_degree] = 2

	Poisson.main(params)
end

function automatic_poisson(n)
	mesh = GT.cartesian_mesh((0,2,0,2,0,2), (n,n,n))

	params = Dict{Symbol,Any}()
	params[:implementation] = :automatic
	params[:mesh] = mesh
	params[:dirichlet_tags] = ["1-face-1","1-face-2","1-face-3","1-face-4"]
	params[:discretization_method] = :continuous_galerkin
	params[:dirichlet_method] = :strong
	params[:integration_degree] = 2

	Poisson.main(params)
end

suite = BenchmarkGroup()
suite["poisson-hand"] = BenchmarkGroup(["Poisson", "handwritten"])
suite["poisson-hand"]["n=10"] = @benchmarkable handwritten_poisson(10)

suite["poisson-auto"] = BenchmarkGroup(["Poisson", "automatic"])
suite["poisson-auto"]["n=10"] = @benchmarkable automatic_poisson(10)

tune!(suite)

results = run(suite, verbose = true)

BenchmarkTools.save("output.json", median(results))

end # module

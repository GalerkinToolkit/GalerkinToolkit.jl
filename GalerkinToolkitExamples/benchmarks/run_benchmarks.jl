module GalerkinToolkitBenchmarkTests

using BenchmarkTools

import GalerkinToolkit as GT
using GalerkinToolkitExamples: Poisson


function handwritten_poisson(n)
	mesh = GT.cartesian_mesh((0,2,0,2),(n,n))

	params = Dict{Symbol,Any}()
	params[:implementation] = :hand_written
	params[:mesh] = mesh
	params[:dirichlet_tags] = ["1-face-1","1-face-2","1-face-3","1-face-4"]
	params[:discretization_method] = :continuous_galerkin
	params[:dirichlet_method] = :strong
	params[:integration_degree] = 2

	Poisson.main(params)
end


suite = BenchmarkGroup()
suite["poisson"] = BenchmarkGroup(["tag1", "tag2"])
suite["poisson"][2] = @benchmarkable handwritten_poisson(2)
suite["poisson"][3] = @benchmarkable handwritten_poisson(3)

tune!(suite)
results = run(suite, verbose = true)

BenchmarkTools.save("output.json", median(results))

end # module

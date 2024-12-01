module GalerkinToolkitBenchmarkTests

using BenchmarkTools

import GalerkinToolkit as GT
using GalerkinToolkitExamples: Poisson


function handwritten_poisson(n)
	"""
	Runs the hand-written Poisson example code, for a 3D
	mesh of dimensions n x n x n.
	"""

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


# Build a benchmark suite for the Poisson example
suite = BenchmarkGroup()
suite["poisson-hand"] = BenchmarkGroup(["Poisson", "handwritten"])
suite["poisson-hand"]["n=10"] = @benchmarkable handwritten_poisson(10)

# Run all benchmarks
tune!(suite)
results = run(suite, verbose = true)

# Save benchmark results for tracking and visualization
BenchmarkTools.save("output.json", median(results))

end # module

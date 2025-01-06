module PoissonCGTests

using Test
using GalerkinToolkitExperiments: GTRunner, FenicsRunner,FenicsxRunner, GridapRunner

display("CG experiment start")
result = Dict()

# solvers = [GridapRunner]
solvers = [GTRunner, GridapRunner]
params = [(3, 1, true), ]


for solver in solvers
    for param in params
        output = solver.poisson_cg_assembly(param...) # TODO: use DrWatson to cache experiments
        display(output)
    end
end

# TODO: result sanity check
@test 1 == 1


# compare


display("CG experiment end")



end # module


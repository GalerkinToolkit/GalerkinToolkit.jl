module PoissonDGTests

using Test
using GalerkinToolkitExperiments: GTRunner, FenicsRunner,FenicsxRunner, GridapRunner

display("DG experiment start")
result = Dict()

# solvers = [GridapRunner]
solvers = [GTRunner]
params = [(3, 3, true), ]


for solver in solvers
    for param in params
        output = solver.poisson_dg_assembly(param...) # TODO: use DrWatson to cache experiments
        display(output)
    end
end

# TODO: result sanity check
@test 1 == 1


# compare


display("DG experiment end")



end # module


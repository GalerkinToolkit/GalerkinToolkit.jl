module StokesTests

using Test
using GalerkinToolkitExperiments: GTRunner, FenicsRunner,FenicsxRunner, GridapRunner

display("Stokes experiment start")
result = Dict()

# solvers = [GridapRunner]
solvers = [GTRunner]
params = [(8, 2, true), ]


for solver in solvers
    for param in params
        output = solver.stokes_assembly(param...) # TODO: use DrWatson to cache experiments
        display(output)
    end
end

# TODO: result sanity check
@test 1 == 1


# compare


display("Stokes experiment end")



end # module


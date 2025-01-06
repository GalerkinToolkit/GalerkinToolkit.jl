module ExperimentsTests

using Test

@testset "Experiments" begin
    # @testset "poisson_cg" begin include("poisson_cg.jl") end
    # @testset "poisson_dg" begin include("poisson_dg.jl") end
    @testset "stokes" begin include("stokes.jl") end
    # TODO: dg, stokes
end

end # module

module GalerkinToolkitExamplesTests

using Test

@testset "GalerkinToolkitExamples" begin
    @testset "poisson" begin
        include("poisson_tests.jl")
    end
end

end # module

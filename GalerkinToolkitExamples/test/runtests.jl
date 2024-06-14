module GalerkinToolkitExamplesTests

using Test

@testset "GalerkinToolkitExamples" begin
    @testset "example001" begin include("example_001_tests.jl") end
end

end # module

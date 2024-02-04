module GalerkinToolkitExamplesTests

using Test

@testset "GalerkinToolkitExamples" begin
    @testset "example001" begin include("example001_tests.jl") end
    @testset "example002" begin include("example002_tests.jl") end
end

end # module

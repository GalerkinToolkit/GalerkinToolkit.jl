module GalerkinToolkitExamplesTests

using Test

@testset "GalerkinToolkitExamples" begin
    @testset "example001" begin include("example001_tests.jl") end
    @testset "example002" begin include("example002_tests.jl") end
end

@test_broken include(joinpath("mpi_array","runtests.jl"))

end # module

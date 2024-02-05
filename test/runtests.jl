module GalerkinToolkitTests

using Test

@testset "GalerkinToolkit" begin
    @testset "Geometry" begin include("geometry_tests.jl") end
    @testset "Integration" begin include("integration_tests.jl") end
    @testset "Interpolation" begin include("interpolation_tests.jl") end
end

end # module

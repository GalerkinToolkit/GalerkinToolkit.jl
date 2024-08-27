module GalerkinToolkitTests

using Test

@testset "GalerkinToolkit" begin
    @testset "Geometry" begin include("geometry_tests.jl") end
    @testset "Symbolics" begin include("symbolics_tests.jl") end
    @testset "Domain" begin include("domain_tests.jl") end
    @testset "Integration" begin include("integration_tests.jl") end
    @testset "Interpolation" begin include("interpolation_tests.jl") end
    #@testset "Assembly" begin include("assembly_tests.jl") end
end

end # module

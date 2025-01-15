module GalerkinToolkitTests

using Test

@testset "GalerkinToolkit" begin
    @testset "Helpers" begin include("new_helpers_tests.jl") end
    @testset "Domain" begin include("new_domain_tests.jl") end
    @testset "Mesh" begin include("new_mesh_tests.jl") end
    @testset "Cartesian mesh" begin include("cartesian_mesh_tests.jl") end
    @testset "Gmsh" begin include("gmsh_tests.jl") end
    @testset "PMesh" begin include("p_mesh_tests.jl") end
    @testset "Topology" begin include("new_topology_tests.jl") end
    @testset "Symbolics" begin include("symbolics_tests.jl") end
    @testset "Quadrature" begin include("new_quadrature_tests.jl") end
    @testset "Quantity" begin include("quantity_tests.jl") end
    @testset "Field" begin include("field_tests.jl") end
    @testset "Integration" begin include("integration_tests.jl") end
    @testset "Space" begin include("new_space_tests.jl") end
    @testset "Visualization" begin include("visualization_tests.jl") end
    @testset "Assembly" begin include("assembly_tests.jl") end
end

end # module

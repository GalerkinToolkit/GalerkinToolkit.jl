module GalerkinToolkitTests

using Test

@testset "GalerkinToolkit" begin
    @testset "Helpers" begin include("helpers_tests.jl") end
    @testset "Domain" begin include("domain_tests.jl") end
    @testset "Mesh" begin include("mesh_tests.jl") end
    @testset "Cartesian mesh" begin include("cartesian_mesh_tests.jl") end
    @testset "Gmsh" begin include("gmsh_tests.jl") end
    @testset "PMesh" begin include("p_mesh_tests.jl") end
    @testset "Topology" begin include("topology_tests.jl") end
    @testset "Symbolics" begin include("symbolics_tests.jl") end
    @testset "Quadrature" begin include("quadrature_tests.jl") end
    #@testset "Quantity" begin include("quantity_tests.jl") end
    @testset "Field" begin include("field_tests.jl") end
    @testset "Integration" begin include("integration_tests.jl") end
    @testset "Space" begin include("space_tests.jl") end
    @testset "Accessors" begin include("accessors_tests.jl") end
    @testset "Visualization" begin include("visualization_tests.jl") end
    @testset "Compiler" begin include("compiler_tests.jl") end
    @testset "Problems" begin include("problems_tests.jl") end
    @testset "ProblemsExt" begin include("problems_ext_tests.jl") end
    @testset "Assembly" begin include("assembly_tests.jl") end
    #@testset "NewIR" begin include("new_ir_tests.jl") end
    @testset "Issue 224" begin include("issue_224.jl") end
    @testset "Issue 230" begin include("issue_230.jl") end
end

end # module

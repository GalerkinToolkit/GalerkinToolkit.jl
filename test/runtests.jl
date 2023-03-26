module GalerkinToolkitTests

using GalerkinToolkit
using Test

@testset "GalerkinToolkit" begin
    @time @testset "meshes" begin include("meshes_tests.jl") end
    @time @testset "write_vtk" begin include("write_vtk_tests.jl") end
    @time @testset "mesh_interface" begin include("mesh_interface_tests.jl") end
    @time @testset "gmsh" begin include("gmsh_tests.jl") end
    @time @testset "p4est" begin include("p4est_tests.jl") end
    @time @testset "integration" begin include("integration_tests.jl") end
end

end # module

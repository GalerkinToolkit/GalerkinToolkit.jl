module GalerkinToolkitTests

using GalerkinToolkit
using Test

@testset "GalerkinToolkit" begin
    @time @testset "meshes" begin include("test_meshes.jl") end
    @time @testset "write_vtk" begin include("test_write_vtk.jl") end
    @time @testset "mesh_interface" begin include("test_mesh_interface.jl") end
    @time @testset "gmsh" begin include("test_gmsh.jl") end
    @time @testset "p4est" begin include("test_p4est.jl") end
end

end # module

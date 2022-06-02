module GalerkinToolkitTests

using GalerkinToolkit
using Test

@time @testset "jagged_array" begin include("jagged_array_tests.jl") end

@time @testset "fe_mesh" begin include("fe_mesh_tests.jl") end

@time @testset "polytopal_complex" begin include("polytopal_complex_tests.jl") end

@time @testset "meshes" begin include("meshes_tests.jl") end

@time @testset "gmsh" begin include("gmsh_tests.jl") end

@time @testset "vtk" begin include("vtk_tests.jl") end

@time @testset "p4est" begin include("p4est_tests.jl") end

end # module

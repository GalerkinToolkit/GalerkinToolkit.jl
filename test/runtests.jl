module GalerkinToolkitTests

using GalerkinToolkit
using Test

@time @testset "jagged_array" begin include("jagged_array_tests.jl") end

@time @testset "fe_mesh" begin include("fe_mesh_tests.jl") end

@time @testset "meshes" begin include("meshes_tests.jl") end

end # module

module GalerkinToolkitTests

using GalerkinToolkit
using Test

@time @testset "JaggedArray" begin include("jagged_array_tests.jl") end

end # module

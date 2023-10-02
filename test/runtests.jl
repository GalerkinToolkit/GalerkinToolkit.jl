module GalerkinToolkitTests

using Test

@time @testset "unittests.jl" begin include("unittests.jl") end

@time @testset "poisson_test.jl" begin include("poisson_test.jl") end

end # module

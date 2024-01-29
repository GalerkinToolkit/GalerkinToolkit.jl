module GalerkinToolkitTests

using Test

@testset "Mesh interface" begin include("mesh_interface_tests.jl") end

@testset "Iso-parametric Poisson" begin include("iso_param_poisson_tests.jl") end

end # module

module GalerkinToolkit
const gk = GalerkinToolkit

using WriteVTK
using FastGaussQuadrature
using StaticArrays
using LinearAlgebra
using ForwardDiff
using Gmsh
using PartitionedArrays
using Combinatorics
using SparseArrays
using Metis

include("geometry.jl")
include("integration.jl")
include("interpolation.jl")
include("domain.jl")

end # module

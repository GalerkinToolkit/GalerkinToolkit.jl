module GalerkinToolkit

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

end # module

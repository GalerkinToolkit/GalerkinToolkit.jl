module GalerkinToolkit
const GT = GalerkinToolkit

using WriteVTK
using FastGaussQuadrature
using StaticArrays
using LinearAlgebra
import LinearAlgebra
import ForwardDiff
using Gmsh
using PartitionedArrays
using Combinatorics
using SparseArrays
using Metis

include("geometry.jl")
include("domain.jl")
include("integration.jl")
include("interpolation.jl")
include("assembly.jl")

end # module

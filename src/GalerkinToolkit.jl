module GalerkinToolkit
const gk = GalerkinToolkit

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
using Metatheory

include("geometry.jl")
include("symbolics.jl")
include("domain.jl")
include("integration.jl")
include("interpolation.jl")
include("assembly.jl")

end # module

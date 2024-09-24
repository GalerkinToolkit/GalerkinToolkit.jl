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
using BlockArrays
using FillArrays
using Combinatorics
using SparseArrays
using Metis
using Metatheory
import Makie

include("geometry.jl")
include("symbolics.jl")
include("domain.jl")
include("integration.jl")
include("interpolation.jl")
include("assembly.jl")
include("visualization.jl")

end # module

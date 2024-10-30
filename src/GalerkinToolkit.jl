module GalerkinToolkit
const GT = GalerkinToolkit

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

include("mesh.jl")
include("p_mesh.jl")
include("symbolics.jl")
include("domain.jl")
include("integration.jl")
include("interpolation.jl")
include("assembly.jl")
include("visualization.jl")

end # module

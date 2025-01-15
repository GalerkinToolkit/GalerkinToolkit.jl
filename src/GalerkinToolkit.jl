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
import MacroTools
import PartitionedSolvers as PS
import AutoHashEquals
import AbstractTrees

include("helpers.jl")
include("abstract_types.jl")
include("domain.jl")
include("mesh.jl")
include("cartesian_mesh.jl")
include("gmsh.jl")
include("p_mesh.jl")
include("topology.jl")
include("symbolics.jl")
include("field.jl")
include("integration.jl")
include("space.jl")
include("quadrature.jl")
include("quantity.jl")
include("visualization.jl")
include("assembly.jl")

#include("mesh.jl")
#include("p_mesh.jl")
#include("symbolics.jl")
#include("domain.jl")
#include("integration.jl")
#include("interpolation.jl")
#include("assembly.jl")
#include("visualization.jl")
#include("helpers.jl")

end # module

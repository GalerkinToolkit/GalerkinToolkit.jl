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

include("new_helpers.jl")
include("new_domain.jl")
include("new_mesh.jl")
include("new_topology.jl")
include("new_space.jl")
include("new_quadrature.jl")
include("cartesian_mesh.jl")
include("gmsh.jl")
include("p_mesh.jl")

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

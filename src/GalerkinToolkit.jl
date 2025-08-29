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
using AutoHashEquals: @auto_hash_equals
import MacroTools
import PartitionedSolvers as PS
import PartitionedSolvers: update
import AutoHashEquals
import AbstractTrees

import ForwardDiff: gradient, jacobian

include("helpers.jl")
include("abstract_types.jl")
include("domain.jl")
include("mesh.jl")
include("cartesian_mesh.jl")
include("gmsh.jl")
include("p_mesh.jl")
include("topology.jl")
#include("symbolics.jl")
include("field.jl")
#include("integration.jl")
include("space.jl")
include("quadrature.jl")
include("accessors.jl")
include("assembly.jl")
include("problems.jl")
include("compiler.jl")
#include("quantity.jl")
include("visualization.jl")
#include("new_ir.jl")

end # module

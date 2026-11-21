module GalerkinToolkitCUDAExt


import GalerkinToolkit as GT
import Adapt
using CUDA
using LinearAlgebra: Transpose
using SparseArrays: nonzeros

include("GalerkinToolkitCUDAExt/cuda_assembly.jl")
include("GalerkinToolkitCUDAExt/cuda_compiler_data.jl")
include("GalerkinToolkitCUDAExt/cuda_compiler.jl")
include("GalerkinToolkitCUDAExt/cuda_execution.jl")

end # modules

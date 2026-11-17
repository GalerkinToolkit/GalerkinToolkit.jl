module GalerkinToolkitCUDAExt


import GalerkinToolkit as GT
using CUDA
using SparseArrays: nonzeros

include("GalerkinToolkitCUDAExt/cuda_assembly.jl")
include("GalerkinToolkitCUDAExt/cuda_compiler_data.jl")

end # modules

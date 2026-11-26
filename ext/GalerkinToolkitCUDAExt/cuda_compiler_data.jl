# Prepared integration data use the backend-neutral Adapt rules in compiler_data.jl.
# Adapting to CuArray preserves the element-major layout, including StaticArrays
# stored inside arrays. Reference tables currently use ordinary device memory.
# CUDA.Const can mark device arrays read-only, but does not place their contents
# in CUDA's separate constant-memory address space.

# Literal node indices let the compiler scalarize short coordinate gathers.
GT.ast_prepare_jacobian_coordinates(ast,::GT.CUDAExecution) =
    GT.ast_gather_jacobian_coordinates(ast,Val(true))

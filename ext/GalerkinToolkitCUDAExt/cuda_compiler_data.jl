# Prepared integration data use the backend-neutral Adapt rules in compiler_data.jl.
# Adapting to CuArray preserves the element-major layout, including StaticArrays
# stored inside arrays. Reference tables currently use ordinary device memory.
# CUDA.Const can mark device arrays read-only, but does not place their contents
# in CUDA's separate constant-memory address space.

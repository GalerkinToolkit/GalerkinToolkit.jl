const CUDACompilerBackend = parentmodule(CUDA.cufunction)

struct CUDAAssemblyCompilerOptions
    opt_level::Int
    fastmath::Bool
end

# CUDA.jl forwards `name` to compiler_config and includes it in the kernel cache.
# Only our option type intercepts that keyword to set GPUCompiler's opt_level.
function Core.kwcall(kwargs::NamedTuple{(:name,),Tuple{CUDAAssemblyCompilerOptions}},
        ::typeof(CUDACompilerBackend.compiler_config),device::CUDA.CuDevice)
    options = kwargs.name
    config = CUDACompilerBackend.compiler_config(device;always_inline=true,fastmath=options.fastmath)
    CUDACompilerBackend.CompilerConfig(config;opt_level=options.opt_level)
end

function GT.validate_matrix_assembly(execution::GT.CUDAExecution)
    CUDA.functional(true) || throw(ArgumentError("CUDA matrix assembly requires a functional GPU"))
    hasfield(CUDACompilerBackend.CompilerConfig,:opt_level) || throw(ArgumentError(
        "CUDA matrix assembly requires a GPUCompiler version with opt_level support"))
    limit = CUDA.attribute(CUDA.device(),CUDA.DEVICE_ATTRIBUTE_MAX_THREADS_PER_BLOCK)
    execution.threads <= limit || throw(ArgumentError("threads exceeds the device's per-block limit"))
    nothing
end

function GT.prepare_matrix_allocation(::GT.CUDAExecution,allocation)
    GT.is_element_csc_allocation(allocation.allocation) || throw(ArgumentError(
        "CUDA matrix assembly currently requires element_csc_assembly()"))
    Adapt.adapt(CuArray,GT.to_layout(allocation,GT.ElementMajorLayout()))
end

struct CUDAAssemblyRunner{K,D,A}
    kernel::K
    data::D
    argument::A
    nelements::Int
    execution::GT.CUDAExecution
end

function GT.assembly_runner(execution::GT.CUDAExecution,prepared::CUDAAssemblyProgram,
        argument,nelements)
    CUDAAssemblyRunner(prepared.kernel,prepared.data,argument,Int(nelements),execution)
end

function GT.assembly_runner(::GT.CUDAExecution,prepare,argument,nelements)
    throw(ArgumentError("CUDA vector assembly is not implemented"))
end

struct CUDAAssemblyElementRange{I<:Integer}
    first::I
    stride::I
    last::I
end

GT.assembly_element_range(elements::CUDAAssemblyElementRange) = elements

@inline function Base.iterate(elements::CUDAAssemblyElementRange,element=elements.first)
    element <= elements.last || return nothing
    element, element + elements.stride
end

function cuda_assembly_kernel!(kernel,argument,data,parameters,nelements::I) where I<:Integer
    threads = I(blockDim().x)
    id = (I(blockIdx().x)-one(I))*threads + I(threadIdx().x)
    stride = threads*I(gridDim().x)
    if id <= nelements
        elements = CUDAAssemblyElementRange(id,stride,nelements)
        @inline kernel(argument,data,elements,parameters)
    end
    nothing
end

function (runner::CUDAAssemblyRunner)(parameters::Tuple)
    runner.nelements == 0 && return nothing
    threads = runner.execution.threads
    blocks = cld(runner.nelements,threads)
    # Keep even the final grid-stride increment representable in the index type.
    index_type = blocks <= div(typemax(Int32)-runner.nelements,threads) ? Int32 : Int
    options = CUDAAssemblyCompilerOptions(runner.execution.opt_level,runner.execution.fastmath)
    @cuda threads=threads blocks=blocks name=options cuda_assembly_kernel!(
        runner.kernel,runner.argument,runner.data,parameters,index_type(runner.nelements))
    nothing
end

function GT.validate_matrix_parameters(::GT.CUDAExecution,parameters)
    parameters isa Tuple || throw(ArgumentError("CUDA assembly parameters must be a tuple"))
    isbitstype(typeof(CUDA.cudaconvert(parameters))) || throw(ArgumentError(
        "CUDA assembly parameters must be scalars or explicitly adapted device-compatible data"))
    nothing
end

function GT.validate_matrix_update(execution::GT.CUDAExecution,parameters,cache)
    GT.validate_matrix_parameters(execution,parameters)
    parameter_count = cache.parameter_count
    if length(parameters) != parameter_count
        throw(ArgumentError("CUDA assembly expects $parameter_count parameters, got $(length(parameters))"))
    end
    nothing
end

# Reset, all contribution launches, and compression/copies use the calling task's
# stream. Complete once at the public assembly/update boundary, not per domain.
GT.finish_matrix_assembly!(::GT.CUDAExecution) = CUDA.synchronize(CUDA.stream())

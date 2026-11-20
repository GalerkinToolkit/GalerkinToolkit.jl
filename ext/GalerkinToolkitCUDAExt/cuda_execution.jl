function GT.validate_matrix_assembly(execution::GT.CUDAExecution)
    CUDA.functional(true) || throw(ArgumentError("CUDA matrix assembly requires a functional GPU"))
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

function cuda_assembly_kernel!(kernel,argument,data,parameters,nelements)
    id = Int((blockIdx().x-1)*blockDim().x+threadIdx().x)
    count = Int(blockDim().x*gridDim().x)
    if id <= nelements
        worker = GT.AssemblyWorker(id,count,nelements,Val(:cyclic))
        @inline kernel(argument,data,worker,parameters)
    end
    nothing
end

function (runner::CUDAAssemblyRunner)(parameters::Tuple)
    runner.nelements == 0 && return nothing
    @cuda threads=runner.execution.threads blocks=cld(runner.nelements,runner.execution.threads) cuda_assembly_kernel!(
        runner.kernel,runner.argument,runner.data,parameters,runner.nelements)
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

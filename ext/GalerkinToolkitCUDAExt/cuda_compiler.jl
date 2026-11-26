# Private, statically sized storage for each generated device invocation.
struct CUDAThreadLocalArray{T,N,S} <: DenseArray{T,N}
    ptr::Core.LLVMPtr{T,0}
end

Base.size(::CUDAThreadLocalArray{T,N,S}) where {T,N,S} = S
Base.IndexStyle(::Type{<:CUDAThreadLocalArray}) = IndexLinear()

@inline @generated function local_alignment(::Type{T}) where T
    :(Val{$(Base.datatype_alignment(T))}())
end

# The custom getindex and setindex! methods provide CUDA-compatible indexing for
# thread-local buffers. Plain MArray indexing can emit host pointer intrinsics.
@inline function Base.getindex(array::CUDAThreadLocalArray{T},index::Int) where T
    unsafe_load(array.ptr,index,local_alignment(T))
end

@inline function Base.setindex!(array::CUDAThreadLocalArray{T},value,index::Int) where T
    unsafe_store!(array.ptr,convert(T,value),index,local_alignment(T))
    array
end

cuda_scratch_shape(dims::Integer...) = Val(Tuple(Int.(dims)))

@inline @generated function local_scratch(::Type{T},::Val{S}) where {T,S}
    isbitstype(T) || return :(throw(ArgumentError("CUDA assembly scratch requires an isbits element type")))
    bytes = Base.aligned_sizeof(T)*prod(S)
    alignment = Base.datatype_alignment(T)
    # Align matrix bases for vectorized accesses to adjacent entries.
    length(S) == 2 && (alignment = max(alignment,16))
    # Raw stack storage avoids the zeroing emitted for MArray allocations. The
    # lifetime starts here and lasts through the inlined caller's use of the buffer.
    pointer_type, lifetime_start = VERSION < v"1.11" ? ("i8*", "llvm.lifetime.start.p0i8") :
        ("ptr", "llvm.lifetime.start.p0")
    ir = """
    declare void @$lifetime_start(i64, $pointer_type)
    define $pointer_type @alloc_scratch() alwaysinline {
        %buffer = alloca i8, i64 $bytes, align $alignment
        call void @$lifetime_start(i64 $bytes, $pointer_type %buffer)
        ret $pointer_type %buffer
    }
    """
    quote
        ptr = Base.llvmcall(($ir,"alloc_scratch"),Core.LLVMPtr{$T,0},Tuple{})
        CUDAThreadLocalArray{$T,$(length(S)),$S}(ptr)
    end
end

@inline function local_zeros(::Type{T},shape::Val) where T
    array = local_scratch(T,shape)
    fill!(array,zero(T))
end

@inline function Base.fill!(array::CUDAThreadLocalArray{T},value) where T
    n = cuda_index_type(length(array))(length(array))
    i = zero(n)
    while i < n
        i += one(i)
        @inbounds array[Int(i)] = value
    end
    array
end

# Internal result of compilation, rather than a user-supplied kernel contract.
struct CUDAAssemblyProgram{K,D}
    kernel::K
    data::D
end

Adapt.@adapt_structure CUDAAssemblyProgram

struct CUDAType{T} end
cuda_capture(::Type{T}) where T = CUDAType{T}()
cuda_capture(value) = value
@inline cuda_input(::CUDAType{T}) where T = T
@inline cuda_input(value) = value

# Initialization placement has already identified invariant data, parameter
# refresh, and private buffers. Host initialization produces explicit numeric
# inputs; the generated device function owns all scratch and parameter refresh.
function GT.ast_assembly_workspaces(::GT.CUDAExecution,captured,parameters,kernel,
        setup,refresh,execution,buffer_names,scopes,guard,guard_node,guard_scope)
    private_names = Set{Symbol}()
    shared_setup,private_setup = Any[],Any[]
    function uses_name(node,names)
        name = GT.ast_leaf_value(node)
        name isa Symbol && return name in names
        GT.ast_head(node) === :quote && return false
        any(child -> uses_name(child,names),GT.ast_children(node))
    end
    for statement in setup
        binding = GT.ast_is_definition(statement) ? GT.ast_leaf_value(GT.ast_lhs(statement)) : nothing
        if binding in buffer_names || uses_name(GT.ast_rhs(statement),private_names)
            push!(private_setup,statement)
            binding isa Symbol && push!(private_names,binding)
        else
            push!(shared_setup,statement)
        end
    end
    # CPU initialization can hoist thousands of individual reference-table
    # entries. Keep their tables in device memory instead of expanding the
    # kernel argument into one scalar/static vector per tabulated entry.
    numeric_getters = (:face_nodes,:node_coordinates,:reference_gradients,
        :reference_tabulator,:domain_faces,:reference_weights)
    function symbols!(names,node)
        name = GT.ast_leaf_value(node)
        if name isa Symbol
            push!(names,name)
        elseif GT.ast_head(node) !== :quote
            foreach(child -> symbols!(names,child),GT.ast_children(node))
        end
        names
    end
    function scratch_allocation(statement)
        children = GT.ast_children(GT.ast_rhs(statement))
        if GT.ast_leaf_value(first(children)) === :alloc_zeros
            # Element matrices are cleared per element and tabulated values are filled before use.
            return local_scratch, children[3], children[4:end]
        end
        local_zeros, children[2], children[3:end]
    end
    is_scratch(statement) = GT.ast_is_definition(statement) &&
        GT.ast_leaf_value(GT.ast_lhs(statement)) in buffer_names
    # Prototype computations determine scratch element types on the CPU. Keep
    # their dependency chain there, even when it reads a reference-table entry.
    host_required = Set{Symbol}()
    for statement in private_setup
        is_scratch(statement) && symbols!(host_required,GT.ast_rhs(statement))
    end
    for statement in refresh
        is_scratch(statement) || continue
        _,_,dims = scratch_allocation(statement)
        foreach(dim -> symbols!(host_required,dim),dims)
    end
    for statement in reverse(shared_setup)
        binding = GT.ast_leaf_value(GT.ast_lhs(statement))
        binding in host_required && symbols!(host_required,GT.ast_rhs(statement))
    end
    numeric_arrays,device_names = Set{Symbol}(),Set{Symbol}()
    host_setup,device_reads = Any[],Any[]
    for statement in shared_setup
        binding = GT.ast_leaf_value(GT.ast_lhs(statement))
        rhs = GT.ast_rhs(statement)
        if GT.ast_is_call(rhs) && GT.ast_leaf_value(first(GT.ast_children(rhs))) in numeric_getters
            push!(numeric_arrays,binding)
        end
        table_read = GT.ast_is_index(rhs) && GT.ast_leaf_value(first(GT.ast_children(rhs))) in numeric_arrays
        if !(binding in host_required) && (table_read || uses_name(rhs,device_names))
            push!(device_reads,statement)
            push!(device_names,binding)
        else
            push!(host_setup,statement)
        end
    end
    shared_setup = host_setup
    # Evaluate scratch dimensions on the host and encode them in the descriptor
    # type. The element type can depend on the types of the runtime parameters.
    host_descriptors = Any[]
    function lower_scratch(statement)
        is_scratch(statement) || return statement
        binding = GT.ast_lhs(statement)
        constructor,element_type,dims = scratch_allocation(statement)
        dim_names = Set{Symbol}()
        foreach(dim -> symbols!(dim_names,dim),dims)
        any(name -> get(scopes,name,0) != 0 || name in private_names,dim_names) &&
            throw(ArgumentError("CUDA assembly requires parameter-independent scratch dimensions"))
        shape = GT.ast_leaf(gensym(:scratch_shape))
        descriptor = GT.ast_call(GT.ast_leaf(cuda_scratch_shape),dims...)
        push!(host_descriptors,GT.ast_definition(shape,descriptor))
        GT.ast_definition(binding,GT.ast_call(GT.ast_leaf(constructor),element_type,shape))
    end
    device_setup = map(lower_scratch,private_setup)
    device_refresh = map(lower_scratch,refresh)

    argument_nodes(args) = GT.ast_leaf_value(args) isa Symbol ? (args,) : GT.ast_children(args)
    allocation,worker = GT.ast_children(first(GT.ast_children(kernel)))
    data,params = GT.ast_leaf(gensym(:data)),GT.ast_leaf(gensym(:parameters))
    parameter_args = first(GT.ast_children(parameters))
    parameter_bindings = [GT.ast_definition(node,GT.ast_index(params,GT.ast_leaf(i)))
        for (i,node) in enumerate(argument_nodes(parameter_args))]
    device_body = GT.ast_block(Any[parameter_bindings...,device_reads...,device_setup...,
        device_refresh...,execution...,GT.ast_leaf(nothing)])
    # Only values actually used by numeric code enter the device argument. Mesh,
    # space, quadrature and other host-only intermediates remain in preparation.
    captured_args = first(GT.ast_children(captured))
    candidates = Any[argument_nodes(captured_args)...]
    append!(candidates,[GT.ast_lhs(s) for s in vcat(shared_setup,host_descriptors) if GT.ast_is_definition(s)])
    used = symbols!(Set{Symbol}(),device_body)
    inputs = unique(GT.ast_leaf_value,filter(node -> GT.ast_leaf_value(node) in used,candidates))
    # Fresh local names are essential: assigning an outer preparation binding
    # inside a nested lambda would capture a Core.Box instead of a device value.
    locals = Dict(name => GT.ast_leaf(gensym(name)) for name in map(GT.ast_leaf_value,inputs))
    function rename(node)
        name = GT.ast_leaf_value(node)
        name isa Symbol && return get(locals,name,node)
        GT.ast_head(node) === :quote && return node
        children = GT.ast_children(node)
        isempty(children) && return node
        GT.ast_replace_children(node,map(rename,children)...)
    end
    input_bindings = map(enumerate(inputs)) do (i,node)
        input = GT.ast_call(GT.ast_leaf(cuda_input),GT.ast_index(data,GT.ast_leaf(i)))
        GT.ast_definition(locals[GT.ast_leaf_value(node)],input)
    end
    device_args = GT.ast_expr(:tuple,allocation,data,worker,params)
    device_kernel = GT.ast_replace_children(kernel,device_args,
        GT.ast_block(Any[input_bindings...,rename(device_body)]))
    captured_inputs = map(node -> GT.ast_call(GT.ast_leaf(cuda_capture),node),inputs)
    captured_data = GT.ast_expr(:tuple,captured_inputs...)
    prepared = GT.ast_expr(:tuple,device_kernel,captured_data)
    host_body = GT.ast_block(Any[shared_setup...,host_descriptors...,prepared])
    GT.ast_replace_children(captured,captured_args,host_body)
end

function GT.generate_assemble_matrix(execution::GT.CUDAExecution,contribution,trial_space,test_space,
        prototype;parameters,optimize_options)
    prepared = generate_cuda_matrix(execution,contribution,trial_space,test_space,prototype;
        parameters,optimize_options)
    device = Adapt.adapt(CuArray,prepared)
    isbitstype(typeof(CUDA.cudaconvert(device))) || throw(ArgumentError(
        "Generated CUDA kernel inputs must be device-compatible numeric data"))
    device
end

# Bake the prepared lookup width and packed dimensions into the device call.
function ast_cuda_contributions(node,index_type)
    GT.ast_is_leaf(node) && return node
    GT.ast_head(node) === :quote && return node
    children = GT.ast_children(node)
    if GT.ast_is_call(node) && length(children) == 7 &&
            GT.ast_leaf_value(first(children)) == :(GT.contribute_element_block_atomic!)
        _,alloc,values,face,base,rows,cols = children
        dims = (GT.ast_leaf_value(rows),GT.ast_leaf_value(cols))
        if all(d -> d isa Integer,dims)
            shape = Val(dims)
            return GT.ast_call(GT.ast_leaf(cuda_contribute_element_block_atomic!),alloc,values,face,base,
                GT.ast_leaf(shape),GT.ast_leaf(index_type))
        end
    end
    GT.ast_replace_children(node,map(child->ast_cuda_contributions(child,index_type),children)...)
end

function generate_cuda_matrix(execution,contribution,trial_space,test_space,prototype;
        parameters,optimize_options)
    GT.num_fields(trial_space) == GT.num_fields(test_space) == 1 ||
        throw(ArgumentError("CUDA matrix generation currently supports single-field CG spaces"))
    domain = GT.domain(contribution)
    GT.num_dims(domain) == GT.num_dims(GT.domain(trial_space)) == GT.num_dims(GT.domain(test_space)) ||
        throw(ArgumentError("CUDA matrix generation currently supports volume assembly"))
    GT.num_faces(domain) == 0 && return CUDAAssemblyProgram(empty_cuda_assembly!,())
    term = GT.optimize(GT.write_assemble_matrix(contribution,trial_space,test_space;
        optimize_options,allocation=prototype,workspace_backend=execution))
    captured,values = GT.capture(GT.parametrize(term,parameters...))
    domain_face_symbol = GT.bindings(term)[2]
    expression = GT.ast_element_range_assembly(GT.expression(captured),domain_face_symbol)
    lowered = GT.ast_initialization_code_motion(expression;workspaces=true,workspace_backend=execution)
    storage = GT.element_csc_matrix_storage(prototype.allocation)
    positions = GT.element_csc_single_domain_cache(storage).positions
    lowered = ast_cuda_contributions(lowered,cuda_index_type(length(positions)))
    kernel,data = GT.evaluate(lowered,values,execution)
    CUDAAssemblyProgram(kernel,data)
end

@inline empty_cuda_assembly!(allocation,data,worker,parameters) = nothing

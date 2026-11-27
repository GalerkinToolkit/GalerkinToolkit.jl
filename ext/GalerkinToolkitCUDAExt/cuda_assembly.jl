# The prepared position table is read-only during device assembly.
@inline function GT.element_csc_position(positions::CUDA.CuDeviceArray{Ti,2},entry,face) where Ti<:Integer
    @inbounds CUDA.Const(positions)[entry,face]
end

const CUDAElementCSCPositions = Union{
    CUDA.CuDeviceArray{<:Integer,2},
    Transpose{<:Integer,<:CUDA.CuDeviceArray{<:Integer,2}},
}
const CUDA_SCATTER_BATCH_WIDTH = 32

@inline cuda_index_type(n) = n <= typemax(Int32) ? Int32 : Int

# Device atomics act on real numbers only, so a complex entry is added as two
# real atomic adds on its (real, imaginary) storage.
@inline function GT.atomic_add!(nzval::CUDA.CuDeviceArray{Complex{T},1,A},slot,v) where {T,A}
    w = convert(Complex{T},v)
    re = reinterpret(Core.LLVMPtr{T,A},pointer(nzval,slot))
    CUDA.atomic_add!(re,real(w))
    CUDA.atomic_add!(re+sizeof(T),imag(w))
    nothing
end

@inline function cuda_scatter_lookup(positions::CUDA.CuDeviceArray,face,base,::Type{I}) where I
    offset = (I(face)-one(I))*I(size(positions,1)) + I(base)
    CUDA.Const(positions),offset,one(I)
end

@inline function cuda_scatter_lookup(positions::Transpose,face,base,::Type{I}) where I
    data = parent(positions)
    stride = I(size(data,1))
    offset = I(face)-one(I) + I(base)*stride
    CUDA.Const(data),offset,stride
end

@inline function GT.contribute_element_csc_block_atomic!(
        nzval,positions::CUDAElementCSCPositions,vals,face,base,shape)
    # A packed table can exceed Int32 even when all its stored slots fit.
    if length(positions) <= typemax(Int32)
        cuda_scatter_packed!(nzval,positions,vals,face,base,shape[1],Int32)
    else
        cuda_scatter_packed!(nzval,positions,vals,face,base,shape[1],Int)
    end
end

# The generated caller supplies the packed shape and lookup width from the host.
@inline function cuda_contribute_element_block_atomic!(malloc,vals,face,base,
        shape::Val{S},::Type{I}) where {S,I}
    base == -1 && return nothing
    alloc = GT.element_csc_matrix_storage(malloc.allocation)
    cache = GT.element_csc_single_domain_cache(alloc)
    if min(S...) <= GT.MAX_UNROLL_SIZE && size(vals) == S
        lookup,offset,stride = cuda_scatter_lookup(cache.positions,face,base,I)
        cuda_scatter_dofs!(alloc.nzval,lookup,vals,offset,stride,shape)
    else
        cuda_scatter_packed!(alloc.nzval,cache.positions,vals,face,base,S[1],I)
    end
    nothing
end

@inline function cuda_scatter_dof!(nzval,lookup,vals,offset::I,stride,entry) where I
    @inbounds begin
        slot = lookup[offset + I(entry)*stride + one(I)]
        if slot != 0
            GT.atomic_add!(nzval,slot,vals[Int(entry+one(entry))])
        end
    end
    nothing
end

@inline @generated function cuda_scatter_dofs!(nzval,lookup,vals,offset,stride,::Val{S}) where S
    rows,cols = S
    index_type = cuda_index_type(prod(S))
    function dof_loop(body,name,n)
        if n <= GT.MAX_UNROLL_SIZE
            return GT.ast_block([body(GT.ast_leaf(index_type(i))) for i in 0:n-1])
        end
        index = GT.ast_leaf(gensym(name))
        quote
            $index = $(zero(index_type))
            while $index < $(index_type(n))
                $(body(index))
                $index += $(one(index_type))
            end
        end
    end
    loops = dof_loop(:col,cols) do col
        dof_loop(:row,rows) do row
            entry = :($col * $(index_type(rows)) + $row)
            GT.ast_call(GT.ast_leaf(:cuda_scatter_dof!),:nzval,:lookup,:vals,:offset,:stride,entry)
        end
    end
    GT.ast_block([loops,GT.ast_leaf(nothing)])
end

@inline function cuda_scatter_packed!(nzval,positions,vals,face,base,nrows,::Type{I}) where I
    lookup,offset,stride = cuda_scatter_lookup(positions,face,base,I)
    n = length(vals)
    if I === Int32 && CUDA_SCATTER_BATCH_WIDTH <= n <= typemax(Int32) && size(vals,1) == nrows
        cuda_scatter_batches!(nzval,lookup,vals,offset,stride,Int32(n))
    elseif n <= typemax(Int32)
        cuda_scatter_atomic!(nzval,lookup,vals,offset,stride,nrows,Int32)
    else
        cuda_scatter_atomic!(nzval,lookup,vals,offset,stride,nrows,Int)
    end
end

@inline function cuda_scatter_batches!(nzval,lookup,vals,offset,stride,n)
    width = Int32(CUDA_SCATTER_BATCH_WIDTH)
    start = Int32(0)
    while start <= n-width
        cuda_scatter_batch!(nzval,lookup,vals,offset+start*stride,start,stride,Val(CUDA_SCATTER_BATCH_WIDTH))
        start += width
    end
    cuda_scatter_tail!(nzval,lookup,vals,offset,stride,start,n)
end

@inline function cuda_scatter_tail!(nzval,lookup,vals,offset,stride,start,n)
    @inbounds while start < n
        slot = lookup[offset + start*stride + Int32(1)]
        start += Int32(1)
        if slot != 0
            GT.atomic_add!(nzval,slot,vals[Int(start)])
        end
    end
    nothing
end

@inline function cuda_scatter_batch!(nzval,lookup,vals,offset,start,stride,::Val{B}) where B
    # Materialize this batch before the slot checks so the compiler can group its loads.
    entries = GT.SVector{B}(ntuple(j -> (@inbounds vals[Int(start+Int32(j))]),Val(B)))
    ntuple(Val(B)) do j
        @inbounds begin
            slot = lookup[offset + (Int32(j)-Int32(1))*stride + Int32(1)]
            if slot != 0
                GT.atomic_add!(nzval,slot,entries[j])
            end
        end
        nothing
    end
    nothing
end

@inline function cuda_scatter_atomic!(nzval,lookup,vals,offset::I,stride::I,nrows,::Type{J}) where {I,J}
    rows,cols = J(size(vals,1)),J(size(vals,2))
    packed_rows = I(nrows)
    col = zero(J)
    @inbounds while col < cols
        lookup_col = offset + I(col)*packed_rows*stride
        tile_col = col*rows
        row = zero(J)
        while row < rows
            slot = lookup[lookup_col + I(row)*stride + one(I)]
            if slot != 0
                GT.atomic_add!(nzval,slot,vals[Int(tile_col + row + one(J))])
            end
            row += one(J)
        end
        col += one(J)
    end
    nothing
end

# The output shares the numeric buffer used by all selected domain arguments.
# Keep sparse construction on the host; only array storage lives on the device.
@static if isdefined(CUDA,:cuSPARSE)
    const CuSparseMatrixCSC = CUDA.cuSPARSE.CuSparseMatrixCSC
else
    const CuSparseMatrixCSC = CUDA.CUSPARSE.CuSparseMatrixCSC
end

# Repeat the index/value bounds so this specializes the CPU fallback's signature.
function GT.compress(alloc::GT.ElementCSCMatrixAllocation{Tv,Ti,C,P,R,V};
        reuse=Val(false)) where {Tv,Ti,C,P<:CuVector{Ti},R<:CuVector{Ti},V<:CuVector{Tv}}
    A = CuSparseMatrixCSC{Tv,Ti}(alloc.colptr,alloc.rowval,alloc.nzval,
        (alloc.nrows,alloc.ncols))
    GT.val_parameter(reuse) ? (A,nothing) : A
end

function GT.compress!(alloc::GT.ElementCSCMatrixAllocation,
        A::CuSparseMatrixCSC,cache)
    values = nonzeros(A)
    values === alloc.nzval || copyto!(values,alloc.nzval)
    A
end

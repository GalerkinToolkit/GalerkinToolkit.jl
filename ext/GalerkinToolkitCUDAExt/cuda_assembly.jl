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

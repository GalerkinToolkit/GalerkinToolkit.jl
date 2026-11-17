

function default_vector_assembly_method()
    monolithic_assembly(coo_assembly())
end

function default_matrix_assembly_method()
    monolithic_assembly(coo_assembly())
end

function vector_assembly_method(;vector=default_vector_assembly_method(),matrix=nothing)
    vector
end

function matrix_assembly_method(;matrix=default_matrix_assembly_method(),vector=nothing)
    matrix
end

function vector_assembly_options(;vector=(;),matrix=nothing)
    vector
end

function matrix_assembly_options(;matrix=(;),vector=nothing)
    matrix
end

function allocate_vector(::Type{T},space_i::AbstractSpace,domains::AbstractDomain...;
        free_or_dirichlet = FREE,
        assembly_method = (;),
        assembly_options = (;),
        kwargs... # hide other unused args
    ) where T
    free_or_diri_i = free_or_dirichlet
    dofs_i = GT.dofs(space_i,free_or_diri_i)
    method = vector_assembly_method(;assembly_method...)
    options = vector_assembly_options(;assembly_options...)
    counter = GT.counter(method,T,dofs_i;options...)
    dof_map_i = free_or_diri_i == FREE ? identity : -
    fields_i = ntuple(identity,num_fields(space_i))
    for domain in domains
        nfaces = num_faces(domain)
        field_face_dofs_i = map(field_i->dofs_accessor(GT.field(space_i,field_i),domain),fields_i)
        field_face_nfaces_around_i = map(field_i->num_faces_around_accesor(GT.domain(GT.field(space_i,field_i)),domain),fields_i)
        allocate_vector_barrier!(counter,nfaces,fields_i,field_face_dofs_i,free_or_diri_i,dof_map_i,field_face_nfaces_around_i)
    end
    alloc = allocate(counter)
    VectorAllocation(alloc,free_or_dirichlet,dof_map_i)
end

function allocate_vector_barrier!(counter,nfaces,fields_i,field_face_dofs_i,free_or_diri_i,map_i,field_face_nfaces_around_i)
    if ! do_loop(counter)
        return counter
    end
    for face in 1:nfaces
        for field_i in fields_i
            nfaces_around_i = field_face_nfaces_around_i[field_i](face)
            for face_around_i in 1:nfaces_around_i
                dofs_i = field_face_dofs_i[field_i](face,face_around_i)
                for dof_i in dofs_i
                    if skip_dof(dof_i,free_or_diri_i)
                        continue
                    end
                    contribute!(counter,nothing,map_i(dof_i),field_i)
                end
            end
        end
    end
    counter
end

# TODO: allocate the matrix on GPUs directly
# TODO: generate code for allocation
function allocate_matrix(::Type{T},space_i::AbstractSpace,space_j::AbstractSpace,domains::AbstractDomain...;
        free_or_dirichlet = (FREE,FREE),
        assembly_method = (;),
        assembly_options = (;),
        kwargs... # hide other unused args
        ) where T
    free_or_diri_i, free_or_diri_j = free_or_dirichlet
    method = matrix_assembly_method(;assembly_method...)
    options = matrix_assembly_options(;assembly_options...)
    dofs_i = GT.dofs(space_i,free_or_diri_i)
    dofs_j = GT.dofs(space_j,free_or_diri_j)
    counter = GT.counter(method,T,dofs_i,dofs_j;options...)
    dof_map_i = free_or_diri_i == FREE ? identity : -
    dof_map_j = free_or_diri_j == FREE ? identity : -
    dof_maps = (dof_map_i,dof_map_j)
    fields_i = ntuple(identity,num_fields(space_i))
    fields_j = ntuple(identity,num_fields(space_j))
    args = (domains, space_i, space_j, fields_i, fields_j, free_or_diri_i, free_or_diri_j, dof_map_i, dof_map_j)
    allocate_matrix_sparsity!(method,counter,args...)
    alloc_0 = allocate(counter)
    alloc = if (alloc_0 isa CachedCSCCounterStep1) || (alloc_0 isa MonolithicAssemblyCounter) # TODO: make it more general
        allocate_matrix_sparsity!(method,alloc_0,args...)
        allocate_matrix_from_sparsity(method,alloc_0,args)
    else 
        alloc_0
    end
    malloc = MatrixAllocation(alloc,free_or_dirichlet,dof_maps)
    prepare_matrix_allocation(method,malloc,space_i,space_j,domains)
end

function allocate_matrix_barrier!(counter, domains, space_i, space_j, fields_i, fields_j, free_or_diri_i, free_or_diri_j, dof_map_i, dof_map_j)
    for domain in domains
        nfaces = num_faces(domain)
        field_face_dofs_i = map(field_i->dofs_accessor(GT.field(space_i,field_i),domain),fields_i)
        field_face_dofs_j = map(field_j->dofs_accessor(GT.field(space_j,field_j),domain),fields_j)
        field_face_nfaces_around_i = map(field_i->num_faces_around_accesor(GT.domain(GT.field(space_i,field_i)),domain),fields_i)
        field_face_nfaces_around_j = map(field_j->num_faces_around_accesor(GT.domain(GT.field(space_j,field_j)),domain),fields_j)
        allocate_matrix_barrier_impl!(counter,nfaces,fields_i,fields_j,field_face_dofs_i,field_face_dofs_j,free_or_diri_i,free_or_diri_j,dof_map_i,dof_map_j,field_face_nfaces_around_i,field_face_nfaces_around_j)
    end
end

function allocate_matrix_barrier_impl!(counter,nfaces,fields_i,fields_j,field_face_dofs_i,field_face_dofs_j,free_or_diri_i,free_or_diri_j,map_i,map_j,field_face_nfaces_around_i,field_face_nfaces_around_j)
    if ! do_loop(counter)
        return counter
    end
    for face in 1:nfaces 
        for field_j in fields_j 
            nfaces_around_j = field_face_nfaces_around_j[field_j](face)
            for field_i in fields_i 
                nfaces_around_i = field_face_nfaces_around_i[field_i](face)
                for face_around_j in 1:nfaces_around_j 
                    dofs_j = field_face_dofs_j[field_j](face,face_around_j)
                    for face_around_i in 1:nfaces_around_i 
                        dofs_i = field_face_dofs_i[field_i](face,face_around_i)
                        for dof_j in dofs_j
                            if skip_dof(dof_j,free_or_diri_j)
                                continue
                            end
                            for dof_i in dofs_i
                                if skip_dof(dof_i,free_or_diri_i)
                                    continue
                                end
                                contribute!(counter,nothing,map_i(dof_i),map_j(dof_j),field_i,field_j)
                            end
                        end
                    end
                end
            end
        end
    end
    counter
end

function skip_dof(i,free_or_dirichlet)
    (i>0 && free_or_dirichlet == DIRICHLET) || (i<0 && free_or_dirichlet == FREE)
end

struct MatrixAllocation{A,B,C} <: AbstractType
    allocation::A
    free_or_dirichlet::B
    dof_maps::C
end

struct VectorAllocation{A,B,C} <: AbstractType
    allocation::A
    free_or_dirichlet::B
    dof_map::C
end

Adapt.@adapt_structure MatrixAllocation

const AssemblyAllocation = Union{VectorAllocation,MatrixAllocation}

Base.eltype(alloc::AssemblyAllocation) = eltype(alloc.allocation)

"""
    contribute!(allocation, values, dofs...)

Accumulate contributions in place. All contribution methods, including sparsity
counters and atomic/element-indexed variants, return `nothing`.
"""
function contribute!(malloc::VectorAllocation,vals,dofs_i,field_i=1)
    free_or_diri_i  = malloc.free_or_dirichlet
    map_i = malloc.dof_map
    alloc = malloc.allocation
    ndofs_i = length(dofs_i)
    for i in 1:ndofs_i
        dof_i = dofs_i[i]
        if skip_dof(dof_i,free_or_diri_i)
            continue
        end
        contribute!(alloc,vals[i],map_i(dof_i),field_i)
    end
    nothing
end

function contribute!(malloc::MatrixAllocation,vals,dofs_i,dofs_j,field_i=1,field_j=1)
    free_or_diri_i, free_or_diri_j = malloc.free_or_dirichlet
    map_i,map_j = malloc.dof_maps
    alloc = malloc.allocation
    ndofs_i = length(dofs_i)
    ndofs_j = length(dofs_j)
    for j in 1:ndofs_j
        dof_j = dofs_j[j]
        if skip_dof(dof_j,free_or_diri_j)
            continue
        end
        for i in 1:ndofs_i
            dof_i = dofs_i[i]
            if skip_dof(dof_i,free_or_diri_i)
                continue
            end
            contribute!(alloc,vals[i,j],map_i(dof_i),map_j(dof_j),field_i,field_j)
        end
    end
    nothing
end

function reset!(malloc::AssemblyAllocation)
    reset!(malloc.allocation)
    malloc
end

function compress(malloc::AssemblyAllocation;reuse=Val(false))
    compress(malloc.allocation;reuse)
end

function compress!(malloc::AssemblyAllocation,A,A_cache)
    compress!(malloc.allocation,A,A_cache)
    A
end

# Monolithic assembly related

function monolithic_assembly(method)
    MonolithicAssembly(method)
end

struct MonolithicAssembly{A} <: AbstractType
    method::A
end

function values_from_vector!(s::MonolithicAssembly,coeffs::BVector,x)
    dofs = axes(coeffs,1)
    nfields = blocklength(dofs)
    for field in 1:nfields
        pend = blocklasts(dofs)[field]
        pini = 1 + pend - length(blocks(dofs)[field])
        copyto!(coeffs[Block(field)],view(x,pini:pend))
    end
    coeffs
end

function values_from_vector(s::MonolithicAssembly,dofs::BRange,x)
    nfields = blocklength(dofs)
    map(1:nfields) do field
        pend = blocklasts(dofs)[field]
        pini = 1 + pend - length(blocks(dofs)[field])
        x[pini:pend]
    end |> BVector
end

function vector_from_values!(s::MonolithicAssembly,x,coeffs::BVector)
    dofs = axes(coeffs,1)
    nfields = blocklength(dofs)
    for field in 1:nfields
        pend = blocklasts(dofs)[field]
        pini = 1 + pend - length(blocks(dofs)[field])
        copyto!(view(x,pini:pend),coeffs[Block(field)])
    end
    x
end

function vector_from_values(s::MonolithicAssembly,coeffs::BVector)
    collect(coeffs)
end

function vector_from_values!(s::MonolithicAssembly,x,coeffs)
    vector_from_values!(s.method,x,coeffs)
end

function vector_from_values(s::MonolithicAssembly,coeffs)
    vector_from_values(s.method,coeffs)
end

function values_from_vector!(s::MonolithicAssembly,coeffs,x)
    values_from_vector!(s.method,coeffs,x)
end

function values_from_vector(s::MonolithicAssembly,dofs,x)
    values_from_vector(s.method,dofs,x)
end

function vector_from_values!(s,x,coeffs)
    if coeffs !== x
        copyto!(x,coeffs)
    end
    x
end

function vector_from_values(s,coeffs)
    coeffs
end

function values_from_vector!(s,coeffs,x)
    if coeffs !== x
        copyto!(coeffs,x)
    end
    coeffs
end

function values_from_vector(s,dofs,x)
    x
end

function counter(s::MonolithicAssembly,::Type{T},dofs_i;block_mask=nothing,kwargs...) where T
    GT.counter(s.method,T,dofs_i;kwargs...)
end

function counter(s::MonolithicAssembly,::Type{T},dofs_i::BRange;
    block_mask=fill(true,blocklength(dofs_i)),
    kwargs...) where T
    counter = GT.counter(s.method,T,dofs_i;kwargs...)
    offsets_i = blocklasts(dofs_i) .- map(length,blocks(dofs_i))
    offsets = (offsets_i,)
    MonolithicAssemblyCounter(counter,offsets,block_mask)
end

function counter(s::MonolithicAssembly,::Type{T},dofs_i,dofs_j;
        block_mask=nothing,kwargs...) where T
    GT.counter(s.method,T,dofs_i,dofs_j;kwargs...)
end

function counter(s::MonolithicAssembly,::Type{T},dofs_i::BRange,dofs_j::BRange;
    block_mask=fill(true,blocklength(dofs_i),blocklength(dofs_j)),
    kwargs...) where T
    counter = GT.counter(s.method,T,dofs_i,dofs_j;kwargs...)
    offsets_i = blocklasts(dofs_i) .- map(length,blocks(dofs_i))
    offsets_j = blocklasts(dofs_j) .- map(length,blocks(dofs_j))
    offsets = (offsets_i,offsets_j)
    MonolithicAssemblyCounter(counter,offsets,block_mask)
end

struct MonolithicAssemblyCounter{A,B,C} <: AbstractType
    counter::A
    offsets::B
    block_mask::C
end

function do_loop(mcounter::MonolithicAssemblyCounter)
    do_loop(mcounter.counter)
end

function reset!(counter::MonolithicAssemblyCounter)
    reset!(counter.counter)
    counter
end

function contribute!(counter::MonolithicAssemblyCounter,v,i,field_i)
    if counter.block_mask[field_i]
        offsets_i, = counter.offsets
        i2 = offsets_i[field_i]+i
        contribute!(counter.counter,v,i2,1)
    end
    nothing
end

function contribute!(counter::MonolithicAssemblyCounter,v,i,j,field_i,field_j)
    if counter.block_mask[field_i,field_j]
        offsets_i,offsets_j = counter.offsets
        i2 = offsets_i[field_i]+i
        j2 = offsets_j[field_j]+j
        contribute!(counter.counter,v,i2,j2,1,1)
    end
    nothing
end

function allocate(mcounter::MonolithicAssemblyCounter)
    (;counter,offsets,block_mask) = mcounter
    alloc = allocate(counter)
    if ((counter isa CSCCounter) && val_parameter(counter.use_cache)) || (counter isa CachedCSCCounterStep1)
        MonolithicAssemblyCounter(alloc, offsets, block_mask)
    else
        reset!(counter)
        MonolithicAssemblyAllocation(alloc,offsets,block_mask)
    end
end

function finalize(mcounter::MonolithicAssemblyCounter)
    (;counter,offsets,block_mask) = mcounter
    alloc = finalize(counter)
    MonolithicAssemblyAllocation(alloc,offsets,block_mask)
end

struct MonolithicAssemblyAllocation{A,B,C} <: AbstractType
    allocation::A
    offsets::B
    block_mask::C
end

Adapt.@adapt_structure MonolithicAssemblyAllocation

Base.eltype(alloc::MonolithicAssemblyAllocation) = eltype(alloc.allocation)

function contribute!(alloc::MonolithicAssemblyAllocation,v,i,field_i)
    if alloc.block_mask[field_i]
        offsets_i, = alloc.offsets
        i2 = offsets_i[field_i]+i
        contribute!(alloc.allocation,v,i2,1)
    end
    nothing
end

function reset!(alloc::MonolithicAssemblyAllocation)
    reset!(alloc.allocation)
    alloc
end

function contribute!(alloc::MonolithicAssemblyAllocation,v,i,j,field_i,field_j)
    if alloc.block_mask[field_i,field_j]
        offsets_i,offsets_j = alloc.offsets
        i2 = offsets_i[field_i]+i
        j2 = offsets_j[field_j]+j
        contribute!(alloc.allocation,v,i2,j2,1,1)
    end
    nothing
end

function compress(malloc::MonolithicAssemblyAllocation;reuse=Val(true))
    compress(malloc.allocation;reuse)
end

function compress!(malloc::MonolithicAssemblyAllocation,A,cache)
    compress!(malloc.allocation,A,cache)
end

# COO related

function coo_assembly()
    COOAssembly()
end

struct COOAssembly <: AbstractType end

function counter(s::COOAssembly,::Type{T},dofs_i;index_type=Int32,eltype=T,vector_type=Vector{eltype}) where T
    @assert vector_type <: Vector "Case not implemented"
    ni = length(dofs_i)
    COOVectorCounter(vector_type,index_type,ni,Ref(0))
end

function counter(s::COOAssembly,::Type{T},dofs_i,dofs_j;
    eltype=T,index_type=Int32,matrix_type=SparseMatrixCSC{eltype,index_type}) where T
    ni = length(dofs_i)
    nj = length(dofs_j)
    COOMatrixCounter(matrix_type,index_type,ni,nj,Ref(0))
end

struct COOVectorCounter{A,B} <: AbstractType
    vector_type::Type{A}
    index_type::Type{B}
    nrows::Int
    nnz::Base.RefValue{Int}
end

struct COOMatrixCounter{A,B} <: AbstractType
    matrix_type::Type{A}
    index_type::Type{B}
    nrows::Int
    ncols::Int
    nnz::Base.RefValue{Int}
end

const COOCounter = Union{COOVectorCounter,COOMatrixCounter}

function do_loop(counter::COOCounter)
    true
end

function reset!(counter::COOCounter)
    counter.nnz[] = 0
    counter
end

function contribute!(counter::COOVectorCounter,v,i,field_i)
    @boundscheck @assert v === nothing
    @boundscheck @assert length(i) == 1
    @boundscheck @assert field_i == 1
    counter.nnz[] += 1
    nothing
end

function contribute!(counter::COOMatrixCounter,v,i,j,field_i,field_j)
    @boundscheck @assert v === nothing
    @boundscheck @assert length(i) == 1
    @boundscheck @assert length(j) == 1
    @boundscheck @assert field_i == 1
    @boundscheck @assert field_j == 1
    counter.nnz[] += 1
    nothing
end

function allocate(counter::COOVectorCounter)
    Ti = counter.index_type
    T = eltype(counter.vector_type)
    n = counter.nnz[]
    I = zeros(Ti,n)
    V = zeros(T,n)
    reset!(counter)
    COOVectorAllocation(I,V,counter)
end

function allocate(counter::COOMatrixCounter)
    Ti = counter.index_type
    T = eltype(counter.matrix_type)
    n = counter.nnz[]
    I = zeros(Ti,n)
    J = zeros(Ti,n)
    V = zeros(T,n)
    reset!(counter)
    COOMatrixAllocation(I,J,V,counter)
end

struct COOVectorAllocation{Tv,Ti,C} <: AbstractType
    I::Vector{Ti}
    V::Vector{Tv}
    counter::C
end


struct COOMatrixAllocation{Tv,Ti,C} <: AbstractType
    I::Vector{Ti}
    J::Vector{Ti}
    V::Vector{Tv}
    counter::C
end

const COOAllocation = Union{COOVectorAllocation,COOMatrixAllocation}

Base.eltype(alloc::COOAllocation) = eltype(alloc.V)

function reset!(alloc::COOAllocation)
    reset!(alloc.counter)
    alloc
end

function contribute!(alloc::COOVectorAllocation,v,i,field_i)
    @boundscheck @assert length(i) == 1
    @boundscheck @assert field_i == 1
    alloc.counter.nnz[] += 1
    n = alloc.counter.nnz[]
    alloc.I[n] = i
    alloc.V[n] = v
    nothing
end

function contribute!(alloc::COOMatrixAllocation,v,i,j,field_i,field_j)
    @boundscheck @assert length(i) == 1
    @boundscheck @assert length(j) == 1
    @boundscheck @assert field_i == 1
    @boundscheck @assert field_j == 1
    alloc.counter.nnz[] += 1
    n = alloc.counter.nnz[]
    alloc.I[n] = i
    alloc.J[n] = j
    alloc.V[n] = v
    nothing
end

function compress(alloc::COOVectorAllocation;reuse=Val(false))
    (;I,V,counter) = alloc
    nrows = counter.nrows
    @assert alloc.counter.vector_type <: Vector
    vec = PartitionedArrays.dense_vector(I,V,nrows)
    cache = nothing
    if val_parameter(reuse)
        (vec, cache)
    else
        vec
    end
end

function compress(alloc::COOMatrixAllocation;reuse=Val(false))
    (;I,J,V,counter) = alloc
    (;nrows,ncols,matrix_type) = counter
    sparse_matrix(matrix_type,I,J,V,nrows,ncols;reuse)
end

function compress!(alloc::COOVectorAllocation,A,cache)
    I = alloc.I
    V = alloc.V
    PartitionedArrays.dense_vector!(A,I,V)
    A
end

function compress!(alloc::COOMatrixAllocation,A,cache)
    V = alloc.V
    sparse_matrix!(A,V,cache)
    A
end

## CSC assembly

function csc_assembly(use_cache=Val(true))
    CSCAssembly(use_cache)
end

struct CSCAssembly{A} <: AbstractType
    use_cache::Val{A}
end

function counter(s::CSCAssembly,::Type{T},dofs_i,dofs_j;
    eltype=T,index_type=Int32,matrix_type=SparseMatrixCSC{eltype,index_type}) where T
    @assert matrix_type <: SparseMatrixCSC
    ni = length(dofs_i)
    nj = length(dofs_j)
    colnnz = zeros(index_type,nj)
    CSCCounter(matrix_type,index_type,ni,nj,colnnz,s.use_cache)
end

struct CSCCounter{A,B,C} <: AbstractType
    matrix_type::Type{A}
    index_type::Type{B}
    nrows::Int
    ncols::Int
    colnnz::Vector{B}
    use_cache::Val{C}
end

# in step 1, we allocate the contribution cache and precompute the colptr & rowval
# in step 2, we precompute the positions of all contributions
struct CachedCSCCounterStep1{A,B} <: AbstractType
    matrix_type::Type{A}
    index_type::Type{B}
    nrows::Int
    ncols::Int
    colptr::Vector{B}
    colnnz::Vector{B}
    rowval::Vector{B}
    cache_count::Base.RefValue{Int}
end

struct CachedCSCCounterStep2{Tv,Ti} <: AbstractType
  nrows::Int
  ncols::Int
  colptr::Vector{Ti}
  rowval::Vector{Ti}
  nzval::Vector{Tv}
  cache_p::Vector{Ti}
  cache_id::Base.RefValue{Int}
end


function do_loop(counter::Union{CSCCounter, CachedCSCCounterStep1, CachedCSCCounterStep2})
    true
end

function reset!(counter::CSCCounter)
    fill!(counter.colnnz,zero(eltype(counter.colnnz)))
    counter
end

function reset!(counter::CachedCSCCounterStep1)
    fill!(counter.colnnz,zero(eltype(counter.colnnz)))
    counter
end

function finalize(counter::CachedCSCCounterStep2)
    CachedCSCAllocation(counter.nrows, counter.ncols, counter.colptr, counter.rowval, counter.nzval, counter.cache_p, Ref(0))
end

function contribute!(counter::CSCCounter,v,i,j,field_i,field_j)
    @boundscheck @assert v === nothing
    @boundscheck @assert length(i) == 1
    @boundscheck @assert length(j) == 1
    @boundscheck @assert field_i == 1
    @boundscheck @assert field_j == 1
    counter.colnnz[j] += 1
    nothing
end

function allocate(counter::CSCCounter)
    Ti = counter.index_type
    T = eltype(counter.matrix_type)
    nrows = counter.nrows
    ncols = counter.ncols
    colnnz = counter.colnnz
    colptr = Vector{Ti}(undef,ncols+1)
    @inbounds for i in 1:ncols
        colptr[i+1] = colnnz[i]
    end
    length_to_ptrs!(colptr)
    ndata = colptr[end] - one(Ti)
    rowval = Vector{Ti}(undef,ndata)
    reset!(counter)

    if val_parameter(counter.use_cache)
        CachedCSCCounterStep1(counter.matrix_type,counter.index_type,nrows,ncols,colptr,colnnz,rowval,Ref(0))
    else
        nzval = zeros(T,ndata)
        CSCAllocation(nrows,ncols,colptr,colnnz,rowval,nzval) # cache_i, cache_j, Ref(0), Ref(0))
    end 
end


function contribute!(counter::CachedCSCCounterStep1,v,i,j,field_i,field_j)
    @boundscheck @assert length(i) == 1
    @boundscheck @assert length(j) == 1
    @boundscheck @assert field_i == 1
    @boundscheck @assert field_j == 1
    @boundscheck @assert v === nothing
    pini = Int(counter.colptr[j])
    pend = pini + Int(counter.colnnz[j]) - 1
    p = searchsortedfirst(counter.rowval,i,pini,pend,Base.Order.Forward)
    if (p>pend)
        # add new entry
        counter.colnnz[j] += 1
        counter.rowval[p] = i
    elseif counter.rowval[p] != i
        # shift one forward from p to pend
        @boundscheck @assert pend+1 < Int(counter.colptr[j+1])
        for k in pend:-1:p
            o = k + 1
            counter.rowval[o] = counter.rowval[k]
        end
        # add new entry
        counter.colnnz[j] += 1
        counter.rowval[p] = i
    else 
        # skip
    end
    counter.cache_count[] += 1
    nothing
end


function contribute!(counter::CachedCSCCounterStep2,v,i,j,field_i,field_j)
    @boundscheck @assert length(i) == 1
    @boundscheck @assert length(j) == 1
    @boundscheck @assert field_i == 1
    @boundscheck @assert field_j == 1
    @boundscheck @assert v === nothing
    pini = Int(counter.colptr[j])
    pend = Int(counter.colptr[j+1]) - 1
    p = searchsortedfirst(counter.rowval,i,pini,pend,Base.Order.Forward)
    cache_id = counter.cache_id[]+1
    counter.cache_p[cache_id] = p
    counter.cache_id[] = cache_id
    nothing
end


function compact_csc_structure!(counter::CachedCSCCounterStep1)
    ncols = counter.ncols
    colnnz = counter.colnnz
    colptr = counter.colptr
    rowval = counter.rowval


    k = 1
    for j in 1:ncols
        pini = Int(colptr[j])
        pend = pini + Int(colnnz[j]) - 1
        for p in pini:pend
            rowval[k] = rowval[p]
            k += 1
        end
    end
    @inbounds for j in 1:ncols
        colptr[j+1] = colnnz[j]
    end
    length_to_ptrs!(colptr)
    nnz = colptr[end]-1
    resize!(rowval,nnz)
    counter
end

function allocate(counter::CachedCSCCounterStep1)
    compact_csc_structure!(counter)
    Ti = counter.index_type
    T = eltype(counter.matrix_type)
    nrows,ncols = counter.nrows,counter.ncols
    colptr,rowval = counter.colptr,counter.rowval
    nnz = length(rowval)
    nzval = zeros(T,nnz)
    cache_p = Vector{Ti}(undef, counter.cache_count[])
    reset!(counter)
    CachedCSCCounterStep2(nrows,ncols,colptr,rowval,nzval,cache_p,Ref(0))
end

struct CSCAllocation{Tv,Ti} <: AbstractType
  nrows::Int
  ncols::Int
  colptr::Vector{Ti}
  colnnz::Vector{Ti}
  rowval::Vector{Ti}
  nzval::Vector{Tv}
end

struct CachedCSCAllocation{Tv,Ti} <: AbstractType
  nrows::Int
  ncols::Int
  colptr::Vector{Ti}
  rowval::Vector{Ti}
  nzval::Vector{Tv}
  cache_p::Vector{Ti}
  cache_id::Base.RefValue{Int}
end


Base.eltype(alloc::Union{CSCAllocation, CachedCSCAllocation}) = eltype(alloc.nzval)

function reset!(alloc::CSCAllocation)
    fill!(alloc.nzval,zero(eltype(alloc.nzval)))
    alloc
end

function reset!(alloc::CachedCSCAllocation)
    fill!(alloc.nzval,zero(eltype(alloc.nzval)))
    alloc.cache_id[] = 0
    alloc
end


function contribute!(alloc::CachedCSCAllocation,v,i,j,field_i,field_j)
    @boundscheck @assert length(i) == 1
    @boundscheck @assert length(j) == 1
    @boundscheck @assert field_i == 1
    @boundscheck @assert field_j == 1
    cache_id = alloc.cache_id[] + 1
    p = alloc.cache_p[cache_id]
    alloc.nzval[p] += v 
    alloc.cache_id[] = cache_id
    nothing
end


function contribute!(alloc::CSCAllocation,v,i,j,field_i,field_j)
    @boundscheck @assert length(i) == 1
    @boundscheck @assert length(j) == 1
    @boundscheck @assert field_i == 1
    @boundscheck @assert field_j == 1
    pini = Int(alloc.colptr[j])
    pend = pini + Int(alloc.colnnz[j]) - 1
    p = searchsortedfirst(alloc.rowval,i,pini,pend,Base.Order.Forward)
    if (p>pend)
        # add new entry
        alloc.colnnz[j] += 1
        alloc.rowval[p] = i
        alloc.nzval[p] = v
    elseif alloc.rowval[p] != i
        # shift one forward from p to pend
        @boundscheck @assert pend+1 < Int(alloc.colptr[j+1])
        for k in pend:-1:p
            o = k + 1
            alloc.rowval[o] = alloc.rowval[k]
            alloc.nzval[o] = alloc.nzval[k]
        end
        # add new entry
        alloc.colnnz[j] += 1
        alloc.rowval[p] = i
        alloc.nzval[p] = v
    else 
        # update existing entry
        alloc.nzval[p] += v 
    end
    nothing
end

function compress(alloc::CachedCSCAllocation{Tv,Ti};reuse=Val(true)) where {Tv,Ti}
    A = SparseMatrixCSC(alloc.nrows,alloc.ncols,alloc.colptr,alloc.rowval,alloc.nzval)
    if val_parameter(reuse)
        (A, nothing)
    else
        A
    end
end

function compress(alloc::CSCAllocation{Tv,Ti};reuse=Val(true)) where {Tv,Ti}
    k = 1
    for j in 1:alloc.ncols
        pini = Int(alloc.colptr[j])
        pend = pini + Int(alloc.colnnz[j]) - 1
        for p in pini:pend
        alloc.nzval[k] = alloc.nzval[p]
        alloc.rowval[k] = alloc.rowval[p]
        k += 1
        end
    end
    @inbounds for j in 1:alloc.ncols
        alloc.colptr[j+1] = alloc.colnnz[j]
    end
    length_to_ptrs!(alloc.colptr)
    nnz = alloc.colptr[end]-1
    resize!(alloc.rowval,nnz)
    resize!(alloc.nzval,nnz)

    A = SparseMatrixCSC(alloc.nrows,alloc.ncols,alloc.colptr,alloc.rowval,alloc.nzval)
    if val_parameter(reuse)
        (A, nothing)
    else
        A
    end
end

function compress!(alloc::Union{CSCAllocation, CachedCSCAllocation},A::SparseMatrixCSC,cache)
    if alloc.nzval !== A.nzval
        copyto!(A.nzval,alloc.nzval)
    end
    A
end

##############################################################################


## Element-indexed CSC assembly and direct vector assembly

"""
    element_csc_assembly()

Sequential matrix assembly with domain-local element-indexed CSC positions. Select
`monolithic_assembly(element_csc_assembly())` through the existing
`assembly_method.matrix` option of `allocate_matrix`.
The existing assembly methods and defaults are kept for comparison.
"""
element_csc_assembly() = ElementCSCAssembly()

struct ElementCSCAssembly <: AbstractType end

struct ElementCSCDomainCache{Ti,P<:AbstractMatrix{Ti},O<:AbstractMatrix{Ti},S<:AbstractArray{Ti,3}} <: AbstractType
    positions::P     # Logical [local entry, domain-local entity], in either layout.
    block_offsets::O # Zero-based [row field, column field]; -1 omits a pair.
    block_sizes::S   # [dimension, row field, column field]
end

struct ElementCSCMatrixAllocation{Tv,Ti,C,P<:AbstractVector{Ti},R<:AbstractVector{Ti},V<:AbstractVector{Tv}} <: AbstractType
    nrows::Int
    ncols::Int
    colptr::P
    rowval::R
    nzval::V
    domain_caches::C # Tuple of ElementCSCDomainCache, in allocation argument order.
end

# A numeric argument for one domain's kernel. Only the selected lookup and the
# shared values are reachable; no domain, mesh, other caches, or CSC graph.
struct ElementCSCDomainAllocation{Tv,Ti,V<:AbstractVector{Tv},C<:ElementCSCDomainCache{Ti}} <: AbstractType
    nzval::V
    cache::C
end

Adapt.@adapt_structure ElementCSCDomainCache
Adapt.@adapt_structure ElementCSCMatrixAllocation
Adapt.@adapt_structure ElementCSCDomainAllocation

"""
    to_layout(allocation, layout)

Change only the element CSC position table's physical layout. Its logical indices
remain `[local_entry, element]`. `EntryMajorLayout()` stores a transposed parent with
contiguous elements for each local entry; `ElementMajorLayout()` stores contiguous local
entries for each element. Other buffers are shared. Layout conversion is separate
from device transfer: use `Adapt.adapt(CuArray, to_layout(allocation, EntryMajorLayout()))`
to upload, then select domain arguments from the adapted parent to share values.
"""
function to_layout(cache::ElementCSCDomainCache,::EntryMajorLayout)
    cache.positions isa Transpose && return cache
    positions = transpose(permutedims(cache.positions))
    ElementCSCDomainCache(positions,cache.block_offsets,cache.block_sizes)
end

function to_layout(cache::ElementCSCDomainCache,::ElementMajorLayout)
    cache.positions isa Transpose || return cache
    ElementCSCDomainCache(permutedims(parent(cache.positions)),cache.block_offsets,cache.block_sizes)
end

function to_layout(alloc::ElementCSCMatrixAllocation,layout::AbstractLayout)
    caches = map(cache->to_layout(cache,layout),alloc.domain_caches)
    caches === alloc.domain_caches && return alloc
    ElementCSCMatrixAllocation(alloc.nrows,alloc.ncols,alloc.colptr,alloc.rowval,alloc.nzval,caches)
end

function to_layout(alloc::ElementCSCDomainAllocation,layout::AbstractLayout)
    cache = to_layout(alloc.cache,layout)
    cache === alloc.cache && return alloc
    ElementCSCDomainAllocation(alloc.nzval,cache)
end

function to_layout(alloc::MonolithicAssemblyAllocation,layout::AbstractLayout)
    allocation = to_layout(alloc.allocation,layout)
    allocation === alloc.allocation && return alloc
    MonolithicAssemblyAllocation(allocation,alloc.offsets,alloc.block_mask)
end

function to_layout(alloc::MatrixAllocation,layout::AbstractLayout)
    allocation = to_layout(alloc.allocation,layout)
    allocation === alloc.allocation && return alloc
    MatrixAllocation(allocation,alloc.free_or_dirichlet,alloc.dof_maps)
end

Base.eltype(alloc::Union{ElementCSCMatrixAllocation,ElementCSCDomainAllocation}) = eltype(alloc.nzval)

# Branch from the shared symbolic CSC builder before creating its flat cache.
function counter(::ElementCSCAssembly,::Type{T},dofs_i,dofs_j;kwargs...) where T
    counter(csc_assembly(),T,dofs_i,dofs_j;kwargs...)
end

allocate_matrix_sparsity!(method,counter,args...) = allocate_matrix_barrier!(counter,args...)

function allocate_matrix_sparsity!(method::MonolithicAssembly,counter,args...)
    allocate_matrix_sparsity!(method.method,counter,args...)
end

# These temporary tuples alias existing numeric buffers. They never become
# part of the matrix allocation, and do not create accessor closures or copy DOFs.
function element_csc_dof_arrays(space,domain)
    dofs = face_dofs(space)
    d,D = num_dims(domain),num_dims(space)
    incidence = d == D ? nothing : face_incidence(topology(mesh(space)),d,D)
    selected_sides = d+1 == D ? faces_around(domain) : nothing
    (;dofs=dofs.data,dof_ptrs=dofs.ptrs,domain_faces=faces(domain),selected_sides,
        cells=incidence === nothing ? nothing : incidence.data,
        cell_ptrs=incidence === nothing ? nothing : incidence.ptrs)
end

function element_csc_num_sides(data,face)
    (data.cell_ptrs === nothing || data.selected_sides !== nothing) && return 1
    mesh_face = data.domain_faces[face]
    Int(data.cell_ptrs[mesh_face+1]-data.cell_ptrs[mesh_face])
end

function element_csc_dof_range(data,face,side)
    mesh_face = data.domain_faces[face]
    cell = if data.cell_ptrs === nothing
        mesh_face
    else
        physical_side = data.selected_sides === nothing ? side : data.selected_sides[face]
        data.cells[Int(data.cell_ptrs[mesh_face])+physical_side-1]
    end
    Int(data.dof_ptrs[cell]):Int(data.dof_ptrs[cell+1])-1
end

function allocate_matrix_sparsity!(::ElementCSCAssembly,counter,domains,space_i,space_j,
        fields_i,fields_j,free_i,free_j,map_i,map_j)
    do_loop(counter) || return counter
    for domain in domains
        rows = map(field->element_csc_dof_arrays(GT.field(space_i,field),domain),fields_i)
        columns = map(field->element_csc_dof_arrays(GT.field(space_j,field),domain),fields_j)
        allocate_element_csc_sparsity!(counter,num_faces(domain),rows,columns,
            fields_i,fields_j,free_i,free_j,map_i,map_j)
    end
    counter
end

function allocate_element_csc_sparsity!(counter,nfaces,rows,columns,fields_i,fields_j,
        free_i,free_j,map_i,map_j)
    for face in 1:nfaces, field_j in fields_j, field_i in fields_i
        row_data,column_data = rows[field_i],columns[field_j]
        for side_j in 1:element_csc_num_sides(column_data,face),
            side_i in 1:element_csc_num_sides(row_data,face)
            row_range = element_csc_dof_range(row_data,face,side_i)
            column_range = element_csc_dof_range(column_data,face,side_j)
            for pj in column_range
                column = column_data.dofs[pj]
                skip_dof(column,free_j) && continue
                for pi in row_range
                    row = row_data.dofs[pi]
                    skip_dof(row,free_i) && continue
                    contribute!(counter,nothing,map_i(row),map_j(column),field_i,field_j)
                end
            end
        end
    end
    counter
end

function allocate_matrix_from_sparsity(method,counter,args)
    alloc = allocate(counter)
    allocate_matrix_barrier!(alloc,args...)
    finalize(alloc)
end

function allocate_matrix_from_sparsity(method::MonolithicAssembly,counter,args)
    allocate_matrix_from_sparsity(method.method,counter,args)
end

function allocate_matrix_from_sparsity(::ElementCSCAssembly,counter,args)
    domains,space_i,space_j = args[1],args[2],args[3]
    entity_counts = map(num_faces,domains)
    element_csc_matrix_allocation(counter,num_fields(space_i),num_fields(space_j),entity_counts)
end

prepare_matrix_allocation(method,malloc,space_i,space_j,domains) = malloc

function prepare_matrix_allocation(method::MonolithicAssembly,malloc,space_i,space_j,domains)
    prepare_matrix_allocation(method.method,malloc,space_i,space_j,domains)
end

function prepare_matrix_allocation(::ElementCSCAssembly,malloc,space_i,space_j,domains)
    storage = element_csc_matrix_storage(malloc.allocation)
    caches = map(storage.domain_caches,domains) do cache,domain
        rows = map(space->element_csc_dof_arrays(space,domain),fields(space_i))
        columns = map(space->element_csc_dof_arrays(space,domain),fields(space_j))
        entries_per_element = prepare_element_csc_domain_blocks!(malloc,cache,space_i,space_j,domain,rows,columns)
        nentities = num_faces(domain)
        checked_element_csc_position_count(eltype(cache.positions),0,(entries_per_element,nentities))
        prepared = ElementCSCDomainCache(zeros(eltype(cache.positions),entries_per_element,nentities),
            cache.block_offsets,cache.block_sizes)
        precompute_element_csc_positions!(malloc,prepared,domain,rows,columns)
        prepared
    end
    allocation = with_element_csc_domain_caches(malloc.allocation,caches)
    MatrixAllocation(allocation,malloc.free_or_dirichlet,malloc.dof_maps)
end

function element_csc_matrix_allocation(csc::CachedCSCCounterStep1,ni,nj,entity_counts::Tuple)
    compact_csc_structure!(csc)
    T = eltype(csc.matrix_type)
    Ti = eltype(csc.rowval)
    caches = map(entity_counts) do nentities
        ElementCSCDomainCache(zeros(Ti,0,nentities),fill(-one(Ti),ni,nj),zeros(Ti,4,ni,nj))
    end
    ElementCSCMatrixAllocation(csc.nrows,csc.ncols,csc.colptr,csc.rowval,
        zeros(T,length(csc.rowval)),caches)
end

function element_csc_matrix_allocation(counter::MonolithicAssemblyCounter,ni,nj,entity_counts::Tuple)
    MonolithicAssemblyAllocation(element_csc_matrix_allocation(counter.counter,ni,nj,entity_counts),
        counter.offsets,counter.block_mask)
end

element_csc_matrix_storage(alloc::Union{ElementCSCMatrixAllocation,ElementCSCDomainAllocation}) = alloc
element_csc_matrix_storage(alloc::MonolithicAssemblyAllocation) = element_csc_matrix_storage(alloc.allocation)

function with_element_csc_domain_caches(alloc::ElementCSCMatrixAllocation,caches)
    ElementCSCMatrixAllocation(alloc.nrows,alloc.ncols,alloc.colptr,alloc.rowval,alloc.nzval,caches)
end

function with_element_csc_domain_caches(alloc::MonolithicAssemblyAllocation,caches)
    MonolithicAssemblyAllocation(with_element_csc_domain_caches(alloc.allocation,caches),
        alloc.offsets,alloc.block_mask)
end

is_element_csc_allocation(alloc) = false
is_element_csc_allocation(::Union{ElementCSCMatrixAllocation,ElementCSCDomainAllocation}) = true
is_element_csc_allocation(alloc::MonolithicAssemblyAllocation) = is_element_csc_allocation(alloc.allocation)

"""
    element_csc_domain(alloc::MatrixAllocation, domain_id)

Select the cache for the `domain_id`-th domain supplied to `allocate_matrix`.
The returned MatrixAllocation is a numeric scatter argument sharing the global
values. Use domain-local entity and side indices with it; reset/compress the
parent allocation. The original six-argument contribution uses the parent.
"""
function element_csc_domain(malloc::MatrixAllocation,domain_id::Integer)
    allocation = element_csc_domain_allocation(malloc.allocation,domain_id)
    MatrixAllocation(allocation,malloc.free_or_dirichlet,malloc.dof_maps)
end

function element_csc_domain_allocation(alloc::ElementCSCMatrixAllocation,domain_id)
    ElementCSCDomainAllocation(alloc.nzval,alloc.domain_caches[domain_id])
end

function element_csc_domain_allocation(alloc::MonolithicAssemblyAllocation,domain_id)
    MonolithicAssemblyAllocation(element_csc_domain_allocation(alloc.allocation,domain_id),
        alloc.offsets,alloc.block_mask)
end

# The legacy methods retain the same runner argument and generated call.
function matrix_assembly_contribution_allocation(malloc::MatrixAllocation,domain_id)
    is_element_csc_allocation(malloc.allocation) ? element_csc_domain(malloc,domain_id) : malloc
end

function element_csc_block_enabled(alloc::ElementCSCMatrixAllocation,field_i,field_j)
    field_i == field_j == 1 ||
        throw(ArgumentError("Wrap element_csc_assembly() in monolithic_assembly() for multiple fields"))
    true
end

element_csc_block_enabled(alloc::MonolithicAssemblyAllocation,field_i,field_j) = alloc.block_mask[field_i,field_j]

function element_csc_has_selected_dofs(data,face,free_or_dirichlet)
    for side in 1:element_csc_num_sides(data,face), p in element_csc_dof_range(data,face,side)
        skip_dof(data.dofs[p],free_or_dirichlet) || return true
    end
    false
end

function checked_element_csc_position_count(::Type{Ti},npositions,shape) where Ti
    nblock = 1
    for dimension in shape
        nblock = Base.Checked.checked_mul(nblock,dimension)
    end
    total = Base.Checked.checked_add(npositions,nblock)
    total <= typemax(Ti) || throw(OverflowError(
        "Element CSC position count $total does not fit index_type=$Ti; use index_type=Int64."))
    total
end

function prepare_element_csc_domain_blocks!(malloc,cache,space_i,space_j,domain,rows,columns)
    alloc = malloc.allocation
    spaces_i,spaces_j = fields(space_i),fields(space_j)
    free_i,free_j = malloc.free_or_dirichlet
    entries_per_element = 0
    for field_j in 1:num_fields(space_j), field_i in 1:num_fields(space_i)
        element_csc_block_enabled(alloc,field_i,field_j) || continue
        nr = max_num_reference_dofs(spaces_i[field_i])
        nc = max_num_reference_dofs(spaces_j[field_j])
        # Counts are local to the integration view. A selected physical side 2
        # has one local side; the incidence buffers resolve the physical side.
        nsi = max_num_faces_around(GT.domain(spaces_i[field_i]),domain)
        nsj = max_num_faces_around(GT.domain(spaces_j[field_j]),domain)
        shape = (nr,nc,nsi,nsj)
        for d in 1:4
            cache.block_sizes[d,field_i,field_j] = shape[d]
        end
        # Omit pairs with no accepted contribution anywhere on this domain.
        # Otherwise every entity has the same layout, including zero slots for
        # constrained DOFs. No per-entity offset or compaction map is retained.
        active = any(1:num_faces(domain)) do face
            element_csc_has_selected_dofs(rows[field_i],face,free_i) &&
                element_csc_has_selected_dofs(columns[field_j],face,free_j)
        end
        active || continue
        total = checked_element_csc_position_count(eltype(cache.positions),entries_per_element,shape)
        cache.block_offsets[field_i,field_j] = entries_per_element
        entries_per_element = total
    end
    entries_per_element
end

function precompute_element_csc_positions!(malloc,cache,domain,rows,columns)
    alloc = malloc.allocation
    free_i,free_j = malloc.free_or_dirichlet
    map_i,map_j = malloc.dof_maps
    for field_j in eachindex(columns), field_i in eachindex(rows)
        offset = cache.block_offsets[field_i,field_j]
        offset == -one(offset) && continue
        row_data,column_data = rows[field_i],columns[field_j]
        shape = ntuple(d->Int(cache.block_sizes[d,field_i,field_j]),Val(4))
        indices = LinearIndices(shape)
        for face in 1:num_faces(domain)
            nsi,nsj = element_csc_num_sides(row_data,face),element_csc_num_sides(column_data,face)
            nsi <= shape[3] && nsj <= shape[4] ||
                throw(DimensionMismatch("Incident side counts exceed the precomputed block dimensions"))
            for side_j in 1:nsj, side_i in 1:nsi
                row_range = element_csc_dof_range(row_data,face,side_i)
                column_range = element_csc_dof_range(column_data,face,side_j)
                for (j,pj) in enumerate(column_range), (i,pi) in enumerate(row_range)
                    row,column = row_data.dofs[pi],column_data.dofs[pj]
                    skip_dof(row,free_i) && continue
                    skip_dof(column,free_j) && continue
                    index = indices[i,j,side_i,side_j]
                    cache.positions[offset+index,face] = element_csc_matrix_position(alloc,
                        map_i(row),map_j(column),field_i,field_j)
                end
            end
        end
    end
    malloc
end

function element_csc_matrix_position(alloc::ElementCSCMatrixAllocation,i,j,field_i,field_j)
    @boundscheck @assert field_i == field_j == 1
    pini,pend = Int(alloc.colptr[j]),Int(alloc.colptr[j+1])-1
    p = searchsortedfirst(alloc.rowval,i,pini,pend,Base.Order.Forward)
    p <= pend && alloc.rowval[p] == i ||
        throw(ArgumentError("Contribution is outside the allocated CSC structure"))
    p
end

function element_csc_matrix_position(alloc::MonolithicAssemblyAllocation,i,j,field_i,field_j)
    offsets_i,offsets_j = alloc.offsets
    element_csc_matrix_position(alloc.allocation,i+offsets_i[field_i],j+offsets_j[field_j],1,1)
end

"""
    contribute_element!(alloc::MatrixAllocation, values, field_i, field_j, element, side_i=1, side_j=1)

Scatter a local matrix through element CSC positions without global DOF IDs.
Element and side indices are local to the selected domain. Constrained entries
are encoded by zero positions; field-pair offsets and dimensions are shared.
Values must be a one-based matrix fitting the selected block, with valid element
and side indices. The entry scatter uses `@inbounds` without runtime validation.
All contribution methods mutate their storage and return `nothing`.
"""
@inline function contribute_element!(malloc::MatrixAllocation,vals,field_i,field_j,
        face::Integer,side_i::Integer=1,side_j::Integer=1)
    contribute_element!(malloc.allocation,vals,field_i,field_j,face,side_i,side_j)
    nothing
end

@inline function contribute_element!(alloc::MonolithicAssemblyAllocation,vals,
        field_i,field_j,face,side_i,side_j)
    if alloc.block_mask[field_i,field_j]
        contribute_element!(alloc.allocation,vals,field_i,field_j,face,side_i,side_j)
    end
    nothing
end

@inline function contribute_element!(alloc::ElementCSCMatrixAllocation,vals,
        field_i,field_j,face,side_i,side_j)
    length(alloc.domain_caches) == 1 || throw(ArgumentError(
        "Select a domain cache with element_csc_domain(alloc,domain_id) before indexed scatter."))
    contribute_element_csc!(alloc.nzval,only(alloc.domain_caches),vals,
        field_i,field_j,face,side_i,side_j)
    nothing
end

@inline function contribute_element!(alloc::ElementCSCDomainAllocation,vals,
        field_i,field_j,face,side_i,side_j)
    contribute_element_csc!(alloc.nzval,alloc.cache,vals,
        field_i,field_j,face,side_i,side_j)
    nothing
end

@inline function contribute_element_csc!(nzval,cache,vals,field_i,field_j,face,side_i,side_j)
    base = element_csc_contribution_base(cache,field_i,field_j,side_i,side_j)
    contribute_element_csc_base!(nzval,cache,vals,field_i,field_j,face,base)
    nothing
end

function element_csc_single_domain_cache(alloc::ElementCSCMatrixAllocation)
    length(alloc.domain_caches) == 1 || throw(ArgumentError(
        "Select a domain cache with element_csc_domain(alloc,domain_id) before indexed scatter."))
    only(alloc.domain_caches)
end

element_csc_single_domain_cache(alloc::ElementCSCDomainAllocation) = alloc.cache

# The field/side part is independent of the element. Code generation evaluates
# it from a known prototype; computed-base calls evaluate it at runtime.
function element_csc_contribution_base(cache::ElementCSCDomainCache,field_i,field_j,side_i,side_j)
    offset = cache.block_offsets[field_i,field_j]
    offset == -one(offset) && return -1
    shape = ntuple(d->Int(cache.block_sizes[d,field_i,field_j]),Val(4))
    Int(offset)+LinearIndices(shape)[1,1,side_i,side_j]-1
end

@inline function contribute_element_csc_base!(nzval,cache,vals,field_i,field_j,face,base)
    base == -1 && return nothing
    shape = ntuple(d->Int(cache.block_sizes[d,field_i,field_j]),Val(2))
    contribute_element_csc_block!(nzval,cache.positions,vals,face,base,shape)
    nothing
end

# Generated literal calls carry their selected field's capacity and base, so
# neither field IDs nor a runtime mask/metadata lookup is needed for scatter.
@inline function contribute_element_block!(malloc::MatrixAllocation,vals,face::Integer,
        base::Integer,nrows::Integer,ncols::Integer)
    base == -1 && return nothing
    alloc = element_csc_matrix_storage(malloc.allocation)
    cache = element_csc_single_domain_cache(alloc)
    contribute_element_csc_block!(alloc.nzval,cache.positions,vals,face,Int(base),(Int(nrows),Int(ncols)))
    nothing
end

# Preserve contiguous linear addressing for ordinary CPU matrices. Device arrays
# and transposed tables use logical indexing, specialized and inlined by type.
@inline element_csc_position_base(positions::AbstractMatrix,face,base) = base
@inline element_csc_position_base(positions::Matrix,face,base) = base+(face-1)*size(positions,1)
@inline element_csc_position(positions::AbstractMatrix,entry,face) = @inbounds positions[entry,face]
@inline element_csc_position(positions::Matrix,entry,face) = @inbounds positions[entry]

@inline function contribute_element_csc_block!(nzval,positions,vals,face,base,shape)
    # Generated callers supply valid one-based values and a precomputed layout.
    position_base = element_csc_position_base(positions,face,base)
    indices = LinearIndices(shape)
    @inbounds for j in axes(vals,2), i in axes(vals,1)
        p = element_csc_position(positions,position_base+indices[i,j],face)
        if p != 0
            nzval[p] += vals[i,j]
        end
    end
    nothing
end

"""
    contribute_element_atomic!(alloc::MatrixAllocation, values, field_i, field_j, element, side_i=1, side_j=1)

Element CSC contribution with atomic additions. Uses the same domain-local
indices, block layout, and valid-input contract as `contribute_element!`.
Complete all concurrent contributions before resetting or compressing the matrix.
"""
@inline function contribute_element_atomic!(malloc::MatrixAllocation,vals,field_i,field_j,
        face::Integer,side_i::Integer=1,side_j::Integer=1)
    contribute_element_atomic!(malloc.allocation,vals,field_i,field_j,face,side_i,side_j)
    nothing
end

@inline function contribute_element_atomic!(alloc::MonolithicAssemblyAllocation,vals,
        field_i,field_j,face,side_i,side_j)
    if alloc.block_mask[field_i,field_j]
        contribute_element_atomic!(alloc.allocation,vals,field_i,field_j,face,side_i,side_j)
    end
    nothing
end

@inline function contribute_element_atomic!(alloc::Union{ElementCSCMatrixAllocation,ElementCSCDomainAllocation},vals,
        field_i,field_j,face,side_i,side_j)
    cache = element_csc_single_domain_cache(alloc)
    base = element_csc_contribution_base(cache,field_i,field_j,side_i,side_j)
    base == -1 && return nothing
    shape = ntuple(d->Int(cache.block_sizes[d,field_i,field_j]),Val(2))
    contribute_element_csc_block_atomic!(alloc.nzval,cache.positions,vals,face,base,shape)
    nothing
end

"""
    contribute_element_block_atomic!(alloc::MatrixAllocation, values, element, base, nrows, ncols)

Atomic counterpart of the prepared element CSC block contribution. The base
and row/column capacities use the existing zero-based block layout.
"""
@inline function contribute_element_block_atomic!(malloc::MatrixAllocation,vals,face::Integer,
        base::Integer,nrows::Integer,ncols::Integer)
    base == -1 && return nothing
    alloc = element_csc_matrix_storage(malloc.allocation)
    cache = element_csc_single_domain_cache(alloc)
    contribute_element_csc_block_atomic!(alloc.nzval,cache.positions,vals,face,
        Int(base),(Int(nrows),Int(ncols)))
    nothing
end

# Buffer-only kernel primitive: no allocation wrapper argument or aggregate
# return slot. Callers handle disabled blocks (base == -1). Inlining also
# avoids materializing the array descriptors at a device call boundary.
@inline function contribute_element_csc_block_atomic!(nzval,positions,vals,face,base,shape)
    position_base = element_csc_position_base(positions,face,base)
    indices = LinearIndices(shape)
    @inbounds for j in axes(vals,2), i in axes(vals,1)
        p = element_csc_position(positions,position_base+indices[i,j],face)
        if p != 0
            Atomix.@atomic :monotonic nzval[p] += vals[i,j]
        end
    end
    nothing
end

# The original six-argument coordinate interface remains on the parent.
function contribute!(alloc::ElementCSCMatrixAllocation,v,i,j,field_i,field_j)
    p = element_csc_matrix_position(alloc,i,j,field_i,field_j)
    alloc.nzval[p] += v
    nothing
end

function reset!(alloc::Union{ElementCSCMatrixAllocation,ElementCSCDomainAllocation})
    fill!(alloc.nzval,zero(eltype(alloc)))
    alloc
end

function compress(alloc::ElementCSCMatrixAllocation;reuse=Val(false))
    A = SparseMatrixCSC(alloc.nrows,alloc.ncols,alloc.colptr,alloc.rowval,alloc.nzval)
    val_parameter(reuse) ? (A,nothing) : A
end

function compress!(alloc::ElementCSCMatrixAllocation,A::SparseMatrixCSC,cache)
    if A.nzval !== alloc.nzval
        copyto!(A.nzval,alloc.nzval)
    end
    A
end

"""
    direct_vector_assembly()

Sequential updates directly into a dense vector. Select
`monolithic_assembly(direct_vector_assembly())` with `assembly_method.vector`.
"""
direct_vector_assembly() = DirectVectorAssembly()

struct DirectVectorAssembly <: AbstractType end

struct DirectVectorAllocation{T} <: AbstractType
    vector::Vector{T}
end

Base.eltype(alloc::DirectVectorAllocation) = eltype(alloc.vector)

is_direct_vector_allocation(::Any) = false
is_direct_vector_allocation(::DirectVectorAllocation) = true
is_direct_vector_allocation(alloc::MonolithicAssemblyAllocation) = is_direct_vector_allocation(alloc.allocation)

function counter(::DirectVectorAssembly,::Type{T},dofs_i;
                 eltype=T,index_type=Int32,vector_type=Vector{eltype}) where T
    vector_type <: Vector || throw(ArgumentError("Direct vector assembly requires a dense Vector"))
    DirectVectorAllocation(zeros(Base.eltype(vector_type),length(dofs_i)))
end

do_loop(::DirectVectorAllocation) = false
allocate(alloc::DirectVectorAllocation) = alloc

function contribute!(alloc::DirectVectorAllocation,v,i,field_i)
    @boundscheck @assert field_i == 1
    alloc.vector[i] += v
    nothing
end

"""
    contribute_atomic!(alloc::VectorAllocation, values, dofs, field=1)

Direct vector contribution with atomic additions, retaining the existing DOF
selection, field offsets, and masks. Reset/compress outside concurrent updates.
"""
function contribute_atomic!(malloc::VectorAllocation,vals,dofs_i,field_i=1)
    free_or_diri_i = malloc.free_or_dirichlet
    map_i = malloc.dof_map
    alloc = malloc.allocation
    @inbounds for i in eachindex(dofs_i)
        dof_i = dofs_i[i]
        skip_dof(dof_i,free_or_diri_i) && continue
        contribute_atomic!(alloc,vals[i],map_i(dof_i),field_i)
    end
    nothing
end

function contribute_atomic!(alloc::MonolithicAssemblyAllocation,v,i,field_i)
    if alloc.block_mask[field_i]
        offsets_i, = alloc.offsets
        contribute_atomic!(alloc.allocation,v,offsets_i[field_i]+i,1)
    end
    nothing
end

function contribute_atomic!(alloc::DirectVectorAllocation,v,i,field_i)
    @inbounds Atomix.@atomic :monotonic alloc.vector[i] += v
    nothing
end

function reset!(alloc::DirectVectorAllocation)
    fill!(alloc.vector,zero(eltype(alloc)))
    alloc
end

function compress(alloc::DirectVectorAllocation;reuse=Val(false))
    val_parameter(reuse) ? (alloc.vector,nothing) : alloc.vector
end

function compress!(alloc::DirectVectorAllocation,b::Vector,cache)
    if b !== alloc.vector
        copyto!(b,alloc.vector)
    end
    b
end





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
    ) where T
    free_or_diri_i = free_or_dirichlet
    dofs_i = GT.dofs(space_i,free_or_diri_i)
    method = vector_assembly_method(;assembly_method...)
    options = vector_assembly_options(;assembly_options...)
    counter = GT.counter(method,T,dofs_i;options...)
    dof_map_i = free_or_diri_i == FREE ? identity : -
    fields_i = ntuple(num_fields(space_i))
    for domain in domains
        nfaces = num_faces(domain)
        field_face_dofs_i = map(field_i->dofs_accessor(GT.fields(space_i,field_i)),fields_i)
        allocate_vector_barrier!(counter,nfaces,fields_i,field_face_dofs_i,dof_map_i)
    end
    alloc = allocate(counter)
    VectorAllocation(alloc,free_or_dirichlet,dof_map_i)
end

function allocate_vector_barrier!(counter,nfaces,fields_i,field_face_dofs_i,map_i)
    if ! do_loop(counter)
        return counter
    end
    for face in 1:nfaces
        field_dofs_i = map(field_i->field_face_dofs_i[field_i](face),fields_i)
        for field_i in fields_i
            dofs_i = field_dofs_i[field_i]
            for dof_i in dofs_i
                if skip_dof(dof_i,free_or_diri_i)
                    continue
                end
                contribute!(counter,nothing,dof_i,field_i,map_i)
            end
        end
    end
    counter
end

function allocate_matrix(::Type{T},space_i::AbstractSpace,space_j::AbstractSpace,domains::AbstractDomain...;
        free_or_dirichlet = (FREE,FREE),
        assembly_method = (;)
        assembly_options = (;)
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
    fields_i = ntuple(num_fields(space_i))
    fields_j = ntuple(num_fields(space_j))
    for domain in domains
        nfaces = num_faces(domain)
        field_face_dofs_i = map(field_i->dofs_accessor(GT.fields(space_i,field_i)),fields_i)
        field_face_dofs_j = map(field_j->dofs_accessor(GT.fields(space_j,field_j)),fields_j)
        allocate_matrix_barrier!(counter,nfaces,fields_i,fields_j,field_face_dofs_i,field_face_dofs_j,dof_map_i,dof_map_j)
    end
    alloc = allocate(counter)
    MatrixAllocation(alloc,free_or_dirichlet,dof_maps)
end

function allocate_matrix_barrier!(counter,nfaces,fields_i,fields_j,field_face_dofs_i,field_face_dofs_j,map_i,map_j)
    if ! do_loop(counter)
        return counter
    end
    for face in 1:nfaces
        field_dofs_i = map(field_i->field_face_dofs_i[field_i](face),fields_i)
        field_dofs_j = map(field_j->field_face_dofs_j[field_j](face),fields_j)
        for field_i in fields_i
            dofs_i = field_dofs_i[field_i]
            for field_j in fields_j
                dofs_j = field_dofs_j[field_j]
                for dof_i in dofs_i
                    if skip_dof(dof_i,free_or_diri_i)
                        continue
                    end
                    for dof_j in dofs_j
                        if skip_dof(dof_j,free_or_diri_j)
                            continue
                        end
                        contribute!(counter,nothing,dof_i,dof_j,field_i,field_j,map_i,map_j)
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
    alloc::A
    free_or_dirichlet::B
    dof_maps::C
end

struct VectorAllocation{A,B,C} <: AbstractType
    alloc::A
    free_or_dirichlet::B
    dof_map::C
end

const AssemblyAllocation = Union{VectorAllocation,MatrixAllocation}

function contribute!(malloc::VectorAllocation,vals,dofs_i,field_i=1,map_i=identity)
    free_or_diri_i  = malloc.free_or_dirichlet
    map_i = malloc.dof_map
    alloc = malloc.alloc
    ndofs_i = length(dofs_i)
    for i in 1:ndofs_i
        dof_i = dofs_i[i]
        if skip_dof(dof_i,free_or_diri_i)
            continue
        end
        contribute!(alloc,vals,dof_i,field_i,map_i)
    end
end

function contribute!(malloc::MatrixAllocation,vals,dofs_i,dofs_j,field_i=1,field_j=1,map_i=identity,map_j=identity)
    free_or_diri_i, free_or_diri_j = malloc.free_or_dirichlet
    map_i,map_j = malloc.dof_maps
    alloc = malloc.alloc
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
            contribute!(alloc,vals,dof_i,dof_j,field_i,field_j,map_i,map_j)
        end
    end
end

function reset!(malloc::AssemblyAllocation)
    reset!(malloc.alloc)
    malloc
end

function compress(malloc::AssemblyAllocation;reuse=Val(false))
    compress(malloc.alloc;reuse)
end

function compress!(malloc::AssemblyAllocation,A,A_cache)
    compress!(malloc.alloc,A,A_cache)
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

function counter(s::MonolithicAssembly,::Type{T},dofs_i;block_mask=nothing,kwargs...)
    GT.counter(s.method,T,dofs_i;kwargs...) where T
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
    offsets_j = blocklasts(dofs_j) .- map(length,blocks(dofs_i))
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

function contribute!(counter::MonolithicAssemblyCounter,v,i,field_i,map_i)
    if counter.block_mask[field_i]
        offsets_i, = counter.offsets
        i2 = offsets_i[field_i]+map_i(i)
        contribute!(counter.counter,v,i2,1,identity)
    end
    counter
end

function contribute!(counter::MonolithicAssemblyCounter,v,i,j,field_i,field_j,map_i,map_j)
    if counter.block_mask[field_i,field_j]
        offsets_i,offsets_j = counter.offsets
        i2 = offsets_i[field_i]+map_i(i)
        j2 = offsets_j[field_j]+map_j(j)
        contribute!(counter.counter,v,i2,j2,1,1,identity,identity)
    end
    counter
end

function allocate(mcounter::MonolithicAssemblyCounter)
    (;counter,offsets,block_mask) = mcounter
    alloc = allocate(counter)
    reset!(counter)
    MonolithicAssemblyAllocation(alloc,offsets,block_mask)
end

struct MonolithicAssemblyAllocation{A,B,C} <: AbstractType
    allocation::A
    offsets::B
    block_mask::C
end

function contribute!(alloc::MonolithicAssemblyAllocation,v,i,field_i,map_i)
    if alloc.block_mask[field_i]
        offsets_i, = alloc.offsets
        i2 = offsets_i[field_i]+map_i(i)
        contribute!(alloc.alloc,v,i2,1,identity)
    end
    alloc
end

function contribute!(alloc::MonolithicAssemblyAllocation,v,i,j,field_i,field_j,map_i,map_j)
    if alloc.block_mask[field_i,field_j]
        offsets_i,offsets_j = alloc.offsets
        i2 = offsets_i[field_i]+map_i(i)
        j2 = offsets_j[field_j]+map_j(j)
        contribute!(alloc.alloc,v,i2,j2,1,1,identity,identity)
    end
    alloc
end

function compress(malloc::MonolithicAssemblyAllocation;reuse=Val(false))
    compress(malloc.alloc;reuse)
end

function comperss!(malloc::MonolithicAssemblyAllocation,A,cache)
    compress!(malloc.alloc,A,cache)
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
    eltype=T,index_type=Int32;matrix_type=SparseMatrixCSC{eltype,Int32}) where T
    ni = length(dofs_i)
    nj = length(dofs_j)
    COOMatrixCounter(matrix_type,index_type,ni,nj,Ref(0))
end

function COOVectorCounter{A,B} <: AbstractType
    vector_type::Type{A}
    index_type::Type{B}
    nrows::Int
    nnz::Ref{Int}
end

function COOMatrixCounter{A,B} <: AbstractType
    matrix_type::Type{A}
    index_type::Type{B}
    nrows::Int
    ncols::Int
    nnz::Ref{Int}
end

const COOCounter = Union{COOVectorCounter,COOMatrixCounter}

function do_loop(counter::COOCounter)
    true
end

function reset!(counter::COOCounter)
    counter.nnz[] = 0
    counter
end

function contribute!(counter::COOVectorCounter,v,i,field_i,map_i)
    @boundscheck @assert v === nothing
    @boundscheck @assert length(i) == 1
    @boundscheck @assert field_i == 1
    counter.nnz[] += 1
    counter
end

function contribute!(counter::COOMatrixCounter,v,i,j,field_i,field_j,map_i,map_j)
    @boundscheck @assert v === nothing
    @boundscheck @assert length(i) == 1
    @boundscheck @assert length(j) == 1
    @boundscheck @assert field_i == 1
    @boundscheck @assert field_j == 1
    counter.nnz[] += 1
    counter
end

function allocate(counter::COOVectorCounter)
    Ti = counter.index_type
    T = eltype(counter.vector_type)
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

function reset!(alloc::COOAllocation)
    reset!(alloc.counter)
    alloc
end

function contribute!(counter::COOVectorAllocation,v,i,field_i,map_i)
    @boundscheck @assert length(i) == 1
    @boundscheck @assert field_i == 1
    alloc.counter.nnz[] += 1
    n = alloc.counter.nnz[]
    alloc.I[n] = map_i(i)
    alloc.V[n] = v
    alloc
end

function contribute!(counter::COOMatrixAllocation,v,i,j,field_i,field_j,map_i,map_j)
    @boundscheck @assert length(i) == 1
    @boundscheck @assert length(j) == 1
    @boundscheck @assert field_i == 1
    @boundscheck @assert field_j == 1
    alloc.counter.nnz[] += 1
    n = alloc.counter.nnz[]
    alloc.I[n] = map_i(i)
    alloc.J[n] = map_j(j)
    alloc.V[n] = v
    alloc
end

function compress(alloc::COOVectorAllocation;reuse=Val(false))
    (;I,V,counter) = alloc
    nrows = counter.nrows
    @assert alloc.counter.vector_type <: Vector
    vec = PartitionedArrays.dense_vector(I,V,nrows)
    cache = nothing
    if val_parameter(reuse)
        (A, cache)
    else
        A
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


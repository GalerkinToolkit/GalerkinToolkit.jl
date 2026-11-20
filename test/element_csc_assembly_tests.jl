module ElementCSCAssemblyTests

using Test
using LinearAlgebra
using SparseArrays
using Adapt
import ForwardDiff
import GalerkinToolkit as GT

const method = GT.element_csc_assembly()
# TODO: Re-enable vector assembly when that scope is supported.
# const vector_method = GT.monolithic_assembly(GT.direct_vector_assembly())
const assembly_method = (;matrix=GT.monolithic_assembly(method))

@testset "Entry-major element CSC lookup layout" begin
    for T in (Float32,Float64), Ti in (Int32,Int64), nfaces in (0,3)
        positions = reshape(Ti[mod(i,7) for i in 1:6nfaces],6,nfaces)
        cache = GT.ElementCSCDomainCache(positions,zeros(Ti,1,1),reshape(Ti[2,3,1,1],4,1,1))
        storage = GT.ElementCSCMatrixAllocation(7,1,Ti[1,8],Ti.(1:7),zeros(T,7),(cache,cache))
        alloc = GT.MatrixAllocation(storage,(GT.FREE,GT.FREE),(identity,identity))
        packed = GT.to_layout(alloc,GT.EntryMajorLayout())
        packed_storage = packed.allocation
        @test packed_storage.nzval === storage.nzval
        @test packed_storage.colptr === storage.colptr
        for selected in packed_storage.domain_caches
            @test selected.positions isa Transpose
            @test size(parent(selected.positions)) == (nfaces,6)
            @test selected.positions == positions
            @test selected.block_offsets === cache.block_offsets
            @test selected.block_sizes === cache.block_sizes
            @test GT.to_layout(selected,GT.EntryMajorLayout()) === selected
        end
        @test GT.to_layout(packed,GT.EntryMajorLayout()) === packed
        restored = GT.to_layout(packed,GT.ElementMajorLayout())
        @test only(unique(c.positions for c in restored.allocation.domain_caches)) == positions
        @test all(c->c.positions isa Matrix{Ti},restored.allocation.domain_caches)
        @test GT.to_layout(alloc,GT.ElementMajorLayout()) === alloc
        adapted = adapt(Array,packed)
        @test adapted.allocation isa GT.ElementCSCMatrixAllocation{T,Ti}
        @test adapted.allocation.domain_caches[1].positions isa Transpose
        values = reshape(T.(1:6),2,3)
        expected = zeros(T,7)
        for face in 1:nfaces
            for j in 1:3, i in 1:2
                p = positions[i+2(j-1),face]
                p == 0 || (expected[p] += values[i,j])
            end
            GT.contribute_element!(GT.element_csc_domain(packed,1),values,1,1,face)
        end
        @test storage.nzval == expected
        GT.reset!(packed)
        for face in 1:nfaces
            GT.contribute_element_block_atomic!(GT.element_csc_domain(packed,2),values,face,0,2,3)
        end
        @test storage.nzval == expected
        selected = GT.element_csc_domain(packed,1)
        @test GT.to_layout(selected,GT.ElementMajorLayout()).allocation.cache.positions == positions
        wrapped = GT.MatrixAllocation(GT.MonolithicAssemblyAllocation(storage,
            (Ti[0],Ti[0]),trues(1,1)),alloc.free_or_dirichlet,alloc.dof_maps)
        @test GT.to_layout(wrapped,GT.EntryMajorLayout()).allocation.allocation.nzval === storage.nzval
        @test adapt(Array,wrapped).allocation.offsets == wrapped.allocation.offsets
    end
end

function index_arrays_have_type(storage,Ti)
    eltype(storage.colptr) === Ti && eltype(storage.rowval) === Ti &&
        all(storage.domain_caches) do cache
            all(name->eltype(getfield(cache,name)) === Ti,(:positions,:block_offsets,:block_sizes))
        end
end

function domain_cache(alloc,domain_id=1)
    storage = GT.element_csc_matrix_storage(alloc.allocation)
    storage isa GT.ElementCSCDomainAllocation ? storage.cache : storage.domain_caches[domain_id]
end

function plain_allocation_data(value)
    if value isa Union{GT.AbstractDomain,GT.AbstractMesh,GT.AbstractSpace}
        return false
    elseif value isa AbstractArray
        return eltype(value) <: Number && isbitstype(eltype(value))
    elseif value isa Number
        return isbitstype(typeof(value))
    elseif value isa Function
        return fieldcount(typeof(value)) == 0
    end
    all(name->plain_allocation_data(getfield(value,name)),fieldnames(typeof(value)))
end

# Inspect a dense element block as a temporary view; the returned allocation stores no
# nested arrays, views, domains, meshes, or spaces.
function csc_position_table(alloc,face=1,domain_id=1,field_i=1,field_j=1)
    cache = domain_cache(alloc,domain_id)
    offset = cache.block_offsets[field_i,field_j]
    offset == -one(offset) && return nothing
    shape = Tuple(Int.(cache.block_sizes[:,field_i,field_j]))
    reshape(view(cache.positions,offset+1:offset+prod(shape),face),shape)
end

function add_reference!(expected,values,rows,columns,free_i=GT.FREE,free_j=GT.FREE,
                        row_offset=0,column_offset=0)
    for j in eachindex(columns), i in eachindex(rows)
        row,column = rows[i],columns[j]
        (free_i == GT.FREE ? row > 0 : row < 0) || continue
        (free_j == GT.FREE ? column > 0 : column < 0) || continue
        expected[abs(row)+row_offset,abs(column)+column_offset] += values[i,j]
    end
    expected
end

function synthetic_block(::Type{T},nr,nc,face,field_i=1,field_j=1) where T
    [T(100face+10field_i+field_j+2i+3j) for i in 1:nr, j in 1:nc]
end

# The new allocator must work even when the legacy accessor is unavailable.
struct ArrayOnlySpace{S} <: GT.AbstractSpace
    space::S
end
GT.domain(space::ArrayOnlySpace) = GT.domain(space.space)
GT.workspace(space::ArrayOnlySpace) = GT.workspace(space.space)
GT.reference_spaces(space::ArrayOnlySpace) = GT.reference_spaces(space.space)
GT.dofs_accessor(::ArrayOnlySpace,::GT.AbstractDomain) = error("Legacy DOF accessor called")

@testset "Element CSC reads existing DOF and incidence arrays" begin
    mesh = GT.cartesian_mesh((0,1,0,1),(2,2);simplexify=false)
    omega = GT.interior(mesh)
    # TODO: Re-enable DG and facet incidence cases.
    # skeleton = GT.skeleton(mesh)
    # V = GT.lagrange_space(omega,2;continuous=Val(false))
    V = GT.lagrange_space(omega,2)
    array_space = ArrayOnlySpace(V)
    # second_side = GT.replace_faces_around(skeleton,fill(2,GT.num_faces(skeleton)))
    # for domain in (omega,skeleton,second_side)
    for domain in (omega,)
        data = GT.element_csc_dof_arrays(array_space,domain)
        @test data.dofs === GT.face_dofs(V).data
        @test data.dof_ptrs === GT.face_dofs(V).ptrs
        @test data.domain_faces === GT.faces(domain)
        @test plain_allocation_data(data)
        #= TODO: Re-enable facet incidence checks with facet assembly.
        if domain !== omega
            incidence = GT.face_incidence(GT.topology(mesh),1,2)
            @test data.cells === incidence.data
            @test data.cell_ptrs === incidence.ptrs
        end
        =#
        legacy = GT.dofs_accessor(V,domain)
        for face in 1:GT.num_faces(domain), side in 1:GT.element_csc_num_sides(data,face)
            @test data.dofs[GT.element_csc_dof_range(data,face,side)] == legacy(face,side)
        end
        alloc = GT.allocate_matrix(Float64,array_space,array_space,domain;assembly_method)
        reference = GT.allocate_matrix(Float64,V,V,domain;assembly_method)
        @test GT.compress(alloc) == GT.compress(reference)
        @test domain_cache(alloc).positions == domain_cache(reference).positions
        @test only(domain_cache(alloc).block_offsets) == 0
        @test plain_allocation_data(alloc)
    end
end

@testset "Element CSC branches from the sparsity builder" begin
    pairs = ((3,2),(1,1),(2,2),(3,1),(1,1),(3,2))
    for T in (Float32,Float64), Ti in (Int32,Int), monolithic in (false,true)
        counter = GT.counter(method,T,1:3,1:2;index_type=Ti)
        for (i,j) in pairs
            GT.contribute!(counter,nothing,i,j,1,1)
        end
        step1 = GT.allocate(counter)
        @test step1 isa GT.CachedCSCCounterStep1
        for (i,j) in pairs
            GT.contribute!(step1,nothing,i,j,1,1)
        end
        # This count belongs exclusively to the legacy flat cache. Invalidating
        # it must not affect the element allocator or trigger Step2 allocation.
        step1.cache_count[] = -1
        builder = monolithic ? GT.MonolithicAssemblyCounter(step1,([0],[0]),fill(true,1,1)) : step1
        allocation = GT.element_csc_matrix_allocation(builder,1,1,(1,))
        storage = GT.element_csc_matrix_storage(allocation)
        @test storage isa GT.ElementCSCMatrixAllocation{T,Ti}
        @test index_arrays_have_type(storage,Ti)
        @test storage.colptr == Ti[1,3,5]
        @test storage.rowval == Ti[1,3,2,3]
        @test all(iszero,storage.nzval)
        @test storage.colptr === step1.colptr && storage.rowval === step1.rowval
        @test isempty(only(storage.domain_caches).positions)
        @test plain_allocation_data(allocation)
        if monolithic
            @test allocation isa GT.MonolithicAssemblyAllocation
            @test allocation.offsets === builder.offsets
            @test allocation.block_mask === builder.block_mask
        end
    end
end

@testset "Element CSC index limits" begin
    limit = Int(typemax(Int32))
    shape = (4,4,1,1)
    @test GT.checked_element_csc_position_count(Int32,limit-16,shape) == limit
    @test_throws OverflowError GT.checked_element_csc_position_count(Int32,limit-15,shape)
    @test GT.checked_element_csc_position_count(Int64,limit,shape) == limit+16
    @test_throws OverflowError GT.checked_element_csc_position_count(Int64,typemax(Int)-8,shape)
    @test_throws OverflowError GT.checked_element_csc_position_count(Int64,0,(typemax(Int),2,1,1))
end

@testset "Sequential element-indexed CSC" begin
    mesh = GT.cartesian_mesh((0,1,0,1,0,1),(2,2,2);simplexify=false)
    omega = GT.interior(mesh)
    for k in 1:3, T in (Float32,Float64), Ti in (Int32,Int)
        V = GT.lagrange_space(omega,k)
        alloc = GT.allocate_matrix(T,V,V,omega;assembly_method,
            assembly_options=(;matrix=(;index_type=Ti)))
        dofs = GT.dofs_accessor(V,omega)
        values = [synthetic_block(T,length(dofs(face)),length(dofs(face)),face)
            for face in 1:GT.num_faces(omega)]
        expected = zeros(T,GT.num_free_dofs(V),GT.num_free_dofs(V))
        for face in reverse(1:GT.num_faces(omega))
            ids = dofs(face)
            GT.contribute_element!(alloc,values[face],1,1,face)
            add_reference!(expected,values[face],ids,ids)
        end
        A,cache = GT.compress(alloc;reuse=Val(true))
        storage = GT.element_csc_matrix_storage(alloc.allocation)
        @test alloc isa GT.MatrixAllocation
        @test plain_allocation_data(alloc)
        @test A isa SparseMatrixCSC{T,Ti}
        @test Matrix(A) == expected
        @test !hasfield(typeof(storage),:cache_id)
        @test size(csc_position_table(alloc)) == ((k+1)^3,(k+1)^3,1,1)
        @test length(domain_cache(alloc).positions) == 8*(k+1)^6
        @test domain_cache(alloc).positions isa Matrix{Ti}
        @test size(domain_cache(alloc).positions) == ((k+1)^6,8)
        @test size(domain_cache(alloc).block_offsets) == (1,1)
        @test only(domain_cache(alloc).block_offsets) == 0
        @test index_arrays_have_type(storage,Ti)
        @test storage.nzval === A.nzval && storage.rowval === A.rowval && storage.colptr === A.colptr
        colptr,rowval = copy(A.colptr),copy(A.rowval)
        GT.reset!(alloc)
        for face in 1:GT.num_faces(omega)
            GT.contribute_element!(alloc,2values[face],1,1,face,1,1)
        end
        @test GT.compress!(alloc,A,cache) === A
        @test Matrix(A) == 2expected
        @test A.colptr == colptr && A.rowval == rowval

        # Calls without element identity remain available for handwritten code.
        GT.reset!(alloc)
        for face in reverse(1:GT.num_faces(omega))
            ids = dofs(face)
            GT.contribute!(alloc,values[face],ids,ids)
        end
        @test Matrix(GT.compress(alloc)) == expected
        other = copy(A)
        fill!(other.nzval,0)
        GT.compress!(alloc,other,nothing)
        @test other == A
    end
end

@testset "Element CSC block capacity and padded stride" begin
    # A partial rectangular block must retain its cached row capacity and the
    # full element-column stride. The second value has a skipped position.
    positions = reshape(Int32.(1:90),30,3)
    positions[5,2] = 0
    cache = GT.ElementCSCDomainCache(positions,zeros(Int32,1,1),
        reshape(Int32[4,5,1,1],4,1,1))
    nzval = zeros(90)
    argument = GT.MatrixAllocation(GT.ElementCSCDomainAllocation(nzval,cache),
        (GT.FREE,GT.FREE),(identity,identity))
    values = reshape(Float64.(1:6),2,3)
    expected = zeros(90)
    expected[[34,38,39,42,43]] = [1,3,4,5,6]
    @test GT.contribute_element_block!(argument,values,2,3,4,5) === nothing
    @test nzval == expected
    for empty_values in (zeros(0,3),zeros(2,0))
        GT.contribute_element_block!(argument,empty_values,2,3,4,5)
        @test nzval == expected
    end
end

@testset "Constrained matrix DOFs" begin
    mesh = GT.cartesian_mesh((0,1,0,1),(2,2);simplexify=false)
    omega = GT.interior(mesh)
    V = GT.lagrange_space(omega,2;dirichlet_boundary=GT.boundary(mesh))
    face_dofs = GT.dofs_accessor(V,omega)
    for free_i in (GT.FREE,GT.DIRICHLET), free_j in (GT.FREE,GT.DIRICHLET)
        alloc = GT.allocate_matrix(Float64,V,V,omega;assembly_method,free_or_dirichlet=(free_i,free_j))
        expected = zeros(length(GT.dofs(V,free_i)),length(GT.dofs(V,free_j)))
        for face in reverse(1:GT.num_faces(omega))
            ids = face_dofs(face)
            values = synthetic_block(Float64,length(ids),length(ids),face)
            GT.contribute_element!(alloc,values,1,1,face)
            add_reference!(expected,values,ids,ids,free_i,free_j)
        end
        @test Matrix(GT.compress(alloc)) == expected
        @test any(iszero,domain_cache(alloc).positions)
    end
    #= TODO: Re-enable direct vector assembly.
    for T in (Float32,Float64), free_i in (GT.FREE,GT.DIRICHLET)
        alloc = GT.allocate_vector(T,V,omega;assembly_method,free_or_dirichlet=free_i)
        expected = zeros(T,length(GT.dofs(V,free_i)))
        for face in reverse(1:GT.num_faces(omega))
            ids = face_dofs(face)
            values = T.(face .+ eachindex(ids))
            GT.contribute!(alloc,values,ids)
            for i in eachindex(ids)
                (free_i == GT.FREE ? ids[i] > 0 : ids[i] < 0) || continue
                expected[abs(ids[i])] += values[i]
            end
        end
        b,cache = GT.compress(alloc;reuse=Val(true))
        @test b == expected
        @test eltype(b) == T
        other = similar(b)
        @test GT.compress!(alloc,other,cache) === other
        @test other == expected
        GT.reset!(alloc)
        @test all(iszero,b)
        GT.contribute!(alloc,T.(1:length(face_dofs(1))),face_dofs(1))
        @test any(!iszero,b)
    end
    =#
end

#= TODO: Re-enable multi-field matrix and vector assembly.
@testset "Field blocks, rectangular matrices, and masks" begin
    mesh = GT.cartesian_mesh((0,1,0,1),(2,1);boundary=false,simplexify=false)
    omega = GT.interior(mesh)
    V1,V2 = GT.lagrange_space(omega,1),GT.lagrange_space(omega,2)
    rows,columns = GT.cartesian_product(V2,V1),GT.cartesian_product(V1,V2)
    mask = Bool[true false; true true]
    alloc = GT.allocate_matrix(Float64,rows,columns,omega;assembly_method,
        assembly_options=(;matrix=(;block_mask=mask)))
    row_spaces,column_spaces = (V2,V1),(V1,V2)
    row_offsets,column_offsets = (0,GT.num_free_dofs(V2)),(0,GT.num_free_dofs(V1))
    expected = zeros(GT.num_free_dofs(rows),GT.num_free_dofs(columns))
    for face in reverse(1:GT.num_faces(omega)), field_i in 2:-1:1, field_j in 2:-1:1
        ids_i = GT.dofs_accessor(row_spaces[field_i],omega)(face)
        ids_j = GT.dofs_accessor(column_spaces[field_j],omega)(face)
        values = synthetic_block(Float64,length(ids_i),length(ids_j),face,field_i,field_j)
        GT.contribute_element!(alloc,values,field_i,field_j,face)
        mask[field_i,field_j] || continue
        add_reference!(expected,values,ids_i,ids_j,GT.FREE,GT.FREE,
            row_offsets[field_i],column_offsets[field_j])
    end
    A = GT.compress(alloc)
    @test alloc isa GT.MatrixAllocation
    @test plain_allocation_data(alloc)
    @test alloc.allocation isa GT.MonolithicAssemblyAllocation
    @test alloc.allocation.allocation isa GT.ElementCSCMatrixAllocation
    @test alloc.allocation.block_mask == mask
    @test !hasfield(typeof(alloc.allocation.allocation),:row_offsets)
    @test !hasfield(typeof(alloc.allocation.allocation),:block_mask)
    @test Matrix(A) == expected
    @test csc_position_table(alloc,1,1,1,2) === nothing
    # Field pairs have different shapes; pack only enabled pairs in each column.
    @test size(domain_cache(alloc).positions) == (9*4+4*4+4*9,2)
    @test domain_cache(alloc).block_offsets == [0 -1; 36 52]
    @test size(domain_cache(alloc).block_sizes) == (4,2,2)
    @test nnz(A[1:GT.num_free_dofs(V2),GT.num_free_dofs(V1)+1:end]) == 0
    for indexed in (false,true)
        GT.reset!(alloc)
        ids_i = GT.dofs_accessor(V2,omega)(1)
        ids_j = GT.dofs_accessor(V1,omega)(1)
        values = ones(length(ids_i),length(ids_j))
        if indexed
            GT.contribute_element!(alloc,values,1,1,1)
        else
            GT.contribute!(alloc,values,ids_i,ids_j,1,1)
        end
        @test sum(A) == length(values)
    end
    vector_alloc = GT.allocate_vector(Float64,rows,omega;assembly_method,
        assembly_options=(;vector=(;block_mask=Bool[false,true])))
    ids_i = GT.dofs_accessor(V2,omega)(1)
    ids_j = GT.dofs_accessor(V1,omega)(1)
    GT.contribute!(vector_alloc,ones(length(ids_i)),ids_i,1)
    GT.contribute!(vector_alloc,ones(length(ids_j)),ids_j,2)
    b = GT.compress(vector_alloc)
    @test vector_alloc isa GT.VectorAllocation
    @test plain_allocation_data(vector_alloc)
    @test vector_alloc.allocation isa GT.MonolithicAssemblyAllocation
    @test vector_alloc.allocation.allocation isa GT.DirectVectorAllocation
    @test all(iszero,b[1:GT.num_free_dofs(V2)])
    @test sum(b) == length(ids_j)
end

=#
@testset "Domain-local CG volume caches" begin
    mesh = GT.cartesian_mesh((0,1,0,1),(2,2);simplexify=false)
    omega = GT.interior(mesh)
    GT.group_faces(mesh,2)["selected elements"] = Int32[4,2]
    GT.group_faces(mesh,2)["empty elements"] = Int32[]
    selected = GT.interior(mesh;group_names=["selected elements"])
    empty_domain = GT.interior(mesh;group_names=["empty elements"])
    V = GT.lagrange_space(omega,1)
    alloc = GT.allocate_matrix(Float64,V,V,selected,omega,empty_domain;assembly_method)
    expected = zeros(GT.num_free_dofs(V),GT.num_free_dofs(V))
    @test plain_allocation_data(alloc)
    for (domain_id,domain) in ((1,selected),(2,omega),(1,selected))
        domain_alloc = GT.element_csc_domain(alloc,domain_id)
        @test plain_allocation_data(domain_alloc)
        @test GT.element_csc_matrix_storage(domain_alloc.allocation).nzval ===
            GT.element_csc_matrix_storage(alloc.allocation).nzval
        face_dofs = GT.dofs_accessor(V,domain)
        for face in reverse(1:GT.num_faces(domain))
            ids = face_dofs(face)
            values = synthetic_block(Float64,length(ids),length(ids),face)
            GT.contribute_element!(domain_alloc,values,1,1,face)
            add_reference!(expected,values,ids,ids)
        end
    end
    @test Matrix(GT.compress(alloc)) == expected
    storage = GT.element_csc_matrix_storage(alloc.allocation)
    @test !hasfield(typeof(storage),:entity_ids)
    @test !hasfield(typeof(storage),:dimension_ptrs)
    @test map(c->length(c.positions),storage.domain_caches) == (2*4^2,4*4^2,0)
    @test domain_cache(alloc,1).block_offsets == fill(Int32(0),1,1)
    @test size(domain_cache(alloc,1).positions) == (16,2)
    @test domain_cache(alloc,2).block_offsets == fill(Int32(0),1,1)
    @test size(domain_cache(alloc,2).positions) == (16,4)
    @test size(domain_cache(alloc,1).block_sizes) == (4,1,1)
    @test_throws ArgumentError GT.contribute_element!(alloc,ones(4,4),1,1,1)
    empty_alloc = GT.allocate_matrix(Float64,V,V,empty_domain;assembly_method)
    @test size(domain_cache(empty_alloc).positions,2) == 0
    @test isempty(domain_cache(empty_alloc).positions)
    @test all(iszero,GT.compress(empty_alloc).nzval)
    @test plain_allocation_data(empty_alloc)

    # Even identical domain arguments own independent lookup buffers while
    # sharing the matrix values. No numerical contribution is deduplicated.
    repeated = GT.allocate_matrix(Float64,V,V,selected,selected;assembly_method)
    @test domain_cache(repeated,1).positions !== domain_cache(repeated,2).positions
    @test domain_cache(repeated,1).positions == domain_cache(repeated,2).positions
    @test GT.element_csc_matrix_storage(GT.element_csc_domain(repeated,1).allocation).nzval ===
        GT.element_csc_matrix_storage(GT.element_csc_domain(repeated,2).allocation).nzval

    #= TODO: Re-enable DG and selected facet sides.
    skeleton = GT.skeleton(mesh)
    Vdg = GT.lagrange_space(omega,1;continuous=Val(false))
    first_side = GT.replace_faces_around(skeleton,fill(1,GT.num_faces(skeleton)))
    second_side = GT.replace_faces_around(skeleton,fill(2,GT.num_faces(skeleton)))
    side_alloc = GT.allocate_matrix(Float64,Vdg,Vdg,second_side;assembly_method)
    side_positions = csc_position_table(side_alloc)
    @test size(side_positions) == (4,4,1,1)
    @test all(>(0),side_positions)
    @test length(domain_cache(side_alloc).positions) == GT.num_faces(skeleton)*4^2

    GT.group_faces(mesh,1)["partial skeleton"] = Int32[first(GT.faces(skeleton))]
    partial = GT.skeleton(mesh;group_names=["partial skeleton"])
    for domains in ((second_side,first_side,skeleton),(second_side,partial))
        alloc = GT.allocate_matrix(Float64,Vdg,Vdg,domains...;assembly_method)
        @test plain_allocation_data(alloc)
        expected = zeros(GT.num_free_dofs(Vdg),GT.num_free_dofs(Vdg))
        # Reverse domain/element/side order to exclude traversal-counter dependence.
        for domain_id in reverse(eachindex(domains))
            domain = domains[domain_id]
            domain_alloc = GT.element_csc_domain(alloc,domain_id)
            numeric = GT.element_csc_matrix_storage(domain_alloc.allocation)
            @test !hasfield(typeof(numeric),:colptr)
            @test !hasfield(typeof(numeric),:domain_caches)
            face_dofs = GT.dofs_accessor(Vdg,domain)
            nsides = GT.max_num_faces_around(omega,domain)
            @test length(domain_cache(alloc,domain_id).positions) == GT.num_faces(domain)*4^2*nsides^2
            for face in reverse(1:GT.num_faces(domain)), side_i in nsides:-1:1, side_j in nsides:-1:1
                ids_i,ids_j = face_dofs(face,side_i),face_dofs(face,side_j)
                physical_i = GT.faces_around(domain) === nothing ? side_i : GT.faces_around(domain)[face]
                physical_j = GT.faces_around(domain) === nothing ? side_j : GT.faces_around(domain)[face]
                values = synthetic_block(Float64,length(ids_i),length(ids_j),face,physical_i,physical_j)
                GT.contribute_element!(domain_alloc,values,1,1,face,side_i,side_j)
                add_reference!(expected,values,ids_i,ids_j)
            end
        end
        @test Matrix(GT.compress(alloc)) == expected
    end
    =#
end

#= TODO: Re-enable multi-field constrained offsets.
@testset "Multi-field constrained offsets" begin
    mesh = GT.cartesian_mesh((0,1,0,1),(2,2);simplexify=false)
    omega = GT.interior(mesh)
    boundary = GT.boundary(mesh)
    V1 = GT.lagrange_space(omega,1;dirichlet_boundary=boundary)
    V2 = GT.lagrange_space(omega,2;dirichlet_boundary=boundary)
    row_spaces,column_spaces = (V1,V2),(V2,V1)
    rows,columns = GT.cartesian_product(row_spaces...),GT.cartesian_product(column_spaces...)
    mask = Bool[true true; false true]
    for free_i in (GT.FREE,GT.DIRICHLET), free_j in (GT.FREE,GT.DIRICHLET)
        alloc = GT.allocate_matrix(Float64,rows,columns,omega;assembly_method,
            free_or_dirichlet=(free_i,free_j),assembly_options=(;matrix=(;block_mask=mask)))
        expected = zeros(length(GT.dofs(rows,free_i)),length(GT.dofs(columns,free_j)))
        row_offsets = (0,length(GT.dofs(V1,free_i)))
        column_offsets = (0,length(GT.dofs(V2,free_j)))
        for face in 1:GT.num_faces(omega), field_j in 1:2, field_i in 1:2
            ids_i = GT.dofs_accessor(row_spaces[field_i],omega)(face)
            ids_j = GT.dofs_accessor(column_spaces[field_j],omega)(face)
            values = synthetic_block(Float64,length(ids_i),length(ids_j),face,field_i,field_j)
            GT.contribute_element!(alloc,values,field_i,field_j,face)
            mask[field_i,field_j] || continue
            add_reference!(expected,values,ids_i,ids_j,free_i,free_j,
                row_offsets[field_i],column_offsets[field_j])
        end
        @test Matrix(GT.compress(alloc)) == expected
        @test plain_allocation_data(alloc)
    end
    for free_i in (GT.FREE,GT.DIRICHLET)
        alloc = GT.allocate_vector(Float64,rows,omega;assembly_method,free_or_dirichlet=free_i)
        expected = zeros(length(GT.dofs(rows,free_i)))
        offsets = (0,length(GT.dofs(V1,free_i)))
        for face in 1:GT.num_faces(omega), field_i in 2:-1:1
            ids = GT.dofs_accessor(row_spaces[field_i],omega)(face)
            values = Float64.(100field_i .+ 10face .+ eachindex(ids))
            GT.contribute!(alloc,values,ids,field_i)
            for i in eachindex(ids)
                (free_i == GT.FREE ? ids[i] > 0 : ids[i] < 0) || continue
                expected[offsets[field_i]+abs(ids[i])] += values[i]
            end
        end
        @test GT.compress(alloc) == expected
        GT.reset!(alloc)
        @test all(iszero,GT.compress(alloc))
    end
end

=#
@testset "Strong Dirichlet volume blocks" begin
    mesh = GT.cartesian_mesh((0,1,0,1,0,1),(4,4,4);simplexify=false)
    omega,boundary = GT.interior(mesh),GT.boundary(mesh)
    V = GT.lagrange_space(omega,1;dirichlet_boundary=boundary)
    face_dofs = GT.dofs_accessor(V,omega)
    alloc = GT.allocate_matrix(Float64,V,V,omega;assembly_method,
        free_or_dirichlet=(GT.FREE,GT.DIRICHLET))
    expected = zeros(GT.num_free_dofs(V),GT.num_dirichlet_dofs(V))
    inactive = Int[]
    for face in 1:GT.num_faces(omega)
        ids = face_dofs(face)
        values = synthetic_block(Float64,length(ids),length(ids),face)
        GT.contribute_element!(alloc,values,1,1,face)
        add_reference!(expected,values,ids,ids,GT.FREE,GT.DIRICHLET)
        all(>(0),ids) && push!(inactive,face)
    end
    @test length(inactive) == 8
    @test size(domain_cache(alloc).positions) == (8^2,4^3)
    @test all(iszero,domain_cache(alloc).positions[:,inactive])
    @test size(domain_cache(alloc).block_offsets) == (1,1)
    @test Matrix(GT.compress(alloc)) == expected
    GT.reset!(alloc)
    for face in inactive
        ids = face_dofs(face)
        GT.contribute_element!(alloc,ones(length(ids),length(ids)),1,1,face)
    end
    @test all(iszero,GT.compress(alloc).nzval)

    # With no strongly constrained DOFs, the free/Dirichlet block has no columns.
    W = GT.lagrange_space(omega,1)
    empty_columns = GT.allocate_matrix(Float64,W,W,omega;assembly_method,
        free_or_dirichlet=(GT.FREE,GT.DIRICHLET))
    @test size(GT.compress(empty_columns)) == (GT.num_free_dofs(W),0)
    @test isempty(GT.compress(empty_columns).nzval)
    @test isempty(domain_cache(empty_columns).positions)
    @test only(domain_cache(empty_columns).block_offsets) == -1

    #= TODO: Re-enable weak boundary scatter.
    # Nitsche uses free DOFs on both domains; only the supplied facets contribute
    # to a boundary-only allocation. Integration itself remains outside this test.
    for domains in ((boundary,),(omega,boundary))
        alloc = GT.allocate_matrix(Float64,W,W,domains...;assembly_method)
        expected = zeros(GT.num_free_dofs(W),GT.num_free_dofs(W))
        for (domain_index,domain) in enumerate(domains)
            domain_alloc = GT.element_csc_domain(alloc,domain_index)
            face_dofs = GT.dofs_accessor(W,domain)
            for face in reverse(1:GT.num_faces(domain))
                ids = face_dofs(face)
                values = synthetic_block(Float64,length(ids),length(ids),face,domain_index)
                GT.contribute_element!(domain_alloc,values,1,1,face)
                add_reference!(expected,values,ids,ids)
            end
            @test size(domain_cache(alloc,domain_index).positions,2) == GT.num_faces(domain)
        end
        @test Matrix(GT.compress(alloc)) == expected
        @test plain_allocation_data(alloc)
    end
    =#
end

@testset "Standard assembly interface and legacy selection" begin
    mesh = GT.cartesian_mesh((0,1,0,1),(2,2);simplexify=false)
    omega = GT.interior(mesh)
    V = GT.lagrange_space(omega,1)
    measure = GT.measure(omega,2)
    a = (u,v) -> GT.∫(x->u(x)*v(x),measure)
    # TODO: Re-enable load vectors.
    # l = v -> GT.∫(x->v(x),measure)
    expected = nothing
    for selected in (GT.element_csc_assembly(),GT.monolithic_assembly(GT.element_csc_assembly()),
                     GT.monolithic_assembly(GT.csc_assembly()),GT.monolithic_assembly(GT.coo_assembly()))
        A,cache = GT.assemble_matrix(a,Float64,V,V;reuse=Val(true),
            assembly_method=(;matrix=selected))
        @test cache.alloc isa GT.MatrixAllocation
        @test sum(A) ≈ 1
        if expected === nothing
            expected = Matrix(A)
        else
            @test Matrix(A) ≈ expected
        end
        fill!(A.nzval,-1)
        @test GT.update_matrix!(A,cache) === A
        @test Matrix(A) ≈ expected
    end
    #= TODO: Re-enable the vector assembly interface.
    for selected in (GT.direct_vector_assembly(),GT.monolithic_assembly(GT.direct_vector_assembly()),
                     GT.monolithic_assembly(GT.coo_assembly()))
        b,cache = GT.assemble_vector(l,Float64,V;reuse=Val(true),
            assembly_method=(;vector=selected))
        @test cache.alloc isa GT.VectorAllocation
        @test sum(b) ≈ 1
        expected_b = copy(b)
        fill!(b,-1)
        @test GT.update_vector!(b,cache) === b
        @test b ≈ expected_b
    end
    =#
end

@testset "Generated domain-local IDs across overlapping domains" begin
    mesh = GT.cartesian_mesh((0,1,0,1),(2,2);simplexify=false)
    omega = GT.interior(mesh)
    GT.group_faces(mesh,2)["reordered cells"] = Int32[4,2]
    selected = GT.interior(mesh;group_names=["reordered cells"])
    V = GT.lagrange_space(omega,2)
    full_measure,selected_measure = GT.measure(omega,4),GT.measure(selected,4)
    a = (u,v) -> GT.∫(x->u(x)*v(x),full_measure)+
                 GT.∫(x->2u(x)*v(x),selected_measure)
    A,cache = GT.assemble_matrix(a,Float64,V,V;assembly_method,reuse=Val(true))
    reference = GT.assemble_matrix(a,Float64,V,V;
        assembly_method=(;matrix=GT.monolithic_assembly(GT.csc_assembly())))
    @test A ≈ reference rtol=1e-12 atol=1e-13
    @test sum(A) ≈ 2
    @test length(cache.loops) == 2
    storage = GT.element_csc_matrix_storage(cache.alloc.allocation)
    @test all(runner->GT.element_csc_matrix_storage(runner.argument.allocation).nzval === A.nzval,cache.loops)
    @test all(i->GT.element_csc_matrix_storage(cache.loops[i].argument.allocation).cache === storage.domain_caches[i],eachindex(cache.loops))
    @test sort(collect(map(c->length(c.positions),storage.domain_caches))) == [2*9^2,4*9^2]
    @test plain_allocation_data(cache.alloc)
    GT.update_matrix!(A,cache)
    @test A ≈ reference rtol=1e-12 atol=1e-13

    # Reassembly must consume the cached positions, not the coordinate fallback.
    saved_positions = map(c->copy(c.positions),storage.domain_caches)
    foreach(c->fill!(c.positions,0),storage.domain_caches)
    GT.update_matrix!(A,cache)
    @test all(iszero,A.nzval)
    foreach((c,p)->copyto!(c.positions,p),storage.domain_caches,saved_positions)
    GT.update_matrix!(A,cache)
    @test A ≈ reference rtol=1e-12 atol=1e-13

    # Geometry/DOFs use mesh IDs 2 and 4; lookup addressing must use local IDs 1 and 2.
    subset_form = (u,v) -> GT.∫(x->u(x)*v(x),selected_measure)
    As,subset_cache = GT.assemble_matrix(subset_form,Float64,V,V;assembly_method,reuse=Val(true))
    expected = GT.assemble_matrix(subset_form,Float64,V,V;
        assembly_method=(;matrix=GT.monolithic_assembly(GT.csc_assembly())))
    @test As ≈ expected rtol=1e-12 atol=1e-13
    @test size(domain_cache(subset_cache.alloc).positions,2) == 2
    @test sum(As) ≈ 0.5
    GT.update_matrix!(As,subset_cache)
    @test As ≈ expected rtol=1e-12 atol=1e-13

    #= TODO: Re-enable generated multi-field DG assembly.
    # The generated path also preserves field offsets and disabled blocks.
    V1 = GT.lagrange_space(omega,1;continuous=Val(false))
    V2 = GT.lagrange_space(omega,2;continuous=Val(false))
    mixed = GT.cartesian_product(V1,V2)
    mixed_form = ((u1,u2),(v1,v2)) -> GT.∫(x->(u1(x)+u2(x))*(v1(x)+v2(x)),full_measure)
    mask_options = (;matrix=(;block_mask=Bool[true false; true true]))
    Am,mixed_cache = GT.assemble_matrix(mixed_form,Float64,mixed,mixed;
        assembly_method,assembly_options=mask_options,reuse=Val(true))
    expected_mixed = GT.assemble_matrix(mixed_form,Float64,mixed,mixed;
        assembly_method=(;matrix=GT.monolithic_assembly(GT.csc_assembly())),assembly_options=mask_options)
    @test Am ≈ expected_mixed rtol=1e-12 atol=1e-13
    @test plain_allocation_data(mixed_cache.alloc)
    GT.update_matrix!(Am,mixed_cache)
    @test Am ≈ expected_mixed rtol=1e-12 atol=1e-13
    =#
end

# Capture completed local matrices from the existing generated integral kernel.
function contribution_calls(expr,name)
    calls = Expr[]
    function visit(node)
        node isa Expr || return
        if node.head == :call && string(first(node.args)) == name
            push!(calls,node)
        end
        foreach(visit,node.args)
    end
    visit(expr)
    calls
end

# TODO: Re-enable the DG contribution-lowering fixture.
# struct OpaqueElementCSCPrototype end
# GT.is_element_csc_allocation(::OpaqueElementCSCPrototype) = true

@testset "Backend-specific contribution expressions" begin
    mesh = GT.cartesian_mesh((0,1,0,1),(2,2);simplexify=false)
    omega = GT.interior(mesh)
    V = GT.lagrange_space(omega,1)
    Q = GT.measure(omega,2)
    a = (u,v)->GT.∫(x->u(x)*v(x),Q)
    contribution = only(GT.contributions(a(GT.form_argument_quantity(V,2),GT.form_argument_quantity(V,1))))
    for backend in (GT.element_csc_assembly(),GT.csc_assembly(),GT.csc_assembly(Val(false)),GT.coo_assembly())
        alloc = GT.allocate_matrix(Float64,V,V,omega;assembly_method=(;matrix=GT.monolithic_assembly(backend)))
        argument = GT.matrix_assembly_contribution_allocation(alloc,1)
        term = GT.optimize(GT.write_assemble_matrix(contribution,V,V;allocation=argument))
        captured,data = GT.capture(GT.parametrize(term))
        expr = GT.expression(captured.body.body)
        if backend isa GT.ElementCSCAssembly
            calls = contribution_calls(expr,"GT.contribute_element_block!")
            @test !isempty(calls)
            @test all(call->length(call.args)==7,calls)
            @test isempty(contribution_calls(expr,"GT.contribute_element!"))
            @test isempty(contribution_calls(expr,"GT.element_csc_contribution_base"))
            @test all(call->all(x->x isa Integer,call.args[end-2:end]),calls)
            @test isempty(contribution_calls(expr,"GT.contribute!"))
            @test !occursin("DOFData",string(expr))
        else
            calls = contribution_calls(expr,"GT.contribute!")
            @test !isempty(calls)
            @test all(call->length(call.args)==7,calls)
            @test isempty(contribution_calls(expr,"GT.contribute_element!"))
            @test isempty(contribution_calls(expr,"GT.element_csc_contribution_base"))
            @test occursin("DOFData",string(expr))
        end
    end

    #= TODO: Re-enable DG facet and multi-field contribution lowering.
    # Test the global DG template with synthetic local values, independently of
    # the unfinished DG integral machinery. Side loops must become four literals.
    skeleton = GT.skeleton(mesh)
    Vdg = GT.lagrange_space(omega,1;continuous=Val(false))
    alloc = GT.allocate_matrix(Float64,Vdg,Vdg,skeleton;assembly_method)
    bindings = (:element,:point,:field_j,:field_i,:side_j,:side_i,:j,:i)
    dependencies = (:U,:V,:Q,:alloc)
    values = :(be[side_i,side_j,field_i,field_j])
    call = GT.matrix_contribution_expr(GT.element_csc_domain(alloc,1),values,
        (Vdg,Vdg,GT.measure(skeleton,2)),dependencies,bindings)
    body = GT.generate_matrix_assembly_template(:(1.0),[2],[2],dependencies,bindings;
        npoints=1,field_ndofs_trial=[4],field_ndofs_test=[4],scalar_type=Float64,contribution=call)
    expr = GT.ast_optimize_with_options(body,Set([:point]),nothing)
    calls = contribution_calls(expr,"GT.contribute_element!")
    @test length(calls) == 4
    @test Set((call.args[end-1],call.args[end]) for call in calls) == Set(((1,1),(2,1),(1,2),(2,2)))
    @test !occursin("DOFData",string(expr))
    # Evaluate the scatter-only generated loop with an allocation argument,
    # ensuring preparation is outside the element loop and remains allocation-free.
    blocks = ntuple(i->fill(Float64(i),4,4),4)
    raw = quote
        for e in 1:ne
            GT.contribute_element!(alloc,blocks[1],1,1,e,1,1)
            GT.contribute_element!(alloc,blocks[2],1,1,e,2,1)
            GT.contribute_element!(alloc,blocks[3],1,1,e,1,2)
            GT.contribute_element!(alloc,blocks[4],1,1,e,2,2)
        end
    end
    argument = GT.element_csc_domain(alloc,1)
    prepared = GT.substitute_element_csc_contribution_bases(raw,argument)
    literal_calls = contribution_calls(prepared,"GT.contribute_element_block!")
    @test length(literal_calls) == 4
    @test [call.args[end-2] for call in literal_calls] == [0,16,32,48]
    @test all(call->call.args[end-1:end] == [4,4],literal_calls)
    @test isempty(contribution_calls(prepared,"GT.element_csc_contribution_base"))
    @test isempty(contribution_calls(prepared,"GT.contribute_element!"))
    @test GT.substitute_element_csc_contribution_bases(raw) == raw
    @test GT.substitute_element_csc_contribution_bases(raw,:alloc) == raw
    opaque = GT.MatrixAllocation(OpaqueElementCSCPrototype(),argument.free_or_dirichlet,argument.dof_maps)
    @test GT.substitute_element_csc_contribution_bases(raw,opaque) == raw
    dynamic = :(GT.contribute_element!(alloc,Ae,field_i,1,e,side_i,1))
    @test GT.substitute_element_csc_contribution_bases(dynamic,argument) == dynamic
    invalid_side = :(GT.contribute_element!(alloc,Ae,1,1,e,3,1))
    @test GT.substitute_element_csc_contribution_bases(invalid_side,argument) == invalid_side
    invalid_field = :(GT.contribute_element!(alloc,Ae,2,1,e,1,1))
    @test GT.substitute_element_csc_contribution_bases(invalid_field,argument) == invalid_field
    unfinished_cache = GT.ElementCSCDomainCache(zeros(Int32,0,GT.num_faces(skeleton)),
        fill(Int32(-1),1,1),zeros(Int32,4,1,1))
    unfinished = GT.MatrixAllocation(GT.ElementCSCDomainAllocation(zeros(Float64,0),unfinished_cache),
        argument.free_or_dirichlet,argument.dof_maps)
    @test GT.substitute_element_csc_contribution_bases(raw,unfinished) == raw
    ambiguous = GT.allocate_matrix(Float64,Vdg,Vdg,skeleton,skeleton;assembly_method)
    @test GT.substitute_element_csc_contribution_bases(raw,ambiguous) == raw
    Vdg2 = GT.lagrange_space(omega,2;continuous=Val(false))
    mixed = GT.cartesian_product(Vdg,Vdg2)
    mixed_alloc = GT.allocate_matrix(Float64,mixed,mixed,skeleton;assembly_method,
        assembly_options=(;matrix=(;block_mask=Bool[true false;true true])))
    mixed_argument = GT.element_csc_domain(mixed_alloc,1)
    disabled = :(GT.contribute_element!(alloc,Ae,1,2,e,1,1))
    @test GT.substitute_element_csc_contribution_bases(disabled,mixed_argument) === nothing
    disabled_dynamic_sides = :(GT.contribute_element!(alloc,Ae,1,2,e,side_i,side_j))
    @test GT.substitute_element_csc_contribution_bases(disabled_dynamic_sides,mixed_argument) === nothing
    enabled_dynamic_sides = :(GT.contribute_element!(alloc,Ae,2,1,e,side_i,side_j))
    @test GT.substitute_element_csc_contribution_bases(enabled_dynamic_sides,mixed_argument) == enabled_dynamic_sides
    rectangular = :(GT.contribute_element!(alloc,Ae,2,1,e,2,1))
    rectangular_call = GT.substitute_element_csc_contribution_bases(rectangular,mixed_argument)
    @test rectangular_call == :(GT.contribute_element_block!(alloc,Ae,e,100,9,4))
    diagonal = :(GT.contribute_element!(alloc,Ae,2,2,e,2,2))
    @test GT.substitute_element_csc_contribution_bases(diagonal,mixed_argument) ==
        :(GT.contribute_element_block!(alloc,Ae,e,451,9,9))
    storage = GT.element_csc_matrix_storage(argument.allocation)
    ne = GT.num_faces(skeleton)
    for e in 1:ne
        GT.contribute_element!(argument,blocks[1],1,1,e,1,1)
        GT.contribute_element!(argument,blocks[2],1,1,e,2,1)
        GT.contribute_element!(argument,blocks[3],1,1,e,1,2)
        GT.contribute_element!(argument,blocks[4],1,1,e,2,2)
    end
    expected = copy(storage.nzval)
    GT.reset!(alloc)
    factory = GT.evaluate(:((blocks,ne)->(()->(alloc->$prepared))),(blocks,ne))
    runner = GT.AssemblyRunner(factory,argument)
    GT.run_assembly!(runner,())
    @test storage.nzval == expected
    GT.reset!(alloc)
    @test (@allocated GT.run_assembly!(runner,())) == 0
    @test storage.nzval == expected
    # Compatible buffers can change while the literal packing remains the same.
    other = GT.allocate_matrix(Float64,Vdg,Vdg,skeleton;assembly_method)
    other_argument = GT.element_csc_domain(other,1)
    GT.run_assembly!(runner,other_argument,())
    @test GT.element_csc_matrix_storage(other_argument.allocation).nzval == expected

    # Unknown prototypes retain the computed-base path and the same result.
    fallback_factory = GT.evaluate(:((blocks,ne)->(()->(alloc->$raw))),(blocks,ne))
    fallback = GT.AssemblyRunner(fallback_factory,argument)
    GT.reset!(alloc)
    GT.run_assembly!(fallback,())
    @test storage.nzval == expected
    GT.reset!(alloc)
    @test (@allocated GT.run_assembly!(fallback,())) == 0
    @test storage.nzval == expected
    =#
end

struct LocalMatrixCollector{T}
    matrices::Vector{Matrix{T}}
end
Base.eltype(::LocalMatrixCollector{T}) where T = T
function GT.contribute!(collector::LocalMatrixCollector,values,rows,columns,field_i,field_j)
    push!(collector.matrices,copy(values))
    nothing
end

@testset "Poisson contributions from generated kernels" begin
    mesh = GT.cartesian_mesh((0,1,0,1,0,1),(2,2,2);simplexify=false)
    omega = GT.interior(mesh)
    for k in 1:3
        V = GT.lagrange_space(omega,k)
        measure = GT.measure(omega,2k)
        a = (u,v) -> GT.∫(x->dot(ForwardDiff.gradient(u,x),ForwardDiff.gradient(v,x)),measure)
        integral = a(GT.form_argument_quantity(V,2),GT.form_argument_quantity(V,1))
        contribution = only(GT.contributions(integral))
        collector = LocalMatrixCollector(Matrix{Float64}[])
        prepare = GT.generate_assemble_matrix_impl(contribution,V,V;allocation=collector)
        GT.run_assembly!(GT.AssemblyRunner(prepare,collector),())
        @test length(collector.matrices) == GT.num_faces(omega)
        alloc = GT.allocate_matrix(Float64,V,V,omega;assembly_method)
        face_dofs = GT.dofs_accessor(V,omega)
        expected = zeros(GT.num_free_dofs(V),GT.num_free_dofs(V))
        for face in reverse(1:GT.num_faces(omega))
            ids = face_dofs(face)
            values = collector.matrices[face]
            GT.contribute_element!(alloc,values,1,1,face)
            add_reference!(expected,values,ids,ids)
        end
        A = GT.compress(alloc)
        @test Matrix(A) ≈ expected rtol=1e-12 atol=1e-13
        @test norm(A-A',Inf) < 1e-11
        @test norm(A*ones(size(A,2)),Inf) < 1e-11
        if k == 1
            generated,generated_cache = GT.assemble_matrix(a,Float64,V,V;assembly_method,reuse=Val(true))
            @test generated ≈ A rtol=1e-12 atol=1e-13
            fill!(generated.nzval,0)
            @test GT.update_matrix!(generated,generated_cache) === generated
            @test generated ≈ A rtol=1e-12 atol=1e-13
        end
        GT.reset!(alloc)
        for face in 1:GT.num_faces(omega)
            GT.contribute_element!(alloc,collector.matrices[face],1,1,face)
        end
        @test Matrix(A) ≈ expected rtol=1e-12 atol=1e-13
    end
end

end

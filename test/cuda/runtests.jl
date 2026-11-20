using Test
using Adapt
using CUDA
using LinearAlgebra
using SparseArrays
import GalerkinToolkit as GT

# Run on a machine with a functional CUDA GPU. Missing hardware is a failure.
CUDA.functional(true) || error("CUDA tests require a functional GPU")
CUDA.allowscalar(false)

function scatter_blocks!(argument,values,nelements,base,nrows,ncols)
    element = (blockIdx().x-1)*blockDim().x+threadIdx().x
    if element <= nelements
        GT.contribute_element_block_atomic!(argument,values,element,base,nrows,ncols)
    end
    nothing
end

function scatter_buffers!(nzval,positions,values,nelements,base,nrows,ncols)
    element = (blockIdx().x-1)*blockDim().x+threadIdx().x
    if element <= nelements && base != -1
        GT.contribute_element_csc_block_atomic!(nzval,positions,values,element,base,(nrows,ncols))
    end
    nothing
end

#= TODO: Re-enable DG sides and multi-field GPU contributions.
function scatter_sides!(argument,values,nelements)
    element = (blockIdx().x-1)*blockDim().x+threadIdx().x
    if element <= nelements
        for side_j in 1:2, side_i in 1:2
            GT.contribute_element_atomic!(argument,values,1,1,element,side_i,side_j)
        end
        # A disabled field pair must not read its block dimensions or positions.
        GT.contribute_element_atomic!(argument,values,1,2,element)
    end
    nothing
end

=#
function launch_scatter!(kernel::F,argument,values,nelements,args...) where F
    nelements == 0 && return nothing
    CUDA.@sync @cuda threads=128 blocks=cld(nelements,128) kernel(argument,values,nelements,args...)
    nothing
end

@testset "CUDA element CSC assembly" begin
    @test Base.get_extension(GT,:GalerkinToolkitCUDAExt) !== nothing
    @testset "Adapt, layouts, and repeated atomic scatter" begin
        for T in (Float32,Float64), Ti in (Int32,Int64), layout in (GT.ElementMajorLayout(),GT.EntryMajorLayout()), nelements in (0,1,257)
            # Rectangular sub-block with a nonzero base and padded row capacity.
            # Many elements share entries; zero encodes a constrained entry.
            positions = reshape(Ti[mod(i,14) for i in 1:30nelements],30,nelements)
            cache = GT.ElementCSCDomainCache(positions,zeros(Ti,1,1),reshape(Ti[4,5,1,1],4,1,1))
            storage = GT.ElementCSCMatrixAllocation(13,1,Ti[1,14],Ti.(1:13),zeros(T,13),(cache,cache))
            host = GT.MatrixAllocation(storage,(GT.FREE,GT.FREE),(identity,identity))
            device = adapt(CuArray,GT.to_layout(host,layout))
            device_storage = device.allocation
            argument = GT.element_csc_domain(device,2)
            @test argument.allocation.nzval === device_storage.nzval
            @test isbitstype(typeof(CUDA.cudaconvert(argument)))
            @test eltype(device_storage.nzval) === T
            @test eltype(device_storage.colptr) === Ti
            lookup = argument.allocation.cache.positions
            if layout isa GT.EntryMajorLayout
                @test lookup isa Transpose
                @test size(parent(lookup)) == (nelements,30)
                @test Array(parent(lookup)) == permutedims(positions)
            else
                @test lookup isa CuMatrix{Ti}
            end
            values = reshape(T.(1:6),2,3)
            device_values = CuArray(values)
            expected = zeros(T,13)
            for element in 1:nelements, j in 1:3, i in 1:2
                p = positions[3+i+4(j-1),element]
                p == 0 || (expected[p] += values[i,j])
            end
            A,matrix_cache = GT.compress(device;reuse=Val(true))
            @test nonzeros(A) === device_storage.nzval
            @test size(A) == (13,1)
            for _ in 1:3
                GT.reset!(device)
                launch_scatter!(scatter_blocks!,argument,device_values,nelements,3,4,5)
                @test GT.compress!(device,A,matrix_cache) === A
                @test Array(nonzeros(A)) == expected
            end
            GT.reset!(device)
            if nelements != 0
                CUDA.@sync @cuda threads=128 blocks=cld(nelements,128) scatter_buffers!(
                    nonzeros(A),lookup,device_values,nelements,3,4,5)
            end
            @test Array(nonzeros(A)) == expected
            # A separate output receives a device-to-device copy during reuse.
            other = copy(A)
            fill!(nonzeros(other),T(-1))
            GT.compress!(device,other,matrix_cache)
            @test Array(nonzeros(other)) == expected
            launch_scatter!(scatter_blocks!,argument,device_values,nelements,-1,4,5)
            @test Array(nonzeros(A)) == expected
            restored = GT.to_layout(adapt(Array,device),GT.ElementMajorLayout())
            @test restored.allocation.domain_caches[2].positions == positions
            @test Matrix(GT.compress(restored)) == reshape(expected,13,1)
        end
    end

#= TODO: Re-enable DG sides and multi-field GPU contributions.
    @testset "Field masks, side indices, and shared domain values" begin
        for T in (Float32,Float64), Ti in (Int32,Int64)
            nelements = 257
            positions = reshape(Ti[mod(i,14) for i in 1:24nelements],24,nelements)
            sizes = zeros(Ti,4,2,2)
            sizes[:,1,1] = Ti[2,3,2,2]
            offsets = fill(-one(Ti),2,2)
            offsets[1,1] = 0
            cache = GT.ElementCSCDomainCache(positions,offsets,sizes)
            storage = GT.ElementCSCMatrixAllocation(13,1,Ti[1,14],Ti.(1:13),zeros(T,13),(cache,cache))
            wrapped = GT.MonolithicAssemblyAllocation(storage,(Ti[0,13],Ti[0,1]),Bool[true false; false false])
            host = GT.MatrixAllocation(wrapped,(GT.FREE,GT.FREE),(identity,identity))
            device = adapt(CuArray,GT.to_layout(host,GT.EntryMajorLayout()))
            first_domain = GT.element_csc_domain(device,1)
            second_domain = GT.element_csc_domain(device,2)
            @test first_domain.allocation.allocation.nzval === second_domain.allocation.allocation.nzval
            @test isbitstype(typeof(CUDA.cudaconvert(first_domain)))
            values = reshape(T.(1:6),2,3)
            device_values = CuArray(values)
            expected = zeros(T,13)
            for element in 1:nelements, sj in 1:2, si in 1:2, j in 1:3, i in 1:2
                p = positions[LinearIndices((2,3,2,2))[i,j,si,sj],element]
                p == 0 || (expected[p] += T(2)*values[i,j])
            end
            A = GT.compress(device)
            for _ in 1:3
                GT.reset!(device)
                launch_scatter!(scatter_sides!,first_domain,device_values,nelements)
                launch_scatter!(scatter_sides!,second_domain,device_values,nelements)
                @test Array(nonzeros(A)) == expected
            end
        end
    end

=#
    @testset "Empty sparse output" begin
        for T in (Float32,Float64), Ti in (Int32,Int64)
            cache = GT.ElementCSCDomainCache(zeros(Ti,0,0),fill(-one(Ti),1,1),zeros(Ti,4,1,1))
            storage = GT.ElementCSCMatrixAllocation(0,0,Ti[1],Ti[],T[],(cache,))
            device = adapt(CuArray,GT.to_layout(storage,GT.EntryMajorLayout()))
            A,cache = GT.compress(device;reuse=Val(true))
            @test size(A) == (0,0)
            @test nnz(A) == 0
            @test nonzeros(A) === device.nzval
            @test GT.reset!(device) === device
            @test GT.compress!(device,A,cache) === A
        end
    end

    @testset "Prepared quadratic FE allocation and reordered subset" begin
        mesh = GT.cartesian_mesh((0,1,0,1),(3,2);simplexify=false)
        omega = GT.interior(mesh)
        GT.group_faces(mesh,2)["cuda subset"] = Int32[6,2]
        subset = GT.interior(mesh;group_names=["cuda subset"])
        V = GT.lagrange_space(omega,2)
        method = (;matrix=GT.monolithic_assembly(GT.element_csc_assembly()))
        host = GT.allocate_matrix(Float64,V,V,omega,subset;assembly_method=method)
        device = adapt(CuArray,GT.to_layout(host,GT.EntryMajorLayout()))
        ndofs = GT.num_dofs(only(GT.reference_spaces(V)))
        values = reshape(Float64.(1:ndofs^2),ndofs,ndofs)
        for (domain_id,domain) in enumerate((omega,subset)), element in 1:GT.num_faces(domain)
            GT.contribute_element!(GT.element_csc_domain(host,domain_id),values,1,1,element)
        end
        expected = GT.compress(host)
        A = GT.compress(device)
        device_values = CuArray(values)
        for _ in 1:3
            GT.reset!(device)
            for (domain_id,domain) in enumerate((omega,subset))
                argument = GT.element_csc_domain(device,domain_id)
                launch_scatter!(scatter_blocks!,argument,device_values,GT.num_faces(domain),0,ndofs,ndofs)
            end
            @test SparseMatrixCSC(A) == expected
        end
    end
end

include("compiler_data_tests.jl")
include("poisson_numeric_tests.jl")
include("execution_tests.jl")

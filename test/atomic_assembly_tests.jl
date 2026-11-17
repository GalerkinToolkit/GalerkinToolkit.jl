module AtomicAssemblyTests

using Test
import GalerkinToolkit as GT

function matrix_argument(nzval,positions,shapes,offsets=zeros(eltype(positions),1,1))
    cache = GT.ElementCSCDomainCache(positions,offsets,shapes)
    GT.MatrixAllocation(GT.ElementCSCDomainAllocation(nzval,cache),
        (GT.FREE,GT.FREE),(identity,identity))
end

function no_warm_allocations(f::F) where F
    for _ in 1:3
        f()
    end
    (@allocated(f()),@allocations(f())) == (0,0)
end

@testset "Atomic element CSC contributions" begin
    for T in (Float32,Float64),Ti in (Int32,Int64)
        positions = reshape(Ti.(1:90),30,3)
        positions[5,2] = 0
        nzval = zeros(T,90)
        argument = matrix_argument(nzval,positions,reshape(Ti[4,5,1,1],4,1,1))
        values = reshape(T.(1:6),2,3)
        expected = zeros(T,90)
        expected[[34,38,39,42,43]] = T[1,3,4,5,6]
        @test GT.contribute_element_block_atomic!(argument,values,2,3,4,5) === argument
        @test nzval == expected
        @test GT.contribute_element_block_atomic!(argument,values,2,-1,4,5) === argument
        @test nzval == expected
        for empty_values in (zeros(T,0,3),zeros(T,2,0))
            GT.contribute_element_block_atomic!(argument,empty_values,2,3,4,5)
            @test nzval == expected
        end
        block_call = let argument=argument,values=values
            ()->GT.contribute_element_block_atomic!(argument,values,2,3,4,5)
        end
        @test no_warm_allocations(block_call)

        # Two test/trial sides with a rectangular block and shared global entries.
        positions = reshape(Ti[mod(i,13)+1 for i in 0:95],24,4)
        positions[3,2] = 0
        shapes = reshape(Ti[2,3,2,2],4,1,1)
        nzval = zeros(T,13)
        atomic = matrix_argument(nzval,positions,shapes)
        sequential = matrix_argument(zeros(T,13),positions,shapes)
        for face in reverse(axes(positions,2)),side_j in 1:2,side_i in 1:2
            values = fill(T(10side_i+side_j),2,3)
            GT.contribute_element!(sequential,values,1,1,face,side_i,side_j)
            @test GT.contribute_element_atomic!(atomic,values,1,1,face,side_i,side_j) === atomic
        end
        @test nzval == sequential.allocation.nzval
        values = ones(T,2,3)
        element_call = let atomic=atomic,values=values
            ()->GT.contribute_element_atomic!(atomic,values,1,1,1,2,2)
        end
        @test no_warm_allocations(element_call)

        # Both domain selectors share numerical storage; no mesh is retained.
        cache = atomic.allocation.cache
        parent_storage = GT.ElementCSCMatrixAllocation(13,1,Ti[1,14],Ti.(1:13),nzval,(cache,cache))
        parent = GT.MatrixAllocation(parent_storage,(GT.FREE,GT.FREE),(identity,identity))
        fill!(nzval,0)
        for domain_id in 1:2
            selected = GT.element_csc_domain(parent,domain_id)
            GT.contribute_element_block_atomic!(selected,values,1,0,2,3)
        end
        expected = zeros(T,13)
        for p in positions[1:6,1]
            p != 0 && (expected[p] += T(2))
        end
        @test nzval == expected
        @test parent_storage.nzval === atomic.allocation.nzval
    end
end

@testset "Atomic field masks and direct vector contributions" begin
    for T in (Float32,Float64)
        positions = reshape(Int32[1,2,3,4],2,2)
        shapes = zeros(Int32,4,2,2)
        shapes[:,1,1] = Int32[1,1,1,1]
        shapes[:,2,2] = Int32[1,1,1,1]
        offsets = Int32[0 -1; -1 1]
        argument = matrix_argument(zeros(T,4),positions,shapes,offsets)
        mask = Bool[true false; false true]
        wrapped = GT.MatrixAllocation(GT.MonolithicAssemblyAllocation(argument.allocation,
            (Int32[0,2],Int32[0,2]),mask),argument.free_or_dirichlet,argument.dof_maps)
        GT.contribute_element_atomic!(wrapped,ones(T,1,1),1,1,1)
        GT.contribute_element_atomic!(wrapped,ones(T,1,1),2,2,2)
        GT.contribute_element_atomic!(wrapped,ones(T,1,1),1,2,1)
        @test argument.allocation.nzval == T[1,0,0,1]

        dofs = Int32[-1,1,-2,2]
        values = T[1,2,3,4]
        for selection in (GT.FREE,GT.DIRICHLET)
            direct = GT.DirectVectorAllocation(zeros(T,4))
            block = GT.MonolithicAssemblyAllocation(direct,(Int32[0,2],),Bool[false,true])
            map_i = selection == GT.FREE ? identity : -
            vector = GT.VectorAllocation(block,selection,map_i)
            @test GT.contribute_atomic!(vector,values,dofs,2) === nothing
            @test direct.vector == (selection == GT.FREE ? T[0,0,2,4] : T[0,0,1,3])
            expected = copy(direct.vector)
            GT.contribute_atomic!(vector,values,dofs,1)
            @test direct.vector == expected
            @test no_warm_allocations(()->GT.contribute_atomic!(vector,values,dofs,2))
        end
    end
end

function concurrent_matrix!(argument,values,repetitions)
    @sync for _ in 1:Threads.nthreads()
        Threads.@spawn begin
            for _ in 1:repetitions,face in 1:4
                GT.contribute_element_block_atomic!(argument,values,face,0,2,2)
            end
        end
    end
end

function concurrent_vector!(argument,values,dofs,repetitions)
    @sync for _ in 1:Threads.nthreads()
        Threads.@spawn begin
            for _ in 1:repetitions
                GT.contribute_atomic!(argument,values,dofs)
            end
        end
    end
end

@testset "Concurrent atomic matrix and vector updates" begin
    for T in (Float32,Float64)
        positions = repeat(reshape(Int32[1,2,0,3],4,1),1,4)
        nzval = zeros(T,3)
        argument = matrix_argument(nzval,positions,reshape(Int32[2,2,1,1],4,1,1))
        values = ones(T,2,2)
        repetitions = 1_000
        for _ in 1:3
            fill!(nzval,0)
            concurrent_matrix!(argument,values,repetitions)
            @test nzval == fill(T(4repetitions*Threads.nthreads()),3)
        end
        storage = GT.DirectVectorAllocation(zeros(T,2))
        vector = GT.VectorAllocation(storage,GT.FREE,identity)
        values,dofs = ones(T,4),Int32[1,2,1,2]
        for _ in 1:3
            GT.reset!(vector)
            concurrent_vector!(vector,values,dofs,repetitions)
            @test storage.vector == fill(T(2repetitions*Threads.nthreads()),2)
        end
    end
end

end # module

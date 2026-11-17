module ContributionTests

using Test
using SparseArrays
import GalerkinToolkit as GT

@testset "COO contributions return nothing" begin
    counter = GT.counter(GT.coo_assembly(),Float64,1:2)
    @test GT.contribute!(counter,nothing,1,1) === nothing
    @test GT.contribute!(counter,nothing,2,1) === nothing
    @test counter.nnz[] == 2
    storage = GT.allocate(counter)
    @test GT.contribute!(storage,4.0,2,1) === nothing
    @test GT.contribute!(storage,-1.0,2,1) === nothing
    @test GT.compress(storage) == [0.0,3.0]

    GT.reset!(storage)
    vector = GT.VectorAllocation(storage,GT.FREE,identity)
    @test GT.contribute!(vector,[1.0,2.0,3.0],[-1,1,2]) === nothing
    @test GT.contribute!(vector,Float64[],Int[]) === nothing
    @test GT.compress(vector) == [2.0,3.0]

    counter = GT.counter(GT.coo_assembly(),Float64,1:2,1:2)
    for j in 1:2,i in 1:2
        @test GT.contribute!(counter,nothing,i,j,1,1) === nothing
    end
    @test counter.nnz[] == 4
    storage = GT.allocate(counter)
    for j in 1:2,i in 1:2
        @test GT.contribute!(storage,Float64(i==j),i,j,1,1) === nothing
    end
    @test Matrix(GT.compress(storage;reuse=Val(false))) == [1.0 0.0; 0.0 1.0]
    GT.reset!(storage)
    matrix = GT.MatrixAllocation(storage,(GT.FREE,GT.FREE),(identity,identity))
    @test GT.contribute!(matrix,[2.0 0.0; 0.0 3.0],1:2,1:2) === nothing
    @test GT.contribute!(matrix,zeros(0,0),Int[],Int[]) === nothing
    @test Matrix(GT.compress(matrix)) == [2.0 0.0; 0.0 3.0]
end

@testset "CSC contributions and cached sparsity stages return nothing" begin
    entries = ((2,2),(1,1),(2,2))
    values = (4.0,2.0,-1.0)
    for cached in (false,true)
        counter = GT.counter(GT.csc_assembly(Val(cached)),Float64,1:2,1:2)
        for (i,j) in entries
            @test GT.contribute!(counter,nothing,i,j,1,1) === nothing
        end
        @test counter.colnnz == [1,2]
        storage = GT.allocate(counter)
        if cached
            @test storage isa GT.CachedCSCCounterStep1
            for (i,j) in entries
                @test GT.contribute!(storage,nothing,i,j,1,1) === nothing
            end
            storage = GT.allocate(storage)
            @test storage isa GT.CachedCSCCounterStep2
            for (i,j) in entries
                @test GT.contribute!(storage,nothing,i,j,1,1) === nothing
            end
            storage = GT.finalize(storage)
        end
        for ((i,j),value) in zip(entries,values)
            @test GT.contribute!(storage,value,i,j,1,1) === nothing
        end
        @test Matrix(GT.compress(storage;reuse=Val(false))) == [2.0 0.0; 0.0 3.0]
    end
end

@testset "Monolithic contribution masks return nothing" begin
    counter = GT.counter(GT.coo_assembly(),Float64,1:2)
    wrapped = GT.MonolithicAssemblyCounter(counter,(Int32[0,1],),Bool[true,false])
    @test GT.contribute!(wrapped,nothing,1,1) === nothing
    @test GT.contribute!(wrapped,nothing,1,2) === nothing
    @test counter.nnz[] == 1

    counter = GT.counter(GT.coo_assembly(),Float64,1:2,1:2)
    mask = Bool[true false; false true]
    wrapped = GT.MonolithicAssemblyCounter(counter,(Int32[0,1],Int32[0,1]),mask)
    @test GT.contribute!(wrapped,nothing,1,1,1,1) === nothing
    @test GT.contribute!(wrapped,nothing,1,1,1,2) === nothing
    @test GT.contribute!(wrapped,nothing,1,1,2,2) === nothing
    @test counter.nnz[] == 2
    storage = GT.allocate(wrapped)
    @test GT.contribute!(storage,2.0,1,1,1,1) === nothing
    @test GT.contribute!(storage,100.0,1,1,1,2) === nothing
    @test GT.contribute!(storage,3.0,1,1,2,2) === nothing
    @test Matrix(GT.compress(storage;reuse=Val(false))) == [2.0 0.0; 0.0 3.0]

    direct = GT.DirectVectorAllocation(zeros(2))
    @test GT.contribute!(direct,4.0,1,1) === nothing
    wrapped = GT.MonolithicAssemblyAllocation(direct,(Int32[0,1],),Bool[true,false])
    @test GT.contribute!(wrapped,1.0,1,1) === nothing
    @test GT.contribute!(wrapped,100.0,1,2) === nothing
    @test GT.contribute_atomic!(direct,2.0,1,1) === nothing
    @test GT.contribute_atomic!(wrapped,3.0,1,1) === nothing
    @test GT.contribute_atomic!(wrapped,100.0,1,2) === nothing
    @test GT.compress(direct) == [10.0,0.0]
end

end # module

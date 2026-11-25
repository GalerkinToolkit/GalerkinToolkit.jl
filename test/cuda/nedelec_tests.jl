module CUDANedelecTests

using Test, CUDA, LinearAlgebra, SparseArrays, StaticArrays
import ForwardDiff
import GalerkinToolkit as GT

const CuSparseCSC = Base.get_extension(GT,:GalerkinToolkitCUDAExt).CuSparseMatrixCSC
const assembly_method = (;matrix=GT.element_csc_assembly())
CUDA.allowscalar(false)

curl(J::SMatrix{3,3}) = SVector(J[3,2]-J[2,3],J[1,3]-J[3,1],J[2,1]-J[1,2])
curl(J::SMatrix{2,2}) = J[2,1]-J[1,2]

@testset "Generated CUDA Nédélec matrices" begin
    for dim in (2,3), dirichlet in (false,true)
        cells = dim == 3 ? (3,3,2) : (4,3)
        mesh = GT.cartesian_mesh(ntuple(i->isodd(i) ? 0 : 1,2dim),cells;simplexify=true)
        mesh = GT.replace_node_coordinates(mesh,
            [SVector{dim}(ntuple(i->1.2x[i]+0.1(sum(x)-x[i]),dim)) for x in GT.node_coordinates(mesh)])
        omega = GT.interior(mesh)
        V = dirichlet ? GT.nedelec_space(omega,0;dirichlet_boundary=GT.boundary(mesh)) :
            GT.nedelec_space(omega,0)
        measure = GT.measure(omega,2)
        scale = 1.25
        a = (u,v)->GT.∫(x->dot(GT.call(curl,ForwardDiff.jacobian(u,x)),GT.call(curl,ForwardDiff.jacobian(v,x)))+
            scale*dot(u(x),v(x)),measure)
        expected = GT.assemble_matrix(a,Float64,V,V;parameters=(scale,),reuse=Val(false),assembly_method)
        execution = GT.CUDAExecution()
        A,cache = GT.assemble_matrix(a,Float64,V,V;execution,assembly_method,parameters=(scale,),reuse=Val(true))
        @test A isa CuSparseCSC
        @test SparseMatrixCSC(A) ≈ expected rtol=2e-11 atol=1e-12
        @test isbitstype(typeof(CUDA.cudaconvert(only(cache.loops).data)))
        for factor in (3.0,0.0)
            @test GT.update_matrix!(A,cache;parameters=(factor,)) === A
            shifted = GT.assemble_matrix(a,Float64,V,V;parameters=(factor,),reuse=Val(false),assembly_method)
            @test SparseMatrixCSC(A) ≈ shifted rtol=2e-11 atol=1e-12
        end
    end
end

@testset "Generated CUDA Nédélec mass matrix in single and double precision" begin
    mesh = GT.cartesian_mesh((0,1,0,1,0,1),(2,2,2);simplexify=true)
    omega = GT.interior(mesh)
    V = GT.nedelec_space(omega,0)
    measure = GT.measure(omega,2)
    a = (u,v)->GT.∫(x->dot(u(x),v(x)),measure)
    for T in (Float32,Float64)
        expected = GT.assemble_matrix(a,T,V,V)
        A,cache = GT.assemble_matrix(a,T,V,V;execution=GT.CUDAExecution(),assembly_method,reuse=Val(true))
        @test eltype(nonzeros(A)) == T
        @test SparseMatrixCSC(A) ≈ expected
        @test GT.update_matrix!(A,cache) === A
        @test SparseMatrixCSC(A) ≈ expected
    end
end

end # module

module CUDAAssemblyExecutionTests

using Test, CUDA, LinearAlgebra, SparseArrays, StaticArrays
import ForwardDiff
import GalerkinToolkit as GT

const CuSparseCSC = Base.get_extension(GT,:GalerkinToolkitCUDAExt).CuSparseMatrixCSC
const assembly_method = (;matrix=GT.element_csc_assembly())
CUDA.allowscalar(false)

@testset "Generated CUDA CG matrix assembly" begin
    @test GT.CUDAExecution().threads == 128
    @test_throws ArgumentError GT.CUDAExecution(threads=0)
    for (dim,simplexify) in ((1,false),(2,false),(2,true),(3,false),(3,true)), k in 1:3
        mesh = GT.cartesian_mesh(ntuple(i->isodd(i) ? 0 : 1,2dim),ntuple(_->2,dim);simplexify)
        mesh = GT.replace_node_coordinates(mesh,
            [SVector{dim}(ntuple(i->1.2x[i]+0.1(sum(x)-x[i]),dim)) for x in GT.node_coordinates(mesh)])
        omega = GT.interior(mesh)
        V = GT.lagrange_space(omega,k)
        measure = GT.measure(omega,2k)
        a = (u,v)->GT.∫(x->dot(ForwardDiff.gradient(u,x),ForwardDiff.gradient(v,x)),measure)
        expected = GT.assemble_matrix(a,Float64,V,V;assembly_method)
        A,cache = GT.assemble_matrix(a,Float64,V,V;
            execution=GT.CUDAExecution(),assembly_method,reuse=Val(true))
        @test A isa CuSparseCSC
        @test SparseMatrixCSC(A) ≈ expected rtol=2e-11 atol=1e-12
        values = nonzeros(A)
        runner = only(cache.loops)
        data = runner.data
        @test isbitstype(typeof(CUDA.cudaconvert(data)))
        @test GT.element_csc_matrix_storage(cache.alloc.allocation).nzval === values
        for _ in 1:2
            @test GT.update_matrix!(A,cache;parameters=()) === A
            @test SparseMatrixCSC(A) ≈ expected rtol=2e-11 atol=1e-12
            @test runner.data === data
            @test nonzeros(A) === values
        end
    end
end

@testset "Generated CUDA parameters, overlapping subsets and constrained DOFs" begin
    for k in (1,2),dirichlet in (false,true)
        mesh = GT.cartesian_mesh((0,1,0,1,0,1),(3,3,3);simplexify=false)
        omega = GT.interior(mesh)
        GT.group_faces(mesh,3)["CUDA generated subset"] = Int32[27,2]
        subset = GT.interior(mesh;group_names=["CUDA generated subset"])
        V = dirichlet ? GT.lagrange_space(omega,k;dirichlet_boundary=GT.boundary(mesh)) :
            GT.lagrange_space(omega,k)
        full_measure,subset_measure = GT.measure(omega,2k),GT.measure(subset,2k)
        scale = 1.25
        a = (u,v)->GT.∫(x->scale*dot(ForwardDiff.gradient(u,x),ForwardDiff.gradient(v,x)),full_measure)+
                  GT.∫(x->scale*u(x)*v(x),subset_measure)
        expected = GT.assemble_matrix(a,Float64,V,V;parameters=(scale,),
            reuse=Val(false),assembly_method)
        execution = GT.CUDAExecution(threads=64)
        A,cache = GT.assemble_matrix(a,Float64,V,V;execution,assembly_method,parameters=(scale,),reuse=Val(true))
        @test SparseMatrixCSC(A) ≈ expected rtol=2e-11 atol=1e-12
        @test length(cache.loops) == 2
        @test cache.execution === execution
        output = nonzeros(A)
        inputs = map(r->r.data,cache.loops)
        for factor in (3.0,0.0,1.0)
            @test GT.update_matrix!(A,cache;parameters=(factor,)) === A
            @test SparseMatrixCSC(A) ≈ (factor/scale)*expected rtol=2e-11 atol=1e-12
            @test all(r.data === data for (r,data) in zip(cache.loops,inputs))
            @test nonzeros(A) === output
        end
        before = Array(nonzeros(A))
        @test_throws ArgumentError GT.update_matrix!(A,cache;parameters=([2.0],))
        @test_throws ArgumentError GT.update_matrix!(A,cache;parameters=())
        @test_throws ArgumentError GT.update_matrix!(A,cache;parameters=(2.0,3.0))
        @test Array(nonzeros(A)) == before
        other = copy(A)
        @test GT.update_matrix!(other,cache;parameters=(2.0,)) === other
        @test SparseMatrixCSC(other) ≈ (2/scale)*expected rtol=2e-11 atol=1e-12
        @test nonzeros(other) !== nonzeros(A)
    end
end

@testset "Generated CUDA mass matrix, storage errors and empty domain" begin
    mesh = GT.cartesian_mesh((0,1,0,1,0,1),(2,2,2);simplexify=false)
    omega = GT.interior(mesh)
    V = GT.lagrange_space(omega,1)
    measure = GT.measure(omega,2)
    a = (u,v)->GT.∫(x->u(x)*v(x),measure)
    execution = GT.CUDAExecution()
    for T in (Float32,Float64)
        expected = GT.assemble_matrix(a,T,V,V)
        A,cache = GT.assemble_matrix(a,T,V,V;execution,assembly_method,reuse=Val(true))
        @test SparseMatrixCSC(A) ≈ expected
        @test GT.update_matrix!(A,cache) === A
        @test SparseMatrixCSC(A) ≈ expected
    end
    @test_throws ArgumentError GT.assemble_matrix(a,Float64,V,V;
        execution=GT.CUDAExecution(threads=10_000),assembly_method)
    @test_throws ArgumentError GT.assemble_matrix(a,Float64,V,V;execution,
        assembly_method=(;matrix=GT.coo_assembly()))
    GT.group_faces(mesh,3)["CUDA generated empty"] = Int32[]
    empty_domain = GT.interior(mesh;group_names=["CUDA generated empty"])
    empty_measure = GT.measure(empty_domain,2)
    empty_form = (u,v)->GT.∫(x->u(x)*v(x),empty_measure)
    A,cache = GT.assemble_matrix(empty_form,Float64,V,V;execution,assembly_method,reuse=Val(true))
    @test nnz(A) == 0
    @test only(cache.loops).nelements == 0
    @test GT.update_matrix!(A,cache) === A
end

@testset "Generated CUDA coordinates and rectangular CG matrix" begin
    mesh = GT.cartesian_mesh((0,1,0,1,0,1),(2,2,2);simplexify=false)
    mesh = GT.replace_node_coordinates(mesh,
        [SVector(x[1]+0.15x[1]*x[2],x[2]+0.1x[1]*x[3],x[3]+0.05x[2]*x[3])
            for x in GT.node_coordinates(mesh)])
    omega = GT.interior(mesh)
    U,V = GT.lagrange_space(omega,2),GT.lagrange_space(omega,1)
    measure = GT.measure(omega,4)
    diffusion = x->1+x[1]^2
    a = (u,v)->GT.∫(x->GT.call(diffusion,x)*dot(ForwardDiff.gradient(u,x),ForwardDiff.gradient(v,x))+
        u(x)*v(x),measure)
    integral = a(GT.form_argument_quantity(U,2),GT.form_argument_quantity(V,1))
    Q = GT.quadrature(only(GT.contributions(integral)))
    coordinates = GT.coordinate_accessor(Q)
    weights = GT.weight_accessor(Q)
    npoints = GT.num_points_accessor(Q)
    trial_dofs, test_dofs = GT.dofs_accessor(U,omega), GT.dofs_accessor(V,omega)
    trial_gradients = GT.shape_function_accessor(ForwardDiff.gradient,U,Q)
    test_gradients = GT.shape_function_accessor(ForwardDiff.gradient,V,Q)
    trial_values = GT.shape_function_accessor(GT.value,U,Q)
    test_values = GT.shape_function_accessor(GT.value,V,Q)
    expected = zeros(GT.num_free_dofs(V),GT.num_free_dofs(U))
    for face in 1:GT.num_faces(omega)
        rows, cols = test_dofs(face,1), trial_dofs(face,1)
        for q in 1:npoints(face)
            scale = diffusion(coordinates(face)(q))
            w = weights(face)(q)
            grad_u, grad_v = trial_gradients(face,1)(q), test_gradients(face,1)(q)
            u, v = trial_values(face,1)(q), test_values(face,1)(q)
            for j in eachindex(cols), i in eachindex(rows)
                expected[rows[i],cols[j]] += w*(scale*dot(grad_u(j),grad_v(i)) + u(j)*v(i))
            end
        end
    end
    cpu = GT.assemble_matrix(a,Float64,U,V;assembly_method)
    @test Matrix(cpu) ≈ expected rtol=2e-11 atol=1e-12
    A,cache = GT.assemble_matrix(a,Float64,U,V;execution=GT.CUDAExecution(),
        assembly_method,reuse=Val(true))
    @test size(A) == size(expected)
    @test size(A,1) != size(A,2)
    @test Matrix(SparseMatrixCSC(A)) ≈ expected rtol=2e-11 atol=1e-12
    @test GT.update_matrix!(A,cache) === A
    @test Matrix(SparseMatrixCSC(A)) ≈ expected rtol=2e-11 atol=1e-12
end

end # module

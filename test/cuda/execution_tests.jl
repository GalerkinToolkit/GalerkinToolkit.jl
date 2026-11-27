module CUDAAssemblyExecutionTests

using Test, CUDA, LinearAlgebra, SparseArrays, StaticArrays
import ForwardDiff
import GalerkinToolkit as GT

const Ext = Base.get_extension(GT,:GalerkinToolkitCUDAExt)
const CuSparseCSC = Ext.CuSparseMatrixCSC
const assembly_method = (;matrix=GT.element_csc_assembly())
CUDA.allowscalar(false)

@testset "Generated CUDA CG matrix assembly" begin
    @test GT.CUDAExecution().threads == 128
    @test GT.CUDAExecution().opt_level == 3
    @test GT.CUDAExecution().fastmath
    @test_throws ArgumentError GT.CUDAExecution(threads=0)
    @test_throws ArgumentError GT.CUDAExecution(opt_level=-1)
    @test_throws ArgumentError GT.CUDAExecution(opt_level=4)
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
    for (T,opt_level) in ((Float32,2),(Float64,3))
        fastmath = T === Float64
        execution = GT.CUDAExecution(;opt_level,fastmath)
        options = Ext.CUDAAssemblyCompilerOptions(execution.opt_level,execution.fastmath)
        config = Ext.CUDACompilerBackend.compiler_config(CUDA.device();name=options)
        @test config.opt_level == opt_level
        @test config.always_inline && config.target.fastmath == fastmath
        expected = GT.assemble_matrix(a,T,V,V;execution=GT.SerialExecution(;fastmath))
        A,cache = GT.assemble_matrix(a,T,V,V;execution,assembly_method,reuse=Val(true))
        @test SparseMatrixCSC(A) ≈ expected
        @test only(cache.loops).execution === execution
        @test GT.update_matrix!(A,cache) === A
        @test SparseMatrixCSC(A) ≈ expected
    end
    execution = GT.CUDAExecution()
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
    for (trial,test,reference) in ((U,V,expected),(V,U,transpose(expected)))
        A,cache = GT.assemble_matrix(a,Float64,trial,test;execution=GT.CUDAExecution(),
            assembly_method,reuse=Val(true))
        @test size(A) == size(reference)
        @test size(A,1) != size(A,2)
        @test Matrix(SparseMatrixCSC(A)) ≈ reference rtol=2e-11 atol=1e-12
        @test GT.update_matrix!(A,cache) === A
        @test Matrix(SparseMatrixCSC(A)) ≈ reference rtol=2e-11 atol=1e-12
    end
end

@testset "Generated CUDA matrix assembly with discrete-field coefficients" begin
    ∇ = ForwardDiff.gradient
    distort(mesh,dim) = GT.replace_node_coordinates(mesh,
        [SVector{dim}(ntuple(i->x[i]+0.08*sin(2x[mod1(i+1,dim)])+0.05*x[i]*x[mod1(i+1,dim)],dim))
         for x in GT.node_coordinates(mesh)])
    g = x -> 1 + x[1]^2 - x[1]*x[end]/2 + 3sin(3x[end])/10
    form(uh,dμ) = (u,v)->GT.∫(x->uh(x)*u(x)*v(x) + dot(∇(uh,x),∇(u,x))*v(x) +
        uh(x)^2*dot(∇(u,x),∇(v,x)),dμ)
    # The field is captured by the form; its arrays become device inputs of the kernel.
    for (dim,simplexify,kf,k) in ((2,false,2,1),(3,true,1,1),(3,true,2,2),(3,false,2,1))
        mesh = distort(GT.cartesian_mesh(ntuple(i->isodd(i) ? 0 : 1,2dim),ntuple(_->2,dim);simplexify),dim)
        ncells = GT.num_faces(mesh,dim)
        GT.group_faces(mesh,dim)["CUDA field subset"] = Int32[1,ncells]
        Ω = GT.interior(mesh)
        Ωs = GT.interior(mesh;group_names=["CUDA field subset"])
        uh = GT.interpolate(GT.analytical_field(g,Ω),
            GT.lagrange_space(Ω,kf;dirichlet_boundary=GT.boundary(mesh)))
        @test GT.num_dirichlet_dofs(GT.space(uh)) > 0
        V = GT.lagrange_space(Ω,k)
        for dμ in (GT.measure(Ω,2kf+k),GT.measure(Ωs,2kf+k))
            a = form(uh,dμ)
            expected = GT.assemble_matrix(a,Float64,V,V;assembly_method)
            A,cache = GT.assemble_matrix(a,Float64,V,V;
                execution=GT.CUDAExecution(),assembly_method,reuse=Val(true))
            @test SparseMatrixCSC(A) ≈ expected rtol=2e-11 atol=1e-12
            runner = only(cache.loops)
            @test isbitstype(typeof(CUDA.cudaconvert(runner.data)))
            @test GT.update_matrix!(A,cache;parameters=()) === A
            @test SparseMatrixCSC(A) ≈ expected rtol=2e-11 atol=1e-12
        end
    end
    # Field values held on the device are read by reference: changing them in place
    # changes the next update_matrix!.
    mesh = distort(GT.cartesian_mesh((0,1,0,1,0,1),(2,2,2);simplexify=true),3)
    Ω = GT.interior(mesh)
    Pf = GT.lagrange_space(Ω,2;dirichlet_boundary=GT.boundary(mesh))
    V = GT.lagrange_space(Ω,1)
    dμ = GT.measure(Ω,5)
    host = GT.interpolate(GT.analytical_field(g,Ω),Pf)
    free,diri = CuArray(GT.free_values(host)),CuArray(GT.dirichlet_values(host))
    uh = GT.discrete_field(Pf,free,diri)
    @test GT.tabulated_lagrange_discrete_field(GT.value,uh,GT.leaf_term(nothing),Ω)
    A,cache = GT.assemble_matrix(form(uh,dμ),Float64,V,V;
        execution=GT.CUDAExecution(),assembly_method,reuse=Val(true))
    @test SparseMatrixCSC(A) ≈ GT.assemble_matrix(form(host,dμ),Float64,V,V;assembly_method) rtol=2e-11 atol=1e-12
    for scale in (2.0,-0.5)
        new_free,new_diri = scale*GT.free_values(host),(1+scale)*GT.dirichlet_values(host)
        copyto!(free,new_free)
        copyto!(diri,new_diri)
        expected = GT.assemble_matrix(form(GT.discrete_field(Pf,new_free,new_diri),dμ),Float64,V,V;
            assembly_method)
        @test GT.update_matrix!(A,cache;parameters=()) === A
        @test SparseMatrixCSC(A) ≈ expected rtol=2e-11 atol=1e-12
    end
end

end # module

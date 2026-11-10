module ProblemsRunnerTests

using Test
import GalerkinToolkit as GT

function warmed_update_allocations(update!::F,output,cache) where F
    for _ in 1:3
        update!(output,cache)
    end
    bytes = @allocated update!(output,cache)
    count = @allocations update!(output,cache)
    bytes,count
end

function runner_tests()
    # Evaluating these functions inside an already compiled caller exercises
    # the latest-world boundary, including parameter-dependent preparation.
    prepare = GT.evaluate(:((input,) -> ((scale,) -> begin
        factor = scale[]
        output -> begin
            for i in eachindex(input)
                output[i] = factor*input[i]
            end
            nothing
        end
    end)),([1.0,2.0],))
    output = zeros(2)
    runner = GT.AssemblyRunner(prepare,output)
    scale = Ref(2.0)
    GT.run_assembly!(runner,(scale,))
    @test output == [2,4]
    scale[] = 3
    GT.run_assembly!(runner,(scale,))
    @test output == [3,6]
    other = zeros(2)
    GT.run_assembly!(runner,other,(scale,))
    @test runner.argument === other
    @test other == [3,6]
    other32 = zeros(Float32,2)
    GT.run_assembly!(runner,other32,(scale,))
    @test other32 == [3,6]

    scalar_prepare = GT.evaluate(:(() -> ((scale,) -> begin
        factor = scale[]
        init -> init+factor
    end)),())
    scalar = GT.ScalarAssemblyRunner(scalar_prepare,1.0)
    @test GT.update_scalar((scalar,);parameters=(scale,)) == 4
    scale[] = 5
    @test GT.update_scalar((scalar,);parameters=(scale,)) == 6
    @test GT.update_scalar((scalar,);parameters=(scale,)) == 6

    fixed = GT.evaluate(:((input,) -> (() -> (output -> copyto!(output,input)))),([1.0,2.0],))
    fixed_runner = GT.AssemblyRunner(fixed,zeros(2))
    for _ in 1:3
        GT.run_assembly!(fixed_runner,())
    end
    @test (@allocated GT.run_assembly!(fixed_runner,())) == 0
    @test (@allocations GT.run_assembly!(fixed_runner,())) == 0

    fixed_scalar = GT.evaluate(:((input,) -> (() -> (init -> init+input[1]))),([2.0],))
    scalar_runner = GT.ScalarAssemblyRunner(fixed_scalar,1.0)
    scalar_loops = (scalar_runner,)
    for _ in 1:3
        GT.update_scalar(scalar_loops)
    end
    @test GT.update_scalar(scalar_loops) == 3
    @test (@allocated GT.update_scalar(scalar_loops)) == 0
    @test (@allocations GT.update_scalar(scalar_loops)) == 0
end

@testset "Generated assembly runners" begin
    runner_tests()
end

function public_runner_tests()
    mesh = GT.cartesian_mesh((0,1,0,1),(1,1);simplexify=false)
    omega = GT.interior(mesh)
    Q = GT.measure(omega,2)
    f = x -> GT.call(sum,x)

    samples,cache = GT.sample(f,Q;reuse=Val(true))
    expected_samples = collect.(samples)
    @test sum(sum,samples) ≈ 4
    foreach(values -> fill!(values,-1),samples)
    @test GT.update_sample!(samples,cache) === samples
    @test collect.(samples) ≈ expected_samples
    other = GT.sample(f,Q)
    foreach(values -> fill!(values,-2),other)
    GT.update_sample!(other,cache)
    @test collect.(other) ≈ expected_samples

    integral = GT.∫(x -> 1,Q)
    value,loops = GT.assemble_scalar(integral;reuse=Val(true))
    @test value ≈ 1
    @test GT.update_scalar(loops) ≈ value
    doubled,double_loops = GT.assemble_scalar(integral+integral;reuse=Val(true))
    @test doubled ≈ 2
    @test GT.update_scalar(double_loops) ≈ doubled
    @test sum(integral) ≈ 1

    contribution = only(GT.contributions(integral))
    faces = zeros(GT.num_faces(omega))
    result,face_cache = GT.assemble_face_contribution!(faces,contribution;reuse=Val(true))
    @test result === faces
    @test sum(faces) ≈ value
    fill!(faces,-1)
    @test GT.update_face_contribution!(faces,face_cache) === faces
    @test sum(faces) ≈ value
    other_faces = zeros(Float32,length(faces))
    GT.update_face_contribution!(other_faces,face_cache)
    @test sum(other_faces) ≈ value
    @test sum(GT.face_contribution(integral,omega)) ≈ value

    V = GT.lagrange_space(omega,1)
    a = (u,v) -> GT.∫(x -> u(x)*v(x),Q)
    l = v -> GT.∫(x -> v(x),Q)
    A,matrix_cache = GT.assemble_matrix(a,Float64,V,V;reuse=Val(true),assembly_method=(;))
    b,vector_cache = GT.assemble_vector(l,Float64,V;reuse=Val(true),assembly_method=(;))
    @test sum(A) ≈ value
    @test sum(b) ≈ value
    @test warmed_update_allocations(GT.update_matrix!,A,matrix_cache) == (0,0)
    @test warmed_update_allocations(GT.update_vector!,b,vector_cache) == (0,0)
    @test sum(A) ≈ value
    @test sum(b) ≈ value
end

@testset "Public generated problem execution" begin
    public_runner_tests()
end

end

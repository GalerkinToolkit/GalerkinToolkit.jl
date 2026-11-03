module CompilerDataTests

using Test
using StaticArrays
import ForwardDiff
import GalerkinToolkit as GT

@testset "Mesh integration layouts" begin
    for dim in (0, 1, 2, 3), nnodes in (0, 5)
        coordinates = SVector{dim,Float32}[SVector{dim,Float32}(10n + c for c in 1:dim) for n in 1:nnodes]
        nodes = [SVector{3,Int32}(10f + n for n in 1:3) for f in 1:4]
        gradients = [SVector{dim,Float64}(100n + 10p + c for c in 1:dim) for n in 1:3, p in 1:4]
        weights = Float64[0.1, 0.2, 0.3, 0.4]
        cpu = GT.MeshIntegrationData(GT.CPULayout(), coordinates, nodes, gradients, weights)
        @test GT.to_layout(cpu, GT.CPULayout()) === cpu
        gpu = GT.to_layout(cpu, GT.GPULayout())
        @test gpu.layout isa GT.GPULayout
        @test gpu.node_coordinates isa Matrix{Float32}
        @test gpu.face_nodes isa Matrix{Int32}
        @test gpu.node_coordinates == [coordinates[n][c] for n in 1:nnodes, c in 1:dim]
        @test gpu.face_nodes == [nodes[f][n] for f in 1:4, n in 1:3]
        @test gpu.reference_gradients === gradients
        @test gpu.reference_weights === weights
        @test GT.to_layout(gpu, GT.GPULayout()) === gpu
        restored = GT.to_layout(gpu, GT.CPULayout())
        @test restored.layout isa GT.CPULayout
        for name in (:node_coordinates, :face_nodes)
            @test typeof(getproperty(restored, name)) === typeof(getproperty(cpu, name))
            @test getproperty(restored, name) == getproperty(cpu, name)
        end
        @test restored.reference_gradients === gradients
        @test restored.reference_weights === weights
    end
end

@testset "Integration data construction" begin
    mesh = GT.cartesian_mesh((0, 1, 0, 1), (2, 3))
    domain = GT.interior(mesh)
    quadrature = GT.quadrature(domain, 4)
    space = GT.lagrange_space(domain, 1)
    reference_space = only(GT.reference_spaces(space))
    reference_quadrature = only(GT.reference_quadratures(quadrature))
    expected_size = (GT.num_dofs(reference_space), length(GT.weights(reference_quadrature)))
    for (f, T) in ((GT.value, Float64), (ForwardDiff.gradient, SVector{2,Float64}),
                   (ForwardDiff.hessian, SMatrix{2,2,Float64,4}))
        data = GT.space_integration_data(f, space, quadrature)
        @test data.tabulate_function === f
        @test data.reference_tabulator isa Matrix{T}
        @test size(data.reference_tabulator) == expected_size
    end
    cpu = GT.mesh_integration_data(mesh, quadrature)
    gpu = GT.to_layout(cpu, GT.GPULayout())
    @test gpu.reference_gradients === cpu.reference_gradients
    @test gpu.reference_weights === cpu.reference_weights
    restored = GT.to_layout(gpu, GT.CPULayout())
    @test restored.node_coordinates == cpu.node_coordinates
    @test restored.face_nodes == cpu.face_nodes
    @test restored.reference_gradients === cpu.reference_gradients
    @test restored.reference_weights === cpu.reference_weights
end

end # module

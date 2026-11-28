module CompilerDataTests

using Test
using StaticArrays
using LinearAlgebra
import ForwardDiff
import GalerkinToolkit as GT

# Active FEM checks cover full-dimensional, scalar CG matrix assembly.
# Preserve unsupported cases in comments for later work.
@testset "Mesh integration layouts" begin
    for dim in (0, 1, 2, 3), nnodes in (0, 5)
        coordinates = SVector{dim,Float32}[SVector{dim,Float32}(10n + c for c in 1:dim) for n in 1:nnodes]
        nodes = [SVector{3,Int32}(10f + n for n in 1:3) for f in 1:4]
        gradients = [SVector{dim,Float64}(100n + 10p + c for c in 1:dim) for n in 1:3, p in 1:4]
        weights = Float64[0.1, 0.2, 0.3, 0.4]
        cpu = GT.MeshIntegrationData(GT.ElementMajorLayout(), coordinates, nodes, gradients, weights)
        @test GT.to_layout(cpu, GT.ElementMajorLayout()) === cpu
        gpu = GT.to_layout(cpu, GT.EntryMajorLayout())
        @test gpu.layout isa GT.EntryMajorLayout
        @test gpu.node_coordinates isa Matrix{Float32}
        @test gpu.face_nodes isa Matrix{Int32}
        @test gpu.node_coordinates == [coordinates[n][c] for n in 1:nnodes, c in 1:dim]
        @test gpu.face_nodes == [nodes[f][n] for f in 1:4, n in 1:3]
        @test gpu.reference_gradients === gradients
        @test gpu.reference_weights === weights
        @test GT.to_layout(gpu, GT.EntryMajorLayout()) === gpu
        restored = GT.to_layout(gpu, GT.ElementMajorLayout())
        @test restored.layout isa GT.ElementMajorLayout
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
    gpu = GT.to_layout(cpu, GT.EntryMajorLayout())
    @test gpu.reference_gradients === cpu.reference_gradients
    @test gpu.reference_weights === cpu.reference_weights
    restored = GT.to_layout(gpu, GT.ElementMajorLayout())
    @test restored.node_coordinates == cpu.node_coordinates
    @test restored.face_nodes == cpu.face_nodes
    @test restored.reference_gradients === cpu.reference_gradients
    @test restored.reference_weights === cpu.reference_weights
end

@testset "Jacobian sum expressions" begin
    # TODO: Re-enable embedded/rectangular Jacobians.
    # for (Di, Dj) in ((1, 1), (2, 1), (3, 2), (2, 3), (3, 3)), T in (Float32, Float64)
    for (Di, Dj) in ((1, 1), (2, 2), (3, 3)), T in (Float32, Float64)
        coordinates = [SVector{Di,T}(10node + i for i in 1:Di) for node in 1:5]
        nodes = Int32[4, 1, 3]
        gradients = [SVector{Dj,T}(node - 2point + j for j in 1:Dj) for node in 1:3, point in 1:4]
        expression = GT.jacobian_sum_function_expr(Di, Dj, zero(T), length(nodes))
        compute = GT._GT_debug_eval(expression)
        for point in axes(gradients, 2)
            expected = sum(coordinates[nodes[node]] * permutedims(gradients[node,point]) for node in eachindex(nodes))
            actual = Base.invokelatest(compute, nodes, coordinates, gradients, point)
            @test actual isa SMatrix{Di,Dj,T,Di*Dj}
            @test actual == expected
        end
    end
end

@testset "Jacobian and weight expressions" begin
    boxes = ((0, 2), (0, 2, 0, 3), (0, 2, 0, 3, 0, 4))
    cells = ((3,), (2, 2), (2, 1, 1))
    warps = (
        x -> SVector(x[1] + 0.1x[1]^2),
        x -> SVector(2x[1] + 0.3x[2] + 0.1x[1]*x[2],
                     -0.2x[1] + 1.5x[2] + 0.05x[1]*x[2]),
        x -> SVector(1.5x[1] + 0.2x[2] + 0.1x[2]*x[3],
                     -0.3x[1] + 1.2x[2] + 0.05x[1]*x[3],
                     0.1x[1] + 0.2x[2] + 1.7x[3] + 0.04x[1]*x[2]),
    )
    meshes = []
    for dim in 1:3, simplex in (false, true)
        mesh = GT.cartesian_mesh(boxes[dim], cells[dim]; simplexify=simplex, boundary=false)
        coordinates = map(warps[dim], GT.node_coordinates(mesh))
        push!(meshes, GT.replace_node_coordinates(mesh, coordinates))
    end
    #= TODO: Re-enable embedded curve and surface geometry.
    mesh = GT.cartesian_mesh(boxes[1], cells[1]; boundary=false)
    coordinates = [SVector(x[1], 0.2x[1]^2) for x in GT.node_coordinates(mesh)]
    push!(meshes, GT.replace_node_coordinates(mesh, coordinates))
    mesh = GT.cartesian_mesh(boxes[2], (2, 1); boundary=false)
    coordinates = [SVector(x[1] + 0.2x[2], -0.3x[1] + 1.2x[2],
                           0.4x[1] + 0.5x[2] + 0.1x[1]*x[2]) for x in GT.node_coordinates(mesh)]
    push!(meshes, GT.replace_node_coordinates(mesh, coordinates))

    =#
    mesh = GT.cartesian_mesh(boxes[2], cells[2]; boundary=false)
    coordinates = [SVector(-x[1], x[2]) for x in GT.node_coordinates(mesh)]
    push!(meshes, GT.replace_node_coordinates(mesh, coordinates))

    for mesh in meshes
        dim = GT.num_dims(mesh)
        quadrature = GT.quadrature(GT.interior(mesh), 4)
        points = GT.coordinates(only(GT.reference_quadratures(quadrature)))
        jacobian = GT.jacobian_accessor(quadrature)
        physical_map = GT.physical_map_accessor(GT.value, mesh, dim)
        expression = GT.jacobian_expr(mesh, quadrature, :mesh_arg, :quadrature_arg, :face_arg, :point_arg)
        compute = GT._GT_debug_eval(:((mesh_arg, quadrature_arg, face_arg, point_arg) -> $expression))
        for face in GT.faces(GT.domain(quadrature))
            expected = jacobian(face)
            phi = physical_map(face)
            for point in eachindex(points)
                actual = Base.invokelatest(compute, mesh, quadrature, face, point)
                @test actual isa SMatrix{GT.num_ambient_dims(mesh),dim}
                @test actual ≈ expected(point)
                @test actual ≈ ForwardDiff.jacobian(phi, points[point])
            end
        end

        # Arguments may be compound expressions rather than bare variable names.
        expression = GT.jacobian_expr(mesh, quadrature, :(meshes_arg[1]), :(quadratures_arg[1]),
            :(face_arg + 1), :(point_arg + 1))
        compute = GT._GT_debug_eval(:((meshes_arg, quadratures_arg, face_arg, point_arg) -> $expression))
        face = last(GT.faces(GT.domain(quadrature)))
        point = lastindex(points)
        @test Base.invokelatest(compute, (mesh,), (quadrature,), face-1, point-1) ≈ jacobian(face)(point)

        @testset "Weight expressions (reference=$is_reference)" for is_reference in (false, true)
            domain = GT.interior(mesh; is_reference_domain=is_reference)
            quadrature = GT.quadrature(domain, 4)
            weights = GT.weight_accessor(quadrature)
            npoints = GT.num_points(only(GT.reference_quadratures(quadrature)))
            expression = GT.weight_expr(mesh, quadrature, :mesh_arg, :quadrature_arg, :face_arg, :point_arg)
            compute = GT._GT_debug_eval(:((mesh_arg, quadrature_arg, face_arg, point_arg) -> $expression))
            for face in GT.faces(domain)
                expected = weights(face)
                for point in 1:npoints
                    actual = Base.invokelatest(compute, mesh, quadrature, face, point)
                    @test actual ≈ expected(point)
                end
            end

            expression = GT.weight_expr(mesh, quadrature, :(meshes_arg[1]), :(quadratures_arg[1]),
                :(face_arg + 1), :(point_arg + 1))
            compute = GT._GT_debug_eval(:((meshes_arg, quadratures_arg, face_arg, point_arg) -> $expression))
            face = last(GT.faces(domain))
            point = npoints
            @test Base.invokelatest(compute, (mesh,), (quadrature,), face-1, point-1) ≈ weights(face)(point)
        end
    end
end

#= TODO: Re-enable rectangular and complex Jacobian pseudoinverses.
@testset "Jacobian pseudoinverses" begin
    for (M, N) in ((2, 1), (3, 1), (3, 2)), T in (Float32, Float64), sample in 1:4
        J = SMatrix{M,N,T}(ntuple(k -> T(sin(k + sample) +
            (mod1(k,M) == cld(k,M) ? 2 : 0)), M*N))
        expected = pinv(J)
        for matrix in (J, MMatrix{M,N,T}(J), SizedMatrix{M,N}(Matrix(J)))
            actual = @inferred GT.jacobian_pinv(matrix)
            @test actual isa SMatrix{N,M,T}
            @test actual ≈ expected
            @test actual * matrix ≈ SMatrix{N,N,T}(I)
            projector = matrix * actual
            @test projector ≈ transpose(projector)
        end
    end

    # Nearly parallel columns whose Gram determinant rounds to zero.
    for T in (Float32, Float64)
        delta = T === Float32 ? T(1e-4) : T(1e-9)
        J = SMatrix{3,2,T}(1, 0, 0, 1, delta, 0)
        actual = GT.jacobian_pinv(J)
        @test actual ≈ pinv(J)
        @test actual * J ≈ SMatrix{2,2,T}(I)
    end

    # Other shapes and complex matrices retain the general pseudoinverse.
    for (M, N, T) in ((4, 2, Float64), (3, 0, Float32), (2, 3, Float64), (3, 2, ComplexF64))
        J = SMatrix{M,N,T}(ntuple(k -> T(sin(k) +
            (mod1(k,M) == cld(k,M) ? 2 : 0)) +
            (T <: Complex ? T(0.2im*cos(k)) : zero(T)), M*N))
        @test GT.jacobian_pinv(J) ≈ pinv(J)
    end
end

=#
@testset "Shape function expressions" begin
    base_mesh = GT.cartesian_mesh((0, 2, 0, 3), (2, 1); boundary=false)
    # TODO: Re-enable embedded geometry.
    # for embedded in (false, true)
    for embedded in (false,)
        coordinates = map(GT.node_coordinates(base_mesh)) do x
            y = SVector(2x[1] + 0.3x[2] + 0.1x[1]*x[2], -0.2x[1] + 1.5x[2])
            embedded ? SVector(y[1], y[2], 0.4x[1] + 0.5x[2]) : y
        end
        mesh = GT.replace_node_coordinates(base_mesh, coordinates)
        quadrature = GT.quadrature(GT.interior(mesh), 4)
        npoints = GT.num_points(only(GT.reference_quadratures(quadrature)))
        for is_reference in (false, true)
            domain = GT.interior(mesh; is_reference_domain=is_reference)
            scalar_space = GT.lagrange_space(domain, 1)
            # TODO: Re-enable vector-valued spaces.
            # vector_space = GT.lagrange_space(domain, 1; tensor_size=Val((2,)))
            for (func, space) in ((GT.value, scalar_space), (ForwardDiff.gradient, scalar_space),
                                  # (GT.value, vector_space), (ForwardDiff.jacobian, vector_space)
                                  )
                accessor = GT.shape_function_accessor(func, space, quadrature)
                ndofs = GT.num_dofs(only(GT.reference_spaces(space)))
                expression = GT.shape_function_expr(func, mesh, space, quadrature,
                    :(meshes_arg[1]), :(spaces_arg[1]), :(quadratures_arg[1]),
                    :(face_arg + 1), :(point_arg + 1), :(dof_arg + 1))
                compute = GT._GT_debug_eval(:((meshes_arg, spaces_arg, quadratures_arg,
                    face_arg, point_arg, dof_arg) -> $expression))
                for face in GT.faces(GT.domain(quadrature))
                    expected_point = accessor(face)
                    for point in 1:npoints
                        expected = expected_point(point)
                        for dof in 1:ndofs
                            actual = Base.invokelatest(compute, (mesh,), (space,), (quadrature,),
                                face-1, point-1, dof-1)
                            @test actual ≈ expected(dof)
                        end
                    end
                end
            end
        end
    end
end

@testset "General reference shape expressions" begin
    for simplex in (false, true)
        mesh = GT.cartesian_mesh((0, 2, 0, 3), (2, 1); simplexify=simplex, boundary=false)
        domain = GT.interior(mesh; is_reference_domain=true)
        # Reference evaluation also accepts a quadrature on the physical domain.
        quadrature = GT.quadrature(GT.interior(mesh), 4)
        points = GT.coordinates(only(GT.reference_quadratures(quadrature)))
        scalar_space = GT.lagrange_space(domain, 2)
        # TODO: Re-enable vector-valued shape operators.
        # vector_space = GT.lagrange_space(domain, 2; tensor_size=Val((2,)))
        square_value = (s, x) -> s(x)^2
        operators = (
            ((s, x) -> ForwardDiff.gradient(s, x), scalar_space),
            (ForwardDiff.hessian, scalar_space),
            ((s, x) -> tr(ForwardDiff.hessian(s, x)), scalar_space),
            (square_value, scalar_space),
            ((s, x) -> s(x + SVector(0.05, -0.03)), scalar_space),
            # ((s, x) -> tr(ForwardDiff.jacobian(s, x)), vector_space),
            # ((s, x) -> sum(abs2, s(x)), vector_space),
        )
        for (func, space) in operators
            # The untabulated accessor applies func to the complete shape.
            accessor = GT.shape_function_accessor_reference(func, space)
            ndofs = GT.num_dofs(only(GT.reference_spaces(space)))
            data = GT.space_integration_data(func, space, quadrature)
            @test data.tabulate_function === func
            @test size(data.reference_tabulator) == (ndofs, length(points))
            for point in eachindex(points), dof in 1:ndofs
                @test data.reference_tabulator[dof,point] ≈ accessor(1)(dof)(points[point])
            end
            expression = GT.shape_function_expr(func, mesh, space, quadrature,
                :(meshes_arg[1]), :(spaces_arg[1]), :(quadratures_arg[1]),
                :(face_arg + 1), :(point_arg + 1), :(dof_arg + 1))
            compute = GT._GT_debug_eval(:((meshes_arg, spaces_arg, quadratures_arg,
                face_arg, point_arg, dof_arg) -> $expression))
            for face in GT.faces(domain), (point, dof) in ((1, 1), (lastindex(points), ndofs))
                actual = Base.invokelatest(compute, (mesh,), (space,), (quadrature,),
                    face-1, point-1, dof-1)
                @test actual ≈ data.reference_tabulator[dof,point]
            end
        end

        # Tabulated values provide an independent reference for nonlinear evaluation.
        reference_space = only(GT.reference_spaces(scalar_space))
        values = GT.tabulator(reference_space)(GT.value, points)
        squared = GT.space_integration_data(square_value, scalar_space, quadrature)
        @test squared.reference_tabulator ≈ permutedims(values .^ 2)

        # A prepared table is a snapshot; changed parameters require rebuilding it.
        scale = Ref(2.0)
        calls = Ref(0)
        func = (s, x) -> begin
            calls[] += 1
            scale[] * s(x)^2
        end
        data = GT.space_integration_data(func, scalar_space, quadrature)
        @test calls[] == length(values)
        @test data.reference_tabulator[1,1] ≈ 2values[1,1]^2
        scale[] = 3.0
        @test data.reference_tabulator[1,1] ≈ 2values[1,1]^2
        refreshed = GT.space_integration_data(func, scalar_space, quadrature)
        @test refreshed.reference_tabulator[1,1] ≈ 3values[1,1]^2
        @test calls[] == 2length(values)
        calls[] = 0
        expression = GT.shape_function_expr(func, mesh, scalar_space, quadrature,
            :mesh_arg, :space_arg, :quadrature_arg, :face_arg, :point_arg, :dof_arg)
        compute = GT._GT_debug_eval(:((mesh_arg, space_arg, quadrature_arg,
            face_arg, point_arg, dof_arg) -> $expression))
        @test calls[] == 0
        @test Base.invokelatest(compute, mesh, scalar_space, quadrature, 1, 1, 1) ≈ 3values[1,1]^2
        @test calls[] == length(values)
    end
end

#= TODO: Re-enable general physical callable/inverse-map operators, vector spaces and embedded geometry.
@testset "General physical shape expressions" begin
    for dim in 1:3, simplex in (false, true)
        box = ntuple(i -> isodd(i) ? 0 : 2, 2dim)
        cells = ntuple(i -> i == 1 ? 2 : 1, dim)
        base_mesh = GT.cartesian_mesh(box, cells; simplexify=simplex, boundary=false)
        A = SMatrix{dim,dim,Float64}(ntuple(dim*dim) do k
            i, j = mod1(k, dim), cld(k, dim)
            i == j ? (i == 1 ? -1.5 : 1.2) : 0.1(i+j)
        end)
        b = SVector{dim,Float64}(0.5i for i in 1:dim)
        # Each cell is affine, with different scaling along the first axis.
        coordinates = map(GT.node_coordinates(base_mesh)) do x
            A * setindex(x, x[1] + 0.15x[1]^2, 1) + b
        end
        mesh = GT.replace_node_coordinates(base_mesh, coordinates)
        domain = GT.interior(mesh)
        # More quadrature points than geometry nodes exercise the tabulation axes.
        quadrature = GT.quadrature(domain, 4)
        points = GT.coordinates(only(GT.reference_quadratures(quadrature)))
        scalar_space = GT.lagrange_space(domain, 2)
        vector_space = GT.lagrange_space(domain, 2; tensor_size=Val((dim,)))
        scalar_values = GT.shape_function_accessor(GT.value, scalar_space, quadrature)
        scalar_gradients = GT.shape_function_accessor(ForwardDiff.gradient, scalar_space, quadrature)
        vector_values = GT.shape_function_accessor(GT.value, vector_space, quadrature)
        vector_jacobians = GT.shape_function_accessor(ForwardDiff.jacobian, vector_space, quadrature)
        jacobians = GT.jacobian_accessor(quadrature)
        physical_map = GT.physical_map_accessor(GT.value, mesh, Val(dim))
        reference_space = only(GT.reference_spaces(scalar_space))
        reference_hessians = GT.tabulator(reference_space)(ForwardDiff.hessian, points)
        reference_shapes = GT.shape_functions(reference_space)
        shift = SVector{dim,Float64}(0.03i for i in 1:dim)
        operators = (
            ((s, x) -> s(x), scalar_space, :value),
            ((s, x) -> ForwardDiff.gradient(s, x), scalar_space, :gradient),
            (ForwardDiff.hessian, scalar_space, :hessian),
            ((s, x) -> tr(ForwardDiff.hessian(s, x)), scalar_space, :laplacian),
            ((s, x) -> s(x)^2 + sum(x), scalar_space, :nonlinear),
            ((s, x) -> s(x + shift), scalar_space, :shifted),
            ((s, x) -> ForwardDiff.jacobian(s, x), vector_space, :jacobian),
            ((s, x) -> tr(ForwardDiff.jacobian(s, x)), vector_space, :divergence),
            ((s, x) -> sum(abs2, s(x)), vector_space, :vector_nonlinear),
        )
        for (func, space, kind) in operators
            ndofs = GT.num_dofs(only(GT.reference_spaces(space)))
            expression = GT.shape_function_expr(func, mesh, space, quadrature,
                :(meshes_arg[1]), :(spaces_arg[1]), :(quadratures_arg[1]),
                :(face_arg + 1), :(point_arg + 1), :(dof_arg + 1))
            compute = GT._GT_debug_eval(:((meshes_arg, spaces_arg, quadratures_arg,
                face_arg, point_arg, dof_arg) -> $expression))
            for face in GT.faces(domain), point in unique((1, lastindex(points)))
                invJ = inv(jacobians(face)(point))
                x = physical_map(face)(points[point])
                for dof in 1:ndofs
                    expected = if kind === :value
                        scalar_values(face)(point)(dof)
                    elseif kind === :gradient
                        scalar_gradients(face)(point)(dof)
                    elseif kind === :hessian || kind === :laplacian
                        H = transpose(invJ) * reference_hessians[point,dof] * invJ
                        kind === :hessian ? H : tr(H)
                    elseif kind === :nonlinear
                        scalar_values(face)(point)(dof)^2 + sum(x)
                    elseif kind === :shifted
                        reference_shapes[dof](points[point] + invJ*shift)
                    elseif kind === :jacobian
                        vector_jacobians(face)(point)(dof)
                    elseif kind === :divergence
                        tr(vector_jacobians(face)(point)(dof))
                    else
                        sum(abs2, vector_values(face)(point)(dof))
                    end
                    actual = Base.invokelatest(compute, (mesh,), (space,), (quadrature,),
                        face-1, point-1, dof-1)
                    @test actual ≈ expected atol=1e-12
                end
            end
        end
    end

    # Face arguments are global mesh IDs, including for a domain subset.
    subset_mesh = GT.cartesian_mesh((0, 2, 0, 1), (2, 1); boundary=false)
    GT.group_faces(subset_mesh, 2)["selected_cell"] = Int32[2]
    subset_domain = GT.interior(subset_mesh; group_names=["selected_cell"])
    subset_quadrature = GT.quadrature(subset_domain, 2)
    subset_space = GT.lagrange_space(subset_domain, 2)
    subset_points = GT.coordinates(only(GT.reference_quadratures(subset_quadrature)))
    subset_values = GT.space_integration_data(GT.value, subset_space, subset_quadrature)
    phi = GT.physical_map_accessor(GT.value, subset_mesh, Val(2))(2)
    func = (s, x) -> s(x)^2 + sum(x)
    expression = GT.shape_function_expr(func, subset_mesh, subset_space, subset_quadrature,
        :mesh_arg, :space_arg, :quadrature_arg, :face_arg, :point_arg, :dof_arg)
    compute = GT._GT_debug_eval(:((mesh_arg, space_arg, quadrature_arg,
        face_arg, point_arg, dof_arg) -> $expression))
    for point in eachindex(subset_points), dof in axes(subset_values.reference_tabulator, 1)
        expected = subset_values.reference_tabulator[dof,point]^2 + sum(phi(subset_points[point]))
        @test Base.invokelatest(compute, subset_mesh, subset_space, subset_quadrature,
            2, point, dof) ≈ expected
    end

    # The general callable fallback currently requires full-dimensional geometry.
    base_mesh = GT.cartesian_mesh((0, 1, 0, 1), (1, 1); boundary=false)
    embedded = GT.replace_node_coordinates(base_mesh,
        [SVector(x[1], x[2], 0.3x[1]) for x in GT.node_coordinates(base_mesh)])
    quadrature = GT.quadrature(GT.interior(embedded), 2)
    space = GT.lagrange_space(GT.interior(embedded), 2)
    @test_throws ArgumentError GT.shape_function_expr(ForwardDiff.hessian,
        embedded, space, quadrature, :mesh_arg, :space_arg, :quadrature_arg,
        :face_arg, :point_arg, :dof_arg)
end

=#
@testset "Change of measure expressions" begin
    # TODO: Re-enable point and embedded measures.
    # shapes = ((0, 0), (1, 0), (3, 0), (1, 1), (2, 1), (2, 2), (3, 1), (3, 2), (3, 3))
    shapes = ((1, 1), (2, 2), (3, 3))
    for (M, N) in shapes
        expression = GT.change_of_measure_expr(:J_arg, M, N)
        compute = GT._GT_debug_eval(:(J_arg -> $expression))
        expression = GT.change_of_measure_expr(:(matrices_arg[1]), M, N)
        compute_compound = GT._GT_debug_eval(:(matrices_arg -> $expression))
        for T in (Float32, Float64)
            J = SMatrix{M,N,T}(ntuple(k -> T(sin(k) +
                (mod1(k,M) == cld(k,M) ? 2 : 0)), M*N))
            @test Base.invokelatest(compute, J) ≈ GT.change_of_measure(J)
            @test Base.invokelatest(compute_compound, (J,)) ≈ GT.change_of_measure(J)
        end
    end
end

# Keep the original formula as an independent numerical baseline.
original_change_of_measure(J) = sqrt(det(transpose(J)*J))

@testset "Change of measure" begin
    # TODO: Re-enable point and embedded measures.
    # shapes = ((0, 0), (1, 0), (3, 0), (1, 1), (2, 1), (2, 2),
    #           (3, 1), (3, 2), (3, 3))
    shapes = ((1, 1), (2, 2), (3, 3))
    for T in (Float32, Float64), (M, N) in shapes, sample in 1:4
        J = SMatrix{M,N,T}(ntuple(k -> T(sin(0.1sample + k) +
            (mod1(k,M) == cld(k,M) ? 2 : 0)), M*N))
        # Singular values give the volume independently of either formula.
        expected = prod(svdvals(Matrix{Float64}(J)))
        for matrix in (J, Matrix(J), MMatrix{M,N,T,M*N}(J), SizedMatrix{M,N}(Matrix(J)))
            actual = @inferred GT.change_of_measure(matrix)
            @test actual ≈ expected rtol=32eps(T)
            @test actual ≈ original_change_of_measure(matrix) rtol=32eps(T)
        end
    end

    for T in (Float32, Float64)
        reversed = SMatrix{2,2,T}(-2, 0, 0, 3)
        @test det(reversed) < 0
        @test GT.change_of_measure(reversed) == T(6)
        @test GT.change_of_measure(Matrix(reversed)) == T(6)
        #= TODO: Re-enable surface measures.
        @test GT.change_of_measure(SMatrix{3,2,T}(1, 0, 1, 0, 1, 1)) ≈ sqrt(T(3))
        @test GT.change_of_measure(SMatrix{3,2,T}(1, 2, 3, 2, 4, 6)) == zero(T)
        =#
        for (M, N) in shapes
            expected = N == 0 ? one(T) : zero(T)
            @test GT.change_of_measure(zeros(SMatrix{M,N,T,M*N})) == expected
        end
        #= TODO: Re-enable nearly parallel surface measures.
        # The Gram determinant loses this small area to cancellation.
        delta = sqrt(eps(T))/4
        nearly_parallel = SMatrix{3,2,T}(1, 0, 0, 1, delta, 0)
        @test GT.change_of_measure(nearly_parallel) ≈ delta
        =#
    end

    for scale in (1.0e-100, 1.0e100)
        expected = scale*scale
        square = SMatrix{2,2}(scale, 0, 0, scale)
        # surface = SMatrix{3,2}(scale, 0, 0, 0, scale, 0)
        @test GT.change_of_measure(square) ≈ expected
        # @test GT.change_of_measure(surface) ≈ expected
    end

    # for (M, N) in ((2, 2), (3, 3), (3, 2)), orientation in (-1, 1)
    for (M, N) in ((2, 2), (3, 3)), orientation in (-1, 1)
        J = SMatrix{M,N}(ntuple(k -> sin(k) + (mod1(k,M) == cld(k,M) ? 2 : 0), M*N))
        coordinates = vec(Matrix(J))
        coordinates[1:M] .*= orientation
        optimized = ForwardDiff.gradient(x -> GT.change_of_measure(SMatrix{M,N}(x)), coordinates)
        original = ForwardDiff.gradient(x -> original_change_of_measure(SMatrix{M,N}(x)), coordinates)
        @test optimized ≈ original
    end
end

@testset "Tabulated Lagrange discrete fields" begin
    ∇ = ForwardDiff.gradient
    no_face_around = GT.leaf_term(nothing)
    forms(f,gf,dμ) = (
        (u,v)->GT.∫(x->f(x)*u(x)*v(x) + dot(gf(x),∇(u,x))*v(x),dμ),
        v->GT.∫(x->f(x)*v(x) + dot(gf(x),∇(v,x)),dμ),
        ()->sum(GT.∫(x->f(x)^2 + dot(gf(x),gf(x)),dμ)))

    # Affine cells: an order-k Lagrange field reproduces a degree-k polynomial, so the
    # array-read field and its physical gradient equal the analytical ones at every point.
    for (dim,simplexify,k,T,dirichlet) in ((2,false,1,Float64,true),(2,false,2,ComplexF64,false),
            (2,true,2,Float64,true),(2,true,1,ComplexF64,true),(3,true,1,Float64,false),
            (3,true,2,ComplexF64,true))
        mesh = GT.cartesian_mesh(ntuple(i->isodd(i) ? 0 : 1,2dim),ntuple(_->2,dim);simplexify)
        mesh = GT.replace_node_coordinates(mesh,
            [SVector{dim}(ntuple(i->1.2x[i]+0.1(sum(x)-x[i]),dim)) for x in GT.node_coordinates(mesh)])
        Ω = GT.interior(mesh)
        c = T <: Complex ? T(1,2) : one(T)
        q = k == 2 ? 1.0 : 0.0
        p = x -> c*(1 + x[1] - x[dim]/2 + q*(x[1]*x[dim] - 3x[1]^2/10))
        gp = x -> c*SVector{dim}(ntuple(i -> i == 1 ? 1 + q*(x[dim] - 3x[1]/5) :
            i == dim ? -1/2 + q*x[1] : 0.0,dim))
        Pf = dirichlet ? GT.lagrange_space(Ω,k;dirichlet_boundary=GT.boundary(mesh)) :
            GT.lagrange_space(Ω,k)
        uh = GT.interpolate(GT.analytical_field(p,Ω),Pf)
        @test eltype(GT.free_values(uh)) == T
        @test (GT.num_dirichlet_dofs(Pf) > 0) == dirichlet
        @test GT.tabulated_lagrange_discrete_field(GT.value,uh,no_face_around,Ω)
        @test GT.tabulated_lagrange_discrete_field(ForwardDiff.gradient,uh,no_face_around,Ω)
        @test !GT.tabulated_lagrange_discrete_field(ForwardDiff.hessian,uh,no_face_around,Ω)
        @test !GT.tabulated_lagrange_discrete_field(GT.value,uh,GT.leaf_term(1),Ω)
        @test !GT.tabulated_lagrange_discrete_field(GT.value,uh,no_face_around,GT.boundary(mesh))
        u_exact,g_exact = GT.analytical_field(p,Ω),GT.analytical_field(gp,Ω)
        V = GT.lagrange_space(Ω,1)
        dμ = GT.measure(Ω,2k+2)
        a,l,s = forms(uh,x->∇(uh,x),dμ)
        a_ref,l_ref,s_ref = forms(u_exact,g_exact,dμ)
        @test GT.assemble_matrix(a,T,V,V) ≈ GT.assemble_matrix(a_ref,T,V,V) rtol=1e-12
        @test GT.assemble_vector(l,T,V) ≈ GT.assemble_vector(l_ref,T,V) rtol=1e-12
        @test s() ≈ s_ref() rtol=1e-12

        # A vector-valued field that is affine in x, with a nonsymmetric Jacobian A, reproduced exactly.
        # The contraction with a nonsymmetric C sees every entry of the Jacobian (a transposed one fails).
        A = SMatrix{dim,dim}([i == j ? 1.0 : i < j ? 0.5i+j : -0.25j for i in 1:dim, j in 1:dim])
        C = SMatrix{dim,dim}([i + 3j/4 for i in 1:dim, j in 1:dim])
        @test abs(dot(c*A,C) - dot(transpose(c*A),C)) > 0.1
        pv = x -> c*(SVector{dim}(ntuple(i->i/2,dim)) + A*x)
        Pv = dirichlet ? GT.lagrange_space(Ω,k;tensor_size=Val((dim,)),dirichlet_boundary=GT.boundary(mesh)) :
            GT.lagrange_space(Ω,k;tensor_size=Val((dim,)))
        uhv = GT.interpolate(GT.analytical_field(pv,Ω),Pv)
        @test GT.tabulated_lagrange_discrete_field(ForwardDiff.jacobian,uhv,no_face_around,Ω)
        κ = dot(c*A,C)
        Cq = GT.uniform_quantity(C)
        lv = v->GT.∫(x->dot(ForwardDiff.jacobian(uhv,x),Cq)*v(x),dμ)
        lv_ref = v->GT.∫(x->κ*v(x),dμ)
        @test GT.assemble_vector(lv,T,V) ≈ GT.assemble_vector(lv_ref,T,V) rtol=1e-12
    end

    # Free and Dirichlet values of different element types: each coefficient is multiplied by the table entry
    # as stored, as in the accessor path. Converting the Int64 free values to Float32 (the promotion of the two
    # value types) would turn 16777217 into 16777216.
    mesh = GT.cartesian_mesh((0,1,0,1),(2,2))
    Ω = GT.interior(mesh)
    V = GT.lagrange_space(Ω,1)
    dμ = GT.measure(Ω,2)
    for Pf in (GT.lagrange_space(Ω,1),GT.lagrange_space(Ω,1;dirichlet_boundary=GT.boundary(mesh)))
        big = 16_777_217
        free = fill(big,GT.num_free_dofs(Pf))
        diri = fill(0.5f0,GT.num_dirichlet_dofs(Pf))
        uh = GT.discrete_field(Pf,free,diri)
        uv = GT.discrete_field(Pf,view(free,:),view(diri,:))
        @test GT.tabulated_lagrange_discrete_field(GT.value,uh,no_face_around,Ω)
        @test !GT.tabulated_lagrange_discrete_field(GT.value,uv,no_face_around,Ω)
        s = f->sum(GT.∫(x->f(x),dμ))
        @test s(uh) ≈ s(uv) rtol=1e-14
        GT.num_dirichlet_dofs(Pf) == 0 && @test s(uh) ≈ big rtol=1e-14
        a = f->((u,v)->GT.∫(x->f(x)*u(x)*v(x),dμ))
        @test GT.assemble_matrix(a(uh),Float64,V,V) ≈ GT.assemble_matrix(a(uv),Float64,V,V) rtol=1e-14
    end

    # Non-affine cells and a subset of the cells: the array-read field (space on all cells)
    # against the same field on a space over the subset, which keeps the accessor path, and
    # against a field whose values are a view, which also keeps the accessor path.
    for (dim,simplexify,k) in ((2,false,2),(3,false,1),(3,true,2))
        mesh = GT.cartesian_mesh(ntuple(i->isodd(i) ? 0 : 1,2dim),ntuple(_->2,dim);simplexify)
        mesh = GT.replace_node_coordinates(mesh,
            [SVector{dim}(ntuple(i->x[i]+0.08*sin(2x[mod1(i+1,dim)])+0.05*x[i]*x[mod1(i+1,dim)],dim))
             for x in GT.node_coordinates(mesh)])
        ncells = GT.num_faces(mesh,dim)
        GT.group_faces(mesh,dim)["tabulated field subset"] = Int32[1,ncells÷2,ncells]
        Ω = GT.interior(mesh)
        Ωs = GT.interior(mesh;group_names=["tabulated field subset"])
        g = x -> 1 + x[1]^2 - x[1]*x[dim]/2 + 3sin(3x[dim])/10
        Pf = GT.lagrange_space(Ω,k;dirichlet_boundary=GT.boundary(mesh))
        uh = GT.interpolate(GT.analytical_field(g,Ω),Pf)
        us = GT.interpolate(GT.analytical_field(g,Ωs),GT.lagrange_space(Ωs,k))
        uv = GT.discrete_field(Pf,view(GT.free_values(uh),:),GT.dirichlet_values(uh))
        @test GT.tabulated_lagrange_discrete_field(GT.value,uh,no_face_around,Ωs)
        @test !GT.tabulated_lagrange_discrete_field(GT.value,us,no_face_around,Ωs)
        @test !GT.tabulated_lagrange_discrete_field(GT.value,uv,no_face_around,Ωs)
        V = GT.lagrange_space(Ω,1)
        dμ = GT.measure(Ωs,2k+2)
        a,l,s = forms(uh,x->∇(uh,x),dμ)
        for f in (us,uv)
            a_ref,l_ref,s_ref = forms(f,x->∇(f,x),dμ)
            @test GT.assemble_matrix(a,Float64,V,V) ≈ GT.assemble_matrix(a_ref,Float64,V,V) rtol=1e-12
            @test GT.assemble_vector(l,Float64,V) ≈ GT.assemble_vector(l_ref,Float64,V) rtol=1e-12
            @test s() ≈ s_ref() rtol=1e-12
        end

        # Free values held in a Fill (dirichlet_field) take the array-read path; the reference holds the same
        # values in views, which keep the accessor path.
        ud = GT.dirichlet_field(Pf,copy(GT.dirichlet_values(uh)))
        @test !(GT.free_values(ud) isa DenseVector)
        udv = GT.discrete_field(Pf,view(collect(GT.free_values(ud)),:),view(copy(GT.dirichlet_values(ud)),:))
        @test GT.tabulated_lagrange_discrete_field(GT.value,ud,no_face_around,Ωs)
        @test !GT.tabulated_lagrange_discrete_field(GT.value,udv,no_face_around,Ωs)
        ad,ld,sd = forms(ud,x->∇(ud,x),dμ)
        ad_ref,ld_ref,sd_ref = forms(udv,x->∇(udv,x),dμ)
        @test GT.assemble_matrix(ad,Float64,V,V) ≈ GT.assemble_matrix(ad_ref,Float64,V,V) rtol=1e-12
        @test GT.assemble_vector(ld,Float64,V) ≈ GT.assemble_vector(ld_ref,Float64,V) rtol=1e-12
        @test sd() ≈ sd_ref() rtol=1e-12

        # Vector-valued fields: value and Jacobian, against the same field with view values.
        J = ForwardDiff.jacobian
        gv = x -> SVector{dim}(ntuple(i->g(x)+i*x[i],dim))
        Pv = GT.lagrange_space(Ω,k;tensor_size=Val((dim,)),dirichlet_boundary=GT.boundary(mesh))
        uhv = GT.interpolate(GT.analytical_field(gv,Ω),Pv)
        uvv = GT.discrete_field(Pv,view(GT.free_values(uhv),:),GT.dirichlet_values(uhv))
        @test GT.tabulated_lagrange_discrete_field(J,uhv,no_face_around,Ωs)
        @test !GT.tabulated_lagrange_discrete_field(J,uvv,no_face_around,Ωs)
        Cq = GT.uniform_quantity(SMatrix{dim,dim}([i + 3j/4 for i in 1:dim, j in 1:dim]))
        vforms(f) = ((u,v)->GT.∫(x->dot(f(x),f(x))*u(x)*v(x) + tr(J(f,x))*dot(∇(u,x),∇(v,x)),dμ),
                     ()->sum(GT.∫(x->tr(J(f,x))^2 + dot(f(x),f(x)) + dot(J(f,x),Cq),dμ)))
        av,sv = vforms(uhv)
        av_ref,sv_ref = vforms(uvv)
        @test GT.assemble_matrix(av,Float64,V,V) ≈ GT.assemble_matrix(av_ref,Float64,V,V) rtol=1e-12
        @test sv() ≈ sv_ref() rtol=1e-12
    end
end

# Timing experiments live outside the correctness suite in an ignored local file.

end # module

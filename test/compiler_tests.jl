module CompilerTests

using Test
using StaticArrays
using LinearAlgebra
import ForwardDiff
import GalerkinToolkit as GT

function compile_term(term, parameters...; arguments=parameters, forbidden_accessors=(), hoist=false)
    term = GT.lambda_term(term,:face,:point,:dof,:field,:face_around)
    term = GT.parametrize(term,parameters...)
    term,data = GT.capture(term)
    expr = GT.expression(term)
    for accessor in forbidden_accessors
        @test !occursin(accessor,string(expr))
    end
    f = GT.evaluate(hoist ? GT.statements_expr(expr) : expr,data)
    Base.invokelatest(f,arguments...)
end

scaled_map_operator(scale) = (phi,xi) -> scale*sum(abs2,phi(xi))

@testset "DOF data and metadata-selected lookups" begin
    mesh = GT.cartesian_mesh((0,2,0,2),(2,2);simplexify=false)
    omega = GT.interior(mesh)
    boundary = GT.boundary(mesh)
    skeleton = GT.skeleton(mesh)
    GT.group_faces(mesh,2)["selected cells"] = Int32[2,4]
    selected = GT.interior(mesh;group_names=["selected cells"])
    second_side = GT.replace_faces_around(skeleton,fill(2,GT.num_faces(skeleton)))
    for continuous in (true,false)
        V = GT.lagrange_space(omega,2;continuous=Val(continuous),
            dirichlet_boundary=continuous ? boundary : nothing)
        data = GT.DOFData(V)
        @test GT.face_dofs(data) === GT.face_dofs(V)
        if continuous
            @test any(dof -> dof < 0,Iterators.flatten(GT.face_dofs(data)))
        end
        for domain in (selected,skeleton,boundary,second_side)
            expr = GT.dofs_expr(V,domain,:V,:mesh,:face,:side)
            @test !occursin("dofs_accessor",string(expr))
            needs_incidence = GT.num_dims(domain) < GT.num_dims(V)
            @test occursin("face_incidence_data",string(expr)) == needs_incidence
            if needs_incidence
                incidence = GT.face_incidence_data(mesh,GT.num_dims(domain),GT.num_dims(V))
                @test GT.face_cells(incidence) === GT.face_incidence(GT.topology(mesh),GT.num_dims(domain),GT.num_dims(V))
            end
            kernel = GT.evaluate(:((V,mesh) -> ((face,side) -> $expr)),(V,mesh))
            Q = GT.quadrature(domain,4)
            domain_data = GT.domain_integration_data(domain)
            @test GT.domain_faces(domain_data) === GT.faces(domain)
            @test fieldnames(typeof(domain_data)) == (:domain_faces,)
            field_expr = GT.field_dofs_expr(V,domain,:V,:Q,1,:face,:side)
            @test !occursin("GT.faces(",string(field_expr))
            field_kernel = GT.evaluate(:((V,Q) -> ((face,side) -> $field_expr)),(V,Q))
            expected = GT.dofs_accessor(V,domain)
            for face in 1:GT.num_faces(domain), side in 1:2
                selected_side = GT.faces_around(domain) === nothing ? side : GT.faces_around(domain)[face]
                mesh_face = GT.faces(domain)[face]
                @test Base.invokelatest(kernel,mesh_face,selected_side) == expected(face,side)
                @test Base.invokelatest(field_kernel,face,side) == expected(face,side)
            end
        end
    end
end

@testset "Metadata-specialized Jacobian callable" begin
    for (Di,Dj) in ((1,1),(2,1),(3,2),(2,3),(3,3)), T in (Float32,Float64),
            nlnodes in (2,3), simplex in (false,true)
        z = zero(T)
        compute = GT.JacobianSumCallable{Di,Dj,z,nlnodes,simplex}()
        @test isbits(compute)
        coordinates = [SVector{Di,T}(10node+i for i in 1:Di) for node in 1:5]
        nodes = Int32[4,1,3]
        gradients = [SVector{Dj,T}(node-2point+j for j in 1:Dj) for node in 1:3, point in 1:4]
        for point in axes(gradients,2)
            # The static node count controls the sum even with a larger input buffer.
            expected = sum(coordinates[nodes[node]]*transpose(gradients[node,point]) for node in 1:nlnodes)
            actual = compute(nodes,coordinates,gradients,point)
            @test actual isa SMatrix{Di,Dj,T,Di*Dj}
            @test actual == expected
        end
    end
end

@testset "Hoisted closure cache isolation" begin
    add = GT.statements_expr(:(() -> (data -> (x -> data+x))))
    multiply = GT.statements_expr(:(() -> (data -> (x -> data*x))))
    @test add.head === :(->)
    @test multiply.head === :(->)
    add_parameters = GT.evaluate(add,())
    multiply_parameters = GT.evaluate(multiply,())
    cached_add_parameters = GT.evaluate(add,())
    add_kernel = Base.invokelatest(add_parameters,2)
    multiply_kernel = Base.invokelatest(multiply_parameters,3)
    cached_add_kernel = Base.invokelatest(cached_add_parameters,4)
    @test Base.invokelatest(add_kernel,1) == 3
    @test Base.invokelatest(multiply_kernel,2) == 6
    @test Base.invokelatest(cached_add_kernel,1) == 5
end

@testset "Integration expression metadata" begin
    for simplex in (false,true), reference in (false,true)
        base = GT.cartesian_mesh((0,2,0,3),(2,1);boundary=false,simplexify=simplex)
        mesh = GT.replace_node_coordinates(base,
            [SVector(x[1]+0.2x[2],x[2]+0.1x[1]*x[2]) for x in GT.node_coordinates(base)])
        GT.group_faces(mesh,2)["selected"] = Int32[2]
        dom = GT.interior(mesh;group_names=["selected"],is_reference_domain=reference)
        Q = GT.quadrature(dom,2)
        V = GT.lagrange_space(dom,1)
        deps = map(GT.leaf_term,(Q,:face,:point))
        coordinate = compile_term(GT.CoordinateTerm(deps),Q)
        weight_term = GT.WeightTerm(deps,2)
        weight = compile_term(weight_term,Q)
        coordinates = GT.coordinate_accessor(Q)(1)
        weights = GT.weight_accessor(Q)(1)
        npoints = GT.num_points(only(GT.reference_quadratures(Q)))
        for point in 1:npoints
            @test Base.invokelatest(coordinate,1,point,1,1,1) ≈ coordinates(point)
            @test Base.invokelatest(weight,1,point,1,1,1) ≈ weights(point)
        end

        coordinates_term = GT.CoordinateTerm(deps)
        map_deps = map(GT.leaf_term,(GT.value,mesh,Val(2),:face))
        physical_map = GT.PhysicalMapTerm(map_deps)
        operators = (GT.value,ForwardDiff.jacobian,
            (phi,xi) -> sum(abs2,phi(xi)),
            (phi,xi) -> ForwardDiff.hessian(q -> phi(q)[2],xi),
            (phi,xi) -> phi(xi+SVector(0.03,-0.02)))
        points = GT.coordinates(only(GT.reference_quadratures(Q)))
        for func in operators
            call = func === GT.value ? GT.call_term(physical_map,coordinates_term) :
                GT.call_term(GT.leaf_term(func;is_compile_constant=true),physical_map,coordinates_term)
            term = GT.optimize(call)
            @test term isa GT.TabulatedTerm
            kernel = compile_term(term,Q,mesh;hoist=true,
                forbidden_accessors=("physical_map_accessor","shape_function_accessor","jacobian_accessor"))
            expected = GT.physical_map_accessor(func,mesh,Val(2))(GT.faces(dom)[1])
            for point in 1:npoints
                @test Base.invokelatest(kernel,1,point,1,1,1) ≈ expected(points[point])
            end
            if func !== GT.value && func !== ForwardDiff.jacobian
                @test GT.prototype(term) ≈ expected(first(points))
            end
        end

        original_operator = scaled_map_operator(2.0)
        call = GT.call_term(GT.leaf_term(original_operator),physical_map,coordinates_term)
        term = GT.optimize(call)
        scaled_map = compile_term(term,Q,mesh,original_operator;hoist=true,
            arguments=(Q,mesh,scaled_map_operator(3.0)))
        @test Base.invokelatest(scaled_map,1,1,1,1,1) ≈ 3sum(abs2,coordinates(1))

        shift = SVector(0.3,-0.1)
        shifted_mesh = GT.replace_node_coordinates(mesh,[x+shift for x in GT.node_coordinates(mesh)])
        term = GT.optimize(GT.call_term(physical_map,coordinates_term))
        shifted_map = compile_term(term,Q,mesh;arguments=(Q,shifted_mesh))
        @test Base.invokelatest(shifted_map,1,1,1,1,1) ≈ coordinates(1)+shift

        for func in (GT.value,ForwardDiff.gradient)
            deps = (GT.leaf_term(func;is_compile_constant=true),GT.leaf_term(V),GT.leaf_term(dom),
                GT.leaf_term(:face),GT.leaf_term(1;is_compile_constant=true),GT.leaf_term(:field),
                GT.leaf_term(:dof),GT.leaf_term(1;is_compile_constant=true),GT.leaf_term(:face_around))
            parent = GT.FormArgumentTerm(deps,2,reference,:interior)
            term = GT.TabulatedTerm(parent,GT.leaf_term(Q),GT.leaf_term(:point))
            shape = compile_term(term,Q,V;
                forbidden_accessors=("form_argument_accessor","shape_function_accessor","jacobian_accessor"))
            expected = GT.shape_function_accessor(func,V,Q)(1)
            for point in 1:npoints, dof in 1:GT.num_dofs(only(GT.reference_spaces(V)))
                actual = Base.invokelatest(shape,1,point,dof,1,1)
                @test actual ≈ expected(point)(dof)
                @test Base.invokelatest(shape,1,point,dof,2,1) == zero(actual)
            end

            # Evaluate the callable argument away from the quadrature points.
            untabulated = compile_term(parent,V,dom;hoist=true,
                forbidden_accessors=("form_argument_accessor","shape_function_accessor"))
            xi = SVector(0.23,0.31)
            mesh_face = GT.faces(dom)[1]
            phi = GT.physical_map_accessor(GT.value,mesh,Val(2))(mesh_face)
            x = reference ? xi : phi(xi)
            J = ForwardDiff.jacobian(phi,xi)
            reference_shapes = GT.shape_functions(only(GT.reference_spaces(V)))
            expected_values = map(reference_shapes) do s
                if func === GT.value
                    s(xi)
                else
                    gradient = ForwardDiff.gradient(s,xi)
                    reference ? gradient : transpose(inv(J))*gradient
                end
            end
            for dof in eachindex(reference_shapes)
                argument = Base.invokelatest(untabulated,1,1,dof,1,1)
                @test Base.invokelatest(argument,x) ≈ expected_values[dof]
                for (field,side) in ((2,1),(1,2))
                    masked = Base.invokelatest(untabulated,1,1,dof,field,side)
                    @test Base.invokelatest(masked,x) == zero(expected_values[dof])
                end
            end
            @test GT.prototype(parent)(x) ≈ first(expected_values)

            call = GT.optimize(GT.call_term(parent,GT.leaf_term(x)))
            @test call isa GT.CallTerm
            @test GT.prototype(call) ≈ first(expected_values)
            direct = compile_term(call,V,dom;hoist=true)
            @test Base.invokelatest(direct,1,1,1,1,1) ≈ first(expected_values)
        end

        # Scalar assembly reuses data expressions derived from the captured quadrature.
        term = GT.lambda_term(weight_term,:face,:point)
        term,data = GT.capture(GT.parametrize(term,Q))
        f = GT.evaluate(GT.statements_expr(GT.expression(term)),data)
        kernel = Base.invokelatest(f,Q)
        @test Base.invokelatest(kernel,1,1) ≈ weights(1)

    end
end

@testset "User matrix assembly and shared integration data" begin
    mesh = GT.cartesian_mesh((0,1,0,1,0,1),(2,2,2);simplexify=false)
    omega = GT.interior(mesh)
    k = 2
    measure = GT.measure(omega,2*k)
    V = GT.lagrange_space(omega,k)
    a = (u,v) -> GT.∫(x -> dot(ForwardDiff.gradient(u,x),ForwardDiff.gradient(v,x)),measure)
    integral = a(GT.form_argument_quantity(V,2),GT.form_argument_quantity(V,1))
    contribution = only(GT.contributions(integral))
    Q = GT.quadrature(contribution)
    A,cache = GT.assemble_matrix(a,Float64,V,V;reuse=Val(true),assembly_method=(;))
    term = GT.optimize(GT.write_assemble_matrix(contribution,V,V;allocation=cache.alloc))
    captured,data = GT.capture(GT.parametrize(term))
    assembly = captured.body.body
    @test GT.prototype(assembly.alloc) === cache.alloc
    @test !any(value -> value === cache.alloc,data)
    expr = GT.expression(assembly.term)
    constructors = Dict(:mesh_integration_data => Expr[], :space_integration_data => Expr[],
        :DOFData => Expr[], :domain_integration_data => Expr[], :face_incidence_data => Expr[])
    callables = GT.JacobianSumCallable[]
    function collect_constructors(node)
        node isa Expr || return
        if node.head === :call && node.args[1] isa Symbol && haskey(constructors,node.args[1])
            push!(constructors[node.args[1]],node)
        elseif node.head === :call && node.args[1] isa GT.JacobianSumCallable
            push!(callables,node.args[1])
        end
        foreach(collect_constructors,node.args)
    end
    collect_constructors(expr)
    @test count(value -> value === V,data) == 1
    @test count(value -> value === Q,data) == 1
    @test GT.expression(assembly.space_trial) == GT.expression(assembly.space_test)
    @test length(unique(constructors[:mesh_integration_data])) == 1
    @test length(unique(constructors[:domain_integration_data])) == 1
    @test isempty(constructors[:face_incidence_data])
    @test length(unique(constructors[:space_integration_data])) == 1
    @test length(unique(callables)) == 1
    @test only(unique(callables)) === GT.JacobianSumCallable{3,3,0.0,8,false}()
    space_arg = only(arg for (arg,value) in zip(captured.args,data) if value === V)
    quadrature_arg = only(arg for (arg,value) in zip(captured.args,data) if value === Q)
    mesh_data_expr = first(constructors[:mesh_integration_data])
    space_data_expr = first(constructors[:space_integration_data])
    @test mesh_data_expr.args[2] == :(GT.mesh(GT.domain($quadrature_arg)))
    @test mesh_data_expr.args[3] == quadrature_arg
    @test space_data_expr.args[3:4] == [space_arg,quadrature_arg]

    empty!(callables)
    assembly_expr = GT.expression(assembly)
    collect_constructors(assembly_expr)
    @test !isempty(callables)
    @test all(callable -> callable === GT.JacobianSumCallable{3,3,0.0,8,false}(),callables)
    @test !occursin("num_points_accessor",string(assembly_expr))
    @test !occursin("max_num_reference_dofs",string(assembly_expr))
    @test !occursin("dofs_accessor",string(assembly_expr))
    @test !occursin("face_incidence_data",string(assembly_expr))
    @test !occursin("eltype(alloc)",string(assembly_expr))
    @test length(unique(constructors[:DOFData])) == 1

    # Assemble the independent reference with the existing accessor kernels.
    expected = zeros(GT.num_free_dofs(V),GT.num_free_dofs(V))
    dofs = GT.dofs_accessor(V,omega)
    gradients = GT.shape_function_accessor(ForwardDiff.gradient,V,Q)
    weights = GT.weight_accessor(Q)
    npoints = GT.num_points_accessor(Q)
    for face in 1:GT.num_faces(omega)
        ids = dofs(face,1)
        point_gradients = gradients(face,1)
        point_weights = weights(face)
        for point in 1:npoints(face)
            dof_gradients = point_gradients(point)
            values = [dof_gradients(i) for i in eachindex(ids)]
            weight = point_weights(point)
            for j in eachindex(ids), i in eachindex(ids)
                expected[ids[i],ids[j]] += dot(values[i],values[j])*weight
            end
        end
    end
    @test size(A) == (125,125)
    @test Matrix(A) ≈ expected rtol=1e-12 atol=1e-13
    @test norm(A-A',Inf) < 1e-11
    @test norm(A*ones(size(A,2)),Inf) < 1e-11
    kernel = Base.invokelatest(only(cache.loops).prepare)
    element_buffers = [getfield(kernel,name) for name in fieldnames(typeof(kernel))
        if getfield(kernel,name) isa Matrix{Float64}]
    @test length(element_buffers) == 1
    GT.update_matrix!(A,cache)
    @test Matrix(A) ≈ expected rtol=1e-12 atol=1e-13
    updated_kernel = Base.invokelatest(only(cache.loops).prepare)
    @test any(name -> getfield(updated_kernel,name) === only(element_buffers),fieldnames(typeof(updated_kernel)))

    A32,cache32 = GT.assemble_matrix(a,Float32,V,V;reuse=Val(true),assembly_method=(;))
    @test eltype(A32) === Float32
    @test Matrix(A32) ≈ expected rtol=2e-6 atol=2e-6
    GT.update_matrix!(A32,cache32)
    @test Matrix(A32) ≈ expected rtol=2e-6 atol=2e-6
end

@testset "Metadata loop bounds for mixed DG fields" begin
    mesh = GT.cartesian_mesh((0,2,0,1),(2,1);boundary=false,simplexify=false)
    omega = GT.interior(mesh)
    measure = GT.measure(omega,4)
    V1 = GT.lagrange_space(omega,1;continuous=Val(false))
    V2 = GT.lagrange_space(omega,2;continuous=Val(false))
    trial = GT.cartesian_product(V1,V2)
    test = trial
    a = ((u1,u2),(v1,v2)) -> GT.∫(x -> (u1(x)+u2(x))*(v1(x)+v2(x)),measure)
    integral = a(GT.form_argument_quantity(trial,2),GT.form_argument_quantity(test,1))
    contribution = only(GT.contributions(integral))
    Q = GT.quadrature(contribution)

    # The bounds must follow each field rather than the largest field.
    expected = zeros(GT.num_free_dofs(test),GT.num_free_dofs(trial))
    trial_shapes = map(V -> GT.shape_function_accessor(GT.value,V,Q),GT.fields(trial))
    test_shapes = map(V -> GT.shape_function_accessor(GT.value,V,Q),GT.fields(test))
    trial_dofs = map(V -> GT.dofs_accessor(V,omega),GT.fields(trial))
    test_dofs = map(V -> GT.dofs_accessor(V,omega),GT.fields(test))
    column_offsets = (0,GT.num_free_dofs(V1))
    row_offsets = (0,GT.num_free_dofs(V1))
    weights = GT.weight_accessor(Q)
    npoints = GT.num_points_accessor(Q)
    for face in 1:GT.num_faces(omega), point in 1:npoints(face)
        weight = weights(face)(point)
        for row_field in 1:2, column_field in 1:2
            rows = test_dofs[row_field](face,1)
            columns = trial_dofs[column_field](face,1)
            row_values = test_shapes[row_field](face,1)(point)
            column_values = trial_shapes[column_field](face,1)(point)
            for j in eachindex(columns), i in eachindex(rows)
                row = row_offsets[row_field]+rows[i]
                column = column_offsets[column_field]+columns[j]
                expected[row,column] += row_values(i)*column_values(j)*weight
            end
        end
    end
    A = GT.assemble_matrix(a,Float64,trial,test;assembly_method=(;))
    @test Matrix(A) ≈ expected rtol=1e-12 atol=1e-13

    # Vector templates use the same metadata bounds without quadrature accessors.
    l = v -> GT.∫(x -> v[1](x)+v[2](x),measure)
    vector_contribution = only(GT.contributions(l(GT.form_argument_quantity(test,1))))
    b,cache = GT.assemble_vector(l,Float32,test;reuse=Val(true),assembly_method=(;))
    vector_term = GT.optimize(GT.write_assemble_vector(vector_contribution,test;allocation=cache.alloc))
    captured,_ = GT.capture(GT.parametrize(vector_term))
    vector_expr = GT.expression(captured)
    @test !occursin("num_points_accessor",string(vector_expr))
    @test !occursin("max_num_reference_dofs",string(vector_expr))
    @test !occursin("dofs_accessor",string(vector_expr))
    @test occursin("DOFData",string(vector_expr))
    @test !occursin("eltype(alloc)",string(vector_expr))
    @test eltype(b) === Float32
    @test sum(b) ≈ 2prod((2,1))
    expected_b = copy(b)
    GT.update_vector!(b,cache)
    @test b ≈ expected_b
end

end

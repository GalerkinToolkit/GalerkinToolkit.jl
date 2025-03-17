abstract type NewAbstractTerm <: AbstractType end

#"""
#This is used to create the argument of functions converting quantities to terms
#"""
#function term_options(domain;domain_face=:domain_face)
#    contents = (;domain,domain_face)
#    TermOptions(contents)
#end
#
#struct TermOptions{A}
#    contents::A
#end
#
#domain(t::TermOptions) = t.contents.domain
#domain_face(t::TermOptions) = t.contents.domain_face


abstract type NewAbstractQuantity <: AbstractType end


# TODO: write test cases
@auto_hash_equals cache=true typearg=true  fields=(value, ) struct LeafTerm{A, B}  <: NewAbstractTerm
    value::A
    prototype::B
end

@auto_hash_equals cache=true fields=(callee, args) struct CallTerm{A, B, C}  <: NewAbstractTerm
    callee::A
    args::B
    prototype::C
end


@auto_hash_equals cache=true fields=(array, index) struct IndexTerm{A, B, C}  <: NewAbstractTerm
    array::A
    index::B
    prototype::C
end

@auto_hash_equals cache=true fields=(body, args) struct LambdaTerm{A, B, C}  <: NewAbstractTerm
    body::A
    args::B
    prototype::C
end

children(t::LeafTerm) = ()
AbstractTrees.children(t::LeafTerm) = (t.value,)

children(t::CallTerm) = (t.callee, t.args...)

children(t::IndexTerm) = (t.array, t.index)

children(t::LambdaTerm) = begin
    if t.args isa NewAbstractTerm
        (t.body, t.args)
    else
        (t.body, t.args...)
    end
end

# TODO: fill prototypes later?
function call_term(callee, args)
    CallTerm(callee, args, nothing)
end

function lambda_term(body, args)
    LambdaTerm(body, args, nothing) 
end

function index_term(array, index)
    IndexTerm(array, index, nothing)
end

AbstractTrees.children(t::NewAbstractTerm) = children(t)


const hashed_terms = Set([LeafTerm, CallTerm, IndexTerm, LambdaTerm])


struct StatementTerm{A, B, C} <: NewAbstractTerm
    lhs::A 
    rhs::B
    prototype::C
end

children(t::StatementTerm) = (t.lhs, t.rhs)

struct BlockTerm{A} <: NewAbstractTerm
    statements::A
end

children(t::BlockTerm) = (t.statements...,)


struct DiscreteFieldTerm{A, B} <: NewAbstractTerm
    discrete_field::A
    prototype::B
    name::Symbol
end
# Maybe we can remove name and use only objectid.
# Be careful with uniform_quantity as two different instances might have the same constant inside.

children(t::DiscreteFieldTerm) = ()

struct TabulatedDiscreteFieldTerm{A, B, C, D, E, F} <: NewAbstractTerm
    linear_operation::A 
    discrete_field::B 
    domain_face::C
    point::D
    measure::E 
    prototype::F
end


children(t::TabulatedDiscreteFieldTerm) = filter(x -> !isnothing(x), (t.linear_operation, t.domain_face, t.point))


prototype(t::BlockTerm) = (length(t.statements) == 0) ? nothing : prototype(t.statements[end])
prototype(t::NewAbstractTerm) = t.prototype


struct NewFormArgumentTerm{A, B, C, D} <: NewAbstractTerm 
    space::A
    domain_face::B 
    dof::C
    prototype::D
    arg::Int
    field::Int
    name::Symbol #TODO: is it needed?
end

children(t::NewFormArgumentTerm) = filter(x -> !isnothing(x), (t.space, t.domain_face, t.dof))

AbstractTrees.children(t::NewFormArgumentTerm) = (t.space, t.domain_face, t.dof, t.arg, t.field)


struct TabulatedFormArgumentTerm{A, B, C, D, E, F} <: NewAbstractTerm 
    linear_operation::A 
    form_argument::NewFormArgumentTerm
    domain_face::B
    point::C
    dof::D
    measure::E
    prototype::F
end


children(t::TabulatedFormArgumentTerm) = filter(x -> !isnothing(x), (t.linear_operation, t.form_argument, t.domain_face, t.point, t.dof))


AbstractTrees.children(t::TabulatedFormArgumentTerm) = (t.linear_operation, t.form_argument.space, t.domain_face, t.point, t.dof, t.measure, t.form_argument.arg, t.form_argument.field)

struct OneFormTerm{A, B, C, D, E, F} <: NewAbstractTerm 
    contribution::A
    domain_face::B
    dof::C
    ndofs::D
    local_vector::E
    prototype::F
end

struct TwoFormTerm{A, B, C, D, E, F, G, H} <: NewAbstractTerm 
    contribution::A
    domain_face::B
    dof_trial::D
    dof_test::C
    ndofs_trial::F
    ndofs_test::E
    local_matrix::G
    prototype::H
end

children(t::OneFormTerm) = filter(x -> !isnothing(x), (t.contribution, t.domain_face, t.dof, t.ndofs, t.local_vector))
children(t::TwoFormTerm) = filter(x -> !isnothing(x), (t.contribution, t.domain_face, t.dof_trial, t.dof_test, t.ndofs_trial, t.ndofs_test, t.local_matrix))

space(t::TabulatedFormArgumentTerm) = space(t.form_argument)
space(t::NewFormArgumentTerm) = t.space 


function build_leaf_term(a, a_sym)
    if a isa NewAbstractTerm
        a 
    else
        LeafTerm(a_sym, a) 
    end
end

macro new_term(expr)
    vars = Set{Symbol}()
    function findvars!(a)
        nothing
    end
    function findvars!(expr::Expr)
        if  expr.head === :(=)
            var = expr.args[1]
            push!(vars,var)
        end
        map(findvars!,expr.args) # always traverse the args. If the user write (a = b = 3) then only a is captured
        nothing
    end

    transform(a::LineNumberNode) = a
    transform(a::Union{Number, Nothing}) = Expr(:call, build_leaf_term, a, a)
    function transform(a::Symbol)
        if a in vars
            a
        else
            Expr(:call, build_leaf_term, a, Expr(:quote, a))
        end
    end

    function transform(expr::Expr)
        
        if  expr.head === :call
            transform_call(expr)
        elseif  expr.head === :ref
            transform_ref(expr)
        elseif  expr.head === :$
            transform_interpolation(expr)
        elseif  expr.head === :(=)
            transform_default(expr)
        elseif  expr.head === :(->)
            transform_lambda(expr)
        # elseif  expr.head === :do
        #     transform_do(expr)
        elseif  expr.head === :block
            transform_default(expr)
        elseif  expr.head === :tuple
            transform_tuple(expr)
        # elseif  expr.head === :.
        #     transform_dot(expr)
        else
            error("Expr with head=$(expr.head) not supported in macro @term")
        end
    end
    function transform_call(expr::Expr)
        args = map(transform,expr.args)
        Expr(:call, call_term, args[1], Expr(:vect, args[2:end]...))
    end
    function transform_lambda(expr::Expr)
        args = map(transform,expr.args)
        Expr(:call, lambda_term, args[2], args[1])
    end
    # function transform_do(expr::Expr)
    #     expr2 = Expr(:call,expr.args[1].args[1],expr.args[2],expr.args[1].args[2])
    #     transform(expr2)
    # end
    function transform_interpolation(expr::Expr)
        expr.args[1]
    end
    function transform_tuple(expr::Expr)
        args = map(transform,expr.args)
        Expr(:tuple,args...)
    end
    # function transform_dot(expr::Expr)
    #     expr
    # end
    function transform_ref(expr::Expr)
        args = map(transform,expr.args)
        quote
            Expr(:call, index_term, args[1], args[2])
        end
    end
    function transform_default(expr::Expr)
        head = expr.head
        args = map(transform,expr.args)
        Expr(head,args...)
    end
    findvars!(expr)
    transform(expr) |> esc
end

function form_argument(space, arg::Int, field::Int) 
    # arg: test or trial
    # field: block of the global array and matrix 
    name = gensym("form-argument-$arg-$field") # TODO: gensym or just a simple Symbol?
    new_quantity(;name) do opts
        NewFormArgumentTerm(space, nothing, nothing, identity, arg, field, name)
    end

end



function lower(t::LeafTerm)
    t.value
end

function lower(t::CallTerm)
    args = map(lower, t.args)
    :($(lower(t.callee))($(args...)))
    # Expr(:call, lower(t.callee), args...)
end


function lower(t::LambdaTerm)
    args = (t.args isa Tuple) ? Expr(:tuple, map(x -> lower(x), t.args)...) : lower(t.args)
    :($args -> ($(lower(t.body))))
    # Expr(:(->), lower(t.args), lower(t.body))
end


function lower(t::IndexTerm)
    array = lower(t.array)
    index = lower(t.index)
    :($array[$index])
    # Expr(:ref, array, index)
end

function lower(t::StatementTerm)
    lhs = lower(t.lhs)
    rhs = lower(t.rhs)
    :($lhs = $rhs)
    # Expr(:(=), lhs, rhs)
end


function lower(t::BlockTerm)
    exprs = map(lower, t.statements)
    Expr(:block, exprs...)
end


function new_quantity(term;name=nothing)
    NewQuantity(term,name)
end

function compile_constant_quantity(v)
    new_quantity(;) do opts 
        LeafTerm(v, v)
    end
end


struct NewDiscreteField{A,B,C,D} <: GT.AbstractField
    mesh::A
    space::B
    free_values::C
    dirichlet_values::D
    name::Symbol
end

free_values(u::NewDiscreteField) = u.free_values
dirichlet_values(u::NewDiscreteField) = u.dirichlet_values
mesh(u::NewDiscreteField) = u.mesh
space(u::NewDiscreteField) = u.space
name(u::NewDiscreteField) = u.name



# TODO remove
function discrete_field_quantity(m::NewDiscreteField)
    name = GT.name(m)
    new_quantity(;name) do opts
        free_vals = GT.free_values(m)
        diri_vals = GT.dirichlet_values(m)
        proto = (length(free_vals) > 0) ? free_vals[1] : diri_vals[1]
        DiscreteFieldTerm(m, x -> proto, name)
    end
end

struct NewQuantity{A,B} <: NewAbstractQuantity
    term::A
    name::B
end


function (f::NewAbstractQuantity)(x::NewAbstractQuantity...)
    new_quantity(;) do opts
        callee = term(f, opts)
        args = map(x) do arg
            term(arg, opts)
        end
        prototype_callee = prototype(callee)
        prototype_args = map(prototype, args)
        CallTerm(callee, args, prototype_callee(prototype_args...))
    end
end

function Base.getindex(a::NewAbstractQuantity, b) # TODO: higher-dimensional indexing
    if !(b isa NewAbstractQuantity)
        index = compile_constant_quantity(b)
    else
        index = b
    end

    new_quantity(;) do opts
        array_term = term(a, opts)
        index_term = term(index, opts)
        proto = prototype(array_term)[1]
        IndexTerm(array_term, index_term, proto)
    end
end


function name(m::NewQuantity)
    @assert m.name !== nothing
    m.name
end

function term(q::NewQuantity,dom)
    q.term(dom)
end

function uniform_quantity(value;name=gensym("uniform"))
    new_quantity(;name) do opts
        UniformTerm(name, value, opts)
    end
end

struct UniformTerm{A, B, C} <: NewAbstractTerm
    name::A
    value::B
    domain::C
end

children(t::UniformTerm) = ()
AbstractTrees.children(t::UniformTerm) = (t.value,)

domain(t::UniformTerm) = t.domain
prototype(t::UniformTerm) = t.value


function value(t::UniformTerm)
    t.value
end

function value(t::DiscreteFieldTerm)
    t.discrete_field
end

struct NumPointsTerm{A, B} <: NewAbstractTerm
    measure::A
    domain_face::B
end

children(t::NumPointsTerm) = filter(x -> !isnothing(x), (t.domain_face))

prototype(a::NumPointsTerm) = 0
domain(t::NumPointsTerm) = domain(t.measure)


struct WeightTerm{A,B,C} <: NewAbstractTerm
    measure::A
    point::B
    domain_face::C
end

children(t::WeightTerm) = filter(x -> !isnothing(x), (t.point, t.domain_face))


domain(t::WeightTerm) = domain(t.measure)
prototype(a::WeightTerm) = 0.0


struct CoordinateTerm{A,B,C} <: NewAbstractTerm
    measure::A
    domain_face::B
    point::C
end

children(t::CoordinateTerm) = filter(x -> !isnothing(x), (t.point, t.domain_face))


domain(t::CoordinateTerm) = domain(t.measure)

prototype(t::CoordinateTerm) = begin
    dom = domain(t.measure)
    mesh = GT.mesh(t.measure.quadrature)
    node_to_x = node_coordinates(mesh)
    node_to_x[1]
end


struct ContributionTerm{A,B,C,D,E} <: NewAbstractTerm
    integrand::A
    weight::B
    point::C # gets reduced
    num_points::D
    prototype::E
end

children(t::ContributionTerm) = filter(x -> !isnothing(x), (t.integrand, t.weight, t.point, t.num_points))

prototype(t::ContributionTerm) = t.prototype

function rewrite_l1(t::OneFormTerm)
    OneFormTerm(rewrite_l1(t.contribution), t.domain_face, t.dof, t.ndofs, t.local_vector, t.prototype)
end

function rewrite_l1(t::TwoFormTerm)
    TwoFormTerm(rewrite_l1(t.contribution), t.domain_face, t.dof_trial, t.dof_test, t.ndofs_trial, t.ndofs_test, t.local_matrix, t.prototype)
end


function rewrite_l1(t::Union{UniformTerm, LeafTerm, NumPointsTerm, WeightTerm, CoordinateTerm, NewFormArgumentTerm})
    t
end

function rewrite_l1(t::IndexTerm)
    IndexTerm(rewrite_l1(t.array), rewrite_l1(t.index), t.prototype)
end

function rewrite_l1(t::LambdaTerm)
    LambdaTerm(rewrite_l1(t.body), t.args, t.prototype)
end


function rewrite_l1(t::CallTerm)
    if (t.callee isa DiscreteFieldTerm)
        if length(t.args) == 1  &&  t.args[1] isa CoordinateTerm
            # TODO: other operators
            TabulatedDiscreteFieldTerm(LeafTerm(:value, value), t.callee, t.args[1].domain_face, t.args[1].point, 
                                    t.args[1].measure, t.callee.prototype) 
        else
            error("arg type not supported for DiscreteField calls!")
        end
    elseif t.callee isa NewFormArgumentTerm
        if length(t.args) == 1  &&  t.args[1] isa CoordinateTerm
            # TODO: other operators
            TabulatedFormArgumentTerm(LeafTerm(:value, value), t.callee, t.args[1].domain_face, t.args[1].point, t.callee.dof,
                                    t.args[1].measure, t.callee.prototype)
        else
            error("arg type not supported for DiscreteField calls!")
        end
    else
        CallTerm(rewrite_l1(t.callee), map(rewrite_l1, t.args), t.prototype)
    end
end


function rewrite_l1(t::ContributionTerm)
    ContributionTerm(rewrite_l1(t.integrand), rewrite_l1(t.weight), t.point, rewrite_l1(t.num_points), t.prototype)
end


function lower_l1_to_l2(t::UniformTerm, name_to_symbol)
    t2 = name_to_symbol[t.name]
    @new_term begin
        value($(LeafTerm(t2, t)))
    end
    # CallTerm(LeafTerm(:value, value), [LeafTerm(t2, t)], prototype(t))
    # :(value($t2))
end


function lower_l1_to_l2(t::LambdaTerm, name_to_symbol)
    body = lower_l1_to_l2(t.body, name_to_symbol)
    @new_term begin
        $(t.args) -> $body
    end
    # LambdaTerm(lower_l1_to_l2(t.body, name_to_symbol), t.args, prototype(t))
end

function lower_l1_to_l2(t::IndexTerm, name_to_symbol)
    IndexTerm(lower_l1_to_l2(t.array, name_to_symbol), lower_l1_to_l2(t.index, name_to_symbol), t.prototype)
end


function lower_l1_to_l2(t::CallTerm, name_to_symbol)
    callee = lower_l1_to_l2(t.callee, name_to_symbol)
    args = map(t.args) do arg 
        lower_l1_to_l2(arg, name_to_symbol)
    end
    CallTerm(callee, args, prototype(t))
end


function lower_l1_to_l2(t::NumPointsTerm, name_to_symbol)
    domain_face_term = t.domain_face
    measure = name_to_symbol[t.measure.name]
    
    
    # face_points = CallTerm(LeafTerm(:num_points_accessor, num_points_accessor), [LeafTerm(measure, t.measure)], x -> prototype(t))
    # CallTerm(face_points, [domain_face_term], prototype(t))
    measure_term = LeafTerm(measure, t.measure)
    @new_term begin
        face_npoints = num_points_accessor($measure_term)
        face_npoints($domain_face_term)
    end
end


function lower_l1_to_l2(t::WeightTerm, name_to_symbol)
    measure = name_to_symbol[t.measure.name]
    point_term = t.point
    domain_face_term = t.domain_face
    measure_term = LeafTerm(measure, t.measure)

    m = GT.mesh(t.measure.quadrature)
    x = node_coordinates(m)[1]
    jacobian_prototype = x * transpose(x)


    # weight_func_0 = CallTerm(LeafTerm(:weight_accessor, weight_accessor), [measure_term], x -> ( (y1, y2) -> 0.0))
    # weight_func = CallTerm(weight_func_0, [domain_face_term], ((y1, y2) -> 0.0))

    # jacobian_func_0 = CallTerm(LeafTerm(:jacobian_accessor, jacobian_accessor), [measure_term], (x -> y -> jacobian_prototype))
    # jacobian_func = CallTerm(jacobian_func_0, [domain_face_term], (y -> jacobian_prototype))
    # jacobian_result = CallTerm(jacobian_func, [point_term], jacobian_prototype)

    # CallTerm(weight_func, [point_term, jacobian_result], prototype(t))

    @new_term  begin
        weight_func = weight_accessor($measure_term)($domain_face_term)
        jacobian_func = jacobian_accessor($measure_term)($domain_face_term)
        J = jacobian_func($point_term)
        weight_func($point_term, J)
    end
end

function new_measure(domain,degree;name=gensym("quadrature"))
    q = quadrature(domain,degree)
    NewMeasure(q,name)
end

struct NewMeasure{A,B}
    quadrature::A
    name::B
end

domain(m::NewMeasure) = domain(m.quadrature)
name(m::NewMeasure) = m.name
mesh(m::NewMeasure) = mesh(m.quadrature)
weight_accessor(m::NewMeasure) = weight_accessor(m.quadrature)
jacobian_accessor(m::NewMeasure) = jacobian_accessor(m.quadrature)
coordinate_accessor(m::NewMeasure) = coordinate_accessor(m.quadrature)
num_points_accessor(m::NewMeasure) = num_points_accessor(m.quadrature)
discrete_field_accessor(f, uh, m::NewMeasure) = begin
    
    discrete_field_accessor(f, uh, m.quadrature)
end

function coordinate_quantity(measure::NewMeasure, point)
    new_quantity(;) do dom
        @assert dom == domain(measure)
        CoordinateTerm(measure, nothing, LeafTerm(point, 1))
    end
end


function lower_l1_to_l2(t::CoordinateTerm, name_to_symbol)
    # coordinate_accessor($measure)($domain_face)($point)

    measure = name_to_symbol[t.measure.name]

    measure_term = LeafTerm(measure, t.measure)
    domain_face_term = t.domain_face
    point_term = t.point

    m = GT.mesh(t.measure.quadrature)
    prototype_point = node_coordinates(m)[1]

    @new_term begin
        point_func = coordinate_accessor($measure_term)($domain_face_term)
        point_func($point_term)
    end

    # point_func_0 = CallTerm(LeafTerm(:coordinate_accessor, coordinate_accessor), [measure_term], x -> y -> prototype_point)
    # point_func = CallTerm(point_func_0, [domain_face_term], y -> prototype_point)
    # CallTerm(point_func, [point_term], prototype_point)
end


function contribution(f, measure::NewMeasure)
    point = :point
    x = coordinate_quantity(measure,point)
    fx = f(x)
    weight = WeightTerm(measure, LeafTerm(point, 1), nothing)
    num_points = NumPointsTerm(measure, nothing)
    domain = GT.domain(measure)
    integrand = term(fx,domain)
    ContributionTerm(
                     integrand,
                     weight,
                     point,
                     num_points, 
                     prototype(integrand))
end

const new_∫ = contribution


function lower_l1_to_l2(t::ContributionTerm, name_to_symbol)
    integrand_expr = lower_l1_to_l2(t.integrand,name_to_symbol)
    weight_expr = lower_l1_to_l2(t.weight,name_to_symbol)
    num_points_expr = lower_l1_to_l2(t.num_points,name_to_symbol)
    point = t.point
    point_term = LeafTerm(point, 1)

    # l_body = CallTerm(LeafTerm(:*, *), [integrand_expr, weight_expr], prototype(integrand_expr))
    # l = LambdaTerm(l_body, LeafTerm(point, 1), x -> prototype(integrand_expr))
    # range = CallTerm(LeafTerm(:(:), :), [LeafTerm(1, 1), num_points_expr], 1:1)
    # CallTerm(LeafTerm(:(sum), sum), [l, range], prototype(t))
    @new_term begin
        sum($point_term -> $integrand_expr * $weight_expr,  1:$num_points_expr)
    end
end



function capture!(a::ContributionTerm, name_to_symbol, name_to_captured_data)
    capture!(a.integrand, name_to_symbol, name_to_captured_data)
    capture!(a.weight, name_to_symbol, name_to_captured_data)
    capture!(a.num_points, name_to_symbol, name_to_captured_data)
end
function capture!(a::Union{CoordinateTerm, WeightTerm, NumPointsTerm}, name_to_symbol, name_to_captured_data)
    name = a.measure.name
    if !haskey(name_to_symbol, name)
        name_to_captured_data[name] = a.measure
    end
end
function capture!(a::UniformTerm, name_to_symbol, name_to_captured_data)
    name = a.name
    if !haskey(name_to_symbol, name)
        name_to_captured_data[name] = a
    end
end

function capture!(a::CallTerm, name_to_symbol, name_to_captured_data)
    capture!(a.callee, name_to_symbol, name_to_captured_data)
    map(a.args) do arg
        capture!(arg, name_to_symbol, name_to_captured_data)
    end
end


function capture!(a::IndexTerm, name_to_symbol, name_to_captured_data)
    capture!(a.array, name_to_symbol, name_to_captured_data)
    capture!(a.index, name_to_symbol, name_to_captured_data)
end

function capture!(a::DiscreteFieldTerm, name_to_symbol, name_to_captured_data)
    name = a.name
    if !haskey(name_to_symbol, name)
        name_to_captured_data[name] = a.discrete_field
    end
end

function capture!(a::LeafTerm, name_to_symbol, name_to_captured_data)
end

function capture!(a::OneFormTerm, name_to_symbol, name_to_captured_data)
    capture!(a.contribution, name_to_symbol, name_to_captured_data)
end

function capture!(a::NewFormArgumentTerm, name_to_symbol, name_to_captured_data)
    name = a.name 
    if !haskey(name_to_symbol, name)
        name_to_captured_data[name] = a
    end
end


function lower_l1_to_l2(t::LeafTerm,name_to_symbol)
    t
end

function lower_l1_to_l2(t::OneFormTerm, name_to_symbol)
    ndofs = lower_l1_to_l2(t.ndofs, name_to_symbol)
    body = lower_l1_to_l2(t.contribution, name_to_symbol)
    result = @new_term begin 
        ($(t.local_vector), $(t.domain_face)) -> begin 
            map!($(t.dof) -> $body, $(t.local_vector), 1:$ndofs)
        end
    end
    
    # dof_lambda_term = LambdaTerm(body, t.dof, x -> prototype(t.body))
    # range = CallTerm(LeafTerm(:(:), :), [LeafTerm(1, 1), lower_l1_to_l2(t.ndofs, name_to_symbol)], 1:1)
    # map!_term = CallTerm(LeafTerm(:map!, map!), [dof_lambda_term, t.local_vector, range], nothing)
    # LambdaTerm(map!_term, (t.local_vector, t.domain_face), x -> nothing)
end


function lower_l1_to_l2(t::TwoFormTerm, name_to_symbol)
    ndofs_trial = lower_l1_to_l2(t.ndofs_trial, name_to_symbol)
    ndofs_test = lower_l1_to_l2(t.ndofs_test, name_to_symbol)
    body = lower_l1_to_l2(t.contribution, name_to_symbol)

    map(x -> map!(y -> x + y, view(a, x, :), 1:5), 1:4)
    # TODO: cartesian indexing of a matrix?
    result = @new_term begin 
        ($(t.local_matrix), $(t.domain_face)) -> begin 
            # two_form_lambda = ($(t.dof_trial), $(t.dof_test)) -> $body # ugly to wrap the cartesian index with another lambda, and maybe incorrect
            # map!(cartesian_arg -> two_form_lambda(cartesian_arg[1], cartesian_arg[2]), $(t.local_matrix), CartesianIndices(($ndofs_trial, $ndofs_test)))
            map($(t.dof_trial) -> map!($(t.dof_test) -> $body, view($(t.local_matrix), :, $(t.dof_trial)), 1:$ndofs_test), 1:$ndofs_trial)
            # do we need a `nothing` here? we do not need to collect the lambda results
        end
    end
end


function lower_l1_to_l2(t::NewFormArgumentTerm, name_to_symbol)
    name = name_to_symbol[t.name]
    LeafTerm(name, t)
end

function lower_l1_to_l2(t::TabulatedFormArgumentTerm,name_to_symbol)
    # NewFormArgumentTerm not needed?
    measure = name_to_symbol[t.measure.name]

    dof_term = t.dof
    point_term = t.point 
    domain_face_term = t.domain_face
    measure_term = LeafTerm(measure, t.measure)
    form_argument_term = lower_l1_to_l2(t.form_argument, name_to_symbol)

    @new_term begin 
        s = space($form_argument_term)
        shape_func = shape_function_accessor($(t.linear_operation), s, $measure_term)($domain_face_term)
        shape_func($point_term, nothing)($dof_term)
    end
    # space_term = CallTerm(LeafTerm(:space, space), [lower_l1_to_l2(t.form_argument, name_to_symbol)], t.form_argument.space) 

    # shape_func_0 = CallTerm(LeafTerm(:shape_function_accessor, shape_function_accessor), [t.linear_operation, space_term, measure_term], (x -> y -> z -> t.prototype))
    # shape_func_1 = CallTerm(shape_func_0, [domain_face_term], y -> z -> t.prototype)
    # shape_func_2 = CallTerm(shape_func_1, [point_term, LeafTerm(:nothing, nothing)], z -> t.prototype)
    # CallTerm(shape_func_2, [dof_term], t.prototype)
end


function lower_l1_to_l2(t::DiscreteFieldTerm,name_to_symbol)
    t_symbol = name_to_symbol[t.name]
    LeafTerm(t_symbol, prototype(t))
    # CallTerm(LeafTerm(:value, value), [LeafTerm(t_symbol, t)], prototype(t))
end


function lower_l1_to_l2(t::TabulatedDiscreteFieldTerm,name_to_symbol)
    measure = name_to_symbol[t.measure.name]
    measure_term = LeafTerm(measure, t.measure)
    domain_face_term = t.domain_face
    point_term = t.point

    m = GT.mesh(t.measure.quadrature)
    x = node_coordinates(m)[1]
    jacobian_prototype = x * transpose(x)
    discrete_field = lower_l1_to_l2(t.discrete_field, name_to_symbol)


    @new_term begin 
        jacobian_func = jacobian_accessor($measure_term)($domain_face_term)
        J = jacobian_func($point_term)
        field_func = discrete_field_accessor($(t.linear_operation), $discrete_field, $measure_term)($domain_face_term)
        field_func($point_term, J)
    end
    # jacobian_func_0 = CallTerm(LeafTerm(:jacobian_accessor, jacobian_accessor), [measure_term], (x -> y -> jacobian_prototype))
    # jacobian_func = CallTerm(jacobian_func_0, [domain_face_term], (y -> jacobian_prototype))
    # jacobian_result = CallTerm(jacobian_func, [point_term], jacobian_prototype)

    # prototype_result = prototype(t)
    
    # field_func_0 = CallTerm(LeafTerm(:discrete_field_accessor, discrete_field_accessor), [t.linear_operation, discrete_field, measure_term], x -> (y1, y2) -> prototype_result) 
    # field_func =  CallTerm(field_func_0, [domain_face_term], (y1, y2) -> prototype_result)
    # CallTerm(field_func, [point_term, jacobian_result], prototype_result)
end


function generate(t::NewAbstractTerm,params...)
    itr = [ name(p)=>Symbol("arg$i") for (i,p) in enumerate(params)  ]
    symbols = map(last,itr)
    name_to_symbol = Dict(itr)
    
    # find all leaf terms 
    name_to_captured_data = Dict()
    capture!(t, name_to_symbol, name_to_captured_data)

    # update name_to_symbol and generate captured data symbols
    captured_data = []
    captured_data_symbols = []
    for (i, (k, v)) in enumerate(name_to_captured_data)
        sym = Symbol("captured_arg_$i")
        name_to_symbol[k] = sym 
        push!(captured_data, v)
        push!(captured_data_symbols, sym)
    end

    # t: L1 (domain specific level) term
    domain_face = :dummy_domain_face
    domain_face_term = LeafTerm(domain_face, 1)
    term_l1_1 = set_free_args(t, domain_face_term) 
    term_l1_2 = LambdaTerm(term_l1_1, domain_face_term, x -> prototype(term_l1_1)) # a wrapper of domain face

    # generate body and optimize
    
    # rewrite L1 term
    term_l1_3 = rewrite_l1(term_l1_2)
    # L2: functional level
    term_l2 = lower_l1_to_l2(term_l1_3, name_to_symbol)
    
    # L3: imperative level
    term_l3 = statements(term_l2)
    
    # L4: imperative level with loops
    
    # L5: julia expression
    expr = lower(term_l3)
    fun_block = expr

    result = quote
        ($(captured_data_symbols...), ) -> begin
            ($(symbols...), ) -> begin
                $fun_block
            end
        end
    end
    (MacroTools.striplines(result), captured_data)
end


# TODO: find a general solution for dofs setting (for accessor terms)
function set_free_args(t::NewAbstractTerm, domain_face, dofs=())
    new_term_args = map(propertynames(t)) do child_name
        child = getproperty(t, child_name)
        if child_name === :domain_face
            # Very important to check that this is indeed a free child
            # out criterion is nothing means free
            if child === nothing
                domain_face
            else 
                child
            end
        elseif child isa NewAbstractTerm
            set_free_args(child, domain_face, dofs)
        elseif applicable(iterate, child)
            map(x -> (x isa NewAbstractTerm) ? set_free_args(x, domain_face, dofs) : x, child)
        else
            child
        end
    end
    new_term_type = typeof(t).name.wrapper
    if new_term_type in hashed_terms
        (new_term_type)(new_term_args[1:end-1]...)
    else 
        (new_term_type)(new_term_args...)
    end
end


function set_free_args(t::OneFormTerm, domain_face, dofs=())
    (contribution_term, nodfs_term, local_vector_term) = map(x -> set_free_args(x, domain_face, dofs), (t.contribution, t.ndofs, t.local_vector))
    dof_term = (t.dof === nothing && length(dofs) >= 1) ? dofs[1] : t.dof
    domain_face_term = (t.domain_face === nothing) ? domain_face : t.domain_face
    OneFormTerm(contribution_term, domain_face_term, dof_term, nodfs_term, local_vector_term, t.prototype)
end 

function set_free_args(t::TwoFormTerm, domain_face, dofs=())
    (contribution_term, ndofs_trial_term, ndofs_test_term, local_matrix_term) = map(x -> set_free_args(x, domain_face, dofs), (t.contribution, t.ndofs_trial, t.ndofs_test, t.local_matrix))
    dof_trial_term = (t.dof_trial === nothing && length(dofs) >= 1) ? dofs[1] : t.dof_trial
    dof_test_term = (t.dof_test === nothing && length(dofs) >= 2) ? dofs[2] : t.dof_test
    domain_face_term = (t.domain_face === nothing) ? domain_face : t.domain_face
    TwoFormTerm(contribution_term, domain_face_term, dof_trial_term, dof_test_term, ndofs_trial_term, ndofs_test_term, local_matrix_term, t.prototype)
end 

function set_free_args(t::NewFormArgumentTerm, domain_face, dofs=())
    dof_term = (t.dof === nothing && length(dofs) >= t.arg) ? dofs[t.arg] : t.dof
    domain_face_term = (t.domain_face === nothing) ? domain_face : t.domain_face
    NewFormArgumentTerm(t.space, domain_face_term, dof_term, t.prototype, t.arg, t.field, t.name)
end


function term_ndofs(t::NewFormArgumentTerm, arg::Int, field::Int)
    if (t.arg == arg && t.field == field) 
        # s = CallTerm(LeafTerm(:space, space), [t], t.space)
        # dom = CallTerm(LeafTerm(:domain, domain), [s], domain(t.space))
        # face_dofs = CallTerm(LeafTerm(:dofs_accessor, dofs_accessor), [s, dom], x -> [1]) 
        # dofs = CallTerm(face_dofs, [t.domain_face], [1])
        # CallTerm(LeafTerm(:length, length), [dofs], 1)
        @new_term begin
            s = space($t)
            dom = domain(s)
            ndofs = num_dofs_accessor(s, dom)($(t.domain_face))
        end
    else
        nothing
    end
end


function term_ndofs(t::NewAbstractTerm, arg::Int, field::Int)
    ndofs = nothing
    for child in children(t)
        ndofs = term_ndofs(child, arg, field)
        if ndofs !== nothing
            break
        end
    end
    ndofs
end

# is it sufficient to check the +/- operators only?
# returns the new term and the field of it
function select_field_impl(field::Vector{Int64}, t::NewAbstractTerm)
    (t, zero(field))
end

function select_field_impl(field::Vector{Int64}, t::NewFormArgumentTerm)
    result_field = zero(field)
    result_field[t.arg] = t.field
    (t, result_field)
end


function isa_subset(field, subfield)
    @assert length(field) == length(subfield)
    for (a, b) in zip(field, subfield)
        if (b != 0) && (a != b)
            return false
        end
    end
    true
end

# TODO: hard to handle input errors, and requires to be (completely) tested
function select_field_impl(field::Vector{Int64}, t::CallTerm)
    (callee, call_field) = select_field_impl(field, t.callee)
    args_with_fields = map(t.args) do arg
        select_field_impl(field, arg)
    end

    fields = (call_field, map(last, args_with_fields)...)

    fields_per_arg = zero(field)
    # -1 represents multipe fields, 0 represents no field, positive int represents the existing field
    for fs in fields 
        for (id, f) in enumerate(fs)
            if fields_per_arg[id] == 0
                fields_per_arg[id] = f 
            elseif fields_per_arg[id] != f && f != 0
                fields_per_arg[id] = -1
            end
        end
    end

    # 1 field: return this field
    if all(x -> x >= 0, fields_per_arg)
        args = map(first, args_with_fields)
        return (CallTerm(callee, args, t.prototype), fields_per_arg)
    end
    
    # TODO: can we process other types?
    # 2+ fields: select the given field, and assume this is a + function
    for id in 1:length(fields_per_arg)
        if fields_per_arg[id] < 0
            fields_per_arg[id] = field[id]
        end
    end
    if callee isa LeafTerm 
        if callee.prototype == +
            args_with_fields_filtered = filter(x-> isa_subset(field, last(x)), args_with_fields)
            args = map(first, args_with_fields_filtered)
            if length(args_with_fields_filtered) == 0
                return (LeafTerm(t.prototype, t.prototype), last(last(args_with_fields))) # To be reduced later
            elseif length(args_with_fields_filtered) == 1
                return args_with_fields_filtered[1] 
            else 
                return (CallTerm(callee, args, t.prototype), fields_per_arg)
            end
        elseif callee.prototype == -
            args = map(first, args_with_fields)
            @assert length(args_with_fields) == 2
            if isa_subset(field, args_with_fields[1][2])
                return (args[1], args_with_fields[1][2])
            elseif isa_subset(field, args_with_fields[2][2])
                return (CallTerm(callee, [args[2]], t.prototype), args_with_fields[2][2])
            else
                return (LeafTerm(t.prototype, t.prototype), last(last(args_with_fields)))  # To be reduced later  
            end
        end
    end
    error("customized linear function to combine different fields not supported!")
    
end


function select_field_impl(field::Vector{Int64}, t::ContributionTerm)
    (int_term, int_field) = select_field_impl(field, t.integrand)
    result = ContributionTerm(int_term, t.weight, t.point, t.num_points, t.prototype)
    (result, int_field)
end


function select_field_impl(field::Vector{Int64}, t::IndexTerm)
    (array_term, array_field) = select_field_impl(field, t.array)
    result = IndexTerm(array_term, t.index, t.prototype) # constant index
    (result, array_field)
end

function select_field_impl(field::Vector{Int64}, t::OneFormTerm)
    (contribution, contribution_field) = select_field_impl(field, t.contribution)
    # TODO: simplify oft without form arguments?
    @assert sum(contribution_field) > 0
    result = OneFormTerm(contribution, t.domain_face, t.dof, t.ndofs, t.local_vector, t.prototype)
    (result, contribution_field)
end


function select_field(field, t::Union{OneFormTerm, TwoFormTerm}) # OneFormTerm or TwoFormTerm
    field_vector = [field...]
    (result, f) = select_field_impl(field_vector, t)
    @assert f == field_vector # TODO: what if f != field?
    result
end

function generate_1_form(field::Int, t::NewAbstractTerm, params...)
    itr = [ name(p)=>Symbol("arg$i") for (i,p) in enumerate(params)  ]
    symbols = map(last,itr)
    name_to_symbol = Dict(itr)
    
    # find all leaf terms 
    name_to_captured_data = Dict()
    capture!(t, name_to_symbol, name_to_captured_data)

    # update name_to_symbol and generate captured data symbols
    captured_data = []
    captured_data_symbols = []
    for (i, (k, v)) in enumerate(name_to_captured_data)
        sym = Symbol("captured_arg_$i")
        name_to_symbol[k] = sym 
        push!(captured_data, v)
        push!(captured_data_symbols, sym)
    end

    dof = :dof 
    dof_term = LeafTerm(dof, 1)
    domain_face = :dummy_domain_face
    domain_face_term = LeafTerm(domain_face, 1)
    local_vector = :b 
    local_vector_term = LeafTerm(local_vector, [prototype(t)])
    
    term_l1_1 = set_free_args(t, domain_face_term, (dof_term, ))
    ndofs_term = term_ndofs(term_l1_1, 1, field)

    term_l1_2 = OneFormTerm(term_l1_1, domain_face_term, dof_term, ndofs_term, local_vector_term, [prototype(term_l1_1)]) # a wrapper of domain face

    # generate body and optimize
    term_l1_3 = select_field((field, ), term_l1_2)
    # rewrite L1 term
    term_l1_4 = rewrite_l1(term_l1_3)
    
    # L2: functional level
    term_l2 = lower_l1_to_l2(term_l1_4, name_to_symbol)
    
    # L3: imperative level
    term_l3 = statements(term_l2)
    
    # L4: imperative level with loops
    
    # L5: julia expression
    expr = lower(term_l3)
    fun_block = expr

    result = quote
        ($(captured_data_symbols...), ) -> begin
            ($(symbols...), ) -> begin
                $fun_block
            end
        end
    end
    (MacroTools.striplines(result), captured_data)
end

# TODO: to be completed and tested
function generate_2_form(field_trial::Int, field_test::Int, t::NewAbstractTerm, params...)
    itr = [ name(p)=>Symbol("arg$i") for (i,p) in enumerate(params)  ]
    symbols = map(last,itr)
    name_to_symbol = Dict(itr)
    
    # find all leaf terms 
    name_to_captured_data = Dict()
    capture!(t, name_to_symbol, name_to_captured_data)

    # update name_to_symbol and generate captured data symbols
    captured_data = []
    captured_data_symbols = []
    for (i, (k, v)) in enumerate(name_to_captured_data)
        sym = Symbol("captured_arg_$i")
        name_to_symbol[k] = sym 
        push!(captured_data, v)
        push!(captured_data_symbols, sym)
    end

    dof_trial = :dof_trial 
    dof_test = :dof_test
    dof_trial_term = LeafTerm(dof_trial, 1)
    dof_test_term = LeafTerm(dof_test, 1)
    domain_face = :dummy_domain_face
    domain_face_term = LeafTerm(domain_face, 1)
    local_matrix = :a
    local_matrix_term = LeafTerm(local_matrix, [[prototype(t)]])
    
    term_l1_1 = set_free_args(t, domain_face_term, (dof_trial_term, dof_test_term))
    ndofs_trial_term = term_ndofs(term_l1_1, 1, field_trial)
    ndofs_test_term = term_ndofs(term_l1_1, 2, field_test)

    term_l1_2 = TwoFormTerm(term_l1_1, domain_face_term, dof_trial_term, dof_test_term, ndofs_trial_term, ndofs_test_term, local_matrix_term, [[prototype(term_l1_1)]]) # a wrapper of domain face

    # generate body and optimize
    term_l1_3 = select_field((field_trial, field_test), term_l1_2)
    # rewrite L1 term
    term_l1_4 = rewrite_l1(term_l1_3)
    
    # L2: functional level
    term_l2 = lower_l1_to_l2(term_l1_4, name_to_symbol)
    
    # L3: imperative level
    term_l3 = statements(term_l2)
    
    # L4: imperative level with loops
    
    # L5: julia expression
    expr = lower(term_l3)
    fun_block = expr

    result = quote
        ($(captured_data_symbols...), ) -> begin
            ($(symbols...), ) -> begin
                $fun_block
            end
        end
    end
    (MacroTools.striplines(result), captured_data)
end


function evaluate(expr_and_captured)
    expr, captured_data = expr_and_captured
    f1 = eval(expr)
    f2 = Base.invokelatest(f1, captured_data...)

    # convert input from quantity to term
    (args..., ) -> begin
        args = map(args) do p
            if p isa NewMeasure 
                p
            elseif p isa NewDiscreteField
                p
            elseif p isa NewQuantity
                term(p, nothing)
            else
                error("param $p is not a uniform quantity or a measure")
            end
        end
        f2(args...)
    end
end

function lambda_args_once!(node::LeafTerm, args)
    return true
end

function lambda_args_once!(node::CallTerm, args)
    all(x -> lambda_args_once!(x, args), (node.callee, node.args...))
end

function lambda_args_once!(node::IndexTerm, args)
    all(x -> lambda_args_once!(x, args), (node.array, node.index))
end

function lambda_args_once!(node::LambdaTerm, args)
    args_lambda = (node.args isa Tuple) ?  map(x -> x.value, node.args) : (node.args.value,)
    for arg in args_lambda
        if arg in args 
            return false
        end
        push!(args, arg)
    end
    lambda_args_once!(node.body, args)
end

function lambda_args_once!(node::OneFormTerm, args)
    return lambda_args_once!(node.contribution, args)
end

function lambda_args_once(node)
    symbols = Set{Symbol}()
    return lambda_args_once!(node, symbols)
end


function statements(node)
    @assert lambda_args_once(node) # keep hash_scope but make an error check. In some cases it will be a bug if we have a lambda function arg name more than once in the term
    root = :root
    scope_level = Dict{Symbol,Int}()
    scope_level[root] = 0
    hash_scope = Dict{UInt,Symbol}()
    scope_rank = Dict{Symbol,Int}()
    scope_rank[root] = 0
    hash_expr = Dict{UInt, NewAbstractTerm}()
    scope_block = Dict(root=>BlockTerm([]))
    function visit(node, depth = 0) # if the terms here are immutable (no terms "appended" after construction) then we can assume there is no graph cycle
        hash = Base.hash(node)
        if haskey(hash_scope, hash)
            return [hash_scope[hash]]
        end
        if node isa LeafTerm
            scopes = [root]
            hash_expr[hash] = node
            hash_scope[hash] = root 
            return scopes
        elseif (node isa CallTerm) || (node isa IndexTerm)
            args = (node isa CallTerm) ? (node.callee, node.args...) : (node.array, node.index)
            scopes_nested = map(args) do arg
                visit(arg, depth)
            end
            scopes = reduce(vcat,scopes_nested) |> unique
            scope = argmax(x->scope_level[x],scopes)
            hash_scope[hash] = scope
            rank = 1 + scope_rank[scope]
            scope_rank[scope] = rank

            var = Symbol("var_$(scope)_$(rank)")
            var_term = LeafTerm(var, prototype(node))
            hash_expr[hash] = var_term 
            args_var = map(args) do arg
                arg_hash = Base.hash(arg)
                hash_expr[arg_hash]
            end

            expr = (node isa CallTerm) ? CallTerm(args_var[1], args_var[2:end], prototype(node)) : IndexTerm(args_var[1], args_var[2], prototype(node))
            assignment = StatementTerm(var_term, expr, prototype(node))
            block = scope_block[scope]
            push!(block.statements, assignment)
            return scopes
            # TODO: check if is sum?
            
        elseif (node isa LambdaTerm) 
            body = node.body
            args = (node.args isa Tuple) ? node.args : (node.args, )
            
            scope_arg = args[end]
            scope = scope_arg.value # take the last arg of lambda function (TODO: do we have 0-arg lambdas?)
            scope_level[scope] = depth + 1
            scope_rank[scope] = 0
            block = BlockTerm([])
            scope_block[scope] = block

            for arg in args 
                arg_hash = Base.hash(arg)
                hash_expr[arg_hash] = arg
                hash_scope[arg_hash] = scope
            end
            
            scopes = visit(body, depth + 1)
            if length(block.statements) == 0 # at least 1 statement
                arg_hash = Base.hash(body)
                arg_var = hash_expr[arg_hash]
                
                push!(block.statements, StatementTerm(LeafTerm(gensym(), prototype(body)), arg_var, prototype(body)))
            end
            scopes = setdiff(scopes,[scope])
            scope = argmax(scope->scope_level[scope],scopes)
            hash_scope[hash] = scope
            rank = 1 + scope_rank[scope]
            scope_rank[scope] = rank

            var = Symbol("var_$(scope)_$(rank)")
            var_term = LeafTerm(var, prototype(node))
            hash_expr[hash] = var_term

            expr = LambdaTerm(block, node.args, x -> prototype(node))

            assignment = StatementTerm(var_term, expr, prototype(node))
            block = scope_block[scope]
            push!(block.statements, assignment)
            return scopes
        else
            error("node type $(typeof(node)) is not implemented!")
        end
    end

    visit(node)
    scope_block[root]
end


function call(f::NewAbstractQuantity, args::NewAbstractQuantity...)
    f(args...)
end

function call(f, args::NewAbstractQuantity...)
    f_quantity = GT.compile_constant_quantity(f)
    f_quantity(args...)
end


function to_quantity(q::NewAbstractQuantity)
    q
end

function to_quantity(q)
    if q isa NewDiscreteField
        discrete_field_quantity(q)
    else
        compile_constant_quantity(q)
    end
    
end

function to_quantities(expr)
    :($to_quantity($expr))
end

function to_quantities(expr::Expr)
    if expr.head == :(->)
        Expr(expr.head, expr.args[1], map(to_quantities, expr.args[2:end])...)
    else
        Expr(expr.head, map(to_quantities, expr.args)...)
    end
end

macro qty(expr)
    expr |> MacroTools.striplines |> to_quantities |> esc
end


# TODO: replace the accessors for the old Measure
function shape_function_accessor(f::typeof(value),space::AbstractSpace,measure::NewMeasure)
    mesh = GT.mesh(measure)
    dom = GT.domain(measure)
    @assert !is_reference_domain(dom)
    d = num_dims(dom)
    @assert num_dims(domain(space)) == d
    # TODO assumes same reference elems for integration and for interpolation
    rid_to_point_to_x = map(coordinates,reference_quadratures(measure.quadrature))
    rid_to_tab = map(rid_to_point_to_x,reference_spaces(space)) do point_to_x, refface
        tabulator(refface)(f,point_to_x)
    end
    face_to_rid = face_reference_id(mesh,d)
    sface_to_face = faces(dom)
    function face_point_dof_s(sface)
        face = sface_to_face[sface]
        rid = face_to_rid[face]
        tab = rid_to_tab[rid]
        function point_dof_s(point,J=nothing)
            function dof_s(dof)
                tab[point,dof]
            end
        end
    end
end


function discrete_field_accessor(f,uh::NewDiscreteField,measure::NewMeasure)
    dom = domain(measure)
    space = GT.space(uh)
    d = num_dims(dom)
    dom2 = domain(space)
    @assert num_dims(dom2) == d
    sface_to_face = faces(dom)
    face_to_rface = inverse_faces(dom)
    rface_to_dofs = face_dofs(space)
    rface_to_rid = face_reference_id(space)
    free_vals = free_values(uh)
    diri_vals = dirichlet_values(uh)
    face_point_dof_s = shape_function_accessor(f,space,measure)
    function face_point_val(sface)
        face = sface_to_face[sface]
        rface = face_to_rface[face]
        dofs = rface_to_dofs[rface]
        rid = rface_to_rid[rface]
        point_dof_s = face_point_dof_s(sface)
        ndofs = length(dofs)
        function point_val(point,J)
            dof_s = point_dof_s(point,J)
            sum(1:ndofs) do i
                dof = dofs[i]
                s = dof_s(i)
                if dof > 0
                    v = free_vals[dof]
                else
                    v = diri_vals[-dof]
                end
                v*s
            end
        end
    end
end


function new_discrete_field(space::AbstractSpace, free_values, dirichlet_values, name=gensym("discrete_field"))
    mesh = space |> GT.mesh
    NewDiscreteField(mesh, space, free_values, dirichlet_values, name)
end


function new_fill_field!(u::NewDiscreteField,z)
    fill!(free_values(u),z)
    fill!(dirichlet_values(u),z)
    u
end

function new_undef_field(::Type{T},space::AbstractSpace) where T
    free_values = allocate_values(T,GT.free_dofs(space))
    dirichlet_values = allocate_values(T,GT.dirichlet_dofs(space))
    new_discrete_field(space,free_values,dirichlet_values)
end

function new_fill_field(z,space::AbstractSpace)
    u = new_undef_field(typeof(z),space)
    new_fill_field!(u,z)
end


function new_zero_field(::Type{T},space::AbstractSpace) where T
    new_fill_field(zero(T),space)
end

struct NewAnalyticalField{A,B,C} <: AbstractField
    definition::A
    quantity::B
    domain::C
end

function new_analytical_field(f,dom::AbstractDomain)
    # TODO: find better way to generate analytical fields and discrete fields
    D = num_dims(dom)
    q = compile_constant_quantity(f)
    NewAnalyticalField(f,q,dom)
end


function interpolate!(f,u::NewDiscreteField)
    interpolate!(f,u,nothing)
end


function interpolate!(f,u::NewDiscreteField,free_or_diri::Union{Nothing,FreeOrDirichlet})
    interpolate_impl!(f,u,free_or_diri)
    # TODO: look into the interpolation to find performance issues
end


function interpolate_impl!(f,u::NewDiscreteField,free_or_diri;location=1)
    
    free_vals = GT.free_values(u)
    diri_vals = GT.dirichlet_values(u)
    space = GT.space(u)
    domain = GT.domain(space)

    d = num_dims(domain)
    dirichlet_dof_location = GT.dirichlet_dof_location(space)
    sface_to_face = faces(domain)
    nfaces = GT.num_faces(domain)
    face_to_dofs = GT.face_dofs(space)
    
    
    # TODO: check correctness and type instability, do code generation with physical maps
    face_to_nodes = face_nodes(space) # 2s
    node_to_x = node_coordinates(space) # 5s
    for sface in 1:nfaces
        face = sface_to_face[sface]
        dofs = face_to_dofs[face]
        ndofs = length(dofs)
        nodes = face_to_nodes[sface]
        for fe_dof in 1:ndofs
            node = nodes[fe_dof]
            x = node_to_x[node]
            gdof = dofs[fe_dof]
            v = f.definition(x)
            if gdof > 0
                if free_or_diri != DIRICHLET
                    free_vals[gdof] = v
                end
            else
                diri_dof = -gdof
                if free_or_diri != FREE && dirichlet_dof_location[diri_dof] == location
                    diri_vals[diri_dof] = v
                end
            end
        end
    end
    u
end


#struct NewIntegral{A}  <: AbstractType
#    contributions::A
#end
#contributions(a::NewIntegral) = a.contributions
#
#
#
#
#"""
#Creates a term representing a placeholder for a (face,point,dof) index.
#"""
#function index_term(name="index";prefix=gensym,int_type::Type{T}=Int) where T
#    IndexTerm(T,prefix(name))
#end
#
#struct IndexTerm{T} <: NewAbstractTerm
#    int_type::Type{T}
#    name::Symbol
#end
#
#function expression!(runtime,t::IndexTerm)
#    t.name
#end
#
#
#abstract type NewAbstractQuantity <: AbstractType end
#
#function new_quantity(term)
#    NewQuantity(term)
#end
#
#struct NewQuantity{A} <: NewAbstractQuantity
#    term::A
#end
#
#function term(q::NewQuantity,opt::TermOptions)
#    q.term(opt)
#end
#
## if it is not a quantity, then process as a constant
#function term(q, opt::TermOptions)
#    NewConstantTerm(q)
#end
#
#function new_constant_quantity(value)
#    new_quantity() do opts
#        NewConstantTerm(value)
#    end
#end
#
#struct NewConstantTerm{A} <: NewAbstractTerm
#    value::A
#end
#
#function expression!(runtime,t::NewConstantTerm)
#    value = compile_symbol!(runtime,t.value)
#    value
#end
#
## function call(f,args::NewAbstractQuantity...)
##    new_quantity() do opts
##        args_term = map(arg->term(arg,opts),args)
##        CallTerm(f,args_term)
##    end
## end
#
## struct CallTerm{A,B} <: NewAbstractTerm
##    callee::A
##    args::B
## end
#
## function (f::NewAbstractQuantity)(x::NewAbstractQuantity)
##    new_quantity() do opts
##        x_term = term(x,opts)
##        f_term = term(f,opts)
##        EvaluateTerm(f_term,x_term)
##    end
## end
##
## struct EvaluateTerm{A,B} <: NewAbstractTerm
##    callee::A
##    arg::B
## end
#
#function new_coordinate_quantity(q::AbstractQuadrature,point::NewAbstractTerm)
#    new_quantity() do opts
#        @assert domain(q) == domain(opts)
#        mesh_face = MeshFaceTerm(domain(opts),domain_face(opts))
#        CoordinateTerm(q,point,mesh_face)
#    end
#end
#
#struct CoordinateTerm{A,B,C} <: NewAbstractTerm
#    quadrature::A
#    point::B
#    mesh_face::C
#end
#
#function expression!(runtime,t::CoordinateTerm)
#    q = t.quadrature
#    rid_point_coordinate_data = map(coordinates,reference_quadratures(q))
#    mesh_face_rid_data = face_reference_id(q)
#    rid_point_coordinate = compile_symbol!(runtime,rid_point_coordinate_data)
#    mesh_face_rid = compile_symbol!(runtime,mesh_face_rid_data)
#    point = expression!(runtime,t.point)
#    mesh_face = expression!(runtime,t.mesh_face)
#    mesh = GT.mesh(GT.domain(q))
#    if is_physical_domain(t.mesh_face.domain)
#        d = num_dims(GT.domain(q))
#        rid_tab_data = map(rid_point_coordinate_data,reference_spaces(mesh,d)) do point_to_x, refface
#            tabulator(refface)(value,point_to_x)
#        end
#        rid_tab = compile_symbol!(runtime,rid_tab_data)
#        node_x = compile_symbol!(runtime,node_coordinates(mesh))
#        mesh_face_nodes = compile_symbol!(runtime,face_nodes(mesh,d))
#        @term begin
#            rid = $mesh_face_reference_id[$mesh_face]
#            tab = $rid_tab[rid]
#            dof_node = $mesh_face_nodes[mesh_face]
#            ndofs = length(dof_node)
#            sum(1:ndofs) do dof
#                node = dof_node[dof]
#                x = $node_x[node]
#                tab[$point,dof]*x
#            end
#        end
#    else
#        @term begin
#            $rid_point_coordinate[$mesh_face_rid[$mesh_face]][$point]
#        end
#    end
#end
#
#struct MeshFaceTerm{A,B} <: NewAbstractTerm
#    domain::A
#    domain_face::B
#end
#
#function expression!(runtime,t::MeshFaceTerm)
#    domain_face_to_mesh_face = compile_symbol!(runtime,faces(t.domain))
#    domain_face = expression!(runtime,t.domain_face)
#    quote
#        $domain_face_to_mesh_face[$domain_face]
#    end
#end
#
#function new_weight_quantity(q::AbstractQuadrature,point::NewAbstractTerm)
#    new_quantity() do opts
#        @assert domain(q) == domain(opts)
#        mesh_face = MeshFaceTerm(domain(opts),domain_face(opts))
#        WeightTerm(q,point,mesh_face)
#    end
#end
#
#struct WeightTerm{A,B,C} <: NewAbstractTerm
#    quadrature::A
#    point::B
#    mesh_face::C
#end
#
#function expression!(runtime,t::WeightTerm)
#    q = t.quadrature
#    rid_point_coordinate_data = map(coordinates,reference_quadratures(q))
#    mesh_face_rid_data = face_reference_id(q)
#    rid_point_weight = compile_symbol!(runtime,map(weights,reference_quadratures(q)))
#    mesh_face_rid = compile_symbol!(runtime,mesh_face_rid_data)
#    point = expression!(runtime,t.point)
#    mesh_face = expression!(runtime,t.mesh_face)
#    mesh = GT.mesh(GT.domain(q))
#    if is_physical_domain(t.mesh_face.domain)
#        d = num_dims(GT.domain(q))
#        rid_tab_data = map(rid_point_coordinate_data,reference_spaces(mesh,d)) do point_to_x, refface
#            tabulator(refface)(ForwardDiff.gradient,point_to_x)
#        end
#        rid_tab = compile_symbol!(runtime,rid_tab_data)
#        node_x = compile_symbol!(runtime,node_coordinates(mesh))
#        mesh_face_nodes = compile_symbol!(runtime,face_nodes(mesh,d))
#        @term begin
#            rid = $mesh_face_rid[$mesh_face]
#            tab = $rid_tab[rid]
#            dof_node = $mesh_face_nodes[$mesh_face]
#            ndofs = length(dof_node)
#            J = sum(1:ndofs) do dof
#                node = dof_node[dof]
#                x = $node_x[node]
#                outer(x,tab[$point,dof]) 
#            end
#            w = $rid_point_weight[rid][$point]
#            w*change_of_measure(J)
#        end
#    else
#        @term begin
#            $rid_point_weight[$mesh_face_rid[$mesh_face]][$point]
#        end
#    end
#end
#
##function new_physical_map_term(domain::AbstractDomain,mesh_face::NewAbstractTerm)
##    mesh = GT.mesh(domain)
##    num_dims = GT.num_dims(domain)
##    NewPhysicalMapTerm(mesh,num_dims,mesh_face)
##end
##
##struct NewPhysicalMapTerm{A,B,C} <: NewAbstractTerm
##    mesh::A
##    num_dims::B
##    mesh_face::C
##end
##
##function expression!(runtime,t::NewPhysicalMapTerm)
##    mesh = t.mesh
##    d = t.num_dims
##    mesh_face_reference_id = compile_symbol!(runtime,face_reference_id(mesh,d))
##    rid_dof_shape_fun = compile_symbol!(runtime,map(shape_functions,reference_spaces(mesh,d)))
##    mesh_face_nodes = compile_symbol!(runtime,face_nodes(mesh,d))
##    node_x = compile_symbol!(runtime,node_coordinates(mesh))
##    @term begin
##        rid = $mesh_face_reference_id[$mesh_face]
##        dof_shape_fun = $rid_dof_shape_fun[rid]
##        dof_node = $mesh_face_nodes[mesh_face]
##        ndofs = length(dof_node)
##        y -> sum(1:ndofs) do dof
##            shape_fun = dof_shape_fun[dof]
##            node = dof_node[dof]
##            x = $node_x[node]
##            shape_fun(y)*x
##        end
##    end
##end
#
## TODO: copied from integration.jl.
#struct NewMeasure{A,B,C,D} <: GT.AbstractType
#    mesh::A
#    quadrature_rule::B
#    domain::C
#    degree::D
#end
#domain(a::NewMeasure) = a.domain
## quadrature_rule(a::NewMeasure) = a.quadrature_rule
#mesh(a::NewMeasure) = a.mesh
## degree(a::NewMeasure) = a.degree
#
#
#function new_measure(dom::AbstractDomain,degree)
#    new_measure(quadrature,dom,degree)
#end
#
#function new_measure(f,dom::AbstractDomain,degree)
#    mesh = GT.mesh(dom)
#    NewMeasure(mesh,f,dom,degree)
#end
#
#function reference_quadratures(measure::NewMeasure)
#    domain = GT.domain(measure)
#    mesh = GT.mesh(domain)
#    d = GT.num_dims(domain)
#    drefid_refdface = GT.reference_spaces(mesh,d)
#    refid_to_quad = map(drefid_refdface) do refdface
#        geo = GT.domain(refdface)
#        measure.quadrature_rule(geo,measure.degree)
#    end
#    refid_to_quad
#end
#
#function face_reference_id(m::NewMeasure)
#    domain = GT.domain(m)
#    mesh = GT.mesh(domain)
#    d = GT.num_dims(domain)
#    face_reference_id(mesh,d)
#end
#
#function quadrature(m::NewMeasure)
#    mesh_quadrature(;
#        domain=domain(m),
#        reference_quadratures = reference_quadratures(m),
#        face_reference_id = face_reference_id(m)
#       )
#end
#
## copy end
#
#struct NumPointsTerm{A, B} <: NewAbstractTerm
#    quadrature::A
#    mesh_face::B
#end
#
#
#function num_points(q::AbstractQuadrature, opts::TermOptions) # TODO: make it a quantity or a term?
#    mesh_face = MeshFaceTerm(domain(opts),domain_face(opts))
#    NumPointsTerm(q, mesh_face)
#end
#
#
#
#struct ParamTerm{A} <: NewAbstractTerm
#    name::Symbol
#    prototype::A
#end
#
#function new_param_quantity(value)
#    name = gensym("param")
#    new_quantity() do opts
#        ParamTerm(name,value) 
#    end
#end
#
#
#function expression!(runtime, t::ParamTerm)
#    get!(runtime,t.name,t.name) # use a unique name for each param term
#    t.name
#end
#
#function expression!(runtime,t::NumPointsTerm)
#    q = t.quadrature
#    rid_point_coordinate_data = map(coordinates,reference_quadratures(q))
#    mesh_face_rid_data = face_reference_id(q)
#    rid_point_weight = compile_symbol!(runtime,map(weights,reference_quadratures(q)))
#    mesh_face_rid = compile_symbol!(runtime,mesh_face_rid_data)
#    # point = expression!(runtime,t.point)
#    mesh_face = expression!(runtime,t.mesh_face)
#    mesh = GT.mesh(GT.domain(q))
#
#    @term begin
#        ws = $rid_point_weight[$mesh_face_rid[$mesh_face]]
#        length(ws)
#    end
#
#end
#
#
#function integrate(f,measure::NewMeasure)
#    point = index_term(:point)
#    q = quadrature(measure)
#    x = new_coordinate_quantity(q,point)
#    dV = new_weight_quantity(q,point)
#    fx = f(x)
#    domain = GT.domain(measure)
#    domain_face = index_term(:domain_face)
#    opts = term_options(domain;domain_face)
#    dV_term = term(dV,opts)
#    fx_term = term(fx,opts)
#    num_points_term = num_points(q, opts)
#    contribution = IntegralTerm(measure,fx_term,dV_term,domain_face,point,num_points_term)
#    contributions = (contribution,)
#    NewIntegral(contributions)
#end
#
#struct IntegralTerm{A,B,C,D,E,F} <: NewAbstractTerm
#    measure::A
#    integrand::B
#    weight::C
#    domain_face::D # remains free
#    point::E # gets reduced
#    num_points::F
#end
#
#measure(t::IntegralTerm) = t.measure
#
#function expression(t::NewAbstractTerm)
#    runtime = IdDict{Any,Symbol}()
#    expr = expression!(runtime,t)
#    expr, runtime
#end
#
#function expression!(runtime,t::IntegralTerm)
#    domain_face = expression!(runtime, t.domain_face)
#    point = expression!(runtime, t.point)
#    fx = expression!(runtime,t.integrand)
#    w = expression!(runtime,t.weight)
#    npoints = expression!(runtime,t.num_points)
#    quote
#        ($domain_face,arg) -> begin
#            $(unpack_runtime_argument(runtime,:arg))
#            sum($point->$fx*$w,1:$npoints)
#        end
#    end
#end
#
#struct NewIntegral{A}  <: AbstractType
#    contributions::A
#end
#contributions(a::NewIntegral) = a.contributions
#
#function Base.sum(int::NewIntegral) # TODO: syntax
#    f = generate_sum(int)
#    f()
#end
#
#function generate_sum(f::Function, params...)
#    param_quantities = map((x -> (typeof(x) <: NewAbstractQuantity) ? x : new_param_quantity(x)), params)
#    int = f(param_quantities...)
#    generate_sum(int, param_quantities...)
#end
#
#function generate_sum(int::NewIntegral,params...)
#    fs = map(c->generate_sum(c,params...),contributions(int))
#    (params2...) -> begin
#        sum(f->f(params2...),fs)
#    end
#end
#
#function generate_sum(ir0::IntegralTerm,params::NewAbstractQuantity...)
#    # @assert length(params) == 0
#    #ir1 = optimize(ir0)
#    #ir2 = lower(ir1)
#    #ir3 = optimize()
#
#    expr, runtime = expression(ir0)
#    
#    dom = GT.domain(ir0.measure)
#    dom_face = ir0.domain_face
#    opts = term_options(dom;domain_face=dom_face)
#    param_terms = map(param -> term(param,opts), params)
#    param_term_exprs = map(x -> x.name, param_terms)
#    arg = runtime_argument(runtime)
#
#    f = eval(expr) # TODO: how can we use the prototype of params?
#
#    nfaces = num_faces(domain(measure(ir0)))
#    (params2...) -> begin 
#        newarg = (;arg..., (zip(param_term_exprs, params2))...)
#
#        sum(domain_face -> invokelatest(f,domain_face,newarg), 1:nfaces)
#    end
#    
#end
#
#function runtime_argument(runtime)
#    (;( key=>val for (val,key) in runtime )...)
#end
#
#function unpack_runtime_argument(runtime,arg)
#    expr = Expr(:block)
#    for k in Base.values(runtime) |> collect |> sort
#        push!(expr.args,:($k = $arg.$k))
#    end
#    expr
#end
#
#function compile_symbol!(runtime,val,name="";prefix=gensym)
#    get!(runtime,val,prefix(name))
#end
#
##int = ∫(α,dΩ)
##f = generate_sum(int,α)
##f(α)
##
##int = ∫(dummy_α,dΩ)
##f = generate_sum(int,dummy_α)
##f(α)
##
## Let us go with this
##int = α -> ∫(α,dΩ)
##f = generate_sum(int,α_prototype)
##f(α)
#
## and this
##int = () -> ∫(α,dΩ)
##f = generate_sum(int)
##f()
#
##int = ∫(α,dΩ)
##f = generate_sum(int)
##f()
#
#
#
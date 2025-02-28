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


const hashed_terms = Set([LeafTerm, CallTerm, IndexTerm, LambdaTerm])


struct StatementTerm{A, B, C} <: NewAbstractTerm
    lhs::A 
    rhs::B
    prototype::C
end

struct BlockTerm{A} <: NewAbstractTerm
    statements::A
end


struct DiscreteFieldTerm{A, B} <: NewAbstractTerm
    discrete_field::A
    prototype::B
    name::Symbol
end

struct TabulatedDiscreteFieldTerm{A, B, C, D, E, F} <: NewAbstractTerm
    linear_operation::A 
    discrete_field::B 
    domain_face::C
    point::D
    measure::E 
    prototype::F
end


prototype(t::BlockTerm) = (length(t.statements) == 0) ? nothing : prototype(t.statements[end])
prototype(t::NewAbstractTerm) = t.prototype


struct NewFormArgumentTerm{A, B, C, D} <: NewAbstractTerm # TODO: rename, check the type
    space::A
    domain_face::B # TODO: is it needed? we have it for coordinate terms
    dof::C
    prototype::D
    arg::Int
    field::Int
    name::Symbol #TODO: is it needed?
end

struct TabulatedFormArgumentTerm{A, B, C, D, E, F} <: NewAbstractTerm # TODO: rename, check the type
    linear_operation::A 
    form_argument::NewFormArgumentTerm
    domain_face::B
    point::C
    dof::D
    measure::E
    prototype::F
end

struct OneFormTerm{A, B, C, D} <: NewAbstractTerm 
    contribution::A
    domain_face::B
    dof::C
    prototype::D
end



space(t::TabulatedFormArgumentTerm) = space(t.form_argument)
space(t::NewFormArgumentTerm) = t.space 

function form_argument(space, arg::Int, field::Int) # TODO: typing
    # arg: test or trial
    # field: block of the global array and matrix 
    @assert arg == 1 && field == 1 # TODO: extend from 1-form with 1 block only
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
    :($(lower(t.args)) -> ($(lower(t.body))))
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

prototype(a::NumPointsTerm) = 0
domain(t::NumPointsTerm) = domain(t.measure)


struct WeightTerm{A,B,C} <: NewAbstractTerm
    measure::A
    point::B
    domain_face::C
end

domain(t::WeightTerm) = domain(t.measure)
prototype(a::WeightTerm) = 0.0


struct CoordinateTerm{A,B,C} <: NewAbstractTerm
    measure::A
    domain_face::B
    point::C
end

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

prototype(t::ContributionTerm) = t.prototype

function rewrite_l1(t::OneFormTerm)
    OneFormTerm(rewrite_l1(t.contribution), t.domain_face, t.dof, t.prototype)
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


function generate_body(t::UniformTerm,name_to_symbol)
    t2 = name_to_symbol[t.name]
    CallTerm(LeafTerm(:value, value), [LeafTerm(t2, t)], prototype(t))
    # :(value($t2))
end


function generate_body(t::LambdaTerm,name_to_symbol)
    LambdaTerm(generate_body(t.body, name_to_symbol), t.args, prototype(t))
end

function generate_body(t::IndexTerm,name_to_symbol)
    IndexTerm(generate_body(t.array, name_to_symbol), generate_body(t.index, name_to_symbol), t.prototype)
end


function generate_body(t::CallTerm,name_to_symbol)
    callee = generate_body(t.callee, name_to_symbol)
    args = map(t.args) do arg 
        generate_body(arg, name_to_symbol)
    end
    CallTerm(callee, args, prototype(t))
end


function generate_body(t::NumPointsTerm,name_to_symbol)
    domain_face_term = t.domain_face
    measure = name_to_symbol[t.measure.name]
    
    face_points = CallTerm(LeafTerm(:num_points_accessor, num_points_accessor), [LeafTerm(measure, t.measure)], x -> prototype(t))
    CallTerm(face_points, [domain_face_term], prototype(t))

    # @term begin
    #     face_npoints = num_points_accessor($measure)
    #     face_npoints($domain_face)
    # end
end


function generate_body(t::WeightTerm,name_to_symbol)
    measure = name_to_symbol[t.measure.name]
    point_term = t.point
    domain_face_term = t.domain_face
    measure_term = LeafTerm(measure, t.measure)

    m = GT.mesh(t.measure.quadrature)
    x = node_coordinates(m)[1]
    jacobian_prototype = x * transpose(x)

    weight_func_0 = CallTerm(LeafTerm(:weight_accessor, weight_accessor), [measure_term], x -> ( (y1, y2) -> 0.0))
    weight_func = CallTerm(weight_func_0, [domain_face_term], ((y1, y2) -> 0.0))

    jacobian_func_0 = CallTerm(LeafTerm(:jacobian_accessor, jacobian_accessor), [measure_term], (x -> y -> jacobian_prototype))
    jacobian_func = CallTerm(jacobian_func_0, [domain_face_term], (y -> jacobian_prototype))
    jacobian_result = CallTerm(jacobian_func, [point_term], jacobian_prototype)

    CallTerm(weight_func, [point_term, jacobian_result], prototype(t))

    # quote
    #     weight_accessor($measure)($domain_face)($point,jacobian_accessor($measure)($domain_face)($point))
    # end
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

function coordinate_quantity(measure::NewMeasure,point)
    new_quantity(;) do dom
        @assert dom == domain(measure)
        CoordinateTerm(measure, nothing, LeafTerm(point, 1))
    end
end


function generate_body(t::CoordinateTerm,name_to_symbol)
    # coordinate_accessor($measure)($domain_face)($point)

    measure = name_to_symbol[t.measure.name]

    measure_term = LeafTerm(measure, t.measure)
    domain_face_term = t.domain_face
    point_term = t.point

    m = GT.mesh(t.measure.quadrature)
    prototype_point = node_coordinates(m)[1]


    point_func_0 = CallTerm(LeafTerm(:coordinate_accessor, coordinate_accessor), [measure_term], x -> y -> prototype_point)
    point_func = CallTerm(point_func_0, [domain_face_term], y -> prototype_point)
    CallTerm(point_func, [point_term], prototype_point)
end


function contribution(f,measure::NewMeasure)
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


function generate_body(t::ContributionTerm,name_to_symbol)
    integrand_expr = generate_body(t.integrand,name_to_symbol)
    weight_expr = generate_body(t.weight,name_to_symbol)
    num_points_expr = generate_body(t.num_points,name_to_symbol)
    point = t.point

    l_body = CallTerm(LeafTerm(:*, *), [integrand_expr, weight_expr], prototype(integrand_expr))
    l = LambdaTerm(l_body, LeafTerm(point, 1), x -> prototype(integrand_expr))
    range = CallTerm(LeafTerm(:(:), :), [LeafTerm(1, 1), num_points_expr], 1:1)
    CallTerm(LeafTerm(:(sum), sum), [l, range], prototype(t))
    # quote
    #     sum($point -> $integrand_expr * $weight_expr,  1:$num_points_expr)
    # end
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
    if !haskey(name_to_symbol, name) # TODO: find the correct form arg
        name_to_captured_data[name] = a
    end
end


function generate_body(t::LeafTerm,name_to_symbol)
    t
end

function generate_body(t::OneFormTerm,name_to_symbol)
    OneFormTerm(generate_body(t.contribution, name_to_symbol), t.domain_face, t.dof, t.prototype) # TODO: array size?
end

function generate_body(t::TabulatedFormArgumentTerm,name_to_symbol)
    # NewFormArgumentTerm not needed?
    measure = name_to_symbol[t.measure.name]

    dof_term = t.dof
    point_term = t.point 
    domain_face_term = t.domain_face
    measure_term = LeafTerm(measure, t.measure)
    
    # space_term = CallTerm(LeafTerm(:space, space), [t.form_argument], t.form_argument.space) # TODO: term?
    space_term = LeafTerm(gensym("space"), t.form_argument.space)

    shape_func_0 = CallTerm(LeafTerm(:shape_function_accessor, shape_function_accessor), [t.linear_operation, space_term, measure_term], (x -> y -> z -> t.prototype))
    shape_func_1 = CallTerm(shape_func_0, [domain_face_term], y -> z -> t.prototype)
    shape_func_2 = CallTerm(shape_func_1, [point_term, LeafTerm(:nothing, nothing)], z -> t.prototype)
    CallTerm(shape_func_2, [dof_term], t.prototype)


end


function generate_body(t::DiscreteFieldTerm,name_to_symbol)
    t_symbol = name_to_symbol[t.name]
    LeafTerm(t_symbol, prototype(t))
    # CallTerm(LeafTerm(:value, value), [LeafTerm(t_symbol, t)], prototype(t))
end


function generate_body(t::TabulatedDiscreteFieldTerm,name_to_symbol)
    measure = name_to_symbol[t.measure.name]
    measure_term = LeafTerm(measure, t.measure)
    domain_face_term = t.domain_face
    point_term = t.point

    m = GT.mesh(t.measure.quadrature)
    x = node_coordinates(m)[1]
    jacobian_prototype = x * transpose(x)


    jacobian_func_0 = CallTerm(LeafTerm(:jacobian_accessor, jacobian_accessor), [measure_term], (x -> y -> jacobian_prototype))
    jacobian_func = CallTerm(jacobian_func_0, [domain_face_term], (y -> jacobian_prototype))
    jacobian_result = CallTerm(jacobian_func, [point_term], jacobian_prototype)

    # point_func_0 = CallTerm(LeafTerm(:coordinate_accessor, coordinate_accessor), [measure_term], x -> y -> prototype_point)
    # point_func = CallTerm(point_func_0, [domain_face_term], y -> prototype_point)
    # CallTerm(point_func, [point_term], prototype_point)

    prototype_result = prototype(t)
    discrete_field = generate_body(t.discrete_field, name_to_symbol)
    field_func_0 = CallTerm(LeafTerm(:discrete_field_accessor, discrete_field_accessor), [t.linear_operation, discrete_field, measure_term], x -> (y1, y2) -> prototype_result) # TODO: proto
    field_func =  CallTerm(field_func_0, [domain_face_term], (y1, y2) -> prototype_result)
    CallTerm(field_func, [point_term, jacobian_result], prototype_result)
end

function generate(t::NewAbstractTerm,params...)
    # t: L1 (domain specific level) term
    domain_face = :dummy_domain_face
    domain_face_term = LeafTerm(domain_face, 1)
    term_l1_1 = set_free_args(t, Dict(:domain_face => domain_face_term))
    term_l1_2 = LambdaTerm(term_l1_1, domain_face_term, x -> prototype(term_l1_1)) # a wrapper of domain face

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

    # generate body and optimize
    
    # rewrite L1 term
    term_l1_3 = rewrite_l1(term_l1_2)
    # L2: functional level
    term_l2 = generate_body(term_l1_3, name_to_symbol)
    
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


  
function set_free_args(t::NewAbstractTerm, args::Dict)
    new_term_args = map(propertynames(t)) do child_name
        child = getproperty(t, child_name)
        if haskey(args, child_name)
            if child === nothing
                args[child_name]
            else 
                child
            end
        elseif child isa NewAbstractTerm
            set_free_args(child, args)
        elseif applicable(iterate, child)
            map(x -> (x isa NewAbstractTerm) ? set_free_args(x, args) : x, child)
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


function generate_1_form(field, int::NewAbstractTerm, params...)
    # TODO: code reuse
    @assert field == 1 # TODO: include more fields

    dof = :dof 
    dof_term = LeafTerm(dof, 1)
    domain_face = :dummy_domain_face
    domain_face_term = LeafTerm(domain_face, 1)

    int2 = set_free_args(int, Dict(:domain_face => domain_face_term, :dof => dof_term))
    t = OneFormTerm(int2, domain_face_term, dof_term, prototype(t))

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

    # generate body and optimize
    # rewrite L1 term
    term_l1 = rewrite_l1(t)

    # L2: functional level
    term_l2 = generate_body(term_l1, name_to_symbol)
    term_l2_wrapped = LambdaTerm(term_l2, domain_face_term, x -> prototype(term_l2))
    
    # L3: imperative level
    term_l3 = statements(term_l2_wrapped)
    
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
    arg = node.args.value
    if arg in args 
        return false
    end
    push!(args, arg)
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
            
        elseif (node isa LambdaTerm) || (node isa OneFormTerm)
            if node isa LambdaTerm
                body = node.body
                args = node.args
            else
                body = node.contribution
                args = node.dof
            end
            scope = args.value # assuming there is only 1 arg in LambdaTerm
            scope_level[scope] = depth + 1
            scope_rank[scope] = 0
            block = BlockTerm([])
            scope_block[scope] = block

            scope_hash = Base.hash(args)
            hash_expr[scope_hash] = args
            hash_scope[scope_hash] = scope

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
            if node isa LambdaTerm
                expr = LambdaTerm(block, node.args, x -> prototype(node))
            else
                expr = OneFormTerm(block, node.domain_face, node.dof, prototype(node))
            end

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
    # TODO: make a quantity of f with @qty 
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
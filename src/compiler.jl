
struct TermInspector{A}
    term::A
end

struct BindingInspector{A}
    binding::A
end

AbstractTrees.children(a::TermInspector) = (map(BindingInspector,bindings(a.term))...,map(TermInspector,dependencies(a.term))...)

function print_term_node(io::IO, node::TermInspector)
    print(io,nameof(typeof(node.term)))
end

function print_term_node(io::IO, node::BindingInspector)
    show(io,node.binding)
end

function print_term_node(io::IO, node)
    show(io,node)
end

function inspect(term::AbstractTerm;maxdepth=5)
    ins = TermInspector(term)
    AbstractTrees.print_tree(print_term_node,stdout,ins;maxdepth)
end

# @term should put every symbol inside
# a LeafNode except the expressions being interpolated,
# which are assumed to to be terms already.
# I have not needed it for the moment.
macro term(expr)
end

const EVALUATED_EXPRS = Dict{Expr,Any}()

function evaluate(expr,captured_data)
    if !haskey(EVALUATED_EXPRS,expr)
        EVALUATED_EXPRS[expr] = eval(expr)
    end
    f1 = EVALUATED_EXPRS[expr]
    f2 = Base.invokelatest(f1, captured_data...)
    f2
end

# No bindings by default
# Bindings are always symbols and dependencies are always terms
function bindings(term::AbstractTerm)
    ()
end

# By default, we just optimize the dependencies
function optimize(term::AbstractTerm)
    dependencies = map(optimize,GT.dependencies(term))
    replace_dependencies(term,dependencies)
end

# We want the leafs to be also AbstractTerm
struct LeafTerm <: AbstractTerm
    value::Any
    prototype::Any
    is_compile_constant::Bool
end

# The prototye is the same as the value by default
leaf_term(v; prototype=v, is_compile_constant = false) = LeafTerm(v, prototype, is_compile_constant)

is_compile_constant(t::LeafTerm) = val_parameter(t.is_compile_constant)

# TODO rename value
value(a::LeafTerm) = a.value
prototype(a::LeafTerm) = a.prototype

AbstractTrees.children(a::TermInspector{<:LeafTerm}) = (a.term.value,)

function dependencies(term::LeafTerm)
    ()
end

function replace_dependencies(term::LeafTerm,dependencies)
    term
end

function expression(term::LeafTerm)
    :( $(term.value) )
end

function lambda_term(body::AbstractTerm,args::Symbol...)
    LambdaTerm(body,args)
end

struct LambdaTerm <: AbstractTerm
    body::AbstractTerm
    args # Tuple of Symbols # all symbols in here get reduced,
end

# TODO this assumes that the prototype
# dos not depend on the arguments
# Maybe it is OK since the arguments
# are typically indices
function prototype(a::LambdaTerm)
    (args...) -> prototype(a.body)
end

# dependencies are always terms
function dependencies(term::LambdaTerm)
    (term.body,)
end

function replace_dependencies(term::LambdaTerm,dependencies)
    (body,) = dependencies
    LambdaTerm(body,term.args)
end

function bindings(term::LambdaTerm)
    term.args
end

function expression(term::LambdaTerm)
    args = Expr(:tuple, term.args...)
    body = expression(term.body)
    :($args -> ($body))
end

function call_term(callee::AbstractTerm,args::AbstractTerm...)
    CallTerm(callee,args)
end

struct CallTerm <: AbstractTerm
    callee::AbstractTerm
    args # Tuple of AbstractTerm
end

function return_prototype(f,args...)
    f(args...)
end

function return_prototype(::typeof(getindex),a,b)
    zero(eltype(a))
end

function return_prototype(::typeof(getindex),a,b::Integer)
    if length(a) != 0
        first(a)
    else
        zero(eltype(a))
    end
end

function prototype(a::CallTerm)
    return_prototype(prototype(a.callee),map(prototype,a.args)...)
end

function dependencies(term::CallTerm)
    (term.callee, term.args...)
end

function replace_dependencies(term::CallTerm,dependencies)
    (callee,args...) = dependencies
    CallTerm(callee,args)
end

function expression(term::CallTerm)
    callee = expression(term.callee)
    args = map(expression,term.args)
    :($(callee)($(args...)))
end

function optimize(term::CallTerm)
    optimize_CallTerm(term.callee,term.args...)
end

function optimize_CallTerm(callee,args...)
    callee2 = optimize(callee)
    args2 = map(optimize,args)
    CallTerm(callee2,args2)
end

function optimize_CallTerm(callee::LeafTerm,args...)
    optimize_CallTerm_LeafTerm(callee.value,callee,args...)
end

function optimize_CallTerm_LeafTerm(value,callee,args...)
    callee2 = optimize(callee)
    args2 = map(optimize,args)
    CallTerm(callee2,args2)
end

struct RefTerm <: AbstractTerm
    container::AbstractTerm
    index::AbstractTerm
end

function optimize(term::RefTerm)
    optimize_RefTerm(term.container,term.index)
end

function optimize_RefTerm(container,index)
    container2 = optimize(container)
    index2 = optimize(index)
    RefTerm(container2,index2)
end

function optimize_CallTerm(callee::RefTerm,args...)
    optimize_CallTerm_RefTerm(callee.container,callee,args...)
end

function optimize_CallTerm_RefTerm(container,callee,args...)
    callee2 = optimize(callee)
    args2 = map(optimize,args)
    CallTerm(callee2,args2)
end

function optimize_CallTerm_LeafTerm(value,callee::LeafTerm,u::RefTerm,x)
    callee2 = optimize(callee)
    u2 = optimize(u)
    x2 = optimize(x)
    term = call_term(callee2,u2,x2)
    optimize(term)
end

function prototype(a::RefTerm)
    #zero(eltype(prototype(a.container)))
    if length(prototype(a.container)) > 0
        first(prototype(a.container))
    else
        zero(eltype(prototype(a.container)))
    end
end

function dependencies(term::RefTerm)
    (term.container,term.index)
end

function replace_dependencies(term::RefTerm,dependencies)
    (container,indices) = dependencies
    RefTerm(container,indices)
end

function expression(term::RefTerm)
    container = expression(term.container)
    index = expression(term.index)
    :($(container)[$(index)])
end

# Parametrizes all values in Leafs if they are in params.
function parametrize(term::AbstractTerm,params...)
    pairs = map(enumerate(params)) do (i,param)
        param=>Symbol("param_$i")
    end
    args = map(last,pairs)
    value_arg = IdDict{Any,Symbol}(pairs...)
    body = parametrize!(value_arg,term)
    lambda_term(body,args...)
end

function parametrize!(value_arg,term::LeafTerm)
    if haskey(value_arg,term.value)
        value = value_arg[term.value]
        leaf_term(value;prototype=prototype(term))
    else
        term
    end
end

function parametrize!(value_arg,term::AbstractTerm)
    dependencies = map(child->parametrize!(value_arg,child),GT.dependencies(term))
    replace_dependencies(term,dependencies)
end

# captures all values in Leafs that are not compile constants
# and have not been reduced
# Make sure that all indices have been reduced
# before using this function!
function capture(term::AbstractTerm)
    value_arg = IdDict{Any,Symbol}()
    reduced_symbols = Base.IdSet{Any}()
    arg = Ref(0)
    body = capture!(value_arg,arg,reduced_symbols,term)
    captured_pairs = collect(value_arg)
    sort!(captured_pairs, by=last)
    captured_values = map(first,captured_pairs)
    captured_args = map(last,captured_pairs)
    captured_term = lambda_term(body,captured_args...)
    captured_term, captured_values
end

function capture!(value_arg,arg,reduced_symbols,term::LeafTerm)
    if is_compile_constant(term)
        return term
    end
    (;value) = term
    if value in reduced_symbols
        return term
    end
    if ! haskey(value_arg,value)
        arg[] += 1
        value_arg[value] = Symbol("captured_arg_$(arg[])")
        leaf_term(value_arg[value];prototype=prototype(term))
    else
        leaf_term(value_arg[value];prototype=prototype(term))
    end
end

function capture!(value_arg,arg,reduced_symbols,term::AbstractTerm)
    bindings = GT.bindings(term)
    if ! isempty(bindings)
        push!(reduced_symbols,bindings...)
    end
    dependencies = map(child->capture!(value_arg,arg,reduced_symbols,child),GT.dependencies(term))
    replace_dependencies(term,dependencies)
end

# Create default index for a form of a given arity.
function index(arity)
    domain_face = :domain_face
    point = :point
    arg_field = ntuple(i->Symbol("field_$i"),arity)
    arg_face_around = ntuple(i->Symbol("face_around_$i"),arity)
    arg_dof = ntuple(i->Symbol("dof_$i"),arity)
    Index(domain_face,point,arg_field,arg_face_around,arg_dof)
end

# All indices that get reduces in a form of arity N
# involving integrals,
# when using numerical integration to compute the integrals
struct Index{N} <: AbstractTerm
    domain_face_index::Symbol
    point_index::Symbol
    arg_field_index::NTuple{N,Symbol}
    arg_face_around_index::NTuple{N,Symbol}
    arg_dof_index::NTuple{N,Symbol}
end

domain_face_index(index::Index) = index.domain_face_index
point_index(index::Index) = index.point_index
field_index(index::Index,arg) = index.arg_field_index[arg]
face_around_index(index::Index,arg) = index.arg_face_around_index[arg]
dof_index(index::Index,arg) = index.arg_dof_index[arg]

# This is the data needed to transform a
# quantity into a term
struct QuantityOptions{A,B} <: AbstractType
    domain::A
    index::B
end

domain(a::QuantityOptions) = a.domain
index(a::QuantityOptions) = a.index

# The objects that are exposed to the user
# inside of an integral

# NB.To define a prototype for a quantity  we need a domain
# It is easier to add the prototype to terms instead
struct Quantity <: AbstractQuantity
    term::Any
end
function quantity(term)
    Quantity(term)
end

function term(qty::Quantity,opts)
    qty.term(opts)
end

function call(callee,args::AbstractQuantity...)
    # TODO: this is a solution for the lambda generated by face_quantity. We need a better call function to simplify this
    is_lambda = (callee isa Function) && !isdefined(Base, Symbol(callee)) && 
            !isdefined(GT, Symbol(callee)) && !isdefined(LinearAlgebra, Symbol(callee)) && !isdefined(ForwardDiff, Symbol(callee))
    f = uniform_quantity(callee;is_compile_constant=!is_lambda)
    call(f,args...)
end

function call(callee::AbstractQuantity,args::AbstractQuantity...)
    quantity() do opts
        f_term = term(callee,opts)
        args_term = map(arg->term(arg,opts),args)
        call_term(f_term,args_term...)
    end
end

function (f::AbstractQuantity)(x::AbstractQuantity)
    call(f,x)
end

function uniform_quantity(v;is_compile_constant=false)
    quantity() do opts
        leaf_term(v; is_compile_constant)
    end
end

function parameter(v)
    ParameterQuantity(v)
end

function parameter(v::DiscreteField)
    v
end

struct ParameterQuantity{A} <: AbstractQuantity
    value::A
end

value(q::ParameterQuantity) = q.value

function term(q::ParameterQuantity,opts)
    t1 = leaf_term(q)
    t2 = call_term(GT.value,t1)
end

function coordinate_quantity(quadrature)
    quantity() do opts
        face = domain_face_index(opts.index)
        point = point_index(opts.index)
        dependencies = map(leaf_term,(quadrature,face,point))
        CoordinateTerm(dependencies)
    end
end

struct CoordinateTerm <: AbstractTerm
    dependencies
end

function prototype(term::CoordinateTerm)
    (quadrature,face,point) = map(prototype,dependencies(term))
    prototype(coordinate_accessor(quadrature))
end

# dependencies are always terms
# Everything needed to compute the value
# including indices
function dependencies(term::CoordinateTerm)
    term.dependencies
end

function replace_dependencies(term::CoordinateTerm,dependencies)
    CoordinateTerm(dependencies)
end

# Not clear if this can contain statements or not
# Not clear if this is even needed for DLS terms
function expression(term::CoordinateTerm)
    (quadrature,face,point) = map(expression,term.dependencies)
    :( coordinate_accessor($quadrature)($face)($point) )
end

function weight_quantity(quadrature)
    quantity() do opts
        face = domain_face_index(opts.index)
        point = point_index(opts.index)
        dependencies = map(leaf_term,(quadrature,face,point))
        D = (num_dims(domain(quadrature)))
        WeightTerm(dependencies, D) # this is required to tabulate jacobian accessors, so that we cannot reuse the compiled code for a different dimension size
    end
end

struct WeightTerm <: AbstractTerm
    dependencies
    D::Int
end

function prototype(term::WeightTerm)
    (quadrature,face,point) = map(prototype,dependencies(term))
    prototype(weight_accessor(quadrature))
end

function dependencies(term::WeightTerm)
    term.dependencies
end

function replace_dependencies(term::WeightTerm,dependencies)
    WeightTerm(dependencies, term.D)
end

function expression(term::WeightTerm)
    (quadrature,face,point) = map(expression,term.dependencies)
    D = term.D
    J = :(jacobian_accessor($quadrature, $(Val(D)))($face)($point))
    :(weight_accessor($quadrature)($face)($point, $J) )
end

function form_argument_quantity(space::AbstractSpace,arg,the_field=1)
    quantity() do opts
        space_domain = GT.domain(space)
        domain = GT.domain(opts)
        D = num_dims(space_domain)
        d = num_dims(domain)
        face_around = GT.face_around(domain)
        index = GT.index(opts)
        face = domain_face_index(index)
        field = field_index(index,arg)
        dof = dof_index(index,arg)

        f = leaf_term(GT.value; is_compile_constant=true)
        the_field_term = leaf_term(the_field; is_compile_constant=true)
        (space_term, domain_term, face_term, field_term, dof_term, face_around_term) = map(leaf_term, (space, domain, face, field, dof, face_around))
        to_leaf_term = x -> (x isa AbstractTerm) ? x : leaf_term(x)

        if D == d
            # the_face_around_term = leaf_term(nothing, is_compile_constant=true)
            the_face_around = 1
            the_face_around_term = leaf_term(the_face_around;is_compile_constant=true)
            # face_around_term = leaf_term(nothing)
            face_around_term = leaf_term(face_around_index(index,arg))
            dependencies = (f,space_term,domain_term,face_term,the_field_term,field_term,dof_term,the_face_around_term,face_around_term)
            FormArgumentTerm(dependencies, D)
        elseif D==d+1 && face_around !== nothing
            # the_face_around_term = leaf_term(face_around, is_compile_constant=true)
            # face_around_term = leaf_term(face_around)
            the_face_around = 1
            the_face_around_term = leaf_term(the_face_around;is_compile_constant=true)
            # face_around_term = leaf_term(nothing)
            face_around_term = leaf_term(face_around_index(index,arg))
            dependencies = (f,space_term,domain_term,face_term,the_field_term,field_term,dof_term,the_face_around_term,face_around_term)
            FormArgumentTerm(dependencies, D)
        else
            face_around_term = leaf_term(face_around_index(index,arg))
            the_face_around = :the_face_around
            the_face_around_term = leaf_term(the_face_around)
            n_faces_around = leaf_term(2;is_compile_constant=true) # Hard coded! But OK in practice.
            dependencies = (f,space_term,domain_term,face_term,the_field_term,field_term,dof_term,the_face_around_term,face_around_term)
            SkeletonTerm(FormArgumentTerm(dependencies, D),n_faces_around,the_face_around)
        end
    end
end

struct FormArgumentTerm <: AbstractTerm
    dependencies
    D::Int
    #d::Int
end

function prototype(term::FormArgumentTerm)
    (f,space,domain,face,the_field,field,dof,the_face_around,face_around) = map(prototype,dependencies(term))
    prototype(form_argument_accessor(space,domain,1))
end

function dependencies(term::FormArgumentTerm)
    term.dependencies
end

function replace_dependencies(term::FormArgumentTerm,dependencies)
    FormArgumentTerm(dependencies, term.D)
end

function replace_the_face_around(term::FormArgumentTerm,the_face_around)
    (f,space,domain,face,the_field,field,dof,_,face_around) = term.dependencies
    dependencies = (f,space,domain,face,the_field,field,dof,the_face_around,face_around)
    FormArgumentTerm(dependencies, term.D)
end

function expression(term::FormArgumentTerm)
    (f,space,domain,face,the_field,field,dof,the_face_around,face_around) = map(expression,term.dependencies)
    # TODO: not tabulated. we need another solution to tabulate it or this cannot be tabulated anyway
    :(form_argument_accessor($f,$space,$domain,$the_field)($face,$the_face_around)($dof,$field,$face_around))
end

#function optimize(term::CallTerm{<:FormArgumentTerm,<:Tuple{<:CoordinateTerm}})
#    coords = optimize(term.args[1])
#    (quadrature,face,point) = map(optimize,dependencies(coords))
#    TabulatedTerm(term.callee,quadrature,point)
#end

function optimize_CallTerm(callee::FormArgumentTerm,coords::CoordinateTerm)
    (quadrature,face,point) = map(optimize,dependencies(coords))
    TabulatedTerm(callee,quadrature,point)
end

for op in [:value,:(ForwardDiff.gradient),:(ForwardDiff.jacobian)]
    @eval begin
        function optimize_CallTerm_LeafTerm(
                ::typeof($op),callee::LeafTerm,v::FormArgumentTerm,coords::CoordinateTerm)
            (quadrature,face,point) = map(optimize,dependencies(coords))
            (f,space,domain,face,the_field,field,dof,the_face_around,face_around) = dependencies(v)
            new_deps = (callee,space,domain,face,the_field,field,dof,the_face_around,face_around)
            v2 = replace_dependencies(v, new_deps)
            TabulatedTerm(v2,quadrature,point)
        end
    end
end

#function optimize(term::CallTerm{<:LeafTerm{typeof(ForwardDiff.gradient)}, <:Tuple{<:FormArgumentTerm, <:CoordinateTerm}})
#    v = optimize(term.args[1]) 
#    coords = optimize(term.args[2])
#    new_deps = (term.callee,space,domain,face,the_field,field,dof,the_face_around,face_around)
#    v2 = replace_dependencies(v, new_deps)
#    TabulatedTerm(v2,quadrature,point)
#end

function term(uh::DiscreteField,opts)
    space = GT.space(uh)
    space_domain = GT.domain(space)
    domain = GT.domain(opts)
    D = num_dims(space_domain)
    d = num_dims(domain)
    the_face_around = GT.face_around(domain)
    index = GT.index(opts)
    face = domain_face_index(index)
    f = value
    if D == d
        face_around = nothing
        dependencies = map(leaf_term,(f,uh,domain,face,face_around))
        DiscreteFieldTerm(dependencies)
    elseif D==d+1 && face_around !== nothing
        dependencies = map(leaf_term,(f,uh,domain,face,face_around))
        DiscreteFieldTerm(dependencies)
    else
        face_around = :the_face_around
        n_faces_around = leaf_term(2) # Hard coded! But OK in practice.
        dependencies = map(leaf_term,(f,uh,domain,face,face_around))
        SkeletonTerm(DiscreteFieldTerm(dependencies),n_faces_around,face_around)
    end
end

struct DiscreteFieldTerm <: AbstractTerm
    dependencies
end

function dependencies(term::DiscreteFieldTerm)
    term.dependencies
end

function replace_dependencies(term::DiscreteFieldTerm,dependencies)
    DiscreteFieldTerm(dependencies)
end

function optimize_CallTerm(callee::DiscreteFieldTerm,coords::CoordinateTerm)
    (quadrature,face,point) = map(optimize,dependencies(coords))
    TabulatedTerm(callee,quadrature,point)
end

#function optimize(term::CallTerm{<:DiscreteFieldTerm,<:Tuple{<:CoordinateTerm}})
#    coords = term.args[1]
#    (quadrature,face,point) = dependencies(coords)
#    TabulatedTerm(term.callee,quadrature,point)
#end

for op in [:value,:(ForwardDiff.gradient),:(ForwardDiff.jacobian)]
    @eval begin
        function optimize_CallTerm_LeafTerm(
                ::typeof($op),callee::LeafTerm,v::DiscreteFieldTerm,coords::CoordinateTerm)
            (quadrature,face,point) = map(optimize,dependencies(coords))
            (f,uh,domain,face,face_around) = dependencies(v)
            new_deps = (callee,uh,domain,face,face_around)
            v2 = replace_dependencies(v, new_deps)
            TabulatedTerm(v2,quadrature,point)
        end
    end
end

#function optimize(term::CallTerm{<:LeafTerm{typeof(ForwardDiff.gradient)}, <:Tuple{<:DiscreteFieldTerm, <:CoordinateTerm}})
#    v = optimize(term.args[1]) 
#    coords = optimize(term.args[2])
#    (quadrature,face,point) = map(optimize,dependencies(coords))
#    (f,uh,domain,face,face_around) = dependencies(v)
#    new_deps = (term.callee,uh,domain,face,face_around)
#    v2 = replace_dependencies(v, new_deps)
#    TabulatedTerm(v2,quadrature,point)
#end

struct TabulatedTerm <: AbstractTerm
    parent::AbstractTerm
    quadrature::AbstractTerm
    point::AbstractTerm
end

function dependencies(term::TabulatedTerm)
    (term.parent,term.quadrature,term.point)
end

function replace_dependencies(term::TabulatedTerm,dependencies)
    parent, quadrature, point = dependencies
    TabulatedTerm(parent,quadrature, point)
end

function prototype(term::TabulatedTerm)
    prototype_TabulatedTerm(term.parent,prototype(term.quadrature),prototype(term.point))
end

function expression(term::TabulatedTerm)
    expression_TabulatedTerm(term.parent,expression(term.quadrature),expression(term.point))
end

function prototype_TabulatedTerm(parent::FormArgumentTerm,quadrature,point)
    (f,space,domain,face,the_field,field,dof,the_face_around,face_around) = map(prototype,dependencies(parent))
    prototype(form_argument_accessor(f,space,quadrature))
end

#function prototype(term::TabulatedTerm{<:FormArgumentTerm})
#    parent = term.parent
#    quadrature = prototype(term.quadrature)
#    (f,space,domain,face,the_field,field,dof,the_face_around,face_around) = map(prototype,dependencies(parent))
#    prototype(form_argument_accessor(space,domain,1))
#end

function prototype_TabulatedTerm(parent::DiscreteFieldTerm,quadrature,point)
    (f,uh,domain,face,face_around) = map(prototype,dependencies(parent))
    prototype(discrete_field_accessor(f,uh,quadrature))
end

#function prototype(term::TabulatedTerm{<:DiscreteFieldTerm})
#    parent = term.parent
#    quadrature = prototype(term.quadrature)
#    (f,uh,domain,face,face_around) = map(prototype,dependencies(parent))
#    prototype(discrete_field_accessor(f,uh,quadrature))
#end

function expression_TabulatedTerm(parent::FormArgumentTerm,quadrature,point)
    form_arg = parent
    (f,space,domain,face,the_field,field,dof,the_face_around,face_around) = map(expression,form_arg.dependencies)
    D = parent.D
    J = :(jacobian_accessor($quadrature, $(Val(D)))($face, $the_face_around)($point))
    # TODO: inline form Accessor. How can we know whether it is physical or reference?
    form_argument_accessor_term(f,space,quadrature,the_field, face,the_face_around, point,J, dof,field,face_around)
    # :(form_argument_accessor($f,$space,$quadrature,$the_field)($face,$the_face_around)($point, $J)($dof,$field,$face_around))
end

#function expression(term::TabulatedTerm{<:FormArgumentTerm})
#    point = expression(term.point)
#    form_arg = term.parent
#    quadrature = expression(term.quadrature)
#    (f,space,domain,face,the_field,field,dof,the_face_around,face_around) = map(expression,form_arg.dependencies)
#    :(form_argument_accessor($f,$space,$quadrature,$the_field)($face,$the_face_around)($point)($dof,$field,$face_around))
#end

function expression_TabulatedTerm(parent::DiscreteFieldTerm,quadrature,point)
    form_arg = parent
    (f,uh,domain,face,face_around) = map(expression,form_arg.dependencies)
    :(discrete_field_accessor($f,$uh,$quadrature)($face,$face_around)($point))
end

#function expression(term::TabulatedTerm{<:DiscreteFieldTerm})
#    point = expression(term.point)
#    form_arg = term.parent
#    quadrature = expression(term.quadrature)
#    (f,uh,domain,face,face_around) = map(expression,form_arg.dependencies)
#    :(discrete_field_accessor($f,$uh,$quadrature)($face,$face_around)($point))
#end

struct SkeletonTerm <: AbstractTerm
    parent::AbstractTerm
    n_faces_around::AbstractTerm
    the_face_around::Symbol # Gets reduced
end

function bindings(term::SkeletonTerm)
    (term.the_face_around,)
end

function dependencies(term::SkeletonTerm)
    (term.parent,term.n_faces_around)
end

function replace_dependencies(term::SkeletonTerm,dependencies)
    (parent,n_faces_around) = dependencies
    SkeletonTerm(parent,n_faces_around,term.the_face_around)
end

function expression(term::SkeletonTerm)
    # Not needed in practice.
    # In theory this is a tuple of length num faces arround (2)
    error()
end

function optimize_RefTerm(container::SkeletonTerm,index)
    the_face_around = index
    replace_the_face_around(container.parent,the_face_around)
end

#function optimize(term::RefTerm{<:SkeletonTerm})
#    the_face_around = term.index
#    replace_the_face_around(term.container.parent,the_face_around)
#end

function optimize_CallTerm_RefTerm(container::SkeletonTerm,callee,args...)
    callee = optimize(callee)
    term2 = CallTerm(callee,args)
    optimize(term2)
end

#function optimize(term::CallTerm{<:RefTerm{<:SkeletonTerm}})
#    callee = optimize(term.callee)
#    term2 = CallTerm(callee,term.args)
#    optimize(term2)
#end

function physical_map(mesh::AbstractMesh,vD)
    quantity() do opts
        D = val_parameter(vD)
        domain = GT.domain(opts)
        d = num_dims(domain)
        face = GT.domain_face_index(opts.index)
        @assert d == D
        f = GT.value
        dependencies = map(leaf_term,(f,mesh,vD,face))
        PhysicalMapTerm(dependencies)
    end
end

struct PhysicalMapTerm <: AbstractTerm
    dependencies
end

function dependencies(term::PhysicalMapTerm)
    term.dependencies
end

function replace_dependencies(term::PhysicalMapTerm,dependencies)
    PhysicalMapTerm(dependencies)
end

function optimize_CallTerm(callee::PhysicalMapTerm,coords::CoordinateTerm)
    (quadrature,face,point) = map(optimize,dependencies(coords))
    TabulatedTerm(callee,quadrature,point)
end

#function optimize(term::CallTerm{<:PhysicalMapTerm,<:Tuple{<:CoordinateTerm}})
#    coords = optimize(term.args[1])
#    (quadrature,face,point) = map(optimize,dependencies(coords))
#    TabulatedTerm(term.callee,quadrature,point)
#end

function optimize_CallTerm_LeafTerm(
    ::typeof(ForwardDiff.jacobian),callee::LeafTerm,v::PhysicalMapTerm,coords::CoordinateTerm)
    (quadrature,face,point) = map(optimize,dependencies(coords))
    (f,mesh,vD,face) = dependencies(v)
    new_deps = (callee,mesh,vD,face)
    v2 = replace_dependencies(v, new_deps)
    TabulatedTerm(v2,quadrature,point)
end

#function optimize(term::CallTerm{<:LeafTerm{typeof(ForwardDiff.jacobian)}, <:Tuple{<:PhysicalMapTerm,<:CoordinateTerm}})
#    v = optimize(term.args[1]) 
#    coords = optimize(term.args[2])
#    (quadrature,face,point) = map(optimize,dependencies(coords))
#    (f,mesh,vD,face) = dependencies(v)
#    new_deps = (term.callee,mesh,vD,face)
#    v2 = replace_dependencies(v, new_deps)
#    TabulatedTerm(v2,quadrature,point)
#end

function prototype_TabulatedTerm(parent::PhysicalMapTerm,quadrature,point)
    (f,mesh,vD,face) = map(prototype,dependencies(parent))
    prototype(physical_map_accessor(f,quadrature,vD))
end

#function prototype(term::TabulatedTerm{<:PhysicalMapTerm})
#    parent = term.parent
#    quadrature = prototype(term.quadrature)
#    (f,mesh,vD,face) = map(prototype,dependencies(parent))
#    prototype(physical_map_accessor(f,quadrature,vD))
#end

function expression_TabulatedTerm(parent::PhysicalMapTerm,quadrature,point)
    (f,mesh,vD,face) = map(expression,parent.dependencies)
    :(physical_map_accessor($f,$quadrature,$vD)($face)($point))
end

#
#function expression(term::TabulatedTerm{<:PhysicalMapTerm})
#    point = expression(term.point)
#    phys_map = term.parent
#    quadrature = expression(term.quadrature)
#    (f,mesh,vD,face) = map(expression,phys_map.dependencies)
#    :(physical_map_accessor($f,$quadrature,$vD)($face)($point))
#end

#function reference_map(mesh::AbstractMesh,vd,vD)
#    quantity() do opts
#        d = val_parameter(vd)
#        D = val_parameter(vD)
#        domain = GT.domain(opts)
#        @assert  d == num_dims(domain)
#        @assert d != D
#        face = GT.domain_face_index(opts.index)
#        the_face_around = :the_face_around
#        f = GT.value
#        dependencies = map(leaf_term,(f,mesh,vd,vD,face,the_face_around))
#        parent = ReferenceMapTerm(dependencies)
#        n_faces_around = leaf_term(2;is_compile_constant=true)
#        SkeletonTerm(parent,n_faces_around,the_face_around)
#    end
#end
#
#struct ReferenceMapTerm{A} <: AbstractTerm
#    dependencies::A
#end
#
#function dependencies(term::ReferenceMapTerm)
#    term.dependencies
#end
#
#function replace_dependencies(term::ReferenceMapTerm,dependencies)
#    ReferenceMapTerm(dependencies)
#end
#
#function replace_the_face_around(term::ReferenceMapTerm,the_face_around)
#    (f,mesh,vd,vD,face,_) = term.dependencies
#    deps = (f,mesh,vd,vD,face,the_face_around)
#    replace_dependencies(term,deps)
#end
#
#function optimize(term::CallTerm{<:ReferenceMapTerm,<:Tuple{<:CoordinateTerm}})
#    coords = optimize(term.args[1])
#    (quadrature,face,point) = map(optimize,dependencies(coords))
#    TabulatedTerm(term.callee,quadrature,point)
#end
#
#function prototype(term::TabulatedTerm{<:ReferenceMapTerm})
#    parent = term.parent
#    quadrature = prototype(term.quadrature)
#    (f,mesh,vD,face) = map(prototype,dependencies(parent))
#    prototype(physical_map_accessor(f,quadrature,vD))
#end
#
#function expression(term::TabulatedTerm{<:ReferenceMapTerm})
#    point = expression(term.point)
#    phys_map = term.parent
#    quadrature = expression(term.quadrature)
#    (f,mesh,vD,face) = map(expression,phys_map.dependencies)
#    :(physical_map_accessor($f,$quadrature,$vD)($face)($point))
#end

function unit_normal(mesh::AbstractMesh,vd)
    quantity() do opts
        d = val_parameter(vd)
        D = num_dims(mesh)
        face = domain_face_index(opts.index)
        @assert d == num_dims(domain(opts))
        if GT.face_around(domain(opts)) === nothing
            the_face_around = :the_face_around
            dependencies = map(leaf_term,(face,the_face_around))
            parent = UnitNormalTerm(dependencies)
            n_faces_around = leaf_term(2;is_compile_constant=true)
            SkeletonTerm(parent,n_faces_around,the_face_around)
        else
            the_face_around = nothing
            dependencies = map(leaf_term,(face,the_face_around))
            UnitNormalTerm(dependencies)
        end
    end
end

struct UnitNormalTerm <: AbstractTerm
    dependencies
end

function dependencies(term::UnitNormalTerm)
    term.dependencies
end

function replace_dependencies(term::UnitNormalTerm,deps)
    UnitNormalTerm(deps)
end

function replace_the_face_around(term::UnitNormalTerm,the_face_around)
    (face,_) = GT.dependencies(term)
    newdeps = (face,the_face_around)
    replace_dependencies(term,newdeps)
end

function optimize_CallTerm(callee::UnitNormalTerm,coords::CoordinateTerm)
    (quadrature,face,point) = map(optimize,dependencies(coords))
    TabulatedTerm(callee,quadrature,point)
end

#function optimize(term::CallTerm{<:UnitNormalTerm,<:Tuple{<:CoordinateTerm}})
#    coords = term.args[1]
#    (quadrature,face,point) = map(optimize,dependencies(coords))
#    TabulatedTerm(term.callee,quadrature,point)
#end

function prototype_TabulatedTerm(parent::UnitNormalTerm,quadrature,point)
    prototype(unit_normal_accessor(quadrature))
end

#function prototype(term::TabulatedTerm{<:UnitNormalTerm})
#    parent = term.parent
#    quadrature = prototype(term.quadrature)
#    #(face,the_face_around) = map(prototype,dependencies(parent))
#    prototype(unit_normal_accessor(quadrature))
#end

function expression_TabulatedTerm(parent::UnitNormalTerm,quadrature,point)
    (face,the_face_around) = map(expression,dependencies(parent))
    :(unit_normal_accessor($quadrature)($face,$the_face_around)($point))
end

#function expression(term::TabulatedTerm{<:UnitNormalTerm})
#    parent = term.parent
#    quadrature = expression(term.quadrature)
#    point = expression(term.point)
#    (face,the_face_around) = map(expression,dependencies(parent))
#    :(unit_normal_accessor($quadrature)($face,$the_face_around)($point))
#end

function dual_basis_quantity(space::AbstractSpace)
    quantity() do opts
        domain_space = GT.domain(space)
        domain = GT.domain(opts)
        d = num_dims(domain)
        index = opts.index
        domain_face = leaf_term(domain_face_index(index))
        dof = leaf_term(dof_index(index,1))
        domain_face_to_face = call_term(map(leaf_term,(GT.faces,domain))...)
        face = RefTerm(domain_face_to_face,domain_face)
        @assert num_dims(domain_space) == d
        face_to_rid = call_term(map(leaf_term,(face_reference_id,space))...)
        rid = RefTerm(face_to_rid,face)
        rid_to_fe = call_term(map(leaf_term,(reference_spaces,space))...)
        fe = RefTerm(rid_to_fe,rid)
        dof_to_sigma = call_term(leaf_term(dual_basis),fe)
        sigma = RefTerm(dof_to_sigma,dof)
        sigma
    end
end

function face_quantity(data,mesh::AbstractMesh,vd;reference=Val(false))
    quantity() do opts
        d = val_parameter(vd)
        domain = GT.domain(opts)
        index = GT.index(opts)
        domain_face = leaf_term(domain_face_index(index))
        @assert num_dims(domain) == d
        domain_face_to_face = call_term(map(leaf_term,(GT.faces,domain))...)
        face = RefTerm(domain_face_to_face,domain_face)
        if val_parameter(reference)
            face_to_rid = call_term(map(leaf_term,face_reference_id,mesh,d)...)
            rid_to_value = leaf_term(data)
            rid = RefTerm(face_to_rid,face)
            value = RefTerm(rid_to_value,rid)
        else
            face_to_value = leaf_term(data)
            value = RefTerm(face_to_value,face)
        end
    end
end

function generate_assemble_scalar(contribution::DomainContribution;parameters=())
    term_0 = write_assemble_scalar(contribution)
    term_1 = optimize(term_0)
    term_2 = parametrize(term_1,parameters...)
    term_3, captured_data = capture(term_2)
    expr_0 = expression(term_3)
    expr_1 = statements_expr(expr_0)
    f = evaluate(expr_1,captured_data)
    f
end

function write_assemble_scalar(contribution::DomainContribution)
    quadrature = GT.quadrature(contribution)
    domain = GT.domain(quadrature)
    arity = Val(0)
    index = GT.index(arity)
    term = GT.term(contribution,index)
    quadrature_term = leaf_term(quadrature)
    init_arg = :init
    init = leaf_term(init_arg)
    ScalarAssemblyTerm(term,quadrature_term,init,init_arg,index)
end


function generate_assemble_vector(contribution::DomainContribution,space::AbstractSpace;parameters=())
    # TODO how to sort the calls to optimize, parametrize, capture, expression, statements, etc
    # needs to be though carefully. We provably will need several calls to optimize and more IR levels
    term_0 = write_assemble_vector(contribution,space)
    term_1 = optimize(term_0)
    term_2 = parametrize(term_1,parameters...)
    term_3, captured_data = capture(term_2)
    expr_0 = expression(term_3)
    f = evaluate(expr_0,captured_data)
    # expr_1 = statements_expr(expr_0)
    # f = evaluate(expr_1,captured_data)
    f
end

function write_assemble_vector(contribution::DomainContribution,space::AbstractSpace)
    quadrature = GT.quadrature(contribution)
    domain = GT.domain(quadrature)
    arity = Val(1)
    index = GT.index(arity)
    term = GT.term(contribution,index)
    space_term = leaf_term(space)
    quadrature_term = leaf_term(quadrature)
    alloc_arg = :alloc
    alloc = leaf_term(alloc_arg)

    field_n_faces_around = map(fields(space)) do field_space
        max_num_faces_around(GT.domain(field_space),domain)
    end

    VectorAssemblyTerm(term,space_term,quadrature_term,alloc,alloc_arg,index,field_n_faces_around)
end


function generate_assemble_matrix(contribution::DomainContribution,space_trial::AbstractSpace,space_test::AbstractSpace;parameters=())
    # TODO how to sort the calls to optimize, parametrize, capture, expression, statements, etc
    # needs to be though carefully. We provably will need several calls to optimize and more IR levels
    term_0 = write_assemble_matrix(contribution,space_trial, space_test)
    term_1 = optimize(term_0)
    term_2 = parametrize(term_1,parameters...)
    term_3, captured_data = capture(term_2)
    expr_0 = expression(term_3)
    f = evaluate(expr_0,captured_data)
    # expr_1 = statements_expr(expr_0)
    # f = evaluate(expr_1,captured_data)
    f
end

function write_assemble_matrix(contribution::DomainContribution,space_trial::AbstractSpace,space_test::AbstractSpace)
    quadrature = GT.quadrature(contribution)
    domain = GT.domain(quadrature)
    arity = Val(2)
    index = GT.index(arity)
    term = GT.term(contribution,index)
    space_trial_term = leaf_term(space_trial)
    space_test_term = leaf_term(space_test)
    quadrature_term = leaf_term(quadrature)
    alloc_arg = :alloc
    alloc = leaf_term(alloc_arg)

    field_n_faces_around_trial = map(fields(space_trial)) do field_space
        max_num_faces_around(GT.domain(field_space),domain)
    end

    field_n_faces_around_test = map(fields(space_test)) do field_space
        max_num_faces_around(GT.domain(field_space),domain)
    end

    MatrixAssemblyTerm(term,space_trial_term,space_test_term,quadrature_term,alloc,alloc_arg,index,field_n_faces_around_trial,field_n_faces_around_test)
end

function generate_sample(f,quadrature::AbstractQuadrature;parameters=())
    term_0 = write_sample(f,quadrature)
    term_1 = optimize(term_0)
    term_2 = parametrize(term_1,parameters...)
    term_3, captured_data = capture(term_2)
    expr_0 = expression(term_3)
    expr_1 = statements_expr(expr_0)
    f = evaluate(expr_1,captured_data)
    f
end

function write_sample(f,quadrature::AbstractQuadrature)
    x = coordinate_quantity(quadrature)
    fx = f(x)
    domain = GT.domain(quadrature)
    arity = Val(0)
    index = GT.index(arity)
    opts = QuantityOptions(domain,index)
    term = GT.term(fx,opts)
    vals_arg = :vals
    vals = leaf_term(vals_arg)
    SampleTerm(term,index,vals,vals_arg)
end

struct SampleTerm <: AbstractTerm
    term::AbstractTerm
    index::Index{0}
    vals::AbstractTerm
    vals_arg::Symbol
end

function bindings(term::SampleTerm)
    index = term.index
    face = domain_face_index(index)
    point = point_index(index)
    (face,point,term.vals_arg)
end

function dependencies(term::SampleTerm)
    (;term,vals) = term
    (term,vals)
end

function replace_dependencies(term::SampleTerm,dependencies)
    (term2,vals) = dependencies
    (;index,vals_arg) = term
    SampleTerm(term2,index,vals,vals_arg)
end

function expression(term::SampleTerm)
    (term2,vals) = map(expression,dependencies(term))
    (face,point,vals_arg) = bindings(term)
    point_v = :($point -> $term2)
    face_point_v = :($face -> $point_v)
    :( $vals_arg ->  sample_loop!($vals,$face_point_v))
end

function sample_loop!(face_point_w,face_point_v)
    nfaces = length(face_point_w)
    for face in 1:nfaces
        point_v = face_point_v(face)
        point_w = face_point_w[face]
        npoints = length(point_w)
        for point in 1:npoints
            point_w[point] = point_v(point)
        end
    end
    face_point_w
end

function generate_assemble_face_contribution(contribution::DomainContribution;parameters=())
    term_0 = write_assemble_face_contribution(contribution)
    term_1 = optimize(term_0)
    term_2 = parametrize(term_1,parameters...)
    term_3, captured_data = capture(term_2)
    expr_0 = expression(term_3)
    expr_1 = statements_expr(expr_0)
    f = evaluate(expr_1,captured_data)
    f
end

function write_assemble_face_contribution(contribution::DomainContribution)
    quadrature = GT.quadrature(contribution)
    domain = GT.domain(quadrature)
    arity = Val(0)
    index = GT.index(arity)
    term = GT.term(contribution,index)
    quadrature_term = leaf_term(quadrature)
    contribs_arg = :init
    contribs = leaf_term(contribs_arg)
    FaceContributionTerm(term,quadrature_term,contribs,contribs_arg,index)
end

# Like ScalarAssemblyTerm but we store intermediate face results into a vector contribs
struct FaceContributionTerm <: AbstractTerm
    term::AbstractTerm
    quadrature::AbstractTerm
    contribs::AbstractTerm
    contribs_arg::Symbol # gets reduced
    index::Index{0} # gets reduced
end

function bindings(term::FaceContributionTerm)
    index = term.index
    face = domain_face_index(index)
    point = point_index(index)
    (term.contribs_arg,face,point)
end

function dependencies(term::FaceContributionTerm)
    (;term,quadrature,contribs) = term
    (term,quadrature,contribs)
end

function replace_dependencies(term::FaceContributionTerm,dependencies)
    (term2,quadrature,contribs) = dependencies
    (;contribs_arg,index) = term
    FaceContributionTerm(term2,quadrature,contribs,contribs_arg,index)
end

function expression(term::FaceContributionTerm)
    (term2,quadrature,contribs) = map(expression,dependencies(term))
    (contribs_arg,face,point) = bindings(term)
    point_v = :($point -> $term2)
    face_point_v = :($face -> $point_v)
    body = :(face_contribution_loop!($contribs,$face_point_v,$quadrature))
    expr = :($contribs_arg->$body)
    expr
end

function face_contribution_loop!(contribs,face_point_v,quadrature)
    domain = GT.domain(quadrature)
    nfaces = num_faces(domain)
    face_npoints = num_points_accessor(quadrature)
    init = zero(eltype(contribs))
    for face in 1:nfaces
        point_v = face_point_v(face)
        npoints = face_npoints(face)
        z = init
        for point in 1:npoints
            z += point_v(point)
        end
        contribs[face] = z
    end
    contribs
end

struct ScalarAssemblyTerm <: AbstractTerm
    term::AbstractTerm
    quadrature::AbstractTerm
    init::AbstractTerm
    init_arg::Symbol
    index::Index{0}
end

struct VectorAssemblyTerm <: AbstractTerm
    term::AbstractTerm
    space::AbstractTerm
    quadrature::AbstractTerm
    alloc::AbstractTerm
    alloc_arg::Symbol
    index::Index{1}
    field_n_faces_around::Tuple
end

struct MatrixAssemblyTerm <: AbstractTerm
    term::AbstractTerm
    space_trial::AbstractTerm
    space_test::AbstractTerm
    quadrature::AbstractTerm
    alloc::AbstractTerm
    alloc_arg::Symbol # gets reduced
    index::Index{2} # gets reduced
    field_n_faces_around_trial::Tuple
    field_n_faces_around_test::Tuple
end

function bindings(term::ScalarAssemblyTerm)
    index = term.index
    face = domain_face_index(index)
    point = point_index(index)
    (term.init_arg,face,point)
end

function dependencies(term::ScalarAssemblyTerm)
    (;term,quadrature,init) = term
    (term,quadrature,init)
end

function bindings(term::VectorAssemblyTerm)
    index = term.index
    face = domain_face_index(index)
    point = point_index(index)
    field = field_index(index,1)
    face_around = face_around_index(index,1)
    dof = dof_index(index,1)
    (term.alloc_arg,face,point,field,face_around,dof)
end

function dependencies(term::VectorAssemblyTerm)
    (;term,space,quadrature,alloc) = term
    (term,space,quadrature,alloc)
end

function bindings(term::MatrixAssemblyTerm)
    index = term.index
    face = domain_face_index(index)
    point = point_index(index)
    field_trial, field_test = field_index(index,1), field_index(index,2)
    face_around_trial, face_around_test = face_around_index(index,1), face_around_index(index,2)
    dof_trial, dof_test = dof_index(index,1), dof_index(index,2)
    (term.alloc_arg,face,point,field_trial,field_test,face_around_trial,face_around_test,dof_trial,dof_test)
end

function dependencies(term::MatrixAssemblyTerm)
    (;term,space_trial,space_test,quadrature,alloc) = term
    (term,space_trial,space_test,quadrature,alloc)
end


function replace_dependencies(term::ScalarAssemblyTerm,dependencies)
    (term2,quadrature,init) = dependencies
    (;init_arg,index) = term
    ScalarAssemblyTerm(term2,quadrature,init,init_arg,index)
end

function replace_dependencies(term::VectorAssemblyTerm,dependencies)
    (term2,space,quadrature,alloc) = dependencies
    (;alloc_arg,index,field_n_faces_around) = term
    VectorAssemblyTerm(term2,space,quadrature,alloc,alloc_arg,index,field_n_faces_around)
end


function replace_dependencies(term::MatrixAssemblyTerm,dependencies)
    (term2,space_trial,space_test,quadrature,alloc) = dependencies
    (;alloc_arg,index,field_n_faces_around_trial,field_n_faces_around_test) = term
    MatrixAssemblyTerm(term2,space_trial,space_test,quadrature,alloc,alloc_arg,index,field_n_faces_around_trial,field_n_faces_around_test)
end


function expression(term::ScalarAssemblyTerm)
    (term2,quadrature,init) = map(expression,dependencies(term))
    (init_arg,face,point) = bindings(term)
    point_v = :($point -> $term2)
    face_point_v = :($face -> $point_v)
    body = :(scalar_assembly_loop($face_point_v,$init,$quadrature))
    expr = :($init_arg->$body)
    expr
end

# We can specialize for particular cases if needed
function expression(term::VectorAssemblyTerm)
    (term2,space,quadrature,alloc) = map(expression,dependencies(term))
    (alloc_arg,face,point,field,face_around,dof) = bindings(term)
    # dof_v = :($dof -> $term2)
    # block_dof_v = :( ($field,$face_around) -> $dof_v)
    # point_block_dof_v = :($point -> $block_dof_v)
    # face_point_block_dof_v = :( $face -> $point_block_dof_v)
    # body = :(vector_assembly_loop!($face_point_block_dof_v,$alloc,$space,$quadrature))
    # expr = :($alloc_arg->$body)
    # expr

    loop_vars = (face, point, field, face_around, dof)
    # expr_L = topological_sort(term2, loop_vars)
    # body = generate_vector_assembly_loop_body(term.field_n_faces_around, expr_L, (space,quadrature,alloc), loop_vars)
    # expr = :($alloc_arg->$body)

    # expr_L_bitmap = topological_sort_bitmap(term2, loop_vars) 
    # body_bitmap = generate_vector_assembly_loop_body_bitmap(term.field_n_faces_around, expr_L_bitmap, (space, quadrature, alloc), loop_vars)
    # expr = :($alloc_arg -> $body_bitmap)

    body = generate_vector_assembly_template(term2, term.field_n_faces_around, (space, quadrature, alloc), loop_vars)

    nfields = length(term.field_n_faces_around)
    ndofs = [:(GT.max_num_reference_dofs(GT.field($space,$i))) for i in 1:nfields]

    loop_var_range = Dict(dof => :(max($(ndofs...),)) )

    body_optimized = ast_optimize(body, loop_var_range)
    expr = :($alloc_arg->$body_optimized)

    expr
end


function expression(term::MatrixAssemblyTerm)
    (term2,space_trial,space_test,quadrature,alloc) = map(expression,dependencies(term))
    (alloc_arg,face,point,field_trial,field_test,face_around_trial,face_around_test,dof_trial,dof_test) = bindings(term)

    # V0: hand-written
    # dof_v = :(($dof_trial) -> ($dof_test) -> $term2)
    # block_dof_v = :( ($field_trial,$field_test,$face_around_trial,$face_around_test) -> $dof_v)
    # point_block_dof_v = :($point -> $block_dof_v)
    # face_point_block_dof_v = :( $face -> $point_block_dof_v)
    # body = :(matrix_assembly_loop!($face_point_block_dof_v,$alloc,$space_trial,$space_test,$quadrature))
    # expr = :($alloc_arg->$body)
    # expr
    loop_vars = (face, point, field_trial, field_test, face_around_trial, face_around_test, dof_trial, dof_test)

    # V1: normal hoisting
    # expr_L = topological_sort(term2, loop_vars) 
    # body = generate_matrix_assembly_loop_body(term.field_n_faces_around_trial, term.field_n_faces_around_test, expr_L, (space_trial, space_test, quadrature, alloc), loop_vars)
    # # display(body)
    # expr = :($alloc_arg->$body)

    # V2: bitmap tabulation
    # expr_L_bitmap = topological_sort_bitmap(term2, loop_vars) 
    # body_bitmap = generate_matrix_assembly_loop_body_bitmap(term.field_n_faces_around_trial, term.field_n_faces_around_test, expr_L_bitmap, (space_trial, space_test, quadrature, alloc), loop_vars)
    # expr = :($alloc_arg -> $body_bitmap)


    # display(body_bitmap)

    # V3: split
    # generate loops for term
    body = generate_matrix_assembly_template(term2, term.field_n_faces_around_trial, term.field_n_faces_around_test, (space_trial, space_test, quadrature, alloc), loop_vars)
    # rewrite

    nfields_trial = length(term.field_n_faces_around_trial)
    nfields_test = length(term.field_n_faces_around_test)
    ndofs_test = [:(GT.max_num_reference_dofs(GT.field($space_test,$i))) for i in 1:nfields_test]
    ndofs_trial = [:(GT.max_num_reference_dofs(GT.field($space_trial,$i))) for i in 1:nfields_trial]

    loop_var_range = Dict(dof_trial => :(max($(ndofs_trial...),)), 
                            dof_test => :(max($(ndofs_test...),)))

    body_optimized = ast_optimize(body, loop_var_range)
    # display(body_optimized)
    expr = :($alloc_arg->$body_optimized)
    
    expr
end

function scalar_assembly_loop(face_point_v,init,quadrature)
    domain = GT.domain(quadrature)
    nfaces = num_faces(domain)
    face_npoints = num_points_accessor(quadrature)

    z = init
    for face in 1:nfaces
        point_v = face_point_v(face)
        npoints = face_npoints(face)
        for point in 1:npoints
            z += point_v(point)
        end
    end
    z
end

function vector_assembly_loop!(face_point_block_dof_v,alloc,space,quadrature)
    domain = GT.domain(quadrature)
    nfaces = num_faces(domain)
    #space_domain = GT.domain(space)
    face_npoints = num_points_accessor(quadrature)
    field_face_n_faces_around = map(fields(space)) do field_space
        num_faces_around_accesor(GT.domain(field_space),domain)
    end
    field_face_dofs = map(fields(space)) do field_space
        dofs_accessor(field_space,domain)
    end
    field_face_around_be = map(fields(space)) do field_space
        n = max_num_faces_around(GT.domain(field_space),domain)
        map(1:n) do _
            m = max_num_reference_dofs(field_space)
            zeros(eltype(alloc),m)
        end
    end
    nfields = num_fields(space)
    z = zero(eltype(alloc))
    for face in 1:nfaces
        point_block_dof_v = face_point_block_dof_v(face)
        # Reset
        for face_around_be in field_face_around_be
            for be in face_around_be
                fill!(be,z)
            end
        end
        # Integrate
        npoints = face_npoints(face)
        for point in 1:npoints
            block_dof_v = point_block_dof_v(point)
            for field in 1:nfields
                face_around_be = field_face_around_be[field]
                n_faces_around = field_face_n_faces_around[field](face)
                for face_around in 1:n_faces_around
                    be = face_around_be[face_around]
                    dof_v = block_dof_v(field,face_around)
                    dofs = field_face_dofs[field](face,face_around)
                    ndofs = length(dofs)
                    for dof in 1:ndofs
                        v = dof_v(dof)
                        be[dof] += v
                    end
                end
            end
        end
        # Contribute
        for field in 1:nfields
            face_around_be = field_face_around_be[field]
            n_faces_around = field_face_n_faces_around[field](face)
            for face_around in 1:n_faces_around
                dofs = field_face_dofs[field](face,face_around)
                be = face_around_be[face_around]
                contribute!(alloc,be,dofs,field)
            end
        end
    end
    alloc
end

function matrix_assembly_loop!(face_point_block_dof_v,alloc,space_trial,space_test,quadrature)
    domain = GT.domain(quadrature)
    nfaces = num_faces(domain)
    # space_domain = GT.domain(space)
    face_npoints = num_points_accessor(quadrature)
    field_face_n_faces_around_trial = map(fields(space_trial)) do field_space
        num_faces_around_accesor(GT.domain(field_space),domain)
    end
    field_face_n_faces_around_test = map(fields(space_test)) do field_space
        num_faces_around_accesor(GT.domain(field_space),domain)
    end

    field_face_dofs_trial = map(fields(space_trial)) do field_space
        dofs_accessor(field_space,domain)
    end

    field_face_dofs_test = map(fields(space_test)) do field_space
        dofs_accessor(field_space,domain)
    end

    field_face_around_be = map(fields(space_trial)) do field_space_trial
        m_trial = max_num_reference_dofs(field_space_trial)
        n_trial = max_num_faces_around(GT.domain(field_space_trial),domain)
        map(fields(space_test)) do field_space_test
            m_test = max_num_reference_dofs(field_space_test)
            n_test = max_num_faces_around(GT.domain(field_space_test),domain)
            # fix the type of tuple(tuple(vector(vector(matrix))))
            n_matrix::Vector{Vector{Matrix{eltype(alloc)}}} = map(1:n_trial) do _
                map(1:n_test) do _
                    zeros(eltype(alloc),m_test,m_trial)
                end
            end
            n_matrix
        end
    end

    nfields_trial = num_fields(space_trial)
    nfields_test = num_fields(space_test)
    z = zero(eltype(alloc))

    for face in 1:nfaces
        point_block_dof_v = face_point_block_dof_v(face)
        # Reset
        for field_test_face_around_be in field_face_around_be
            for face_around_be in field_test_face_around_be
                for face_around_test_be in face_around_be
                    for be in face_around_test_be
                        fill!(be,z)
                    end
                end
            end
        end
        # Integrate
        npoints = face_npoints(face)
        for point in 1:npoints
            block_dof_v = point_block_dof_v(point)
            for field_trial in 1:nfields_trial
                n_faces_around_trial = field_face_n_faces_around_trial[field_trial](face)
                for field_test in 1:nfields_test
                    face_around_be = field_face_around_be[field_trial][field_test]
                    n_faces_around_test = field_face_n_faces_around_test[field_test](face)
                    
                    for face_around_trial in 1:n_faces_around_trial
                        for face_around_test in 1:n_faces_around_test
                            be = face_around_be[face_around_trial][face_around_test]
                            dof_v = block_dof_v(field_trial,field_test,face_around_trial,face_around_test)

                            dofs_trial = field_face_dofs_trial[field_trial](face,face_around_trial)
                            dofs_test = field_face_dofs_test[field_test](face,face_around_test)
                            ndofs_trial, ndofs_test = length(dofs_trial), length(dofs_test)
                            for dof_trial in 1:ndofs_trial
                                dof_test_v = dof_v(dof_trial)
                                for dof_test in 1:ndofs_test
                                    # v = dof_v(dof_trial, dof_test)
                                    v = dof_test_v(dof_test)
                                    be[dof_test, dof_trial] += v 
                                end
                            end
                        end
                    end
                end
            end
        end

        # Contribute
        for field_trial in 1:nfields_trial
            n_faces_around_trial = field_face_n_faces_around_trial[field_trial](face)
            for field_test in 1:nfields_test
                face_around_be = field_face_around_be[field_trial][field_test]
                n_faces_around_test = field_face_n_faces_around_test[field_test](face)
                for face_around_trial in 1:n_faces_around_trial
                    for face_around_test in 1:n_faces_around_test
                        dofs_trial = field_face_dofs_trial[field_trial](face,face_around_trial)
                        dofs_test = field_face_dofs_test[field_test](face,face_around_test)
                        be = face_around_be[face_around_trial][face_around_test]
                        contribute!(alloc,be,dofs_test, dofs_trial,field_test, field_trial)
                    end
                end
            end
        end
    end
    alloc
end

function last_symbol(expr_L)
    for block in reverse(expr_L)
        statements = block.args
        if length(statements) > 0
            return statements[end].args[1]
        end
    end
    return nothing # if not found
end

function alloc_zeros(symbol, args...)
    if args[1] === Any
        Array{args[1]}(undef, args[2:end])
    else
        zeros(args...)
    end
end

function generate_vector_assembly_loop_body(field_n_faces_around, expr_L, dependencies, bindings)
    # expr_L: a list of block in the order of (global, face, point, field, face_around, dof)
    # dependencies: (space, quadrature, alloc)
    # bindings: (face, point, field, face_around, dof)
    (space, quadrature, alloc) = dependencies
    (face, point, field_symbol, face_around_symbol, dof) = bindings # TODO: replace field_symbol and face_around_symbol
    v = last_symbol(expr_L)
    block = Expr(:block)
    nfields = length(field_n_faces_around)

    assignment = :(domain = GT.domain($quadrature))
    push!(block.args,assignment)

    assignment = :(T = eltype($alloc))
    push!(block.args, assignment)

    assignment = :(z = zero(T))
    push!(block.args, assignment)

    assignment = :(nfaces = GT.num_faces(domain))
    push!(block.args,assignment)

    assignment = :(face_npoints = GT.num_points_accessor($quadrature))
    push!(block.args,assignment)

    # Get single field spaces
    for field in 1:nfields
        var = Symbol("space_for_$(field)")
        expr = :(GT.field($space,$field))
        assignment = :($var = $expr)
        push!(block.args, assignment)
    end

    # Get face_dofs
    for field in 1:nfields
        space_field = Symbol("space_for_$(field)")
        var_field = Symbol("face_dofs_for_$(field)")
        expr = :(GT.dofs_accessor($space_field,domain))
        assignment = :($var_field = $expr)
        push!(block.args, assignment)
    end

    # Allocate face vectors
    for field in 1:nfields
        n_faces_around = field_n_faces_around[field]
        for face_around in 1:n_faces_around
            space_field = Symbol("space_for_$(field)")
            var_str = "be_for_$(field)_$(face_around)"
            var_face_around = Symbol(var_str)
            expr = :(alloc_zeros($var_str, T,GT.max_num_reference_dofs($space_field)))
            assignment = :($var_face_around = $expr)
            push!(block.args, assignment)
        end
    end

    # statements without deps
    push!(block.args, expr_L[1].args...)

    # Face loop
    face_loop_head = :($face = 1:nfaces)
    face_loop_body = Expr(:block)
    face_loop = Expr(:for,face_loop_head,face_loop_body)
    push!(block.args,face_loop)

    assignment = :(npoints = face_npoints($face))
    push!(face_loop_body.args,assignment)

    # statements depend on face
    push!(face_loop_body.args, expr_L[2].args...)

    for field in 1:nfields
        n_faces_around = field_n_faces_around[field]
        for face_around in 1:n_faces_around
            dofs = Symbol("dofs_for_$(field)_$(face_around)")
            accessor = Symbol("face_dofs_for_$(field)")
            assignment = :($dofs = $(accessor)($face,$face_around))
            push!(face_loop_body.args,assignment)
        end
    end

    for field in 1:nfields
        n_faces_around = field_n_faces_around[field]
        for face_around in 1:n_faces_around
            dofs = Symbol("dofs_for_$(field)_$(face_around)")
            ndofs = Symbol("n_dofs_for_$(field)_$face_around")
            assignment = :($ndofs = length($(dofs)))
            push!(face_loop_body.args,assignment)
        end
    end

    for field in 1:nfields
        n_faces_around = field_n_faces_around[field]
        for face_around in 1:n_faces_around
            be = Symbol("be_for_$(field)_$(face_around)")
            assignment = :(fill!($be,z))
            push!(face_loop_body.args,assignment)
        end
    end

    # Point loop
    point_loop_head = :($point = 1:npoints)
    point_loop_body = Expr(:block)
    point_loop = Expr(:for,point_loop_head,point_loop_body)
    push!(face_loop_body.args,point_loop)

    # assignment = :(block_dof_v = point_block_dof_v(point))
    # push!(point_loop_body.args,assignment)
    # statements depend on point
    push!(point_loop_body.args, expr_L[3].args...)

    for field in 1:nfields
        n_faces_around = field_n_faces_around[field]
        # statements depend on field
        push!(point_loop_body.args, :($field_symbol = $field))
        push!(point_loop_body.args, expr_L[4].args...)
        
        for face_around in 1:n_faces_around
            # accessor = Symbol("dofs_v_for_$(field)_$(face_around)")
            # assignment = :( $accessor = block_dof_v($field,$face_around)  )
            # push!(point_loop_body.args,assignment)
            # statements depend on face around
            push!(point_loop_body.args, :($face_around_symbol = $face_around))
            push!(point_loop_body.args, expr_L[5].args...)
            
            # DOF loop
            # do just after the field and face_around computation, or otherwise we need to replace the variable names in expr_L 
            ndofs = Symbol("n_dofs_for_$(field)_$face_around")
            dof_loop_head = :($dof = 1:$ndofs)
            dof_loop_body = Expr(:block)
            dof_loop = Expr(:for,dof_loop_head,dof_loop_body)
            push!(point_loop_body.args,dof_loop)
            be = Symbol("be_for_$(field)_$(face_around)")
            # accessor = Symbol("dofs_v_for_$(field)_$(face_around)")
            # statements depend on dof
            push!(dof_loop_body.args, expr_L[6].args...)
            assignment = :($be[$dof] += $v)
            push!(dof_loop_body.args,assignment)
        end
    end

    # for field in 1:nfields
    #     n_faces_around = field_n_faces_around[field]
    #     for face_around in 1:n_faces_around
    #         ndofs = Symbol("n_dofs_for_$(field)_$face_around")
    #         dof_loop_head = :(dof = 1:$ndofs)
    #         dof_loop_body = Expr(:block)
    #         dof_loop = Expr(:for,dof_loop_head,dof_loop_body)
    #         push!(point_loop_body.args,dof_loop)
    #         be = Symbol("be_for_$(field)_$(face_around)")
    #         accessor = Symbol("dofs_v_for_$(field)_$(face_around)")
    #         assignment = :($be[dof] += $accessor(dof))
    #         push!(dof_loop_body.args,assignment)
    #     end
    # end

    for field in 1:nfields
        n_faces_around = field_n_faces_around[field]
        for face_around in 1:n_faces_around
            dofs = Symbol("dofs_for_$(field)_$(face_around)")
            be = Symbol("be_for_$(field)_$(face_around)")
            assignment = :(GT.contribute!(alloc,$be,$dofs,$field))
            push!(face_loop_body.args,assignment)
        end
    end

    block
end




function generate_matrix_assembly_loop_body(field_n_faces_around_trial, field_n_faces_around_test, expr_L, dependencies, bindings)
    # expr_L: a list of block in the order of (global, face, point, field_trial, field_test, face_around_trial, face_around_test, dof_trial, dof_test)
    # dependencies: (space_trial, space_test, quadrature, alloc)
    # bindings: (face, point, field_trial, field_test, face_around_trial, face_around_test, dof_trial, dof_test)
    (space_trial, space_test, quadrature, alloc) = dependencies
    (face, point, field_trial_symbol, field_test_symbol, face_around_trial_symbol, face_around_test_symbol, dof_trial, dof_test) = bindings # TODO: replace field_symbol and face_around_symbol
    v = last_symbol(expr_L)
    block = Expr(:block)
    nfields_trial = length(field_n_faces_around_trial)
    nfields_test = length(field_n_faces_around_test)

    assignment = :(domain = GT.domain($quadrature))
    push!(block.args,assignment)

    assignment = :(T = eltype($alloc))
    push!(block.args, assignment)

    assignment = :(z = zero(T))
    push!(block.args, assignment)

    assignment = :(nfaces = GT.num_faces(domain))
    push!(block.args,assignment)

    assignment = :(face_npoints = GT.num_points_accessor($quadrature))
    push!(block.args,assignment)

    # Get single field spaces 
    for field in 1:nfields_trial
        var = Symbol("space_for_$(field)_trial")
        expr = :(GT.field($space_trial,$field))
        assignment = :($var = $expr)
        push!(block.args, assignment)
    end

    for field in 1:nfields_test
        var = Symbol("space_for_$(field)_test")
        expr = :(GT.field($space_test,$field))
        assignment = :($var = $expr)
        push!(block.args, assignment)
    end

    # Get face_dofs
    for field in 1:nfields_trial
        space_field = Symbol("space_for_$(field)_trial")
        var_field = Symbol("face_dofs_for_$(field)_trial")
        expr = :(GT.dofs_accessor($space_field,domain))
        assignment = :($var_field = $expr)
        push!(block.args, assignment)
    end

    for field in 1:nfields_test
        space_field = Symbol("space_for_$(field)_test")
        var_field = Symbol("face_dofs_for_$(field)_test")
        expr = :(GT.dofs_accessor($space_field,domain))
        assignment = :($var_field = $expr)
        push!(block.args, assignment)
    end

    
    # Allocate face vectors
    for field_trial in 1:nfields_trial
        n_faces_around_trial = field_n_faces_around_trial[field_trial]
        for field_test in 1:nfields_test
            n_faces_around_test = field_n_faces_around_test[field_test]
            for face_around_trial in 1:n_faces_around_trial
                for face_around_test in 1:n_faces_around_test
                    space_field_trial = Symbol("space_for_$(field_trial)_trial")
                    space_field_test = Symbol("space_for_$(field_test)_test")
                    var_str = "be_for_$(field_trial)_$(field_test)_$(face_around_trial)_$(face_around_test)"
                    var_face_around = Symbol(var_str)
                    expr = :(alloc_zeros($var_str,T,GT.max_num_reference_dofs($space_field_test),GT.max_num_reference_dofs($space_field_trial)))
                    assignment = :($var_face_around = $expr)
                    push!(block.args, assignment)
                end
            end
        end
    end

    # statements without deps
    push!(block.args, expr_L[1].args...)

    # Face loop
    face_loop_head = :($face = 1:nfaces)
    face_loop_body = Expr(:block)
    face_loop = Expr(:for,face_loop_head,face_loop_body)
    push!(block.args,face_loop)

    assignment = :(npoints = face_npoints($face))
    push!(face_loop_body.args,assignment)


    # statements depend on face
    push!(face_loop_body.args, expr_L[2].args...)

    # trial
    for field in 1:nfields_trial
        n_faces_around = field_n_faces_around_trial[field]
        for face_around in 1:n_faces_around
            dofs = Symbol("dofs_for_$(field)_$(face_around)_trial")
            accessor = Symbol("face_dofs_for_$(field)_trial")
            assignment = :($dofs = $(accessor)($face,$face_around))
            push!(face_loop_body.args,assignment)
        end
    end

    # test
    for field in 1:nfields_test
        n_faces_around = field_n_faces_around_test[field]
        for face_around in 1:n_faces_around
            dofs = Symbol("dofs_for_$(field)_$(face_around)_test")
            accessor = Symbol("face_dofs_for_$(field)_test")
            assignment = :($dofs = $(accessor)($face,$face_around))
            push!(face_loop_body.args,assignment)
        end
    end

    # trial
    for field in 1:nfields_trial
        n_faces_around = field_n_faces_around_trial[field]
        for face_around in 1:n_faces_around
            dofs = Symbol("dofs_for_$(field)_$(face_around)_trial")
            ndofs = Symbol("n_dofs_for_$(field)_$(face_around)_trial")
            assignment = :($ndofs = length($(dofs)))
            push!(face_loop_body.args,assignment)
        end
    end

    # test
    for field in 1:nfields_test
        n_faces_around = field_n_faces_around_test[field]
        for face_around in 1:n_faces_around
            dofs = Symbol("dofs_for_$(field)_$(face_around)_test")
            ndofs = Symbol("n_dofs_for_$(field)_$(face_around)_test")
            assignment = :($ndofs = length($(dofs)))
            push!(face_loop_body.args,assignment)
        end
    end


    for field_trial in 1:nfields_trial
        n_faces_around_trial = field_n_faces_around_trial[field_trial]
        for field_test in 1:nfields_test
            n_faces_around_test = field_n_faces_around_test[field_test]
            for face_around_trial in 1:n_faces_around_trial
                for face_around_test in 1:n_faces_around_test
                    var_str = "be_for_$(field_trial)_$(field_test)_$(face_around_trial)_$(face_around_test)"
                    be = Symbol(var_str)
                    assignment = :(fill!($be,z))
                    push!(face_loop_body.args,assignment)
                end
            end
        end
    end


    # Point loop
    point_loop_head = :($point = 1:npoints)
    point_loop_body = Expr(:block)
    point_loop = Expr(:for,point_loop_head,point_loop_body)
    push!(face_loop_body.args,point_loop)


    # statements depend on point
    push!(point_loop_body.args, expr_L[3].args...)



    for field_trial in 1:nfields_trial
        n_faces_around_trial = field_n_faces_around_trial[field_trial]
        # statements depend on field
        push!(point_loop_body.args, :($field_trial_symbol = $field_trial))
        push!(point_loop_body.args, expr_L[4].args...)
        for field_test in 1:nfields_test
            n_faces_around_test = field_n_faces_around_test[field_test]
            # statements depend on field
            push!(point_loop_body.args, :($field_test_symbol = $field_test))
            push!(point_loop_body.args, expr_L[5].args...)

            for face_around_trial in 1:n_faces_around_trial
                # statements depend on face around
                push!(point_loop_body.args, :($face_around_trial_symbol = $face_around_trial))
                push!(point_loop_body.args, expr_L[6].args...)
                for face_around_test in 1:n_faces_around_test
                    # statements depend on face around
                    push!(point_loop_body.args, :($face_around_test_symbol = $face_around_test))
                    push!(point_loop_body.args, expr_L[7].args...)


                    # DOF loop
                    # do just after the field and face_around computation, or otherwise we need to replace the variable names in expr_L 
                    ndofs_trial = Symbol("n_dofs_for_$(field_trial)_$(face_around_trial)_trial")
                    dof_trial_loop_head = :($dof_trial = 1:$ndofs_trial)
                    dof_trial_loop_body = Expr(:block)
                    dof_trial_loop = Expr(:for,dof_trial_loop_head,dof_trial_loop_body)
                    push!(point_loop_body.args,dof_trial_loop)
                    push!(dof_trial_loop_body.args, expr_L[8].args...)


                    ndofs_test = Symbol("n_dofs_for_$(field_test)_$(face_around_test)_test")
                    dof_test_loop_head = :($dof_test = 1:$ndofs_test)
                    dof_test_loop_body = Expr(:block)
                    dof_test_loop = Expr(:for,dof_test_loop_head,dof_test_loop_body)
                    push!(dof_trial_loop_body.args,dof_test_loop)
                    push!(dof_test_loop_body.args, expr_L[9].args...)



                    be = Symbol("be_for_$(field_trial)_$(field_test)_$(face_around_trial)_$(face_around_test)")
                    # statements depend on dof
                    
                    assignment = :($be[$dof_test, $dof_trial] += $v)
                    push!(dof_test_loop_body.args,assignment)
                end
            end
        end
    end

    for field_trial in 1:nfields_trial
        n_faces_around_trial = field_n_faces_around_trial[field_trial]
        for field_test in 1:nfields_test
            n_faces_around_test = field_n_faces_around_test[field_test]
            for face_around_trial in 1:n_faces_around_trial
                for face_around_test in 1:n_faces_around_test
                    dofs_test = Symbol("dofs_for_$(field_test)_$(face_around_test)_test")
                    dofs_trial = Symbol("dofs_for_$(field_trial)_$(face_around_trial)_trial")
                    be = Symbol("be_for_$(field_trial)_$(field_test)_$(face_around_trial)_$(face_around_test)")
                    assignment = :(GT.contribute!(alloc,$be,$dofs_test,$dofs_trial,$field_test,$field_trial))
                    push!(face_loop_body.args,assignment)
                end
            end
        end
    end
    block
end

function lowbit(a::Int)
    a & (-a) # trick to get the lowest bit of an integer.
end

function lca_deps(a::Int, b::Int, len::Int)
    lca = 0
    if a == b
        lca = a
    else
        diff = a ⊻ b
        lca = (lowbit(diff) - 1) & a
    end
    result = a & ~(lca)
    return result # bitmap a - lca
end

function binomial_tree_tabulation(expr_L, unroll_loops)
    result = Dict()
    used_vars_to_node = Dict()
    len_expr = length(expr_L)
    len_deps = trailing_zeros(len_expr)
    loop_bits = 0
    for (i, is_unrolled) in enumerate(unroll_loops)
        if !is_unrolled
            loop_bits += 2^(i-1)
        end
    end

    function binomial_tree_tabulation_impl(position, last_dep)
        function check_used(a) # TODO: many functions defined here. check performance
        end

        function check_used(a::Union{Array, Tuple})
            map(check_used, a)
        end
        function check_used(a::Expr)
            map(check_used, a.args)
        end
        function check_used(a::Symbol)
            if haskey(used_vars_to_node, a)
                deps = lca_deps(used_vars_to_node[a], position, len_deps)
                if  (deps & loop_bits) != 0   # &&  (position - used_vars_to_node[a]) & loop_bits != 0
                    result[a] = deps
                end
            end
        end

        function insert_variables!(a)
        end
        function insert_variables!(a::Union{Array, Tuple})
            map(insert_variables!, a)
        end
        function insert_variables!(a::Expr) 
            if a.head === :(=) && a.args[1] isa Symbol
                var = a.args[1]
                used_vars_to_node[var] = position
            end
            map(insert_variables!, a.args)
        end

        # step 1: Check if the variables used in this block are in the used var dict. 
        # For all variables used but exist, insert them into the result $r$. 
        # the mem alloc size is updated as the loop path from the Least Common Ancestor LCA(a, s[v]) to v. 

        check_used(expr_L[position + 1])
        # step 2: Traverse over the binomial tree, DFS in the left-most order and call the same step 2.  
        for i in len_deps:-1:(last_dep + 1)
            binomial_tree_tabulation_impl(position | (2 ^ (i - 1)), i)
        end
        # step 3: Insert all variables defined in this block
        insert_variables!(expr_L[position+1])
    end
    
    binomial_tree_tabulation_impl(0, 0)
    return result
end

# TODO: loops in lambda functions? we need to tabulate the computation of that


# input expr_L: expression in bitmap
# output: an expression (bitmap) array for each unrolled subset. 
# unrolled_expr_L has a length of 2^l where l is the number of vars
# each of them is a nested array, representing the dependency of each unrolled var
function unroll_expr_L(expr_L, bindings, unrolled_deps_and_length)
    # TODO: block extraction here. be careful with the dependencies. The simplified code may not depend on a loop var but it is implicitly dependent bacause of block extraction.
    unrolled_expr_L::Vector{Any} = [nothing for i in expr_L]
    status = [0 for _ in bindings] # status of unrolled variables
    len_deps = length(bindings)

    var_deps = Dict() # deps of unrolled variables. this is used to generate new names for each unrolled loop
    unrolled_var_deps = Dict()
    
    function replace_statement_deps_unroll(a::Symbol)
        if haskey(var_deps, a)
            suffix_status = map(x -> status[x], var_deps[a])
            suffix = map(x -> "_$(x)", suffix_status)
            result_str = string(a, suffix...)
            result = Symbol(result_str)
            unrolled_var_deps[result] = (var_deps[a], suffix_status)
            result
        else
            a
        end
    end
        
    function replace_statement_deps_unroll(a)
        a
    end
    function replace_statement_deps_unroll(a::Expr)
        Expr(a.head, map(replace_statement_deps_unroll, a.args)...)
    end

    function generate_similar_block_array(statement, unrolled_index)
        if length(unrolled_index) == 0
            replace_statement_deps_unroll(statement)
        else 
            new_unrolled_index = unrolled_index[2:end]
            last_dep = unrolled_index[1]
            dep = unrolled_deps_and_length[last_dep] # must be an unrolled var
            len_unroll = if dep[1] === nothing
                dep[2]
            else
                len_unroll = if dep[1] isa Tuple
                    result = dep[2]
                    for i in dep[1]
                        result = result[status[i]]
                    end
                    result
                else
                    index = status[dep[1]]
                    if index == 0
                        index = 1
                    end
                    dep[2][index]
                end
            end
            
            map(1:len_unroll) do x
                status[last_dep] = x
                result = generate_similar_block_array(statement, new_unrolled_index)
                status[last_dep] = 0
                result
            end

        end
    end

    function unroll_expr_L_impl!(position, last_dep, last_unrolled_index = ())
        # alloc array or use a single block
        is_unrolled = (last_dep > 0) && (unrolled_deps_and_length[last_dep] !== nothing)
        new_unrolled_index = is_unrolled ? (last_unrolled_index..., last_dep) : last_unrolled_index
        
        for statement in expr_L[position+1].args # update deps
            if statement.head === :(=) && statement.args[1] isa Symbol
                var = statement.args[1]
                var_deps[var] = new_unrolled_index
            end
        end
        unrolled_expr_L[position+1] = generate_similar_block_array(expr_L[position+1], new_unrolled_index)
        for i in len_deps:-1:(last_dep+1)
            unroll_expr_L_impl!(position | (2 ^ (i - 1)), i, new_unrolled_index)
        end
    end
    
    unroll_expr_L_impl!(0, 0, ())
    return unrolled_expr_L, unrolled_var_deps
end

function expr_L_protos(expr_L, bindings, unroll_loop_length)
    proto_block = Expr(:block)
    var_proto::Dict{Symbol, Symbol} = Dict()
    len_expr = length(expr_L)
    len_deps = trailing_zeros(len_expr)
    status = [1 for _ in bindings] # status of unrolled variables
    binding_index = Dict([(binding => i)  for (i, binding) in enumerate(bindings)])
    proto_var_count = 0

    function replace_vars_proto(a)
        a
    end

    function replace_vars_proto(a::Symbol)
        if haskey(var_proto, a)
            var_proto[a]
        elseif haskey(binding_index, a)
            index = binding_index[a]
            status[index]
        else
            a
        end
    end

    function replace_vars_proto(a::Expr)
        if a.head === :(=) && a.args[1] isa Symbol && !haskey(var_proto, a.args[1])
            var = a.args[1]
            proto_var_count += 1
            var_proto[var] = Symbol("proto_$proto_var_count")
        end
        Expr(a.head, map(replace_vars_proto, a.args)...)
    end

    function replace_vars_proto_unrolled!(statements, position)
        if statements isa Expr
            for statement in statements.args
                new_statement = replace_vars_proto(statement)
                push!(proto_block.args, new_statement) 
            end
        else
            next_unroll = 0
            new_position = position
            while (next_unroll == 0) || (unroll_loop_length[next_unroll] === nothing)
                next_unroll = trailing_zeros(new_position) + 1
                new_position -= 2 ^ (next_unroll - 1)
            end
            dep = unroll_loop_length[next_unroll]
            len_unroll = (dep[1] === nothing) ? dep[2] : dep[2][status[dep[1]]] # TODO: handle tuples
            @assert len_unroll == length(statements)
            for i in 1:len_unroll
                status[next_unroll] = i
                replace_vars_proto_unrolled!(statements[i], new_position)
            end
            status[next_unroll] = 1
        end
    end

    for position in 0:(len_expr - 1)
        statements = expr_L[position+1]
        replace_vars_proto_unrolled!(statements, position)
    end
    return proto_block, var_proto
end


function expr_replace_tabulates(expr_L, unrolled_var_deps, tabulated_variables, bindings, var_proto, unroll_loop_length, unroll_deps_length, max_lengths_symbol = :array_lengths)
    
    alloc_block = Expr(:block)
    len_expr = length(expr_L)
    len_deps = trailing_zeros(len_expr)
    status = [0 for i in 1:len_deps]

    for (k, v) in tabulated_variables
        deps = :()
        fill!(status, 0)
        if haskey(unrolled_var_deps, k)
            for (i, j) in zip(unrolled_var_deps[k]...)
                status[i] = j
            end
        end
        for i in len_deps:-1:1
            if v & (2 ^ (i-1)) != 0 && unroll_loop_length[i] === nothing
                # get exact length with unrolled deps
                dep = if unroll_deps_length[i] === nothing
                    :($max_lengths_symbol[$i])
                elseif status[unroll_deps_length[i]] == 0
                    :(reduce(max, $max_lengths_symbol[$i]))
                else
                    :($max_lengths_symbol[$i][$(status[unroll_deps_length[i]])]) # TODO: tuples?
                end
                push!(deps.args, dep)
            end
        end

        proto = var_proto[k]
        statement = :( $k = zeros(typeof($proto), ($deps)...))
        push!(alloc_block.args, statement)
    end

    function replace_vars(a)
        a
    end
    function replace_vars(a::Symbol)
        if haskey(tabulated_variables, a)
            result = Expr(:ref, a)
            v = tabulated_variables[a]
            for i in len_deps:-1:1
                if v & (2 ^ (i-1)) != 0 && unroll_loop_length[i] === nothing
                    push!(result.args, bindings[i])
                end
            end
            result
        else
            a
        end
    end
    function replace_vars(a::Expr)
        Expr(a.head, map(replace_vars, a.args)...)
    end

    function replace_vars(a::Union{Tuple, Array})
        map(replace_vars, a)
    end

    for i in 1:length(expr_L)
        expr_L[i] = replace_vars(expr_L[i])
    end
    alloc_block, expr_L
end


# TODO: include heuristic loop fusion?
function get_side_loops(position, current_dep, bindings, unroll_loops, unrolled_expr_L, unroll_loop_length,
        loop_range, unroll_deps_length, unroll_deps_length_symbol, current_unrolled_vars_status = [])
    result = Expr(:block)
    len = length(bindings)
    unrolled_vars = zeros(Int, length(unroll_loops))
    for (i, j) in current_unrolled_vars_status
        unrolled_vars[i] = j
    end
    function get_side_loops_impl!(pos, last_dep)
        n_statements = 0
        loop_body = Expr(:block)
        # select the correct block
        statements = unrolled_expr_L[pos+1]
        for i in 1:len
            if unrolled_vars[i] != 0
                statements = statements[unrolled_vars[i]]
            end
        end
        if unroll_loops[last_dep] == false
            pre_loop = loop_range[last_dep][1:end - 1]
            push!(loop_body.args, pre_loop...)

            loop_head = loop_range[last_dep][end]
            loop = Expr(:for,loop_head,loop_body)
            
            n_statements += length(statements.args)
            push!(loop_body.args, statements.args...)

            for i in len:-1:last_dep + 1
                loop_expr, n_statements_child = get_side_loops_impl!(pos | (2 ^ (i - 1)), i)
                if n_statements_child > 0
                    push!(loop_body.args, loop_expr)
                    n_statements += n_statements_child
                end
            end
            
            # TODO: this is ugly. use unroll_deps_length for these branches
            if unroll_deps_length[last_dep] !== nothing
                assignment = unroll_deps_length_symbol(unrolled_vars, last_dep) 
                loop = Expr(:block, assignment, loop)
            end
            return  loop, n_statements
        else
            loop_length = 0 
            loop_dep = unroll_loop_length[last_dep]
            if loop_dep[1] === nothing
                loop_length = loop_dep[2]
            elseif unrolled_vars[loop_dep[1]] == 0
                loop_length = reduce(max, loop_dep[2])
            else
                loop_length = loop_dep[2][unrolled_vars[loop_dep[1]]]
            end
            loop_var = loop_range[last_dep][end][1]
            for i in 1:loop_length
                unrolled_vars[last_dep] = i
                assignment = :($(loop_var) = $i)
                push!(loop_body.args, assignment)
                local_statements = statements[i]
                n_statements += length(local_statements.args)
                push!(loop_body.args, local_statements.args...)
                # update_ndofs
                for j in len:-1:(last_dep + 1)
                    loop_expr, n_statements_child = get_side_loops_impl!(pos | (2 ^ (j - 1)), j)
                    if n_statements_child > 0
                        push!(loop_body.args, loop_expr)
                        n_statements += n_statements_child
                    end
                end
            end
            unrolled_vars[last_dep] = 0

            return  loop_body, n_statements
        end
    end

    for i in len:-1:current_dep + 2
        loop_expr, n_statements = get_side_loops_impl!(position | (2 ^ (i - 1)), i)
        if n_statements > 0
            push!(result.args, loop_expr)
        end
    end
    return result
end




function generate_vector_assembly_loop_body_bitmap(field_n_faces_around, expr_L, dependencies, bindings)
    # expr_L: a list of block in the order of (global, face, point, field, face_around, dof)
    # dependencies: (space, quadrature, alloc)
    # bindings: (face, point, field, face_around, dof)
    (space, quadrature, alloc) = dependencies
    (face, point, field_symbol, face_around_symbol, dof) = bindings # TODO: replace field_symbol and face_around_symbol
    unroll_loops = (false, false, true, true, false)
    unroll_loop_length  = (nothing, nothing, (nothing, length(field_n_faces_around)), (3, field_n_faces_around), nothing)
    unroll_deps_length = (nothing, nothing, nothing, nothing, 3)
    unroll_deps_length_symbol = (status, last_dep) -> begin
                 if last_dep == 5  # dof trial
                    field = max(status[3], 1)
                    face_around = max(status[4], 1)
                    ndofs = Symbol("n_dofs_for_$(field)_$(face_around)")
                    assignment = :(ndofs_local = $ndofs)
                end  
    end

    v = last_symbol(expr_L)
    new_v = Symbol("vector_contribution_result")
    push!(expr_L[end].args, :($new_v = $v)) 
    v = new_v

    unrolled_expr_L, unrolled_var_deps = unroll_expr_L(expr_L, bindings, unroll_loop_length)
    tabulated_variables = binomial_tree_tabulation(unrolled_expr_L, unroll_loops)
    proto_block, var_proto = expr_L_protos(unrolled_expr_L, bindings, unroll_loop_length)
    expr_alloc, unrolled_expr_L = expr_replace_tabulates(unrolled_expr_L, unrolled_var_deps, tabulated_variables, bindings, var_proto, unroll_loop_length, unroll_deps_length)


    block = Expr(:block)
    nfields = length(field_n_faces_around)

    assignment = :(domain = GT.domain($quadrature))
    push!(block.args,assignment)

    assignment = :(T = eltype($alloc))
    push!(block.args, assignment)

    assignment = :(z = zero(T))
    push!(block.args, assignment)

    assignment = :(nfaces = GT.num_faces(domain))
    push!(block.args,assignment)

    assignment = :(face_npoints = GT.num_points_accessor($quadrature))
    push!(block.args,assignment)

    # Get single field spaces
    for field in 1:nfields
        var = Symbol("space_for_$(field)")
        expr = :(GT.field($space,$field))
        assignment = :($var = $expr)
        push!(block.args, assignment)
    end

    # Get face_dofs
    for field in 1:nfields
        space_field = Symbol("space_for_$(field)")
        var_field = Symbol("face_dofs_for_$(field)")
        expr = :(GT.dofs_accessor($space_field,domain))
        assignment = :($var_field = $expr)
        push!(block.args, assignment)
    end

    max_num_dofs = :(())
    for field in 1:nfields
        space_field = Symbol("space_for_$(field)")
        push!(max_num_dofs.args, :(GT.max_num_reference_dofs($space_field)))
    end

    # Allocate face vectors
    for field in 1:nfields
        n_faces_around = field_n_faces_around[field]
        for face_around in 1:n_faces_around
            space_field = Symbol("space_for_$(field)")
            var_str = "be_for_$(field)_$(face_around)"
            var_face_around = Symbol(var_str)
            expr = :(alloc_zeros($var_str, T,GT.max_num_reference_dofs($space_field)))
            assignment = :($var_face_around = $expr)
            push!(block.args, assignment)
        end
    end

    # max_length_deps
    assignment = :(array_lengths = (nfaces, 
                GT.max_num_points($quadrature), 
                $nfields, 
                $field_n_faces_around, 
                $max_num_dofs)
    )
    
    loop_range = ((:($face = 1:nfaces),), 
                    (:(npoints = face_npoints($face)), :($point = 1:npoints)), 
                    ((field_symbol, nfields), ),
                    ((face_around_symbol, nothing), ),
                    (:($dof = 1:ndofs_local), ),  # this should depend on face 
    )

    function get_side_loops_local(position, current_dep, current_unrolled_vars_status = [])
        get_side_loops(position, current_dep, bindings, unroll_loops, unrolled_expr_L, unroll_loop_length, loop_range, unroll_deps_length, unroll_deps_length_symbol, current_unrolled_vars_status)
    end

    
    push!(block.args, assignment)
    push!(block.args, proto_block.args...)
    push!(block.args, expr_alloc.args...)

    # statements without deps
    push!(block.args, unrolled_expr_L[1].args...)

    loops = get_side_loops_local(0, 0)
    push!(block.args, loops.args...)

    # Face loop
    face_loop_head = :($face = 1:nfaces)
    face_loop_body = Expr(:block)
    face_loop = Expr(:for,face_loop_head,face_loop_body)
    push!(block.args,face_loop)

    assignment = :(npoints = face_npoints($face))
    push!(face_loop_body.args,assignment)

    # statements depend on face
    push!(face_loop_body.args, unrolled_expr_L[2].args...)
    loops = get_side_loops_local(1, 1)
    push!(face_loop_body.args, loops.args...)

    for field in 1:nfields
        n_faces_around = field_n_faces_around[field]
        for face_around in 1:n_faces_around
            dofs = Symbol("dofs_for_$(field)_$(face_around)")
            accessor = Symbol("face_dofs_for_$(field)")
            assignment = :($dofs = $(accessor)($face,$face_around))
            push!(face_loop_body.args,assignment)
        end
    end

    for field in 1:nfields
        n_faces_around = field_n_faces_around[field]
        for face_around in 1:n_faces_around
            dofs = Symbol("dofs_for_$(field)_$(face_around)")
            ndofs = Symbol("n_dofs_for_$(field)_$face_around")
            assignment = :($ndofs = length($(dofs)))
            push!(face_loop_body.args,assignment)
        end
    end

    for field in 1:nfields
        n_faces_around = field_n_faces_around[field]
        for face_around in 1:n_faces_around
            be = Symbol("be_for_$(field)_$(face_around)")
            assignment = :(fill!($be,z))
            push!(face_loop_body.args,assignment)
        end
    end

    # Point loop
    point_loop_head = :($point = 1:npoints)
    point_loop_body = Expr(:block)
    point_loop = Expr(:for,point_loop_head,point_loop_body)
    push!(face_loop_body.args,point_loop)

    # statements depend on point
    ndofs = Symbol("n_dofs_for_1_1")
    assignment = :(ndofs_local = $ndofs)
    push!(point_loop_body.args, assignment)

    # statements depend on point
    push!(point_loop_body.args, unrolled_expr_L[4].args...)
    loops = get_side_loops_local(3, 2)
    push!(point_loop_body.args, loops.args...)


    for field in 1:nfields
        n_faces_around = field_n_faces_around[field]
        # statements depend on field
        push!(point_loop_body.args, :($field_symbol = $field))
        push!(point_loop_body.args, unrolled_expr_L[8][field].args...)
        unrolled_var_status = ((3, field),)
        loops = get_side_loops_local(7, 3, unrolled_var_status)
        push!(point_loop_body.args, loops.args...)

        for face_around in 1:n_faces_around
            # statements depend on face around
            push!(point_loop_body.args, :($face_around_symbol = $face_around))
            push!(point_loop_body.args, unrolled_expr_L[16][field][face_around].args...)
            unrolled_var_status = ((3, field), (4, face_around))
            loops = get_side_loops_local(15, 4, unrolled_var_status)
            push!(point_loop_body.args, loops.args...)

            # DOF loop
            ndofs = Symbol("n_dofs_for_$(field)_$(face_around)")
            dof_loop_head = :($dof = 1:$ndofs)
            dof_loop_body = Expr(:block)
            dof_loop = Expr(:for,dof_loop_head,dof_loop_body)
            push!(point_loop_body.args,dof_loop)
            be = Symbol("be_for_$(field)_$(face_around)")
            # statements depend on dof
            suffix = map(x -> "_$x", [field, face_around])
            v_string = string(v, suffix...)
            push!(dof_loop_body.args, unrolled_expr_L[32][field][face_around].args...)
            assignment = :($be[$dof] += $(Symbol(v_string)))
            push!(dof_loop_body.args,assignment)
        end
    end

    for field in 1:nfields
        n_faces_around = field_n_faces_around[field]
        for face_around in 1:n_faces_around
            dofs = Symbol("dofs_for_$(field)_$(face_around)")
            be = Symbol("be_for_$(field)_$(face_around)")
            assignment = :(GT.contribute!(alloc,$be,$dofs,$field))
            push!(face_loop_body.args,assignment)
        end
    end

    block
end




function generate_matrix_assembly_loop_body_bitmap(field_n_faces_around_trial, field_n_faces_around_test, expr_L, dependencies, bindings)
    # expr_L: a list of block in the order of (global, face, point, field_trial, field_test, face_around_trial, face_around_test, dof_trial, dof_test)
    # dependencies: (space_trial, space_test, quadrature, alloc)
    # bindings: (face, point, field_trial, field_test, face_around_trial, face_around_test, dof_trial, dof_test)
    (space_trial, space_test, quadrature, alloc) = dependencies
    (face, point, field_trial_symbol, field_test_symbol, face_around_trial_symbol, face_around_test_symbol, dof_trial, dof_test) = bindings 
    unroll_loops = (false, false, true, true, true, true, false, false)
    unroll_loop_length  = (nothing, nothing, (nothing, length(field_n_faces_around_trial)), (nothing, length(field_n_faces_around_test)), (3, field_n_faces_around_trial), (4, field_n_faces_around_test), nothing, nothing)
    unroll_deps_length = [nothing, nothing, nothing, nothing, nothing, nothing, 3, 4]
    unroll_deps_length_symbol = (status, last_dep) -> begin
                 if last_dep == 7  # dof trial
                    field = max(status[3], 1)
                    face_around = max(status[5], 1)
                    ndofs = Symbol("n_dofs_for_$(field)_$(face_around)_trial")
                    assignment = :(ndofs_trial_local = $ndofs)
                elseif last_dep == 8 # dof test
                    field = max(status[4], 1)
                    face_around = max(status[6], 1)
                    ndofs = Symbol("n_dofs_for_$(field)_$(face_around)_test")
                    assignment = :(ndofs_test_local = $ndofs)
                end  
    end

    v = last_symbol(expr_L)
    new_v = Symbol("matrix_contribution_result")
    push!(expr_L[end].args, :($new_v = $v))
    v = new_v

    unrolled_expr_L, unrolled_var_deps = unroll_expr_L(expr_L, bindings, unroll_loop_length)
    tabulated_variables = binomial_tree_tabulation(unrolled_expr_L, unroll_loops)
    proto_block, var_proto = expr_L_protos(unrolled_expr_L, bindings, unroll_loop_length)
    expr_alloc, unrolled_expr_L = expr_replace_tabulates(unrolled_expr_L, unrolled_var_deps, tabulated_variables, bindings, var_proto, unroll_loop_length, unroll_deps_length)


    block = Expr(:block)
    nfields_trial = length(field_n_faces_around_trial)
    nfields_test = length(field_n_faces_around_test)


    assignment = :(domain = GT.domain($quadrature))
    push!(block.args,assignment)

    assignment = :(T = eltype($alloc))
    push!(block.args, assignment)

    assignment = :(z = zero(T))
    push!(block.args, assignment)
    
    assignment = :(nfaces = GT.num_faces(domain))
    push!(block.args,assignment)

    assignment = :(face_npoints = GT.num_points_accessor($quadrature))
    push!(block.args,assignment)

    # Get single field spaces 
    for field in 1:nfields_trial
        var = Symbol("space_for_$(field)_trial")
        expr = :(GT.field($space_trial,$field))
        assignment = :($var = $expr)
        push!(block.args, assignment)
    end

    for field in 1:nfields_test
        var = Symbol("space_for_$(field)_test")
        expr = :(GT.field($space_test,$field))
        assignment = :($var = $expr)
        push!(block.args, assignment)
    end

    # Get face_dofs
    for field in 1:nfields_trial
        space_field = Symbol("space_for_$(field)_trial")
        var_field = Symbol("face_dofs_for_$(field)_trial")
        expr = :(GT.dofs_accessor($space_field,domain))
        assignment = :($var_field = $expr)
        push!(block.args, assignment)
    end

    for field in 1:nfields_test
        space_field = Symbol("space_for_$(field)_test")
        var_field = Symbol("face_dofs_for_$(field)_test")
        expr = :(GT.dofs_accessor($space_field,domain))
        assignment = :($var_field = $expr)
        push!(block.args, assignment)
    end

    
    max_num_dofs_trial = :(())
    max_num_dofs_test = :(())
    # Allocate face vectors

    for field_trial in 1:nfields_trial 
        space_field_trial = Symbol("space_for_$(field_trial)_trial")
        push!(max_num_dofs_trial.args, :(GT.max_num_reference_dofs($space_field_trial)))
    end
    for field_test in 1:nfields_test 
        space_field_test = Symbol("space_for_$(field_test)_test")
        push!(max_num_dofs_test.args, :(GT.max_num_reference_dofs($space_field_test)))
    end

    for field_trial in 1:nfields_trial # TODO: insert it into expr_L
        n_faces_around_trial = field_n_faces_around_trial[field_trial]
        space_field_trial = Symbol("space_for_$(field_trial)_trial")
        for field_test in 1:nfields_test
            n_faces_around_test = field_n_faces_around_test[field_test]
            space_field_test = Symbol("space_for_$(field_test)_test")
            for face_around_trial in 1:n_faces_around_trial
                for face_around_test in 1:n_faces_around_test
                    var_str = "be_for_$(field_trial)_$(field_test)_$(face_around_trial)_$(face_around_test)"
                    var_face_around = Symbol(var_str)
                    expr = :(alloc_zeros($var_str,T,GT.max_num_reference_dofs($space_field_test),GT.max_num_reference_dofs($space_field_trial)))
                    assignment = :($var_face_around = $expr)
                    push!(block.args, assignment)
                end
            end
        end
    end
    
    # max_length_deps

    assignment = :(array_lengths = (nfaces, 
                GT.max_num_points($quadrature), 
                $nfields_trial, 
                $nfields_test, 
                $field_n_faces_around_trial, 
                $field_n_faces_around_test, 
                $max_num_dofs_trial,
                $max_num_dofs_test)
    )
    
    loop_range = ((:($face = 1:nfaces),), 
                    (:(npoints = face_npoints($face)), :($point = 1:npoints)), 
                    ((field_trial_symbol, nfields_trial), ),
                    ((field_test_symbol, nfields_test), ),
                    ((face_around_trial_symbol, nothing), ),
                    ((face_around_test_symbol, nothing), ),
                    (:($dof_trial = 1:ndofs_trial_local), ),  # this should depend on face 
                    (:($dof_test = 1:ndofs_test_local), ),    # this should depend on face
    )

    function get_side_loops_local(position, current_dep, current_unrolled_vars_status = [])
        get_side_loops(position, current_dep, bindings, unroll_loops, unrolled_expr_L, unroll_loop_length, loop_range, unroll_deps_length, unroll_deps_length_symbol, current_unrolled_vars_status)
    end

    push!(block.args, assignment)

    push!(block.args, proto_block.args...)
    push!(block.args, expr_alloc.args...)


    # statements without deps
    push!(block.args, unrolled_expr_L[1].args...)

    loops = get_side_loops_local(0, 0)
    push!(block.args, loops.args...)

    # Face loop
    face_loop_head = :($face = 1:nfaces)
    face_loop_body = Expr(:block)
    face_loop = Expr(:for,face_loop_head,face_loop_body)
    push!(block.args,face_loop)

    assignment = :(npoints = face_npoints($face))
    push!(face_loop_body.args,assignment)

    # statements depend on face
    push!(face_loop_body.args, unrolled_expr_L[2].args...)
    loops = get_side_loops_local(1, 1)
    push!(face_loop_body.args, loops.args...)

    # trial
    for field in 1:nfields_trial
        n_faces_around = field_n_faces_around_trial[field]
        for face_around in 1:n_faces_around
            dofs = Symbol("dofs_for_$(field)_$(face_around)_trial")
            accessor = Symbol("face_dofs_for_$(field)_trial")
            assignment = :($dofs = $(accessor)($face,$face_around))
            push!(face_loop_body.args,assignment)
        end
    end

    # test
    for field in 1:nfields_test
        n_faces_around = field_n_faces_around_test[field]
        for face_around in 1:n_faces_around
            dofs = Symbol("dofs_for_$(field)_$(face_around)_test")
            accessor = Symbol("face_dofs_for_$(field)_test")
            assignment = :($dofs = $(accessor)($face,$face_around))
            push!(face_loop_body.args,assignment)
        end
    end

    # trial
    for field in 1:nfields_trial
        n_faces_around = field_n_faces_around_trial[field]
        for face_around in 1:n_faces_around
            dofs = Symbol("dofs_for_$(field)_$(face_around)_trial")
            ndofs = Symbol("n_dofs_for_$(field)_$(face_around)_trial")
            assignment = :($ndofs = length($(dofs)))
            push!(face_loop_body.args,assignment)
        end
    end

    # test
    for field in 1:nfields_test
        n_faces_around = field_n_faces_around_test[field]
        for face_around in 1:n_faces_around
            dofs = Symbol("dofs_for_$(field)_$(face_around)_test")
            ndofs = Symbol("n_dofs_for_$(field)_$(face_around)_test")
            assignment = :($ndofs = length($(dofs)))
            push!(face_loop_body.args,assignment)
        end
    end


    for field_trial in 1:nfields_trial
        n_faces_around_trial = field_n_faces_around_trial[field_trial]
        for field_test in 1:nfields_test
            n_faces_around_test = field_n_faces_around_test[field_test]
            for face_around_trial in 1:n_faces_around_trial
                for face_around_test in 1:n_faces_around_test
                    var_str = "be_for_$(field_trial)_$(field_test)_$(face_around_trial)_$(face_around_test)"
                    be = Symbol(var_str)
                    assignment = :(fill!($be,z))
                    push!(face_loop_body.args,assignment)
                end
            end
        end
    end


    # Point loop
    point_loop_head = :($point = 1:npoints)
    point_loop_body = Expr(:block)
    point_loop = Expr(:for,point_loop_head,point_loop_body)
    push!(face_loop_body.args,point_loop)


    ndofs = Symbol("n_dofs_for_1_1_trial")
    assignment = :(ndofs_trial_local = $ndofs)
    push!(point_loop_body.args, assignment)

    ndofs = Symbol("n_dofs_for_1_1_test")
    assignment = :(ndofs_test_local = $ndofs)
    push!(point_loop_body.args, assignment)



    # statements depend on point
    push!(point_loop_body.args, unrolled_expr_L[4].args...)
    loops = get_side_loops_local(3, 2)
    push!(point_loop_body.args, loops.args...)


    for field_trial in 1:nfields_trial
        n_faces_around_trial = field_n_faces_around_trial[field_trial]
        # statements depend on field
        push!(point_loop_body.args, :($field_trial_symbol = $field_trial))
        push!(point_loop_body.args, unrolled_expr_L[8][field_trial].args...)
        unrolled_var_status = ((3, field_trial),)
        loops = get_side_loops_local(7, 3, unrolled_var_status)
        push!(point_loop_body.args, loops.args...)


        for field_test in 1:nfields_test
            n_faces_around_test = field_n_faces_around_test[field_test]
            # statements depend on field
            push!(point_loop_body.args, :($field_test_symbol = $field_test))
            push!(point_loop_body.args, unrolled_expr_L[16][field_trial][field_test].args...)
            unrolled_var_status = ((3, field_trial), (4, field_test))
            loops = get_side_loops_local(15, 4, unrolled_var_status)
            push!(point_loop_body.args, loops.args...)

            for face_around_trial in 1:n_faces_around_trial
                # statements depend on face around
                push!(point_loop_body.args, :($face_around_trial_symbol = $face_around_trial))
                push!(point_loop_body.args, unrolled_expr_L[32][field_trial][field_test][face_around_trial].args...)
                unrolled_var_status = ((3, field_trial), (4, field_test), (5, face_around_trial))
                loops = get_side_loops_local(31, 5, unrolled_var_status)
                push!(point_loop_body.args, loops.args...)
                for face_around_test in 1:n_faces_around_test
                    # statements depend on face around
                    push!(point_loop_body.args, :($face_around_test_symbol = $face_around_test))
                    push!(point_loop_body.args, unrolled_expr_L[64][field_trial][field_test][face_around_trial][face_around_test].args...)
                    unrolled_var_status = ((3, field_trial), (4, field_test), (5, face_around_trial), (6, face_around_test))
                    loops = get_side_loops_local(63, 6, unrolled_var_status)
                    push!(point_loop_body.args, loops.args...)

                    # DOF loop
                    # do just after the field and face_around computation, or otherwise we need to replace the variable names in expr_L 
                    ndofs_trial = Symbol("n_dofs_for_$(field_trial)_$(face_around_trial)_trial")
                    dof_trial_loop_head = :($dof_trial = 1:$ndofs_trial)
                    dof_trial_loop_body = Expr(:block)
                    dof_trial_loop = Expr(:for,dof_trial_loop_head,dof_trial_loop_body)
                    push!(point_loop_body.args,dof_trial_loop)
                    
                    push!(dof_trial_loop_body.args, unrolled_expr_L[128][field_trial][field_test][face_around_trial][face_around_test].args...)
                    loops = get_side_loops_local(127, 7, unrolled_var_status)
                    push!(dof_trial_loop_body.args, loops.args...)

                    ndofs_test = Symbol("n_dofs_for_$(field_test)_$(face_around_test)_test")
                    dof_test_loop_head = :($dof_test = 1:$ndofs_test)
                    dof_test_loop_body = Expr(:block)
                    dof_test_loop = Expr(:for,dof_test_loop_head,dof_test_loop_body)
                    push!(dof_trial_loop_body.args,dof_test_loop)
                    
                    push!(dof_test_loop_body.args, unrolled_expr_L[256][field_trial][field_test][face_around_trial][face_around_test].args...) # assuming that the last block is not empty
                    # loops = get_side_loops(255, 9)
                    # push!(dof_test_loop_body.args, loops.args...)


                    be = Symbol("be_for_$(field_trial)_$(field_test)_$(face_around_trial)_$(face_around_test)")
                    # statements depend on dof
                    suffix = map(x -> "_$x", [field_trial, field_test, face_around_trial, face_around_test])
                    v_string = string(v, suffix...)
                    assignment = :($be[$dof_test, $dof_trial] += $(Symbol(v_string)))
                    push!(dof_test_loop_body.args,assignment)
                end
            end
        end
    end

    for field_trial in 1:nfields_trial
        n_faces_around_trial = field_n_faces_around_trial[field_trial]
        for field_test in 1:nfields_test
            n_faces_around_test = field_n_faces_around_test[field_test]
            for face_around_trial in 1:n_faces_around_trial
                for face_around_test in 1:n_faces_around_test
                    dofs_test = Symbol("dofs_for_$(field_test)_$(face_around_test)_test")
                    dofs_trial = Symbol("dofs_for_$(field_trial)_$(face_around_trial)_trial")
                    be = Symbol("be_for_$(field_trial)_$(field_test)_$(face_around_trial)_$(face_around_test)")
                    assignment = :(GT.contribute!(alloc,$be,$dofs_test,$dofs_trial,$field_test,$field_trial))
                    push!(face_loop_body.args,assignment)
                end
            end
        end
    end
    block
end


#for op in (:+,:-,:*,:/,:\,:^)
#    @eval begin
#        (Base.$op)(a::AbstractQuantity,b::AbstractQuantity) = call(Base.$op,a,b)
#        (Base.$op)(a::Number,b::AbstractQuantity) = call(Base.$op,GT.compile_constant_quantity(a),b)
#        (Base.$op)(a::AbstractQuantity,b::Number) = call(Base.$op,a,GT.compile_constant_quantity(b))
#    end
#end

# Operations

# Base

#for op in (:+,:-,:sqrt,:abs,:abs2,:real,:imag,:conj,:transpose,:adjoint,:*,:/,:\,:^,:getindex)
#    @eval begin
#        function get_symbol!(index,val::typeof(Base.$op),name="";prefix=gensym)
#            $( Expr(:quote,op) )
#        end
#    end
#end

for op in (:+,:-,:sqrt,:abs,:abs2,:real,:imag,:conj,:transpose,:adjoint)
  @eval begin
      (Base.$op)(a::AbstractQuantity) = call(Base.$op,a)
  end
end

function Base.getindex(a::AbstractQuantity,b::AbstractQuantity)
    quantity() do opts
        term_a = term(a,opts)
        term_b = term(b,opts)
        RefTerm(term_a,term_b)
    end
end

 function Base.getindex(a::AbstractQuantity,b::Integer)
    f = uniform_quantity(b;is_compile_constant=true)
    a[f]
 end

for op in (:+,:-,:*,:/,:\,:^)
  @eval begin
      (Base.$op)(a::AbstractQuantity,b::AbstractQuantity) = call(Base.$op,a,b)
      #(Base.$op)(a::Number,b::AbstractQuantity) = call(Base.$op,GT.uniform_quantity(a),b)
      #(Base.$op)(a::AbstractQuantity,b::Number) = call(Base.$op,a,GT.uniform_quantity(b))
  end
end

for op in (:*,:/,:\,:^)
    @eval begin
        (Base.$op)(a::Number,b::AbstractQuantity) = call(Base.$op,GT.uniform_quantity(a),b)
        (Base.$op)(a::AbstractQuantity,b::Number) = call(Base.$op,a,GT.uniform_quantity(b))
    end
end


# LinearAlgebra

for op in (:inv,:det,:norm,:tr)
  @eval begin
      #function get_symbol!(index,val::typeof(LinearAlgebra.$op),name="";prefix=gensym)
      #    $( Expr(:quote,op) )
      #end
    (LinearAlgebra.$op)(a::AbstractQuantity) = call(LinearAlgebra.$op,a)
  end
end

for op in (:dot,:cross)
  @eval begin
      #function get_symbol!(index,val::typeof(LinearAlgebra.$op),name="";prefix=gensym)
      #    $( Expr(:quote,op) )
      #end
      (LinearAlgebra.$op)(a::AbstractQuantity,b::AbstractQuantity) = call(LinearAlgebra.$op,a,b)
    #   (LinearAlgebra.$op)(a::Number,b::AbstractQuantity) = call(LinearAlgebra.$op,GT.uniform_quantity(a),b)
    #   (LinearAlgebra.$op)(a::AbstractQuantity,b::Number) = call(LinearAlgebra.$op,a,GT.uniform_quantity(b))
  end
end


for op in (:(ForwardDiff.gradient),:(ForwardDiff.jacobian),:(ForwardDiff.hessian))
    @eval begin
        ($op)(a::AbstractQuantity,b::AbstractQuantity) = call($op,a,b)
    end
end


# a copy of the older version. need to be compared with statements over terms
function topological_sort_bitmap(expr,deps)
    temporary = gensym()
    expr_L = [Expr(:block) for _ in 1:2^length(deps)]
    marks = Dict{UInt,Any}()
    marks_deps = Dict{UInt,Int}()
    variable_count = 0

    function visit(expr_n::Union{Symbol, Function, Number, Nothing, Val})
        id_n = hash(expr_n)
        marks[id_n] = expr_n
        i = findfirst(e->expr_n===e,deps)
        j = i === nothing ? 0 : (1<<(Int(i)-1))
        marks_deps[id_n] = j
        expr_n , j
    end

    function setup_expr_default(expr_n)
        args = expr_n.args
        r = map(visit,args)
        args_var = map(first,r)
        j = reduce(|, map(last,r))
        expr_n_new = Expr(expr_n.head,args_var...)
        expr_n_new, j
    end
    function setup_expr_call(expr_n)
        args = expr_n.args
        r = map(visit,args)
        args_var = map(first,r)
        j = reduce(|, map(last,r))
        expr_n_new = Expr(expr_n.head,args_var...)
        expr_n_new, j
    end
    function setup_expr_lambda(expr_n) # do topological sort with deps & new deps, get a list of exprs and merge
        body = expr_n.args[2]
        args = (expr_n.args[1] isa Symbol) ? (expr_n.args[1], ) : expr_n.args[1].args
        new_deps = (deps..., args...)
        body_sorted = topological_sort_bitmap(body, new_deps) # TODO: is it efficient?
        # TODO: last var if it is a constant, or replace the function calls as a constant

        for i in 1:2^length(deps) # hoisting
            push!(expr_L[i].args, body_sorted[i].args...)
        end

        # compute dependencies
        j = reduce(|, filter(i -> length(body_sorted[i]) > 0, 1:2^length(new_deps) ))
        j &= 2^length(deps) - 1 # remove lambda arg deps

        merged_body_sorted = vcat(map(x -> x.args, body_sorted[2^length(deps) + 1: 2^length(new_deps)])...)
        expr_n_new = Expr(expr_n.head,expr_n.args[1],merged_body_sorted)
        expr_n_new, j
    end
    function visit(expr_n)
        id_n = hash(expr_n)
        if haskey(marks,id_n)
            if marks[id_n] !== temporary
                return marks[id_n], marks_deps[id_n]
            else
                error("Graph has cycles! This is not possible.")
            end
        end
        marks[id_n] = temporary
        variable_count += 1
        var = Symbol("var_$variable_count")
        if isa(expr_n,Expr)
            if expr_n.head === :call # || expr_n.head === :if  # || expr_n.head === :ref
                expr_n_new, j = setup_expr_call(expr_n)
            elseif expr_n.head === :(->)
                expr_n_new, j = setup_expr_lambda(expr_n)
            else
                expr_n_new, j = setup_expr_default(expr_n)
            end
            assignment = :($var = $expr_n_new)
        else
            j = 0
            assignment = :($var = $expr_n)
        end
        marks[id_n] = var
        marks_deps[id_n] = j
        push!(expr_L[j+1].args,assignment)
        var, j
    end
    visit(expr)
    expr_L
end

function topological_sort(expr,deps)
    temporary = gensym()
    expr_L = [Expr(:block) for _ in 0:length(deps)]
    marks = Dict{UInt,Any}()
    marks_deps = Dict{UInt,Int}()
    function visit(expr_n::Union{Symbol, Function, Number, Nothing, Val})
        id_n = hash(expr_n)
        marks[id_n] = expr_n
        i = findfirst(e->expr_n===e,deps)
        j = i === nothing ? 0 : Int(i)
        marks_deps[id_n] = j
        expr_n , j
    end

    function setup_expr_default(expr_n)
        args = expr_n.args
        r = map(visit,args)
        args_var = map(first,r)
        j = maximum(map(last,r))
        expr_n_new = Expr(expr_n.head,args_var...)
        expr_n_new, j
    end
    function setup_expr_call(expr_n)
        args = expr_n.args
        r = map(visit,args)
        args_var = map(first,r)
        j = maximum(map(last,r))
        expr_n_new = Expr(expr_n.head,args_var...)
        expr_n_new, j
    end
    function setup_expr_lambda(expr_n)
        body = expr_n.args[2]
        body_sorted = topological_sort(body,())[1]
        j = length(deps)
        expr_n_new = Expr(expr_n.head,expr_n.args[1],body_sorted)
        expr_n_new, j
    end
    function visit(expr_n)
        id_n = hash(expr_n)
        if haskey(marks,id_n)
            if marks[id_n] !== temporary
                return marks[id_n], marks_deps[id_n]
            else
                error("Graph has cycles! This is not possible.")
            end
        end
        marks[id_n] = temporary
        var = gensym()
        if isa(expr_n,Expr)
            if expr_n.head === :call # || expr_n.head === :ref
                expr_n_new, j = setup_expr_call(expr_n)
            elseif expr_n.head === :(->)
                expr_n_new, j = setup_expr_lambda(expr_n)
            else
                expr_n_new, j = setup_expr_default(expr_n)
            end
            assignment = :($var = $expr_n_new)
        else
            j = 0
            assignment = :($var = $expr_n)
        end
        marks[id_n] = var
        marks_deps[id_n] = j
        push!(expr_L[j+1].args,assignment)
        var, j
    end
    visit(expr)
    expr_L
end


function statements_expr(node)
    # @assert lambda_args_once(node) # keep hash_scope but make an error check. In some cases it will be a bug if we have a lambda function arg name more than once in the term
    root = :root
    scope_level = Dict{Symbol,Int}()
    scope_level[root] = 0
    hash_scope = Dict{UInt,Symbol}()
    scope_rank = Dict{Symbol,Int}()
    scope_rank[root] = 0
    hash_expr = Dict{UInt, Any}()
    scope_block = Dict(root=>Expr(:block))
    function visit(node, depth = 0) # if the terms here are immutable (no terms "appended" after construction) then we can assume there is no graph cycle
        hash = Base.hash(node)
        if haskey(hash_scope, hash)
            return [hash_scope[hash]]
        end
        if node isa Symbol || node isa Number || node isa Function || node isa Nothing || node isa Val
            scopes = [root]
            hash_expr[hash] = node
            hash_scope[hash] = root 
            return scopes
        elseif node isa Expr
            if  node.head === :call || node.head ===  :ref || node.head ===  :block
                
                # args = (node.head === :call) ? (node.callee, node.args...) : (node.array, node.index)
                args = filter(x->!(x isa LineNumberNode), node.args)
                scopes_nested = map(args) do arg
                    visit(arg, depth)
                end
                scopes = reduce(vcat,scopes_nested) |> unique
                scope = argmax(x->scope_level[x],scopes)
                hash_scope[hash] = scope
                rank = 1 + scope_rank[scope]
                scope_rank[scope] = rank

                args_var = map(args) do arg
                    arg_hash = Base.hash(arg)
                    hash_expr[arg_hash]
                end
                
                if node.head === :block && length(args_var) === 1
                    hash_expr[hash] = args_var[1]
                else
                    var = Symbol("var_$(scope)_$(rank)")
                    hash_expr[hash] = var
                    expr = Expr(node.head,args_var...)
                    assignment = :($var = $expr)
                    block = scope_block[scope]
                    push!(block.args, assignment)
                end
                return scopes
            elseif node.head === :(->)
                body = node.args[2]
                args = (node.args[1] isa Expr) ? node.args[1].args : (node.args[1], )
                block = Expr(:block)

                scope = (length(args) > 0) ? args[1] : gensym()
                scope_level[scope] = depth + 1
                scope_rank[scope] = 0                
                scope_block[scope] = block
    
                for arg in args 
                    arg_hash = Base.hash(arg)
                    hash_expr[arg_hash] = arg
                    hash_scope[arg_hash] = scope
                end
                
                scopes = visit(body, depth + 1)
                if length(block.args) == 0 # at least 1 statement
                    arg_hash = Base.hash(body)
                    arg_var = hash_expr[arg_hash]
                    push!(block.args, arg_var)
                    # push!(block.args, :($(gensym()) = $arg_var))
                end

                scopes = setdiff(scopes,[scope])
                scope = argmax(scope->scope_level[scope],scopes)
                hash_scope[hash] = scope
                rank = 1 + scope_rank[scope]
                scope_rank[scope] = rank
    
                var = Symbol("var_$(scope)_$(rank)")
                hash_expr[hash] = var
    
                expr = Expr(:(->), node.args[1], block) #  LambdaTerm(block, node.args, x -> prototype(node))
    
                assignment = :($var = $expr) #  StatementTerm(var_term, expr, prototype(node))
                block = scope_block[scope]
                push!(block.args, assignment)
                return scopes
            else
                @show node
                dump(node)
                @show typeof(node)
                error("A")
            end
        else
            @show node
            @show typeof(node)
            error("B")
        end


    end

    visit(node)
    scope_block[root]
end


# TODO: just a backup. maybe not needed any more.
function statements_expr_with_loops(node)
    # @assert lambda_args_once(node) # keep hash_scope but make an error check. In some cases it will be a bug if we have a lambda function arg name more than once in the term
    root = :root
    scope_level = Dict{Symbol,Int}()
    scope_level[root] = 0
    hash_scope = Dict{UInt,Symbol}()
    scope_rank = Dict{Symbol,Int}()
    scope_rank[root] = 0
    hash_expr = Dict{UInt, Any}()
    scope_block = Dict(root=>Expr(:block))
    function visit(node, depth = 0, last_scope = :root; has_value = true) # if the terms here are immutable (no terms "appended" after construction) then we can assume there is no graph cycle
        hash = Base.hash(node)
        if haskey(hash_scope, hash)
            return [hash_scope[hash]]
        end
        if node isa Symbol || node isa Number || node isa Function || node isa Nothing || node isa Val
            scopes = [root]
            hash_expr[hash] = node
            hash_scope[hash] = root 
            return scopes
        elseif node isa Expr
            if  node.head === :call || node.head ===  :ref || node.head ===  :block || node.head === :tuple
                
                # args = (node.head === :call) ? (node.callee, node.args...) : (node.array, node.index)
                args = filter(x->!(x isa LineNumberNode), node.args)
                scopes_nested = map(args) do arg
                    visit(arg, depth, last_scope)
                end
                scopes = reduce(vcat,scopes_nested) |> unique
                scope = argmax(x->scope_level[x],scopes)
                hash_scope[hash] = scope
                rank = 1 + scope_rank[scope]
                scope_rank[scope] = rank

                args_var = map(args) do arg
                    arg_hash = Base.hash(arg)
                    hash_expr[arg_hash]
                end
                
                if node.head === :block && length(args_var) == 1
                    hash_expr[hash] = args_var[1]
                elseif has_value
                    var = Symbol("var_$(scope)_$(rank)")
                    hash_expr[hash] = var
                    expr = Expr(node.head,args_var...)
                    assignment = :($var = $expr)
                    block = scope_block[scope]
                    push!(block.args, assignment)
                end
                return scopes
            elseif node.head === :(->)
                body = node.args[2]
                args = (node.args[1] isa Expr) ? node.args[1].args : (node.args[1], )
                block = Expr(:block)

                scope = (length(args) > 0) ? args[1] : gensym()
                scope_level[scope] = depth + 1
                scope_rank[scope] = 0                
                scope_block[scope] = block
    
                for arg in args 
                    arg_hash = Base.hash(arg)
                    hash_expr[arg_hash] = arg
                    hash_scope[arg_hash] = scope
                end
                
                scopes = visit(body, depth + 1, scope)
                if length(block.args) == 0 # at least 1 statement
                    arg_hash = Base.hash(body)
                    arg_var = hash_expr[arg_hash]
                    push!(block.args, arg_var)
                    # push!(block.args, :($(gensym()) = $arg_var))
                end

                scopes = setdiff(scopes,[scope])
                scope = argmax(scope->scope_level[scope],scopes)
                hash_scope[hash] = scope
                rank = 1 + scope_rank[scope]
                scope_rank[scope] = rank
    
                var = Symbol("var_$(scope)_$(rank)")
                hash_expr[hash] = var
    
                expr = Expr(:(->), node.args[1], block) #  LambdaTerm(block, node.args, x -> prototype(node))
    
                assignment = :($var = $expr) #  StatementTerm(var_term, expr, prototype(node))
                block = scope_block[scope]
                push!(block.args, assignment)
                return scopes
            elseif node.head === :for
                body = node.args[2]
                for_variables = node.args[1].args[1]
                args = (for_variables isa Expr) ? for_variables : (for_variables, )
                block = Expr(:block)
                # TODO: optimize loop ranges?

                scope = (length(args) > 0) ? args[1] : gensym()
                scope_level[scope] = depth + 1
                scope_rank[scope] = 0                
                scope_block[scope] = block
    
                for arg in args 
                    arg_hash = Base.hash(arg)
                    hash_expr[arg_hash] = arg
                    hash_scope[arg_hash] = scope
                end
                
                scopes = visit(body, depth + 1, scope; has_value = false)
                # if length(block.args) == 0 # at least 1 statement
                #     arg_hash = Base.hash(body)
                #     arg_var = hash_expr[arg_hash]
                #     push!(block.args, arg_var)
                #     # push!(block.args, :($(gensym()) = $arg_var))
                # end

                scope = last_scope
                scopes = [scope]
                hash_scope[hash] = scope
                hash_expr[hash] = :nothing
                
                expr = Expr(:for, node.args[1], block)
                block = scope_block[scope]
                push!(block.args, expr)
                return scopes
            elseif node.head in [:(=), :(+=), :(-=), :(*=), :(/=)] # TODO: do we need to make it complete? other operators: \=  ÷=  %=  ^=  &=  |=  ⊻=  >>>=  >>=  <<=
                args = node.args[2]
                var = node.args[1]
                visit(args, depth, last_scope)
                scope = last_scope
                scopes = [scope]
                hash_scope[hash] = scope
                
                rank = 1 + scope_rank[scope]
                scope_rank[scope] = rank

                hash_expr[hash] = :nothing
                hash_expr[Base.hash(var)] = var
                expr = hash_expr[Base.hash(args)]

                assignment = Expr(node.head, var, expr)
                # assignment = :($var $(node.head) $expr)
                block = scope_block[scope]
                push!(block.args, assignment)
                return scopes
            else
                @show node
                dump(node)
                @show typeof(node)
                error("A")
            end
        else
            @show node
            @show typeof(node)
            error("B")
        end


    end

    visit(node; has_value = false)
    scope_block[root]
end


function ast_optimize(expr, loop_var_range)

    # TODO: ugly. find a better way to simplify it
    unrolled_vars = Set([:field_1, :field_2])
    expr2 = ast_loop_unroll(expr, unrolled_vars) |> ast_constant_folding

    unrolled_vars_2 = Set([:face_around_1 , :face_around_2])
    expr3 = ast_loop_unroll(expr2, unrolled_vars_2) |> ast_constant_folding


    expr4, var_count = ast_flatten(expr3, 0)
    expr4

    expr5 = ast_array_unroll(expr4)

          
    expr6, _ = ast_tabulate(expr5, 0, loop_var_range) 
    
    # expr6

    expr7 = ast_remove_dead_code(expr6)

    expr8 = ast_topological_sort(expr7)
    expr8, var_count = ast_flatten(expr8, var_count)
    expr8 = ast_topological_sort(expr8)
    expr8
end

function generate_matrix_assembly_template(term, field_n_faces_around_trial, field_n_faces_around_test, dependencies, bindings)
    # expr_L: a list of block in the order of (global, face, point, field_trial, field_test, face_around_trial, face_around_test, dof_trial, dof_test)
    # dependencies: (space_trial, space_test, quadrature, alloc)
    # bindings: (face, point, field_trial, field_test, face_around_trial, face_around_test, dof_trial, dof_test)
    (space_trial, space_test, quadrature, alloc) = dependencies
    (face, point, field_trial_symbol, field_test_symbol, face_around_trial_symbol, face_around_test_symbol, dof_trial, dof_test) = bindings # TODO: replace field_symbol and face_around_symbol
    # v = last_symbol(expr_L)

    nfields_trial = length(field_n_faces_around_trial)
    nfields_test = length(field_n_faces_around_test)

    max_face_around_trial = max(field_n_faces_around_trial...)
    max_face_around_test = max(field_n_faces_around_test...)

    field_n_faces_around_trial = :(($(field_n_faces_around_trial...), ))
    field_n_faces_around_test = :(($(field_n_faces_around_test...), ))

    assignment = quote
        domain = $(GT.domain)($quadrature)
        T = eltype($alloc)
        z = zero(T)
        nfaces = $(GT.num_faces)(domain)
        face_npoints = $(GT.num_points_accessor)($quadrature)

        be = alloc_zeros("be",Any, $max_face_around_test, $max_face_around_trial, $nfields_test, $nfields_trial)
        for $field_trial_symbol in 1:$nfields_trial
            n_faces_around_trial_init = $field_n_faces_around_trial[$field_trial_symbol]
            for $field_test_symbol in 1:$nfields_test
                n_faces_around_test_init = $field_n_faces_around_test[$field_test_symbol]
                for $face_around_trial_symbol in 1:n_faces_around_trial_init
                    for $face_around_test_symbol in 1:n_faces_around_test_init
                        var_str_init = "be_for_$($field_trial_symbol)_$($field_test_symbol)_$($face_around_trial_symbol)_$($face_around_test_symbol)"
                        be[$face_around_test_symbol, $face_around_trial_symbol, $field_test_symbol, $field_trial_symbol] = 
                            alloc_zeros(var_str_init, T, $(GT.max_num_reference_dofs)($(GT.field)($space_test,$field_test_symbol)),$(GT.max_num_reference_dofs)($(GT.field)($space_trial,$field_trial_symbol)))
                    end
                end
            end
        end

        for $face = 1:nfaces
            npoints = face_npoints($face)
            for $field_trial_symbol in 1:$nfields_trial
                n_faces_around_trial_zero = $field_n_faces_around_trial[$field_trial_symbol]
                for $field_test_symbol in 1:$nfields_test
                    n_faces_around_test_zero = $field_n_faces_around_test[$field_test_symbol]
                    for $face_around_trial_symbol in 1:n_faces_around_trial_zero
                        for $face_around_test_symbol in 1:n_faces_around_test_zero
                            fill!(be[$face_around_test_symbol, $face_around_trial_symbol, $field_test_symbol, $field_trial_symbol], z)
                        end
                    end
                end
            end

            for $point in 1:npoints
                for $field_trial_symbol in 1:$nfields_trial
                    n_faces_around_trial = $field_n_faces_around_trial[$field_trial_symbol]
                    for $field_test_symbol in 1:$nfields_test
                        n_faces_around_test = $field_n_faces_around_test[$field_test_symbol]
                        for $face_around_trial_symbol in 1:n_faces_around_trial
                            dofs_trial = $(GT.dofs_accessor)($(GT.field)($space_trial,$field_trial_symbol),domain)($face,$face_around_trial_symbol)
                            ndofs_trial = length(dofs_trial)
                            for $face_around_test_symbol in 1:n_faces_around_test
                                dofs_test = $(GT.dofs_accessor)($(GT.field)($space_test,$field_test_symbol),domain)($face,$face_around_test_symbol)
                                ndofs_test = length(dofs_test)
                                for $dof_trial in 1:ndofs_trial
                                    for $dof_test in 1:ndofs_test
                                        v = $term
                                        be[$face_around_test_symbol, $face_around_trial_symbol, $field_test_symbol, $field_trial_symbol][$dof_test, $dof_trial] += v
                                    end
                                end
                            end
                        end
                    end
                end
            end

            for $field_trial_symbol in 1:$nfields_trial
                n_faces_around_trial_contribute = $field_n_faces_around_trial[$field_trial_symbol]
                for $field_test_symbol in 1:$nfields_test
                    n_faces_around_test_contribute = $field_n_faces_around_test[$field_test_symbol]
                    for $face_around_trial_symbol in 1:n_faces_around_trial_contribute
                        for $face_around_test_symbol in 1:n_faces_around_test_contribute
                            $(GT.contribute!)($alloc,be[$face_around_test_symbol, $face_around_trial_symbol, $field_test_symbol, $field_trial_symbol],
                                                $(GT.dofs_accessor)($(GT.field)($space_test,$field_test_symbol),domain)($face,$face_around_test_symbol),
                                                $(GT.dofs_accessor)($(GT.field)($space_trial,$field_trial_symbol),domain)($face,$face_around_trial_symbol),
                                                $field_test_symbol,$field_trial_symbol)
                        end
                    end
                end
            end

        end
    end

    block = assignment
    ast_clean_up!(block)
    block
end



function generate_vector_assembly_template(term, field_n_faces_around, dependencies, bindings)
    # expr_L: a list of block in the order of (global, face, point, field_trial, field_test, face_around_trial, face_around_test, dof_trial, dof_test)
    # dependencies: (space_trial, space_test, quadrature, alloc)
    # bindings: (face, point, field_trial, field_test, face_around_trial, face_around_test, dof_trial, dof_test)
    (space, quadrature, alloc) = dependencies
    (face, point, field_symbol, face_around_symbol, dof) = bindings
    # v = last_symbol(expr_L)

    nfields = length(field_n_faces_around)

    max_face_around = max(field_n_faces_around...)

    field_n_faces_around = :(($(field_n_faces_around...), ))

    assignment = quote
        domain = $(GT.domain)($quadrature)
        T = eltype($alloc)
        z = zero(T)
        nfaces = $(GT.num_faces)(domain)
        face_npoints = $(GT.num_points_accessor)($quadrature)

        be = alloc_zeros("be",Any, $max_face_around, $nfields)
        for $field_symbol in 1:$nfields
            n_faces_around_init = $field_n_faces_around[$field_symbol]
            for $face_around_symbol in 1:n_faces_around_init
                var_str_init = "be_for_$($field_symbol)_$($face_around_symbol)"
                be[$face_around_symbol, $field_symbol] = 
                    alloc_zeros(var_str_init, T, $(GT.max_num_reference_dofs)($(GT.field)($space,$field_symbol)) )
            end
        end

        for $face = 1:nfaces
            npoints = face_npoints($face)
            for $field_symbol in 1:$nfields
                n_faces_around_zero = $field_n_faces_around[$field_symbol]
                for $face_around_symbol in 1:n_faces_around_zero
                    fill!(be[$face_around_symbol, $field_symbol], z)
                end
            end

            for $point in 1:npoints
                for $field_symbol in 1:$nfields
                    n_faces_around = $field_n_faces_around[$field_symbol]
                    for $face_around_symbol in 1:n_faces_around
                        dofs = $(GT.dofs_accessor)($(GT.field)($space,$field_symbol),domain)($face,$face_around_symbol)
                        ndofs = length(dofs)
                        for $dof in 1:ndofs
                            v = $term
                            be[$face_around_symbol, $field_symbol][$dof] += v
                        end
                    end
                end
            end

            for $field_symbol in 1:$nfields
                n_faces_around_contribute = $field_n_faces_around[$field_symbol]
                for $face_around_symbol in 1:n_faces_around_contribute
                    $(GT.contribute!)($alloc,be[$face_around_symbol, $field_symbol],
                                        $(GT.dofs_accessor)($(GT.field)($space,$field_symbol),domain)($face,$face_around_symbol),
                                        $field_symbol)
                end

            end

        end
    end

    block = assignment
    ast_clean_up!(block)
    block
end

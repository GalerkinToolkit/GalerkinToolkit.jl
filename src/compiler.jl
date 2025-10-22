
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
    coordinate(at_any_index(quadrature_accessor(quadrature)))
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
    :( coordinate(at_point(at_face(quadrature_accessor($quadrature),$face),$point)) )
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
    weight(at_any_index(quadrature_accessor(quadrature)))
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
    :(weight(at_point(at_face(quadrature_accessor($quadrature),$face),$point)) )
end

function form_argument_quantity(space::AbstractSpace,arg,the_field=1)
    quantity() do opts
        space_domain = GT.domain(space)
        domain = GT.domain(opts)
        D = num_dims(space_domain)
        d = num_dims(domain)
        faces_around = GT.faces_around(domain)
        index = GT.index(opts)
        face = domain_face_index(index)
        field = field_index(index,arg)
        dof = dof_index(index,arg)

        f = leaf_term(GT.value; is_compile_constant=true)
        the_field_term = leaf_term(the_field; is_compile_constant=true)
        (space_term, domain_term, face_term, field_term, dof_term, face_around_term) = map(leaf_term, (space, domain, face, field, dof, face_around))
        to_leaf_term = x -> (x isa AbstractTerm) ? x : leaf_term(x)
        is_reference = is_reference_domain(space_domain)

        if D == d
            # the_face_around_term = leaf_term(nothing, is_compile_constant=true)
            the_face_around = 1
            the_face_around_term = leaf_term(the_face_around;is_compile_constant=true)
            # face_around_term = leaf_term(nothing)
            face_around_term = leaf_term(face_around_index(index,arg))
            dependencies = (f,space_term,domain_term,face_term,the_field_term,field_term,dof_term,the_face_around_term,face_around_term)
            FormArgumentTerm(dependencies, D, is_reference, :interior)
        elseif D==d+1 && faces_around !== nothing
            #NB this is not constant anymore
            #the_face_around_term = leaf_term(face_around, is_compile_constant=true)
            #face_around_term = leaf_term(face_around)
            the_face_around_term = RefTerm(leaf_term(GT.faces_around(domain)),leaf_term(face))
            face_around_term = the_face_around_term
            # face_around_term = leaf_term(face_around_index(index,arg))
            dependencies = (f,space_term,domain_term,face_term,the_field_term,field_term,dof_term,the_face_around_term,face_around_term)
            FormArgumentTerm(dependencies, D, is_reference, :boundary)
        else
            face_around_term = leaf_term(face_around_index(index,arg))
            the_face_around = :the_face_around
            the_face_around_term = leaf_term(the_face_around)
            n_faces_around = leaf_term(2;is_compile_constant=true) # Hard coded! But OK in practice.
            dependencies = (f,space_term,domain_term,face_term,the_field_term,field_term,dof_term,the_face_around_term,face_around_term)
            SkeletonTerm(FormArgumentTerm(dependencies, D, is_reference, :skeleton), n_faces_around,the_face_around)
        end
    end
end

struct FormArgumentTerm <: AbstractTerm
    dependencies
    D::Int
    is_reference::Bool
    integral_type::Symbol
    #d::Int
end

function prototype(term::FormArgumentTerm)
    @error "We dont want untabulated shape functions"
end

function dependencies(term::FormArgumentTerm)
    term.dependencies
end

function replace_dependencies(term::FormArgumentTerm,dependencies)
    FormArgumentTerm(dependencies, term.D, term.is_reference, term.integral_type)
end

function replace_the_face_around(term::FormArgumentTerm,the_face_around)
    (f,space,domain,face,the_field,field,dof,_,face_around) = term.dependencies
    dependencies = (f,space,domain,face,the_field,field,dof,the_face_around,face_around)
    FormArgumentTerm(dependencies, term.D, term.is_reference, term.integral_type)
end

function expression(term::FormArgumentTerm)
    @error "We dont want untabulated shape functions"
    (f,space,domain,face,the_field,field,dof,the_face_around,face_around) = map(expression,term.dependencies)
    # TODO: not tabulated. we need another solution to tabulate it or this cannot be tabulated anyway
    #:(form_argument_accessor($f,$space,$domain,$the_field)($face,$the_face_around)($dof,$field,$face_around))
    :(ifelse($face_around == $the_face_around && $field == $the_field,shape_functions($f,at_point(at_face_around(at_face(space_accessor($space,$domain),$face),$the_face_around),$point))[$dof],zero(eltype(shape_functions($f,at_any_index(space_accessor($space,$domain)))))))
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
    faces_around = GT.faces_around(domain)
    index = GT.index(opts)
    face = domain_face_index(index)
    f = value
    if D == d
        face_around = nothing
        dependencies = map(leaf_term,(f,uh,domain,face,face_around))
        DiscreteFieldTerm(dependencies)
    elseif D==d+1 && faces_around !== nothing
        dependencies = map(leaf_term,(f,uh,domain,face,))
        the_face_around_term = face_around_term(index,d,D)
        DiscreteFieldTerm((dependencies...,the_face_around_term))
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

function replace_the_face_around(term::DiscreteFieldTerm,the_face_around)
    (f,uh,domain,face,_) = term.dependencies
    dependencies = (f,uh,domain,face,the_face_around)
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
    uh_acc = field_accessor(uh,quadrature;tabulate=(f,))
    uh_point = at_any_index(uh_acc)
    field(f,uh_point)
    #prototype(discrete_field_accessor(f,uh,quadrature))
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
    mesh_acc = :(quadrature_accessor($quadrature,Val($D)))
    mesh_ready = :(at_point(at_face_around(at_face($mesh_acc,$face),$the_face_around),$point))
    z = :(zero(eltype(shape_functions($f,at_any_index(tabulate($f,space_accessor($space,$mesh_acc)))))))
    mask = :($face_around == $the_face_around && $field == $the_field)
    space_acc = :(tabulate($f,space_accessor($space,$mesh_acc)))
    sfun = :(shape_functions($f,at_point(at_face_around(at_face($space_acc,$face),$the_face_around),$mesh_ready))[$dof])
    :(ifelse($mask,$sfun,$z))
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
    :(field($f,at_point(at_face_around(at_face(tabulate($f,field_accessor($uh,$quadrature)),$face),$face_around),$point)))
    #:(discrete_field_accessor($f,$uh,$quadrature)($face,$face_around)($point))
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
    acc = quadrature_accessor(quadrature,vD)
    coordinate(f,at_any_index(acc))
    #prototype(physical_map_accessor(f,quadrature,vD))
end

#function prototype(term::TabulatedTerm{<:PhysicalMapTerm})
#    parent = term.parent
#    quadrature = prototype(term.quadrature)
#    (f,mesh,vD,face) = map(prototype,dependencies(parent))
#    prototype(physical_map_accessor(f,quadrature,vD))
#end

function expression_TabulatedTerm(parent::PhysicalMapTerm,quadrature,point)
    (f,mesh,vD,face) = map(expression,parent.dependencies)
    :(coordinate($f,at_point(at_face(quadrature_accessor($quadrature,$vD),$face),$point)))
    #:(physical_map_accessor($f,$quadrature,$vD)($face)($point))
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
        if GT.faces_around(domain(opts)) === nothing
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
    mesh = GT.mesh(domain(quadrature))
    vD = Val(num_dims(mesh))
    acc = compute(unit_normal,quadrature_accessor(quadrature,vD))
    unit_normal(at_any_index(acc))
    #prototype(unit_normal_accessor(quadrature))
end

#function prototype(term::TabulatedTerm{<:UnitNormalTerm})
#    parent = term.parent
#    quadrature = prototype(term.quadrature)
#    #(face,the_face_around) = map(prototype,dependencies(parent))
#    prototype(unit_normal_accessor(quadrature))
#end

function expression_TabulatedTerm(parent::UnitNormalTerm,quadrature,point)
    (face,the_face_around) = map(expression,dependencies(parent))
    mesh = :(mesh(domain($quadrature)))
    :(unit_normal(at_point(at_face_around(at_face(compute(unit_normal,quadrature_accessor($quadrature,Val(num_dims($mesh)))),$face),$the_face_around),$point)))
    #:(unit_normal_accessor($quadrature)($face,$the_face_around)($point))
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

#TODO: code generation
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


function generate_assemble_vector(contribution::DomainContribution,space::AbstractSpace;parameters=(),optimize_options=nothing)
    term_0 = write_assemble_vector(contribution,space; optimize_options)
    term_1 = optimize(term_0)
    term_2 = parametrize(term_1,parameters...)
    term_3, captured_data = capture(term_2)
    expr_0 = expression(term_3)
    f = evaluate(expr_0,captured_data)
    # expr_1 = statements_expr(expr_0)
    # f = evaluate(expr_1,captured_data)
    f
end

function write_assemble_vector(contribution::DomainContribution,space::AbstractSpace;optimize_options=nothing)
    quadrature = GT.quadrature(contribution)
    domain = GT.domain(quadrature)
    arity = Val(1)
    index = GT.index(arity)
    term = GT.term(contribution,index)
    space_term = leaf_term(space)
    quadrature_term = leaf_term(quadrature)
    alloc_arg = :alloc
    alloc = leaf_term(alloc_arg)

    field_n_faces_around = [max_num_faces_around(GT.domain(field_space),domain) for field_space in fields(space)]

    VectorAssemblyTerm(term,space_term,quadrature_term,alloc,alloc_arg,index,field_n_faces_around, optimize_options)
end


function generate_assemble_matrix(contribution::DomainContribution,space_trial::AbstractSpace,space_test::AbstractSpace;parameters=(),optimize_options=nothing)
    term_0 = write_assemble_matrix(contribution,space_trial, space_test;optimize_options)
    term_1 = optimize(term_0)
    term_2 = parametrize(term_1,parameters...)
    term_3, captured_data = capture(term_2)
    expr_0 = expression(term_3)
    f = evaluate(expr_0,captured_data)
    # expr_1 = statements_expr(expr_0)
    # f = evaluate(expr_1,captured_data)
    f
end

function write_assemble_matrix(contribution::DomainContribution,space_trial::AbstractSpace,space_test::AbstractSpace;optimize_options=nothing)
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

    field_n_faces_around_trial = [max_num_faces_around(GT.domain(field_space),domain) for field_space in fields(space_trial)]
    field_n_faces_around_test = [max_num_faces_around(GT.domain(field_space),domain) for field_space in fields(space_test)]
    

    MatrixAssemblyTerm(term,space_trial_term,space_test_term,quadrature_term,alloc,alloc_arg,index,field_n_faces_around_trial,field_n_faces_around_test,optimize_options)
end

#TODO: code generation
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
    q_acc = quadrature_accessor(quadrature)
    init = zero(eltype(contribs))
    for face in 1:nfaces
        point_v = face_point_v(face)
        npoints = num_points(at_face(q_acc,face))
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
    field_n_faces_around::Vector{Int}
    optimize_options
end

struct MatrixAssemblyTerm <: AbstractTerm
    term::AbstractTerm
    space_trial::AbstractTerm
    space_test::AbstractTerm
    quadrature::AbstractTerm
    alloc::AbstractTerm
    alloc_arg::Symbol # gets reduced
    index::Index{2} # gets reduced
    field_n_faces_around_trial::Vector{Int}
    field_n_faces_around_test::Vector{Int}
    optimize_options
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
    (;alloc_arg,index,field_n_faces_around,optimize_options) = term
    VectorAssemblyTerm(term2,space,quadrature,alloc,alloc_arg,index,field_n_faces_around,optimize_options)
end


function replace_dependencies(term::MatrixAssemblyTerm,dependencies)
    (term2,space_trial,space_test,quadrature,alloc) = dependencies
    (;alloc_arg,index,field_n_faces_around_trial,field_n_faces_around_test,optimize_options) = term
    MatrixAssemblyTerm(term2,space_trial,space_test,quadrature,alloc,alloc_arg,index,field_n_faces_around_trial,field_n_faces_around_test,optimize_options)
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
    options = term.optimize_options
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

    body_optimized = ast_optimize_with_options(body, loop_var_range, options)
    expr = :($alloc_arg->$body_optimized)

    expr
end


function expression(term::MatrixAssemblyTerm)
    (term2,space_trial,space_test,quadrature,alloc) = map(expression,dependencies(term))
    (alloc_arg,face,point,field_trial,field_test,face_around_trial,face_around_test,dof_trial,dof_test) = bindings(term)
    options = term.optimize_options
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

    body_optimized = ast_optimize_with_options(body, loop_var_range, options)
    # display(body_optimized)
    expr = quote 
        $alloc_arg -> begin
            if GT.num_faces(GT.domain($quadrature)) == 0
                return
            end
            $body_optimized
        end
    end
    expr
end

function scalar_assembly_loop(face_point_v,init,quadrature)
    domain = GT.domain(quadrature)
    nfaces = num_faces(domain)
    z = init
    q_faces = each_face(quadrature)
    for face in 1:nfaces
        point_v = face_point_v(face)
        npoints = num_points(q_faces[face])
        for point in 1:npoints
            z += point_v(point)
        end
    end
    z
end


function alloc_zeros(symbol, args...)
    if args[1] === Any
        Array{args[1]}(undef, args[2:end])
    else
        zeros(args...)
    end
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





# TODO: this is used in scalar and sample assembly. remove it and use code generation
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




function ast_optimize(expr, loop_var_range, options = nothing)
    # unroll -> flatten -> array_unroll -> tabulate -> topological_sort -> array_aliasing 
    # do not remove: unroll -> array_unroll,  remove_dead_code, flatten
    expr2 = ast_loop_unroll(expr) |> ast_constant_folding

    expr3 = ast_loop_unroll(expr2) |> ast_constant_folding
    var_count = 0

    expr4, var_count = ast_flatten(expr3, var_count)
    expr4

    expr5 = ast_array_unroll(expr4)

          
    expr6, var_count = ast_tabulate(expr5, var_count, loop_var_range) 
    
    # expr6

    expr7 = ast_remove_dead_code(expr6)

    if options === nothing || options.topological_sort == true
        expr8 = ast_topological_sort(expr7)
        expr8, var_count = ast_flatten(expr8, var_count)
        expr8 = ast_topological_sort(expr8)
    else
        expr8 = expr7
    end
    
    if options === nothing || options.array_aliasing == true
        expr9 = ast_array_aliasing(expr8) 
    else
        expr9 = expr8
    end

    expr9 = ast_remove_dead_code(expr9)

    expr9
end



function ast_optimize_2(expr, loop_var_range, options = nothing)
    var_count = 0
    
    expr2, var_count = ast_flatten(expr, var_count)

    expr3, var_count = ast_tabulate(expr2, var_count, loop_var_range) 
    
    expr4 = ast_loop_unroll(expr3)  |> ast_array_unroll |> ast_constant_folding

    expr5 = ast_loop_unroll(expr4)  |> ast_array_unroll |> ast_constant_folding
    
    # remove dead code twice. this is important because we have ifelse statements flattened
    expr6 = ast_remove_dead_code(expr5) |> ast_constant_folding |> ast_remove_dead_code 


    expr7 = if options === nothing || options.topological_sort == true
        ast_topological_sort(expr6)
    else
        expr6
    end

    expr8 = if options === nothing || options.array_aliasing == true
        ast_array_aliasing(expr7)
    else
        expr7
    end

    expr8 = ast_remove_dead_code(expr8)

    expr9, var_count = ast_flatten(expr8, var_count) 
    # expr9, var_count = ast_tabulate(expr9, var_count, loop_var_range) 

    expr10 = expr9 |> ast_constant_folding 
    expr10 = if options === nothing || options.topological_sort == true
        ast_topological_sort(expr10)
    else
        expr10
    end

    expr10 = ast_remove_dead_code(expr10)

    expr10
end


function ast_optimize_3(expr, loop_var_range, options = nothing) # optimize steps with hoist only, but not tabulate
    # TODO: implement
    expr
end

function ast_optimize_with_options(expr, loop_var_range, options)
    if options === nothing || options.order == 1
        ast_optimize(expr, loop_var_range, options)
    elseif options.order == 2 
        ast_optimize_2(expr, loop_var_range, options)
    else # no optimize
        expr
    end
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
        domain = GT.domain($quadrature)
        T = eltype($alloc)
        z = zero(T)
        nfaces = GT.num_faces(domain)
        q_acc = GT.quadrature_accessor($quadrature)

        be = alloc_zeros("be",Any, $max_face_around_test, $max_face_around_trial, $nfields_test, $nfields_trial)
        for $field_trial_symbol in 1:$nfields_trial 
            n_faces_around_trial_init = $field_n_faces_around_trial[$field_trial_symbol]
            for $field_test_symbol in 1:$nfields_test
                n_faces_around_test_init = $field_n_faces_around_test[$field_test_symbol]
                for $face_around_trial_symbol in 1:n_faces_around_trial_init
                    for $face_around_test_symbol in 1:n_faces_around_test_init
                        var_str_init = "be_for_$($field_test_symbol)_$($field_trial_symbol)_$($face_around_test_symbol)_$($face_around_trial_symbol)"
                        be[$face_around_test_symbol, $face_around_trial_symbol, $field_test_symbol, $field_trial_symbol] = 
                            alloc_zeros(var_str_init, T, GT.max_num_reference_dofs(GT.field($space_test,$field_test_symbol)),GT.max_num_reference_dofs(GT.field($space_trial,$field_trial_symbol)))
                    end
                end
            end
        end

        for $face = 1:nfaces
            npoints = num_points(at_face(q_acc,$face))
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
                            dofs_trial = GT.dofs_accessor(GT.field($space_trial,$field_trial_symbol),domain)($face,$face_around_trial_symbol)
                            ndofs_trial = length(dofs_trial)
                            for $face_around_test_symbol in 1:n_faces_around_test
                                dofs_test = GT.dofs_accessor(GT.field($space_test,$field_test_symbol),domain)($face,$face_around_test_symbol)
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
                            GT.contribute!($alloc,be[$face_around_test_symbol, $face_around_trial_symbol, $field_test_symbol, $field_trial_symbol],
                                                GT.dofs_accessor(GT.field($space_test,$field_test_symbol),domain)($face,$face_around_test_symbol),
                                                GT.dofs_accessor(GT.field($space_trial,$field_trial_symbol),domain)($face,$face_around_trial_symbol),
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
        domain = GT.domain($quadrature)
        T = eltype($alloc)
        z = zero(T)
        nfaces = GT.num_faces(domain)
        q_acc = GT.quadrature_accessor($quadrature)

        be = alloc_zeros("be",Any, $max_face_around, $nfields)
        for $field_symbol in 1:$nfields
            n_faces_around_init = $field_n_faces_around[$field_symbol]
            for $face_around_symbol in 1:n_faces_around_init
                var_str_init = "be_for_$($field_symbol)_$($face_around_symbol)"
                be[$face_around_symbol, $field_symbol] = 
                    alloc_zeros(var_str_init, T, GT.max_num_reference_dofs(GT.field($space,$field_symbol)) )
            end
        end

        for $face = 1:nfaces
            npoints = num_points(at_face(q_acc,$face))
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
                        dofs = GT.dofs_accessor(GT.field($space,$field_symbol),domain)($face,$face_around_symbol)
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
                    GT.contribute!($alloc,be[$face_around_symbol, $field_symbol],
                                        GT.dofs_accessor(GT.field($space,$field_symbol),domain)($face,$face_around_symbol),
                                        $field_symbol)
                end

            end

        end
    end

    block = assignment
    ast_clean_up!(block)
    block
end

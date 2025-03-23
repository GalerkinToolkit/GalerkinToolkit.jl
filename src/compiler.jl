
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
struct LeafTerm{A,B,C} <: AbstractTerm
    value::A
    prototype::B
    is_compile_constant::Val{C}
end

# The prototye is the same as the value by default
leaf_term(v; prototype=v, is_compile_constant::Val=Val(false)) = LeafTerm(v, prototype, is_compile_constant)

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

struct LambdaTerm{A,B} <: AbstractTerm
    body::A
    args::B # all symbols in here get reduced,
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

struct CallTerm{A,B} <: AbstractTerm
    callee::A
    args::B
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

struct RefTerm{A,B} <: AbstractTerm
    container::A
    index::B
end

function prototype(a::RefTerm)
    zero(eltype(prototype(a.container)))
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
    :($(container)[$(args...)])
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
struct Quantity{A} <: AbstractQuantity
    term::A
end
function quantity(term)
    Quantity(term)
end

function term(qty::Quantity,opts)
    qty.term(opts)
end

function call(callee,args::AbstractQuantity...)
    f = uniform_quantity(callee;is_compile_constant=Val(true))
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

function uniform_quantity(v;is_compile_constant=Val(false))
    quantity() do opts
        leaf_term(v; is_compile_constant)
    end
end

function coordinate_quantity(quadrature)
    quantity() do opts
        face = domain_face_index(opts.index)
        point = point_index(opts.index)
        dependencies = map(leaf_term,(quadrature,face,point))
        CoordinateTerm(dependencies)
    end
end

struct CoordinateTerm{A} <: AbstractTerm
    dependencies::A
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
        WeightTerm(dependencies)
    end
end

struct WeightTerm{A} <: AbstractTerm
    dependencies::A
end

function prototype(term::WeightTerm)
    (quadrature,face,point) = map(prototype,dependencies(term))
    prototype(weight_accessor(quadrature))
end

function dependencies(term::WeightTerm)
    term.dependencies
end

function replace_dependencies(term::WeightTerm,dependencies)
    WeightTerm(dependencies)
end

function expression(term::WeightTerm)
    (quadrature,face,point) = map(expression,term.dependencies)
    :( weight_accessor($quadrature)($face)($point) )
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

        f = leaf_term(GT.value; is_compile_constant=Val(true))
        the_field_term = leaf_term(the_field; is_compile_constant=Val(true))
        (space_term, domain_term, face_term, field_term, dof_term, face_around_term) = map(leaf_term, (space, domain, face, field, dof, face_around))
        to_leaf_term = x -> (x isa AbstractTerm) ? x : leaf_term(x)

        if D == d
            the_face_around_term = leaf_term(nothing, is_compile_constant=Val(true))
            face_around_term = leaf_term(nothing)
            dependencies = (f,space_term,domain_term,face_term,the_field_term,field_term,dof_term,the_face_around_term,face_around_term)
            FormArgumentTerm(dependencies)
        elseif D==d+1 && face_around !== nothing
            the_face_around_term = leaf_term(face_around, is_compile_constant=Val(true))
            face_around_term = leaf_term(face_around)
            dependencies = (f,space_term,domain_term,face_term,the_field_term,field_term,dof_term,the_face_around_term,face_around_term)
            FormArgumentTerm(dependencies)
        else
            face_around_term = leaf_term(face_around_index(index,arg))
            the_face_around_term = leaf_term(:the_face_around)
            # the_face_around = leaf_term(:the_face_around, is_compile_constant=Val(true))
            n_faces_around = leaf_term(2) # Hard coded! But OK in practice.
            dependencies = (f,space_term,domain_term,face_term,the_field_term,field_term,dof_term,the_face_around_term,face_around_term)
            SkeletonTerm(FormArgumentTerm(dependencies),n_faces_around,the_face_around)
        end
    end
end

struct FormArgumentTerm{A} <: AbstractTerm
    dependencies::A
end

function prototype(term::FormArgumentTerm)
    (f,space,domain,face,the_field,field,dof,the_face_around,face_around) = map(prototype,dependencies(term))
    prototype(form_argument_accessor(space,domain,1))
end

function dependencies(term::FormArgumentTerm)
    term.dependencies
end

function replace_dependencies(term::FormArgumentTerm,dependencies)
    FormArgumentTerm(dependencies)
end

function replace_the_face_around(term::FormArgumentTerm,the_face_around)
    (f,space,domain,face,the_field,field,dof,_,face_around) = term.dependencies
    dependencies = (f,space,domain,face,the_field,field,dof,the_face_around,face_around)
    FormArgumentTerm(dependencies)
end

function expression(term::FormArgumentTerm)
    (f,space,domain,face,the_field,field,dof,the_face_around,face_around) = map(expression,term.dependencies)
    :(form_argument_accessor($f,$space,$domain,$the_field)($face,$the_face_around)($dof,$field,$face_around))
end

function optimize(term::CallTerm{<:FormArgumentTerm,<:Tuple{<:CoordinateTerm}})
    coords = optimize(term.args[1])
    (quadrature,face,point) = map(optimize,dependencies(coords))
    TabulatedTerm(term.callee,quadrature,point)
end

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
        DiscreteTerm(dependencies)
    elseif D==d+1 && face_around !== nothing
        dependencies = map(leaf_term,(f,uh,domain,face,face_around))
        DiscreteTerm(dependencies)
    else
        face_around = :the_face_around
        n_faces_around = leaf_term(2) # Hard coded! But OK in practice.
        dependencies = map(leaf_term,(f,uh,domain,face,face_around))
        SkeletonTerm(DiscreteTerm(dependencies),n_faces_around,face_around)
    end
end

struct DiscreteFieldTerm{A} <: AbstractTerm
    dependencies::A
end

function optimize(term::CallTerm{<:DiscreteFieldTerm,Tuple{<:CoordinateTerm}})
    coords = term.args[1]
    (quadrature,face,point) = dependencies(coords)
    TabulatedTerm(term.calee,quadrature,point)
end

struct TabulatedTerm{A,B,C} <: AbstractTerm
    parent::A
    quadrature::B
    point::C
end

function dependencies(term::TabulatedTerm)
    (term.parent,term.quadrature,term.point)
end

function replace_dependencies(term::TabulatedTerm,dependencies)
    parent, quadrature, point = dependencies
    TabulatedTerm(parent,quadrature, point)
end

function prototype(term::TabulatedTerm{<:FormArgumentTerm})
    parent = term.parent
    quadrature = prototype(term.quadrature)
    (space,domain,face,the_field,field,dof,the_face_around,face_around) = map(prototype,dependencies(parent))
    prototype(form_argument_accessor(space,domain,1))
end

function prototype(term::TabulatedTerm{<:DiscreteFieldTerm})
    parent = term.parent
    quadrature = prototype(term.quadrature)
    (f,uh,domain,face,face_around) = map(prototype,dependencies(parent))
    prototype(discrete_field_accessor(f,uh,quadrature))
end

function expression(term::TabulatedTerm{<:FormArgumentTerm})
    point = expression(term.point)
    form_arg = term.parent
    quadrature = expression(term.quadrature)
    (f,space,domain,face,the_field,field,dof,the_face_around,face_around) = map(expression,form_arg.dependencies)
    :(form_argument_accessor($f,$space,$quadrature,$the_field)($face,$the_face_around)($point)($dof,$field,$face_around))
end

function expression(term::TabulatedTerm{<:DiscreteFieldTerm})
    point = expression(term.point)
    form_arg = term.parent
    quadrature = expression(term.quadrature)
    (f,uh,domain,face,face_around) = map(expression,form_arg.dependencies)
    :(discrete_field_accessor($f,$space,$quadrature)($face,$face_around)($point))
end

struct SkeletonTerm{A,B,C} <: AbstractTerm
    parent::A
    n_faces_around::B
    the_face_around::C # Gets reduced
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

function optimize(term::RefTerm{<:SkeletonTerm})
    the_face_around = term.index
    replace_the_face_around(term.container.parent,the_face_around)
end

function generate_assemble_scalar(contribution::DomainContribution;parameters=())
    # TODO how to sort the calls to optimize, parametrize, capture, expression, statements, etc
    # needs to be though carefully. We provably will need several calls to optimize and more IR levels
    term_0 = write_assemble_scalar(contribution,space)
    term_1 = optimize(term_0)
    term_2 = parametrize(term_1,parameters...)
    term_3, captured_data = capture(term_2)
    term_4 = statements(term_3)
    expr_0 = expression(term_4)

    f = evaluate(expr_0,captured_data)
    f
end

function generate_assemble_vector(contribution::DomainContribution,space::AbstractSpace;parameters=())
    # TODO how to sort the calls to optimize, parametrize, capture, expression, statements, etc
    # needs to be though carefully. We provably will need several calls to optimize and more IR levels
    term_0 = write_assemble_vector(contribution,space)
    term_1 = optimize(term_0)
    term_2 = parametrize(term_1,parameters...)
    term_3, captured_data = capture(term_2)
    expr_0 = expression(term_3)
    #TODO
    #expr_1 = statements(expr_0)
    expr_1 = expr_0
    f = evaluate(expr_1,captured_data)
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
    VectorAssemblyTerm(term,space_term,quadrature_term,alloc,alloc_arg,index)
end


function generate_assemble_matrix(contribution::DomainContribution,space_trial::AbstractSpace,space_test::AbstractSpace;parameters=())
    # TODO how to sort the calls to optimize, parametrize, capture, expression, statements, etc
    # needs to be though carefully. We provably will need several calls to optimize and more IR levels
    term_0 = write_assemble_matrix(contribution,space_trial, space_test)
    term_1 = optimize(term_0)
    term_2 = parametrize(term_1,parameters...)
    term_3, captured_data = capture(term_2)
    expr_0 = expression(term_3)
    #TODO
    #expr_1 = statements(expr_0)
    expr_1 = expr_0
    f = evaluate(expr_1,captured_data)
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
    MatrixAssemblyTerm(term,space_trial_term,space_test_term,quadrature_term,alloc,alloc_arg,index)
end
struct VectorAssemblyTerm{A,B,C,D,E,F} <: AbstractTerm
    term::A # <: Term
    space::B # <: Term
    quadrature::C # <: Term
    alloc::D
    alloc_arg::E # gets reduced
    index::F # gets reduced
end

struct MatrixAssemblyTerm{A,B,C,D,E,F,G} <: AbstractTerm
    term::A # <: Term
    space_trial::B # <: Term
    space_test::C
    quadrature::D # <: Term
    alloc::E
    alloc_arg::F # gets reduced
    index::G # gets reduced
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

function replace_dependencies(term::VectorAssemblyTerm,dependencies)
    (term2,space,quadrature,alloc) = dependencies
    (;alloc_arg,index) = term
    VectorAssemblyTerm(term2,space,quadrature,alloc,alloc_arg,index)
end


function replace_dependencies(term::MatrixAssemblyTerm,dependencies)
    (term2,space_trial,space_test,quadrature,alloc) = dependencies
    (;alloc_arg,index) = term
    MatrixAssemblyTerm(term2,space_trial,space_test,quadrature,alloc,alloc_arg,index)
end

# We can specialize for particular cases if needed
function expression(term::VectorAssemblyTerm)
    (term2,space,quadrature,alloc) = map(expression,dependencies(term))
    (alloc_arg,face,point,field,face_around,dof) = bindings(term)
    dof_v = :($dof -> $term2)
    block_dof_v = :( ($field,$face_around) -> $dof_v)
    point_block_dof_v = :($point -> $block_dof_v)
    face_point_block_dof_v = :( $face -> $point_block_dof_v)
    nfaces = :(num_faces($domain))
    body = :(vector_assembly_loop!($face_point_block_dof_v,$alloc,$space,$quadrature))
    expr = :($alloc_arg->$body)
    expr
end

# TODO: fix it
function expression(term::MatrixAssemblyTerm)
    (term2,space_trial,space_test,quadrature,alloc) = map(expression,dependencies(term))
    (alloc_arg,face,point,field_trial,field_test,face_around_trial,face_around_test,dof_trial,dof_test) = bindings(term)
    dof_v = :(($dof_trial, $dof_test) -> $term2)
    block_dof_v = :( ($field_trial,$field_test,$face_around_trial,$face_around_test) -> $dof_v)
    point_block_dof_v = :($point -> $block_dof_v)
    face_point_block_dof_v = :( $face -> $point_block_dof_v)
    nfaces = :(num_faces($domain))
    body = :(matrix_assembly_loop!($face_point_block_dof_v,$alloc,$space_trial,$space_test,$quadrature))
    expr = :($alloc_arg->$body)
    expr
end


function vector_assembly_loop!(face_point_block_dof_v,alloc,space,quadrature)
    domain = GT.domain(quadrature)
    nfaces = num_faces(domain)
    space_domain = GT.domain(space)
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

# TODO:fix
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

    # TODO: set matrix 
    field_face_around_be = map(fields(space_trial)) do field_space_trial
        n_trial = max_num_faces_around(GT.domain(field_space_trial),domain)
        map(fields(space_test)) do field_space_test
            n_test = max_num_faces_around(GT.domain(field_space_test),domain)
            map(1:n_trial) do _
                map(1:n_test) do _
                    m_trial = max_num_reference_dofs(field_space_trial)
                    m_test = max_num_reference_dofs(field_space_test)
                    zeros(eltype(alloc),m_trial, m_test)
                end
            end
        end
    end


    nfields_trial = num_fields(space_trial)
    nfields_test = num_fields(space_test)
    z = zero(eltype(alloc))

    for face in 1:nfaces
        point_block_dof_v = face_point_block_dof_v(face)
        # TODO: Reset
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
                                for dof_test in 1:ndofs_test
                                    v = dof_v(dof_trial, dof_test)
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
                    # function contribute!(alloc::MatrixAllocation,b,dofs_test,dofs_trial,field_test=1,field_trial=1)
                end
            end
        end
    end
    alloc
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

for op in (:+,:-,:sqrt,:abs,:abs2,:real,:imag,:conj,:transpose,:adjoint,:*,:/,:\,:^,:getindex)
    @eval begin
        function get_symbol!(index,val::typeof(Base.$op),name="";prefix=gensym)
            $( Expr(:quote,op) )
        end
    end
end

for op in (:+,:-,:sqrt,:abs,:abs2,:real,:imag,:conj,:transpose,:adjoint)
  @eval begin
      (Base.$op)(a::AbstractQuantity) = call(Base.$op,a)
  end
end

function Base.getindex(a::AbstractQuantity,b::AbstractQuantity)
    call(getindex,a,b)
end

#function Base.getindex(a::AbstractQuantity,b::Number)
#    Base.getindex(a,compile_constant_quantity(b))
#end

for op in (:+,:-,:*,:/,:\,:^)
  @eval begin
      (Base.$op)(a::AbstractQuantity,b::AbstractQuantity) = call(Base.$op,a,b)
      #(Base.$op)(a::Number,b::AbstractQuantity) = call(Base.$op,GT.uniform_quantity(a),b)
      #(Base.$op)(a::AbstractQuantity,b::Number) = call(Base.$op,a,GT.uniform_quantity(b))
  end
end

# LinearAlgebra

for op in (:inv,:det,:norm,:tr)
  @eval begin
      function get_symbol!(index,val::typeof(LinearAlgebra.$op),name="";prefix=gensym)
          $( Expr(:quote,op) )
      end
    (LinearAlgebra.$op)(a::AbstractQuantity) = call(LinearAlgebra.$op,a)
  end
end

for op in (:dot,:cross)
  @eval begin
      function get_symbol!(index,val::typeof(LinearAlgebra.$op),name="";prefix=gensym)
          $( Expr(:quote,op) )
      end
      (LinearAlgebra.$op)(a::AbstractQuantity,b::AbstractQuantity) = call(LinearAlgebra.$op,a,b)
      #(LinearAlgebra.$op)(a::Number,b::AbstractQuantity) = call(LinearAlgebra.$op,GT.uniform_quantity(a),b)
      #(LinearAlgebra.$op)(a::AbstractQuantity,b::Number) = call(LinearAlgebra.$op,a,GT.uniform_quantity(b))
  end
end



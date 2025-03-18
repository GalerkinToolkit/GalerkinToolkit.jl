
# @term should put every symbol inside
# a LeafNode except the expressions being interpolated,
# which are assumed to to be terms already.
# I have not needed it for the moment.
macro term(expr)
end

function evaluate(expr,captured_data)
    f1 = eval(expr)
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
    dependencies = map(optimize,dependencies(term))
    replace_dependencies(term,dependencies)
end

# We want the leafs to be also AbstractTerm
struct LeafTerm{A,B,C} <: AbstractTerm
    value::A
    prototype::B
    is_compile_constant::Val{C}
end

# The prototye is the same as the value by default
leaf_term(v,prototype=v;is_compile_constant::Val) = LeafTerm(v,v,is_compile_constant)

is_compile_constant(t::LeafTerm) = val_parameter(t.is_compile_constant)

# TODO rename value
value(a::LeafTerm) = a.value
prototype(a::LeafTerm) = a.prototype

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

function prototype(a::CallTerm)
    return_prototype(prototype(a.callee),map(prototype,a.args)...)
end

function dependencies(term::CallTerm)
    (term.callee,calee.args...)
end

function replace_dependencies(term::CallTerm,dependencies)
    (callee,args...) = dependencies
    CallTerm(callee,args)
end

function expression(term::CallTerm)
    callee = expression(term.body)
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
    if haskey(value_arg,value)
        value = value_arg[term.value]
    else
        value = term.value
    end
end

function parametrize!(value_arg,term::AbstractTerm)
    dependencies = map(child->parametrize!(value_arg,child),dependencies(term))
    replace_dependencies(term,dependencies)
end

# captures all values in Leafs that are not compile constants
# and have not been reduced
# Make sure that all indices have been reduced
# before using this function!
function capture(term::AbstractTerm)
    value_arg = IdDict{Any,Symbol}
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
    end
    value_arg[value]
    term
end

function capture!(value_arg,arg,reduced_symbols,term::AbstractTerm)
    bindings = GT.bindings(term)
    append!(reduced_symbols,bindings)
    dependencies = map(child->capture!(value_arg,arg,child),dependencies(term))
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

# The objects that are exposed to the user
# inside of an integral
abstract type AbstractQuantity <: AbstractType end

function quantity(term)
    Quantity(term)
end

# NB.To define a prototype for a quantity  we need a domain
# It is easier to add the prototype to terms instead
struct Quantity{A,B}
    term::A
end

function term(qty::Quantity,opts)
    qyt.term(opts)
end

function call(callee,args::AbstractQuantity...)
    f = compile_constant_quantity(callee)
    call(f,args...)
end

function call(callee::AbstractQuantity,args::AbstractQuantity...)
    quantity() do opts
        f_term = term(f,opts)
        args_term = map(arg->term(arg,opts),args)
        call_term(f_term,args_term...)
    end
end

function compile_constant_quantity(v)
    quantity() do opts
        LeafTerm(v)
    end
end

function coordinate_quantity(quadrature)
    quantity() do opts
        face = domain_face_index(term.index)
        point = point_index(term.index)
        dependencies = map(LeafTerm,(quadrature,face,point))
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
        face = domain_face_index(term.index)
        point = point_index(term.index)
        dependencies = map(LeafTerm,(quadrature,face,point))
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
        the_face_around = GT.face_around(domain)
        index = GT.index(opts)
        face = domain_face_index(index)
        field = field_index(index,arg)
        dof = dof_index(index,arg)
        if D == d
            the_face_around = nothing
            face_around = nothing
            dependencies = map(LeafTerm,(space,domain,face,the_field,field,dof,the_face_around,face_around))
            FormArgumentTerm(dependencies)
        elseif D==d+1 && face_around !== nothing
            the_face_around = face_around
            dependencies = map(LeafTerm,(space,domain,face,the_field,field,dof,the_face_around,face_around))
            FormArgumentTerm(dependencies)
        else
            face_around = face_around_index(index,arg)
            the_face_around = :the_face_around
            dependencies = map(LeafTerm,(space,domain,face,the_field,field,dof,the_face_around,face_around))
            n_faces_around = LeafTerm(2) # Hard coded! But OK in practice.
            SkeletonTerm(FormArgumentTerm(dependencies),n_faces_around,the_face_around)
        end
    end
end

struct FormArgumentTerm{A} <: AbstractTerm
    dependencies::A
end

function prototype(term::FormArgumentTerm)
    (space,domain,face,the_field,field,dof,the_face_around,face_around) = map(prototype,dependencies(term))
    prototype(form_argument_accessor(space,domain,1))
end

function dependencies(term::FormArgumentTerm)
    term.dependencies
end

function replace_dependencies(term::FormArgumentTerm,dependencies)
    FormArgumentTerm(dependencies)
end

function replace_the_face_around(term::FormArgumentTerm,the_face_around)
    (space,domain,face,the_field,field,dof,_,face_around) = term.dependencies
    dependencies = (space,domain,face,the_field,field,dof,the_face_around,face_around)
    FormArgumentTerm(dependencies)
end

function expression(term::FormArgumentTerm)
    (space,domain,face,the_field,field,dof,the_face_around,face_around) = map(expression,term.dependencies)
    :(form_argument_accessor($space,$domain,$the_field)($face,$the_face_around)($dof,$field,$face_around))
end

function optimize(term::CallTerm{<:FormArgumentTerm,Tuple{<:CoordinateTerm}})
    coords = term.args[1]
    (quadrature,face,point) = dependencies(coords)
    TabulatedTerm(term.calee,quadrature,point)
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
        dependencies = map(LeafTerm,(f,uh,domain,face,face_around))
        DiscreteTerm(dependencies)
    elseif D==d+1 && face_around !== nothing
        dependencies = map(LeafTerm,(f,uh,domain,face,face_around))
        DiscreteTerm(dependencies)
    else
        face_around = :the_face_around
        n_faces_around = LeafTerm(2) # Hard coded! But OK in practice.
        dependencies = map(LeafTerm,(f,uh,domain,face,face_around))
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
    (space,domain,face,the_field,field,dof,the_face_around,face_around) = map(expression,form_arg.dependencies)
    :(form_argument_accessor($space,$quadrature,$the_field)($face,$the_face_around)($point)($dof,$field,$face_around))
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

function generate_assemble_vector(contribution::DomainContribution,space::AbstractSpace;parameters=())
    # TODO how to sort the calls to optimize, parametrize, capture, expression, statements, etc
    # needs to be though carefully. We provably will need several calls to optimize and more IR levels
    term_0 = write_assemble_vector(contribution,space)
    term_1 = optimize(term_0)
    term_2 = parametrize(term_1,parameters...)
    term_3, captured_data = capture(term_2)
    expr_0 = expression(term_3)
    expr_1 = statements(expr_0)
    f = evaluate(expr_1,captured_data)
    f
end

function write_assemble_vector(contribution::DomainContribution,space::AbstractSpace)
    quadrature = GT.quadrature(contribution)
    domain = GT.domain(quadrature)
    arity = Val(1)
    index = GT.index(arity)
    term = GT.term(contribution,index)
    space_term = LeafTerm(space)
    quadrature_term = LeafTerm(quadrature)
    arg = :alloc
    alloc = LeafTerm(alloc)
    VectorAssemblyTerm(term,space_term,quadrature_term,alloc,arg,index)
end

struct VectorAssemblyTerm{A,B,C,D,E,F} <: AbstractTerm
    term::A # <: Term
    space::B # <: Term
    quadrature::C # <: Term
    alloc::D
    arg::E # gets reduced
    index::F # gets reduced
end

function bindings(term::VectorAssemblyTerm)
    index = term.index
    face = domain_face_index(index)
    point = point_index(index)
    field = field_index(index,1)
    face_around = face_around_index(index,1)
    dof = dof_index(index,1)
    (term.arg,face,point,field,face_around,dof)
end

function dependencies(term::VectorAssemblyTerm)
    (;term,space,quadrature,alloc) = term
    (term,space,quadrature,alloc)
end

function replace_dependencies(term::VectorAssemblyTerm,dependencies)
    (term,space,quadrature,alloc) = dependencies
    (;arg,index) = term
    VectorAssemblyTerm(term,space,quadrature,alloc,arg,index)
end

# We can specialize for particular cases if needed
function expression(term::VectorAssemblyTerm)
    (term,space,quadrature,alloc) = map(expression,dependencies(term))
    (arg,face,point,field,face_around,dof) = bindings(term)
    dof_v = :($dof -> $term)
    block_dof_v = :( ($field,$face_around) -> $dof_v)
    point_block_dof_v = :($point -> $block_dof_v)
    face_point_block_dof_v = :( $face -> $point_block_dof_v)
    nfaces = :(num_faces($domain))
    body = :(vector_assembly_loop!($face_point_block_dof_v,$alloc,$space,$quadrature))
    expr = :($arg->$body)
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
            zeros(eltype(alloc),m,n)
        end
    end
    z = zero(eltype(alloc))
    for face in 1:nfaces
        point_block_dof_v = face_point_block_dof_v(face)
        field_n_faces_around = face_field_n_faces_around(face)
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
                be = face_around_be[face_around]
                contribute!(alloc,be,dofs,field)
            end
        end
    end
    alloc
end


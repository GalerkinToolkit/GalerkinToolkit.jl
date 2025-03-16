
# @term should put every symbol inside
# a LeafNode except the expressions being interpolated,
# which are assumed to to be terms already.
macro term(expr)
end

abstract type AbstractTerm <: AbstractType end

# By default, we just optimize the children
function optimize(term::AbstractTerm)
    children = map(optimize,children(term))
    replace_children(term,children)
end

# We want the leafs to be also AbstractTerm
struct LeafTerm{A} <: AbstractTerm
    value::A
end

function children(term::LeafTerm)
    (,)
end

function replace_children(term::LeafTerm,children)
    term
end

function expression(term::LeafTerm)
    :( $(term.value) )
end

function lambda_term(body,args...)
    LambdaTerm(body,args)
end

# Should "dummy"  be LeafTerms, or just symbols?
# The answer depends if they should be children or not.
# Children should always be terms
# "the children of a Term are always terms"
# A LeafTerm has no children
struct LambdaTerm{A,B<:Tuple} <: AbstractTerm
    body::A
    args::B # all symbols in here get reduced,
end

# Children are always terms
function children(term::LambdaTerm)
    (term.body,)
end

function replace_children(term::LambdaTerm,children)
    (body,) = children
    LambdaTerm(body,term.args)
end

function expression(term::LambdaTerm)
    args = Expr(:tuple, term.args...)
    body = expression(term.body)
    :($args -> ($body))
end

function call_term(callee,args...)
    CallTerm(callee,args)
end

struct CallTerm{A,B<:Tuple} <: AbstractTerm
    callee::A
    args::B
end

function children(term::CallTerm)
    (callee.body,calee.args...)
end

function replace_children(term::CallTerm,children)
    (callee,args...) = children
    CallTerm(callee,args)
end

function expression(term::CallTerm)
    callee = expression(term.body)
    args = map(expression,term.args)
    :($(callee)($(args...)))
end

# Parametrizes all values in Leafs if they are in params.
function parametrize(term::AbstractTerm,params...)
    pairs = map(enumerate(params)) do (i,param)
        param=>Symbol("param_$(arg[])")
    end
    args = map(last,pairs)
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
    children = map(child->parametrize!(value_arg,child),children(term))
    replace_children(term,children)
end

# captures all values in Leafs
# Make sure that all indices have been reduced
# before using this function!
function capture(term::AbstractTerm)
    value_arg = IdDict{Any,Symbol}
    arg = Ref(0)
    body = capture!(value_arg,arg,term)
    captured_pairs = collect(value_arg)
    sort!(captured_pairs, by=last)
    captured_values = map(first,captured_pairs)
    captured_args = map(last,captured_pairs)
    captured_term = lambda_term(body,captured_args...)
    captured_term, captured_values
end

function capture!(value_arg,arg,term::LeafTerm)
    (;value) = term
    if ! haskey(value_arg,value)
        arg[] += 1
        value_arg[value] = Symbol("captured_arg_$(arg[])")
    end
    value_arg[value]
end

function capture!(value_arg,arg,term::AbstractTerm)
    children = map(child->capture!(value_arg,arg,child),children(term))
    replace_children(term,children)
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
    domain_face_symbol::Symbol
    point_symbol::Symbol
    arg_field_symbol::NTuple{N,Symbol}
    arg_face_around_symbol::NTuple{N,Symbol}
    arg_dof_symbol::NTuple{N,Symbol}
end

domain_face_symbol(index::Index) = index.domain_face_symbol
point_symbol(index::Index) = index.point_symbol
field_symbol(index::Index,arg) = index.arg_field_symbol[arg]
face_around_symbol(index::Index,arg) = index.arg_face_around_symbol[arg]
dof_symbol(index::Index,arg) = index.arg_dof_symbol[arg]

domain_face_term(index::Index) = LeafTerm(index.domain_face_symbol)
point_term(index::Index) = LeafTerm(index.point_symbol)
field_term(index::Index,arg) = LeafTerm(index.arg_field_symbol[arg])
face_around_term(index::Index,arg) = LeafTerm(index.arg_face_around_symbol[arg])
dof_term(index::Index,arg) = LeafTerm(index.arg_dof_symbol[arg])

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

struct Quantity{A}
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
        CoordinateTerm(quadrature,index(opts))
    end
end

struct CoordinateTerm{A,B} <: AbstractTerm
    quadrature::A
    index::B
end

# Children are always terms
function children(term::CoordinateTerm)
    (LeafTerm(term.quadrature),)
end

function replace_children(term::CoordinateTerm,children)
    (quadrature_term,) = children
    quadrature = quadrature_term.value
    CoordinateTerm(quadrature,term.index)
end

# Not clear if this can contain statements or not
# Not clear if this is even needed for DLS terms
function expression(term::CoordinateTerm)
    quadrature = term.quadrature
    face = domain_face_symbol(term.index)
    point = point_symbol(term.index)
    :( coordinate_accessor($quadrature)($face)($point) )
end

function weight_quantity(quadrature)
    quantity() do opts
        WeightTerm(LeafTerm(quadrature),index(opts))
    end
end

# Implement like for CoordinateTerm
struct WeightTermTerm{A,B} <: AbstractTerm
    quadrature::A
    index::B
end

function form_argument_quantity(space::AbstractSpace,arg,the_field=1)
    quantity() do opts
        the_face_around = nothing
        FormArgumentTerm(space,domain(opts),arg,the_field,the_face_around,index(opts))
    end
end

struct FormArgumentTerm{A,B,C,D,E} <: AbstractTerm
    space::A
    domain::B
    arg::C
    the_field::D
    the_face_around::E
end

function children(term::FormArgumentTerm)
    # Dispatch whether or not face_around has been reduced
end

function children_reduced(term::FormArgumentTerm)
    (;space,domain,the_field) = term
    map(LeafTerm,(space,domain,the_field))
end

function children_skeleton(term::FormArgumentTerm)
    (;space,domain,the_field,the_face_around) = term
    map(LeafTerm,(space,domain,the_field))
end


function replace_children(term::FormArgumentTerm,children)
    space, the_field, the_face_around = map(value,children)
    # easy to implement
    error()
end

function expression(term::FormArgumentTerm)
    (;index,arg,space,the_field,the_face_around) = term
    face = domain_face_symbol(index)
    point = point_symbol(index)
    dof = dof_symbol(index,arg)
    field = field_symbol(index,arg)
    face_around = face_around_symbol(index,arg)
    :( form_argument_accessor($space,$the_field)($face,$the_face_around)($dof,$face_around) )
end


function generate_vector_assembly_loop(contribution::DomainContribution,space::AbstractSpace,params...)
    # TODO how to sort the calls to optimize, parametrize, capture, expression, statements, etc
    # needs to be though carefully. We provably will need several calls to optimize and more IR levels
    term_0 = write_vector_assembly_loop(contribution,space)
    term_1 = optimize(term_0)
    term_2 = parametrize(term_1,params...)
    term_3, captured_data = capture(term_2)
    expr_0 = expression(term_3)
    expr_1 = statements(expr_0)
    f = evaluate(expr_1,captured_data)
    f
end

function write_vector_assembly_loop(contribution::DomainContribution,space::AbstractSpace)
    quadrature = GT.quadrature(contribution)
    domain = GT.domain(quadrature)
    arity = Val(1)
    index = GT.index(arity)
    term = GT.term(contribution,index)
    space_term = LeafTerm(space)
    quadrature_term = LeafTerm(quadrature)
    VectorAssemblyLoopTerm(term,index,space_term,quadrature_term)
end

struct VectorAssemblyLoopTerm{A,B,C,D} <: AbstractTerm
    point_value::A # <: Term
    index::B # gets reduced
    space::C # <: Term
    quadrature::D # <: Term
end

function children(term::VectorAssemblyLoopTerm)
    (;point_value,space,quadrature) = term
    (point_value,space,quadrature)
end

function replace_children(term::VectorAssemblyLoopTerm,children)
    (point_value,space,quadrature) = children
    (;index) = term
    VectorAssemblyLoopTerm(term,index,space_term,domain_term)
end


function expression(term::VectorAssemblyLoopTerm)
    index = GT.index(term)
    domain_term = GT.domain_term(term)
    space_term = GT.space_term(term)
    term_expr = expression(term)
    domain_expr = expression(domain_term)
    space_expr = expression(space_term)
    face = index.domain_face
    point = index.point
    field = index.arg_field[1]
    face_around = index.arg_face_around[1]
    dof = index.arg_dof[1]
    dof_v_expr = :($dof -> $term_expr)
    block_dof_v_expr = :( ($field,$face_around) -> $dof_v_expr)
    point_block_dof_v_expr = :($point -> $block_dof_v_expr)
    face_point_block_dof_v_expr = :( $face -> $point_block_dof_v_expr)
    nfaces_expr = :(num_faces($domain_expr))
    alloc_expr = :alloc
    block_be_expr = :block_be
    body = :(vector_assembly_loop!($face_point_block_dof_v_expr,$alloc_expr,$block_be_expr,$space_expr,$quadrature_expr))
    expr = :(($alloc_expr,$be_expr)->$body)
    expr
end

function vector_assembly_loop!(face_point_block_dof_v,alloc,block_be,space,quadrature)
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
    z = zero(eltype(eltype(block_be)))
    for face in 1:nfaces
        point_block_dof_v = face_point_block_dof_v(face)
        field_n_faces_around = face_field_n_faces_around(face)
        # Reset
        for be in block_be
            fill!(be,z)
        end
        # Integrate
        npoints = face_npoints(face)
        for point in 1:npoints
            block_dof_v = point_block_dof_v(point)
            for field in 1:nfields
                n_faces_around = field_face_n_faces_around[field](face)
                for face_around in 1:n_faces_around
                    be = block_be[field,face_around]
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
            n_faces_around = field_face_n_faces_around[field](face)
            for face_around in 1:n_faces_around
                be = block_be[field,face_around]
                contribute!(alloc,be,dofs,field)
            end
        end
    end
    alloc
end


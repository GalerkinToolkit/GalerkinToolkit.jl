#macro term(expr)
#    transform(a) = a # :($(Expr(:quote,a)))
#    function transform(expr::Expr)
#        if  expr.head === :call
#            transform_call(expr)
#        elseif  expr.head === :ref
#            transform_ref(expr)
#        else
#            transform_default(expr)
#        end
#    end
#    function transform_call(expr::Expr)
#        args = map(transform,expr.args)
#        quote
#            Expr(:call,$(args...))
#        end
#    end
#    function transform_ref(expr::Expr)
#        args = map(transform,expr.args)
#        quote
#            Expr(:ref,$(args...))
#        end
#    end
#    function transform_default(expr::Expr)
#        head = expr.head
#        args = map(transform,expr.args)
#        Expr(head,args...)
#    end
#    transform(expr) |> esc
#end

macro term(expr)
    vars = Set{Symbol}()
    function findvars!(a)
        nothing
    end
    function findvars!(expr::Expr)
        if  expr.head === :(=)
            var = expr.args[1]
            push!(vars,var)
        # else
        #     map(findvars!,expr.args)
        end
        map(findvars!,expr.args) # always traverse the args. If the user write (a = b = 3) then only a is captured
        nothing
    end
    transform(a) = a
    function transform(a::Symbol)
        if a in vars
            a
        else
            :($(Expr(:quote,a)))
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
        elseif  expr.head === :do
            transform_do(expr)
        elseif  expr.head === :block
            transform_default(expr)
        elseif  expr.head === :tuple
            transform_tuple(expr)
        elseif  expr.head === :.
            transform_dot(expr)
        else
            error("Expr with head=$(expr.head) not supported in macro @term")
        end
    end
    function transform_call(expr::Expr)
        args = map(transform,expr.args)
        quote
            Expr(:call,$(args...))
        end
    end
    function transform_ref(expr::Expr)
        args = map(transform,expr.args)
        quote
            Expr(:ref,$(args...))
        end
    end
    function transform_lambda(expr::Expr)
        args = map(transform,expr.args)
        ret = quote
            Expr(:(->),$(args...)) # TODO: check the behavior of transform_lambda. it doesn't work when there is a lambda in terms if written as $(args...)
        end
    end
    function transform_do(expr::Expr)
        expr2 = Expr(:call,expr.args[1].args[1],expr.args[2],expr.args[1].args[2])
        transform(expr2)
    end
    function transform_interpolation(expr::Expr)
        expr.args[1]
    end
    function transform_tuple(expr::Expr)
        args = map(transform,expr.args)
        quote
            Expr(:tuple,$(args...))
        end
    end
    function transform_dot(expr::Expr)
        expr
    end
    function transform_default(expr::Expr)
        head = expr.head
        args = map(transform,expr.args)
        Expr(head,args...)
    end
    findvars!(expr)
    transform(expr) |> esc
end

index(t::AbstractTerm) = t.index
#num_dims(t::AbstractTerm) = t.dim
free_dims(t::AbstractTerm) = t.free_dims
prototype(t::AbstractTerm) = t.prototype

function call(f,args::AbstractTerm...)
    if length(args) == 1
        unary_call_term(f,args...)
    elseif length(args) == 2
        binary_call_term(f,args...)
    else
        multivar_call_term(f,args)
    end
end

function (f::AbstractTerm)(x::AbstractTerm)
    call(call,f,x)
end

for op in (:+,:-,:*,:/,:\,:^)
  @eval begin
      (Base.$op)(a::AbstractTerm,b::AbstractTerm) = call(Base.$op,a,b)
      (Base.$op)(a::Number,b::AbstractTerm) = call(Base.$op,GT.constant_term(a),b)
      (Base.$op)(a::AbstractTerm,b::Number) = call(Base.$op,a,GT.constant_term(b))
  end
end

constant_term(args...;compile_constant=false) = ConstantTerm(args...,compile_constant)

struct ConstantTerm{A,B,C} <: AbstractTerm
    value::A
    index::B
    compile_constant::C
end

AbstractTrees.children(a::ConstantTerm) = (a.value,)

function expression(a::ConstantTerm)
    if a.compile_constant
        :($(a.value))
    else
        get_symbol!(index(a),a.value,"constant_quantity_value")
    end
end

free_dims(a::ConstantTerm) = Int[]

prototype(t::ConstantTerm) = t.value

reference_face_term(args...) = ReferenceFaceTerm(args...)

#struct ReferenceFaceTerm{A,B,C} <: AbstractTerm
#    dim::A
#    rid_to_value::B
#    index::C
#end
#
#physical_face_term(args...) = PhysicalFaceTerm(args...)
#
#struct PhysicalFaceTerm{A,B,C,D} <: AbstractTerm
#    dim::A
#    sface_to_value::B
#    face_to_sface::C
#    index::D
#end
#
#AbstractTrees.children(a::PhysicalFaceTerm) =
#(
# "dim = $(a.dim)",
# "sface_to_value = $(summary(a.sface_to_value))",
# "face_to_sface = $(summary(a.face_to_sface))",
# "face = $(face_index(a.index,a.dim))",
#)
#
#function expression(a::PhysicalFaceTerm)
#    face = face_index(a.index,a.dim)
#    sface_to_value = get_symbol!(a.index,a.sface_to_value,"sface_to_value_$(a.dim)d")
#    face_to_sface = get_symbol!(a.index,a.face_to_sface,"face_to_sface_$(a.dim)d")
#    @term $sface_to_value[$face_to_sface[$face]]
#end
#
#function prototype(a::PhysicalFaceTerm)
#    zero(eltype(a.sface_to_value))
#end

unary_call_term(args...) = UnaryCallTerm(args...)

struct UnaryCallTerm{A,B} <: AbstractTerm
   callee::A
   arg::B
end

AbstractTrees.children(a::UnaryCallTerm) = (a.callee,a.arg)

index(a::UnaryCallTerm) = index(a.arg)
function prototype(a::UnaryCallTerm)
    return_prototype(a.callee,prototype(a.arg))
end

free_dims(a::UnaryCallTerm) = free_dims(a.arg)


function expression(c::UnaryCallTerm)
    f = c.callee
    a = c.arg
    arg = expression(a)
    # ndofs = expression(a.ndofs)
    # dof = a.dof
    expr = @term begin
        $f($arg)
        # fun = $dof -> $c*$s
        # sum(fun,1:$ndofs)
    end
end

binary_call_term(args...) = BinaryCallTerm(args...)

struct BinaryCallTerm{A,B,C} <: AbstractTerm
    callee::A
    arg1::B
    arg2::C
end

AbstractTrees.children(a::BinaryCallTerm) = (a.callee,a.arg1,a.arg2)

index(a::BinaryCallTerm) = index(a.arg1)

function expression(a::BinaryCallTerm)
    expr1 = expression(a.arg1)
    expr2 = expression(a.arg2)
    g = get_symbol!(index(a),a.callee,"callee")
    :($g($expr1,$expr2))
end

function prototype(a::BinaryCallTerm)
    return_prototype(a.callee,prototype(a.arg1),prototype(a.arg2))
end

free_dims(a::BinaryCallTerm) = union(free_dims(a.arg1),free_dims(a.arg2))

#function num_dims(a::BinaryCallTerm)
#    dim = promote_dim(num_dims(a.arg1),num_dims(a.arg2))
#    if dim !== nothing
#        return dim
#    end
#    min(num_dims(a.arg1),num_dims(a.arg2))
#end

multivar_call_term(args...) = MultivarCallTerm(args...)

struct MultivarCallTerm{A,B} <: AbstractTerm
    callee::A
    args::B
end

AbstractTrees.children(a::MultivarCallTerm) = (a.callee,a.args...)

index(a::MultivarCallTerm) = index(a.args[1])

function expression(a::MultivarCallTerm)
    exprs = map(expression,a.args)
    g = get_symbol!(index(a),a.callee,"callee")
    :($g($(exprs...)))
end

function prototype(a::MultivarCallTerm)
    ps = map(prototype,a.args)
    return_prototype(a.callee,ps...)
end

free_dims(a::MultivarCallTerm) = union(map(free_dims,a.args)...)


boundary_term(args...) = BoundaryTerm(args...)

struct BoundaryTerm{A,B,C,D} <: AbstractTerm
    dim::A
    dim_around::B
    term::C
    face_around::D
end

AbstractTrees.children(a::BoundaryTerm) =
(
 "dim = $(a.dim)","dim_around = $(a.dim_around)",a.face_around,a.term)

index(a::BoundaryTerm) = index(a.term)

function expression(a::BoundaryTerm)
    expr = expression(a.term)
    topo = topology(mesh(index(a)))
    dim = a.dim
    dim2 = a.dim_around
    face_to_cells = get_symbol!(index(a),face_incidence(topo,dim,dim2),"face_to_cells")
    cell_around = expression(a.face_around)
    face = face_index(index(a),dim)
    actual = :($face_to_cells[$face][$cell_around])
    dummy = face_index(index(a),dim2)
    expr = substitute(expr,dummy=>actual)
    dummy = face_around_dummy_index(index(a),dim,dim2)
    substitute(expr,dummy=>cell_around)
end

free_dims(a::BoundaryTerm) = setdiff(free_dims(a.term),[a.dim_around])

prototype(a::BoundaryTerm) = prototype(a.term)

function binary_call_term(f,a,b::BoundaryTerm)
    t = binary_call_term(f,a,b.term)
    boundary_term(b.dim,b.dim_around,t,b.face_around)
end

function binary_call_term(f,a::BoundaryTerm,b)
    t = binary_call_term(f,a.term,b)
    boundary_term(a.dim,a.dim_around,t,a.face_around)
end

function binary_call_term(f,a::BoundaryTerm,b::BoundaryTerm)
    BinaryCallTerm(f,a,b)
end

skeleton_term(args...) = SkeletonTerm(args...)

function skeleton_term(dim,dim_around,term::BoundaryTerm)
    if dim == term.dim && dim_around == term.dim_around
        term
    else
        SkeletonTerm(dim,dim_around,term)
    end
end

struct SkeletonTerm{A,B,C} <: AbstractTerm
    dim::A
    dim_around::B
    term::C
end

AbstractTrees.children(a::SkeletonTerm) = ("dim = $(a.dim)","dim_around = $(a.dim_around)",a.term)

index(a::SkeletonTerm) = index(a.term)

prototype(a::SkeletonTerm) = [prototype(a.term),prototype(a.term)]

free_dims(a::SkeletonTerm) = setdiff(free_dims(a.term),[a.dim_around])

function expression(a::SkeletonTerm)
    expr = expression(a.term)
    topo = topology(mesh(index(a)))
    dim = a.dim
    dim2 = a.dim_around
    face_to_cells = get_symbol!(index(a),face_incidence(topo,dim,dim2),"face_to_cells")
    cell_dummy = face_index(index(a),dim2)
    cell_around = face_around_dummy_index(index(a),dim,dim2)
    face = face_index(index(a),dim)
    cells = :($face_to_cells[$face])
    cell_actuall = :($face_to_cells[$face][$cell_around])
    expr = substitute(expr,cell_dummy=>cell_actuall)
    @term begin
        fun = $cell_around -> $expr
        map(fun,1:length($cells))
    end
    #axis_to_face_around = face_around_index(index(a))
    #if length(axis_to_face_around) == 0
    #    return @term begin
    #        fun = $cell_around -> $expr
    #        map(fun,1:length($cells))
    #    end
    #end
    #exprs = map(axis_to_face_around) do face_around
    #    :(LinearAlgebra.I[$i,$face_around])
    #end
    #if length(exprs) > 1
    #    deltas = Expr(:call,:tuple,exprs...)
    #    mask = :(prod($deltas))
    #else
    #    mask = exprs[1]
    #end
    #@term begin
    #    fun = $cell_around -> $expr*$mask
    #    map(fun,1:length($cells))
    #end
end

function binary_call_term(::typeof(getindex),a::SkeletonTerm,b::ConstantTerm)
    t = a.term
    boundary_term(a.dim,a.dim_around,t,b)
    #if form_arity(index(a)) == 0
    #    t = a.term
    #elseif form_arity(index(a)) == 1
    #    face_around = face_around_index(index(a))[1]
    #    c = expr_term(Int[],face_around,0,index(a))
    #    mask = call(==,c,b)
    #    t = mask*a.term
    #elseif form_arity(index(a)) == 2
    #    face_around_1, face_around_2 = face_around_index(index(a))
    #    c1 = expr_term(Int[],face_around_1,0,index(a))
    #    c2 = expr_term(Int[],face_around_2,0,index(a))
    #    mask1 = call(==,c1,b)
    #    mask2 = call(==,c2,b)
    #    t = (mask1*mask2)*a.term
    #else
    #    error("case not implemented, but possible to implement it")
    #end
    #boundary_term(a.dim,a.dim_around,t,b)
end

expr_term(args...;dof=nothing) = ExprTerm(args...,dof)

struct ExprTerm{A,B,C,D,E} <: AbstractTerm
    free_dims::A
    expr::B
    prototype::C
    index::D
    dof::E # TODO ugly. Hack for dual_basis(a::AbstractSpace) + sign flip in RT space
end

dof_index(a::ExprTerm) = a.dof

AbstractTrees.children(a::ExprTerm) = ("dim = $(a.free_dims)",)

expression(a::ExprTerm) = a.expr

#user_getindex_term(args...) = UserGetindexTerm(args...)
#
#function user_getindex_term(a::SkeletonTerm,array_index)
#    boundary_term(a.dim,a.dim_around,a.term,array_index)
#end

#struct UserGetindexTerm{A,B} <: AbstractTerm
#    parent::A
#    array_index::B
#end
#
#num_dims(a::UserGetindexTerm) = num_dims(a.parent)
#
#function expression(a::UserGetindexTerm)
#    expr = expression(a.parent)
#    i = expression(a.array_index)
#    :($expr[$i])
#end
#
#AbstractTrees.children(a::UserGetindexTerm) = (a.parent,a.array_index)

reference_point_term(args...) = ReferencePointTerm(args...)

struct ReferencePointTerm{A,B,C,D} <: AbstractTerm
    dim::A
    rid_to_point_to_value::B
    face_to_rid::C
    index::D
end

free_dims(a::ReferencePointTerm) = [a.dim]

AbstractTrees.children(a::ReferencePointTerm) =
(
 "dim = $(a.dim)",
 "rid_to_point_to_value = $(summary(a.rid_to_point_to_value))",
 "face_to_rid = $(summary(a.face_to_rid))",
)

function expression(a::ReferencePointTerm)
    face = face_index(a.index,a.dim)
    point = point_index(a.index)
    rid_to_point_to_value = get_symbol!(a.index,a.rid_to_point_to_value,"rid_to_point_to_value_$(a.dim)d")
    face_to_rid = get_symbol!(a.index,a.face_to_rid,"face_to_rid_$(a.dim)d")
    @term $rid_to_point_to_value[$face_to_rid[$face]][$point]
end

function prototype(a::ReferencePointTerm)
    zero(eltype(eltype(a.rid_to_point_to_value)))
end

physical_point_term(args...) = PhysicalPointTerm(args...)

struct PhysicalPointTerm{A,B,C,D} <: AbstractTerm
    dim::A
    sface_to_point_to_value::B
    face_to_sface::C
    index::D
end

AbstractTrees.children(a::PhysicalPointTerm) =
(
 "dim = $(a.dim)",
 "sface_to_point_to_value = $(summary(a.sface_to_point_to_value))",
 "face_to_sface = $(summary(a.face_to_sface))",
)

function expression(a::PhysicalPointTerm)
    face = face_index(a.index,a.dim)
    point = point_index(a.index)
    sface_to_point_to_value = get_symbol!(a.index,a.sface_to_point_to_value,"sface_to_point_to_value_$(a.dim)d")
    face_to_sface = get_symbol!(a.index,a.face_to_sface,"face_to_sface_$(a.dim)d")
    @term $sface_to_point_to_value[$face_to_sface[$face]][$point]
end

function prototype(a::PhysicalPointTerm)
    zero(eltype(eltype(a.sface_to_point_to_value)))
end

#struct PhysicalMapTerm{A,B} <: AbstractTerm
#    dim::A
#    index::B
#end
#
#reference_map_term(args...) = ReferenceMapTerm(args...)
#
#struct ReferenceMapTerm{A,B,C} <: AbstractTerm
#    dim::A
#    dim_around::B
#    index::C
#end
#
#inverse_map_term(args...) = InverseMapTerm(args...)
#
#struct InverseMapTerm{A} <: AbstractTerm
#    direct_map::A
#end

#reference_shape_function_term(args...) = ReferenceShapeFunctionTerm(args...)

#struct ReferenceShapeFunctionTerm{A,B,C,D,E,F} <: AbstractTerm
#    dim::A
#    rid_to_reffe::B
#    face_to_rid::C
#    caller::D
#    dof::E # dummy
#    index::F
#end
#
#AbstractTrees.children(a::ReferenceShapeFunctionTerm) =
#(
# "dim = $(a.dim)",
# "rid_to_reffe = $(summary(a.rid_to_reffe))",
# "face_to_rid = $(summary(a.face_to_rid))",
# "caller =$(a.caller)",
# "face = $(face_index(a.index,a.dim))",
#)
#
#function expression(a::ReferenceShapeFunctionTerm)
#    face = face_index(a.index,a.dim)
#    dof = a.dof
#    rid_to_dof_to_fun = get_symbol!(a.index,map(shape_functions,a.rid_to_reffe),"rid_to_dof_to_fun")
#    face_to_rid = get_symbol!(a.index,a.face_to_rid,"face_to_rid_$(a.dim)d")
#    @term x -> $caller($rid_to_dof_to_fun[$face_to_rid[$face]][$dof],x)
#end


reference_shape_function_term(args...) = ReferenceShapeFunctionTerm(args...)

struct ReferenceShapeFunctionTerm{A,B,C,D,E} <: AbstractTerm
    dim::A
    rid_to_dof_to_value::B
    face_to_rid::C
    dof::D
    index::E
end

free_dims(a::ReferenceShapeFunctionTerm) = [a.dim]

AbstractTrees.children(a::ReferenceShapeFunctionTerm) =
(
 "dim = $(a.dim)",
 "rid_to_dof_to_value = $(summary(a.rid_to_dof_to_value))",
 "face_to_rid = $(summary(a.face_to_rid))",
)

function expression(a::ReferenceShapeFunctionTerm)
    face = face_index(a.index,a.dim)
    dof = a.dof
    rid_to_dof_to_value = get_symbol!(a.index,a.rid_to_dof_to_value,"rid_to_point_to_value_$(a.dim)d")
    face_to_rid = get_symbol!(a.index,a.face_to_rid,"face_to_rid_$(a.dim)d")
    @term $rid_to_dof_to_value[$face_to_rid[$face]][$dof]
end

function prototype(a::ReferenceShapeFunctionTerm)
    a.rid_to_dof_to_value |> first |> first
end

function binary_call_term(f,a::ReferenceShapeFunctionTerm,b::ReferencePointTerm)
    tab = map(b.rid_to_point_to_value) do point_to_x
        map(a.rid_to_dof_to_value) do dof_to_f
            f.(dof_to_f,permutedims(point_to_x))
        end
    end
    tabulator_term(f,tab,a,b)
    #ridb_rida_tab = get_symbol!(index(a),tab,"tabulator")
    ##a_rid_to_dof_to_value = get_symbol!(a.index,a.rid_to_dof_to_value,"rid_to_dof_to_value_$(a.dim)d")
    #a_face_to_rid = get_symbol!(a.index,a.face_to_rid,"face_to_rid_$(a.dim)d")
    #b_rid_to_point_to_value = get_symbol!(b.index,b.rid_to_point_to_value,"rid_to_point_to_value_$(b.dim)d")
    #b_face_to_rid = get_symbol!(b.index,b.face_to_rid,"face_to_rid_$(b.dim)d")
    #face_a = face_index(index(a),a.dim)
    #face_b = face_index(index(b),b.dim)
    #dof = a.dof
    #point = point_index(index(b))
    #expr = @term begin
    #    rid_a = $a_face_to_rid[$face_a]
    #    rid_b = $b_face_to_rid[$face_b]
    #    $ridb_rida_tab[rid_b][rid_a][$dof,$point]
    #end
    #dims = union(free_dims(a),free_dims(b))
    #p = prototype(a)(prototype(b))
    #expr_term(dims,expr,p,index)
end

tabulator_term(args...) = TabulatorTerm(args...)

struct TabulatorTerm{A,B,C,D} <: AbstractTerm
    callee::A
    ridb_rida_tabulator::B
    arg1::C
    arg2::D
end

AbstractTrees.children(a::TabulatorTerm) =
(
 a.callee,
 a.arg1,
 a.arg2,
)

free_dims(t::TabulatorTerm) = union(free_dims(t.arg1),free_dims(t.arg2))

index(t::TabulatorTerm) = index(t.arg1)

prototype(a::TabulatorTerm) = a.ridb_rida_tabulator |> eltype |> eltype |> zero

function expression(t::TabulatorTerm)
    a = t.arg1
    b = t.arg2
    ridb_rida_tab = get_symbol!(index(a),t.ridb_rida_tabulator,"tabulator")
    #a_rid_to_dof_to_value = get_symbol!(a.index,a.rid_to_dof_to_value,"rid_to_dof_to_value_$(a.dim)d")
    a_face_to_rid = get_symbol!(a.index,a.face_to_rid,"face_to_rid_$(a.dim)d")
    b_rid_to_point_to_value = get_symbol!(b.index,b.rid_to_point_to_value,"rid_to_point_to_value_$(b.dim)d")
    b_face_to_rid = get_symbol!(b.index,b.face_to_rid,"face_to_rid_$(b.dim)d")
    face_a = face_index(index(a),a.dim)
    face_b = face_index(index(b),b.dim)
    dof = a.dof
    point = point_index(index(b))
    expr = @term begin
        rid_a = $a_face_to_rid[$face_a]
        rid_b = $b_face_to_rid[$face_b]
        $ridb_rida_tab[rid_b][rid_a][$dof,$point]
    end
end

form_argument_term(args...) = FormArgumentTerm(args...)

struct FormArgumentTerm{A,B,C,D} <: AbstractTerm
    axis::A
    field::B
    functions::C
    evaluated::D
end

function dof_index(t::FormArgumentTerm)
    dof_index(t.functions)
end

function dof_index(t::TabulatorTerm)
    dof_index(t.arg1)
end

function dof_index(t::ReferenceShapeFunctionTerm)
    t.dof
end

AbstractTrees.children(a::FormArgumentTerm) =
(
 "axis = $(a.axis)",
 "field = $(a.field)",
 a.functions
)

index(a::FormArgumentTerm) = index(a.functions)

free_dims(a::FormArgumentTerm) = free_dims(a.functions)

prototype(a::FormArgumentTerm) = prototype(a.functions)

function expression(c::FormArgumentTerm)
    @assert c.evaluated
    form_argument_expression(c,c.functions)
end

function form_argument_expression(c,::TabulatorTerm)
    field = field_index(index(c),c.axis)
    dof = dof_index(index(c),c.axis)
    t = c.functions
    a = t.arg1
    b = t.arg2
    ridb_rida_tab = get_symbol!(index(a),t.ridb_rida_tabulator,"tabulator")
    #a_rid_to_dof_to_value = get_symbol!(a.index,a.rid_to_dof_to_value,"rid_to_dof_to_value_$(a.dim)d")
    a_face_to_rid = get_symbol!(a.index,a.face_to_rid,"face_to_rid_$(a.dim)d")
    b_rid_to_point_to_value = get_symbol!(b.index,b.rid_to_point_to_value,"rid_to_point_to_value_$(b.dim)d")
    b_face_to_rid = get_symbol!(b.index,b.face_to_rid,"face_to_rid_$(b.dim)d")
    face_a = face_index(index(a),a.dim)
    face_b = face_index(index(b),b.dim)
    dof = a.dof
    point = point_index(index(b))
    # TODO if field is a integer instead of a symbol we can return a much more optimized expression
    # Idem for face_around
    D = a.dim
    d = target_dim(index(a))
    if d == D || is_boundary(domain(index(a)))
        expr = @term begin
            rid_a = $a_face_to_rid[$face_a]
            rid_b = $b_face_to_rid[$face_b]
            tab = $ridb_rida_tab[rid_b][rid_a]
            field_bool = $(c.field) == $field
            GalerkinToolkit.getindex_if(field_bool,tab,$dof,$point)
        end
    else
        face_around = face_around_index(index(a),c.axis)
        face_around_dummy = face_around_dummy_index(index(a),d,D)
        expr = @term begin
            rid_a = $a_face_to_rid[$face_a]
            rid_b = $b_face_to_rid[$face_b]
            tab = $ridb_rida_tab[rid_b][rid_a]
            field_bool = $(c.field) == $field
            geom_bool = $face_around == $face_around_dummy
            bool = field_bool & geom_bool
            GalerkinToolkit.getindex_if(bool,tab,$dof,$point)
        end
    end
end

function getindex_if(b,a,i...;otherwise=zero(eltype(a)))
    if b
        a[i...]
    else
        otherwise
    end
end

function binary_call_term(f,a::FormArgumentTerm,b::ReferencePointTerm)
    t = binary_call_term(f,a.functions,b)
    evaluated = true
    form_argument_term(a.axis,a.field,t,evaluated)
end

discrete_function_term(args...) = DiscreteFunctionTerm(args...)

struct DiscreteFunctionTerm{A,B,C,D} <: AbstractTerm
    coefficients::A
    functions::B
    dof::C
    ndofs::D
end

AbstractTrees.children(a::DiscreteFunctionTerm) =
(
 a.coefficients,
 a.functions,
 a.dof,
 a.ndofs,
)

function free_dims(a::DiscreteFunctionTerm)
    union(free_dims(a.coefficients),free_dims(a.functions),free_dims(a.ndofs))
end

function expression(a::DiscreteFunctionTerm)
    s = expression(a.functions)
    c = expression(a.coefficients)
    dof = a.dof
    ndofs = expression(a.ndofs)
    @term begin
        x -> begin
            g = $dof -> $c*$s(x)
            n = $ndofs
            sum(g,1:n)
        end
    end
end

function prototype(a::DiscreteFunctionTerm)
    s = prototype(a.functions)
    c = prototype(a.coefficients)
    x -> c*s(x) + c*s(x)
end

index(a::DiscreteFunctionTerm) = index(a.functions)

const LinearOperators = Union{typeof(call),typeof(ForwardDiff.gradient),typeof(ForwardDiff.jacobian)}

function sum_init(f,r,s)
    sum(f,r;init=s)
end

function expression(c::BinaryCallTerm{<:LinearOperators,<:DiscreteFunctionTerm,<:Any})
    f = c.callee
    a = c.arg1
    b = c.arg2
    s0 = zero(prototype(c))
    s0_sym = get_symbol!(index(c),s0,"s0")
    s = expression(call(f,a.functions,b))
    v = expression(a.coefficients)
    ndofs = expression(a.ndofs)
    dof = a.dof
    expr = :(GT.sum_init($dof -> $v*$s,1:$ndofs,$(s0_sym)))
end

#function binary_call_term(f::LinearOperators,a::DiscreteFunctionTerm,b::ReferencePointTerm)
#    s = expression(call(f,a.functions,b))
#    c = expression(a.coefficients)
#    ndofs = expression(a.ndofs)
#    dof = a.dof
#    expr = @term begin
#        fun = $dof -> $c*$s
#        sum(fun,1:$ndofs)
#    end
#    dims = union(free_dims(a),free_dims(b))
#    px = prototype(b)
#    pf = prototype(a)
#    p = pf(px)
#    expr_term(dims,expr,p,index(a))
#end

discrete_function_call_term(args...) = DiscreteFunctionCallTerm(args...)

physical_map_term(args...) = PhysicalMapTerm(args...)

struct PhysicalMapTerm{A,B} <: AbstractTerm
    dim::A
    index::B
end

free_dims(a::PhysicalMapTerm) = [a.dim]


function prototype(a::PhysicalMapTerm)
    D = num_ambient_dims(mesh(index(a)))
    x-> zero(SVector{D,Float64})
end


function Base.:(==)(a::PhysicalMapTerm,b::PhysicalMapTerm)
    a.dim == b.dim
end

function discrete_function_term(a::PhysicalMapTerm)
    dof = gensym("dummy-physical-dof")
    d = a.dim
    face_to_rid = face_reference_id(mesh(index(a)),d)
    rid_to_dof_to_fun = map(shape_functions,reference_spaces(mesh(index(a)),d))
    functions = reference_shape_function_term(d,rid_to_dof_to_fun,face_to_rid,dof,index(a))
    g_to_v = node_coordinates(mesh(index(a)))
    g_to_value = get_symbol!(index(a),g_to_v,"g_to_value")
    face_to_dof_to_g = get_symbol!(index(a),face_nodes(mesh(index(a)),d),"face_to_dof_to_g")
    face = GT.face_index(index(a),d)
    expr = :($g_to_value[$face_to_dof_to_g[$face][$dof]])
    p = zero(eltype(g_to_v))
    coeffs = GT.expr_term([d],expr,p,index(a))
    expr = :(length($face_to_dof_to_g[$face]))
    ndofs = expr_term([d],expr,0,index(a))
    discrete_function_term(coeffs,functions,dof,ndofs)
end

function expression(a::PhysicalMapTerm)
    expression(discrete_function_term(a))
end

#function binary_call_term(f::typeof(call),a::PhysicalMapTerm,b::ReferencePointTerm)
#    a2 = discrete_function_term(a)
#    call(f,a2,b)
#end
#
#function binary_call_term(f::typeof(ForwardDiff.jacobian),a0::PhysicalMapTerm,b::ReferencePointTerm)
#    a = discrete_function_term(a0)
#    s = expression(call(ForwardDiff.gradient,a.functions,b))
#    c = expression(a.coefficients)
#    ndofs = expression(a.ndofs)
#    dof = a.dof
#    expr = @term begin
#        fun = $dof -> $c*transpose($s)
#        sum(fun,1:$ndofs)
#    end
#    dims = union(free_dims(a),free_dims(b))
#    px = prototype(b)
#    pf = prototype(a)
#    p = pf(px)
#    expr_term(dims,expr,p,index(a))
#end

function expression(c::BinaryCallTerm{typeof(call),<:PhysicalMapTerm,<:Any})
    f = c.callee
    a = discrete_function_term(c.arg1)
    b = c.arg2
    fab = call(f,a,b)
    expression(fab)
end

function expression(c::BinaryCallTerm{typeof(ForwardDiff.jacobian),<:PhysicalMapTerm,<:Any})
    f = c.callee
    a = discrete_function_term(c.arg1)
    b = c.arg2
    s = expression(call(ForwardDiff.gradient,a.functions,b))
    c = expression(a.coefficients)
    ndofs = expression(a.ndofs)
    dof = a.dof
    expr = @term begin
        fun = $dof -> $c*transpose($s)
        sum(fun,1:$ndofs)
    end
end

function inv_map(f,x0)
    function pseudo_inverse_if_not_square(J)
        m,n = size(J)
        if m != n
            pinv(J)
        else
            inv(J)
        end
    end
    function invf(fx)
        x = x0
        tol = 1.0e-12
        J = nothing
        niters = 100
        for _ in 1:niters
            J = ForwardDiff.jacobian(f,x)
            Jinv = pseudo_inverse_if_not_square(J)
            dx = Jinv*(fx-f(x))
            x += dx
            if norm(dx) < tol
                return x
            end
        end
        error("Max iterations reached")
        x
    end
end

function return_prototype(::typeof(inv_map),f,x0)
    fx -> x0
end

function get_symbol!(index,val::typeof(inv_map),name="";prefix=index.data.prefix)
    :(GalerkinToolkit.inv_map)
end

const ComposedWithInverseTerm{A,B,C} = BinaryCallTerm{typeof(∘),A,BinaryCallTerm{typeof(inv_map),B,C}}

const FunctionCallTerm{A,B} = BinaryCallTerm{typeof(call),A,B}

function binary_call_term(::typeof(call),a::ComposedWithInverseTerm,b::FunctionCallTerm)
    f = a.arg1
    x = b.arg2
    phi1 = a.arg2.arg1
    phi2 = b.arg1
    if phi1 != phi2
        return binary_call_term_physical_maps(call,a,b,phi1,phi2)
    end
    f(x)
end

function binary_call_term(d::typeof(ForwardDiff.gradient),a::ComposedWithInverseTerm,b::FunctionCallTerm)
    phi1 = a.arg2.arg1
    phi2 = b.arg1
    if phi1 != phi2
        return binary_call_term_physical_maps(d,a,b,phi1,phi2)
    end
    phi = phi1
    f = a.arg1
    x = b.arg2
    J = call(ForwardDiff.jacobian,phi,x)
    Jt = call(transpose,J)
    gradf = call(d,f,x)
    call(\,Jt,gradf)
end

function binary_call_term(d::typeof(ForwardDiff.jacobian),a::ComposedWithInverseTerm,b::FunctionCallTerm)
    phi1 = a.arg2.arg1
    phi2 = b.arg1
    if phi1 != phi2
        return binary_call_term_physical_maps(d,a,b,phi1,phi2)
    end
    phi = phi1
    f = a.arg1
    x = b.arg2
    J = call(ForwardDiff.jacobian,phi,x)
    Jf = call(d,f,x)
    call(/,Jf,J)
end

#  w(q) = u(ϕ(q))
#  ∇(w,q) = Jt(ϕ,q)*∇(u,ϕ(q))
#  Jt(ϕ,q)\∇(w,q) = ∇(u,ϕ(q))
#
#
#  J(w,q) = J(u,ϕ(q))*J(ϕ,q)
#  J(w,q)/J(ϕ,q) = J(u,ϕ(q))

function binary_call_term_physical_maps(op,a,b,phi1,phi2)
    return BinaryCallTerm(op,a,b)
end

function binary_call_term_physical_maps(::typeof(call),a,b,phi1::PhysicalMapTerm,phi2::PhysicalMapTerm)
    @assert phi1.dim != phi2.dim
    phi = reference_map_term(phi2.dim,phi1.dim,index(a))
    f = a.arg1
    x = b.arg2
    f(phi(x))
end

function binary_call_term_physical_maps(d::typeof(ForwardDiff.gradient),a,b,phi1::PhysicalMapTerm,phi2::PhysicalMapTerm)
    @assert phi1.dim != phi2.dim
    phi = reference_map_term(phi2.dim,phi1.dim,index(a))
    f = a.arg1
    x = b.arg2
    J = call(ForwardDiff.jacobian,phi1,phi(x))
    Jt = call(transpose,J)
    gradf = call(d,f,phi(x))
    call(\,Jt,gradf)
end

reference_map_term(args...) = ReferenceMapTerm(args...)

struct ReferenceMapTerm{A,B,C} <: AbstractTerm
    dim::A
    dim_around::B
    index::C
end

AbstractTrees.children(a::ReferenceMapTerm) =
(
 "dim = $(a.dim)",
 "dim_around = $(a.dim_around)",
)


free_dims(a::ReferenceMapTerm) = ([a.dim,a.dim_around])

prototype(a::ReferenceMapTerm) = x-> zero(SVector{a.dim_around,Float64})

#function binary_call_term(f::typeof(call),a::ReferenceMapTerm,b::ReferencePointTerm)
#    a2 = discrete_function_term(a)
#    call(f,a2,b)
#end

function expression(c::BinaryCallTerm{typeof(call),<:ReferenceMapTerm,<:Any})
    f = c.callee
    a = discrete_function_term(c.arg1)
    b = c.arg2
    fab = call(f,a,b)
    expression(fab)
end

function discrete_function_term(a::ReferenceMapTerm)
    dof = gensym("dummy-reference-dof")
    d = a.dim
    D = a.dim_around
    dface_to_drid = face_reference_id(mesh(index(a)),d)
    drid_to_dof_to_fun = map(shape_functions,reference_spaces(mesh(index(a)),d))
    functions = reference_shape_function_term(d,drid_to_dof_to_fun,dface_to_drid,dof,index(a))
    dface_to_drid = get_symbol!(index(a),face_reference_id(mesh(index(a)),d),"dface_to_drid")
    dface = face_index(index(a),d)
    drid_to_ndofs = get_symbol!(index(a),map(length,drid_to_dof_to_fun),"drid_to_ndofs")
    expr = :($drid_to_ndofs[$dface_to_drid[$dface]])
    ndofs = expr_term([d],expr,0,index(a))
    Dface = face_index(index(a),D)
    topo = topology(mesh(index(a)))
    dface_to_Dfaces_data, dface_to_ldfaces_data = GT.face_incidence_ext(topo,d,D)
    dface_to_Dfaces = get_symbol!(index(a),dface_to_Dfaces_data,"dface_to_Dfaces")
    dface_to_ldfaces = get_symbol!(index(a),dface_to_ldfaces_data,"dface_to_ldfaces")
    Dface_to_ldface_to_perm = get_symbol!(index(a),GT.face_permutation_ids(topo,D,d),"Dface_to_ldface_to_perm")
    Dface_to_Drid = get_symbol!(index(a),face_reference_id(mesh(index(a)),D),"Dface_to_Drid")
    Drid_to_refDface = reference_spaces(mesh(index(a)),D)
    Drid_to_ldface_perm_dof_to_x_data = map(refDface->face_node_coordinates(refDface,d),Drid_to_refDface)
    Drid_to_ldface_perm_dof_to_x = get_symbol!(index(a),Drid_to_ldface_perm_dof_to_x_data,"Drid_to_ldface_perm_dof_to_x")
    Dface_around = face_around_dummy_index(index(a),d,D)
    expr = @term begin
        drid = $dface_to_drid[$dface]
        Drid = $Dface_to_Drid[$Dface]
        ldface = $dface_to_ldfaces[$dface][$Dface_around]
        perm = $Dface_to_ldface_to_perm[$Dface][ldface]
        $Drid_to_ldface_perm_dof_to_x[Drid][ldface][perm][$dof]
    end
    p = Drid_to_ldface_perm_dof_to_x_data |> eltype |> eltype |> eltype |> eltype |> zero
    coeffs = expr_term([d,D],expr,p,index)
    discrete_function_term(coeffs,functions,dof,ndofs)
end

function binary_call_term(f,a::ReferenceShapeFunctionTerm,b::BinaryCallTerm{typeof(call),<:ReferenceMapTerm,<:ReferencePointTerm})
    phi = b.arg1
    x = b.arg2
    d = phi.dim
    D = phi.dim_around
    m = mesh(index(a))
    topo = topology(m)
    drid_to_refdface = reference_spaces(m,d)
    Drid_to_refDface = reference_spaces(m,D)
    # NB the TODOs below can be solved by introducing two extra nesting levels
    # TODO this assumes the same reffes for mesh and quadrature
    drid_Drid_ldface_perm_to_dof_and_point_data = map(x.rid_to_point_to_value,drid_to_refdface) do point_to_x,refdface
        # TODO this assumes the same reffes for mesh and interpolation
        map(a.rid_to_dof_to_value,Drid_to_refDface) do dof_to_f,refDface
            ldface_perm_varphi = reference_map(refdface,refDface)
            map(ldface_perm_varphi) do perm_varphi
                map(perm_varphi) do varphi
                    point_to_q = varphi.(point_to_x)
                    f.(dof_to_f,permutedims(point_to_q))
                end
            end
        end
    end
    face_tabulator_term(f,drid_Drid_ldface_perm_to_dof_and_point_data,a,b)
    #drid_Drid_ldface_perm_to_dof_and_point = get_symbol!(index(a),drid_Drid_ldface_perm_to_dof_and_point_data,"face_tabulator")
    #Drid_to_dof_to_value = get_symbol!(a.index,a.rid_to_dof_to_value,"Drid_to_dof_to_f")
    #Dface_to_rid = get_symbol!(a.index,a.face_to_rid,"Dface_to_Drid")
    #dface_to_Dfaces_data, dface_to_ldfaces_data = GT.face_incidence_ext(topo,d,D)
    #dface_to_Dfaces = get_symbol!(index(a),dface_to_Dfaces_data,"dface_to_Dfaces")
    #dface_to_ldfaces = get_symbol!(index(a),dface_to_ldfaces_data,"dface_to_ldfaces")
    #Dface_to_ldface_to_perm = get_symbol!(index(a),GT.face_permutation_ids(topo,D,d),"Dface_to_ldface_to_perm")
    #Dface_to_Drid = get_symbol!(index(a),face_reference_id(mesh(index(a)),D),"Dface_to_Drid")
    #dface_to_drid = get_symbol!(index(a),face_reference_id(mesh(index(a)),d),"dface_to_drid")
    #dface = face_index(index(a),d)
    #Dface = face_index(index(b),D)
    #dof = a.dof
    #point = point_index(index(b))
    #Dface_around = face_around_dummy_index(index(a),d,D)
    #expr = @term begin
    #    drid = $dface_to_drid[$dface]
    #    Drid = $Dface_to_Drid[$Dface]
    #    Dfaces = $dface_to_Dfaces[$dface]
    #    ldface = $dface_to_ldfaces[$dface][$Dface_around]
    #    perm = $Dface_to_ldface_to_perm[$Dface][ldface]
    #    $drid_Drid_ldface_perm_to_dof_and_point[drid][Drid][ldface][perm][$dof,$point]
    #end
    #dims = union(free_dims(a),free_dims(b))
    #p = f(prototype(a),prototype(b))
    #expr_term(dims,expr,p,index(a))
end

function binary_call_term(f,a::FormArgumentTerm,b::BinaryCallTerm{typeof(call),<:ReferenceMapTerm,<:ReferencePointTerm})
    t = binary_call_term(f,a.functions,b)
    form_argument_term(a.axis,a.field,t,true)
end

face_tabulator_term(args...) = FaceTabulatorTerm(args...)

struct FaceTabulatorTerm{A,B,C,D} <: AbstractTerm
    callee::A
    drid_Drid_ldface_perm_to_dof_and_point_data::B
    arg1::C
    arg2::D
end

AbstractTrees.children(a::FaceTabulatorTerm) =
(
 a.callee,
 a.arg1,
 a.arg2,
)

free_dims(t::FaceTabulatorTerm) = union(free_dims(t.arg1),free_dims(t.arg2))

index(t::FaceTabulatorTerm) = index(t.arg1)

prototype(a::FaceTabulatorTerm) = a.drid_Drid_ldface_perm_to_dof_and_point_data |> eltype |> eltype |> eltype |> eltype |> zero

function expression(t::FaceTabulatorTerm)
    a = t.arg1
    b = t.arg2
    phi = b.arg1
    x = b.arg2
    d = phi.dim
    D = phi.dim_around
    m = mesh(index(a))
    topo = topology(m)
    drid_to_refdface = reference_spaces(m,d)
    Drid_to_refDface = reference_spaces(m,D)
    drid_Drid_ldface_perm_to_dof_and_point = get_symbol!(index(a),t.drid_Drid_ldface_perm_to_dof_and_point_data,"face_tabulator")
    Drid_to_dof_to_value = get_symbol!(a.index,a.rid_to_dof_to_value,"Drid_to_dof_to_f")
    Dface_to_rid = get_symbol!(a.index,a.face_to_rid,"Dface_to_Drid")
    dface_to_Dfaces_data, dface_to_ldfaces_data = GT.face_incidence_ext(topo,d,D)
    dface_to_Dfaces = get_symbol!(index(a),dface_to_Dfaces_data,"dface_to_Dfaces")
    dface_to_ldfaces = get_symbol!(index(a),dface_to_ldfaces_data,"dface_to_ldfaces")
    Dface_to_ldface_to_perm = get_symbol!(index(a),GT.face_permutation_ids(topo,D,d),"Dface_to_ldface_to_perm")
    Dface_to_Drid = get_symbol!(index(a),face_reference_id(mesh(index(a)),D),"Dface_to_Drid")
    dface_to_drid = get_symbol!(index(a),face_reference_id(mesh(index(a)),d),"dface_to_drid")
    dface = face_index(index(a),d)
    Dface = face_index(index(b),D)
    dof = a.dof
    point = point_index(index(b))
    Dface_around = face_around_dummy_index(index(a),d,D)
    @term begin
        drid = $dface_to_drid[$dface]
        Drid = $Dface_to_Drid[$Dface]
        ldface = $dface_to_ldfaces[$dface][$Dface_around]
        perm = $Dface_to_ldface_to_perm[$Dface][ldface]
        $drid_Drid_ldface_perm_to_dof_and_point[drid][Drid][ldface][perm][$dof,$point]
    end
end

function form_argument_expression(c,::FaceTabulatorTerm)
    field = field_index(index(c),c.axis)
    dof = dof_index(index(c),c.axis)
    t = c.functions
    a = t.arg1
    b = t.arg2
    phi = b.arg1
    x = b.arg2
    d = phi.dim
    D = phi.dim_around
    m = mesh(index(a))
    topo = topology(m)
    drid_to_refdface = reference_spaces(m,d)
    Drid_to_refDface = reference_spaces(m,D)
    drid_Drid_ldface_perm_to_dof_and_point = get_symbol!(index(a),t.drid_Drid_ldface_perm_to_dof_and_point_data,"face_tabulator")
    Drid_to_dof_to_value = get_symbol!(a.index,a.rid_to_dof_to_value,"Drid_to_dof_to_f")
    Dface_to_rid = get_symbol!(a.index,a.face_to_rid,"Dface_to_Drid")
    dface_to_Dfaces_data, dface_to_ldfaces_data = GT.face_incidence_ext(topo,d,D)
    dface_to_Dfaces = get_symbol!(index(a),dface_to_Dfaces_data,"dface_to_Dfaces")
    dface_to_ldfaces = get_symbol!(index(a),dface_to_ldfaces_data,"dface_to_ldfaces")
    Dface_to_ldface_to_perm = get_symbol!(index(a),GT.face_permutation_ids(topo,D,d),"Dface_to_ldface_to_perm")
    Dface_to_Drid = get_symbol!(index(a),face_reference_id(mesh(index(a)),D),"Dface_to_Drid")
    dface_to_drid = get_symbol!(index(a),face_reference_id(mesh(index(a)),d),"dface_to_drid")
    dface = face_index(index(a),d)
    Dface = face_index(index(b),D)
    dof = a.dof
    point = point_index(index(b))
    Dface_around = face_around_dummy_index(index(a),d,D)
    D = a.dim
    d = target_dim(index(a))
    if is_boundary(domain(index(a)))
        expr = @term begin
            drid = $dface_to_drid[$dface]
            Drid = $Dface_to_Drid[$Dface]
            Dfaces = $dface_to_Dfaces[$dface]
            ldface = $dface_to_ldfaces[$dface][$Dface_around]
            perm = $Dface_to_ldface_to_perm[$Dface][ldface]
            tab = $drid_Drid_ldface_perm_to_dof_and_point[drid][Drid][ldface][perm]
            field_bool = $(c.field) == $field
            GalerkinToolkit.getindex_if(field_bool,tab,$dof,$point)
        end
    else
        face_around = face_around_index(index(a),c.axis)
        face_around_dummy = face_around_dummy_index(index(a),d,D)
        expr = @term begin
            drid = $dface_to_drid[$dface]
            Drid = $Dface_to_Drid[$Dface]
            Dfaces = $dface_to_Dfaces[$dface]
            ldface = $dface_to_ldfaces[$dface][$Dface_around]
            perm = $Dface_to_ldface_to_perm[$Dface][ldface]
            tab = $drid_Drid_ldface_perm_to_dof_and_point[drid][Drid][ldface][perm]
            field_bool = $(c.field) == $field
            geom_bool = $face_around == $face_around_dummy
            bool = field_bool & geom_bool
            GalerkinToolkit.getindex_if(bool,tab,$dof,$point)
        end
    end
end


#function find_face_adound(face,faces)
#    findfirst(i->i==face,faces)
#end

unit_normal_term(args...) = UnitNormalTerm(args...)

struct UnitNormalTerm{A} <: AbstractTerm
    outwards_normals::A
end

free_dims(a::UnitNormalTerm) = free_dims(a.outwards_normals)
prototype(a::UnitNormalTerm) = x -> prototype(a.outwards_normals)
index(a::UnitNormalTerm) = index(a.outwards_normals)

function expression(c::BinaryCallTerm{typeof(call),<:UnitNormalTerm,<:Any})
    f = c.callee
    D = num_dims(mesh(index(c)))
    a = c.arg1
    x = c.arg2
    phiD = physical_map_term(D,index(c))
    J = call(ForwardDiff.jacobian,phiD,x)
    n = call(map_unit_normal,J,a.outwards_normals)
    expression(n)
end

function map_unit_normal(J,n)
    Jt = transpose(J)
    pinvJt = transpose(inv(Jt*J)*Jt)
    v = pinvJt*n
    m = sqrt(v⋅v)
    if m < eps()
        return zero(v)
    else
        return v/m
    end
end

function get_symbol!(index,val::typeof(map_unit_normal),name="";prefix=index.data.prefix)
    :(GalerkinToolkit.map_unit_normal)
end

#discrete_function_term(args...) = DiscreteFunctionTerm(args...)

#struct DiscreteFunctionTerm{A,B} <: AbstractTerm
#    coefficients::A
#    functions::B
#    ndofs::C
#end
#
#function expression(a::DiscreteFunctionTerm)
#    s = expression(a.functions)
#    c = expression(a.coefficients)
#    f = c*s
#    dof = a.functions.dof # TODO it assumes that it matches with coefficents
#    ndofs = num_dofs(a.functions)
#    @term begin
#        g = $dof -> $f
#        x -> sum(g,$ndof)
#    end
#end
#
#face_dof_term(args...) = FaceDofTerm(args...)
#
#struct FaceDofTerm{A,B,C,D,E} <: AbstractTerm
#    dim::A
#    dof_to_value::B
#    face_to_dofs::C
#    dof::D
#    index::E
#end


function topological_sort(expr,deps)
    temporary = gensym()
    expr_L = [Expr(:block) for _ in 0:length(deps)]
    marks = Dict{UInt,Symbol}()
    marks_deps = Dict{UInt,Int}()
    function visit(expr_n::Symbol)
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
        args = expr_n.args[2:end]
        r = map(visit,args)
        args_var = map(first,r)
        j = maximum(map(last,r))
        expr_n_new = Expr(expr_n.head,expr_n.args[1],args_var...)
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

#const gradient = ForwardDiff.gradient
#const jacobian = ForwardDiff.jacobian
#const hessian = ForwardDiff.hessian
#function simplify(expr::Expr)
#    # r_invmap needs to be processed before r_tabulator
#    # We use call since @rule is not able to identify function calls on variable slots (phi in this case)
#    r_invmap = @slots a b c phi @rule call( a ∘ inverse_map_impl(phi,b) , call(phi,c) ) --> call(a,c)
#
#    # $(\nabla (f\circ\varphi^{-1}))(\varphi(x))$ -> $(J(\varphi(x)))^{-T} (\nabla f)(x)  $
#    # matrix inversion and jacobian may introduce floating point error
#    r_invmap_gradient = @slots a b c phi @rule gradient(a ∘ inverse_map_impl(phi,b), call(phi,c))  -->  inv(jacobian(phi, c))' * gradient(a, c)
#
#    # we assume to have the same type of reference face (b[e] == h[i])
#    r_tabulator = @slots a b c d e f g h i @rule call(face_function(a,b,c,d,e),reference_value(f,h,i)[g]) --> face_function_value(map(reference_tabulator,a,f),b,c,d,e,g)
#    r_jacobian_tabulator = @slots a b c d e f g h i @rule jacobian(face_function(a,b,c,d,e),reference_value(f,h,i)[g]) --> jacobian_face_function_value(map(gradient_reference_tabulator,a,f),b,c,d,e,g)
#    
#
#    r_shape_gradient_tabulator = @slots rid_to_fs face_to_rid face dof rid_to_coords face_to_rid2 sface point @rule gradient(face_shape_function(rid_to_fs,face_to_rid,face,dof), reference_value(rid_to_coords, face_to_rid2, sface)[point]) --> 
#        face_shape_function_value(map(gradient_reference_tabulator, rid_to_fs, rid_to_coords), face_to_rid, face, dof, face_to_rid2, sface, point)
#
#    r_shape_tabulator = @slots rid_to_fs face_to_rid face dof rid_to_coords face_to_rid2 sface point @rule call(face_shape_function(rid_to_fs,face_to_rid,face,dof), reference_value(rid_to_coords, face_to_rid2, sface)[point]) --> 
#        face_shape_function_value(map(reference_tabulator, rid_to_fs, rid_to_coords), face_to_rid, face, dof, face_to_rid2, sface, point)
#
#    
#    rules = [
#             r_invmap_gradient,
#             r_invmap,
#             r_tabulator,
#             r_jacobian_tabulator,
#             r_shape_tabulator,
#             r_shape_gradient_tabulator,
#            ]
#
#    expr2 = nothing
#    for rule in rules
#        expr2 = rule(expr)
#        if expr2 !== nothing
#            break
#        end
#    end
#    if expr2 === nothing
#        args = map(simplify,expr.args)
#        Expr(expr.head,args...)
#    else
#        simplify(expr2)
#    end
#end

function simplify(expr)
    expr
end
#
#function simplify(expr::Expr)
#    expr = rewrite(expr)
#    expr = inline_lambdas(expr)
#end

function rewrite(expr)
    expr
end

function rewrite(expr::Expr)
    #r_map_getindex_enumerate = @slots f r i @rule getindex(map(f,enumerate(r)),i) --> call(f,(i,r[i]))
    r_map_getindex_1 = @slots f r i @rule getindex(map(f,r),i) --> call(f,r[i])
    r_map_getindex_2 = @slots f r s i @rule getindex(map(f,r,s),i) --> call(f,r[i],s[i])
    rules = [
             #r_map_getindex_enumerate,
             r_map_getindex_1,
             r_map_getindex_2,
            ]

    expr2 = nothing
    for rule in rules
        expr2 = rule(expr)
        if expr2 !== nothing
            break
        end
    end
    if expr2 === nothing
        args = map(rewrite,expr.args)
        Expr(expr.head,args...)
    else
        rewrite(expr2)
    end
end

function inline_lambdas(expr)
    MacroTools.postwalk(expr) do ex
        if !isa(ex,Expr)
            return ex
        end
        if ex.head !== :call
            return ex
        end
        if ex.args[1] !== :call && ex.args[1] != :(GalerkinToolkit.call)
            return ex
        end
        lambda = ex.args[2]
        if !isa(lambda,Expr)
            return ex
        end
        if lambda.head !== :(->)
            return ex
        end
        signature, body = lambda.args
        if isa(signature,Symbol)
            x = ex.args[3]
            return substitute(body,signature=>x)
        end
        if !isa(signature,Expr)
            return ex
        end
        if signature.head === :tuple
            actual_x = ex.args[3:end]
            dummy_x = signature.args
            for i in 1:length(dummy_x)
                if isa(dummy_x[i],Symbol)
                    body = substitute(body,dummy_x[i]=>actual_x[i])
                end
            end
            return body
        end
        ex
    end
end

function substitute(expr,old_new)
    old,new = old_new
    MacroTools.postwalk(ex -> ex === old ? new : ex ,expr)
end


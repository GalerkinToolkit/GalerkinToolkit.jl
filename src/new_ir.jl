abstract type NewAbstractTerm <: AbstractType end

"""
This is used to create the argument of functions converting quantities to terms
"""
function term_options(domain;domain_face=:domain_face)
    contents = (;domain,domain_face)
    TermOptions(contents)
end

struct TermOptions{A}
    contents::A
end

domain(t::TermOptions) = t.contents.domain
domain_face(t::TermOptions) = t.contents.domain_face

abstract type NewAbstractQuantity <: AbstractType end

function new_quantity(term)
    NewQuantity(term)
end

struct NewQuantity{A} <: NewAbstractQuantity
    term::A
end

function term(q::NewQuantity,opt::TermOptions)
    q.term(opt)
end

function uniform_quantity(value)
    new_quantity() do opts
        UniformTerm(value)
    end
end

struct UniformTerm{A} <: NewAbstractTerm
    value::A
end

function generate_body(t::UniformTerm,param_to_symbol,domain_face)
    param_to_symbol[t.value]
end

struct NumPointsTerm{A,B} <: NewAbstractTerm
    quadrature::A
    domain_face::B
end

function generate_body(t::NumPointsTerm,param_to_symbol,domain_face)
    quadrature = param_to_symbol[t.quadrature]
    @term begin
        face_npoints = GT.num_points_accessor($quadrature)
        face_npoints($domain_face)
    end
end

struct WeightTerm{A,B,C} <: NewAbstractTerm
    quadrature::A
    domain_face::B
    point::C
end

function generate_body(t::WeightTerm,param_to_symbol,domain_face)
    quadrature = param_to_symbol[t.quadrature]
    point = t.point
    @term begin
        face_point_J = GT.jacobian_accessor($quadrature)
        face_point_w = GT.weight_accessor($quadrature)
        J = face_point_J($domain_face)($point)
        face_point_w($domain_face)($point,J)
    end
end

function coordinate_quantity(quadrature::AbstractQuadrature,point)
    new_quantity() do opts
        domain_face = GT.domain_face(opts)
        CoordinateTerm(quadrature,domain_face,point)
    end
end

struct CoordinateTerm{A,B,C} <: NewAbstractTerm
    quadrature::A
    domain_face::B
    point::C
end

function contribution(f,quadrature::AbstractQuadrature)
    domain_face = :domain_face
    point = :point
    x = coordinate_quantity(quadrature,point)
    fx = f(x)
    weight = WeightTerm(quadrature,domain_face,point)
    num_points = NumPointsTerm(quadrature,domain_face)
    domain = GT.domain(weight.quadrature)
    opts = term_options(domain;domain_face)
    integrand = term(fx,opts)
    ContributionTerm(
                     integrand,
                     weight,
                     domain_face,
                     point,
                     num_points, )
end

struct ContributionTerm{A,B,C,D,E} <: NewAbstractTerm
    integrand::A
    weight::B
    domain_face::C # remains free. Not sure if we need it
    point::D # gets reduced
    num_points::E
end

function generate_body(t::ContributionTerm,param_to_symbol,domain_face)
    integrand_expr = generate_body(t.integrand,param_to_symbol,domain_face)
    weight_expr = generate_body(t.weight,param_to_symbol,domain_face)
    num_points_expr = generate_body(t.num_points,param_to_symbol,domain_face)
    point = t.point
    quote
        sum($point -> $integrand_expr * $weight_expr,  1:$num_points_expr)
    end
end

function generate(t::NewAbstractTerm,params...)
    itr = [ p=>Symbol("arg$i") for (i,p) in enumerate(params)  ]
    @show symbols = map(last,itr)
    param_to_symbol = IdDict(itr )
    display(param_to_symbol)
    domain_face = :dummy_domain_face
    expr = generate_body(t,param_to_symbol,domain_face)
    quote
        ($(symbols...),) -> begin
            $domain_face -> begin
                $expr
            end
        end
    end
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

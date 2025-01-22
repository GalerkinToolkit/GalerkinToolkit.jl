
abstract type NewAbstractTerm <: AbstractType end

"""
Creates a term representing a placeholder for a (face,point,dof) index.
"""
function index_term(name="index";prefix=gensym,int_type::Type{T}) where T
    IndexTerm(T,prefix(name))
end

struct IndexTerm{T} <: NewAbstractTerm
    int_type::Type{T}
    name::Symbol
end

function expression!(runtime,t::IndexTerm)
    t.name
end

"""

This is used to create the argument of functions converting quantities to terms
"""
function term_options(;domain,domain_face=index_term(:domain_face))
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

function new_constant_quantity(value)
    new_quantity() do opts
        NewConstantTerm(value)
    end
end

struct NewConstantTerm{A} <: NewAbstractTerm
    value::A
end

function expression!(runtime,t::NewConstantTerm)
    value = compile_symbol!(runtime,t.value)
    value
end

#function call(f,args::NewAbstractQuantity...)
#    new_quantity() do opts
#        args_term = map(arg->term(arg,opts),args)
#        CallTerm(f,args_term)
#    end
#end
#
#struct CallTerm{A,B} <: NewAbstractTerm
#    callee::A
#    args::B
#end

#function (f::NewAbstractQuantity)(x::NewAbstractQuantity)
#    new_quantity() do opts
#        x_term = term(x,opts)
#        f_term = term(f,opts)
#        EvaluateTerm(f_term,x_term)
#    end
#end
#
#struct EvaluateTerm{A,B} <: NewAbstractTerm
#    callee::A
#    arg::B
#end

function new_coordinate_quantity(q::AbstractQuadrature,point::NewAbstractTerm)
    new_quantity() do opts
        @assert domain(q) == domain(opts)
        mesh_face = MeshFaceTerm(domain,domain_face)
        CoordinateTerm(q,point,mesh_face)
    end
end

struct CoordinateTerm{A,B,C} <: NewAbstractTerm
    quadrature::A
    point::B
    mesh_face::C
end

function expression!(runtime,t::CoordinateTerm)
    q = t.quadrature
    rid_point_coordinate_data = map(coordinates,reference_quadratures(q))
    mesh_face_rid_data = face_reference_id(q)
    rid_point_coordinate = compile_symbol!(runtime,rid_point_coordinate_data)
    mesh_face_rid = compile_symbol!(runtime,mesh_face_rid_data)
    point = expression!(runtime,t.point)
    mesh_face = expression!(runtime,t.mesh_face)
    mesh = GT.mesh(GT.domain(q))
    if is_physical_domain(t.domain)
        d = num_dims(GT.domain(q))
        rid_tab_data = map(rid_point_coordinate_data,reference_spaces(mesh,d)) do point_to_x, refface
            tabulator(refface)(value,point_to_x)
        end
        rid_tab = compile_symbol!(runtime,rid_tab_data)
        node_x = compile_symbol!(runtime,node_coordinates(mesh))
        mesh_face_nodes = compile_symbol!(runtime,face_nodes(mesh,d))
        @term begin
            rid = $mesh_face_reference_id[$mesh_face]
            tab = $rid_tab[rid]
            dof_node = $mesh_face_nodes[mesh_face]
            ndofs = length(dof_node)
            sum(1:ndofs) do dof
                node = dof_node[dof]
                x = $node_x[node]
                tab[$point,dof]*x
            end
        end
    else
        @term begin
            $rid_point_coordinate[$mesh_face_rid[$mesh_face]]
        end
    end
end

struct MeshFaceTerm{A,B} <: NewAbstractTerm
    domain::A
    domain_face::B
end

function expression!(runtime,t::MeshFaceTerm)
    domain_face_to_mesh_face = compile_symbol!(runtime,faces(t.domain))
    domain_face = expression!(runtime,t.domain_face)
    quote
        $domain_face_to_mesh_face[$domain_face]
    end
end

function new_weight_quantity(q::AbstractQuadrature,point::NewAbstractTerm)
    new_quantity() do opts
        @assert domain(q) == domain(opts)
        mesh_face = MeshFaceTerm(domain,domain_face)
        WeightTerm(q,point,mesh_face)
    end
end

struct WeightTerm{A,B,C} <: NewAbstractTerm
    quadrature::A
    point::B
    mesh_face::C
end

function expression!(runtime,t::WeightTerm)
    q = t.quadrature
    rid_point_coordinate_data = map(coordinates,reference_quadratures(q))
    mesh_face_rid_data = face_reference_id(q)
    rid_point_weight = compile_symbol!(runtime,map(weights,reference_quadratures(q)))
    mesh_face_rid = compile_symbol!(runtime,mesh_face_rid_data)
    point = expression!(runtime,t.point)
    mesh_face = expression!(runtime,t.mesh_face)
    mesh = GT.mesh(GT.domain(q))
    if is_physical_domain(t.domain)
        d = num_dims(GT.domain(q))
        rid_tab_data = map(rid_point_coordinate_data,reference_spaces(mesh,d)) do point_to_x, refface
            tabulator(refface)(ForwardDiff.gradient,point_to_x)
        end
        rid_tab = compile_symbol!(runtime,rid_tab_data)
        node_x = compile_symbol!(runtime,node_coordinates(mesh))
        mesh_face_nodes = compile_symbol!(runtime,face_nodes(mesh,d))
        @term begin
            rid = $mesh_face_reference_id[$mesh_face]
            tab = $rid_tab[rid]
            dof_node = $mesh_face_nodes[mesh_face]
            ndofs = length(dof_node)
            J = sum(1:ndofs) do dof
                node = dof_node[dof]
                x = $node_x[node]
                outer(x,tab[$point,dof])
            end
            w = $rid_point_weight[rid][$point]
            w*change_of_measure(J)
        end
    else
        @term begin
            $rid_point_weight[$mesh_face_rid[$mesh_face]]
        end
    end
end

#function new_physical_map_term(domain::AbstractDomain,mesh_face::NewAbstractTerm)
#    mesh = GT.mesh(domain)
#    num_dims = GT.num_dims(domain)
#    NewPhysicalMapTerm(mesh,num_dims,mesh_face)
#end
#
#struct NewPhysicalMapTerm{A,B,C} <: NewAbstractTerm
#    mesh::A
#    num_dims::B
#    mesh_face::C
#end
#
#function expression!(runtime,t::NewPhysicalMapTerm)
#    mesh = t.mesh
#    d = t.num_dims
#    mesh_face_reference_id = compile_symbol!(runtime,face_reference_id(mesh,d))
#    rid_dof_shape_fun = compile_symbol!(runtime,map(shape_functions,reference_spaces(mesh,d)))
#    mesh_face_nodes = compile_symbol!(runtime,face_nodes(mesh,d))
#    node_x = compile_symbol!(runtime,node_coordinates(mesh))
#    @term begin
#        rid = $mesh_face_reference_id[$mesh_face]
#        dof_shape_fun = $rid_dof_shape_fun[rid]
#        dof_node = $mesh_face_nodes[mesh_face]
#        ndofs = length(dof_node)
#        y -> sum(1:ndofs) do dof
#            shape_fun = dof_shape_fun[dof]
#            node = dof_node[dof]
#            x = $node_x[node]
#            shape_fun(y)*x
#        end
#    end
#end

function integrate(f,measure::NewMeasure)
    point = index_term(:point)
    x = new_coordinate_quantity(measure,point)
    dV = new_weight_quantity(measure,point)
    fx = f(x)
    domain = GT.domain(measure)
    domain_face = index_term(:domain_face)
    opts = term_options(domain;domain_face)
    dV_term = term(dV,opts)
    fx_term = term(fx,opts)
    contribution = IntegralTerm(fx_term,dV_term,domain_face,point)
    contributions = (contribution,)
    NewIntegral(contributions)
end

struct IntegralTerm{A,B,C,D,E,F} <: NewAbstractTerm
    measure::A
    integrand::B
    weight::C
    domain_face::D # remains free
    point::E # gets reduced
    num_points::F
end

function expression(t::NewAbstractTerm)
    runtime = IdDict{Any,Symbol}()
    expr = expression!(runtime,t)
    expr, runtime
end

function expression!(runtime,t::IntegralTerm)
    domain_face = t.domain_face
    point = t.point
    fx = expression!(runtime,t.integrand)
    w = expression!(runtime,t.weight)
    npoints = expression!(runtime,t.num_points)
    quote
        ($domain_face,arg) -> begin
            $(unpack_runtime_argument(runtime,:arg))
            sum($point->$fx*$w,1:$npoints)
        end
    end
end

struct NewIntegral{A}  <: AbstractType
    contributions::A
end
contributions(a::NewIntegral) = a.contributions

function sum(int::NewIntegral)
    f = generate_sum(int)
    f()
end

function generate_sum(int::NewIntegral,params...)
    fs = map(c->generate_sum(c,params...),contributions(int))
    params2... -> begin
        sum(f->f(params2...),fs)
    end
end

function generate_sum(ir0::IntegralTerm,params::NewAbstractQuantity...)
    @assert length(params) == 0
    #ir1 = optimize(ir0)
    #ir2 = lower(ir1)
    #ir3 = optimize()
    #domain = GT.domain(ir0)
    #domain_face = ir0.domain_face
    #opts = term_options(domain;domain_face)
    #params_terms = map(param -> term(params,opts),params)
    expr, runtime = expression(ir0)
    f = eval(expr)
    arg = runtime_argument(runtime)
    nfaces = num_faces(domain(measure(ir0)))
    () -> sum(domain_face -> invokelatest(f,domain_face,arg), 1:nfaces)
end

function runtime_argument(runtime)
    (;( key=>val for (val,key) in runtime )...)
end

function unpack_runtime_argument(runtime,arg)
    expr = Expr(:block)
    for k in Base.values(index.data.dict) |> collect |> sort
        push!(expr.args,:($k = $arg.$k))
    end
    expr
end

function compile_symbol!(runtime,val,name="";prefix=gensym)
    get!(runtime,val,prefix(name))
end

#int = ∫(α,dΩ)
#f = generate_sum(int,α)
#f(α)
#
#int = ∫(dummy_α,dΩ)
#f = generate_sum(int,dummy_α)
#f(α)
#
# Let us go with this
#int = α -> ∫(α,dΩ)
#f = generate_sum(int,α_prototype)
#f(α)

# and this
#int = () -> ∫(α,dΩ)
#f = generate_sum(int)
#f()

#int = ∫(α,dΩ)
#f = generate_sum(int)
#f()




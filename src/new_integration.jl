
function measure(dom::AbstractDomain,degree)
    measure(quadrature,dom,degree)
end

function measure(f,dom::AbstractDomain,degree)
    mesh = GT.mesh(dom)
    Measure(mesh,f,dom,degree)
end

struct Measure{A,B,C,D} <: GT.AbstractType
    mesh::A
    quadrature_rule::B
    domain::C
    degree::D
end

domain(a::Measure) = a.domain
function PartitionedArrays.partition(a::Measure{<:PMesh})
    map(partition(a.domain)) do domain
        GT.measure(a.quadrature_rule,domain,a.degree)
    end
end

function reference_quadratures(measure::Measure)
    domain = GT.domain(measure)
    mesh = GT.mesh(domain)
    d = GT.num_dims(domain)
    drefid_refdface = GT.reference_spaces(mesh,d)
    refid_to_quad = map(drefid_refdface) do refdface
        geo = GT.domain(refdface)
        measure.quadrature_rule(geo,measure.degree)
    end
    refid_to_quad
end

function face_reference_id(m::Measure)
    domain = GT.domain(m)
    mesh = GT.mesh(domain)
    d = GT.num_dims(domain)
    face_reference_id(mesh,d)
end

function coordinates(measure::Measure)
    domain = GT.domain(measure)
    coordinates(measure,domain)
end

#function coordinates(measure::Measure{<:PMesh})
#    q = map(GT.coordinates,partition(measure))
#    term = map(GT.term,q)
#    GT.quantity(term) 
#end

function coordinates(measure::Measure,domain::ReferenceDomain)
    mesh = GT.mesh(domain)
    d = GT.num_dims(domain)
    face_to_refid = GT.face_reference_id(mesh,d)
    domface_to_face = GT.faces(domain)
    refid_to_quad = reference_quadratures(measure)
    refid_to_coords = map(GT.coordinates,refid_to_quad)
    prototype = first(GT.coordinates(first(refid_to_quad)))

    point_quantity(refid_to_coords,domain;reference=true)
end

function coordinates(measure::Measure,domain::PhysicalDomain)
    Ω = domain
    Ωref = GT.reference_domain(Ω)
    ϕ = GT.physical_map(mesh(Ω),num_dims(Ω))
    x = GT.coordinates(measure,Ωref)
    ϕ(x)
end

function weights(measure::Measure)
    domain = GT.domain(measure)
    weights(measure,domain)
end

#function weights(measure::Measure{<:PMesh})
#    q = map(GT.weights,partition(measure))
#    term = map(GT.term,q)
#    prototype = map(GT.prototype,q) |> PartitionedArrays.getany
#    domain = measure.domain
#    GT.quantity(term,prototype,domain)
#end

function weights(measure::Measure,domain::ReferenceDomain)
    mesh = GT.mesh(domain)
    d = GT.num_dims(domain)
    face_to_refid = GT.face_reference_id(mesh,d)
    domface_to_face = GT.faces(domain)
    refid_to_quad = reference_quadratures(measure)
    refid_to_ws = map(GT.weights,refid_to_quad)
    prototype = first(GT.weights(first(refid_to_quad)))
    point_quantity(refid_to_ws,domain;reference=true)
    #GT.quantity() do index
    #    domface = face_index(index,d)
    #    point = point_index(index)
    #    domface_to_face_sym = get_symbol!(index,domface_to_face,gensym("domface_to_face"))
    #    face_to_refid_sym = get_symbol!(index,face_to_refid,gensym("face_to_refid"))
    #    refid_to_ws_sym = get_symbol!(index,refid_to_ws,gensym("refid_to_ws"))
    #    expr = @term begin
    #        face = $domface_to_face_sym[$domface]
    #        ws = $refid_to_ws_sym[$face_to_refid_sym[face]]
    #        ws[$point]
    #    end
    #    expr_term(d,expr,prototype,index)
    #end
end

function weights(measure::Measure,domain::PhysicalDomain)
    Ω = domain
    Ωref = GT.reference_domain(Ω)
    ϕ = GT.physical_map(mesh(Ω),num_dims(Ω))
    w = GT.weights(measure,Ωref)
    x = GT.coordinates(measure,Ωref)
    J = ForwardDiff.jacobian(ϕ,x)
    # J = GT.call(ForwardDiff.jacobian,ϕ,x)
    dV = GT.call(change_of_measure,J)
    GT.call(*,w,dV)
end

function num_points(measure::Measure)
    domain = GT.domain(measure)
    mesh = GT.mesh(domain)
    d = GT.num_dims(domain)
    face_to_refid = GT.face_reference_id(mesh,d)
    domface_to_face = GT.faces(domain)
    refid_to_quad = reference_quadratures(measure)
    refid_to_ws = map(GT.weights,refid_to_quad)
    index -> begin
        face = face_index(index, d)
        #domface_to_face_sym = get_symbol!(index,domface_to_face,gensym("domface_to_face"))
        face_to_refid_sym = get_symbol!(index,face_to_refid,gensym("face_to_refid"))
        refid_to_ws_sym = get_symbol!(index,refid_to_ws,gensym("refid_to_ws"))
            #face = $domface_to_face_sym[$domface]
        expr = @term begin
            ws = $refid_to_ws_sym[$face_to_refid_sym[$face]]
            length(ws)
        end
        expr_term(d,expr,0,index)
    end
end

function integrate(f,measure::Measure)
    domain = GT.domain(measure)
    x = GT.coordinates(measure)
    fx = f(x)
    contrib = integrate_impl(fx,measure)
    integral((measure=>contrib,))
end
const ∫ = integrate

function integrate_impl(fx,measure)
    w = GT.weights(measure)
    GT.call(*,fx,w)
    #p_fx = GT.prototype(fx)
    #p_w = GT.prototype(w)
    #t_fx = GT.term(fx)
    #t_w = GT.term(w)
    #prototype = p_fx*p_w + p_fx*p_w
    #GT.quantity(prototype,GT.domain(measure)) do index
    #    GT.call(*,t_fx(index),t_w(index))
    #end
end

function integrate(f,measure::Measure{<:PMesh})
    x = GT.coordinates(measure)
    fx = f(x)
    q = map(GT.integrate_impl,partition(fx),partition(measure))
    term = map(GT.term,q)
    contrib = GT.quantity(term)
    integral((measure=>contrib,))
end

function integral(contribs)
    Integral(contribs)
end

# TODO Rename contributions with measure_contribution ?
struct Integral{A}  <: GT.AbstractType
    contributions::A
end
contributions(a::Integral) = a.contributions

function Base.sum(int::Integral)
    sum(map(sum_contribution,contributions(int)))
end

function sum_contribution(contrib)
    measure, qty = contrib
    sum_contribution(measure,qty)
end

function sum_contribution(measure,qty)
    domain = GT.domain(measure)
    nfaces = GT.num_faces(domain)
    facemask = fill(true,nfaces)
    sum_contribution_impl(qty,measure,facemask)
end

#function sum_contribution(measure::Measure{<:PMesh},qty::AbstractQuantity{<:PMesh})
#    domain = GT.domain(measure)
#    mesh = domain |> GT.mesh
#    d = GT.num_dims(domain)
#    # TODO allow the user to skip or not to skip ghosts
#    map(partition(measure),partition(qty),GT.face_partition(mesh,d)) do mymeasure,myqty,myfaces
#        mydom = GT.domain(mymeasure)
#        facemask = (part_id(myfaces) .== local_to_owner(myfaces))[GT.faces(mydom)]
#        sum_contribution_impl(myqty,mymeasure,facemask)
#    end |> sum
#end

function sum_contribution_impl(qty,measure,facemask)
    # TODO some code duplication with face_contribution_impl
    # Yes, but this allocates less memory
    term_npoints = GT.num_points(measure) # TODO: term?
    dom = domain(measure)
    d = GT.num_dims(dom)
    index = GT.generate_index(dom)
    sface_to_face = get_symbol!(index,faces(dom),"sface_to_face")
    t = term(qty, index)
    face = face_index(index,d)
    point = point_index(index)
    expr_qty = t |> expression |> simplify
    expr_npoints = term_npoints(index) |> expression |> simplify
    # TODO merge statements
    s_qty = GT.topological_sort(expr_qty,(face,point))
    s_npoints = GT.topological_sort(expr_npoints,(face,))
    expr = quote
        (z,facemask,storage) -> begin
            s = zero(z)
            $(unpack_index_storage(index,:storage))
            $(s_qty[1])
            $(s_npoints[1])
            nfaces = length(facemask)
            for sface in 1:nfaces
                $face = $sface_to_face[sface]
                if ! facemask[sface]
                    continue
                end
                $(s_qty[2])
                npoints = $(s_npoints[2])
                for $point in 1:npoints
                    s += $(s_qty[3])
                end
            end
            s
        end
    end
    loop = eval(expr)
    storage = GT.index_storage(index)
    z = zero(GT.prototype(t))
    s = invokelatest(loop,z,facemask,storage)
    s
end

function face_contribution(int::Integral,domain)
    for (measure,qty) in int.contributions
        domain2 = GT.domain(measure)
        if domain2 == domain
            return face_contribution(measure,qty)
        end
    end
    error("Domain not in integral")
end

function face_contribution(measure,qty)
    domain = GT.domain(measure)
    nfaces = GT.num_faces(domain)
    facemask = fill(true,nfaces)
    face_contribution_impl(qty,measure,facemask)
end

function face_contribution_impl(qty,measure,facemask)
    term_npoints = GT.num_points(measure)
    dom = domain(measure)
    d = GT.num_dims(dom)
    index = GT.generate_index(dom)
    t = term(qty, index)
    face = face_index(index,d)
    point = point_index(index)

    expr_qty = t |> expression |> simplify
    expr_npoints = term_npoints(index) |> expression |> simplify
    # TODO merge statements
    s_qty = GT.topological_sort(expr_qty,(face,point))
    s_npoints = GT.topological_sort(expr_npoints,(face,))
    expr = quote
        (facevals,facemask,storage) -> begin
            $(unpack_index_storage(index,:storage))
            $(s_qty[1])
            $(s_npoints[1])
            nfaces = length(facemask)
            z = zero(eltype(facevals))
            for $face in 1:nfaces
                if ! facemask[$face]
                    continue
                end
                $(s_qty[2])
                npoints = $(s_npoints[2])
                s = z
                for $point in 1:npoints
                    s += $(s_qty[3])
                end
                facevals[$face] = s
            end
            facevals
        end
    end
    loop! = eval(expr)
    storage = GT.index_storage(index)
    z = zero(GT.prototype(t))
    nfaces = length(facemask)
    facevals = fill(z,nfaces)
    invokelatest(loop!,facevals,facemask,storage)
    facevals
end

function face_diameter(Ω)
    dΩ = GT.measure(Ω,1)
    d = GT.num_dims(Ω)
    u = GT.analytical_field(x->1,Ω)
    int = ∫(u,dΩ)
    r = face_contribution(int,Ω)
    r .= r .^ (1/d)
    r
end

function face_diameter_field(Ω)
    dims = GT.face_diameter(Ω)
    face_constant_field(dims,Ω)
end

function Base.:+(int1::Integral,int2::Integral)
    # TODO merge contributions on the same domain?
    contribs = (GT.contributions(int1)...,GT.contributions(int2)...)
    GT.integral(contribs)
end

function Base.:-(int1::Integral,int2::Integral)
    int1 + (-1)*int2
end

function Base.:-(int1::Real,int2::Integral)
    @assert int1 == 0
    (-1)*int2
end

function Base.:*(v::Number,int::Integral)
    contribs = map(GT.contributions(int)) do domain_and_contribution
        domain, contribution = domain_and_contribution
        domain => v*contribution
    end
    integral(contribs)
end

function Base.:*(int::Integral,v::Number)
    v*int
end

function Base.:/(int::Integral,v::Number)
    (1/v)*int
end

function quadrature(m::Measure)
    mesh_quadrature(;
        domain=domain(m),
        reference_quadratures = reference_quadratures(m),
        face_reference_id = face_reference_id(m)
       )
end

function num_points_accessor(measure::Measure)
    num_points_accessor(quadrature(measure))
end

function coordinate_accessor(measure::Measure)
    coordinate_accessor(quadrature(measure))
end
function jacobian_accessor(measure::Measure)
    jacobian_accessor(quadrature(measure))
end
function weight_accessor(measure::Measure)
    weight_accessor(quadrature(measure))
end

"""
    abstract type AbstractQuadrature

# Basic queries

- [`coordinates`](@ref)
- [`weights`](@ref)

# Basic constructors

- [`default_quadrature`](@ref)
- [`duffy_quadrature`](@ref)
- [`tensor_product_quadrature`](@ref)

# Supertype hierarchy

    AbstractQuadrature <: GT.AbstractType
"""
abstract type AbstractQuadrature <: GT.AbstractType end

struct GenericCuadrature{A,B} <: AbstractQuadrature
    coordinates::A
    weights::B
end
struct Cuadrature{D,T} <: GT.AbstractType
    coordinates::Vector{SVector{D,T}}
    weights::Vector{T}
end
function quadrature(coordinates,weights)
    GenericCuadrature(coordinates,weights)
end
function quadrature(coordinates::Vector{SVector{D,T}},weights::Vector{T}) where {D,T}
    Cuadrature(coordinates,weights)
end

"""
"""
function duffy_quadrature(geo,degree)
    @assert is_unit_simplex(geo)
    D = num_dims(geo)
    if D == 0
        x = zeros(SVector{0,real_type},1)
        w = ones(real_type,1)
        return quadrature(x,w)
    end
    function map_to(a,b,(points,weights))
      points_ab = similar(points)
      weights_ab = similar(weights)
      points_ab .= 0.5*(b-a)*points .+ 0.5*(a+b)
      weights_ab .= 0.5*(b-a)*weights
      (points_ab, weights_ab)
    end
    function duffy_map(q)
        D = length(q)
        a = 1.0
        m = ntuple(Val(D)) do i
            if i == 1
                q[i]
            else
                a *= (1-q[i-1])
                a*q[i]
            end
        end
        typeof(q)(m)
    end
    n = ceil(Int, (degree + 1.0) / 2.0 )
    beta = 0
    dim_to_quad_1d = map(1:(D-1)) do d
        alpha = (D-1)-(d-1)
        map_to(0,1,gaussjacobi(n,alpha,beta))
    end
    quad_1d = map_to(0,1,gausslegendre(n))
    push!(dim_to_quad_1d,quad_1d)
    coords_per_dir = map(first,dim_to_quad_1d)
    weights_per_dir =  map(last,dim_to_quad_1d)
    a = 0.5
    for d in (D-1):-1:1
        ws_1d = weights_per_dir[d]
        ws_1d[:] *= a
        a *= 0.5
    end
    m = prod(map(length,weights_per_dir))
    Tv = real_type(geo)
    w = zeros(Tv,m)
    x = zeros(SVector{D,Tv},m)
    tensor_product!(identity,x,coords_per_dir)
    tensor_product!(prod,w,weights_per_dir)
    x .= duffy_map.(x)
    quadrature(x,w)
end

"""
"""
function tensor_product_quadrature(geo,degree_per_dir)
    @assert is_n_cube(geo)
    @assert is_axis_aligned(geo)
    my_bounding_box = bounding_box(geo)
    D = num_dims(geo)
    limits_per_dir = ntuple(i->(my_bounding_box[1][i],my_bounding_box[2][i]),Val(D))
    n_per_dir = map(d->ceil(Int,(d+1)/2),degree_per_dir)
    function quadrature_1d(n,limits)
        x,w = FastGaussQuadrature.gausslegendre(n)
        a,b = limits
        x .= (0.5*(b-a)) .*x .+ (0.5*(b+a))
        w .*= 0.5*(b-a)
        quadrature(x,w)
    end
    quad_per_dir = map(quadrature_1d,n_per_dir,limits_per_dir)
    coords_per_dir = map(coordinates,quad_per_dir)
    weights_per_dir = map(weights,quad_per_dir)
    m = prod(map(length,weights_per_dir))
    Tv = real_type(geo)
    w = zeros(Tv,m)
    x = zeros(SVector{D,Tv},m)
    tensor_product!(identity,x,coords_per_dir)
    tensor_product!(prod,w,weights_per_dir)
    quadrature(x,w)
end

function tensor_product!(f,result,values_per_dir)
    shape = Tuple(map(length,values_per_dir))
    cis = CartesianIndices(shape)
    lis = LinearIndices(cis)
    for ci in cis
        li = lis[ci]
        result[li] = f(map((q,i)->q[i],values_per_dir,Tuple(ci)))
    end
    result
end

"""
"""
function default_quadrature(geo,degree)
    if is_n_cube(geo) && is_axis_aligned(geo)
        D = num_dims(geo)
        degree_per_dir = repeat_per_dir(geo,degree)
        tensor_product_quadrature(geo,degree_per_dir)
    elseif is_unit_simplex(geo)
        duffy_quadrature(geo,degree)
    else
        error("Not implemented")
    end
end

function measure(dom::AbstractDomain,degree)
    measure(default_quadrature,dom,degree)
end

function measure(f,dom::AbstractDomain,degree)
    mesh = GT.mesh(dom)
    Measure(mesh,f,dom,degree)
end

struct Measure{A,B,C,D}
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
    d = GT.face_dim(domain)
    drefid_refdface = GT.reference_faces(mesh,d)
    refid_to_quad = map(drefid_refdface) do refdface
        geo = GT.geometry(refdface)
        measure.quadrature_rule(geo,measure.degree)
    end
    refid_to_quad
end

function coordinates(measure::Measure)
    domain = GT.domain(measure)
    coordinates(measure,domain)
end

function coordinates(measure::Measure{<:PMesh})
    q = map(GT.coordinates,partition(measure))
    term = map(GT.term,q)
    prototype = map(GT.prototype,q) |> PartitionedArrays.getany
    domain = measure.domain
    GT.quantity(term,prototype,domain)
end

function coordinates(measure::Measure,domain::ReferenceDomain)
    mesh = GT.mesh(domain)
    d = GT.face_dim(domain)
    face_to_refid = GT.face_reference_id(mesh,d)
    domface_to_face = GT.faces(domain)
    refid_to_quad = reference_quadratures(measure)
    refid_to_coords = map(GT.coordinates,refid_to_quad)
    prototype = first(GT.coordinates(first(refid_to_quad)))
    GT.quantity(prototype,domain) do index
        domface = index.face
        point = index.point
        domface_to_face_sym = get!(index.dict,domface_to_face,gensym("domface_to_face"))
        face_to_refid_sym = get!(index.dict,face_to_refid,gensym("face_to_refid"))
        refid_to_coords_sym = get!(index.dict,refid_to_coords,gensym("refid_to_coords"))
        @term begin
            face = $domface_to_face_sym[$domface]
            coords = reference_value($refid_to_coords_sym,$face_to_refid_sym,face)
            coords[$point]
        end
    end
end

function coordinates(measure::Measure,domain::PhysicalDomain)
    Ω = domain
    Ωref = GT.reference_domain(Ω)
    ϕ = GT.domain_map(Ωref,Ω)
    x = GT.coordinates(measure,Ωref)
    ϕ(x)
end

function weights(measure::Measure)
    domain = GT.domain(measure)
    weights(measure,domain)
end

function weights(measure::Measure{<:PMesh})
    q = map(GT.weights,partition(measure))
    term = map(GT.term,q)
    prototype = map(GT.prototype,q) |> PartitionedArrays.getany
    domain = measure.domain
    GT.quantity(term,prototype,domain)
end

function weights(measure::Measure,domain::ReferenceDomain)
    mesh = GT.mesh(domain)
    d = GT.face_dim(domain)
    face_to_refid = GT.face_reference_id(mesh,d)
    domface_to_face = GT.faces(domain)
    refid_to_quad = reference_quadratures(measure)
    refid_to_ws = map(GT.weights,refid_to_quad)
    prototype = first(GT.weights(first(refid_to_quad)))
    GT.quantity(prototype,domain) do index
        domface = index.face
        point = index.point
        domface_to_face_sym = get!(index.dict,domface_to_face,gensym("domface_to_face"))
        face_to_refid_sym = get!(index.dict,face_to_refid,gensym("face_to_refid"))
        refid_to_ws_sym = get!(index.dict,refid_to_ws,gensym("refid_to_ws"))
        @term begin
            face = $domface_to_face_sym[$domface]
            ws = reference_value($refid_to_ws_sym,$face_to_refid_sym,face)
            ws[$point]
        end
    end
end

function change_of_measure(J)
    sqrt(det(transpose(J)*J))
end

function weights(measure::Measure,domain::PhysicalDomain)
    Ω = domain
    Ωref = GT.reference_domain(Ω)
    ϕ = GT.domain_map(Ωref,Ω)
    w = GT.weights(measure,Ωref)
    x = GT.coordinates(measure,Ωref)
    J = GT.call(ForwardDiff.jacobian,ϕ,x)
    dV = GT.call(change_of_measure,J)
    GT.call(*,w,dV)
end

function num_points(measure::Measure)
    domain = GT.domain(measure)
    mesh = GT.mesh(domain)
    d = GT.face_dim(domain)
    face_to_refid = GT.face_reference_id(mesh,d)
    domface_to_face = GT.faces(domain)
    refid_to_quad = reference_quadratures(measure)
    refid_to_ws = map(GT.weights,refid_to_quad)
    index -> begin
        domface = index.face
        domface_to_face_sym = get!(index.dict,domface_to_face,gensym("domface_to_face"))
        face_to_refid_sym = get!(index.dict,face_to_refid,gensym("face_to_refid"))
        refid_to_ws_sym = get!(index.dict,refid_to_ws,gensym("refid_to_ws"))
        @term begin
            face = $domface_to_face_sym[$domface]
            ws = reference_value($refid_to_ws_sym,$face_to_refid_sym,face)
            length(ws)
        end
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
    prototype = map(GT.prototype,q) |> PartitionedArrays.getany
    domain = measure |> GT.domain
    contrib = GT.quantity(term,prototype,domain)
    integral((measure=>contrib,))
end

function integral(contribs)
    Integral(contribs)
end

# TODO Rename contributions with measure_contribution ?
struct Integral{A}
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

function sum_contribution(measure::Measure{<:PMesh},qty::AbstractQuantity{<:PMesh})
    domain = GT.domain(measure)
    mesh = domain |> GT.mesh
    d = GT.num_dims(domain)
    # TODO allow the user to skip or not to skip ghosts
    map(partition(measure),partition(qty),GT.face_partition(mesh,d)) do mymeasure,myqty,myfaces
        mydom = GT.domain(mymeasure)
        facemask = (part_id(myfaces) .== local_to_owner(myfaces))[GT.faces(mydom)]
        sum_contribution_impl(myqty,mymeasure,facemask)
    end |> sum
end

function sum_contribution_impl(qty,measure,facemask)
    # TODO some code duplication with face_contribution_impl
    # Yes, but this allocates less memory
    term_qty = GT.term(qty)
    term_npoints = GT.num_points(measure)
    face = :face
    point = :point
    index = GT.index(;face,point)
    expr_qty = term_qty(index) |> simplify
    expr_npoints = term_npoints(index) |> simplify
    # TODO merge statements
    s_qty = GT.topological_sort(expr_qty,(face,point))
    s_npoints = GT.topological_sort(expr_npoints,(face,))
    expr = quote
        function loop(z,facemask,storage)
            s = zero(z)
            $(unpack_storage(index.dict,:storage))
            $(s_qty[1])
            $(s_npoints[1])
            nfaces = length(facemask)
            for $face in 1:nfaces
                if ! facemask[$face]
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
    storage = GT.storage(index)
    z = zero(GT.prototype(qty))
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
    term_qty = GT.term(qty)
    term_npoints = GT.num_points(measure)
    face = :face
    point = :point
    index = GT.index(;face,point)
    expr_qty = term_qty(index) |> simplify
    expr_npoints = term_npoints(index) |> simplify
    # TODO merge statements
    s_qty = GT.topological_sort(expr_qty,(face,point))
    s_npoints = GT.topological_sort(expr_npoints,(face,))
    expr = quote
        function loop!(facevals,facemask,storage)
            $(unpack_storage(index.dict,:storage))
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
    storage = GT.storage(index)
    z = zero(GT.prototype(qty))
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


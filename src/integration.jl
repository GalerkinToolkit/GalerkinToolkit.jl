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

    AbstractQuadrature <: gk.AbstractType
"""
abstract type AbstractQuadrature <: gk.AbstractType end

struct GenericCuadrature{A,B} <: AbstractQuadrature
    coordinates::A
    weights::B
end
struct Cuadrature{D,T} <: gk.AbstractType
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
    Measure(f,dom,degree)
end

struct Measure{A,B,C}
    quadrature_rule::A
    domain::B
    degree::C
end

domain(a::Measure) = a.domain

function reference_quadratures(measure::Measure)
    domain = gk.domain(measure)
    @assert isa(gk.domain_style(domain),GlobalDomain)
    mesh = gk.mesh(domain)
    d = gk.face_dim(domain)
    drefid_refdface = gk.reference_faces(mesh,d)
    refid_to_quad = map(drefid_refdface) do refdface
        geo = gk.geometry(refdface)
        measure.quadrature_rule(geo,measure.degree)
    end
    refid_to_quad
end

function coordinates(measure::Measure)
    domain = gk.domain(measure)
    coordinates(measure,gk.domain_style(domain))
end

function coordinates(measure::Measure,stype::GlobalDomain{true})
    domain = gk.domain(measure)
    mesh = gk.mesh(domain)
    d = gk.face_dim(domain)
    face_to_refid = gk.face_reference_id(mesh,d)
    domface_to_face = gk.faces(domain)
    refid_to_quad = reference_quadratures(measure)
    refid_to_coords = map(gk.coordinates,refid_to_quad)
    prototype = first(gk.coordinates(first(refid_to_quad)))
    gk.quantity(prototype,domain) do index
        domface = index.face
        face = domface_to_face[domface]
        refid = face_to_refid[face]
        coords = refid_to_coords[refid]
        point = index.point
        coords[point]
    end
end

function weights(measure::Measure)
    domain = gk.domain(measure)
    weights(measure,gk.domain_style(domain))
end

function weights(measure::Measure,stype::GlobalDomain{true})
    domain = gk.domain(measure)
    mesh = gk.mesh(domain)
    d = gk.face_dim(domain)
    face_to_refid = gk.face_reference_id(mesh,d)
    domface_to_face = gk.faces(domain)
    refid_to_quad = reference_quadratures(measure)
    refid_to_w = map(gk.weights,refid_to_quad)
    prototype = first(gk.weights(first(refid_to_quad)))
    gk.quantity(prototype,domain) do index
        domface = index.face
        face = domface_to_face[domface]
        refid = face_to_refid[face]
        w = refid_to_w[refid]
        point = index.point
        w[point]
    end
end

function num_points(measure::Measure)
    domain = gk.domain(measure)
    @assert isa(gk.domain_style(domain),gk.GlobalDomain)
    mesh = gk.mesh(domain)
    d = gk.face_dim(domain)
    face_to_refid = gk.face_reference_id(mesh,d)
    domface_to_face = gk.faces(domain)
    refid_to_quad = reference_quadratures(measure)
    refid_to_w = map(gk.weights,refid_to_quad)
    index -> begin
        domface = index.face
        face = domface_to_face[domface]
        refid = face_to_refid[face]
        w = refid_to_w[refid]
        length(w)
    end
end

function integrate(f,measure::Measure)
    domain = gk.domain(measure)
    x = gk.coordinates(measure)
    w = gk.weights(measure)
    fx = f(x)
    p_fx = gk.prototype(fx)
    p_w = gk.prototype(w)
    t_fx = gk.term(fx)
    t_w = gk.term(w)
    prototype = p_fx*p_w + p_fx*p_w
    num_points = gk.num_points(measure)
    contrib = gk.quantity(prototype,domain) do index
        np = num_points(index)
        sum(1:np) do point
            index2 = replace_point(index,point)
            t_fx(index2)*t_w(index2)
        end
    end
    integral((domain=>contrib,))
end
const ∫ = integrate

integral(contribs) = Integral(contribs)
struct Integral{A}
    contributions::A
end
contributions(a::Integral) = a.contributions

function sum_contribution(contrib)
    domain, qty = contrib
    z = zero(gk.prototype(qty))
    nfaces = gk.num_faces(domain)
    term = gk.term(qty)
    for face in 1:nfaces
        index = gk.index(;face)
        z += term(index)
    end
    z
end

function Base.sum(int::Integral)
    sum(map(sum_contribution,contributions(int)))
end


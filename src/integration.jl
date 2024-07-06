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
        face = domface_to_face[domface]
        refid = face_to_refid[face]
        coords = refid_to_coords[refid]
        point = index.point
        coords[point]
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
    refid_to_w = map(GT.weights,refid_to_quad)
    prototype = first(GT.weights(first(refid_to_quad)))
    GT.quantity(prototype,domain) do index
        domface = index.face
        face = domface_to_face[domface]
        refid = face_to_refid[face]
        w = refid_to_w[refid]
        point = index.point
        w[point]
    end
end

function weights(measure::Measure,domain::PhysicalDomain)
    Ω = domain
    Ωref = GT.reference_domain(Ω)
    ϕ = GT.domain_map(Ωref,Ω)
    w = GT.weights(measure,Ωref)
    x = GT.coordinates(measure,Ωref)
    f(J) = sqrt(det(transpose(J)*J))
    J = GT.call(ForwardDiff.jacobian,ϕ,x)
    dV = GT.call(f,J)
    GT.call(*,w,dV)
end

function num_points(measure::Measure)
    domain = GT.domain(measure)
    mesh = GT.mesh(domain)
    d = GT.face_dim(domain)
    face_to_refid = GT.face_reference_id(mesh,d)
    domface_to_face = GT.faces(domain)
    refid_to_quad = reference_quadratures(measure)
    refid_to_w = map(GT.weights,refid_to_quad)
    index -> begin
        domface = index.face
        face = domface_to_face[domface]
        refid = face_to_refid[face]
        w = refid_to_w[refid]
        length(w)
    end
end

function integrate(f,measure::Measure)
    domain = GT.domain(measure)
    x = GT.coordinates(measure)
    fx = f(x)
    contrib = integrate_impl(fx,measure)
    integral((domain=>contrib,))
end
const ∫ = integrate

function integrate_impl(fx,measure)
    w = GT.weights(measure)
    p_fx = GT.prototype(fx)
    p_w = GT.prototype(w)
    t_fx = GT.term(fx)
    t_w = GT.term(w)
    prototype = p_fx*p_w + p_fx*p_w
    num_points = GT.num_points(measure)
    contrib = GT.quantity(prototype,GT.domain(measure)) do index
        np = num_points(index)
        sum(1:np) do point
            index2 = replace_point(index,point)
            t_fx(index2)*t_w(index2)
        end
    end
end

function integrate(f,measure::Measure{<:PMesh})
    x = GT.coordinates(measure)
    fx = f(x)
    q = map(GT.integrate_impl,partition(fx),partition(measure))
    term = map(GT.term,q)
    prototype = map(GT.prototype,q) |> PartitionedArrays.getany
    domain = measure |> GT.domain
    contrib = GT.quantity(term,prototype,domain)
    integral((domain=>contrib,))
end

function integral(contribs)
    Integral(contribs)
end

# TODO Rename contributions with domain_contribution
struct Integral{A}
    contributions::A
end
contributions(a::Integral) = a.contributions

function Base.sum(int::Integral)
    sum(map(sum_contribution,contributions(int)))
end

function sum_contribution(contrib)
    domain, qty = contrib
    sum_contribution(domain,qty)
end

function sum_contribution(domain,qty)
    nfaces = GT.num_faces(domain)
    facemask = fill(true,nfaces)
    sum_contribution_impl(qty,facemask)
end

function sum_contribution(domain::AbstractDomain{<:PMesh},qty::AbstractQuantity{<:PMesh})
    mesh = domain |> GT.mesh
    d = GT.num_dims(domain)
    # TODO allow the user to skip or not to skip ghosts
    map(partition(domain),partition(qty),GT.face_partition(mesh,d)) do mydom,myqty,myfaces
        facemask = (part_id(myfaces) .== local_to_owner(myfaces))[GT.faces(mydom)]
        sum_contribution_impl(myqty,facemask)
    end |> sum
end

function sum_contribution_impl(qty,facemask)
    # TODO some code duplication with face_contribution_impl
    z = zero(GT.prototype(qty))
    nfaces = length(facemask)
    term = GT.term(qty)
    for face in 1:nfaces
        if ! facemask[face]
            continue
        end
        index = GT.index(;face)
        z += term(index)
    end
    z
end

function face_contribution(int::Integral,domain)
    for (domain2,qty) in int.contributions
        if domain2 == domain
            return face_contribution(domain2,qty)
        end
    end
    error("Domain not in integral")
end

function face_contribution(domain,qty)
    nfaces = GT.num_faces(domain)
    facemask = fill(true,nfaces)
    face_contribution_impl(qty,facemask)
end

function face_contribution_impl(qty,facemask)
    z = zero(GT.prototype(qty))
    nfaces = length(facemask)
    r = fill(z,nfaces)
    term = GT.term(qty)
    for face in 1:nfaces
        if ! facemask[face]
            continue
        end
        index = GT.index(;face)
        r[face] = term(index)
    end
    r
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


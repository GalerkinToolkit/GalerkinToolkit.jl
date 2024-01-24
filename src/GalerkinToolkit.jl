module GalerkinToolkit

using WriteVTK
using FastGaussQuadrature
using StaticArrays
using LinearAlgebra
using ForwardDiff
using Gmsh
using PartitionedArrays
using Combinatorics
using SparseArrays

abstract type GalerkinToolkitDataType end
function Base.show(io::IO,data::GalerkinToolkitDataType)
    print(io,"GalerkinToolkit.$(nameof(typeof(data)))(…)")
end

val_parameter(a) = a
val_parameter(::Val{a}) where a = a
num_dims(a) = val_parameter(a.num_dims)
coordinates(a) = a.coordinates
weights(a) = a.weights

struct UnitSimplex{D} <: GalerkinToolkitDataType
    num_dims::Val{D}
end
function Base.show(io::IO,k::MIME"text/plain",data::UnitSimplex)
    D = num_dims(data)
    print(io,"Unit simplex of dimension $D")
end

struct UnitNCube{D} <: GalerkinToolkitDataType
    num_dims::Val{D}
end
function Base.show(io::IO,k::MIME"text/plain",data::UnitNCube)
    D = num_dims(data)
    print(io,"Unit cube of dimension $D")
end

function unit_simplex(num_dims)
    D = val_parameter(num_dims)
    UnitSimplex(Val(D))
end

function unit_n_cube(num_dims)
    D = val_parameter(num_dims)
    UnitNCube(Val(D))
end

struct GenericCuadrature{A,B} <: GalerkinToolkitDataType
    coordinates::A
    weights::B
end
struct Cuadrature{D,T} <: GalerkinToolkitDataType
    coordinates::Vector{SVector{D,T}}
    weights::Vector{T}
end
function quadrature(coordinates,weights)
    GenericCuadrature(coordinates,weights)
end
function quadrature(coordinates::Vector{SVector{D,T}},weights::Vector{T}) where {D,T}
    Cuadrature(coordinates,weights)
end

function default_quadrature(geo::UnitSimplex;
    degree,
    real_type=Float64)
    duffy_quadrature(geo;degree,real_type)
end

function default_quadrature(geo::UnitNCube;
    degree,
    real_type=Float64)
    degree_per_dir = ntuple(i->degree,Val(num_dims(geo)))
    tensor_product_quadrature(geo;degree_per_dir,real_type)
end

function duffy_quadrature(geo::UnitSimplex;
    degree,
    real_type = Float64,
    )
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
    w = zeros(real_type,m)
    x = zeros(SVector{D,real_type},m)
    tensor_product!(identity,x,coords_per_dir)
    tensor_product!(prod,w,weights_per_dir)
    x .= duffy_map.(x)
    quadrature(x,w)
end

function tensor_product_quadrature(geo::UnitNCube;
    degree_per_dir,
    real_type = Float64,
    )

    D = num_dims(geo)
    limits_per_dir = ntuple(i->(0,1),Val(D))
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
    w = zeros(real_type,m)
    x = zeros(SVector{D,real_type},m)
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

end # module

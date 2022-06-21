
call(f,x...) = f(x...)
evaluate(g,f,x) = broadcast(g,f,fill(x))
evaluate(f,x) = evaluate(call,f,x)

function sample(g,f,x)
  axesf = axes(Base.broadcastable(f))
  x′= reshape_sampling_points(x,axesf)
  broadcast(g,f,x′)
end
sample(f,x) = sample(call,f,x)
reshape_sampling_points(x::AbstractVector,::Tuple{}) = x
reshape_sampling_points(x::AbstractVector,::Tuple{Any}) = permutedims(x)

function inverse_map end

struct Monomial{N,T} <: Function
  exps::NTuple{N,T}
  Monomial(exps::NTuple{N,T}) where {N,T<:Integer} = new{N,T}(exps)
  Monomial(exps::Vararg{T,N}) where {N,T<:Integer} = new{N,T}(exps)
end
(f::Monomial)(x) = prod(x.^f.exps)
Base.show(io::IO,m::Monomial) = print(io,"Monomial",m.exps)
Base.broadcastable(f::Monomial) = fill(f)

struct NodalValue{T} <: Function
  node::T
end
(x::NodalValue)(f) = f(x.node)
Base.show(io::IO,m::NodalValue) = print(io,"NodalValue",m.node)
Base.broadcastable(f::NodalValue) = fill(f)

scale(p,coeff,f) = ScaledFunction(p,coeff,f)
struct ScaledFunction{P,A,F} <: Function
  p::P
  coeff::A
  f::F
end
(g::ScaledFunction)(x) = g.p(g.coeff,g.f(x))

linear_combination(p,coeff::AbstractVector,b::AbstractVector) = FunctionCombination(p,coeff,b)
linear_combination(p,coeff::AbstractMatrix,b::AbstractVector) = [ linear_combination(p,view(coeff,i,:),b) for i in 1:size(coeff,1) ]
linear_combination(a,b) = linear_combination(*,a,b)
struct FunctionCombination{P,C<:AbstractVector,F<:AbstractVector} <: Function
  p::P
  coeff::C
  f::F
end
(g::FunctionCombination)(x) = linear_combination_value(g.p,g.coeff,value(g.f,x))

function typeof_linear_combination(f,a,b)
  za = zero(eltype(a))
  zb = zero(eltype(b))
  typeof(f(za,zb)+f(za,zb))
end
function linear_combination_value(f,a::AbstractVector,b::AbstractVector)
  z = zero(typeof_linear_combination(f,a,b))
  @assert length(a) == length(b)
  for i in 1:length(a)
    z += f(a[i],b[i])
  end
  z
end

struct Operator{T}
  f::T
end
(o::Operator)(args...) = FunctionOperation(o.f,args)
struct FunctionOperation{T,A} <: Function
  f::T
  args::A
end
(g::FunctionOperation)(x) = g.f(map(i->i(x),g.args)...)

struct AffineMap{A,B} <: Function
  jacobian::A
  offset::B
end
(f::AffineMap)(x) = f.jacobian*x+f.offset
function inverse_map(f::AffineMap)
  A = f.jacobian
  b = f.offset
  Ainv = inv(A)
  AffineMap(Ainv,-(Ainv*b))
end

Q_basis(k::Integer,::Val{d}) where d = Q_basis(ntuple(i->k,Val(d)))
Q_basis(orders::Integer...) = Q_basis(orders)
function Q_basis(orders::Tuple)
  cis = CartesianIndices(orders.+1)
  exps = [ Tuple(cis[i]).-1  for i in 1:length(cis) ]
  map(Monomial,exps)
end

function P_basis(k::Integer,::Val{d}) where d
  m = Q_basis(k,Val(d))
  filter(mi->sum(mi.exps)<=k,m)
end

function P̃_basis(k::Integer,::Val{d}) where d
  m = P_basis(k,Val(d))
  filter(i->sum(i.exps)>(k-1),m)
end

function S̃_basis(k::Integer,::Val{2})
  e1 = SVector(1,0)
  e2 = SVector(0,1)
  function f(α)
    m1 = Monomial(α-1,k-α+1)
    m2 = Monomial(α,k-α)
    linear_combination([-e1,e2],[m1,m2])
  end
  [ f(α) for α in 1:k ]
end

function S̃_basis(k::Integer,::Val{3})
  e1 = SVector(1,0,0)
  e2 = SVector(0,1,0)
  e3 = SVector(0,0,1)
  function a(α,β)
    m1 = Monomial(α-1,k-α-β-2,β-1)
    m2 = Monomial(α,k-α-β-1,β-1)
    linear_combination([-e1,e2],[m1,m2])
  end
  function b(α,β)
    m1 = Monomial(k-α-β+1,β-1,α)
    m3 = Monomial(k-α-β+2,β-1,α-1)
    linear_combination([-e1,e3],[m1,m3])
  end
  function c(α)
    m2 = Monomial(0,α-1,k-α+1)
    m3 = Monomial(0,α,k-α)
    linear_combination([-e2,e3],[m2,m3])
  end
  m = [ c(α) for α in 1:k]
  for β in 1:k
    for α in 1:(k+1-β)
      push!(m,a(α,β))
      push!(m,b(α,β))
    end
  end
  m
end

function cartesian_product(
  args...;
  inner=SVector,
  outer=ms->vcat(ms...))

  d = length(args)
  ms = map(args,ntuple(z->z,Val(d))) do s,i
    coeff = inner(ntuple(k-> (k==i ? 1 : 0),Val(d)))
    map(si->scale(*,coeff,si),basis(s))
  end
  outer(ms)
end

function direct_sum(a...)
  # TODO improve type stability of this one
  vcat(args...)
end




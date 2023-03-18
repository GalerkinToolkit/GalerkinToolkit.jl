
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

function cartesian_product(
  args...;
  inner=MathTuple,
  outer=ms->BlockArrays.mortar(collect(ms)))

  d = length(args)
  ms = map(args,ntuple(z->z,Val(d))) do s,i
    coeff = inner(ntuple(k-> (k==i ? 1 : 0),Val(d)))
    map(si->scale(*,coeff,si),s)
  end
  outer(ms)
end

function direct_sum(a...)
  # TODO improve type stability of this one
  vcat(args...)
end

function vector_valued_basis(basis,::Val{d}) where d
  args = ntuple(i->basis,Val(d))
  cartesian_product(args...,inner=SVector,outer=ms->vcat(ms...))
end





function prefix!(a)
  a[1] = one(eltype(a))
  n = length(a)
  @inbounds for i in 1:(n-1)
    a[i+1] += a[i]
  end
  a
end

function rewind!(a)
  n = length(a)
  @inbounds for i in (n-1):-1:1
    a[i+1] = a[i]
  end
  a[1] = one(eltype(a))
  a
end

struct JaggedArray{T,A,B} <: AbstractVector{T}
  data::A
  ptrs::B
end

const JArray{T} = JaggedArray{SubArray{T,1,Vector{T},Tuple{UnitRange{Int32}},true},Vector{T},Vector{Int32}}

function JaggedArray(data,ptrs)
  Tp = eltype(ptrs)
  T = typeof(view(data,Tp(1):Tp(0)))
  A = typeof(data)
  B = typeof(ptrs)
  JaggedArray{T,A,B}(data,ptrs)
end

function JaggedArray(a::AbstractArray{<:AbstractArray{T}}) where T
  n = length(a)
  ptrs = Vector{Int32}(undef,n+1)
  u = one(eltype(ptrs))
  @inbounds for i in 1:n
    ai = a[i]
    ptrs[i+1] = length(ai)
  end
  prefix!(ptrs)
  ndata = ptrs[end]-u
  data = Vector{T}(undef,ndata)
  p = 1
  @inbounds for i in 1:n
    ai = a[i]
    for j in 1:length(ai)
      aij = ai[j]
      data[p] = aij
      p += 1
    end
  end
  JaggedArray(data,ptrs)
end

JaggedArray(a::JaggedArray) = a

Base.size(a::JaggedArray) = (length(a.ptrs)-1,)

function Base.getindex(a::JaggedArray,i::Integer)
  u = one(eltype(a.ptrs))
  pini = a.ptrs[i]
  pend = a.ptrs[i+1]-u
  view(a.data,pini:pend)
end

function Base.convert(::Type{J},vv) where J <: JaggedArray
  a = JaggedArray(vv)
  J(a.data,a.ptrs)
end

function Base.convert(::Type{J},vv::J) where J <: JaggedArray
  vv
end


struct MathTuple{A<:Tuple}
  tuple::A
end
math_tuple(args...) = MathTuple(args)
Base.:+(a::MathTuple,b::MathTuple) = MathTuple(a.tuple .+ b.tuple)
Base.:-(a::MathTuple,b::MathTuple) = MathTuple(a.tuple .- b.tuple)
Base.:*(a,b::MathTuple) = MathTuple(map(i->a*b,b.tuple))
Base.:*(a::MathTuple,b) = MathTuple(map(i->i*b,a.tuple))
Base.:/(a,b::MathTuple) = MathTuple(map(i->a/i,b.tuple))
Base.:/(a::MathTuple,b) = MathTuple(map(i->i/b,a.tuple))
(f::MathTuple)(x) = MathTuple(i->i(x),f.tuple)
Base.length(a::MathTuple) = length(a.tuple)
Base.length(a::MathTuple{A}) where A = length(A)
@inline Base.getindex(a::MathTuple,i::Integer) = a.tuple[i]

module JaggedArrayTests

using GalerkinToolkit
using Test

a = [[1,2],[3,4,5],Int[],[3,4]]
b = JaggedArray(a)
@test a == b
@test b === JaggedArray(b)

T = typeof(b)
c = T(b.data,b.ptrs)
@test c == b

d = collect(c)


end


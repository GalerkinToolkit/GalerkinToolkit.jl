module SpaceTests

import GalerkinToolkit as GT
using Test
using LinearAlgebra

using InteractiveUtils

#@code_warntype GT.unit_simplex(Val(3))

spx0 = GT.unit_simplex(0)
spx1 = GT.unit_simplex(1)
spx2 = GT.unit_simplex(2)
spx3 = GT.unit_simplex(3)
@test GT.num_dims(spx3) == 3
display(spx3)

cube0 = GT.unit_n_cube(0)
cube1 = GT.unit_n_cube(1)
cube2 = GT.unit_n_cube(2)
cube3 = GT.unit_n_cube(3)
@test GT.num_dims(cube3) == 3
display(cube3)

#@code_warntype GT.options(cube0)
#@code_warntype GT.options()

options = GT.options(cube0)
@show GT.real_type(options)
@show GT.int_type(options)

#@code_warntype GT.real_type(options)

degree = 2
qua = GT.quadrature(cube2,degree)

#@code_warntype GT.quadrature(spx2,2)

@test cube2 === GT.domain(qua)
x = GT.coordinates(qua)
@test sum(GT.weights(qua)) ≈ 1

degree = 2
qua = GT.quadrature(spx2,degree)
@test spx2 === GT.domain(qua)
x = GT.coordinates(qua)
@test sum(GT.weights(qua)) ≈ 0.5

fe = GT.lagrange_space(cube2)
display(fe)

fe2 = GT.lagrange_space(cube2)

@test isequal(fe,fe2)
@test hash(fe) == hash(fe2)

x = GT.monomial_exponents(fe)

#@code_warntype GT.monomial_exponents(fe)

display(x)

x = GT.node_coordinates(fe)
#@code_warntype GT.node_coordinates(fe)
display(x)

@test GT.num_dofs(fe) == GT.num_nodes(fe)

A = GT.tabulator(fe)(GT.value,x)
n = GT.num_dofs(fe)
@test A ≈ Matrix{Float64}(I,n,n)

t = GT.tabulator(fe)
#@code_warntype t(GT.value,x)

qua = GT.node_quadrature(fe)
#@code_warntype GT.node_quadrature(fe)

fe = GT.lagrange_space(cube2;tensor_size=Val((2,)))

@test GT.num_dofs(fe) == 2*GT.num_nodes(fe)

A = GT.tabulator(fe)(GT.value,x)

@show GT.node_dofs(fe)
@show GT.dof_node(fe)

fe = GT.lagrange_space(cube2;order=0)
x = GT.node_coordinates(fe)
@test x[1] ≈ [0.5,0.5]

fe = GT.lagrange_space(spx0)
x = GT.monomial_exponents(fe)
display(x)

x = GT.node_coordinates(fe)
display(x)

@test GT.order(fe) == 0

A = GT.tabulator(fe)(GT.value,x)
n = GT.num_dofs(fe)
@test A ≈ Matrix{Float64}(I,n,n)

end # module

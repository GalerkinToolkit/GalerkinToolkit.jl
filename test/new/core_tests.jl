module CoreTests

import GalerkinToolkit as GT
import StaticArrays as SA
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

Tv = GT.real_type(GT.options(cube2))
Ti = GT.int_type(GT.options(cube2))
Tr = GT.reference_int_type(GT.options(cube2))
space0 = GT.lagrange_space(cube0)
space1 = GT.lagrange_space(cube1)
space2 = GT.lagrange_space(cube2)
node_coordinates = SA.SVector{2,Tv}[(0,0),(1,0),(0,1),(1,1)]
face_nodes = Vector{Vector{Ti}}[[[1],[2],[3],[4]],[[1,2],[3,4],[1,3],[2,4]],[[1,2,3,4]]]
face_reference_id = Vector{Tr}[[1,1,1,1],[1,1,1,1],[1]]
outward_normals = SA.SVector{2,Tv}[(0,-1),(0,1),(-1,0),(1,0)]
reference_spaces = ((space0,),(space1,),(space2,))
mesh = GT.mesh(;
     node_coordinates,
     face_nodes,
     face_reference_id,
     reference_spaces,
     outward_normals
    )

d = 1
domain = GT.domain(mesh,d)

degree = 2
q = GT.quadrature(domain,degree)
q = GT.node_quadrature(domain)

face_point_x = GT.coordinate_accessor(q)
#@code_warntype GT.coordinate_accessor(q)
face_point_J = GT.jacobian_accessor(q)
face_point_dV = GT.weight_accessor(q)

face = 2
point = 1
x = face_point_x(face)(point)
J = face_point_J(face)(point)
dV = face_point_dV(face)(point,J)

face_nodes = Vector{Ti}[[1,2,3,4]]
face_reference_id = Tr[1]
reference_spaces = (space2,)
chain = GT.chain(;
     node_coordinates,
     face_nodes,
     face_reference_id,
     reference_spaces,
     outward_normals
    )

mesh = GT.mesh(chain)

chain2 = GT.chain(mesh)

#@code_warntype face_point_x(face)
#f = face_point_x(face)
#@code_warntype f(point)


end # module


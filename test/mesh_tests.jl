module MeshTests

import GalerkinToolkit as GT
import StaticArrays as SA
using Test
using LinearAlgebra

#using InteractiveUtils

#@code_warntype GT.unit_simplex(Val(3))

spx0 = GT.unit_simplex(0)
spx1 = GT.unit_simplex(1)
spx2 = GT.unit_simplex(2)
spx3 = GT.unit_simplex(3)

cube0 = GT.unit_n_cube(0)
cube1 = GT.unit_n_cube(1)
cube2 = GT.unit_n_cube(2)
cube3 = GT.unit_n_cube(3)

Tv = GT.real_type(GT.options(cube2))
Ti = GT.int_type(GT.options(cube2))
Tr = GT.reference_int_type(GT.options(cube2))
space0 = GT.lagrange_space(cube0,1)
space1 = GT.lagrange_space(cube1,1)
space2 = GT.lagrange_space(cube2,1)
node_coordinates = SA.SVector{2,Tv}[(0,0),(1,0),(0,1),(1,1)]
face_nodes = Vector{Vector{Ti}}[[[1],[2],[3],[4]],[[1,2],[3,4],[1,3],[2,4]],[[1,2,3,4]]]
face_reference_id = Vector{Tr}[[1,1,1,1],[1,1,1,1],[1]]
normals = SA.SVector{2,Tv}[(0,-1),(0,1),(-1,0),(1,0)]
reference_spaces = ((space0,),(space1,),(space2,))
mesh = GT.mesh(;
     node_coordinates,
     face_nodes,
     face_reference_id,
     reference_spaces,
     normals
    )

d = 1
domain = GT.domain(mesh,d)

degree = 2
q = GT.quadrature(domain,degree)
q = GT.node_quadrature(domain)

face = 2
point = 1
q_faces = GT.each_face(q)
q_points = GT.each_point(q_faces[face])
q_point = q_points[point]

x = GT.coordinate(q_point)
J = GT.jacobian(q_point)
dV = GT.weight(q_point)

face_nodes = Vector{Ti}[[1,2,3,4]]
face_reference_id = Tr[1]
reference_spaces = (space2,)
chain = GT.chain(;
     node_coordinates,
     face_nodes,
     face_reference_id,
     reference_spaces,
     normals
    )

mesh = GT.mesh(chain)

chain2 = GT.chain(mesh)

end # module

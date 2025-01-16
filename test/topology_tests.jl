module TopologyTests

import GalerkinToolkit as GT
import StaticArrays as SA
using Test

#using InteractiveUtils

spx0 = GT.unit_simplex(0)
spx1 = GT.unit_simplex(1)
spx2 = GT.unit_simplex(2)
spx3 = GT.unit_simplex(3)

cube0 = GT.unit_n_cube(0)
cube1 = GT.unit_n_cube(1)
cube2 = GT.unit_n_cube(2)
cube3 = GT.unit_n_cube(3)

mesh = GT.mesh(spx0)
mesh = GT.mesh(spx1)
mesh = GT.mesh(spx2)
mesh = GT.mesh(spx3)

mesh = GT.mesh(cube0)
mesh = GT.mesh(cube1)
mesh = GT.mesh(cube2)
mesh = GT.mesh(cube3)

order = 1
fe = GT.lagrange_space(cube2,order)
mesh = GT.complexify(fe)

order = 2
fe = GT.lagrange_space(cube2,order)
mesh = GT.complexify(fe)

order = 3
fe = GT.lagrange_space(cube2,order)
mesh = GT.complexify(fe)

display(GT.face_nodes(mesh,0))
display(GT.face_nodes(mesh,1))
display(GT.face_nodes(mesh,2))

#@code_warntype GT.mesh(cube0)

#@code_warntype GT.complexify(fe)

topo = GT.topology(cube0)
@test GT.num_dims(topo) == 0
@test GT.num_faces(topo,0) == 1
GT.face_incidence(topo,0,0)
GT.face_permutation_ids(topo,0,0)
@test GT.face_reference_id(topo,0) == [1]
@show GT.reference_topologies(topo,0)

mesh = GT.mesh(cube0)
mesh = GT.mesh(cube1)
topo = GT.topology(mesh)

mesh = GT.mesh(cube2)
topo = GT.topology(mesh)
mesh = GT.mesh(cube3)
topo = GT.topology(mesh)

mesh = GT.mesh(spx0)
mesh = GT.mesh(spx1)
topo = GT.topology(mesh)
mesh = GT.mesh(spx2)
topo = GT.topology(mesh)
mesh = GT.mesh(spx3)
topo = GT.topology(mesh)

topo = GT.topology(spx0)
topo = GT.topology(spx1)
topo = GT.topology(spx2)
topo = GT.topology(spx3)

@show GT.reference_topologies(topo,0)
@show GT.reference_topologies(topo,1)
@show GT.reference_topologies(topo,2)
@show GT.reference_topologies(topo,3)

Tv = GT.real_type(GT.options(spx2))
Ti = GT.int_type(GT.options(spx2))
Tr = GT.reference_int_type(GT.options(spx2))
space2 = GT.lagrange_space(spx2,1)
node_coordinates = SA.SVector{2,Tv}[(0,0),(1,0),(0,1),(1,1)]
face_nodes = Vector{Ti}[[1,2,3],[2,3,4]]
face_reference_id = Tr[1,1]
reference_spaces = (space2,)
chain = GT.chain(;
     node_coordinates,
     face_nodes,
     face_reference_id,
     reference_spaces,
    )

mesh = GT.mesh(chain)

mesh2 = GT.complexify(mesh)

@show GT.workspace(mesh2)

topo = GT.topology(mesh2)
topo2 = GT.topology(mesh2)

@test topo === topo2


end # module

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
@test topo == deepcopy(topo)

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
@test GT.is_cell_complex(mesh2)

@test GT.node_coordinates(mesh2) === GT.node_coordinates(mesh)

@show GT.num_faces(mesh2)

@test length(GT.reference_spaces(mesh2,2)) == 1
@test length(GT.reference_spaces(mesh2,1)) == 1
@test length(GT.reference_spaces(mesh2,0)) == 1

@show GT.face_reference_id(mesh2,0)
@show GT.face_reference_id(mesh2,1)
@show GT.face_reference_id(mesh2,2)
@show GT.face_nodes(mesh2,1)
@show GT.face_nodes(mesh2,0)
@show GT.face_nodes(mesh2,2)

@show GT.normals(mesh2)
@show GT.periodic_nodes(mesh2)

@show GT.group_faces(mesh2,2)
@show GT.group_faces(mesh2,1)
@show GT.group_faces(mesh2,0)

@show GT.workspace(mesh2)

topo = GT.topology(mesh2)
topo2 = GT.topology(mesh2)

@test topo === topo2

glue = GT.complexify_glue(GT.workspace(topo))

@show GT.periodic_faces_permutation_id(topo,0)
@show GT.periodic_faces_permutation_id(topo,1)
@show GT.periodic_faces_permutation_id(topo,2)

@show GT.periodic_faces(topo,0)
@show GT.periodic_faces(topo,1)
@show GT.periodic_faces(topo,2)

@show GT.face_permutation_ids(topo,2,0)
@show GT.face_permutation_ids(topo,2,1)
@show GT.face_permutation_ids(topo,1,0)

@show GT.face_reference_id(topo,0)
@show GT.face_reference_id(topo,1)
@show GT.face_reference_id(topo,2)

@show GT.face_incidence(topo,0,1)
@show GT.face_incidence(topo,1,2)
@show GT.face_incidence(topo,1,0)
@show GT.face_incidence(topo,2,1)
@show GT.face_incidence(topo,0,0)
@show GT.face_incidence(topo,2,2)
@show GT.face_incidence(topo,0,2)
@show GT.face_incidence(topo,2,0)

@show GT.parent_face_face(topo,1)
@show GT.parent_face_face(topo,2)
@show GT.parent_face_face(topo,0)
@show GT.face_nodes(mesh2,2)
@show GT.vertex_parent_faces(glue,0)
@show GT.vertex_parent_faces(glue,1)
@show GT.vertex_parent_faces(glue,2)
@show GT.parent_face_vertices(glue,0)
@show GT.parent_face_vertices(glue,1)
@show GT.parent_face_vertices(glue,2)
@show GT.vertex_node(glue)
@show GT.node_vertex(glue)
@show GT.num_faces(topo,1)
@show GT.num_faces(topo,0)
@show GT.num_faces(topo,2)


#@show isassigned(topo.face_incidence,0+1,0+1)
#@show isassigned(topo.face_incidence,1+1,0+1)
#@show isassigned(topo.face_incidence,2+1,0+1)
#
#@show isassigned(topo.face_incidence,0+1,1+1)
#@show isassigned(topo.face_incidence,1+1,1+1)
#@show isassigned(topo.face_incidence,2+1,1+1)
#
#@show isassigned(topo.face_incidence,0+1,2+1)
#@show isassigned(topo.face_incidence,1+1,2+1)
#@show isassigned(topo.face_incidence,2+1,2+1)

domain = (0,1,0,1)
cells = (2,2)
mesh = GT.cartesian_mesh(domain,cells)
topo = GT.topology(mesh)

@show GT.face_reference_id(mesh,0)
@show GT.face_reference_id(mesh,1)
@show GT.face_reference_id(mesh,2)
@show GT.face_nodes(mesh,1)
@show GT.face_nodes(mesh,0)
@show GT.face_nodes(mesh,2)
@show GT.normals(mesh)
@show GT.periodic_nodes(mesh)
@show GT.group_faces(mesh,2)
@show GT.group_faces(mesh,1)
@show GT.group_faces(mesh,0)


GT.face_incidence(topo,0,0)
GT.face_incidence(topo,1,0)
GT.face_incidence(topo,2,0)

GT.face_incidence(topo,0,1)
GT.face_incidence(topo,1,1)
GT.face_incidence(topo,2,1)

GT.face_incidence(topo,0,2)
GT.face_incidence(topo,1,2)
GT.face_incidence(topo,2,2)

GT.face_permutation_ids(topo,2,0)
GT.face_permutation_ids(topo,2,1)
GT.face_permutation_ids(topo,1,0)

domain = (0,1,0,1,0,1)
cells = (10,10,10)
mesh = GT.cartesian_mesh(domain,cells)
topo = GT.topology(mesh)


#@show "column 0"
#@show isassigned(topo.face_incidence,0+1,0+1)
#@show isassigned(topo.face_incidence,1+1,0+1)
#@show isassigned(topo.face_incidence,2+1,0+1)
#@show isassigned(topo.face_incidence,3+1,0+1)
#
#@show "column 1"
#@show isassigned(topo.face_incidence,0+1,1+1)
#@show isassigned(topo.face_incidence,1+1,1+1)
#@show isassigned(topo.face_incidence,2+1,1+1)
#@show isassigned(topo.face_incidence,3+1,1+1)
#
#@show "column 2"
#@show isassigned(topo.face_incidence,0+1,2+1)
#@show isassigned(topo.face_incidence,1+1,2+1)
#@show isassigned(topo.face_incidence,2+1,2+1)
#@show isassigned(topo.face_incidence,3+1,2+1)
#
#@show "column 3"
#@show isassigned(topo.face_incidence,0+1,3+1)
#@show isassigned(topo.face_incidence,1+1,3+1)
#@show isassigned(topo.face_incidence,2+1,3+1)
#@show isassigned(topo.face_incidence,3+1,3+1)


GT.face_incidence(topo,0,0)
GT.face_incidence(topo,1,0)
GT.face_incidence(topo,2,0)
GT.face_incidence(topo,3,0)

GT.face_incidence(topo,0,1)
GT.face_incidence(topo,1,1)
GT.face_incidence(topo,2,1)
GT.face_incidence(topo,3,1)

GT.face_incidence(topo,0,2)
GT.face_incidence(topo,1,2)
GT.face_incidence(topo,2,2)
GT.face_incidence(topo,3,2)

GT.face_incidence(topo,0,3)
GT.face_incidence(topo,1,3)
GT.face_incidence(topo,2,3)
GT.face_incidence(topo,3,3)

end # module

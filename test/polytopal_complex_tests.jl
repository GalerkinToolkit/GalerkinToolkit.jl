module PolytopalComplexTests

using GalerkinToolkit
using Meshes
using Test

grid = CartesianGrid(2,2)
mesh = fe_mesh(grid)
poly = polytopal_complex(mesh)
@test polytopal_complex(mesh) === poly

face_incidence(poly,0,0)
face_incidence(poly,0,1)
face_incidence(poly,0,2)
face_incidence(poly,1,0)
face_incidence(poly,1,1)
face_incidence(poly,1,2)
face_incidence(poly,2,0)
face_incidence(poly,2,1)
face_incidence(poly,2,2)

groups = physical_groups(poly)

@test node_vertex(poly) == 1:num_nodes(poly)
@test vertex_node(poly) == 1:num_nodes(poly)





end # module

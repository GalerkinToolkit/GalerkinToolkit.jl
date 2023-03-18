module PolytopalComplexTests

using GalerkinToolkit
using Meshes
using Test
using WriteVTK
using MappedArrays

file = msh_file(@__DIR__,"gmsh","higher_order_2D.msh")
poly = polytopal_complex(file)

face_faces(poly,0,0)
face_faces(poly,0,1)
face_faces(poly,0,2)
face_faces(poly,1,0)
face_faces(poly,1,1)
face_faces(poly,1,2)
face_faces(poly,2,0)
face_faces(poly,2,1)
face_faces(poly,2,2)

grid = CartesianGrid(2,2)
mesh = fe_mesh(grid)
poly = polytopal_complex(mesh)
@test polytopal_complex(mesh) === poly

grid = CartesianGrid(2,2)
poly = polytopal_complex(grid)
face_faces(poly,0,0)
face_faces(poly,0,1)
face_faces(poly,0,2)
face_faces(poly,1,0)
face_faces(poly,1,1)
face_faces(poly,1,2)
face_faces(poly,2,0)
face_faces(poly,2,1)
face_faces(poly,2,2)

groups = physical_groups(poly)
@test node_vertex(poly) == 1:num_nodes(poly)
@test vertex_node(poly) == 1:num_nodes(poly)

grid = CartesianGrid(3,3,3)

face = first(grid)
face_faces(face,3,0)
face_faces(face,3,1)
face_faces(face,3,2)
face_faces(face,3,3)
face_faces(face,2,3)
face_faces(face,1,3)
face_faces(face,0,3)
face_faces(face,0,0)
face_faces(face,0,1)
face_faces(face,0,2)
face_faces(face,1,0)
face_faces(face,1,1)
face_faces(face,1,2)
face_faces(face,2,0)
face_faces(face,2,1)
face_faces(face,2,2)

poly = polytopal_complex(grid)

face_faces(poly,3,0)
face_faces(poly,3,1)
face_faces(poly,3,2)
face_faces(poly,3,3)
face_faces(poly,2,3)
face_faces(poly,1,3)
face_faces(poly,0,3)
face_faces(poly,0,0)
face_faces(poly,0,1)
face_faces(poly,0,2)
face_faces(poly,1,0)
face_faces(poly,1,1)
face_faces(poly,1,2)
face_faces(poly,2,0)
face_faces(poly,2,1)
face_faces(poly,2,2)

file = msh_file(@__DIR__,"gmsh","demo.msh")
poly = polytopal_complex(file)

face_faces(poly,3,0)
face_faces(poly,3,1)
face_faces(poly,3,2)
face_faces(poly,3,3)
face_faces(poly,2,3)
face_faces(poly,1,3)
face_faces(poly,0,3)
face_faces(poly,0,0)
face_faces(poly,0,1)
face_faces(poly,0,2)
face_faces(poly,1,0)
face_faces(poly,1,1)
face_faces(poly,1,2)
face_faces(poly,2,0)
face_faces(poly,2,1)
face_faces(poly,2,2)


end # module

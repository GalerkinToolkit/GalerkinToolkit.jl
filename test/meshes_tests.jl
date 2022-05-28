module MeshesTests

using GalerkinToolkit
using Meshes
using Test

grid = CartesianGrid(2,2)

mesh = fe_mesh(grid)

groups = physical_groups(mesh)
@show group_ids(groups)
@show group_names(groups)
@show face_nodes(mesh,0)
@show face_nodes(mesh,1)
@show face_nodes(mesh,2)
@show ref_faces(mesh,0)
@show ref_faces(mesh,1)
@show ref_faces(mesh,2)
@show face_ref_id(mesh,0)
@show face_ref_id(mesh,1)
@show face_ref_id(mesh,2)



end # module

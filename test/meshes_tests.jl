module MeshesTests

using GalerkinToolkit
using Meshes
using Test

grid = CartesianGrid(3,3)

mesh = fe_mesh(grid)

groups = physical_groups(mesh)
group_ids(groups)
group_names(groups)
@test group_faces(groups,3) == [4]
@test group_faces(groups,5) == [1,3,4]
@test face_nodes(mesh,0) == [[1],[4],[13],[16]]
@test face_nodes(mesh,1) == [[1,2],[5,1],[2,3],[3,4],[4,8],[9,5],[8,12],[14,13],[13,9],[15,14],[12,16],[16,15]]
@test face_nodes(mesh,2) == [[1,2,6,5],[2,3,7,6],[3,4,8,7],[5,6,10,9],[6,7,11,10],[7,8,12,11],[9,10,14,13],[10,11,15,14],[11,12,16,15]]
ref_faces(mesh,0)
ref_faces(mesh,1)
ref_faces(mesh,2)
face_ref_id(mesh,0)
face_ref_id(mesh,1)
face_ref_id(mesh,2)


end # module
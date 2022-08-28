module FEMeshTests

using GalerkinToolkit
using StaticArrays
using Test

mesh = FEMesh{3}(2)

@test ambient_dim(mesh) == 3
@test domain_dim(mesh) == 2
@test num_faces(mesh,0) == 0
@test num_faces(mesh,1) == 0
@test num_faces(mesh,2) == 0
@test has_hanging_nodes(mesh) == false
@test has_periodic_nodes(mesh) == false
@test node_coordinates(mesh) == SVector{3,Float64}[]
face_ref_id(mesh,0)
face_ref_id(mesh,1)
face_ref_id(mesh,2)
ref_faces(mesh,0)
ref_faces(mesh,1)
ref_faces(mesh,2)

groups = physical_groups(mesh)
physical_group!(groups,2,"foo")
physical_group!(groups,0,"bar")
physical_group_faces(groups,2,"foo")
@test physical_group_id(groups,2,"foo") == 1
@test physical_group_name(groups,2,1) == "foo"
@test physical_group_names(groups,1) == String[]
@test physical_group_names(groups,0) == ["bar"]
@test physical_group_names(groups,2) == ["foo"]
@test physical_group_ids(groups,1) == Int[]
@test physical_group_ids(groups,0) == [2]
@test physical_group_ids(groups,2) == [1]


end

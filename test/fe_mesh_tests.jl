module FEMeshTests

using GalerkinToolkit
using StaticArrays
using Test

mesh = SimpleFEMesh{SVector{3,Float64}}(VOID,2)

@test ambient_dim(mesh) == 3
@test domain_dim(mesh) == 2
@test num_faces(mesh,0) == 0
@test num_faces(mesh,1) == 0
@test num_faces(mesh,2) == 0
@test is_hanging(mesh) == false
@test is_periodic(mesh) == false
@test node_coordinates(mesh) == SVector{3,Float64}[]
face_ref_id(mesh,0)
face_ref_id(mesh,1)
face_ref_id(mesh,2)
ref_faces(mesh,0)
ref_faces(mesh,1)
ref_faces(mesh,2)

groups = physical_groups(mesh)
add_group!(groups,2,"foo")
add_group!(groups,0,"bar")
group_faces(groups,"foo")
@test group_id(groups,"foo") == 1
@test group_name(groups,1) == "foo"
@test group_names(groups) == ["bar","foo"]
@test group_ids(groups) == [1,2]
@test group_names(groups,1) == String[]
@test group_names(groups,0) == ["bar"]
@test group_names(groups,2) == ["foo"]
@test group_ids(groups,1) == Int[]
@test group_ids(groups,0) == [2]
@test group_ids(groups,2) == [1]


end
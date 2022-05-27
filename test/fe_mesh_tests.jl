module FEMeshTests

using GalerkinToolkit
using StaticArrays
using Test

mesh = SimpleFEMesh{SVector{3,Float64}}(void,2)

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
add_group!(groups,"foo")
group_faces(groups,"foo",0)
group_faces(groups,"foo",1)
group_faces(groups,"foo",2)
@test group_id(groups,"foo") == 1
@test group_name(groups,1) == "foo"
@test group_names(groups) == ["foo"]
@test group_ids(groups) == [1]


end

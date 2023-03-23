module TestMeshInterface

using GalerkinToolkit
using Meshes
using Test
using StaticArrays
using WriteVTK

ref_tri = Triangle(Point.([(0,0),(1,0),(0,1)]))
ref_qua = Quadrangle(Point.([(0,0),(1,0),(1,1),(0,1)]))

x = SVector{2,Float64}[(0,0),(1,0),(0,1),(1,1),(0,2),(1,2)]
c = [[1,2,3],[2,3,4],[3,4,6,5]]
i = [1,1,2]
r = (ref_tri,ref_qua)

mesh_node_coordinates = x
mesh_face_nodes = [[Int[]],[Int[]],c]
mesh_face_reference_id = [Int[],Int[],i]
mesh_reference_faces = ((),(),r)

mesh = GenericMesh(
                mesh_node_coordinates,
                mesh_face_nodes,
                mesh_face_reference_id,
                mesh_reference_faces
               )

@test dimension(mesh) == 2
@test embedded_dimension(mesh) == 2

groups_2 = Dict(1=>physical_group([3,1],"region"))
groups_1 = typeof(groups_2)()
groups_0 = typeof(groups_2)()
groups = [groups_0,groups_1,groups_2]
mesh_with_groups = set_phyisical_groups(mesh,groups)

d = dimension(mesh)
vtk_grid("mesh",vtk_args(mesh_with_groups,d)...) do vtk
    vtk_physical_groups!(vtk,d,mesh_with_groups)
end

end # module

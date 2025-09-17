module CartesianMeshTests

using Test

import GalerkinToolkit as GT
import PartitionedArrays as PA

#using InteractiveUtils

domain = (0,1,0,1,0,1)
cells = (2,2,2)
mesh = GT.cartesian_mesh(domain,cells)
smesh = GT.simplexify(mesh)

domain = (0,1,0,1)
cells = (2,2)

mesh = GT.cartesian_mesh(domain,cells,boundary=false)
mesh = GT.cartesian_mesh(domain,cells,simplexify=true)
mesh = GT.cartesian_mesh(domain,cells,boundary=false,simplexify=true)
mesh = GT.cartesian_mesh(domain,cells)

mesh = GT.cartesian_mesh(domain,cells)
GT.geometry_names(mesh)
doms = GT.geometries(mesh,1)

mesh_faces = GT.each_face(mesh)
group_faces1 = findall(mesh_faces) do mesh_face
    lnode_x = GT.node_coordinates(mesh_face)
    xm = sum(lnode_x) / length(lnode_x)
    xm[2] < 0.5
end
D = GT.num_dims(mesh)
GT.group_faces(mesh,D)["1"] = group_faces1

Ω = GT.interior(mesh)
∂Ω = GT.boundary(Ω)

Ω1 = GT.interior(mesh;group_names=["1"])
∂Ω1 = GT.boundary(Ω1)

parts = PA.DebugArray(LinearIndices((2,)))
mesh = GT.with_mesh_partitioner(;parts) do
    domain = (0,1,0,1)
    cells = (2,2)
    GT.cartesian_mesh(domain,cells)
end

@test GT.is_partitioned(mesh)

Ω = GT.interior(mesh)
@test GT.is_partitioned(Ω)



end #module

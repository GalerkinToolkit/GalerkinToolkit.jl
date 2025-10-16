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

display(GT.periodic_nodes(mesh))

println("===")
domain = (0,1,0,1)
cells = (2,2)
mesh = GT.cartesian_mesh(domain,cells;periodic=(true,false))
display(GT.periodic_nodes(mesh))
topo = GT.topology(mesh)
display(collect(enumerate(GT.face_nodes(mesh,0))))
display(collect(enumerate(GT.face_nodes(mesh,1))))
display(collect(enumerate(GT.face_nodes(mesh,2))))
#display(collect(enumerate(GT.face_incidence(topo,1,0))))
display(collect(enumerate(GT.face_incidence(topo,2,0))))


#display(GT.periodic_faces(topo,0))
#display(GT.periodic_faces(topo,1))

xxx

domain = (0,1,0,1,0,1)
cells = (2,2,2)
mesh = GT.cartesian_mesh(domain,cells;periodic=(false,true,false))
display(GT.periodic_nodes(mesh))

parts = PA.DebugArray(LinearIndices((2,)))
mesh = GT.with_mesh_partitioner(;parts) do
    domain = (0,1,0,1)
    cells = (2,2)
    GT.cartesian_mesh(domain,cells)
end

@test GT.is_partitioned(mesh)

Ω = GT.interior(mesh)
@test GT.is_partitioned(Ω)

mesh2 = GT.centralize(mesh)


end #module

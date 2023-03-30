module MeshTopologyTests

using GalerkinToolkit
using Meshes

p = Point(0,1,2)

topo = mesh_topology(p)

@show face_incidence(topo,0,0)
@show reference_faces(topo,0)

s = Segment(Meshes.Point.([(0,0),(1,0)]))

topo = mesh_topology(s)

t = Triangle(Meshes.Point.([(0,0),(1,0),(0,1)]))

topo = mesh_topology(t)

t = Tetrahedron(Meshes.Point.([(0,0,0),(1,0,0),(0,1,0),(0,0,1)]))

∂t = face_boundary(t)

#display(face_nodes(∂t,2))
#display(face_nodes(∂t,1))
topo = mesh_topology(∂t)

@show typeof(topo)

end # module

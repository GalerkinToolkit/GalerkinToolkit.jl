module TestMeshes

using GalerkinToolkit
using Meshes
using Test

p1 = Point(0,1,2)

@test dimension(p1) == 0
@test embedded_dimension(p1) == 3

t1 = Triangle(Meshes.Point.([(0,0),(1,0),(0,1)]))
@test dimension(t1) == 2
@test embedded_dimension(t1) == 2
@test face_nodes(t1,1) == [[1, 2], [2, 3], [3, 1]]

end # module

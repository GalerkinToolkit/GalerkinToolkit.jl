module GmshTests

using GalerkinToolkit
using Test

file = msh_file(@__DIR__,"gmsh","periodic.msh")
mesh = fe_mesh(file)

groups = physical_groups(mesh)
@test group_names(groups) == ["PL", "PR", "P", "R", "L", "S"]
@test group_ids(groups) == [1, 2, 5, 6, 8, 9]
@test group_dims(groups) == [0, 0, 1, 1, 1, 2]
@test group_names(groups,0) == ["PL", "PR"]
@test group_ids(groups,0) == [1, 2]
@test group_dims(groups,0) == [0, 0]
@test group_faces(groups,1) == [1,4]

end # module

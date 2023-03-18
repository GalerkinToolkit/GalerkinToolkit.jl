module GmshTests

using GalerkinToolkit
using Test

file = msh_file(@__DIR__,"gmsh","higher_order_2D.msh")
mesh = fe_mesh(file)

file = msh_file(@__DIR__,"gmsh","twoTetraeder.msh")
mesh = fe_mesh(file)

file = msh_file(@__DIR__,"gmsh","t1.msh")
mesh = fe_mesh(file)

file = msh_file(@__DIR__,"gmsh","periodic.msh")
mesh = fe_mesh(file)

groups = physical_groups(mesh)
@test physical_group_names(groups,0) == ["PL", "PR"]
@test physical_group_ids(groups,0) == [1, 2]
@test physical_group_faces(groups,0,1) == [1,4]


end # module

module P4estTests

using GalerkinToolkit
using Meshes
using WriteVTK
using LinearAlgebra

dir = mkpath(joinpath(@__DIR__,"vtk"))

grid = CartesianGrid(1,1)
amr = p4est_mesh_refiner(grid,initial_level=1)
D = domain_dim(grid)

flags = [true,false,false,false]
refine!(amr,flags,num_levels=4)


mesh = fe_mesh(amr)

vtk_grid(joinpath(dir,"p4est0"),vtk_args(mesh,0)...) do vtk
  physical_groups!(vtk,mesh,0)
end

vtk_grid(joinpath(dir,"p4est1"),vtk_args(mesh,1)...) do vtk
  physical_groups!(vtk,mesh,1)
end

vtk_grid(joinpath(dir,"p4est2"),vtk_args(mesh,2)...) do vtk
  physical_groups!(vtk,mesh,2)
end

flags = fill(false,7)
flags[4] = true
refine!(amr,flags)
mesh = fe_mesh(amr)
vtk_grid("tmp",vtk_args(mesh,2)...) do vtk end

amr = p4est_mesh_refiner(grid,initial_level=4)
let mesh = fe_mesh(amr)
  for i in 1:5
    node_to_x = node_coordinates(mesh)
    face_to_nodes = face_nodes(mesh,D)
    flags = map(face_to_nodes) do nodes
      xs = view(node_to_x,nodes)
      xm = sum(xs)/length(xs)
      R = 0.5
      norm(xm) < R
    end
    if i%2 == 0
      coarsen!(amr,flags,num_levels=2)
    else
      refine!(amr,flags)
    end
    mesh = fe_mesh(amr)
  end
  vtk_grid(joinpath(dir,"p4est"),vtk_args(mesh,D)...) do vtk end
end







end # module

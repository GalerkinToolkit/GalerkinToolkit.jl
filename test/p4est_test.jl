module P4estTests

using GalerkinToolkit
using Meshes
using WriteVTK
using LinearAlgebra

grid = CartesianGrid(1,1)
amr = p4est_amr(grid,initial_level=1)
D = domain_dim(grid)

flags = [true,false,false,false]
new_to_old = refine!(amr,flags)
@show new_to_old
flags = fill(false,7)
flags[4] = true

f2c = refine!(amr,flags)
@show f2c

f2c = balance!(amr)
@show f2c
aaaaa

new_to_old = new_to_old[]
mesh = fe_mesh(amr)
vtk_grid("tmp",vtk_args(mesh,2)...) do vtk
  vtk["old"] = new_to_old
end


aaaaa
flags = fill(false,7)
flags[4] = true
new_to_old = refine!(amr,flags)
mesh = fe_mesh(amr)
vtk_grid("tmp",vtk_args(mesh,2)...) do vtk
  vtk["old"] = new_to_old
end

aaaa


amr = p4est_amr(grid,initial_level=4)

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
      coarsen!(amr,flags)
    else
      refine!(amr,flags)
    end
    mesh = fe_mesh(amr)
  end
  vtk_grid("tmp",vtk_args(mesh,D)...) do vtk end
end







end # module

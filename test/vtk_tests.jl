module VTKTests

using GalerkinToolkit
using Meshes
using WriteVTK
using Test
using MappedArrays

grid = CartesianGrid(4,4)
mesh = fe_mesh(grid)

dir = mkpath(joinpath(@__DIR__,"vtk"))
fn = joinpath(dir,"ex")

vtk_grid(fn,vtk_args(mesh,2)...) do vtk end
vtk_grid(fn,vtk_args(mesh,1)...) do vtk end
vtk_grid(fn,vtk_args(mesh,0)...) do vtk end

vtk_grid(fn,vtk_args(mesh,1)...) do vtk
  physical_groups!(vtk,mesh,1)
end

vtk_grid(fn,vtk_args(mesh,2)...) do vtk
  ids = [3,5,9,2]
  vtk["group"] = ids[classify_nodes(mesh,ids)]
end

sphere = Sphere((0.,0.,0.), 1.)
grid = discretize(sphere, RegularDiscretization(50,50))
poly = polytopal_complex(grid)

x = node_coordinates(poly)
face_to_mask = mappedarray(face_nodes(poly,2)) do nodes
  xn = view(x,nodes)
  xm = sum(xn)/length(xn)
  xm[1]+xm[2]<0.7
end
groups = physical_groups(poly)
add_group!(groups,2,"foo")
group_faces!(groups,findall(face_to_mask),2,"foo")

vtk_grid(fn,vtk_args(poly,2)...) do vtk
  physical_groups!(vtk,poly,2)
end

grid = CartesianGrid(4,4,4)
vtk_grid(fn,vtk_args(grid,3)...) do vtk
end


end # module

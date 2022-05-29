module VTKTests

using GalerkinToolkit
using Meshes
using WriteVTK
using Test

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

end # module

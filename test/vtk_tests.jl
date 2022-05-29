module VTKTests

using GalerkinToolkit
using Meshes
using WriteVTK

grid = CartesianGrid(4,4)
mesh = fe_mesh(grid)

dir = mkpath(joinpath(@__DIR__,"vtk"))
fn = joinpath(dir,"ex")

vtk_grid(fn,vtk_args(mesh,2)...) do vtk end
vtk_grid(fn,vtk_args(mesh,1)...) do vtk end
vtk_grid(fn,vtk_args(mesh,0)...) do vtk end

end # module

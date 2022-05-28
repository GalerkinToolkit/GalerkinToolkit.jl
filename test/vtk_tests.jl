module VTKTests

using GalerkinToolkit
using Meshes
using WriteVTK

grid = CartesianGrid(4,4)
mesh = fe_mesh(grid)

vtk_grid("filename", vtk_args(mesh,1)...) do vtk
end

end # module

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

# To think about this
@test_broken begin
groups = physical_groups(mesh)
vtk_grid(fn,vtk_args(mesh,1)...) do vtk
  physical_groups!(vtk,groups,1)
end
true
end

end # module

module TestGmsh

using GalerkinToolkit
using WriteVTK
using Test

# TODO: compare generated files with a reference

function write_mesh(mesh,filebase)
    for d in 0:dimension(mesh)
        fn = "$(filebase)_$(d)"
        vtk_grid(fn,vtk_args(mesh,d)...) do vtk
            vtk_physical_groups!(vtk,d,mesh)
        end
    end
end

file = joinpath(@__DIR__,"..","assets","twoTetraeder.msh")
mesh = gmsh_mesh(file)
write_mesh(mesh,"twoTetraeder")

file = joinpath(@__DIR__,"..","assets","t1.msh")
mesh = gmsh_mesh(file)
write_mesh(mesh,"t1")

file = joinpath(@__DIR__,"..","assets","periodic.msh")
mesh = gmsh_mesh(file)
@test length(periodic_nodes(mesh)) != 0
write_mesh(mesh,"periodic")

file = joinpath(@__DIR__,"..","assets","higher_order_2D.msh")
mesh = gmsh_mesh(file)
write_mesh(mesh,"higher_order_2D")

end # module

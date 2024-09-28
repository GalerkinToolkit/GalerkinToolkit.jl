module VisualizationTests

import GalerkinToolkit as GT
using GLMakie
using WriteVTK
using PartitionedArrays

domain = (0,1,0,1,0,1)
cells = (2,2,2)
mesh = GT.cartesian_mesh(domain,cells,simplexify=true)
np = 2
parts = DebugArray(LinearIndices((np,)))
pmesh = GT.partition_mesh(mesh,np;parts)

plt = GT.plot(pmesh)

fig = Makie.plot(plt;color=GT.FaceData("__OWNER__"),edgecolor=:black)
display(fig)

xxx

for mesh2 in (mesh,pmesh)

    plt = GT.plot(mesh2)
    vtk_grid("mesh",plt) |> close
    vtk_grid("mesh",mesh2) |> close

    plt = GT.shrink(plt,scale=0.7)
    vtk_grid("shrink",plt) |> close

    Ω = GT.interior(mesh2)
    u = GT.analytical_field(sum,Ω)

    plt = GT.plot(Ω)
    GT.plot!(plt,u;label="u")
    vtk_grid("domain",plt) |> close

    vtk_grid("domain",Ω) do plt
        GT.plot!(plt,u;label="u")
    end

    paraview_collection("domain",plt) do pvd
        vtk_grid("domain1",plt) do plt
            GT.plot!(plt,u;label="u")
            pvd[1] = plt
        end
    end

    paraview_collection("domain",Ω) do pvd
        vtk_grid("domain1",Ω) do plt
            GT.plot!(plt,u;label="u")
            pvd[1] = plt
        end
    end

end

Ω = GT.interior(mesh)
u = GT.analytical_field(sum,Ω)
plt = GT.plot(Ω)
GT.plot!(plt,u;label="u")
fig = Makie.plot(plt,color=GT.NodeData("u"))
display(fig)

fig = Makie.plot(Ω;color=u)
display(fig)

plt = GT.plot(mesh)
fig = GT.makie3d(plt;shading=Makie.NoShading,color=:blue)
GT.makie3d1d!(plt,color=:red)
GT.makie2d!(plt;shading=Makie.NoShading,color=:pink)
GT.makie2d1d!(plt,color=:cyan)
GT.makie1d!(plt,color=:green)
GT.makie0d!(plt,color=:black)
display(fig)

color = GT.FaceData("3-face-1")
fig = GT.makie3d1d(plt;color)
GT.makie3d!(plt;shading=Makie.NoShading,color=:blue)
GT.makie2d!(plt;shading=Makie.NoShading,color=:pink)
GT.makie2d1d!(plt;color)
GT.makie1d!(plt;color)
GT.makie0d!(plt;color)
display(fig)

fig = GT.makieplot(plt;dim=3)
display(fig)

plt = GT.restrict_to_dim(plt,3)
fig = GT.makieplot(plt)
display(fig)

fig = Makie.plot(plt;dim=3)
display(fig)

fig = Makie.plot(mesh;dim=2:3,shrink=0.6,edgecolor=:darkblue)
display(fig)



end # module

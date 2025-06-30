module VisualizationTests

import GalerkinToolkit as GT
using GLMakie
using WriteVTK
using PartitionedArrays

using StaticArrays

S = SVector{2,Float64}
vertices = S[(0,0),(1,0),(0,1),(1,1),(0,1),(1,1),(0,2),(1,2)]
conn = [1 2 3; 2 4 3; 5 7 6; 6 7 8]
fig = Makie.mesh(vertices,conn;color=:pink)
display(fig)

S = SVector{3,Float64}
vertices1 = S[(0,0,0),(1,0,0),(0,0,1),(1,0,1),(0,0,1),(1,0,1),(0,0,2),(1,0,2)]
vertices2 = vertices1 .+ Ref(S(0,0.1,0))
conn1 = [1 2 3; 2 4 3; 5 7 6; 6 7 8]
conn2 = [1 3 2; 2 3 4; 5 6 7; 6 8 7]
vertices = vcat(vertices1,vertices2)
conn = vcat(conn1,conn2 .+ length(vertices1))
fig = Makie.mesh(vertices,conn;color=:pink)
display(fig)

domain = (0,1,0,1,0,1)
n = 4
cells = (n,n,n)
mesh = GT.cartesian_mesh(domain,cells)

fig = Makie.plot(mesh;dim=0:3,shrink=0.6,strokecolor=:darkblue)
display(fig)
xxx

plt = GT.plot(mesh)


#plt = GT.shrink(plt;scale=0.7)
#plt = GT.restrict_to_dim(plt,2)
#vtk_grid("kk",plt) |> close

fig = GT.makie3d(plt)
GT.makie3d1d!(plt,color=:red)
GT.makie2d!(plt;color=:pink)
GT.makie2d1d!(plt,color=:cyan)
GT.makie1d!(plt,color=:green)
GT.makie0d!(plt,color=:black)
display(fig)

xxx




plt = GT.plot(mesh)
fig = Makie.plot(plt)
display(fig)

Ω = GT.interior(mesh)
fig = Makie.plot(Ω;color=:pink)
display(fig)

plt = GT.plot(Ω)
fig = Makie.plot(plt)
display(fig)

for s in  (false,true)

    domain = (0,1,0,1,0,1)
    cells = (4,4,4)
    mesh = GT.cartesian_mesh(domain,cells;simplexify=s)

    plt = GT.plot(mesh)

    fig = Makie.plot(plt;color=nothing,strokecolor=:black)
    display(fig)

    plt = GT.simplexify(plt)
    fig = Makie.plot(plt;color=:pink,strokecolor=:black)
    display(fig)

    np = 2
    parts = DebugArray(LinearIndices((np,)))
    pmesh = GT.partition_mesh(mesh,np;parts,renumber=true)
    vtk_grid("pmesh",pmesh) |> close

    plt = GT.plot(pmesh)

    #fig = Makie.plot(plt;color=GT.FaceData("__OWNER__"),strokecolor=:black)
    #display(fig)

    for mesh2 in (mesh,)#pmesh)

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

    #domain = (0,1,0,1)
    #cells = (4,4)
    #mesh = GT.cartesian_mesh(domain,cells;simplexify=s)


    Ω = GT.interior(mesh)
    u = GT.analytical_field(sum,Ω)
    v = GT.analytical_field(identity,Ω)
    plt = GT.plot(Ω)
    GT.plot!(plt,u;label="u")
    GT.plot!(plt,v;label="v")
    fig = Makie.plot(plt,color=GT.NodeData("u"))
    Makie.plot!(plt,color=nothing,strokecolor=:black,warp_by_vector=GT.NodeData("v"),warp_scale=0.1)
    #Makie.arrows2d!(plt,GT.NodeData("v"),lengthscale=0.1,color=GT.NodeData("u"))
    #Makie.arrows2d!(v;lengthscale=0.1,color=u)
    #display(fig)

    fig = Makie.plot(Ω;color=u)
    Makie.plot!(Ω;color=nothing,strokecolor=:black,warp_by_vector=v,warp_scale=0.1)
    display(fig)

    plt = GT.plot(mesh)
    fig = GT.makie3d(plt)
    GT.makie3d1d!(plt,color=:red)
    fig = GT.makie2d(plt;shading=Makie.NoShading)
    GT.makie2d1d!(plt,color=:cyan)
    GT.makie1d!(plt,color=:green)
    GT.makie0d!(plt,color=:black)
    display(fig)

    color = GT.FaceData("1-face-1")
    fig = GT.makie3d1d(plt;color)
    GT.makie3d!(plt;color=:blue)
    GT.makie2d!(plt;color=:pink)
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

    fig = Makie.plot(mesh;dim=2:3,shrink=0.6,strokecolor=:darkblue)
    display(fig)

    fig = Makie.plot(mesh;dim=0:3,shrink=0.6,strokecolor=:darkblue)
    display(fig)
    xxx

end


end # module

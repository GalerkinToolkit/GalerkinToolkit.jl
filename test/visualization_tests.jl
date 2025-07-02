module VisualizationTests

import GalerkinToolkit as GT
using GLMakie
using WriteVTK
using PartitionedArrays

using StaticArrays
using LinearAlgebra

# Setup we like in 2d
S = SVector{2,Float64}
vertices = S[(0,0),(1,0),(0,1),(1,1),(0,1),(1,1),(0,2),(1,2)]
conn = [1 2 3; 2 4 3; 5 7 6; 6 7 8]
edges = [1 2; 2 4; 3 4; 1 3; 4 8; 3 7; 7 8]
segments = zeros(S,2*size(edges,1))
segments[1:2:end-1] = vertices[edges[:,1]]
segments[2:2:end] = vertices[edges[:,2]]
directions = [ SVector(rand(),rand()) for _ in 1:length(vertices)]
fig = Figure()
ax = Axis(fig[1,1];aspect=DataAspect())
hidespines!(ax)
hidedecorations!(ax)
Makie.mesh!(vertices,conn)
Makie.scatter!(vertices;markersize=13,color=:black)
Makie.linesegments!(segments;linewidth=2,color=:black)
Makie.arrows2d!(vertices,directions,lengthscale = 0.5,color=norm.(directions))
display(fig)

# Setup we like in 3d
S = SVector{3,Float64}
vertices = S[(0,0,0),(1,0,0),(0,1,0),(1,1,0),(0,1,0),(1,1,0),(0,2,1),(1,2,1)]
directions = [ SVector(rand(),rand(),rand()) for _ in 1:length(vertices)]
faces = [1 2 3; 2 4 3; 5 6 7; 6 8 7]
edges = [1 2; 2 4; 3 4; 1 3; 4 8; 3 7; 7 8]
segments = zeros(S,2*size(edges,1))
segments[1:2:end-1] = vertices[edges[:,1]]
segments[2:2:end] = vertices[edges[:,2]]
fig = Figure()
#ax = Axis3(fig[1,1],aspect=:data)
#hidespines!(ax)
#hidedecorations!(ax)
ax = LScene(fig[1,1],show_axis=false)
Makie.mesh!(vertices,faces)
Makie.scatter!(vertices;markersize=13,color=:black)
Makie.linesegments!(segments;linewidth=2,color=:black)
Makie.arrows3d!(vertices,directions,lengthscale=0.5,color=norm.(directions))
display(fig)

# Visualizing plot objects
domain = (0,1,0,1,0,1)
n = 3
cells = (n,n,n)
mesh = GT.cartesian_mesh(domain,cells)
plt = GT.plot(mesh)
plt = GT.shrink(plt)
fig = Figure()
ax = Axis3(fig[1,1],aspect=:data)

GT.makie0d!(plt;dim=3,color=:red)
GT.makie0d!(plt;dim=2,color=:blue)
GT.makie0d!(plt;dim=1,color=:black)
GT.makie0d!(plt;dim=0,color=:green)
GT.makie1d!(plt;dim=3,color=:red)
GT.makie1d!(plt;dim=2,color=:blue)
GT.makie1d!(plt;dim=1,color=:black)
GT.makie2d!(plt;dim=2,color=:blue)
GT.makie2d!(plt;dim=3,color=:red)
display(fig)

fig = Figure()
ax = Axis3(fig[1,1],aspect=:data)
color = GT.FaceData("2-face-1")
GT.makie2d!(plt;dim=2,color)
GT.makie1d!(plt;dim=1,color)
GT.makie0d!(plt;dim=0,color)
display(fig)

# Visualizing plot objects
domain = (0,1,0,1)
n = 3
cells = (n,n)
mesh = GT.cartesian_mesh(domain,cells)
plt = GT.plot(mesh)
plt = GT.shrink(plt)
fig = Figure()
ax = Axis(fig[1,1],aspect=DataAspect())
GT.makie0d!(plt;dim=2,color=:blue)
GT.makie0d!(plt;dim=1,color=:black)
GT.makie0d!(plt;dim=0,color=:green)
GT.makie1d!(plt;dim=2,color=:blue)
GT.makie1d!(plt;dim=1,color=:black)
GT.makie2d!(plt;dim=2,color=:blue)
display(fig)

fig = Figure()
ax = Axis(fig[1,1],aspect=DataAspect())
color = GT.FaceData("1-face-1")
GT.makie2d!(plt;dim=2,color)
GT.makie1d!(plt;dim=1,color)
GT.makie0d!(plt;dim=0,color)
display(fig)

fig = Figure()
ax = Axis(fig[1,1],aspect=DataAspect())
warp_by_vector=fill(SVector(0.1,0.1),GT.num_nodes(plt.mesh))
warp_scale = 1
GT.makie2d!(plt)
GT.makie2d!(plt;warp_by_vector,warp_scale,color=:red)
GT.makie1d!(plt;warp_by_vector,warp_scale,color=:black)
GT.makie0d!(plt;warp_by_vector,warp_scale,color=:black)
display(fig)

fig = Figure()
ax = Axis3(fig[1,1],aspect=:data)
warp_by_scalar=fill(0.1,GT.num_nodes(plt.mesh))
warp_scale = 1
GT.makie2d!(plt)
GT.makie2d!(plt;warp_by_scalar,warp_scale,color=:red)
GT.makie1d!(plt;warp_by_scalar,warp_scale,color=:black)
GT.makie0d!(plt;warp_by_scalar,warp_scale,color=:black)
display(fig)

plt = GT.plot(mesh)
fig = Figure()
ax = Axis(fig[1,1],aspect=DataAspect())
shrink = 0.6
GT.makie0d!(plt;dim=2,shrink,color=:blue)
GT.makie0d!(plt;dim=1,shrink,color=:black)
GT.makie0d!(plt;dim=0,shrink,color=:green)
GT.makie1d!(plt;dim=2,shrink,color=:blue)
GT.makie1d!(plt;dim=1,shrink,color=:black)
GT.makie2d!(plt;dim=2,shrink,color=:blue)
display(fig)

GT.node_data(plt)["aux"] = rand(GT.num_nodes(plt.mesh))
fig = Figure()
ax = Axis(fig[1,1],aspect=DataAspect())
color = GT.NodeData("aux")
GT.makie2d!(plt;dim=2,color)
GT.makie1d!(plt;dim=1,color)
GT.makie0d!(plt;dim=0,color)
display(fig)

mesh = GT.mesh_from_msh(joinpath(@__DIR__,"..","assets","solid.msh"))
plt = GT.plot(mesh)
plt = GT.shrink(plt)
fig = Figure()
color = color=GT.FaceData("surface_1")
Axis3(fig[1,1],aspect=:data)
GT.makie2d!(plt;dim=2,color,colormap=:bluesreds)
GT.makie1d!(plt;dim=2,color,colormap=:bluesreds)
GT.makie0d!(plt;dim=2,color,colormap=:bluesreds)
display(fig)

plt = GT.plot(mesh)
plt = GT.skin(plt)
fig = Figure()
ax = Axis3(fig[1,1],aspect=:data)
hidespines!(ax)
hidedecorations!(ax)
#ax = LScene(fig[1,1],show_axis=false)
GT.makie2d!(plt;color=:pink)
GT.makie1d!(plt;color=:black,linewidth=2)
display(fig)

domain = (0,1,0,1,0,1)
n = 10
cells = (n,n,n)
mesh = GT.cartesian_mesh(domain,cells)

Ω = GT.interior(mesh)
u = GT.analytical_field(sum,Ω)
plt = GT.plot(Ω)
GT.plot!(plt,u;label="u")
fig = Figure()
ax = Axis3(fig[1,1],aspect=:data)
hidedecorations!(ax)
GT.makie2d!(plt;color=GT.NodeData("u"))
GT.makie1d!(plt)
GT.makie0d!(plt)
display(fig)

Γ = GT.boundary(mesh,physical_names=["2-face-1","2-face-3"])
plt = GT.plot(Γ)
u = GT.analytical_field(sum,Γ)
GT.plot!(plt,u;label="u")
fig = Figure()
ax = Axis3(fig[1,1],aspect=:data)
hidedecorations!(ax)
GT.makie2d!(plt;color=GT.NodeData("u"))
GT.makie1d!(plt)
GT.makie0d!(plt)
display(fig)


#xxx
#
#
#
#fig = Makie.plot(mesh;dim=3,axis=(aspect=:data,))
#
#Makie.plot!(mesh;draw=:face_edges)
#display(fig)
#xxx
#
#fig = Figure(size=(600,600))
#ax = Axis3(fig[1,1],aspect = :data)
#
#plt = GT.plot(mesh)
#Makie.plot!(plt;strokecolor=:black)
#
##Makie.plot!(mesh;fxaa=true,transparency=false,overdraw=false,color=:cyan,strokecolor=:darkblue,strokewidth=1.5)
#display(fig)
#xxx
#
#fig = Figure()
#ax = Axis3(fig[1,1])
#Makie.plot!(mesh;dim=0:3,shrink=0.6,strokecolor=:darkblue,strokewidth=2)
#display(fig)
#xxx
#
#plt = GT.plot(mesh)
#
#
##plt = GT.shrink(plt;scale=0.7)
##plt = GT.restrict_to_dim(plt,2)
##vtk_grid("kk",plt) |> close
#
#fig = GT.makie3d(plt)
#GT.makie3d1d!(plt,color=:red)
#GT.makie2d!(plt;color=:pink)
#GT.makie2d1d!(plt,color=:cyan)
#GT.makie1d!(plt,color=:green)
#GT.makie0d!(plt,color=:black)
#display(fig)
#
#xxx
#
#
#
#
#plt = GT.plot(mesh)
#fig = Makie.plot(plt)
#display(fig)
#
#Ω = GT.interior(mesh)
#fig = Makie.plot(Ω;color=:pink)
#display(fig)
#
#plt = GT.plot(Ω)
#fig = Makie.plot(plt)
#display(fig)
#
#for s in  (false,true)
#
#    domain = (0,1,0,1,0,1)
#    cells = (4,4,4)
#    mesh = GT.cartesian_mesh(domain,cells;simplexify=s)
#
#    plt = GT.plot(mesh)
#
#    fig = Makie.plot(plt;color=nothing,strokecolor=:black)
#    display(fig)
#
#    plt = GT.simplexify(plt)
#    fig = Makie.plot(plt;color=:pink,strokecolor=:black)
#    display(fig)
#
#    np = 2
#    parts = DebugArray(LinearIndices((np,)))
#    pmesh = GT.partition_mesh(mesh,np;parts,renumber=true)
#    vtk_grid("pmesh",pmesh) |> close
#
#    plt = GT.plot(pmesh)
#
#    #fig = Makie.plot(plt;color=GT.FaceData("__OWNER__"),strokecolor=:black)
#    #display(fig)
#
#    for mesh2 in (mesh,)#pmesh)
#
#        plt = GT.plot(mesh2)
#        vtk_grid("mesh",plt) |> close
#        vtk_grid("mesh",mesh2) |> close
#
#        plt = GT.shrink(plt,scale=0.7)
#        vtk_grid("shrink",plt) |> close
#
#        Ω = GT.interior(mesh2)
#        u = GT.analytical_field(sum,Ω)
#
#        plt = GT.plot(Ω)
#        GT.plot!(plt,u;label="u")
#        vtk_grid("domain",plt) |> close
#
#        vtk_grid("domain",Ω) do plt
#            GT.plot!(plt,u;label="u")
#        end
#
#        paraview_collection("domain",plt) do pvd
#            vtk_grid("domain1",plt) do plt
#                GT.plot!(plt,u;label="u")
#                pvd[1] = plt
#            end
#        end
#
#        paraview_collection("domain",Ω) do pvd
#            vtk_grid("domain1",Ω) do plt
#                GT.plot!(plt,u;label="u")
#                pvd[1] = plt
#            end
#        end
#
#    end
#
#    #domain = (0,1,0,1)
#    #cells = (4,4)
#    #mesh = GT.cartesian_mesh(domain,cells;simplexify=s)
#
#
#    Ω = GT.interior(mesh)
#    u = GT.analytical_field(sum,Ω)
#    v = GT.analytical_field(identity,Ω)
#    plt = GT.plot(Ω)
#    GT.plot!(plt,u;label="u")
#    GT.plot!(plt,v;label="v")
#    fig = Makie.plot(plt,color=GT.NodeData("u"))
#    Makie.plot!(plt,color=nothing,strokecolor=:black,warp_by_vector=GT.NodeData("v"),warp_scale=0.1)
#    #Makie.arrows2d!(plt,GT.NodeData("v"),lengthscale=0.1,color=GT.NodeData("u"))
#    #Makie.arrows2d!(v;lengthscale=0.1,color=u)
#    #display(fig)
#
#    fig = Makie.plot(Ω;color=u)
#    Makie.plot!(Ω;color=nothing,strokecolor=:black,warp_by_vector=v,warp_scale=0.1)
#    display(fig)
#
#    plt = GT.plot(mesh)
#    fig = GT.makie3d(plt)
#    GT.makie3d1d!(plt,color=:red)
#    fig = GT.makie2d(plt;shading=Makie.NoShading)
#    GT.makie2d1d!(plt,color=:cyan)
#    GT.makie1d!(plt,color=:green)
#    GT.makie0d!(plt,color=:black)
#    display(fig)
#
#    color = GT.FaceData("1-face-1")
#    fig = GT.makie3d1d(plt;color)
#    GT.makie3d!(plt;color=:blue)
#    GT.makie2d!(plt;color=:pink)
#    GT.makie2d1d!(plt;color)
#    GT.makie1d!(plt;color)
#    GT.makie0d!(plt;color)
#    display(fig)
#
#    fig = GT.makieplot(plt;dim=3)
#    display(fig)
#
#    plt = GT.restrict_to_dim(plt,3)
#    fig = GT.makieplot(plt)
#    display(fig)
#
#    fig = Makie.plot(plt;dim=3)
#    display(fig)
#
#    fig = Makie.plot(mesh;dim=2:3,shrink=0.6,strokecolor=:darkblue)
#    display(fig)
#
#    fig = Makie.plot(mesh;dim=0:3,shrink=0.6,strokecolor=:darkblue)
#    display(fig)
#    xxx
#
#end


end # module

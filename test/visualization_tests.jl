module VisualizationTests

import GalerkinToolkit as GT
using GLMakie
using WriteVTK
using PartitionedArrays

using StaticArrays
using LinearAlgebra


# mesh to plt

domain = (0,2,0,2)
cells = (4,4)
mesh = GT.cartesian_mesh(domain,cells)
parts = DebugArray(LinearIndices((4,)))
pmesh = GT.partitioned(mesh,parts)

plt = GT.plot(mesh)
pplt = GT.plot(pmesh)
#plts = GT.subdomains(pplt)
plt = GT.centralize(pplt)


# VTK

domain = (0,1,0,1,0,1)
cells = (4,4,4)
mesh = GT.cartesian_mesh(domain,cells)

plt = GT.plot(mesh)
vtk_grid("mesh",plt) |> close
vtk_grid("mesh",mesh) |> close

plt = GT.shrink(plt,scale=0.7)
vtk_grid("shrink",plt) |> close

Ω = GT.interior(mesh)
u = GT.analytical_field(sum,Ω)

plt = GT.plot(Ω)
GT.plot!(plt,u;label="u")
vtk_grid("domain",plt) |> close

vtk_grid("domain",Ω) do plt
    GT.plot!(plt,u;label="u")
end

paraview_collection("domain") do pvd
    vtk_grid("domain1",plt) do plt
        GT.plot!(plt,u;label="u")
        pvd[1] = plt
    end
end

paraview_collection("domain") do pvd
    vtk_grid("domain1",Ω) do plt
        GT.plot!(plt,u;label="u")
        pvd[1] = plt
    end
end

### VTK parallel

domain = (0,1,0,1)
cells = (10,10)
mesh = GT.cartesian_mesh(domain,cells)
parts = DebugArray(LinearIndices((4,)))
mesh = GT.partitioned(mesh,parts)

plt = GT.plot(mesh)
pvtk_grid("pmesh",plt) |> close
pvtk_grid("pmesh",mesh) |> close

paraview_collection("pmesh") do pvd
    pvtk_grid("pmesh",mesh) do plt
        pvd[1] = plt
    end
end

# Makie

counter = 0
dir = mkpath(joinpath(@__DIR__,"..","pictures"))

## Setup we like in 2d
#S = SVector{2,Float64}
#vertices = S[(0,0),(1,0),(0,1)]
#vertces = [0.0 0.0; 1.0 0.0; 0.0 1.0]
##conn = [1 2 3]
#conn = Int32[1 2 3]
#Makie.mesh(vertices,conn)
#counter += 1; Makie.save(joinpath(dir,"fig_$counter.png"),Makie.current_figure())

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
counter += 1; Makie.save(joinpath(dir,"fig_$counter.png"),Makie.current_figure())

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
counter += 1; Makie.save(joinpath(dir,"fig_$counter.png"),Makie.current_figure())

# Single-face domain
cube = GT.unit_n_cube(3)
axis = (;aspect=:data)
GT.makie_surfaces(cube;axis)
GT.makie_edges!(cube;color=:black)
counter += 1; Makie.save(joinpath(dir,"fig_$counter.png"),Makie.current_figure())

triangle = GT.unit_simplex(Val(2))
axis = (;aspect=Makie.DataAspect())
GT.makie_surfaces(triangle;axis)
GT.makie_edges!(triangle;color=:black)
counter += 1; Makie.save(joinpath(dir,"fig_$counter.png"),Makie.current_figure())

# Parallel 

domain = (0,1,0,1)
n = 10
cells = (n,n)
mesh = GT.cartesian_mesh(domain,cells)
np = 4
parts = DebugArray(LinearIndices((np,)))
mesh = GT.partitioned(mesh,parts)

plt = GT.plot(mesh)

fig = Figure()
ax = Axis3(fig[1,1],aspect=:data)
hidespines!(ax)
hidedecorations!(ax)
GT.makie_surfaces!(plt;color=GT.FaceColor("__OWNER__"))
GT.makie_edges!(plt;color=:black)
counter += 1; Makie.save(joinpath(dir,"fig_$counter.png"),Makie.current_figure())

fig = Figure()
ax = Axis3(fig[1,1],aspect=:data)
hidespines!(ax)
hidedecorations!(ax)
GT.makie_surfaces!(mesh;color=GT.FaceColor("__OWNER__"))
counter += 1; Makie.save(joinpath(dir,"fig_$counter.png"),Makie.current_figure())


# Parallel (old)

domain = (0,1,0,1,0,1)
n = 10
cells = (n,n,n)
mesh = GT.cartesian_mesh(domain,cells)
np = 5
parts = DebugArray(LinearIndices((np,)))
#pmesh = GT.partition_mesh(mesh,np;parts,renumber=true)
#
#pplt = GT.plot(pmesh)
#plt = GT.centralize(pplt)

fig = Figure()
ax = Axis3(fig[1,1],aspect=:data)
hidespines!(ax)
hidedecorations!(ax)
GT.makie_surfaces!(plt;color=GT.FaceColor("__OWNER__"))
GT.makie_edges!(plt;color=:black)
counter += 1; Makie.save(joinpath(dir,"fig_$counter.png"),Makie.current_figure())

#fig = Figure()
#ax = Axis3(fig[1,1],aspect=:data)
#hidespines!(ax)
#hidedecorations!(ax)
#GT.makie_surfaces!(pmesh;color=GT.FaceColor("__OWNER__"))
#counter += 1; Makie.save(joinpath(dir,"fig_$counter.png"),Makie.current_figure())

# Visualizing objects directly
domain = (0,1,0,1,0,1)
n = 3
cells = (n,n,n)
mesh = GT.cartesian_mesh(domain,cells)
fig = GT.makie_surfaces(mesh;axis=(;aspect=:data))
GT.makie_edges!(mesh;color=:black)
GT.makie_vertices!(mesh;color=:black)
counter += 1; Makie.save(joinpath(dir,"fig_$counter.png"),Makie.current_figure())

fig = Figure()
ax = Axis3(fig[1,1],aspect=:data)
hidespines!(ax)
hidedecorations!(ax)
Ω = GT.interior(mesh)
u = GT.analytical_field(sum,Ω)
v = GT.analytical_field(identity,Ω)
w = GT.analytical_field(x->-x,Ω)
GT.makie_surfaces!(Ω;color=u)
GT.makie_edges!(Ω;color=u)
GT.makie_vertices!(Ω;color=u)
GT.makie_arrows2d!(Ω,v;color=u)
GT.makie_arrows3d!(Ω,w;color=u)
counter += 1; Makie.save(joinpath(dir,"fig_$counter.png"),Makie.current_figure())

fig = GT.makie_surfaces(Ω;color=u)
GT.makie_edges!(Ω;warp_by_vector=v,color=u)
counter += 1; Makie.save(joinpath(dir,"fig_$counter.png"),Makie.current_figure())

# Visualizing plot objects
domain = (0,1,0,1,0,1)
n = 3
cells = (n,n,n)
mesh = GT.cartesian_mesh(domain,cells)
plt = GT.plot(mesh)
plt = GT.shrink(plt)

fig = Figure()
ax = Axis3(fig[1,1],aspect=:data)
GT.makie_vertices!(plt;dim=3,color=:red)
GT.makie_vertices!(plt;dim=2,color=:blue)
GT.makie_vertices!(plt;dim=1,color=:black)
GT.makie_vertices!(plt;dim=0,color=:green)
GT.makie_edges!(plt;dim=3,color=:red)
GT.makie_edges!(plt;dim=2,color=:blue)
GT.makie_edges!(plt;dim=1,color=:black)
GT.makie_surfaces!(plt;dim=2,color=:blue)
GT.makie_surfaces!(plt;dim=3,color=:red)
counter += 1; Makie.save(joinpath(dir,"fig_$counter.png"),Makie.current_figure())

fig = Figure()
ax = Axis3(fig[1,1],aspect=:data)
color = GT.FaceColor("2-face-1")
GT.makie_surfaces!(plt;dim=2,color)
GT.makie_edges!(plt;dim=1,color)
GT.makie_vertices!(plt;dim=0,color)
counter += 1; Makie.save(joinpath(dir,"fig_$counter.png"),Makie.current_figure())

# Visualizing plot objects
domain = (0,1,0,1)
n = 3
cells = (n,n)
mesh = GT.cartesian_mesh(domain,cells)
plt = GT.plot(mesh)
plt = GT.shrink(plt)
fig = Figure()
ax = Axis(fig[1,1],aspect=DataAspect())
GT.makie_vertices!(plt;dim=2,color=:blue)
GT.makie_vertices!(plt;dim=1,color=:black)
GT.makie_vertices!(plt;dim=0,color=:green)
GT.makie_edges!(plt;dim=2,color=:blue)
GT.makie_edges!(plt;dim=1,color=:black)
GT.makie_surfaces!(plt;dim=2,color=:blue)
counter += 1; Makie.save(joinpath(dir,"fig_$counter.png"),Makie.current_figure())

fig = Figure()
ax = Axis(fig[1,1],aspect=DataAspect())
color = GT.FaceColor("1-face-1")
GT.makie_surfaces!(plt;dim=2,color)
GT.makie_edges!(plt;dim=1,color)
GT.makie_vertices!(plt;dim=0,color)
counter += 1; Makie.save(joinpath(dir,"fig_$counter.png"),Makie.current_figure())

fig = Figure()
ax = Axis(fig[1,1],aspect=DataAspect())
warp_by_vector=fill(SVector(0.1,0.1),GT.num_nodes(plt.mesh))
warp_scale = 1
GT.makie_surfaces!(plt)
GT.makie_surfaces!(plt;warp_by_vector,warp_scale,color=:red)
GT.makie_edges!(plt;warp_by_vector,warp_scale,color=:black)
GT.makie_vertices!(plt;warp_by_vector,warp_scale,color=:black)
counter += 1; Makie.save(joinpath(dir,"fig_$counter.png"),Makie.current_figure())

fig = Figure()
ax = Axis3(fig[1,1],aspect=:data)
warp_by_scalar=fill(0.1,GT.num_nodes(plt.mesh))
warp_scale = 1
GT.makie_surfaces!(plt)
GT.makie_surfaces!(plt;warp_by_scalar,warp_scale,color=:red)
GT.makie_edges!(plt;warp_by_scalar,warp_scale,color=:black)
GT.makie_vertices!(plt;warp_by_scalar,warp_scale,color=:black)
counter += 1; Makie.save(joinpath(dir,"fig_$counter.png"),Makie.current_figure())

plt = GT.plot(mesh)
fig = Figure()
ax = Axis(fig[1,1],aspect=DataAspect())
shrink = 0.6
GT.makie_vertices!(plt;dim=2,shrink,color=:blue)
GT.makie_vertices!(plt;dim=1,shrink,color=:black)
GT.makie_vertices!(plt;dim=0,shrink,color=:green)
GT.makie_edges!(plt;dim=2,shrink,color=:blue)
GT.makie_edges!(plt;dim=1,shrink,color=:black)
GT.makie_surfaces!(plt;dim=2,shrink,color=:blue)
counter += 1; Makie.save(joinpath(dir,"fig_$counter.png"),Makie.current_figure())

GT.node_data(plt)["aux"] = rand(GT.num_nodes(plt.mesh))
fig = Figure()
ax = Axis(fig[1,1],aspect=DataAspect())
color = GT.NodeColor("aux")
GT.makie_surfaces!(plt;dim=2,color)
GT.makie_edges!(plt;dim=1,color)
GT.makie_vertices!(plt;dim=0,color)
counter += 1; Makie.save(joinpath(dir,"fig_$counter.png"),Makie.current_figure())

mesh = GT.mesh_from_msh(joinpath(@__DIR__,"..","assets","solid.msh"))
plt = GT.plot(mesh)
plt = GT.shrink(plt)
fig = Figure()
color = color=GT.FaceColor("surface_1")
Axis3(fig[1,1],aspect=:data)
GT.makie_surfaces!(plt;dim=2,color,colormap=:bluesreds)
GT.makie_edges!(plt;dim=2,color,colormap=:bluesreds)
GT.makie_vertices!(plt;dim=2,color,colormap=:bluesreds)
counter += 1; Makie.save(joinpath(dir,"fig_$counter.png"),Makie.current_figure())

plt = GT.plot(mesh)
plt = GT.skin(plt)
fig = Figure()
ax = Axis3(fig[1,1],aspect=:data)
hidespines!(ax)
hidedecorations!(ax)
#ax = LScene(fig[1,1],show_axis=false)
GT.makie_surfaces!(plt;color=:pink)
GT.makie_edges!(plt;color=:black,linewidth=2)
counter += 1; Makie.save(joinpath(dir,"fig_$counter.png"),Makie.current_figure())

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
GT.makie_surfaces!(plt;color=GT.NodeColor("u"))
GT.makie_edges!(plt)
GT.makie_vertices!(plt)
counter += 1; Makie.save(joinpath(dir,"fig_$counter.png"),Makie.current_figure())

Γ = GT.boundary(mesh,group_names=["2-face-1","2-face-3"])
plt = GT.plot(Γ)
u = GT.analytical_field(sum,Γ)
GT.plot!(plt,u;label="u")
fig = Figure()
ax = Axis3(fig[1,1],aspect=:data)
hidedecorations!(ax)
GT.makie_surfaces!(plt;color=GT.NodeColor("u"))
GT.makie_edges!(plt)
GT.makie_vertices!(plt)
counter += 1; Makie.save(joinpath(dir,"fig_$counter.png"),Makie.current_figure())

# Arrows 2d
v = GT.analytical_field(identity,Ω)
GT.plot!(plt,v;label="v")
fig = Figure()
ax = Axis3(fig[1,1],aspect=:data)
hidedecorations!(ax)
GT.makie_surfaces!(plt)
GT.makie_arrows2d!(plt,GT.NodeColor("v");color=GT.NodeColor("u"),lengthscale=0.1)
counter += 1; Makie.save(joinpath(dir,"fig_$counter.png"),Makie.current_figure())

# Arrows 3d
v = GT.analytical_field(identity,Ω)
GT.plot!(plt,v;label="v")
fig = Figure()
ax = Axis3(fig[1,1],aspect=:data)
hidedecorations!(ax)
GT.makie_surfaces!(plt)
GT.makie_arrows3d!(plt,GT.NodeColor("v");color=GT.NodeColor("u"),lengthscale=0.1)
counter += 1; Makie.save(joinpath(dir,"fig_$counter.png"),Makie.current_figure())


end # module

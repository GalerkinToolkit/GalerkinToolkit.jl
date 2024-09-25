module VisualizationTests

import GalerkinToolkit as GT
using GLMakie
using PartitionedArrays

#domain = (0,1,0,1)
#cells = (4,4)
#np = 4
#parts = identity(LinearIndices((np,)))
#parts = DebugArray(LinearIndices((np,)))
#mesh = GT.cartesian_mesh(domain,cells,simplexify=true)
#pmesh = GT.partition_mesh(mesh,np;parts)
#
#plt = GT.plot(pmesh)
#GT.save_vtk("pmesh",plt)
#GT.save_vtk("pmesh",pmesh)
#
#plt = GT.shrink(plt,scale=0.7)
#GT.save_vtk("shrink",plt)

#fig = GT.makie2d(plt;shading=Makie.NoShading,color=GT.FaceColor("__OWNER__"))
#display(fig)

domain = (0,1,0,1,0,1)
cells = (2,2,2)
mesh = GT.cartesian_mesh(domain,cells,simplexify=true)

plt = GT.plot(mesh)
GT.save_vtk("mesh",plt)
GT.save_vtk("mesh",mesh)

plt = GT.shrink(plt,scale=0.7)
GT.save_vtk("shrink",plt)

Ω = GT.interior(mesh)
u = GT.analytical_field(sum,Ω)

plt = GT.plotnew(Ω;fields=(;u))
GT.save_vtk("domain",plt)
fig = Makie.plot(plt,color=GT.NodeColor("u"))
display(fig)

fig = Makie.plot(Ω;color=u)
display(fig)


#using StaticArrays
#x = SVector{2,Float64}[[0,0],[1,1]]
#fig = Makie.linesegments(x,color=[1,2])
#display(fig)
#
#xxx

plt = GT.plot(mesh)
#fig = Figure()
#ax = Axis(fig[1,1], aspect=DataAspect())
fig = GT.makie3d(plt;shading=Makie.NoShading,color=:blue)
GT.makie3d1d!(plt,color=:red)
GT.makie2d!(plt;shading=Makie.NoShading,color=:pink)
GT.makie2d1d!(plt,color=:cyan)
GT.makie1d!(plt,color=:green)
GT.makie0d!(plt,color=:black)
display(fig)

color = GT.FaceColor("3-face-1")
fig = GT.makie3d1d(plt;color)
GT.makie3d!(plt;shading=Makie.NoShading,color=:blue)
GT.makie2d!(plt;shading=Makie.NoShading,color=:pink)
GT.makie2d1d!(plt;color)
GT.makie1d!(plt;color)
GT.makie0d!(plt;color)
display(fig)

#plt = GT.restrict_to_dim(plt,3)
fig = GT.makieplot(plt;dim=3)
display(fig)

plt = GT.restrict_to_dim(plt,3)
fig = GT.makieplot(plt)
display(fig)

fig = Makie.plot(plt;dim=3)
display(fig)

fig = Makie.plot(mesh;dim=3,shrink=0.6)
display(fig)




#fig = Figure()
#ax = Axis(fig[1,1], aspect=DataAspect())
#GT.render_with_makie!(ax,plt;color=GT.FaceColor("1-face-1")#)
#hidedecorations!(ax)
#hidespines!(ax)
#display(fig)

#fig = GT.makieplot(plt;color=GT.FaceColor("2-face-1"),alpha=0.1)

#display(fig)



#plt = GT.plot(Ω)
#GT.plot!(plt,uh;label="uh")
#GT.plot!(plt,uh;label="uh")
#GT.plot!(plt,uh;label="uh")
#GT.plot!(plt,uh;label="uh")
#save_vtk(plt)


#makie_render(plt,render=:surface,showedges=true,color=GT.NodeColor("0-face-1"),points3d=true,edges3d=true)

#plt = visualize(mtf,render=:surface,showedges=true,color=GT.NodeColor("0-face-1"),points3d=true,edges3d=true)

end # module

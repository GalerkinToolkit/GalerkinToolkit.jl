module VisualizationTests

import GalerkinToolkit as GT
using GLMakie

domain = (0,1,0,1,0,1)
cells = (2,2,2)
mesh = GT.cartesian_mesh(domain,cells,simplexify=true)

plt = GT.plot(mesh)
GT.plot!(plt,GT.physical_faces)
GT.save_vtk("mesh",plt)
GT.save_vtk("mesh",mesh)

plt = GT.shrink(plt,scale=0.7)
#GT.save_vtk("shrink",plt)

#fig = Figure()
#ax = Axis(fig[1,1], aspect=DataAspect())
fig = GT.makievolumes(plt;shading=Makie.NoShading)
GT.makievolumeedges!(plt,color=:red)
#GT.makiefaceedges!(plt,color=:red)
#GT.makieedges!(plt)
#GT.makievertices!(plt)
display(fig)


#fig = Figure()
#ax = Axis(fig[1,1], aspect=DataAspect())
#GT.render_with_makie!(ax,plt;color=GT.FaceColor("1-face-1"))
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

module VisualizationTests

import GalerkinToolkit as GT
import Makie
using GLMakie

domain = (0,1,0,1)
cells = (3,3)
mesh = GT.cartesian_mesh(domain,cells,simplexify=true,complexify=false)

plt = GT.plot(mesh)
GT.plot!(plt,GT.physical_faces)
GT.save_vtk("mesh",plt)
GT.save_vtk("mesh",mesh)

plt = GT.shrink(plt)

fig = Figure()
ax = Axis(fig, aspect=DataAspect())
GT.render_with_makie!(ax,plt;color=GT.FaceColor("1-face-1"))
display(fig)



#plt = GT.plot(Ω)
#GT.plot!(plt,uh;label="uh")
#GT.plot!(plt,uh;label="uh")
#GT.plot!(plt,uh;label="uh")
#GT.plot!(plt,uh;label="uh")
#save_vtk(plt)


#makie_render(plt,render=:surface,showedges=true,color=GT.NodeColor("0-face-1"),points3d=true,edges3d=true)

#plt = visualize(mtf,render=:surface,showedges=true,color=GT.NodeColor("0-face-1"),points3d=true,edges3d=true)

end # module

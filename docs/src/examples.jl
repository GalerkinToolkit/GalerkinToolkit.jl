# # Examples

# ## Poisson
#
# Import dependencies.

import GalerkinToolkit as GT
import Gmsh: gmsh
import GLMakie
import ForwardDiff
using LinearAlgebra

# Define geometry via gmsh API.

#if gmsh.isInitialized() == 1
#    gmsh.finalize()
#end
#gmsh.initialize()
#gmsh.option.setNumber("General.Terminal",0)
#v1 = gmsh.model.occ.addBox(0,0,0,1,1,1)
#c1 = gmsh.model.occ.addCylinder(0,0.5,0.5,1,0,0,0.4)
#c2 = gmsh.model.occ.addCylinder(0.5,0,0.5,0,1,0,0.4)
#v2, = gmsh.model.occ.fuse([(3,c1)],[(3,c2)])
#v3, = gmsh.model.occ.cut([(3,v1)],v2)
#gmsh.model.occ.synchronize()
#gmsh.model.addPhysicalGroup(3,map(last,v3),-1,"volume")
#gmsh.model.addPhysicalGroup(2,[14],-1,"top")
#gmsh.model.addPhysicalGroup(2,[15],-1,"bottom")
#gmsh.model.mesh.generate(3)
#gmsh.model.mesh.renumberElements()
#gmsh.model.mesh.renumberNodes()
#
## Create a mesh from the gmsh-generated mesh.
#
#mesh =  GT.mesh_from_gmsh_module()
#GLMakie.plot(mesh,strokecolor=:black)


mesh = GT.mesh_from_gmsh("assets/mesh1.msh")

# Solve the Poisson equation with GalerkinToolkit.

k = 1
Ω = GT.interior(mesh)
Γd1 = GT.boundary(mesh;physical_names=["top"])
Γd2 = GT.boundary(mesh;physical_names=["bottom"])


g1 = GT.analytical_field(x->1,Ω)
g2 = GT.analytical_field(x->2,Ω)
g = GT.piecewiese_field(g1,g2)
Γd = GT.piecewiese_domain(Γd1,Γd2)
V = GT.lagrange_space(Ω,k;dirichlet_boundary=Γd)
dΩ = GT.measure(Ω,2*k)
uhd = GT.dirichlet_field(Float64,V)
GT.interpolate_dirichlet!(g,uhd)
gradient(u) = x->ForwardDiff.gradient(u,x)
∇(u,x) = GT.call(gradient,u)(x)
a(u,v) = GT.∫( x->∇(u,x)⋅∇(v,x), dΩ)
l(v) = 0
x,A,b = GT.linear_problem(uhd,a,l)
x .= A\b
uh = GT.solution_field(uhd,x)
GLMakie.plot(Ω,color=uh,strokecolor=:black)



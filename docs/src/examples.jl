# # Examples

# ## Poisson
#
# Import dependencies.

import GalerkinToolkit as GT
import GLMakie
import ForwardDiff
using LinearAlgebra

# Load a mesh from a msh file.

## TODO a more elegant way of getting the mesh
mesh = GT.mesh_from_gmsh("../../assets/mesh1.msh")
GLMakie.plot(mesh,strokecolor=:black)

# This mesh has defined to physical groups for surfaces, "top" and "bottom".

GT.physical_names(mesh,2)


# Solve the Poisson equation imposing Dirichlet boundary conditions at the top and bottom faces.

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


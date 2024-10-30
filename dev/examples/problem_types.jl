# # Problem types
#
# ## Linear, steady-state, single-field

# Solve the following Poisson equation on the unit square,
#
# ```math
# \left\lbrace
# \begin{aligned}
# -\Delta u = f  \ &\text{in} \ \Omega,\\
# u = g \ &\text{on}\ \partial\Omega,\\
# \end{aligned}
# \right.
# ```
# with $f=1$ and $g=0$.
#
# Solve it with a piece-wise linear Lagrange interpolation, and visualize the result.

import GalerkinToolkit as GT
import ForwardDiff
import GLMakie
using LinearAlgebra
domain = (0,1,0,1)
cells = (10,10)
mesh = GT.cartesian_mesh(domain,cells;simplexify=true)
dirichlet_tag = "dirichlet"
GT.label_boundary_faces!(mesh;physical_name=dirichlet_tag)
Ω = GT.interior(mesh)
Γd = GT.boundary(mesh;physical_names=[dirichlet_tag])
k = 1
V = GT.lagrange_space(Ω,k;dirichlet_boundary=Γd)
uhd = GT.dirichlet_field(Float64,V)
g = GT.analytical_field(x->0,Ω)
f = GT.analytical_field(x->1,Ω)
GT.interpolate_dirichlet!(g,uhd)
dΩ = GT.measure(Ω,2*k)
gradient(u) = x->ForwardDiff.gradient(u,x)
∇(u,x) = GT.call(gradient,u)(x)
a(u,v) = GT.∫( x->∇(u,x)⋅∇(v,x), dΩ)
l(v) = GT.∫( x->v(x)*f(x), dΩ)
x,A,b = GT.linear_problem(uhd,a,l)
x .= A\b
uh = GT.solution_field(uhd,x)
GLMakie.plot(Ω;color=uh,strokecolor=:black)

# !!! warning
#     TODO:
#     - 2D domains should be visualized as 2D plots by default
#     - Transparent background so that figures look good in dark mode.

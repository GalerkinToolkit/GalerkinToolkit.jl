# # p-Laplacian
#
# ## Problem statement
#
# Solve the following p-Laplacian equation on the unit square,
#
# ```math
# \left\lbrace
# \begin{aligned}
# -\nabla \cdot \left( |\nabla u|^{q-2} \ \nabla u \right) = f\ &\text{in}\ \Omega,\\
# u = g \ &\text{on}\ \partial\Omega,\\
# \end{aligned}
# \right.
# ```
# with $f=1$ and $g=0$ and $q=3$.
#
# Solve it with a piece-wise bi-linear Lagrange interpolation, and visualize the result.
#

# ## Implementation

# Load dependencies form Julia stdlib.

using LinearAlgebra

# Import other dependencies

import GalerkinToolkit as GT
import PartitionedSolvers as PS
import ForwardDiff
import GLMakie as Makie
import FileIO # hide

# Code

domain = (0,1,0,1)
cells = (10,10)
mesh = GT.cartesian_mesh(domain,cells)
dirichlet_tag = "dirichlet"
GT.label_boundary_faces!(mesh;physical_name=dirichlet_tag)
Ω = GT.interior(mesh)
Γd = GT.boundary(mesh;physical_names=[dirichlet_tag])
k = 1
V = GT.lagrange_space(Ω,k;dirichlet_boundary=Γd)
g = GT.analytical_field(x->0,Ω)
f = GT.analytical_field(x->1,Ω)
const q = 3
flux(∇u) = norm(∇u)^(q-2) * ∇u
dflux(∇du,∇u) = (q-2)*norm(∇u)^(q-4)*(∇u⋅∇du)*∇u+norm(∇u)^(q-2)*∇du
uh = GT.rand_field(Float64,V)
GT.interpolate_dirichlet!(g,uh)
dΩ = GT.measure(Ω,2*k)
∇ = ForwardDiff.gradient
res(u) = v -> GT.∫( x-> ∇(v,x)⋅GT.call(flux,∇(u,x)) - f(x)*v(x) , dΩ)
jac(u) = (du,v) -> GT.∫( x-> ∇(v,x)⋅GT.call(dflux,∇(du,x),∇(u,x)) , dΩ)
p = GT.nonlinear_problem(uh,res,jac)
s = PS.newton_raphson(p,verbose=true)
s = PS.solve(s)
uh = GT.solution_field(uh,s)
Makie.plot(Ω;color=uh,strokecolor=:black)
FileIO.save(joinpath(@__DIR__,"fig_pt_plaplacian.png"),Makie.current_figure()) # hide

# ![](fig_pt_plaplacian.png)

# Now, by showing the intermediate results in the iteration process

uh = GT.rand_field(Float64,V)
GT.interpolate_dirichlet!(g,uh)
p = GT.nonlinear_problem(uh,res,jac)
s = PS.newton_raphson(p)
color = Makie.Observable(uh)
fig = Makie.plot(Ω;color,strokecolor=:black)
fn = joinpath(@__DIR__,"fig_pt_plaplacian.gif")
Makie.record(fig,fn,PS.history(s);framerate=2) do s
    color[] = GT.solution_field(uh,s)
end
nothing # hide

# ![](fig_pt_plaplacian.gif)


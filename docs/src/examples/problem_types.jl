# # Problem types
#

import GalerkinToolkit as GT
import PartitionedSolvers as PS
import ForwardDiff
import GLMakie as Makie
using LinearAlgebra
import FileIO # hide

# ## Poisson
#
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
# Solve it with a piece-wise bi-linear Lagrange interpolation, and visualize the result.
#

domain = (0,1,0,1)
cells = (10,10)
mesh = GT.cartesian_mesh(domain,cells)
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
∇ = ForwardDiff.gradient
a(u,v) = GT.∫( x->∇(u,x)⋅∇(v,x), dΩ)
l(v) = GT.∫( x->v(x)*f(x), dΩ)
p = GT.linear_problem(uhd,a,l)
s = PS.LinearAlgebra_lu(p)
s = PS.solve(s)
uh = GT.solution_field(uhd,s)
Makie.plot(Ω;color=uh,strokecolor=:black)
FileIO.save(joinpath(@__DIR__,"fig_pt_poisson.png"),Makie.current_figure()) # hide

# ![](fig_pt_poisson.png)

# !!! warning
#     * TODO The unit square should look square not rectangular
#     * TODO Use a more complex 2d geometry. The map of the Netherlands?
#     

# ## p-Laplacian
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

# ## Heat equation
#
#
# !!! warning
#     TODO Key missing things
#     * Implement theta method in PartitionedSolvers (easy)
#     * function `GT.ode_problem` (challenging)

# ## Wave equation
#
# As a 2nd order ODE
#
# !!! warning
#     TODO Key missing things
#     * Implement Newmark method in PartitionedSolvers (easy)
#
# Reducing to a 1st order ODE
#
# !!! warning
#     TODO This should be quite straight forward once the functions above have been implemented
#
# ## Helmholtz equation
#
# !!! warning
#     TODO This will illustrate how to use complex numbers.
#     It should be possible to implement the example with the current state of the code.
#     (except for periodic boundary conditions)
#     Help wanted.
#
# ## Stokes equation
#
# !!! warning
#     TODO Lid cavity problem. This will illustrate how to use vector-valued spaces and multifield.
#     The key ingredient missing is to develop a strategy to impose a zero mean constrain on the pressure.
#
#

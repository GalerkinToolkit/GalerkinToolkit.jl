# # Hello, world!
#
# In this example, we show how to solve the "Hello, world" PDE example:
# the Poisson equation on the unit square with Dirichlet boundary conditions.
#
# ```math
# \left\lbrace
# \begin{aligned}
# -\Delta u = f  \ &\text{in} \ \Omega,\\
# u = g \ &\text{on}\ \partial\Omega,\\
# \end{aligned}
# \right.
# ```
# with $f=0$ and $g(x)=\text{sum}(x)$.
#
#

# ## Automatic assembly
#
# Load dependencies form Julia stdlib.

using LinearAlgebra

# Import other dependencies

import GalerkinToolkit as GT
import PartitionedSolvers as PS
import ForwardDiff
import GLMakie as Makie
import FileIO # hide

# Generate the computational mesh.

domain = (0,1,0,1)
cells = (10,10)
mesh = GT.cartesian_mesh(domain,cells)
nothing # hide

# Visualize the mesh.

Makie.plot(mesh,color=:pink,strokecolor=:blue)
FileIO.save(joinpath(@__DIR__,"fig_poisson_1.png"),Makie.current_figure()) # hide

# ![](fig_poisson_1.png)

# Define the Dirichlet boundary.

dirichlet_tag = "dirichlet"
GT.label_boundary_faces!(mesh;physical_name=dirichlet_tag)
nothing # hide

# Defile computational domains.

Ω = GT.interior(mesh)
Γd = GT.boundary(mesh;physical_names=[dirichlet_tag])
nothing # hide

# Define manufactured fields

g = GT.analytical_field(sum,Ω)
f = GT.analytical_field(x->0,Ω)
nothing # hide

# Define the interpolation space.

k = 1
V = GT.lagrange_space(Ω,k;dirichlet_boundary=Γd)
nothing # hide

# Interpolate Dirichlet values.

uhd = GT.dirichlet_field(Float64,V)
GT.interpolate_dirichlet!(g,uhd)
nothing # hide

# Visualize the Dirichlet field.

Makie.plot(Ω,color=uhd,strokecolor=:blue)
FileIO.save(joinpath(@__DIR__,"fig_poisson_2.png"),Makie.current_figure()) # hide

# ![](fig_poisson_2.png)

# Define numerical integration.

degree = 2*k
dΩ = GT.measure(Ω,degree)
nothing # hide

# Define weak form.

∇ = ForwardDiff.gradient
a = (u,v) -> GT.∫( x->∇(u,x)⋅∇(v,x), dΩ)
l = v -> GT.∫( x->v(x)*f(x), dΩ)
nothing # hide

# Assemble the problem use the automatic assembly loop generator

p = GT.linear_problem(uhd,a,l)
nothing # hide

# Solve the problem

s = PS.LinearAlgebra_lu(p)
s = PS.solve(s)
nothing # hide

# Build the FE solution.

uh = GT.solution_field(uhd,s)
nothing # hide

# Visualize the solution.

Makie.plot(Ω;color=uh,strokecolor=:black)
FileIO.save(joinpath(@__DIR__,"fig_poisson_3.png"),Makie.current_figure()) # hide

# ![](fig_poisson_3.png)

# Compute the L2 norm of the discretization error.

eh = x -> uh(x) - g(x)
el2 = GT.∫( x->abs2(eh(x)), dΩ) |> sum |> sqrt



# ## Hand-written assembly
#
#

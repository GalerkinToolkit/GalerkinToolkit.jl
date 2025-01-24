# # Poisson equation
#

# Load dependencies form Julia stdlib.

using LinearAlgebra

# Import other dependencies

import GalerkinToolkit as GT
import PartitionedSolvers as PS
import ForwardDiff
import GLMakie as Makie
import FileIO # hide

# Read and visualize the mesh

assets_dir = normpath(joinpath(@__DIR__,"..","..","..","assets"))
msh_file = joinpath(assets_dir,"model.msh")
mesh = GT.mesh_from_gmsh(msh_file)
nothing # hide

Makie.plot(mesh,color=:pink,strokecolor=:blue)
FileIO.save(joinpath(@__DIR__,"fig_poisson_eq_0.png"),Makie.current_figure()) # hide

# ![](fig_poisson_eq_0.png)

# Define domains and visualize them

Ω = GT.interior(mesh)
Makie.plot(Ω,color=:pink)
FileIO.save(joinpath(@__DIR__,"fig_poisson_eq_1.png"),Makie.current_figure()) # hide

# ![](fig_poisson_eq_1.png)

dirichlet_names = ["sides"]
Γd = GT.boundary(mesh;physical_names=dirichlet_names)
Makie.plot(Γd,color=:pink)
FileIO.save(joinpath(@__DIR__,"fig_poisson_eq_2.png"),Makie.current_figure()) # hide

# ![](fig_poisson_eq_2.png)

neumann_names = ["circle", "triangle", "square"]
Γn = GT.boundary(mesh;physical_names=neumann_names)
Makie.plot(Γn,color=:pink)
FileIO.save(joinpath(@__DIR__,"fig_poisson_eq_3.png"),Makie.current_figure()) # hide

# ![](fig_poisson_eq_3.png)


# Define forcing data

f = GT.analytical_field(x->1.0,Ω)
g = GT.analytical_field(x->2.0,Ω)
h = GT.analytical_field(x->3.0,Ω)

# Define the interpolation space.

k = 1
V = GT.lagrange_space(Ω,k;dirichlet_boundary=Γd)
nothing # hide

# Interpolate Dirichlet values.

T = Float64
uhd = GT.dirichlet_field(T,V)
GT.interpolate_dirichlet!(g,uhd)
nothing # hide

# Visualize the Dirichlet field.

Makie.plot(Ω,color=uhd,strokecolor=:blue)
FileIO.save(joinpath(@__DIR__,"fig_poisson_eq_4.png"),Makie.current_figure()) # hide

# ![](fig_poisson_eq_4.png)

# Define numerical integration.

degree = 2*k
dΩ = GT.measure(Ω,degree)
dΓn = GT.measure(Γn,degree)
nothing # hide

# Define weak form.

const ∇ = ForwardDiff.gradient
a = (u,v) -> GT.∫( x->∇(u,x)⋅∇(v,x), dΩ)
l = v -> GT.∫( x->v(x)*f(x), dΩ) + GT.∫( x->v(x)*h(x), dΓn)
nothing # hide

# Assemble the problem using the automatic assembly loop generator

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
FileIO.save(joinpath(@__DIR__,"fig_poisson_eq_5.png"),Makie.current_figure()) # hide

# ![](fig_poisson_eq_5.png)



# # Hello, World!
#
# ![](fig_poisson_3.png)
#
# ## Problem statement
#
# In this example, we show how to solve the "Hello, world" PDE example:
# the Poisson equation on the unit hyper-cube $\Omega  =[0,1]^d$, $d\in\{2,3\}$, with Dirichlet boundary conditions.
#
# ```math
# \left\lbrace
# \begin{aligned}
# -\Delta u = f  \ &\text{in} \ \Omega,\\
# u = g \ &\text{on}\ \partial\Omega,\\
# \end{aligned}
# \right.
# ```
# with $f=0$ and $g(x)=\text{sum}(x)$. In this case, we know that the solution is $u=g$ which allows us to check that we solve the
# problem correctly, by integration an error norm.

#  To solve this PDE, we use a conventional Galerkin finite element (FE) method with conforming Lagrangian FE spaces (see, e.g., [1] for specific details on this formulation).

# ## Implementation
#
# Load dependencies form Julia stdlib.

using LinearAlgebra

# Import other dependencies

import GalerkinToolkit as GT
import PartitionedSolvers as PS
import ForwardDiff
import GLMakie as Makie
import LinearSolve
import FileIO # hide

# Generate the computational mesh.

domain = (0,1,0,1)
cells = (10,10)
mesh = GT.cartesian_mesh(domain,cells)
nothing # hide

# Visualize the mesh.

axis = (aspect = Makie.DataAspect(),)
Makie.plot(mesh;color=:pink,strokecolor=:blue,axis)
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

# Define differential operators

const ∇ = ForwardDiff.gradient
Δ(f,x) = tr(ForwardDiff.jacobian(y->∇(f,y),x))

# Define manufactured fields.

g = GT.analytical_field(sum,Ω)
f = GT.analytical_field(x->-Δ(g.definition,x),Ω)
nothing # hide

# Define the interpolation space.

k = 1
V = GT.lagrange_space(Ω,k;dirichlet_boundary=Γd)
nothing # hide

# Interpolate Dirichlet values.

T = Float64
uhd = GT.zero_dirichlet_field(T,V)
GT.interpolate_dirichlet!(g,uhd)
nothing # hide

# Visualize the Dirichlet field.

Makie.plot(Ω;color=uhd,strokecolor=:blue,axis)
FileIO.save(joinpath(@__DIR__,"fig_poisson_2.png"),Makie.current_figure()) # hide

# ![](fig_poisson_2.png)

# Define numerical integration.

degree = 2*k
dΩ = GT.measure(Ω,degree)
nothing # hide

# Define weak form.

a = (u,v) -> GT.∫( x->∇(u,x)⋅∇(v,x), dΩ)
l = v -> GT.∫( x->v(x)*f(x), dΩ)
nothing # hide

# Assemble the problem using the automatic assembly loop generator.
# We generate a `SciMLBase.LinearProblem` object that can be solved with `LinearSolve`.

p = GT.SciMLBase_LinearProblem(uhd,a,l)
sol = LinearSolve.solve(p)
nothing # hide

# Transform the solution object into the FE solution object.

uh = GT.solution_field(uhd,sol)
nothing # hide

# Visualize the solution.

Makie.plot(Ω;color=uh,strokecolor=:black,axis)
FileIO.save(joinpath(@__DIR__,"fig_poisson_2-1.png"),Makie.current_figure()) # hide

# ![](fig_poisson_2-1.png)

# We can also create a linear problem object to be solved with `PartitionedSolvers`.
# This is specially useful for the distributed case, but it also works for sequential runs.

p = GT.PartitionedSolvers_linear_problem(uhd,a,l)
nothing # hide

# Solve the problem

s = PS.LinearAlgebra_lu(p)
s = PS.solve(s)
nothing # hide

# Build the FE solution.

uh = GT.solution_field(uhd,s)
nothing # hide

# Visualize the solution.

Makie.plot(Ω;color=uh,strokecolor=:black,axis)
FileIO.save(joinpath(@__DIR__,"fig_poisson_3.png"),Makie.current_figure()) # hide

# ![](fig_poisson_3.png)

# Compute the L2 norm of the discretization error.

eh = x -> uh(x) - g(x)
el2 = GT.∫( x->abs2(eh(x)), dΩ) |> sum |> sqrt

# ## Final program

module Program 

using LinearAlgebra
import GalerkinToolkit as GT
import PartitionedSolvers as PS
import ForwardDiff
import GLMakie as Makie

function main(;domain,cells)
    mesh = GT.cartesian_mesh(domain,cells)
    dirichlet_tag = "dirichlet"
    GT.label_boundary_faces!(mesh;physical_name=dirichlet_tag)
    Ω = GT.interior(mesh)
    Γd = GT.boundary(mesh;physical_names=[dirichlet_tag])
    ∇ = ForwardDiff.gradient
    Δ(f,x) = tr(ForwardDiff.jacobian(y->∇(f,y),x))
    g = GT.analytical_field(sum,Ω)
    f = GT.analytical_field(x->-Δ(g.definition,x),Ω)
    k = 1
    V = GT.lagrange_space(Ω,k;dirichlet_boundary=Γd)
    T = Float64
    uhd = GT.zero_dirichlet_field(T,V)
    GT.interpolate_dirichlet!(g,uhd)
    degree = 2*k
    dΩ = GT.measure(Ω,degree)
    a = (u,v) -> GT.∫( x->∇(u,x)⋅∇(v,x), dΩ)
    l = v -> GT.∫( x->v(x)*f(x), dΩ)
    p = GT.PartitionedSolvers_linear_problem(uhd,a,l)
    s = PS.LinearAlgebra_lu(p)
    s = PS.solve(s)
    uh = GT.solution_field(uhd,s)
    eh = x -> uh(x) - g(x)
    el2 = GT.∫( x->abs2(eh(x)), dΩ) |> sum |> sqrt
end

end # module

# Run it for a 2d case.

Program.main(domain=(0,1,0,1),cells=(10,10))

# Run it for a 3d case.

Program.main(domain=(0,1,0,1,0,1),cells=(10,10,10))

# ## References
#
# [1] C. Johnson. *Numerical Solution of Partial Differential Equations by the Finite Element Method*. Dover Publications, 2009.

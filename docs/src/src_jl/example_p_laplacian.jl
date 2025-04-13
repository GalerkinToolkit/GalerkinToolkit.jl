# # p-Laplacian
#
# ![](fig_p_laplacian.gif)

# ## Problem statement

# Find the scalar-field $u$ such that
#
# ```math
# \left\lbrace
# \begin{aligned}
# -\nabla \cdot \left( |\nabla u|^{p-2} \ \nabla u \right) = f\ &\text{in}\ \Omega,\\
# u = -1 \ &\text{on} \ \Gamma_0,\\
# u = 1 \ &\text{on} \ \Gamma_1,\\
# \left( |\nabla u|^{p-2}\ \nabla u \right)\cdot n = 0 \ &\text{elsewhere on} \ \partial\Omega,
# \end{aligned}
# \right.
# ```
# with $p>2$.
# The vector field $n$ is the outwards unit normal vector to $\partial\Omega$. The computational domains are defined in the mesh file `model.msh`. The domain $\Omega$ is represented by the 3D faces in this mesh. The domain $\Gamma_0$ is represented by the physical group named `"sides"` and $\Gamma_1$ is the union of the physical groups named `"circle"`, `"triangle"`, and `"square"`.  
#  To solve this PDE, we use a conventional Galerkin finite element method with conforming Lagrangian FE spaces.

# ## Implementation


# Load dependencies form Julia stdlib.

using LinearAlgebra

# Import other dependencies.

import GalerkinToolkit as GT
import PartitionedSolvers as PS
import NonlinearSolve
import ForwardDiff
import GLMakie as Makie
import FileIO # hide

# Read and visualize the mesh.

assets_dir = normpath(joinpath(@__DIR__,"..","..","..","assets"))
msh_file = joinpath(assets_dir,"model.msh")
mesh = GT.mesh_from_msh(msh_file)
nothing # hide

Makie.plot(mesh,color=:pink,strokecolor=:blue)
FileIO.save(joinpath(@__DIR__,"fig_p_laplacian_1.png"),Makie.current_figure()) # hide

# ![](fig_p_laplacian_1.png)

# Define domains.

dirichlet_0_names = ["sides"]
dirichlet_1_names = ["circle", "triangle", "square"]
Ω = GT.interior(mesh)
Γ0 = GT.boundary(mesh;physical_names=dirichlet_0_names)
Γ1 = GT.boundary(mesh;physical_names=dirichlet_1_names)
Γd = GT.piecewise_domain(Γ0,Γ1)
nothing # hide

# Define forcing data.

g0 = GT.analytical_field(x->-1.0,Ω)
g1 = GT.analytical_field(x->1.0,Ω)
g = GT.piecewise_field(g0,g1)
nothing # hide

# Define the interpolation space.

k = 1
V = GT.lagrange_space(Ω,k;dirichlet_boundary=Γd)
nothing # hide

# Interpolate Dirichlet values.

T = Float64
uh = GT.rand_field(T,V)
GT.interpolate_dirichlet!(g,uh)
nothing # hide

# Visualize the Dirichlet field.

Makie.plot(Ω,color=uh,strokecolor=:blue)
FileIO.save(joinpath(@__DIR__,"fig_p_laplacian_2.png"),Makie.current_figure()) # hide

# ![](fig_p_laplacian_2.png)

# Define numerical integration.

degree = 2*k
dΩ = GT.measure(Ω,degree)
nothing # hide

# Define weak form.

const ∇ = ForwardDiff.gradient
const q = 3
flux(∇u) = norm(∇u)^(q-2) * ∇u
dflux(∇du,∇u) = (q-2)*norm(∇u)^(q-4)*(∇u⋅∇du)*∇u+norm(∇u)^(q-2)*∇du
res = u -> v -> GT.∫( x-> ∇(v,x)⋅GT.call(flux,∇(u,x)), dΩ)
jac = u -> (du,v) -> GT.∫( x-> ∇(v,x)⋅GT.call(dflux,∇(du,x),∇(u,x)) , dΩ)
nothing # hide

# Define non-linear problem using the automatic assembly loop generator.
# We can define a problem object from `SciMLBase` that can be solved with `NonlinearSolve.jl`.

p = GT.SciMLBase_NonlinearProblem(uh,res,jac)
nothing # hide

# Solve it with `NonlinearSolve.jl`.

sol = NonlinearSolve.solve(p;show_trace=Val(true))
@assert sol.retcode == NonlinearSolve.ReturnCode.Success
nothing # hide

# Get the FE solution object

uh = GT.solution_field(uh,sol)
nothing # hide

# Visualize the solution

Makie.plot(Ω;color=uh,strokecolor=:black)
FileIO.save(joinpath(@__DIR__,"fig_p_laplacian_2-1.png"),Makie.current_figure()) # hide

# ![](fig_p_laplacian_2-1.png)

# We can also create a nonlinear problem object to be solved with `PartitionedSolvers`.
# This is specially useful for the distributed case, but it also works for sequential runs.

uh = GT.rand_field(Float64,V)
p = GT.PartitionedSolvers_nonlinear_problem(uh,res,jac)
nothing # hide

# Define a nonlinear solver with and solve the problem `PartitionedSolvers.jl`

s = PS.newton_raphson(p,verbose=true)
s = PS.solve(s)
uh = GT.solution_field(uh,s)

# Visualize the solution

Makie.plot(Ω;color=uh,strokecolor=:black)
FileIO.save(joinpath(@__DIR__,"fig_p_laplacian_3.png"),Makie.current_figure()) # hide

# ![](fig_p_laplacian_3.png)

# Now, solve while showing the intermediate results in the iteration process.

uh = GT.rand_field(Float64,V)
GT.interpolate_dirichlet!(g,uh)
p = GT.PartitionedSolvers_nonlinear_problem(uh,res,jac)
s = PS.newton_raphson(p)
color = Makie.Observable(uh)
fig = Makie.plot(Ω;color,strokecolor=:black)
fn = joinpath(@__DIR__,"fig_p_laplacian.gif")
Makie.record(fig,fn,PS.history(s);framerate=2) do s
    color[] = GT.solution_field(uh,s)
end
nothing # hide

# ![](fig_p_laplacian.gif)
#

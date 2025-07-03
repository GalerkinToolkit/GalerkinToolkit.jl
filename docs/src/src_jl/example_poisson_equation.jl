# # Poisson equation
#
# ![](fig_poisson_1.png)
#
# ## Problem statement

#  We want to solve a simple Poisson equation, but this time on a more complex domain and also including Neumann boundary conditions.
#  Dirichlet boundary conditions are applied on $\Gamma_{\rm D}$, non-homogeneous Neumann conditions are applied on $\Gamma_{\rm N}$,   and homogeneous Neumann boundary conditions are applied on the remaining portion of the boundary. The computational domains are defined in the mesh file `model.msh`. The domain $\Omega$ is represented by the 3D faces in this mesh. The domain $\Gamma_{\rm D}$ is represented by the physical group named `"sides"` and $\Gamma_{\rm N}$ is the union of the physical groups named `"circle"`, `"triangle"`, and `"square"`.
#
#  Formally, the problem to solve is: find the scalar field $u$ such that
#
# ```math
# \left\lbrace
# \begin{aligned}
# -\Delta u = f  \ &\text{in} \ \Omega,\\
# u = g \ &\text{on}\ \Gamma_{\rm D},\\
# \nabla u\cdot n = h \ &\text{on}\  \Gamma_{\rm N},\\
# \nabla u\cdot n = 0 \ &\text{elsewhere on}\  \partial \Omega,\\
# \end{aligned}
# \right.
# ```
#  being $n$ the outwards unit normal vector to $\partial\Omega$. In this example, we chose $f(x) = 1$, $g(x) = 2$, and $h(x)=3$. The variable $x$ is the position vector $x=(x_1,x_2,x_3)$.
#
# ## Numerical scheme
#
#  We use a conventional Galerkin finite element (FE) method with conforming Lagrangian FE spaces (see, e.g., [Johnson2009](@cite)).  The weak form equation we solve consists in finding $u_h\in V_g$ such that $a(u_h,v) = \ell(v)$ for all $v\in V_0$. To this end we build a space $V$ spanned by continuous and piece-wise Lagrangian basis functions. The auxiliary spaces $V_g$ and $V_0$ are the subsets of $V$ that fulfill the Dirichlet boundary condition $g$ and $0$ on $\partial\Omega$ respectively. The bilinear and linear forms are
# ```math
#   a(u,v) \doteq \int_{\Omega} \nabla v \cdot \nabla u \ {\rm d}\Omega, \quad b(v) \doteq \int_{\Omega} v\ f  \ {\rm  d}\Omega + \int_{\Gamma_{\rm N}} v\ h \ {\rm d}\Gamma_{\rm N}.
# ```
#
# This equation results in a system of linear algebraic equations that is solved using an external linear solver from [`LinearSolve.jl`](https://github.com/SciML/LinearSolve.jl).

# ## Implementation

import FileIO # hide
using LinearAlgebra
import GalerkinToolkit as GT
import LinearSolve
import ForwardDiff
import GLMakie as Makie

#Read the mesh file
assets_dir = normpath(joinpath(@__DIR__,"..","..","..","assets"))
msh_file = joinpath(assets_dir,"model.msh")
mesh = GT.mesh_from_msh(msh_file)

#Geometry
dirichlet_names = ["sides"]
neumann_names = ["circle", "triangle", "square"]
Ω = GT.interior(mesh)
Γd = GT.boundary(mesh;physical_names=dirichlet_names)
Γn = GT.boundary(mesh;physical_names=neumann_names)

#Define forcing data
f = GT.analytical_field(x->1.0,Ω)
g = GT.analytical_field(x->2.0,Ω)
h = GT.analytical_field(x->3.0,Ω)

#Define the interpolation space.
k = 1
V = GT.lagrange_space(Ω,k;dirichlet_boundary=Γd)

#Interpolate Dirichlet values.
T = Float64
uhd = GT.zero_dirichlet_field(T,V)
GT.interpolate_dirichlet!(g,uhd)

#Define numerical integration.
degree = 2*k
dΩ = GT.measure(Ω,degree)
dΓn = GT.measure(Γn,degree)
nothing # hide

#Define weak form.
∇ = ForwardDiff.gradient
a = (u,v) -> GT.∫( x->∇(u,x)⋅∇(v,x), dΩ)
l = v -> GT.∫( x->v(x)*f(x), dΩ) + GT.∫( x->v(x)*h(x), dΓn)

#Assemble the problem and solve it
p = GT.SciMLBase_LinearProblem(uhd,a,l)
sol = LinearSolve.solve(p)

#Build the FE solution.
uh = GT.solution_field(uhd,sol)

#Visualize the solution.
fig = Makie.Figure()
elevation = 0.24π
azimuth = -0.55π
aspect = :data
ax = Makie.Axis3(fig[1,1];aspect,elevation,azimuth)
Makie.hidespines!(ax)
Makie.hidedecorations!(ax)
GT.makie_surface!(Ω;color=uh)
FileIO.save(joinpath(@__DIR__,"fig_poisson_1.png"),Makie.current_figure()) # hide
nothing # hide



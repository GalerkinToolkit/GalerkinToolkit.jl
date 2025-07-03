# # Hello, World!
#
# ![](fig_hello_world_1.png)
#
#
# ## Problem statement
#
# In this example, we show how to solve the "Hello, world" PDE example:
# the Laplace equation on the unit hyper-cube $\Omega  =[0,1]^d$, $d=3$, with Dirichlet boundary conditions.
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

#  ## Numerical scheme
#  We use a conventional Galerkin finite element (FE) method with conforming Lagrangian FE spaces (see, e.g., [Johnson2009](@cite)).  The weak form equation we solve consists in finding $u_h\in V_g$ such that $a(u_h,v) = \ell(v)$ for all $v\in V_0$. To this end we build a space $V$ spanned by continuous and piece-wise Lagrangian basis functions. The auxiliary spaces $V_g$ and $V_0$ are the subsets of $V$ that fulfill the Dirichlet boundary condition $g$ and $0$ on $\partial\Omega$ respectively. The bilinear and linear forms are
# ```math
#   a(u,v) \doteq \int_{\Omega} \nabla v \cdot \nabla u \ {\rm d}\Omega, \quad b(v) \doteq \int_{\Omega} v\ f  \ {\rm  d}\Omega.
# ```
#
# This equation results in a system of linear algebraic equations that is solved using an external linear solver from [`LinearSolve.jl`](https://github.com/SciML/LinearSolve.jl).

# ## Implementation


import FileIO # hide
using LinearAlgebra
import GalerkinToolkit as GT
import PartitionedSolvers as PS
import ForwardDiff
import GLMakie as Makie
import LinearSolve
domain = (0,1,0,1,0,1)
cells = (10,10,10)
mesh = GT.cartesian_mesh(domain,cells)
Ω = GT.interior(mesh)
Γd = GT.boundary(mesh)
∇ = ForwardDiff.gradient
g = GT.analytical_field(sum,Ω)
f = GT.analytical_field(x->0,Ω)
k = 1
V = GT.lagrange_space(Ω,k;dirichlet_boundary=Γd)
T = Float64
uhd = GT.zero_dirichlet_field(T,V)
GT.interpolate_dirichlet!(g,uhd)
degree = 2*k
dΩ = GT.measure(Ω,degree)
a = (u,v) -> GT.∫( x->∇(u,x)⋅∇(v,x), dΩ)
l = v -> GT.∫( x->v(x)*f(x), dΩ)
p = GT.SciMLBase_LinearProblem(uhd,a,l)
sol = LinearSolve.solve(p)
uh = GT.solution_field(uhd,sol)
eh = x -> uh(x) - g(x)
el2 = GT.∫( x->abs2(eh(x)), dΩ) |> sum |> sqrt
@assert el2 < 1.0e-9
fig = Makie.Figure()
ax = Makie.Axis3(fig[1,1],aspect=:data)
Makie.hidespines!(ax)
Makie.hidedecorations!(ax)
GT.makie_surface!(Ω;color=uh)
FileIO.save(joinpath(@__DIR__,"fig_hello_world_1.png"),Makie.current_figure()) # hide
nothing # hide



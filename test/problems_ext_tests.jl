module ProblemsExt

import GalerkinToolkit as GT
import LinearSolve
import NonlinearSolve
import ForwardDiff
using Test
using LinearAlgebra

# Linear problem

n = 10
domain = (0,1,0,1)
cells = (n,n)
mesh = GT.cartesian_mesh(domain,cells)
D = GT.num_dims(mesh)
Ω = GT.interior(mesh)
Γ = GT.boundary(mesh)
order = 1
degree = 2*order
V = GT.lagrange_space(Ω,order;dirichlet_boundary=Γ)
dΩ = GT.measure(Ω,degree)
u = GT.analytical_field(sum,Ω)
uhd = GT.interpolate_dirichlet(u,V)
∇ = (u,q) -> ForwardDiff.gradient(u,q)
a = (u,v) -> GT.∫( q->∇(u,q)⋅∇(v,q), dΩ)
l = v -> 0

prob =  GT.SciMLBase_LinearProblem(uhd,a,l)
sol = LinearSolve.solve(prob)
uh = GT.solution_field(uhd,sol)

eh = q -> u(q) - uh(q)
∇eh =q -> ∇(u,q) - ∇(uh,q)
tol = 1.0e-10
el2 = GT.∫( q->abs2(eh(q)), dΩ) |> sum |> sqrt
eh1 = GT.∫( q->∇eh(q)⋅∇eh(q), dΩ) |> sum |> sqrt
@test el2 < tol
@test eh1 < tol

# Nonlinear problem

n = 10
domain = (0,1,0,1)
cells = (n,n)
mesh = GT.cartesian_mesh(domain,cells)
D = GT.num_dims(mesh)
Ω = GT.interior(mesh)
order = 1
degree = 2*order
V = GT.lagrange_space(Ω,order;dirichlet_boundary=Γ)
dΩ = GT.measure(Ω,degree)
uh = GT.rand_field(Float64,V)

∇ = (u,q) -> ForwardDiff.gradient(u,q)
const q = 3
flux(∇u) = norm(∇u)^(q-2) * ∇u
dflux(∇du,∇u) = (q-2)*norm(∇u)^(q-4)*(∇u⋅∇du)*∇u+norm(∇u)^(q-2)*∇du
res = u -> v -> GT.∫( x-> ∇(v,x)⋅GT.call(flux,∇(u,x)) - v(x), dΩ)
jac = u -> (du,v) -> GT.∫( x-> ∇(v,x)⋅GT.call(dflux,∇(du,x),∇(u,x)) , dΩ)

prob = GT.SciMLBase_NonlinearProblem(uh,res,jac)
sol = NonlinearSolve.solve(prob;show_trace=Val(true))
@test sol.retcode == NonlinearSolve.ReturnCode.Success
uh = GT.solution_field(uh,sol)
uhl2 = GT.∫( q->abs2(uh(q)), dΩ) |> sum |> sqrt
@test abs(uhl2-0.09133166701839236)<tol


end # module

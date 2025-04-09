module ProblemsExt

import GalerkinToolkit as GT
import LinearSolve
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
Γ = GT.boundary(mesh)
order = 1
degree = 2*order
V = GT.lagrange_space(Ω,order;dirichlet_boundary=Γ)
dΩ = GT.measure(Ω,degree)
u = GT.analytical_field(sum,Ω)
uhd = GT.interpolate_dirichlet(u,V)


end # module

# # Discretization methods

import GalerkinToolkit as GT
import PartitionedSolvers as PS
import GLMakie as Makie
import ForwardDiff
using LinearAlgebra
using Test
import FileIO # hide

domain = (0,1,0,1)
cells = (10,10)
mesh = GT.cartesian_mesh(domain,cells)
D = GT.num_dims(mesh)
GT.label_boundary_faces!(mesh;physical_name="Γd")
Ω = GT.interior(mesh)
Γd = GT.boundary(mesh;physical_names=["Γd"])
n_Γd = GT.unit_normal(mesh,D-1)
u = GT.analytical_field(x->sum(x),Ω)
#TODO these two definitions should be enough
∇ = ForwardDiff.gradient
Δ = (u,x) -> tr(ForwardDiff.jacobian(y->∇(u,y),x))
#gradient(u) = x->ForwardDiff.gradient(u,x)
#jacobian(u) = x->ForwardDiff.jacobian(u,x)
#laplacian(u) = x-> tr(ForwardDiff.jacobian(y->ForwardDiff.gradient(u,y),x))
#Δ(u) = GT.call(laplacian,u)
#∇(u) = GT.call(gradient,u)
#Δ(u,x) = Δ(u)(x)
#∇(u,x) = ∇(u)(x)
tol = 1e-10
order = 2
γ = order*(order+1)/10

# ## Continuous Galerkin

V = GT.lagrange_space(Ω,order;dirichlet_boundary=Γd)
uhd = GT.dirichlet_field(Float64,V)
GT.interpolate_dirichlet!(u,uhd)
dΩ = GT.measure(Ω,2*order)
a_Ω = (u,v) -> GT.∫( x->∇(u,x)⋅∇(v,x), dΩ)
l_Ω = (v) -> GT.∫( x-> -v(x)*GT.call(Δ,u,x), dΩ)
p = GT.linear_problem(uhd,a_Ω,l_Ω)
s = PS.solve(p)
uh = GT.solution_field(uhd,s)
eh = x -> u(x) - uh(x)
el2 = sqrt(sum(GT.∫( x->abs2(eh(x)), dΩ)))
@test el2 < tol

# ## Static condensation
#
# !!! warning
#     TODO
#

# ## Nitche Method

V = GT.lagrange_space(Ω,order)
n_Γd = GT.unit_normal(mesh,D-1)
h_Γd = GT.face_diameter_field(Γd)
dΓd = GT.measure(Γd,2*order)
a_Γd = (u,v) -> GT.∫(
    x-> (γ/h_Γd(x))*v(x)*u(x)-
    v(x)*n_Γd(x)⋅∇(u,x)-
    n_Γd(x)⋅∇(v,x)*u(x), dΓd)
l_Γd = v -> GT.∫(
    x-> (γ/h_Γd(x))*v(x)*u(x)-
    n_Γd(x)⋅∇(v,x)*u(x), dΓd)
a_ni = (u,v) -> a_Ω(u,v) + a_Γd(u,v)
l_ni = v -> l_Ω(v) + l_Γd(v)
p = GT.linear_problem(Float64,V,a_ni,l_ni)
s = PS.solve(p)
uh = GT.solution_field(V,s)
eh = x -> u(x) - uh(x)
el2 = sqrt(sum(GT.∫( x->abs2(eh(x)), dΩ)))
@test el2 < tol

# ## Lagrange multipliers

Q = GT.lagrange_space(Γd,order-1;conformity=:L2)
VxQ = V × Q
A = ((u,p),(v,q)) -> a_Ω(u,v) + GT.∫(
    x->(u(x)+p(x))*(v(x)+q(x))-
    u(x)*v(x)-p(x)*q(x), dΓd)
L = ((v,q),) -> l_Ω(v) + GT.∫(x->u(x)*q(x), dΓd)
p = GT.linear_problem(Float64,VxQ,A,L)
s = PS.solve(p)
uh, = GT.solution_field(VxQ,s)
eh = x -> u(x) - uh(x)
el2 = sqrt(sum(GT.∫( x->abs2(eh(x)), dΩ)))
@test el2 < tol

# ## Interior Penalty

mean(u) = 0.5*(u[1]+u[2])
jump(u,n) = u[2]*n[2] + u[1]*n[1]
GT.label_interior_faces!(mesh;physical_name="Λ")
Λ = GT.skeleton(mesh;physical_names=["Λ"])
dΛ = GT.measure(Λ,2*order)
n_Λ = GT.unit_normal(mesh,D-1)
h_Λ = GT.face_diameter_field(Λ)
V = GT.lagrange_space(Ω,order;conformity=:L2)
a_Λ = (u,v) -> GT.∫(
                    x-> (γ/h_Λ(x))*jump(v(x),n_Λ(x))⋅jump(u(x),n_Λ(x))-
                    jump(v(x),n_Λ(x))⋅mean(∇(u,x))-
                    mean(∇(v,x))⋅jump(u(x),n_Λ(x)), dΛ)
a_ip = (u,v) -> a_Ω(u,v) + a_Γd(u,v) + a_Λ(u,v)
l_ip = l_ni
p = GT.linear_problem(Float64,V,a_ip,l_ip)
s = PS.solve(p)
uh = GT.solution_field(V,s)
eh = x -> u(x) - uh(x)
el2 = sqrt(sum(GT.∫( x->abs2(eh(x)), dΩ)))
@test el2 < tol

# We can also combine the interior penalty the impose inter-element continuity with the auxiliary Lagrange multiplier
# field to impose Dirichlet boundary conditions.

VxQ = V × Q
A_ip = ((u,p),(v,q)) -> A((u,p),(v,q)) + a_Λ(u,v)
L_ip = L
p = GT.linear_problem(Float64,VxQ,A_ip,L_ip)
s = PS.solve(p)
uh, = GT.solution_field(VxQ,s)
eh = x -> u(x) - uh(x)
@test el2 < tol

# ## Hybrid method
#
# !!! warning
#     TODO. Raviart-Thomas interpolation needed for this case
#

# ## HDG
#
# !!! warning
#     TODO. This one is challenging
#

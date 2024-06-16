module InterpolationTests

import GalerkinToolkit as gk
using GalerkinToolkit: ∫, ×
using Test
import ForwardDiff
using LinearAlgebra

D = 2
order = 3
fe = gk.lagrangian_fe(gk.unit_n_cube(Val(D)),order)

gk.face_nodes(fe,0)
gk.face_nodes(fe,1)
gk.face_nodes(fe,2)

gk.face_interior_nodes(fe,0)
gk.face_interior_nodes(fe,1)
gk.face_interior_nodes(fe,2)

gk.face_interior_node_permutations(fe,0)
gk.face_interior_node_permutations(fe,1)
gk.face_interior_node_permutations(fe,2)

gk.face_dofs(fe,0)
gk.face_dofs(fe,1)
gk.face_dofs(fe,2)

gk.face_own_dofs(fe,0)
gk.face_own_dofs(fe,1)
gk.face_own_dofs(fe,2)

gk.face_own_dof_permutations(fe,0)
gk.face_own_dof_permutations(fe,1)
gk.face_own_dof_permutations(fe,2)

domain = (0,1,0,1)
cells = (3,3)
mesh = gk.cartesian_mesh(domain,cells)
gk.label_boundary_faces!(mesh;physical_name="boundary_faces")

topo = gk.topology(mesh)

gk.face_permutation_ids(topo,2,0)
gk.face_permutation_ids(topo,2,1)
gk.face_permutation_ids(topo,2,2)

Ω = gk.domain(mesh)
Ωref = gk.domain(mesh;is_reference_domain=true)
ϕ = gk.domain_map(Ωref,Ω)

D = gk.num_dims(mesh)
Γdiri = gk.domain(mesh;face_dim=D-1,physical_names=["boundary_faces"])

V = gk.iso_parametric_space(Ωref;dirichlet_boundary=Γdiri)

v = gk.zero_field(Float64,V)
v2 = gk.zero_field(Float64,V)

u = gk.analytical_field(x->sum(x),Ω)
uref = u∘ϕ

gk.interpolate!(uref,v)
gk.interpolate_dirichlet!(uref,v2)

Y = V×V
y = gk.zero_field(Float64,Y)
y1, y2 = y

order = 3
V = gk.lagrange_space(Ωref,order)
gk.face_dofs(V)

V = gk.lagrange_space(Ωref,order;dirichlet_boundary=Γdiri)
gk.face_dofs(V)

w = gk.zero_field(Float64,V)
w2 = gk.zero_field(Float64,V)
w3 = gk.zero_field(Float64,V)

gk.interpolate!(uref,w)
gk.interpolate_dirichlet!(uref,w2)
gk.interpolate_dirichlet!(uref,w3)

x = gk.free_values(w3)
x .= rand(length(x))

uh = gk.zero_field(Float64,V)
gk.interpolate!(uref,uh)
eh(q) = u(ϕ(q)) - uh(q)
tol = 1.e-12
degree = 2
dΩref = gk.measure(Ωref,degree)
function dV(q)
    J = ForwardDiff.jacobian(ϕ,q)
    abs(det(J))
end
el2 = ∫( q->abs2(eh(q))*dV(q), dΩref) |> sum |> sqrt
@test el2 < tol

# TODO
#udiri = gk.analytical_field(sum,Γdiri)
#uh = gk.zero_field(Float64,V)
#gk.interpolate_free!(uref,uh)
#gk.interpolate_dirichlet!(udiri,uh)
#eh(q) = u(ϕ(q)) - uh(q)
#el2 = ∫( q->abs2(eh(q))*dV(q), dΩref) |> sum |> sqrt
#@test el2 < tol

V = gk.lagrange_space(Γdiri,order)
gk.face_dofs(V)

Γ1 = gk.domain(mesh;face_dim=D-1,physical_names=["1-face-1"])
Γ2 = gk.domain(mesh;face_dim=D-1,physical_names=["1-face-3"])
Γ3 = gk.domain(mesh;face_dim=D-2,physical_names=["0-face-1"])

u1 = gk.analytical_field(x->1.0,Ωref)
u2 = gk.analytical_field(x->2.0,Ωref)
u3 = gk.analytical_field(x->3.0,Ωref)

# TODO better names than piecewiese_field and piecewiese_domain?
udiri = gk.piecewiese_field(u1,u2,u3)
Γdiri = gk.piecewiese_domain(Γ1,Γ2,Γ3)

order = 1
V = gk.lagrange_space(Ωref,order;dirichlet_boundary=Γdiri)

gk.face_dofs(V)

V |> gk.dirichlet_dof_location

uh2 = gk.zero_field(Float64,V)
gk.interpolate_dirichlet!(udiri,uh2)

order = 1
V = gk.lagrange_space(Ωref,order;dirichlet_boundary=Γdiri,conformity=:L2)

outdir = mkpath(joinpath(@__DIR__,"..","output"))
gk.vtk_plot(joinpath(outdir,"omega_ref"),Ωref;refinement=40) do plt
    gk.plot!(plt,v;label="v")
    gk.plot!(plt,v2;label="v2")
    gk.plot!(plt,y1;label="y1")
    gk.plot!(plt,y2;label="y2")
    gk.plot!(plt,w;label="w")
    gk.plot!(plt,w2;label="w2")
    gk.plot!(plt,w3;label="w3")
    gk.plot!(plt,uh2;label="uh2")
end

end #module

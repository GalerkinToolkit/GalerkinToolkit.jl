module InterpolationTests

import GalerkinToolkit as GT
using GalerkinToolkit: ∫, ×
using Test
import ForwardDiff
using LinearAlgebra
using BlockArrays
using WriteVTK

D = 2
order = 3
fe = GT.lagrangian_fe(GT.unit_n_cube(Val(D)),order)

GT.face_nodes(fe,0)
GT.face_nodes(fe,1)
GT.face_nodes(fe,2)

GT.face_interior_nodes(fe,0)
GT.face_interior_nodes(fe,1)
GT.face_interior_nodes(fe,2)

GT.face_interior_node_permutations(fe,0)
GT.face_interior_node_permutations(fe,1)
GT.face_interior_node_permutations(fe,2)

GT.face_dofs(fe,0)
GT.face_dofs(fe,1)
GT.face_dofs(fe,2)

GT.face_own_dofs(fe,0)
GT.face_own_dofs(fe,1)
GT.face_own_dofs(fe,2)

GT.face_own_dof_permutations(fe,0)
GT.face_own_dof_permutations(fe,1)
GT.face_own_dof_permutations(fe,2)

domain = (0,1,0,1)
cells = (3,3)
mesh = GT.cartesian_mesh(domain,cells)
GT.label_boundary_faces!(mesh;physical_name="boundary_faces")

topo = GT.topology(mesh)

GT.face_permutation_ids(topo,2,0)
GT.face_permutation_ids(topo,2,1)
GT.face_permutation_ids(topo,2,2)

Ω = GT.interior(mesh)
Ωref = GT.interior(mesh;is_reference_domain=true)
ϕ = GT.domain_map(Ωref,Ω)

D = GT.num_dims(mesh)
Γdiri = GT.boundary(mesh;physical_names=["boundary_faces"])

V = GT.iso_parametric_space(Ωref;dirichlet_boundary=Γdiri)

v = GT.zero_field(Float64,V)
v2 = GT.zero_field(Float64,V)

u = GT.analytical_field(x->sum(x),Ω)
uref = u∘ϕ

GT.interpolate!(uref,v)
GT.interpolate_dirichlet!(uref,v2)

Y = V×V
y = GT.zero_field(Float64,Y)
display(GT.free_values(y))
display(GT.free_dofs(Y))
y1, y2 = y

@test GT.free_values(y1) === blocks(GT.free_values(y))[1]
@test GT.free_values(y2) === blocks(GT.free_values(y))[2]

GT.interpolate!(uref,y1)
GT.interpolate!(uref,y,1)

GT.interpolate_free!(uref,y1)
GT.interpolate_free!(uref,y,1)

GT.interpolate_dirichlet!(uref,y1)
GT.interpolate_dirichlet!(uref,y,1)

order = 3
V = GT.lagrange_space(Ωref,order)
GT.face_dofs(V)

V = GT.lagrange_space(Ωref,order;dirichlet_boundary=Γdiri)
GT.face_dofs(V)

w = GT.zero_field(Float64,V)
w2 = GT.zero_field(Float64,V)
w3 = GT.zero_field(Float64,V)

GT.interpolate!(uref,w)
GT.interpolate_dirichlet!(uref,w2)
GT.interpolate_dirichlet!(uref,w3)

x = GT.free_values(w3)
x .= rand(length(x))

uh = GT.zero_field(Float64,V)
GT.interpolate!(uref,uh)
eh(q) = u(ϕ(q)) - uh(q)
tol = 1.e-12
degree = 2
dΩref = GT.measure(Ωref,degree)
function dV(q)
    J = ForwardDiff.jacobian(ϕ,q)
    abs(det(J))
end
el2 = ∫( q->abs2(eh(q))*dV(q), dΩref) |> sum |> sqrt
@test el2 < tol

uhd = GT.dirichlet_field(Float64,V)
GT.interpolate_dirichlet!(uref,uh)

# TODO
#udiri = GT.analytical_field(sum,Γdiri)
#uh = GT.zero_field(Float64,V)
#GT.interpolate_free!(uref,uh)
#GT.interpolate_dirichlet!(udiri,uh)
#eh(q) = u(ϕ(q)) - uh(q)
#el2 = ∫( q->abs2(eh(q))*dV(q), dΩref) |> sum |> sqrt
#@test el2 < tol

V = GT.lagrange_space(Γdiri,order)
GT.face_dofs(V)

Γ1 = GT.boundary(mesh;physical_names=["1-face-1"])
Γ2 = GT.boundary(mesh;physical_names=["1-face-3"])
Γ3 = GT.boundary(mesh;physical_names=["0-face-1"])

u1 = GT.analytical_field(x->1.0,Ωref)
u2 = GT.analytical_field(x->2.0,Ωref)
u3 = GT.analytical_field(x->3.0,Ωref)

# TODO better names than piecewiese_field and piecewiese_domain?
udiri = GT.piecewiese_field(u1,u2,u3)
Γdiri = GT.piecewiese_domain(Γ1,Γ2,Γ3)

order = 1
V = GT.lagrange_space(Ωref,order;dirichlet_boundary=Γdiri)

GT.face_dofs(V)

V |> GT.dirichlet_dof_location

uh2 = GT.zero_field(Float64,V)
GT.interpolate_dirichlet!(udiri,uh2)

order = 1
V = GT.lagrange_space(Ωref,order;dirichlet_boundary=Γdiri,conformity=:L2)
GT.face_dofs(V)
V = GT.lagrange_space(Ωref,order-1;dirichlet_boundary=Γdiri,conformity=:L2)
GT.face_dofs(V)

outdir = mkpath(joinpath(@__DIR__,"..","output"))
vtk_grid(joinpath(outdir,"omega_ref"),Ωref;plot_params=(;refinement=40)) do plt
    GT.plot!(plt,v;label="v")
    GT.plot!(plt,v2;label="v2")
    GT.plot!(plt,y1;label="y1")
    GT.plot!(plt,y2;label="y2")
    GT.plot!(plt,w;label="w")
    GT.plot!(plt,w2;label="w2")
    GT.plot!(plt,w3;label="w3")
    GT.plot!(plt,uh2;label="uh2")
end

Γ = GT.boundary(mesh;physical_names=["boundary_faces"])

V = GT.lagrange_space(Γ,order;conformity=:L2)
GT.face_dofs(V)
V = GT.lagrange_space(Γ,order-1;conformity=:L2)
GT.reference_fes(V) # TODO why 2 reference fes?
GT.face_dofs(V)

end #module

module InterpolationTests

import GalerkinToolkit as gk
using GalerkinToolkit: ×
using Test

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

outdir = mkpath(joinpath(@__DIR__,"..","output"))

domain = (0,1,0,1)
cells = (4,4)
mesh = gk.cartesian_mesh(domain,cells)

topo = gk.topology(mesh)

gk.face_permutation_ids(topo,2,0)
gk.face_permutation_ids(topo,2,1)
gk.face_permutation_ids(topo,2,2)

Ω = gk.domain(mesh)
Ωref = gk.domain(mesh;is_reference_domain=true)
ϕ = gk.domain_map(Ωref,Ω)

D = gk.num_dims(mesh)
Γdiri = gk.domain(mesh;face_dim=D-1)

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

gk.interpolate!(uref,w)
gk.interpolate_dirichlet!(uref,w2)

gk.vtk_plot(joinpath(outdir,"omega_ref"),Ωref;refinement=4) do plt
    gk.plot!(plt,v;label="v")
    gk.plot!(plt,v2;label="v2")
    gk.plot!(plt,y1;label="y1")
    gk.plot!(plt,y2;label="y2")
    gk.plot!(plt,v;label="w")
    gk.plot!(plt,v2;label="w2")
end

V = gk.lagrange_space(Γdiri,order)
gk.face_dofs(V)

end #module

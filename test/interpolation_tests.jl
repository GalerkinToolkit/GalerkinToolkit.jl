module InterpolationTests

import GalerkinToolkit as GT
using GalerkinToolkit: ∫, ×
using Test
import ForwardDiff
using LinearAlgebra
using BlockArrays
using WriteVTK
using StaticArrays

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

for major in (:component,:node)

    fe = GT.lagrangian_fe(GT.unit_n_cube(Val(D)),order;shape=(2,),major)

    GT.face_nodes(fe,0)
    GT.face_nodes(fe,1)
    GT.face_nodes(fe,2)

    GT.face_dofs(fe,0)
    GT.face_dofs(fe,1)
    GT.face_dofs(fe,2)

    GT.face_own_dofs(fe,0)
    GT.face_own_dofs(fe,1)
    GT.face_own_dofs(fe,2)

    GT.face_own_dof_permutations(fe,0)
    GT.face_own_dof_permutations(fe,1)
    GT.face_own_dof_permutations(fe,2)

end

fe = GT.lagrangian_fe(GT.unit_n_cube(Val(D)),order,space=:P)

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

for geo in (GT.unit_n_cube(Val(2)),GT.unit_simplex(Val(2)))
    for order in (0,1)
        fe = GT.raviart_thomas_fe(geo,order)
        GT.face_dofs(fe,0)
        GT.face_dofs(fe,1)
        GT.face_dofs(fe,2)
        @show GT.face_own_dofs(fe,0)
        @show GT.face_own_dofs(fe,1)
        @show GT.face_own_dofs(fe,2)
        GT.face_own_dof_permutations(fe,0)
        GT.face_own_dof_permutations(fe,1)
        GT.face_own_dof_permutations(fe,2)
    end
end

domain = (0,1,0,1)
cells = (2,2)
mesh = GT.cartesian_mesh(domain,cells;simplexify=true)
GT.group_boundary_faces!(mesh;group_name="boundary_faces")
GT.group_interior_faces!(mesh;group_name="interior_faces")

order = 1

D = GT.num_dims(mesh)
Ω = GT.interior(mesh)
Λ = GT.skeleton(mesh;group_names=["interior_faces"])
dΛ = GT.measure(Λ,2*order)
n = GT.unit_normal(mesh,D-1)
V = GT.raviart_thomas_space(Ω,order)
#V = GT.lagrange_space(Ω,order,shape=(D,))
uh = GT.zero_field(Float64,V)
u = GT.analytical_field(identity,Ω)
GT.interpolate!(u,uh)
rh = GT.rand_field(Float64,V)

#plt = GT.plot(Ω,refinement=10)
#GT.plot!(plt,uh;label="uh")
#GT.plot!(plt,rh;label="rh")
#
##foreach(1:GT.num_free_dofs(V)) do i
##    sh = GT.zero_field(Float64,V)
##    free_vals = GT.free_values(sh)
##    free_vals[i] = 1
##    GT.plot!(plt,sh;label="sh-$i")
##end
#vtk_grid("rt",plt) |> close
#
#xxxx

domain = (0,1,0,1)
cells = (3,3)
mesh = GT.cartesian_mesh(domain,cells)
GT.group_boundary_faces!(mesh;group_name="boundary_faces")

topo = GT.topology(mesh)

GT.face_permutation_ids(topo,2,0)
GT.face_permutation_ids(topo,2,1)
GT.face_permutation_ids(topo,2,2)

Ω = GT.interior(mesh)
Ωref = GT.interior(mesh;is_reference_domain=true)
D = GT.num_dims(mesh)
ϕ = GT.physical_map(mesh,D)

D = GT.num_dims(mesh)
Γdiri = GT.boundary(mesh;group_names=["boundary_faces"])

V = GT.iso_parametric_space(Ωref;dirichlet_boundary=Γdiri)

@test V.cache.face_dofs == GT.face_dofs(V)

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

@test V.cache.face_dofs == GT.face_dofs(V)

V = GT.lagrange_space(Ωref,order;dirichlet_boundary=Γdiri)

@test V.cache.face_dofs == GT.face_dofs(V)

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

uhd = GT.zero_dirichlet_field(Float64,V)
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

@test V.cache.face_dofs == GT.face_dofs(V)

Γdiri = GT.boundary(mesh;group_names=["1-face-1","1-face-3","0-face-1"])
udiri = analytical_field(Γdiri;piecewise=true) do x,name
    if name === "1-face-1"
        1.0
    elseif name == "1-face-2"
        2.0
    elseif name == "1-face-2"
        3.0
    else
        0.0
    end
end

order = 1
V = GT.lagrange_space(Ωref,order;dirichlet_boundary=Γdiri)

@test V.cache.face_dofs == GT.face_dofs(V)

V |> GT.dirichlet_dof_location

uh2 = GT.zero_field(Float64,V)
GT.interpolate_dirichlet!(udiri,uh2)

order = 1
V = GT.lagrange_space(Ωref,order;dirichlet_boundary=Γdiri,continuous=false)

@test V.cache.face_dofs == GT.face_dofs(V)

V = GT.lagrange_space(Ωref,order-1;dirichlet_boundary=Γdiri,continuous=false)

@test V.cache.face_dofs == GT.face_dofs(V)

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

Γ = GT.boundary(mesh;group_names=["boundary_faces"])

V = GT.lagrange_space(Γ,order;continuous=false)

@test V.cache.face_dofs == GT.face_dofs(V)

V = GT.lagrange_space(Γ,order-1;continuous=false)
GT.reference_fes(V) # TODO why 2 reference fes?

@test V.cache.face_dofs == GT.face_dofs(V)


order = 3
m = GT.analytical_field(x->SVector(false,true),Γ)
V = GT.lagrange_space(Ω,order;shape=(2,),dirichlet_boundary=m)

@test V.cache.face_dofs == GT.face_dofs(V)

uh = GT.rand_field(Float64,V)

vtk_grid(joinpath(outdir,"Vvec"),Ω;plot_params=(;refinement=10)) do plt
    GT.plot!(plt,uh;label="uh")
    GT.plot!(plt,x->uh(x)[1];label="uh1")
    GT.plot!(plt,x->uh(x)[2];label="uh2")
end

Γ1 = GT.boundary(mesh;group_names=["1-face-1"])
Γ2 = GT.boundary(mesh;group_names=["1-face-3"])

order = 3
m1 = GT.analytical_field(x->SVector(false,true),Γ1)
m2 = GT.analytical_field(x->SVector(true,false),Γ2)
m = GT.piecewise_field(m1,m2)
V = GT.lagrange_space(Ω,order;shape=(2,),dirichlet_boundary=m)

@test V.cache.face_dofs == GT.face_dofs(V)

uh = GT.rand_field(Float64,V)
vtk_grid(joinpath(outdir,"Vvec2"),Ω;plot_params=(;refinement=10)) do plt
    GT.plot!(plt,uh;label="uh")
    GT.plot!(plt,x->uh(x)[1];label="uh1")
    GT.plot!(plt,x->uh(x)[2];label="uh2")
end

order = 0
V = GT.lagrange_space(Ω,order,space=:P,dirichlet_boundary=GT.last_dof())

@test V.cache.face_dofs == GT.face_dofs(V)

uh = GT.rand_field(Float64,V)

vtk_grid(joinpath(outdir,"Vpdisc"),Ω;plot_params=(;refinement=10)) do plt
    GT.plot!(plt,uh;label="uh")
end

end #module

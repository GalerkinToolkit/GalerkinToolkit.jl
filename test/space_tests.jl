module SpaceTests

import GalerkinToolkit as GT
using GalerkinToolkit: ∫, ×
using Test
import ForwardDiff
using LinearAlgebra
using BlockArrays
using WriteVTK
using StaticArrays
using PartitionedArrays

#using InteractiveUtils

domain = (0,1,0,1)
cells_per_dir = (4,4)
parts_per_dir = (2,2)
mesh = GT.cartesian_mesh(domain,cells_per_dir)

Ω = GT.interior(mesh)
Γ = GT.boundary(mesh)
Λ = GT.skeleton(mesh)
F = GT.domain(mesh,1)

order = 2

V = GT.lagrange_space(Ω,order)
V = GT.lagrange_space(Ω,order;dirichlet_boundary=Γ)

V = GT.lagrange_space(Γ,order)

uh = GT.rand_field(Float64,V)
plt = GT.plot(F)
GT.plot!(plt,uh,label="uh")
vtk_grid("plt",plt) |> close


V = GT.lagrange_space(Λ,order)

domain = (0,1,0,1)
cells_per_dir = (4,4)
parts_per_dir = (2,2)
np = prod(parts_per_dir)
parts = DebugArray(LinearIndices((np,)))
mesh = GT.cartesian_pmesh(domain,cells_per_dir,parts,parts_per_dir)

Ω = GT.interior(mesh)
Γ = GT.boundary(mesh)
Λ = GT.skeleton(mesh)

#vtk_grid("pgamma",Γ) |> close
#vtk_grid("plam",Λ) |> close

# TODO
#map(partition(Ω)) do dom
#    imesh = GT.mesh(dom)
#    ids = GT.face_local_indices(imesh,2)
#    local_to_owner(ids)
#end |> display
#vtk_grid("pdom",Ω) |> close

order = 1

V = GT.lagrange_space(Ω,order)
V = GT.lagrange_space(Ω,order;dirichlet_boundary=Γ)

@test GT.face_dofs(V) isa PVector

#map(partition(GT.free_dofs(V))) do ids
#    display(local_to_owner(ids))
#end

map(partition(V)) do space
    display(GT.face_dofs(space))
end

#u = GT.analytical_field(sum,Ω)
#T = Float64
#uh = GT.undef_field(T,V)
#GT.interpolate!(u,uh)


#@code_warntype GT.unit_simplex(Val(3))

spx0 = GT.unit_simplex(0)
spx1 = GT.unit_simplex(1)
spx2 = GT.unit_simplex(2)
spx3 = GT.unit_simplex(3)
@test GT.num_dims(spx3) == 3
display(spx3)

cube0 = GT.unit_n_cube(0)
cube1 = GT.unit_n_cube(1)
cube2 = GT.unit_n_cube(2)
cube3 = GT.unit_n_cube(3)
@test GT.num_dims(cube3) == 3
display(cube3)

#@code_warntype GT.options(cube0)
#@code_warntype GT.options()

options = GT.options(cube0)
@show GT.real_type(options)
@show GT.int_type(options)

#@code_warntype GT.real_type(options)

degree = 2
qua = GT.quadrature(cube2,degree)

#@code_warntype GT.quadrature(spx2,2)

@test cube2 === GT.domain(qua)
x = GT.coordinates(qua)
@test sum(GT.weights(qua)) ≈ 1

degree = 2
qua = GT.quadrature(spx2,degree)
@test spx2 === GT.domain(qua)
x = GT.coordinates(qua)
@test sum(GT.weights(qua)) ≈ 0.5

fe = GT.lagrange_space(cube2,1)
display(fe)

fe2 = GT.lagrange_space(cube2,1)

@test isequal(fe,fe2)
@test hash(fe) == hash(fe2)

x = GT.monomial_exponents(fe)

#@code_warntype GT.monomial_exponents(fe)

display(x)

x = GT.node_coordinates(fe)
#@code_warntype GT.node_coordinates(fe)
display(x)

@test GT.num_dofs(fe) == GT.num_nodes(fe)

A = GT.tabulator(fe)(GT.value,x)
n = GT.num_dofs(fe)
@test A ≈ Matrix{Float64}(I,n,n)

t = GT.tabulator(fe)
#@code_warntype t(GT.value,x)

qua = GT.node_quadrature(fe)
#@code_warntype GT.node_quadrature(fe)

fe = GT.lagrange_space(cube2,1;tensor_size=Val((2,)))

@test GT.num_dofs(fe) == 2*GT.num_nodes(fe)

A = GT.tabulator(fe)(GT.value,x)

@show GT.node_dofs(fe)
@show GT.dof_node(fe)

fe = GT.lagrange_space(cube2,0)
x = GT.node_coordinates(fe)
@test x[1] ≈ [0.5,0.5]

fe = GT.lagrange_space(spx0,1)
x = GT.monomial_exponents(fe)
display(x)

x = GT.node_coordinates(fe)
display(x)

@test GT.order(fe) == 0

A = GT.tabulator(fe)(GT.value,x)
n = GT.num_dofs(fe)
@test A ≈ Matrix{Float64}(I,n,n)

mesh = GT.simplexify(fe)
mesh = GT.complexify(fe)
mesh = GT.simplexify(cube2)
mesh = GT.complexify(cube2)

D = 0
order = 1
fe = GT.lagrange_space(GT.unit_n_cube(Val(D)),order)
@test GT.interior_nodes(fe) == [1]

D = 1
order = 1
fe = GT.lagrange_space(GT.unit_n_cube(Val(D)),order)

@test GT.num_dims(fe) == 1
@test GT.num_nodes(fe) == 2
@test GT.order(fe) == 1
@test GT.conforming(fe) == true
@test GT.space_type(fe) == :Q

@test GT.interior_nodes(fe) == Int[]

D = 2
order = 1
fe = GT.lagrange_space(GT.unit_n_cube(Val(D)),order)
@show GT.interior_nodes(fe)

@test GT.face_nodes(fe,0) == [[1], [2], [3], [4]]
@test GT.face_nodes(fe,1) == [[1, 2], [3, 4], [1, 3], [2, 4]]
@test GT.face_nodes(fe,2) == [[1, 2, 3, 4]]

@test GT.face_interior_nodes(fe,0) == [[1], [2], [3], [4]]
@test GT.face_interior_nodes(fe,1) == [[], [], [], []]
@test GT.face_interior_nodes(fe,2) == [Int64[]]


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

D = 2
order = 3
fe = GT.lagrange_space(GT.unit_n_cube(Val(D)),order)

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

    fe = GT.lagrange_space(GT.unit_n_cube(Val(D)),order;tensor_size=Val((2,)),major=Val(major))

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

fe = GT.lagrange_space(GT.unit_n_cube(Val(D)),order;space_type=:P)

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
        fe = GT.raviart_thomas_space(geo,order)
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
mesh = GT.cartesian_mesh(domain,cells;simplexify=false)
GT.label_boundary_faces!(mesh;physical_name="boundary_faces")
GT.label_interior_faces!(mesh;physical_name="interior_faces")

order = 1

D = GT.num_dims(mesh)
Ω = GT.interior(mesh)
Λ = GT.skeleton(mesh;physical_names=["interior_faces"])
dΛ = GT.measure(Λ,2*order)
n = GT.unit_normal(mesh,D-1)


V = GT.lagrange_space(Ω,1)
GT.reference_face_own_dofs(V,0) |> display
GT.reference_face_own_dofs(V,1) |> display
GT.reference_face_own_dofs(V,2) |> display

display(GT.face_dofs(V))


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
GT.label_boundary_faces!(mesh;physical_name="boundary_faces")

topo = GT.topology(mesh)

GT.face_permutation_ids(topo,2,0)
GT.face_permutation_ids(topo,2,1)
GT.face_permutation_ids(topo,2,2)

Ω = GT.interior(mesh)
Ωref = GT.interior(mesh;is_reference_domain=true)
D = GT.num_dims(mesh)
ϕ = GT.physical_map(mesh,D)

D = GT.num_dims(mesh)
Γdiri = GT.boundary(mesh;physical_names=["boundary_faces"])

V = GT.lagrange_space(Ωref,1;dirichlet_boundary=Γdiri)

#@test V.cache.face_dofs == GT.face_dofs(V)

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
#GT.interpolate!(uref,y,1)

GT.interpolate_free!(uref,y1)
#GT.interpolate_free!(uref,y,1)

GT.interpolate_dirichlet!(uref,y1)
#GT.interpolate_dirichlet!(uref,y,1)

order = 3
V = GT.lagrange_space(Ωref,order)

@test GT.workspace(V).face_dofs == GT.face_dofs(V)

V = GT.lagrange_space(Ωref,order;dirichlet_boundary=Γdiri)

@test GT.workspace(V).face_dofs == GT.face_dofs(V)

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

@test GT.workspace(V).face_dofs == GT.face_dofs(V)

Γ1 = GT.boundary(mesh;physical_names=["1-face-1"])
Γ2 = GT.boundary(mesh;physical_names=["1-face-3"])
Γ3 = GT.boundary(mesh;physical_names=["0-face-1"])

u1 = GT.analytical_field(x->1.0,Ω)
u2 = GT.analytical_field(x->2.0,Ω)
u3 = GT.analytical_field(x->3.0,Ω)

# TODO better names than piecewise_field and piecewise_domain?
udiri = GT.piecewise_field(u1,u2,u3)
Γdiri = GT.piecewise_domain(Γ1,Γ2,Γ3)

order = 1
V = GT.lagrange_space(Ω,order;dirichlet_boundary=Γdiri)

@test GT.workspace(V).face_dofs == GT.face_dofs(V)

V |> GT.dirichlet_dof_location

uh2 = GT.zero_field(Float64,V)
GT.interpolate_dirichlet!(udiri,uh2)

order = 1
V = GT.lagrange_space(Ωref,order;dirichlet_boundary=Γdiri,conformity=:L2)

@test GT.workspace(V).face_dofs == GT.face_dofs(V)

V = GT.lagrange_space(Ωref,order-1;dirichlet_boundary=Γdiri,conformity=:L2)

@test GT.workspace(V).face_dofs == GT.face_dofs(V)

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

@test GT.workspace(V).face_dofs == GT.face_dofs(V)

V = GT.lagrange_space(Γ,order-1;conformity=:L2)
GT.reference_spaces(V) # TODO why 2 reference fes?

@test GT.workspace(V).face_dofs == GT.face_dofs(V)


order = 3
m = GT.analytical_field(x->SVector(false,true),Γ)
V = GT.lagrange_space(Ω,order;tensor_size=Val((2,)),dirichlet_boundary=m)

@test GT.workspace(V).face_dofs == GT.face_dofs(V)

uh = GT.rand_field(Float64,V)

vtk_grid(joinpath(outdir,"Vvec"),Ω;plot_params=(;refinement=10)) do plt
    GT.plot!(plt,uh;label="uh")
    GT.plot!(plt,x->uh(x)[1];label="uh1")
    GT.plot!(plt,x->uh(x)[2];label="uh2")
end

Γ1 = GT.boundary(mesh;physical_names=["1-face-1"])
Γ2 = GT.boundary(mesh;physical_names=["1-face-3"])

order = 3
m1 = GT.analytical_field(x->SVector(false,true),Γ1)
m2 = GT.analytical_field(x->SVector(true,false),Γ2)
m = GT.piecewise_field(m1,m2)
V = GT.lagrange_space(Ω,order;tensor_size=Val((2,)),dirichlet_boundary=m)

@test GT.workspace(V).face_dofs == GT.face_dofs(V)

uh = GT.rand_field(Float64,V)
vtk_grid(joinpath(outdir,"Vvec2"),Ω;plot_params=(;refinement=10)) do plt
    GT.plot!(plt,uh;label="uh")
    GT.plot!(plt,x->uh(x)[1];label="uh1")
    GT.plot!(plt,x->uh(x)[2];label="uh2")
end

order = 0
V = GT.lagrange_space(Ω,order,space_type=:P,dirichlet_boundary=GT.last_dof())

@test GT.workspace(V).face_dofs == GT.face_dofs(V)

uh = GT.rand_field(Float64,V)

vtk_grid(joinpath(outdir,"Vpdisc"),Ω;plot_params=(;refinement=10)) do plt
    GT.plot!(plt,uh;label="uh")
end
end # module

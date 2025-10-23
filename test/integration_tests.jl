module IntegrationTests

using Test
import GalerkinToolkit as GT
import PartitionedArrays as pa
using GalerkinToolkit: ∫
using LinearAlgebra
import ForwardDiff
using AbstractTrees

spx0 = GT.unit_simplex(0)
spx1 = GT.unit_simplex(1)
spx2 = GT.unit_simplex(2)
spx3 = GT.unit_simplex(3)

cube0 = GT.unit_n_cube(0)
cube1 = GT.unit_n_cube(1)
cube2 = GT.unit_n_cube(2)
cube3 = GT.unit_n_cube(3)

degree = 4
quad = GT.quadrature(spx0,degree)
quad = GT.quadrature(spx1,degree)
quad = GT.quadrature(spx2,degree)
quad = GT.quadrature(spx3,degree)

quad = GT.quadrature(cube0,degree)
@test sum(GT.weights(quad)) ≈ 1
quad = GT.quadrature(cube1,degree)
@test sum(GT.weights(quad)) ≈ 1
quad = GT.quadrature(cube2,degree)
@test sum(GT.weights(quad)) ≈ 1
quad = GT.quadrature(cube3,degree)
@test sum(GT.weights(quad)) ≈ 1

outdir = mkpath(joinpath(@__DIR__,"..","output"))

domain = (0,2,0,2)
cells = (8,8)
mesh = GT.cartesian_mesh(domain,cells)

D = GT.num_dims(mesh)
Ω = GT.interior(mesh)
Ωref = GT.interior(mesh;is_reference_domain=true)
ϕ = GT.physical_map(mesh,D)
u = GT.analytical_field(x->sum(x),Ω)

degree = 2
dΩref = GT.measure(Ωref,degree)
int = ∫(dΩref) do q
    x = ϕ(q)
    J = ForwardDiff.jacobian(ϕ,q)
    dV = abs(det(J))
    u(x)*dV
end
sum(int)
@test_broken sum(int) ≈ 8


int = ∫(dΩref) do q
    J = ForwardDiff.jacobian(ϕ,q)
    dV = abs(det(J))
    dV
end

@test_broken sum(int) ≈ 4

dΩ = GT.measure(Ω,degree)
int = ∫(dΩ) do x
    u(x)
end

@test sum(int) ≈ 8

u = GT.analytical_field(x->1,Ω)
int = ∫(u,dΩ)
@test sum(int) ≈ 4

D = GT.num_dims(mesh)
Γref = GT.boundary(mesh;
                 is_reference_domain=true,
                 group_names=["1-face-2","1-face-4"])

Γ = GT.physical_domain(Γref)
h = GT.face_diameter_field(Γ)
n = GT.unit_normal(mesh,D-1)

function dS(J)
    Jt = transpose(J)
    sqrt(det(Jt*J))
end

dΓref = GT.measure(Γref,degree)
α = GT.physical_map(mesh,D-1)
#β = GT.reference_map(mesh,D-1,D)

int = ∫(dΓref) do p
    J = ForwardDiff.jacobian(α,p)
    dS(J)
end
@test_broken sum(int) ≈ 4

dΓ = GT.measure(Γ,degree)
int = ∫(x->norm(n(x)),dΓ)
sum_int = sum(int)
@test sum_int ≈ 4

uref = GT.analytical_field(x->1,Γ)
#int = ∫(dΓref) do p
#    p
#    q = β(p)
#    J = ForwardDiff.jacobian(α,p)
#    uref(q)*dS(J)
#end
#sum(int) ≈ 4

int = ∫(dΩref) do q
    J = ForwardDiff.jacobian(ϕ,q)
    dV = abs(det(J))
    dV
end +
∫(dΓref) do p
    J = ForwardDiff.jacobian(α,p)
    dS(J)
end

@test_broken sum(int) ≈ 8

#fd = GT.face_diameter(Γ)
#@test sum(fd) ≈ 4

h = GT.face_diameter_field(Γ)

Λref = GT.skeleton(mesh;is_reference_domain=true)

Λ = GT.physical_domain(Λref)
dΛref = GT.measure(Λref,degree)
dΛ = GT.measure(Λ,degree)
ϕ_Λref_Λ = GT.physical_map(mesh,D-1)
#ϕ_Λref_Ωref = GT.reference_map(mesh,D-1,D)

#jump(u,ϕ,q) = u(ϕ(q)[2])-u(ϕ(q)[1])

#int = 10*∫(dΛref) do p
#    #q = ϕ_Λref_Ωref(p)[2]
#    #index = GT.generate_index(Γref)
#    #t = GT.term(q,index)
#    #print_tree(t)
#    #q
#    #xxx
#    J = ForwardDiff.jacobian(ϕ_Λref_Λ,p)
#    3*jump(uref,ϕ_Λref_Ωref,p)*dS(J)
#end
#s = sum(int*1)
#@test s + 1 ≈ 1
#@test sum(int/1) + 1 ≈ 1

int = ∫(x->norm(n[1](x)+n[2](x)),dΛ)
sum_int = sum(int)
@test sum_int  + 1 ≈ 1

int = 10*∫(dΛref) do p
    J = ForwardDiff.jacobian(ϕ_Λref_Λ,p)
    h(p)*dS(J)
    #3*h(p)*dS(J)
end
s = sum(int*1)

domain = (0,1,0,1)
cells = (4,4)
mesh = GT.cartesian_mesh(domain,cells)
GT.group_interior_faces!(mesh;group_name="interior_faces")
GT.group_boundary_faces!(mesh;group_name="boundary_faces")
Ω = GT.interior(mesh)
Γ = GT.boundary(mesh;group_names=["boundary_faces"])
Λ = GT.skeleton(mesh;group_names=["interior_faces"])

order = 2
degree = 2*order
dΩ = GT.measure(Ω,degree)
dΓ = GT.measure(Γ,degree)
dΛ = GT.measure(Λ,degree)

V = GT.lagrange_space(Ω,order;dirichlet_boundary=Γ)

#face_point_x = GT.coordinate_accessor(dΩ)
#face_point_J = GT.jacobian_accessor(dΩ)
#face_point_dV = GT.weight_accessor(dΩ)
#face_point_dof_s = GT.shape_function_accessor(GT.value,V,dΩ)
#face_point_dof_∇s = GT.shape_function_accessor(ForwardDiff.gradient,V,dΩ)
#
#face = 3
#point = 4
#dof = 2
#
#point_x = face_point_x(face)
#point_J = face_point_J(face)
#point_dV = face_point_dV(face)
#point_dof_s = face_point_dof_s(face)
#point_dof_∇s = face_point_dof_∇s(face)
#
#x = point_x(point)
#J = point_J(point)
#dV = point_dV(point,J)


# Parallel

#domain = (0,2,0,2)
#cells_per_dir = (4,4)
#parts_per_dir = (2,2)
#np = prod(parts_per_dir)
#parts = pa.DebugArray(LinearIndices((np,)))
#partition_strategy = GT.partition_strategy(graph_nodes=:cells,graph_edges=:nodes,ghost_layers=0)
#mesh = GT.cartesian_mesh(domain,cells_per_dir;parts_per_dir,parts,partition_strategy)
#
#GT.group_boundary_faces!(mesh;group_name="boundary_faces")
#Ω = GT.interior(mesh)
#Ωref = GT.interior(mesh;is_reference_domain=true)
#u = GT.analytical_field(x->sum(x),Ω)
#ϕ = GT.domain_map(Ωref,Ω)
#uref = u∘ϕ
#
#degree = 2
#dΩref = GT.measure(Ωref,degree)
#int = ∫(dΩref) do q
#    x = ϕ(q)
#    J = ForwardDiff.jacobian(ϕ,q)
#    dV = abs(det(J))
#    u(x)*dV
#end
#
#s = sum(int)
#@test s ≈ 8
#
#int = ∫(dΩref) do q
#    J = ForwardDiff.jacobian(ϕ,q)
#    dV = abs(det(J))
#    dV
#end
#
#@test sum(int) ≈ 4
#
#dΩ = GT.measure(Ω,degree)
#int = ∫(dΩ) do x
#    u(x)
#end
#
#@test sum(int) ≈ 8
#
#u = GT.analytical_field(x->1,Ω)
#int = ∫(u,dΩ)
#@test sum(int) ≈ 4
#
#D = GT.num_dims(mesh)
#Γref = GT.boundary(mesh;
#                 is_reference_domain=true,
#                 group_names=["1-face-2","1-face-4"])
#
#Γ = GT.physical_domain(Γref)
#
#function dS(J)
#    Jt = transpose(J)
#    sqrt(det(Jt*J))
#end
#
#dΓref = GT.measure(Γref,degree)
#α = GT.domain_map(Γref,Γ)
#
#β = GT.domain_map(Γref,Ωref)
#
#int = ∫(dΓref) do p
#    J = ForwardDiff.jacobian(α,p)
#    dS(J)
#end
#@test sum(int) ≈ 4
#
#uref = GT.analytical_field(x->1,Ωref)
#int = ∫(dΓref) do p
#    q = β(p)
#    J = ForwardDiff.jacobian(α,p)
#    uref(q)*dS(J)
#end
#sum(int) ≈ 4
#
#int = ∫(dΩref) do q
#    J = ForwardDiff.jacobian(ϕ,q)
#    dV = abs(det(J))
#    dV
#end +
#∫(dΓref) do p
#    J = ForwardDiff.jacobian(α,p)
#    dS(J)
#end
#
#@test sum(int) ≈ 8


end # module

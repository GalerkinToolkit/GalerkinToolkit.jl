module IntegrationTests

using Test
import GalerkinToolkit as gk
import PartitionedArrays as pa
using GalerkinToolkit: ∫
using LinearAlgebra
import ForwardDiff

spx0 = gk.unit_simplex(0)
spx1 = gk.unit_simplex(1)
spx2 = gk.unit_simplex(2)
spx3 = gk.unit_simplex(3)

cube0 = gk.unit_n_cube(0)
cube1 = gk.unit_n_cube(1)
cube2 = gk.unit_n_cube(2)
cube3 = gk.unit_n_cube(3)

degree = 4
quad = gk.default_quadrature(spx0,degree)
quad = gk.default_quadrature(spx1,degree)
quad = gk.default_quadrature(spx2,degree)
quad = gk.default_quadrature(spx3,degree)

quad = gk.default_quadrature(cube0,degree)
@test sum(gk.weights(quad)) ≈ 1
quad = gk.default_quadrature(cube1,degree)
@test sum(gk.weights(quad)) ≈ 1
quad = gk.default_quadrature(cube2,degree)
@test sum(gk.weights(quad)) ≈ 1
quad = gk.default_quadrature(cube3,degree)
@test sum(gk.weights(quad)) ≈ 1

outdir = mkpath(joinpath(@__DIR__,"..","output"))

domain = (0,2,0,2)
cells = (8,8)
mesh = gk.cartesian_mesh(domain,cells)
gk.label_boundary_faces!(mesh;physical_name="boundary_faces")
gk.label_interior_faces!(mesh;physical_name="interior_faces")

Ω = gk.domain(mesh)
Ωref = gk.domain(mesh;is_reference_domain=true)
ϕ = gk.domain_map(Ωref,Ω)
u = gk.analytical_field(x->sum(x),Ω)

degree = 2
dΩref = gk.measure(Ωref,degree)
int = ∫(dΩref) do q
    x = ϕ(q)
    J = ForwardDiff.jacobian(ϕ,q)
    dV = abs(det(J))
    u(x)*dV
end

@test sum(int) ≈ 8

int = ∫(dΩref) do q
    J = ForwardDiff.jacobian(ϕ,q)
    dV = abs(det(J))
    dV
end

@test sum(int) ≈ 4

dΩ = gk.measure(Ω,degree)
int = ∫(dΩ) do x
    u(x)
end

@test sum(int) ≈ 8

u = gk.analytical_field(x->1,Ω)
int = ∫(u,dΩ)
@test sum(int) ≈ 4


D = gk.num_dims(mesh)
Γref = gk.domain(mesh;
                 face_dim=D-1,
                 is_reference_domain=true,
                 physical_names=["1-face-2","1-face-4"])

Γ = gk.physical_domain(Γref)

function dS(J)
    Jt = transpose(J)
    sqrt(det(Jt*J))
end

dΓref = gk.measure(Γref,degree)
α = gk.domain_map(Γref,Γ)

β = gk.domain_map(Γref,Ωref;face_around=1)

int = ∫(dΓref) do p
    J = ForwardDiff.jacobian(α,p)
    dS(J)
end
@test sum(int) ≈ 4

uref = gk.analytical_field(x->1,Ωref)
int = ∫(dΓref) do p
    q = β(p)
    J = ForwardDiff.jacobian(α,p)
    uref(q)*dS(J)
end
sum(int) ≈ 4

int = ∫(dΩref) do q
    J = ForwardDiff.jacobian(ϕ,q)
    dV = abs(det(J))
    dV
end +
∫(dΓref) do p
    J = ForwardDiff.jacobian(α,p)
    dS(J)
end

@test sum(int) ≈ 8

@test sum(gk.face_diameter(Γ)) ≈ 4

h = gk.face_diameter_field(Γ)

Λref = gk.domain(mesh;
                 face_dim=D-1,
                 is_reference_domain=true,
                 physical_names=["interior_faces"])

Λ = gk.physical_domain(Λref)
dΛref = gk.measure(Λref,degree)
ϕ_Λref_Λ = gk.domain_map(Λref,Λ)
ϕ_Λref_Ωref = gk.domain_map(Λref,Ωref)

jump(u) = u[2]-u[1]

int = 10*∫(dΛref) do p
    q = ϕ_Λref_Ωref(p)
    J = ForwardDiff.jacobian(ϕ_Λref_Λ,p)
    3*jump(uref(q))*dS(J)
end

@test sum(int*1) + 1 ≈ 1
@test sum(int/1) + 1 ≈ 1

# Parallel

domain = (0,2,0,2)
cells_per_dir = (4,4)
parts_per_dir = (2,2)
np = prod(parts_per_dir)
parts = pa.DebugArray(LinearIndices((np,)))
partition_strategy = gk.partition_strategy(graph_nodes=:cells,graph_edges=:nodes,ghost_layers=0)
mesh = gk.cartesian_mesh(domain,cells_per_dir;parts_per_dir,parts,partition_strategy)

gk.label_boundary_faces!(mesh;physical_name="boundary_faces")
Ω = gk.domain(mesh)
Ωref = gk.domain(mesh;is_reference_domain=true)
u = gk.analytical_field(x->sum(x),Ω)
ϕ = gk.domain_map(Ωref,Ω)
uref = u∘ϕ

degree = 2
dΩref = gk.measure(Ωref,degree)
int = ∫(dΩref) do q
    x = ϕ(q)
    J = ForwardDiff.jacobian(ϕ,q)
    dV = abs(det(J))
    u(x)*dV
end

@test sum(int) ≈ 8

int = ∫(dΩref) do q
    J = ForwardDiff.jacobian(ϕ,q)
    dV = abs(det(J))
    dV
end

@test sum(int) ≈ 4

dΩ = gk.measure(Ω,degree)
int = ∫(dΩ) do x
    u(x)
end

@test sum(int) ≈ 8

u = gk.analytical_field(x->1,Ω)
int = ∫(u,dΩ)
@test sum(int) ≈ 4

D = gk.num_dims(mesh)
Γref = gk.domain(mesh;
                 face_dim=D-1,
                 is_reference_domain=true,
                 physical_names=["1-face-2","1-face-4"])

Γ = gk.physical_domain(Γref)

function dS(J)
    Jt = transpose(J)
    sqrt(det(Jt*J))
end

dΓref = gk.measure(Γref,degree)
α = gk.domain_map(Γref,Γ)

β = gk.domain_map(Γref,Ωref;face_around=1)

int = ∫(dΓref) do p
    J = ForwardDiff.jacobian(α,p)
    dS(J)
end
@test sum(int) ≈ 4

uref = gk.analytical_field(x->1,Ωref)
int = ∫(dΓref) do p
    q = β(p)
    J = ForwardDiff.jacobian(α,p)
    uref(q)*dS(J)
end
sum(int) ≈ 4

int = ∫(dΩref) do q
    J = ForwardDiff.jacobian(ϕ,q)
    dV = abs(det(J))
    dV
end +
∫(dΓref) do p
    J = ForwardDiff.jacobian(α,p)
    dS(J)
end

@test sum(int) ≈ 8


end # module

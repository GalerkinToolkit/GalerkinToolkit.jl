module AssemblyTests

import GalerkinToolkit as gk
using GalerkinToolkit: ∫, ×
using Test
import ForwardDiff
using LinearAlgebra

outdir = mkpath(joinpath(@__DIR__,"..","output"))

domain = (0,1,0,1)
cells = (4,4)
mesh = gk.cartesian_mesh(domain,cells)
gk.label_interior_faces!(mesh;physical_name="interior_faces")
gk.label_boundary_faces!(mesh;physical_name="boundary_faces")

Ω = gk.domain(mesh)
Ωref = gk.domain(mesh;is_reference_domain=true)
ϕ = gk.domain_map(Ωref,Ω)

D = gk.num_dims(mesh)
Γdiri = gk.domain(mesh;face_dim=D-1,physical_names=["1-face-1","1-face-3"])

Γref = gk.domain(mesh;
                 face_dim=D-1,
                 is_reference_domain=true,
                 physical_names=["1-face-2","1-face-4"])

Γ = gk.physical_domain(Γref)

V = gk.iso_parametric_space(Ωref;dirichlet_boundary=Γdiri)

degree = 2
dΩref = gk.measure(Ωref,degree)
ϕ = gk.domain_map(Ωref,Ω)

dΓref = gk.measure(Γref,degree)
α = gk.domain_map(Γref,Γ)
β = gk.domain_map(Γref,Ωref;face_around=1)

Λref = gk.domain(mesh;
                 face_dim=D-1,
                 is_reference_domain=true,
                 physical_names=["interior_faces"])

Λ = gk.physical_domain(Λref)
dΛref = gk.measure(Λref,degree)
ϕ_Λref_Λ = gk.domain_map(Λref,Λ)
ϕ_Λref_Ωref = gk.domain_map(Λref,Ωref)

function dV(J)
    abs(det(J))
end

function dS(J)
    Jt = transpose(J)
    sqrt(det(Jt*J))
end

jump(u) = u[2]-u[1]

function l(v)
    ∫(dΩref) do q
        J = ForwardDiff.jacobian(ϕ,q)
        v(q)*dV(J)
    end +
    ∫(dΓref) do p
        q = β(p)
        J = ForwardDiff.jacobian(α,p)
        v(q)*dS(J)
    end +
    ∫(dΛref) do p
        q = ϕ_Λref_Ωref(p)
        J = ForwardDiff.jacobian(ϕ_Λref_Λ,p)
        jump(v(q))*dS(J)
    end
end

b = gk.assemble_vector(l,V)

function l(v)
    ∫(dΛref) do p
        q = ϕ_Λref_Ωref(p)
        J = ForwardDiff.jacobian(ϕ_Λref_Λ,p)
        jump(v(q))*dS(J)
    end
end

b = gk.assemble_vector(l,V)
@test sum(b)+1 ≈ 1

V² = V × V

function l((v1,v2))
    ∫(dΩref) do q
        J = ForwardDiff.jacobian(ϕ,q)
        v1(q)*v2(q)*dV(J)
    end +
    ∫(dΛref) do p
        q = ϕ_Λref_Ωref(p)
        J = ForwardDiff.jacobian(ϕ_Λref_Λ,p)
        jump(v1(q))*jump(v2(q))*dS(J)
    end
end

b = gk.assemble_vector(l,V²)

function a(u,v)
    ∫(dΩref) do q
        J = ForwardDiff.jacobian(ϕ,q)
        u(q)*v(q)*dV(J)
    end +
    ∫(dΓref) do p
        q = β(p)
        J = ForwardDiff.jacobian(α,p)
        u(q)*v(q)*dS(J)
    end +
    ∫(dΛref) do p
        q = ϕ_Λref_Ωref(p)
        J = ForwardDiff.jacobian(ϕ_Λref_Λ,p)
        jump(v(q))*jump(u(q))*dS(J)
    end
end

A = gk.assemble_matrix(a,V,V)

function a((u1,u2),(v1,v2))
    ∫(dΩref) do q
        J = ForwardDiff.jacobian(ϕ,q)
        v1(q)*v2(q)*u1(q)*u2(q)*dV(J)
    end +
    ∫(dΛref) do p
        q = ϕ_Λref_Ωref(p)
        J = ForwardDiff.jacobian(ϕ_Λref_Λ,p)
        jump(v1(q))*jump(v2(q))*jump(u1(q))*jump(u2(q))*dS(J)
    end
end

A = gk.assemble_matrix(a,V²,V²)

function dV(ϕ,q)
    J = ForwardDiff.jacobian(ϕ,q)
    abs(det(J))
end

a(u,v) = ∫( q->u(q)*v(q)*dV(ϕ,q), dΩref)

f = gk.analytical_field(sum,Ω)

l(v) = ∫( q->f(ϕ(q))*v(q)*dV(ϕ,q), dΩref)

V = gk.iso_parametric_space(Ωref)
uh = gk.zero_field(Float64,V)

x,A,b = gk.linear_problem(uh,a,l)
x .= A\b

function ∇(u,phi,q)
   J = ForwardDiff.jacobian(phi,q)
   g = ForwardDiff.gradient(u,q)
   J\g
end

a(u,v) = ∫( q->∇(u,ϕ,q)⋅∇(v,ϕ,q)*dV(ϕ,q), dΩref)
l(v) = 0

x,A,b = gk.linear_problem(uh,a,l)
x .= A\b

# Poisson solve

domain = (0,1,0,1)
cells = (4,4)
mesh = gk.cartesian_mesh(domain,cells)
gk.label_boundary_faces!(mesh;physical_name="boundary_faces")

Ω = gk.domain(mesh)
Ωref = gk.reference_domain(Ω)
ϕ = gk.domain_map(Ωref,Ω)

D = gk.num_dims(mesh)
Γdiri = gk.domain(mesh;face_dim=D-1,physical_names=["boundary_faces"])

#V = gk.iso_parametric_space(Ωref;dirichlet_boundary=Γdiri)

# TODO not working for order > 2
order = 2
V = gk.lagrange_space(Ωref,order;dirichlet_boundary=Γdiri)

u = gk.analytical_field(sum,Ω)
uh = gk.zero_field(Float64,V)
# TODO
#gk.interpolate_dirichlet!(q->u(ϕ(q)),uh)
gk.interpolate_dirichlet!(u∘ϕ,uh)

function ∇(u,q)
   J = ForwardDiff.jacobian(ϕ,q)
   g = ForwardDiff.gradient(u,q)
   J\g
end

function dV(q)
    J = ForwardDiff.jacobian(ϕ,q)
    abs(det(J))
end

degree = 2
dΩref = gk.measure(Ωref,degree)

a(u,v) = ∫( q->∇(u,q)⋅∇(v,q)*dV(q), dΩref)
l(v) = 0

x,A,b = gk.linear_problem(uh,a,l)
x .= A\b


# TODO
# Functions like this ones should
# work as AbstractQuantities?
eh(q) = u(ϕ(q)) - uh(q)
∇eh(q) = ForwardDiff.gradient(u,ϕ(q)) - ∇(uh,q)

tol = 1.e-12
el2 = ∫( q->abs2(eh(q))*dV(q), dΩref) |> sum |> sqrt
@test el2 < tol

eh1 = ∫( q->∇eh(q)⋅∇eh(q)*dV(q), dΩref) |> sum |> sqrt
@test el2 < tol

end # module

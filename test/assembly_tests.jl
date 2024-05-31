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
    1
end

jump(u) = u[2]-u[1]

function l(v)
    ∫(dΩref) do q
        J = gk.call(ForwardDiff.jacobian,ϕ,q)
        dVq = gk.call(dV,J)
        gk.call(*,v(q),dVq)
    end +
    ∫(dΓref) do p
        q = β(p)
        J = gk.call(ForwardDiff.jacobian,α,p)
        dSq = gk.call(dS,J)
        gk.call(*,v(q),dSq)
        v(q)
    end +
    ∫(dΛref) do p
        q = ϕ_Λref_Ωref(p)
        J = gk.call(ForwardDiff.jacobian,ϕ_Λref_Λ,p)
        dSq = gk.call(dS,J)
        jvq = gk.call(jump,v(q))
        gk.call(*,jvq,dSq)
    end
end

b = gk.assemble_vector(l,V)

function l(v)
    ∫(dΛref) do p
        q = ϕ_Λref_Ωref(p)
        J = gk.call(ForwardDiff.jacobian,ϕ_Λref_Λ,p)
        dSq = gk.call(dS,J)
        jvq = gk.call(jump,v(q))
        gk.call(*,jvq,dSq)
    end
end

b = gk.assemble_vector(l,V)
@test sum(b)+1 ≈ 1

V² = V × V

function l((v1,v2))
    ∫(dΩref) do q
        J = gk.call(ForwardDiff.jacobian,ϕ,q)
        dVq = gk.call(dV,J)
        v = gk.call(*,v1(q),v2(q))
        gk.call(*,v,dVq)
    end
end

b = gk.assemble_vector(l,V²)

end # module

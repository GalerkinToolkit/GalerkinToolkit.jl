module AssemblyTests

import GalerkinToolkit as gk
using GalerkinToolkit: ∫
using Test
import ForwardDiff
using LinearAlgebra

outdir = mkpath(joinpath(@__DIR__,"..","output"))

domain = (0,1,0,1)
cells = (4,4)
mesh = gk.cartesian_mesh(domain,cells)

Ω = gk.domain(mesh)
Ωref = gk.domain(mesh;is_reference_domain=true)
ϕ = gk.domain_map(Ωref,Ω)

D = gk.num_dims(mesh)
Γdiri = gk.domain(mesh;face_dim=D-1)

V = gk.iso_parametric_space(Ωref;dirichlet_boundary=Γdiri)

degree = 2
dΩref = gk.measure(Ωref,degree)
ϕ = gk.domain_map(Ωref,Ω)

l((v,)) = ∫(dΩref) do q
    J = gk.call(ForwardDiff.jacobian,ϕ,q)
    dV = gk.call(J->abs(det(J)),J)
    gk.call(*,v(q),dV)
end

b = gk.assemble_vector(l,V)

display(b)






end # module

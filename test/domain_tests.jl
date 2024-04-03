module DomainTests

import GalerkinToolkit as gk
using Test

outdir = mkpath(joinpath(@__DIR__,"..","output"))

domain = (0,1,0,1)
cells = (4,4)
mesh = gk.cartesian_mesh(domain,cells)

Ω = gk.domain(mesh)
Ωref = gk.domain(mesh;is_reference_domain=true)

@test Ω == Ω
@test Ω != Ωref

u = gk.analytical_field(x->sum(x),Ω)

gk.vtk_plot(joinpath(outdir,"omega"),Ω;refinement=4) do plt
    gk.plot!(plt,u;label="u")
end

ϕ = gk.domain_map(Ωref,Ω)
uref = u∘ϕ

gk.vtk_plot(joinpath(outdir,"omega_ref"),Ωref;refinement=4) do plt
    gk.plot!(plt,uref;label="u")
end

D = gk.num_dims(mesh)
Γref = gk.domain(mesh;face_dim=D-1,is_reference_domain=true)

ϕ = gk.domain_map(Γref,Ωref;face_around=1)
g = uref∘ϕ

gk.vtk_plot(joinpath(outdir,"gamma_ref"),Γref) do plt
    gk.plot!(plt,g;label="u")
    gk.plot!(plt;label="u2") do q
        x = ϕ(q)
        uref(x)
    end
end


end # module

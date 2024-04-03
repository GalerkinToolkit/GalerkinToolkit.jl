module DomainTests

import GalerkinToolkit as gk
using Test

outdir = mkpath(joinpath(@__DIR__,"..","output"))

domain = (0,1,0,1)
cells = (4,4)
mesh = gk.cartesian_mesh(domain,cells)
gk.label_interior_faces!(mesh;physical_name="interior_faces")
gk.label_boundary_faces!(mesh;physical_name="boundary_faces")

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
Γref = gk.domain(mesh;
                 face_dim=D-1,
                 is_reference_domain=true,
                 physical_names=["boundary_faces"])

ϕ = gk.domain_map(Γref,Ωref;face_around=1)
g = uref∘ϕ

gk.vtk_plot(joinpath(outdir,"gamma_ref"),Γref) do plt
    gk.plot!(plt,g;label="u")
    gk.plot!(plt;label="u2") do q
        x = ϕ(q)
        uref(x)
    end
end

Λref = gk.domain(mesh;
                 face_dim=D-1,
                 is_reference_domain=true,
                 physical_names=["interior_faces"])

ϕ⁺,ϕ⁻ = gk.domain_map(Λref,Ωref)

gk.vtk_plot(joinpath(outdir,"lambda_ref"),Λref) do plt
    gk.plot!(plt;label="jump_u") do q
        x⁺ = ϕ⁺(q)
        x⁻ = ϕ⁻(q)
        gk.call(-,uref(x⁺),uref(x⁻))
    end
end

end # module

module InterpolationTests

import GalerkinToolkit as gk
using GalerkinToolkit: ×
using Test

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

v = gk.zero_field(Float64,V)
v2 = gk.zero_field(Float64,V)

u = gk.analytical_field(x->sum(x),Ω)
uref = u∘ϕ

gk.interpolate!(uref,v)
gk.interpolate_dirichlet!(uref,v2)

Y = V×V
y = gk.zero_field(Float64,Y)
y1, y2 = y

gk.vtk_plot(joinpath(outdir,"omega_ref"),Ωref;refinement=4) do plt
    gk.plot!(plt,v;label="v")
    gk.plot!(plt,v2;label="v2")
    gk.plot!(plt,y1;label="y1")
    gk.plot!(plt,y2;label="y2")
end

end #module

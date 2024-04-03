module DomainTests

import GalerkinToolkit as gk

outdir = mkpath(joinpath(@__DIR__,"..","output"))

domain = (0,1,0,1)
cells = (4,4)
mesh = gk.cartesian_mesh(domain,cells)
Ω = gk.domain(mesh)
u = gk.analytical_field(x->sum(x),Ω)

gk.vtk_plot(joinpath(outdir,"omega"),Ω;refinement=4) do plt
    gk.plot!(plt,u;label="u")
end


end # module

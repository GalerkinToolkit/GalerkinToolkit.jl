module GmshTests

import GalerkinToolkit as GT

outdir = mkpath(joinpath(@__DIR__,"..","output"))
msh =  joinpath(@__DIR__,"..","assets","quad.msh")
mesh = GT.mesh_from_msh(msh)

mesh = GT.mesh_from_msh(msh,complexify=false)
mesh = GT.complexify(mesh)

msh =  joinpath(@__DIR__,"..","assets","demo.msh")
mesh = GT.mesh_from_msh(msh)

end # module

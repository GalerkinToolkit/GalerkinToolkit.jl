module CartesianMeshTests

import GalerkinToolkit as GT

#using InteractiveUtils

domain = (0,1,0,1,0,1)
cells = (2,2,2)
mesh = GT.cartesian_mesh(domain,cells)
smesh = GT.simplexify(mesh)

domain = (0,1,0,1)
cells = (2,2)

mesh = GT.cartesian_mesh(domain,cells,boundary=false)
mesh = GT.cartesian_mesh(domain,cells,simplexify=true)
mesh = GT.cartesian_mesh(domain,cells,boundary=false,simplexify=true)
mesh = GT.cartesian_mesh(domain,cells)

end #module

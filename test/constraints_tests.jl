module ConstraintsTests

import GalerkinToolkit as GT

domain = (0,1,0,1)
cells = (2,2)
mesh = GT.cartesian_mesh(domain,cells;periodic=(true,false))
Ω = GT.interior(mesh)
Γ = GT.boundary(mesh)

order = 1

V0 = GT.lagrange_space(Ω,order)
C = GT.periodic_constraints(V0)

display(C)

V = GT.constrain(V0,C)

end # module

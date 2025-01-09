module TopologyTests

import GalerkinToolkit as GT
using Test

using InteractiveUtils

spx0 = GT.unit_simplex(0)
spx1 = GT.unit_simplex(1)
spx2 = GT.unit_simplex(2)
spx3 = GT.unit_simplex(3)

cube0 = GT.unit_n_cube(0)
cube1 = GT.unit_n_cube(1)
cube2 = GT.unit_n_cube(2)
cube3 = GT.unit_n_cube(3)

mesh = GT.complexify(spx0)
mesh = GT.complexify(spx1)
mesh = GT.complexify(spx2)
mesh = GT.complexify(spx3)

mesh = GT.complexify(cube0)
mesh = GT.complexify(cube1)
mesh = GT.complexify(cube2)
mesh = GT.complexify(cube3)

order = 1
fe = GT.lagrange_space(cube2;order)
mesh = GT.complexify(fe)

order = 2
fe = GT.lagrange_space(cube2;order)
mesh = GT.complexify(fe)

order = 3
fe = GT.lagrange_space(cube2;order)
mesh = GT.complexify(fe)

display(GT.face_nodes(mesh,0))
display(GT.face_nodes(mesh,1))
display(GT.face_nodes(mesh,2))

#@code_warntype GT.complexify(cube0)

#@code_warntype GT.complexify(fe)



end # module

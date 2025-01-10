module DomainTests

import GalerkinToolkit as GT
import StaticArrays as SA
using Test
using LinearAlgebra

using InteractiveUtils

#@code_warntype GT.unit_simplex(Val(3))

spx0 = GT.unit_simplex(0)
spx1 = GT.unit_simplex(1)
spx2 = GT.unit_simplex(2)
spx3 = GT.unit_simplex(3)
@test GT.num_dims(spx3) == 3
display(spx3)

cube0 = GT.unit_n_cube(0)
cube1 = GT.unit_n_cube(1)
cube2 = GT.unit_n_cube(2)
cube3 = GT.unit_n_cube(3)
@test GT.num_dims(cube3) == 3
display(cube3)

#@code_warntype GT.options(cube0)
#@code_warntype GT.options()

options = GT.options(cube0)
@show GT.real_type(options)
@show GT.int_type(options)

mesh = GT.mesh(spx0)
mesh = GT.mesh(spx1)
mesh = GT.mesh(spx2)
mesh = GT.mesh(spx3)

mesh = GT.mesh(cube0)
mesh = GT.mesh(cube1)
mesh = GT.mesh(cube2)
mesh = GT.mesh(cube3)

@show GT.faces(cube3)
@show GT.inverse_faces(cube3)
@show GT.geometries(cube3)
@show GT.geometries(cube3,2)
@show GT.num_geometries(cube3,2)

perms = GT.vertex_permutations(cube2)

#@code_warntype GT.mesh(cube3)

#@code_warntype GT.vertex_permutations(cube2)

mesh = GT.mesh(cube3)

domain = GT.domain(mesh,2)

@test mesh === GT.mesh(domain)

GT.faces(domain)
GT.inverse_faces(domain)

mesh = GT.simplexify(cube0)
mesh = GT.simplexify(cube1)
mesh = GT.simplexify(cube2)
mesh = GT.simplexify(cube3)

mesh = GT.simplexify(spx0)
mesh = GT.simplexify(spx1)
mesh = GT.simplexify(spx2)
mesh = GT.simplexify(spx3)

end # module

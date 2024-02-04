module IntegrationTests

using Test
import GalerkinToolkit as gk

spx0 = gk.unit_simplex(0)
spx1 = gk.unit_simplex(1)
spx2 = gk.unit_simplex(2)
spx3 = gk.unit_simplex(3)

cube0 = gk.unit_n_cube(0)
cube1 = gk.unit_n_cube(1)
cube2 = gk.unit_n_cube(2)
cube3 = gk.unit_n_cube(3)

degree = 4
quad = gk.default_quadrature(spx0,degree)
quad = gk.default_quadrature(spx1,degree)
quad = gk.default_quadrature(spx2,degree)
quad = gk.default_quadrature(spx3,degree)

quad = gk.default_quadrature(cube0,degree)
@test sum(gk.weights(quad)) ≈ 1
quad = gk.default_quadrature(cube1,degree)
@test sum(gk.weights(quad)) ≈ 1
quad = gk.default_quadrature(cube2,degree)
@test sum(gk.weights(quad)) ≈ 1
quad = gk.default_quadrature(cube3,degree)
@test sum(gk.weights(quad)) ≈ 1


end # module

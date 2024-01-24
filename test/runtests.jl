module GalerkinToolkitTests

using Test
import GalerkinToolkit as gt

spx0 = gt.unit_simplex(0)
spx1 = gt.unit_simplex(1)
spx2 = gt.unit_simplex(2)
spx3 = gt.unit_simplex(3)
display(spx3)

cube0 = gt.unit_n_cube(0)
cube1 = gt.unit_n_cube(1)
cube2 = gt.unit_n_cube(2)
cube3 = gt.unit_n_cube(3)
display(cube3)

quad = gt.default_quadrature(spx0;degree=2)
quad = gt.default_quadrature(spx1;degree=2)
quad = gt.default_quadrature(spx2;degree=2)
quad = gt.default_quadrature(spx3;degree=2)

quad = gt.default_quadrature(cube0;degree=2)
quad = gt.default_quadrature(cube1;degree=2)
quad = gt.default_quadrature(cube2;degree=2)
quad = gt.default_quadrature(cube3;degree=2)

quad = gt.default_quadrature(cube1;degree=4,real_type=Float32)
quad = gt.default_quadrature(spx1;degree=4,real_type=Float32)



end # module

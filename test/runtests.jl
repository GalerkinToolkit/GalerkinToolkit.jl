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

@show typeof(spx0)

fe = gt.lagrangian_fe(spx0,order=1)
fe = gt.lagrangian_fe(spx1,order=1)
fe = gt.lagrangian_fe(spx2,order=1)
fe = gt.lagrangian_fe(spx3,order=1)
display(fe)

fe = gt.lagrangian_fe(cube0,order=1)
fe = gt.lagrangian_fe(cube1,order=1)
fe = gt.lagrangian_fe(cube2,order=1)
fe = gt.lagrangian_fe(cube3,order=1)
display(fe)

fe = gt.lagrangian_fe(cube2,order=3)
@show gt.node_coordinates(fe)

#∂spx0 = gt.boundary(spx0)
#∂spx0 = gt.boundary(spx1)
#∂cube0 = gt.boundary(cube0)



end # module

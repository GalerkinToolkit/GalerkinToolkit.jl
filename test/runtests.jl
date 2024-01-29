module GalerkinToolkitTests

using Test
import GalerkinToolkit as gt
using WriteVTK

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

@show typeof(spx0)
@show typeof(cube0)
@show typeof(spx1)
@show typeof(cube1)

degree = 4
quad = gt.default_quadrature(spx0,degree)
quad = gt.default_quadrature(spx1,degree)
quad = gt.default_quadrature(spx2,degree)
quad = gt.default_quadrature(spx3,degree)

quad = gt.default_quadrature(cube0,degree)
quad = gt.default_quadrature(cube1,degree)
quad = gt.default_quadrature(cube2,degree)
quad = gt.default_quadrature(cube3,degree)


order = 1
fe = gt.lagrangian_fe(spx0,order)
fe = gt.lagrangian_fe(spx1,order)
fe = gt.lagrangian_fe(spx2,order)
fe = gt.lagrangian_fe(spx3,order)
display(fe)

fe = gt.lagrangian_fe(cube0,order)
fe = gt.lagrangian_fe(cube1,order)
fe = gt.lagrangian_fe(cube2,order)
fe = gt.lagrangian_fe(cube3,order)
display(fe)

fe = gt.lagrangian_fe(cube0,order)
@show gt.monomial_exponents(fe)
@show gt.node_coordinates(fe)
fe = gt.lagrangian_fe(cube2,order)
@show gt.node_coordinates(fe)

spx2 = gt.unit_simplex(2)
quad = gt.default_quadrature(spx2,degree)
fe = gt.lagrangian_fe(spx2,order)
funs = gt.shape_functions(fe)
x = gt.coordinates(quad)
B = broadcast(gt.value,permutedims(funs),x)
display(B)
tabulator = gt.tabulator(fe)
A = tabulator(gt.value,x)
@test A≈B
x = gt.node_coordinates(fe)
A = tabulator(gt.value,x)

fe = gt.lagrangian_fe(spx2,order;shape=(3,))
funs = gt.shape_functions(fe)
x = gt.coordinates(quad)
B = broadcast(gt.value,permutedims(funs),x)
display(B)
tabulator = gt.tabulator(fe)
A = tabulator(gt.value,x)
@test A≈B
x = gt.node_coordinates(fe)
A = tabulator(gt.value,x)

fe = gt.lagrangian_fe(spx2,order;shape=())
funs = gt.shape_functions(fe)
x = gt.coordinates(quad)
B = broadcast(gt.value,permutedims(funs),x)
display(B)
tabulator = gt.tabulator(fe)
A = tabulator(gt.value,x)
@test A≈B
x = gt.node_coordinates(fe)
A = tabulator(gt.value,x)

outdir = mkpath(joinpath(@__DIR__,"..","output"))

msh =  joinpath(@__DIR__,"..","assets","demo.msh")
mesh = gt.mesh_from_gmsh(msh;complexify=false)
vtk_grid(joinpath(outdir,"demo"),gt.vtk_args(mesh)...) do vtk
    gt.vtk_physical_faces!(vtk,mesh)
end

@show gt.unit_simplex(0) |> gt.boundary
@show gt.unit_simplex(1) |> gt.boundary
@show gt.unit_simplex(2) |> gt.boundary
@show gt.unit_simplex(3) |> gt.boundary

@show gt.unit_n_cube(0) |> gt.boundary
@show gt.unit_n_cube(1) |> gt.boundary
@show gt.unit_n_cube(2) |> gt.boundary
@show gt.unit_n_cube(3) |> gt.boundary

order = 2
gt.lagrangian_fe(spx1,order) |> gt.boundary |> gt.topology

mesh = gt.mesh_from_gmsh(msh;complexify=false)

new_mesh, old_to_new = gt.complexify(mesh)

mesh = gt.mesh_from_gmsh(msh)

domain = (0,1,0,1)
cells = (2,2)
mesh = gt.cartesian_mesh(domain,cells)
mesh = gt.cartesian_mesh(domain,cells,boundary=false)
mesh = gt.cartesian_mesh(domain,cells,simplexify=true)
mesh = gt.cartesian_mesh(domain,cells,boundary=false,simplexify=true)

mesh = gt.mesh_from_gmsh(msh)
face_groups = gt.physical_faces(mesh)
group_names = gt.physical_names(mesh,2)
group_names = gt.physical_names(mesh)
group_names = gt.physical_names(mesh;merge_dims=true)
node_groups = gt.physical_nodes(mesh;merge_dims=true)
node_groups = gt.physical_nodes(mesh;merge_dims=true,disjoint=true)

vmesh, vglue = gt.visualization_mesh(mesh)

#∂spx0 = gt.boundary(spx0)
#∂spx0 = gt.boundary(spx1)
#∂cube0 = gt.boundary(cube0)



end # module

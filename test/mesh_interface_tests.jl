module MeshInterfaceTests

using Test
import GalerkinToolkit as gk
using WriteVTK

spx0 = gk.unit_simplex(0)
spx1 = gk.unit_simplex(1)
spx2 = gk.unit_simplex(2)
spx3 = gk.unit_simplex(3)
display(spx3)

cube0 = gk.unit_n_cube(0)
cube1 = gk.unit_n_cube(1)
cube2 = gk.unit_n_cube(2)
cube3 = gk.unit_n_cube(3)
display(cube3)

@show typeof(spx0)
@show typeof(cube0)
@show typeof(spx1)
@show typeof(cube1)

degree = 4
quad = gk.default_quadrature(spx0,degree)
quad = gk.default_quadrature(spx1,degree)
quad = gk.default_quadrature(spx2,degree)
quad = gk.default_quadrature(spx3,degree)

quad = gk.default_quadrature(cube0,degree)
quad = gk.default_quadrature(cube1,degree)
quad = gk.default_quadrature(cube2,degree)
quad = gk.default_quadrature(cube3,degree)


order = 1
fe = gk.lagrangian_fe(spx0,order)
fe = gk.lagrangian_fe(spx1,order)
fe = gk.lagrangian_fe(spx2,order)
fe = gk.lagrangian_fe(spx3,order)
display(fe)

fe = gk.lagrangian_fe(cube0,order)
fe = gk.lagrangian_fe(cube1,order)
fe = gk.lagrangian_fe(cube2,order)
fe = gk.lagrangian_fe(cube3,order)
display(fe)

fe = gk.lagrangian_fe(cube0,order)
@show gk.monomial_exponents(fe)
@show gk.node_coordinates(fe)
fe = gk.lagrangian_fe(cube2,order)
@show gk.node_coordinates(fe)

spx2 = gk.unit_simplex(2)
quad = gk.default_quadrature(spx2,degree)
fe = gk.lagrangian_fe(spx2,order)
funs = gk.shape_functions(fe)
x = gk.coordinates(quad)
B = broadcast(gk.value,permutedims(funs),x)
display(B)
tabulator = gk.tabulator(fe)
A = tabulator(gk.value,x)
@test A≈B
x = gk.node_coordinates(fe)
A = tabulator(gk.value,x)

fe = gk.lagrangian_fe(spx2,order;shape=(3,))
funs = gk.shape_functions(fe)
x = gk.coordinates(quad)
B = broadcast(gk.value,permutedims(funs),x)
display(B)
tabulator = gk.tabulator(fe)
A = tabulator(gk.value,x)
@test A≈B
x = gk.node_coordinates(fe)
A = tabulator(gk.value,x)

fe = gk.lagrangian_fe(spx2,order;shape=())
funs = gk.shape_functions(fe)
x = gk.coordinates(quad)
B = broadcast(gk.value,permutedims(funs),x)
display(B)
tabulator = gk.tabulator(fe)
A = tabulator(gk.value,x)
@test A≈B
x = gk.node_coordinates(fe)
A = tabulator(gk.value,x)

outdir = mkpath(joinpath(@__DIR__,"..","output"))

msh =  joinpath(@__DIR__,"..","assets","demo.msh")
mesh = gk.mesh_from_gmsh(msh;complexify=false)
vtk_grid(joinpath(outdir,"demo"),gk.vtk_args(mesh)...) do vtk
    gk.vtk_physical_faces!(vtk,mesh)
end

@show gk.unit_simplex(0) |> gk.boundary
@show gk.unit_simplex(1) |> gk.boundary
@show gk.unit_simplex(2) |> gk.boundary
@show gk.unit_simplex(3) |> gk.boundary

@show gk.unit_n_cube(0) |> gk.boundary
@show gk.unit_n_cube(1) |> gk.boundary
@show gk.unit_n_cube(2) |> gk.boundary
@show gk.unit_n_cube(3) |> gk.boundary

order = 2
gk.lagrangian_fe(spx1,order) |> gk.boundary |> gk.topology

mesh = gk.mesh_from_gmsh(msh;complexify=false)

new_mesh, old_to_new = gk.complexify(mesh)

mesh = gk.mesh_from_gmsh(msh)

domain = (0,1,0,1)
cells = (2,2)
mesh = gk.cartesian_mesh(domain,cells)
mesh = gk.cartesian_mesh(domain,cells,boundary=false)
mesh = gk.cartesian_mesh(domain,cells,simplexify=true)
mesh = gk.cartesian_mesh(domain,cells,boundary=false,simplexify=true)

mesh = gk.mesh_from_gmsh(msh)
face_groups = gk.physical_faces(mesh)
group_names = gk.physical_names(mesh,2)
group_names = gk.physical_names(mesh)
group_names = gk.physical_names(mesh;merge_dims=true)
node_groups = gk.physical_nodes(mesh;merge_dims=true)
node_groups = gk.physical_nodes(mesh;merge_dims=true,disjoint=true)

vmesh, vglue = gk.visualization_mesh(mesh)

#∂spx0 = gk.boundary(spx0)
#∂spx0 = gk.boundary(spx1)
#∂cube0 = gk.boundary(cube0)



end # module

module MeshInterfaceTests

using Test
import GalerkinToolkit as gk
using WriteVTK
using PartitionedArrays
using Metis

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
fe = gk.lagrange_mesh_face(spx0,order)
fe = gk.lagrange_mesh_face(spx1,order)
fe = gk.lagrange_mesh_face(spx2,order)
fe = gk.lagrange_mesh_face(spx3,order)
display(fe)

fe = gk.lagrange_mesh_face(cube0,order)
fe = gk.lagrange_mesh_face(cube1,order)
fe = gk.lagrange_mesh_face(cube2,order)
fe = gk.lagrange_mesh_face(cube3,order)
display(fe)

fe = gk.lagrange_mesh_face(cube0,order)
@show gk.monomial_exponents(fe)
@show gk.node_coordinates(fe)
fe = gk.lagrange_mesh_face(cube2,order)
@show gk.node_coordinates(fe)

spx2 = gk.unit_simplex(2)
quad = gk.default_quadrature(spx2,degree)
fe = gk.lagrange_mesh_face(spx2,order)
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
    gk.vtk_physical_nodes!(vtk,mesh)
end
for d in 0:gk.num_dims(mesh)
    vtk_grid(joinpath(outdir,"demo_$d"),gk.vtk_args(mesh,d)...) do vtk
        gk.vtk_physical_faces!(vtk,mesh,d)
        gk.vtk_physical_nodes!(vtk,mesh,d)
    end
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
gk.lagrange_mesh_face(spx1,order) |> gk.boundary |> gk.topology

mesh = gk.mesh_from_gmsh(msh;complexify=false)

new_mesh, old_to_new = gk.complexify(mesh)

mesh = gk.mesh_from_gmsh(msh)

domain = (0,1,0,1)
cells = (2,2)
mesh = gk.cartesian_mesh(domain,cells)
mesh = gk.cartesian_mesh(domain,cells,boundary=false)
mesh = gk.cartesian_mesh(domain,cells,simplexify=true)
mesh = gk.cartesian_mesh(domain,cells,boundary=false,simplexify=true)

mesh = gk.cartesian_mesh(domain,cells)
vtk_grid(joinpath(outdir,"cartesian"),gk.vtk_args(mesh)...) do vtk
    gk.vtk_physical_faces!(vtk,mesh)
    gk.vtk_physical_nodes!(vtk,mesh)
end
for d in 0:gk.num_dims(mesh)
    vtk_grid(joinpath(outdir,"cartesian_$d"),gk.vtk_args(mesh,d)...) do vtk
        gk.vtk_physical_faces!(vtk,mesh,d)
        gk.vtk_physical_nodes!(vtk,mesh,d)
    end
end

mesh = gk.mesh_from_gmsh(msh)
face_groups = gk.physical_faces(mesh)
group_names = gk.physical_names(mesh,2)
group_names = gk.physical_names(mesh)
group_names = gk.physical_names(mesh;merge_dims=true)
node_groups = gk.physical_nodes(mesh;merge_dims=true)
node_groups = gk.physical_nodes(mesh;merge_dims=true,disjoint=true)

vmesh, vglue = gk.visualization_mesh(mesh)

np = 4
ranks = DebugArray(LinearIndices((np,)))
mesh = gk.mesh_from_gmsh(msh)
pmesh = gk.partition_mesh(Metis.partition,ranks,mesh,via=:nodes)
function setup(mesh,ids,rank)
    face_to_owner = zeros(Int,sum(gk.num_faces(mesh)))
    D = gk.num_dims(mesh)
    for d in 0:D
        face_to_owner[gk.face_range(mesh,d)] = local_to_owner(gk.face_indices(ids,d))
    end
    pvtk_grid(joinpath(outdir,"pmesh"),gk.vtk_args(mesh)...;part=rank,nparts=np) do vtk
        gk.vtk_physical_faces!(vtk,mesh)
        gk.vtk_physical_nodes!(vtk,mesh)
        vtk["piece"] = fill(rank,sum(gk.num_faces(mesh)))
        vtk["owner"] = local_to_owner(gk.node_indices(ids))
        vtk["owner"] = face_to_owner
    end
end
map(setup,partition(pmesh),gk.index_partition(pmesh),ranks)

end # module

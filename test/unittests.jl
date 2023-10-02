module UnitTests

import GalerkinToolkit as glk
using StaticArrays
using WriteVTK
using SparseArrays
using LinearAlgebra
using ForwardDiff
using Test
using PartitionedArrays
using Metis


# TODO
# before more cleanup:
# generate a high order mesh from a linear non oriented mesh
#
# visualization mesh [done]
#
# remove interpolation [done]
# remove gradient! and value! [done]
# swap order in tabulation matrix. Important!!
# topology reference_faces[end][1].boundary
# cleanups
# parametrize with Ti and Tf
# think about order if it needs to be a tuple

#lagrangian_reference_face(geo)
#lagrangian_reference_element(geo)
#lagrangian_reference_element(geo,shape=:scaler)
#lagrangian_reference_element(geo,shape=())
#lagrangian_reference_element(geo,shape=(2,))

#segment
#num_dims(segment)
#face_reference_id(segment,d)
#face_reference_id(segment.boundary,d)
#face_reference_id(segment.boundary.topology,d)
#shape_functions(segment,d)
#node_coordinates(segment,d)
#
#mesh
#num_dims(mesh)
#face_reference_id(mesh,d)
#face_reference_id(mesh.topology,d)
#physical_groups(mesh)
#
#mesh2 = set_data(mesh;physical_groups,topology)

outdir = mkpath(joinpath(@__DIR__,"..","output"))

msh =  joinpath(@__DIR__,"..","assets","demo.msh")
mesh = glk.mesh_from_gmsh(msh;complexify=false)

np = 4
ranks = DebugArray(LinearIndices((np,)))
# TODO inclure global "glue" outside pmesh?
pmesh = glk.partition_mesh(Metis.partition,ranks,mesh,via=:nodes)
pmesh = glk.partition_mesh(Metis.partition,ranks,mesh,via=:cells)

map(pmesh,ranks) do mesh,rank
    pvtk_grid(joinpath(outdir,"pdebug"),glk.vtk_args(mesh)...;part=rank,nparts=np) do vtk
        glk.vtk_physical_groups!(vtk,mesh)
        vtk["piece"] = fill(rank,sum(glk.num_faces(mesh)))
        vtk["owner"] = local_to_owner(mesh.local_nodes)
        vtk["interface"] = map(colors->Int(length(colors)!=1),mesh.local_node_colors)
    end
end

domain = (1,2,1,2)
cells = (2,2)
mesh = glk.cartesian_mesh(domain,cells)

D = glk.num_dims(mesh)
geo = glk.reference_faces(mesh,D)|>first|>glk.geometry
reffe = glk.lagrangian_reference_element(geo,order=3)

local_dofs = glk.Object(;
    num_dims=Val(2),
    face_reference_id = [Int[],Int[],fill(1,glk.num_faces(mesh,2))],
    reference_faces = ((),(),[reffe])
)

dof_glue = glk.dof_glue_from_mesh_and_local_dofs(mesh,local_dofs)

geo = glk.unit_n_cube(Val(1))
reffe = glk.lagrangian_reference_element(geo)
A = reffe.shape_functions.tabulation_matrix(glk.value,reffe.node_coordinates)
B = reffe.shape_functions.tabulation_matrix(ForwardDiff.gradient,reffe.node_coordinates)

reffe = glk.lagrangian_reference_element(geo,shape=(2,),order=2,major=:node)
A = reffe.shape_functions.tabulation_matrix(glk.value,reffe.node_coordinates)
A = reffe.shape_functions.tabulation_matrix(ForwardDiff.jacobian,reffe.node_coordinates)

reffe = glk.lagrangian_reference_element(geo,shape=(2,),major=:component)
A = reffe.shape_functions.tabulation_matrix(glk.value,reffe.node_coordinates)
A = reffe.shape_functions.tabulation_matrix(ForwardDiff.jacobian,reffe.node_coordinates)
display(reffe.node_to_dofs)

reffe = glk.lagrangian_reference_element(geo,shape=(),major=:component)
A = reffe.shape_functions.tabulation_matrix(glk.value,reffe.node_coordinates)
A = reffe.shape_functions.tabulation_matrix(ForwardDiff.jacobian,reffe.node_coordinates)
display(reffe.node_to_dofs)

geo = glk.unit_n_cube(Val(2))
reffe = glk.lagrangian_reference_element(geo,order=3,shape=(2,))
display(reffe.face_own_dofs)
display(reffe.face_own_dof_permutations)

reffe = glk.lagrangian_reference_element(geo,shape=(2,3))
A = reffe.shape_functions.tabulation_matrix(glk.value,reffe.node_coordinates)

geo = glk.unit_simplex(Val(1))
perms = glk.vertex_permutations_from_geometry(geo)
@test length(perms) == 2
glk.vertex_permutations(geo) == perms

geo = glk.unit_simplex(Val(3))
reffe = glk.lagrangian_reference_face(geo)

@show reffe
display(reffe)

perms = glk.vertex_permutations_from_geometry(geo)
glk.vertex_permutations(geo) == perms

geo = glk.unit_n_cube(Val(2))
perms = glk.vertex_permutations_from_geometry(geo)
@test length(perms) == 8
glk.vertex_permutations(geo) == perms

quad1 = glk.lagrangian_reference_face(geo)
node_perms = glk.interior_node_permutations_from_reference_face(quad1)
@test all(map(i->all(i .!= 0),node_perms))

quad2 = glk.lagrangian_reference_face(geo,order=2)
node_perms = glk.interior_node_permutations_from_reference_face(quad2)
@test all(map(i->all(i .!= 0),node_perms))

quad3 = glk.lagrangian_reference_face(geo,order=3)
node_perms = glk.interior_node_permutations_from_reference_face(quad3)
@test all(map(i->all(i .!= 0),node_perms))

domain = (1,2,1,2,1,2)
cells = (2,2,2)
mesh = glk.structured_simplex_mesh_with_boundary(domain,cells)
vtk_grid(joinpath(outdir,"debug"),glk.vtk_args(mesh)...) do vtk
    glk.vtk_physical_groups!(vtk,mesh)
end

tet = glk.unit_simplex(Val(3))
refface = glk.lagrangian_reference_face(tet,order=2)
mesh = glk.mesh_from_reference_face(refface)
vtk_grid(joinpath(outdir,"debug"),glk.vtk_args(mesh)...) do vtk
    vtk["nodeid"] = 1:6
    glk.vtk_physical_groups!(vtk,mesh)
end

vtk_grid(joinpath(outdir,"debug_boundary"),glk.vtk_args(glk.boundary(refface))...) do vtk 
    vtk["nodeid"] = 1:6
end

hex = glk.unit_n_cube(Val(3))
tri_hex = glk.simplexify_reference_geometry(hex)
vtk_grid(joinpath(outdir,"debug"),glk.vtk_args(tri_hex)...) do vtk
    glk.vtk_physical_groups!(vtk,tri_hex)
end

order = 2
refface = glk.lagrangian_reference_face(hex;order)
mesh = glk.simplexify_reference_face(refface)

vtk_grid(joinpath(outdir,"debug"),glk.vtk_args(mesh)...) |> vtk_save

domain = (1,2,1,2,1,2)
cells = (2,2,2)
mesh = glk.cartesian_mesh(domain,cells,boundary=false,complexify=false)

vtk_grid(joinpath(outdir,"cartesian"),glk.vtk_args(mesh)...) do vtk
    glk.vtk_physical_groups!(vtk,mesh)
end

mesh = glk.cartesian_mesh(domain,cells,complexify=false)

vtk_grid(joinpath(outdir,"cartesian"),glk.vtk_args(mesh)...) do vtk
    glk.vtk_physical_groups!(vtk,mesh)
end

domain = (1,2,1,2)
cells = (2,2)
mesh = glk.cartesian_mesh(domain,cells,simplexify=false,boundary=false,complexify=false)
vismesh, visglue = glk.visualization_mesh_from_mesh(mesh)
vismesh, visglue = glk.visualization_mesh_from_mesh(mesh,order=2)
vismesh, visglue = glk.visualization_mesh_from_mesh(mesh,resolution=4)

pointdata = zeros(glk.num_nodes(vismesh))
for face in 1:glk.num_faces(mesh,2)
    fine_nodes = visglue.face_fine_nodes[face]
    for fine_node in fine_nodes
        pointdata[fine_node] = face
    end
end

vtk_grid(joinpath(outdir,"vismesh"),glk.vtk_args(vismesh,2)...) do vtk
    vtk["cellid"] = visglue.parent_face
    vtk["pointdata"] = pointdata
end

mesh = glk.cartesian_mesh(domain,cells,simplexify=true,boundary=false,complexify=false)
vismesh, visglue = glk.visualization_mesh_from_mesh(mesh)
vismesh, visglue = glk.visualization_mesh_from_mesh(mesh,order=2)
vismesh, visglue = glk.visualization_mesh_from_mesh(mesh,resolution=4)

domain = (1,2,1,2,1,2)
cells = (2,2,2)
mesh = glk.cartesian_mesh(domain,cells,boundary=false,complexify=false)

mesh,_ = glk.complexify_mesh(mesh)

vtk_grid(joinpath(outdir,"cartesian"),glk.vtk_args(mesh)...) do vtk
    glk.vtk_physical_groups!(vtk,mesh)
end

mesh = glk.cartesian_mesh(domain,cells)

vtk_grid(joinpath(outdir,"cartesian"),glk.vtk_args(mesh)...) do vtk
    glk.vtk_physical_groups!(vtk,mesh)
end

mesh = glk.cartesian_mesh(domain,cells,simplexify=true)

vtk_grid(joinpath(outdir,"cartesian"),glk.vtk_args(mesh)...) do vtk
    glk.vtk_physical_groups!(vtk,mesh)
end

mesh = glk.cartesian_mesh(domain,cells,simplexify=true,boundary=false)

vtk_grid(joinpath(outdir,"cartesian"),glk.vtk_args(mesh)...) do vtk
    glk.vtk_physical_groups!(vtk,mesh)
end

vertex = glk.unit_n_cube(Val(0))
vertex1 = glk.lagrangian_reference_face(vertex)

segment = glk.unit_n_cube(Val(1))
vtk_grid(joinpath(outdir,"segment"),glk.vtk_args(segment.boundary)...) |> vtk_save

segment2 = glk.lagrangian_reference_face(segment)
segment3 = glk.lagrangian_reference_face(segment,order=2)
segment4 = glk.lagrangian_reference_face(segment,order=4)

quad = glk.unit_n_cube(Val(2))
vtk_grid(joinpath(outdir,"quad"),glk.vtk_args(quad.boundary)...) |> vtk_save

quad4 = glk.lagrangian_reference_face(quad,order=1)
quad9 = glk.lagrangian_reference_face(quad,order=2)
vtk_grid(joinpath(outdir,"quad9"),glk.vtk_args(quad9.boundary)...) |> vtk_save

mesh_quad4 = glk.mesh_from_reference_face(quad4)

vtk_grid(joinpath(outdir,"mesh_quad4"),glk.vtk_args(mesh_quad4)...) |> vtk_save

order = 1
segment = glk.unit_n_cube(Val(1))
segment2 = glk.lagrangian_reference_face(segment;order)
quad = glk.unit_n_cube(Val(2))
quad4 = glk.lagrangian_reference_face(quad;order)

node_coordinates = SVector{2,Float64}[(0,0),(1,0),(2,0),(0,1),(1,1),(2,1),(0,2),(1,2),(2,2)]
face_nodes = [
   Vector{Int}[],
   [[1,2],[2,3],[1,4],[4,7],[7,8],[8,9],[3,6],[6,9]],
   [[1,2,4,5],[2,3,5,6],[4,5,7,8],[5,6,8,9]]
  ]
face_reference_id = [Int[],Int[1,1,1,1,1,1,1,1],[1,1,1,1]]
reference_faces = ([],[segment2],[quad4])
physical_groups = [
  Dict([]),
  Dict(["face_1"=>[1,2],"face_2"=>[3,4],"boundary"=>[1,2,3,4,5,6,7,8]]),
  Dict(["domain"=>[1,2,3,4]])]
mesh = (;
    num_dims=Val(2),node_coordinates,
    face_nodes,face_reference_id,
    reference_faces,physical_groups)

d = 2
vtk_grid(joinpath(outdir,"mesh2"),glk.vtk_args(mesh,d)...) do vtk
    glk.vtk_physical_groups!(vtk,mesh,d)
end

d = 1
vtk_grid(joinpath(outdir,"mesh1"),glk.vtk_args(mesh,d)...) do vtk
    glk.vtk_physical_groups!(vtk,mesh,d)
end

vtk_grid(joinpath(outdir,"mesh"),glk.vtk_args(mesh)...) do vtk
    glk.vtk_physical_groups!(vtk,mesh)
end

new_mesh, old_to_new = glk.complexify_mesh(mesh)

vtk_grid(joinpath(outdir,"new_mesh"),glk.vtk_args(new_mesh)...) do vtk
    glk.vtk_physical_groups!(vtk,new_mesh)
end

msh =  joinpath(@__DIR__,"..","assets","demo.msh")
mesh = glk.mesh_from_gmsh(msh;complexify=true)
@test glk.periodic_nodes(mesh) == ([] => [])

vtk_grid(joinpath(outdir,"mesh"),glk.vtk_args(mesh)...) do vtk
    glk.vtk_physical_groups!(vtk,mesh)
end

for d in 1:3
    vtk_grid(joinpath(outdir,"mesh_$d"),glk.vtk_args(mesh,d)...) do vtk
        glk.vtk_physical_groups!(vtk,mesh,d)
    end
end

end # module

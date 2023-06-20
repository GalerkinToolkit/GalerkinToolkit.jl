module TestGmsh

using GalerkinToolkit
using WriteVTK
using Test

# TODO: compare generated files with a reference

function write_mesh(mesh,filebase)
    for d in 0:dimension(mesh)
        fn = "$(filebase)_$(d)"
        vtk_grid(fn,vtk_args(mesh,d)...) do vtk
            vtk_physical_groups!(vtk,mesh,d)
        end
    end
end

file = joinpath(@__DIR__,"..","assets","twoTetraeder.msh")
mesh = gmsh_mesh(file)
write_mesh(mesh,"twoTetraeder")

file = joinpath(@__DIR__,"..","assets","t1.msh")
mesh = gmsh_mesh(file)
write_mesh(mesh,"t1")

file = joinpath(@__DIR__,"..","assets","periodic.msh")
mesh = gmsh_mesh(file)
constraints = periodic_node_constraints(mesh)
display(constraints)
@test length(periodic_nodes(constraints)) != 0
write_mesh(mesh,"periodic")

file = joinpath(@__DIR__,"..","assets","higher_order_2D.msh")
mesh = gmsh_mesh(file)
write_mesh(mesh,"higher_order_2D")

complex, new_faces = face_complex(mesh)
vtk_grid("complex",vtk_args(complex)...) do vtk
    #vtk_physical_groups!(vtk,complex)
end

vtk_grid("higher_order_2D",vtk_args(mesh)...) do vtk
    vtk_physical_groups!(vtk,mesh)
end


file = joinpath(@__DIR__,"..","assets","solid.msh")
mesh = gmsh_mesh(file)

vtk_grid("solid",vtk_args(mesh)...) do vtk
    vtk_physical_groups!(vtk,mesh)
end

D = dimension(mesh)
domain_tags = ["material_2"]
dirichlet_tags = ["surface_2"]
neumann_tags = ["surface_1"]

D = dimension(mesh)-1
domain_tags = ["surface_2"]
dirichlet_tags = ["surface_2_c"]
neumann_tags = String[]

node_to_tag = classify_mesh_nodes(mesh,dirichlet_tags,D)
cell_to_tag = classify_mesh_faces(mesh,domain_tags,D)
face_to_tag = classify_mesh_faces(mesh,neumann_tags,D-1)

vtk_grid("solid_2",vtk_args(mesh,D-1)...) do vtk
    vtk_physical_groups!(vtk,mesh,D-1)
    vtk["neumann_tag",WriteVTK.VTKCellData()] = face_to_tag
end

nnodes = num_nodes(mesh)
cell_to_nodes = face_nodes(mesh,D)
cells_in_domain = collect(Int32,findall(i->i>0,cell_to_tag))

domain_cell_to_dofs, dofs_and_nondofs = restrict_face_dofs(nnodes,cell_to_nodes,cells_in_domain)

faces_in_neumann = collect(Int32,findall(i->i>0,face_to_tag))

dof_to_tag = view(node_to_tag,first(dofs_and_nondofs))
free_and_dirichlet = partition_from_mask(i->i==0,dof_to_tag)

n_domain_cells = length(cells_in_domain)

u(x) = sum(x)
tag_to_u = fill(x->1.,length(dirichlet_tags))

free_values = zeros(length(first(free_and_dirichlet)))
dirichlet_values = zeros(length(last(free_and_dirichlet)))

dof_permutation = permutation(free_and_dirichlet)

dof_to_node = first(dofs_and_nondofs)
node_to_coordinates = node_coordinates(mesh)
n_free_dofs = length(first(free_and_dirichlet))
for domain_cell in 1:n_domain_cells
    dofs = domain_cell_to_dofs[domain_cell]
    for dof in dofs
        node = dof_to_node[dof]
        x = node_to_coordinates[node]
        p = dof_permutation[dof]
        if p > n_free_dofs
            tag = dof_to_tag[dof]
            dirichlet_values[p-n_free_dofs] = tag_to_u[tag](x)
        else
            free_values[p] = u(x)
        end
    end
end

node_values = zeros(nnodes)
node_values[view(dof_to_node,first(free_and_dirichlet))] = free_values
node_values[view(dof_to_node,last(free_and_dirichlet))] = dirichlet_values

vtk_grid("solid_3",vtk_args(mesh,D)...) do vtk
    vtk_physical_groups!(vtk,mesh,D)
    vtk["domain_tag",WriteVTK.VTKCellData()] = cell_to_tag
    vtk["dirichlet_tag",WriteVTK.VTKPointData()] = node_to_tag
    vtk["values",WriteVTK.VTKPointData()] = node_values
end

end # module

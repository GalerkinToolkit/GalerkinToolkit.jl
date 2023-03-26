module IntegrationTests

using GalerkinToolkit
using WriteVTK

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

#nnodes = num_nodes(mesh)
#node_to_tag = zeros(Int32,nnodes)
#for d in D:-1:0
#    face_to_nodes = face_nodes(mesh,d)
#    face_groups = physical_groups(mesh,d)
#    for (tag,name) in enumerate(dirichlet_tags)
#        for (name2,faces) in face_groups
#            if name != name2
#                continue
#            end
#            for face in faces
#                nodes = face_to_nodes[face]
#                node_to_tag[nodes] .= tag
#            end
#        end
#    end
#end

cell_to_tag = classify_mesh_faces(mesh,domain_tags,D)

#ncells = num_faces(mesh,D)
#cell_to_tag = zeros(Int32,ncells)
#cell_groups = physical_groups(mesh,D)
#for (tag,name) in enumerate(domain_tags)
#    for (name2,cells) in cell_groups
#        if name != name2
#            continue
#        end
#        cell_to_tag[cells] .= tag
#    end
#end

face_to_tag = classify_mesh_faces(mesh,neumann_tags,D-1)

#nfaces = num_faces(mesh,D-1)
#face_to_tag = zeros(Int32,nfaces)
#face_groups = physical_groups(mesh,D-1)
#for (tag,name) in enumerate(neumann_tags)
#    for (name2,faces) in face_groups
#        if name != name2
#            continue
#        end
#        face_to_tag[faces] .= tag
#    end
#end

vtk_grid("solid_2",vtk_args(mesh,D-1)...) do vtk
    vtk_physical_groups!(vtk,mesh,D-1)
    vtk["neumann_tag",WriteVTK.VTKCellData()] = face_to_tag
end

cells_in_domain = collect(Int32,findall(i->i>0,cell_to_tag))
cell_to_nodes = face_nodes(mesh,D)
nnodes = num_nodes(mesh)
node_to_touched = fill(false,nnodes)
for nodes in view(cell_to_nodes,cells_in_domain)
    node_to_touched[nodes] .= true
end

dofs_and_nondofs = partition_from_mask(node_to_touched)
domain_cell_to_dofs = JaggedArray{Int32,Int32}(view(cell_to_nodes,cells_in_domain))
node_permutation = permutation(dofs_and_nondofs)
domain_cell_to_dofs.data[:] = view(node_permutation,domain_cell_to_dofs.data)

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

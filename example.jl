module TMP

import GalerkinToolkit as gk
using WriteVTK

function test_two_level_mesh_with_nonperiodic_square_unit_cell()
    # tests a 4 x 4 unit cell in a 1 x 1 coarse mesh
    domain = (0,10,0,10)
    cells = (1,1)
    coarse_mesh = gk.cartesian_mesh(domain,cells)
    unit_cell_mesh = gk.cartesian_mesh((0,1,0,1), (4, 4)) # non periodic example 
    final_mesh, final_glue = gk.two_level_mesh(coarse_mesh, unit_cell_mesh)

    # periodicity info -- shouldn't apply  
    periodic_nodes = gk.periodic_nodes(unit_cell_mesh)
    fine_pnode_to_fine_node = periodic_nodes.first 
    fine_pnode_to_master_fine_node = periodic_nodes.second 
    node_ids = collect(1:gk.num_nodes(unit_cell_mesh))
    fine_node_to_master_fine_node = copy(node_ids)
    fine_node_to_master_fine_node[
        fine_pnode_to_fine_node] = fine_pnode_to_master_fine_node 

    # visualize unit cell 
    vtk_grid(joinpath("output","20240305-square-gmsh"),gk.vtk_args(unit_cell_mesh)...) do vtk
        gk.vtk_physical_faces!(vtk,unit_cell_mesh)
        gk.vtk_physical_nodes!(vtk,unit_cell_mesh)
        vtk["periodic_master_id"] = fine_node_to_master_fine_node
        vtk["node_id"] = node_ids
    end

    # visualize final mesh
    vtk_grid(joinpath("output","20240305-final-gmsh"),gk.vtk_args(final_mesh)...) do vtk
        gk.vtk_physical_faces!(vtk,final_mesh)
        gk.vtk_physical_nodes!(vtk,final_mesh)
    end
end 

function test_two_level_mesh_with_square_periodic_unit_cell()

    # Load periodic fine (unit cell) mesh and get periodic info
    unit_cell_mesh_fpath = joinpath(
        @__DIR__, "assets", "coarse_periodic_right_left_top_bottom.msh")
    unit_cell_mesh = gk.mesh_from_gmsh(unit_cell_mesh_fpath)
    periodic_nodes = gk.periodic_nodes(unit_cell_mesh)
    fine_pnode_to_fine_node = periodic_nodes.first 
    fine_pnode_to_master_fine_node = periodic_nodes.second 

    # make two level mesh using periodic unit cell 
    domain = (0,10,0,10)
    cells = (1,1)
    coarse_mesh = gk.cartesian_mesh(domain,cells)
    periodic_final_mesh, periodic_final_glue = gk.two_level_mesh(
        coarse_mesh, unit_cell_mesh)

    # labeling periodic nodes 
    node_ids = collect(1:gk.num_nodes(unit_cell_mesh))
    fine_node_to_master_fine_node = copy(node_ids)
    fine_node_to_master_fine_node[
        fine_pnode_to_fine_node] = fine_pnode_to_master_fine_node 

    # Fixing vertices indirection
    # assumes one level of indirection for periodicity (e.g., vertex -> master -> master)
    n_fine_nodes = length(node_ids)
    for fine_node in 1:n_fine_nodes
        master_fine_node = fine_node_to_master_fine_node[fine_node]
        if fine_node == master_fine_node
            continue  
        end

        master_master_fine_node = fine_node_to_master_fine_node[master_fine_node]
        fine_node_to_master_fine_node[fine_node] = master_master_fine_node 
    end 

    # visualize unit cell 
    vtk_grid(
        joinpath("output","20240305-periodic-square-gmsh"),
        gk.vtk_args(unit_cell_mesh)...) do vtk
        gk.vtk_physical_faces!(vtk,unit_cell_mesh)
        gk.vtk_physical_nodes!(vtk,unit_cell_mesh)
        vtk["periodic_master_id"] = fine_node_to_master_fine_node
        vtk["node_id"] = node_ids
    end

    # visualize final mesh
    vtk_grid(
        joinpath("output","20240305-periodic-final-gmsh"),
        gk.vtk_args(periodic_final_mesh)...) do vtk
        gk.vtk_physical_faces!(vtk,periodic_final_mesh)
        gk.vtk_physical_nodes!(vtk,periodic_final_mesh)
    end
end

end # module TMP

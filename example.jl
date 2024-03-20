"""
TODO:
3D cartesian and naming and gmsh

geometry.jl
geometry_test.jl # clean up example.jl
descriptive names for assets and geo files
2D square-<square, triangl>-4x4
3D  

for testing correctness of outputs:
   given a julia object, store the integer for correct object contents (find this func, not ash)   
    take some example node ids and corresponding coordinates of interface nodes against
        the hardcoded physical coordinates that would be expected for those points

mesh generation should be of order of multigrid solver example02

loading of sequential total mesh 

parallel mesh support in gmsh?

Go through assets/ and remove unneeded geos/meshes 

Test using puzzle piece mesh??
"""
module TMP

import GalerkinToolkit as gk
using WriteVTK

function test_two_level_mesh_with_nonperiodic_square_unit_cell()
    # tests a 4 x 4 unit cell in a 1 x 1 coarse mesh
    coarse_domain = (0,10,0,10)
    coarse_mesh_dims = (1,1)
    coarse_mesh = gk.cartesian_mesh(coarse_domain,coarse_mesh_dims)

    unit_cell_domain = (0,1,0,1)
    unit_cell_dims = (4,4)
    unit_cell_mesh = gk.cartesian_mesh(unit_cell_domain, unit_cell_dims) 
    final_mesh, final_glue = gk.two_level_mesh(coarse_mesh, unit_cell_mesh)

    # visualize the glk unit cell 
    visualize_unit_cell_mesh(
        unit_cell_mesh, joinpath(
            "output", "nonperiodic-square-4x4-unit-cell"))

    # visualize final mesh with 1x1 coarse cell
    vtk_grid(joinpath(
        "output",
        "final-mesh-with-nonperiodic-square-unit-cell-and-coarse-1x1-mesh"),
        gk.vtk_args(final_mesh)...) do vtk
        gk.vtk_physical_faces!(vtk,final_mesh)
        gk.vtk_physical_nodes!(vtk,final_mesh)
    end

    # tests a 4 x 4 unit cell in a 2 x 2 coarse mesh
    coarse_domain = (0,10,0,10)
    coarse_mesh_dims = (2,2)
    coarse_mesh = gk.cartesian_mesh(coarse_domain,coarse_mesh_dims)
    final_mesh, final_glue = gk.two_level_mesh(coarse_mesh, unit_cell_mesh)
 
    # visualize final mesh with 2x2 coarse cell
    vtk_grid(joinpath(
        "output",
        "final-mesh-with-nonperiodic-square-unit-cell-and-2x2-coarse-mesh"),
        gk.vtk_args(final_mesh)...) do vtk
        gk.vtk_physical_faces!(vtk,final_mesh)
        gk.vtk_physical_nodes!(vtk,final_mesh)
    end
end 

function test_two_level_mesh_with_periodic_square_unit_cell()

    # Load periodic fine (unit cell) mesh and get periodic info
    unit_cell_mesh_fpath = joinpath(
        @__DIR__, "assets", "coarse_periodic_right_left_top_bottom_triangle.msh")
    unit_cell_mesh = gk.mesh_from_gmsh(unit_cell_mesh_fpath)

    # visualize the periodic gmsh unit cell 
    visualize_unit_cell_mesh(unit_cell_mesh, joinpath(
        "output", "periodic-square-4x4-unit-cell-gmsh"))

    # visualize final mesh with 2x2 coarse mesh and 4x4 unit cell
    periodic_final_mesh, periodic_final_glue = gk.two_level_mesh(
        coarse_mesh, unit_cell_mesh)
    n_nodes = gk.num_nodes(periodic_final_mesh)
    vtk_grid(
        joinpath("output",
        "final-mesh-with-periodic-square-4x4-unit-cell-gmsh-and-2x2-coarse-mesh"),
        gk.vtk_args(periodic_final_mesh)...) do vtk
        gk.vtk_physical_faces!(vtk,periodic_final_mesh)
        gk.vtk_physical_nodes!(vtk,periodic_final_mesh)
        vtk["node_ids"] = collect(1:n_nodes)
    end
end

function visualize_unit_cell_mesh(unit_cell_mesh, outpath)

    # Get periodic node info about unit cell
    periodic_nodes = gk.periodic_nodes(unit_cell_mesh)
    fine_pnode_to_fine_node = periodic_nodes.first 
    fine_pnode_to_master_fine_node = periodic_nodes.second 
    
    # Labeling periodic nodes 
    node_ids = collect(1:gk.num_nodes(unit_cell_mesh))
    fine_node_to_master_fine_node = copy(node_ids)
    fine_node_to_master_fine_node[
        fine_pnode_to_fine_node] = fine_pnode_to_master_fine_node

    # Handles periodic vertices whose master node is the master of another periodic 
    # vertex assuming one level of indirection for periodicity 
    # (e.g., vertex -> master -> master)
    n_fine_nodes = length(node_ids)
    for fine_node in 1:n_fine_nodes
        master_fine_node = fine_node_to_master_fine_node[fine_node]
        if fine_node == master_fine_node
            continue  
        end

        master_master_fine_node = fine_node_to_master_fine_node[master_fine_node]
        fine_node_to_master_fine_node[fine_node] = master_master_fine_node 
    end 

    # Visualize in paraview
    vtk_grid(
        outpath,
        gk.vtk_args(unit_cell_mesh)...) do vtk
        gk.vtk_physical_faces!(vtk,unit_cell_mesh)
        gk.vtk_physical_nodes!(vtk,unit_cell_mesh)
        vtk["periodic_master_id"] = fine_node_to_master_fine_node
        vtk["node_id"] = node_ids
    end
end 

TMP.test_two_level_mesh_with_periodic_square_unit_cell()

end # module TMP

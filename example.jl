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

    # periodicity info -- shouldn't apply  
    periodic_nodes = gk.periodic_nodes(unit_cell_mesh)
    fine_pnode_to_fine_node = periodic_nodes.first 
    fine_pnode_to_master_fine_node = periodic_nodes.second 
    node_ids = collect(1:gk.num_nodes(unit_cell_mesh))
    fine_node_to_master_fine_node = copy(node_ids)
    fine_node_to_master_fine_node[
        fine_pnode_to_fine_node] = fine_pnode_to_master_fine_node 

    # visualize unit cell 
    vtk_grid(joinpath(
        "output",
        "nonperiodic-square-unit-cell"),gk.vtk_args(unit_cell_mesh)...) do vtk
        gk.vtk_physical_faces!(vtk,unit_cell_mesh)
        gk.vtk_physical_nodes!(vtk,unit_cell_mesh)
        vtk["periodic_master_id"] = fine_node_to_master_fine_node
        vtk["node_id"] = node_ids
    end

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
    # TODO: 
    unit_cell_mesh_fpath_4x4 = joinpath(
        @__DIR__, "assets", "coarse_periodic_right_left_top_bottom_triangle.msh")
    unit_cell_mesh_4x4 = gk.mesh_from_gmsh(unit_cell_mesh_fpath_4x4)

    unit_cell_mesh_fpath_2x2 = joinpath(
        @__DIR__, "assets", "coarse_periodic_right_left_top_bottom_2x2.msh")
    unit_cell_mesh_2x2 = gk.mesh_from_gmsh(unit_cell_mesh_fpath_2x2)

    periodic_nodes_4x4 = gk.periodic_nodes(unit_cell_mesh_4x4)
    fine_pnode_to_fine_node_4x4 = periodic_nodes_4x4.first 
    fine_pnode_to_master_fine_node_4x4 = periodic_nodes_4x4.second 

    # TODO: visualize the coarse mesh in gmsh to see what the numbering
    # of nodes is
    # @assert enumeration of coarse faces and final nodes should match 

    ## NOTE: For debugging purposes, during the meeting i do comment out 
    # the visualization of the 1 x 1 coarse mesh
    # # coarse 1x1 mesh
    # coarse_domain = (0,10,0,10)
    # coarse_mesh_dims = (1,1)
    # coarse_mesh = gk.cartesian_mesh(coarse_domain,coarse_mesh_dims)

    # # make two level mesh using periodic unit cell and 1x1 coarse mesh
    # periodic_final_mesh, periodic_final_glue = gk.two_level_mesh(
    #     coarse_mesh, unit_cell_mesh)

    # # labeling periodic nodes 
     node_ids = collect(1:gk.num_nodes(unit_cell_mesh_4x4))
     fine_node_to_master_fine_node_4x4 = copy(node_ids)
     fine_node_to_master_fine_node_4x4[
          fine_pnode_to_fine_node_4x4] = fine_pnode_to_master_fine_node_4x4
    # 
    # ## Visualize unit cell 
    # # Must fix vertices indirection
    # # assumes one level of indirection for periodicity (e.g., vertex -> master -> master)
     n_fine_nodes = length(node_ids)
     for fine_node in 1:n_fine_nodes
          master_fine_node = fine_node_to_master_fine_node_4x4[fine_node]
          if fine_node == master_fine_node
              continue  
          end

          master_master_fine_node = fine_node_to_master_fine_node_4x4[master_fine_node]
          fine_node_to_master_fine_node_4x4[fine_node] = master_master_fine_node 
      end 

      @show vtk_grid(
          joinpath("output","periodic-square-4x4-unit-cell-gmsh"),
          gk.vtk_args(unit_cell_mesh_4x4)...) do vtk
          gk.vtk_physical_faces!(vtk,unit_cell_mesh_4x4)
          gk.vtk_physical_nodes!(vtk,unit_cell_mesh_4x4)
          vtk["periodic_master_id"] = fine_node_to_master_fine_node_4x4
          vtk["node_id"] = node_ids
      end

    # ## visualize final mesh with 1x1 coarse mesh
    # vtk_grid(
    #     joinpath("output",
    #     "final-mesh-with-periodic-square-unit-cell-gmsh-and-1x1-coarse-mesh"),
    #     gk.vtk_args(periodic_final_mesh)...) do vtk
    #     gk.vtk_physical_faces!(vtk,periodic_final_mesh)
    #     gk.vtk_physical_nodes!(vtk,periodic_final_mesh)
    # end

    ## coarse 2x2 mesh
    # should be equidistant
    coarse_domain = (0,4,0,4)
    coarse_mesh_dims = (2,2)
    coarse_mesh = gk.cartesian_mesh(coarse_domain,coarse_mesh_dims)
     
    ## make two level mesh using periodic unit cell and coarse 2x2 mesh
    periodic_final_mesh, periodic_final_glue = gk.two_level_mesh(
        coarse_mesh, unit_cell_mesh_2x2)

    ## visualize final mesh with 2x2 coarse mesh
    vtk_grid(
        joinpath("output",
        "final-mesh-with-periodic-square-2x2-unit-cell-gmsh-and-2x2-coarse-mesh"),
        gk.vtk_args(periodic_final_mesh)...) do vtk
        gk.vtk_physical_faces!(vtk,periodic_final_mesh)
        gk.vtk_physical_nodes!(vtk,periodic_final_mesh)
    end

    ## visualize final mesh with 2x2 coarse mesh and 4x4 unit cell
     
    periodic_final_mesh, periodic_final_glue = gk.two_level_mesh(
        coarse_mesh, unit_cell_mesh_4x4)
    n_nodes = gk.num_nodes(periodic_final_mesh)
    @show n_nodes
    vtk_grid(
        joinpath("output",
        "final-mesh-with-periodic-square-4x4-unit-cell-gmsh-and-2x2-coarse-mesh"),
        gk.vtk_args(periodic_final_mesh)...) do vtk
        gk.vtk_physical_faces!(vtk,periodic_final_mesh)
        gk.vtk_physical_nodes!(vtk,periodic_final_mesh)
        vtk["node_ids"] = collect(1:n_nodes)
    end
end

TMP.test_two_level_mesh_with_periodic_square_unit_cell()

end # module TMP

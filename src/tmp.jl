import GalerkinToolkit as gk
# TODO: load easier quad mesh
function two_level_mesh(outpath, coarse_mesh,fine_mesh;boundary_names=nothing)

    D = gk.num_dims(fine_mesh)
    n_fine_cells = gk.num_faces(fine_mesh,D)
    n_fine_nodes = gk.num_nodes(fine_mesh)
    fine_refcell = first(gk.reference_faces(fine_mesh,D))
    coarse_refcell = first(gk.reference_faces(coarse_mesh,D))
    n_coarse_cells = gk.num_faces(coarse_mesh,D)
    fine_cell_local_node_to_fine_node = gk.face_nodes(fine_mesh,D)
    coarse_cell_lnode_to_coarse_node = gk.face_nodes(coarse_mesh,D)
    coarse_node_to_x = gk.node_coordinates(coarse_mesh)
    topo = gk.topology(coarse_mesh)

    # The user can provide custom names for physical groups on the boundary
    # but we give some default value
    if boundary_names === nothing
        boundary_names = [
            [ "$d-face-$face" for face in 1:gk.num_faces(gk.boundary(coarse_refcell),d)] 
            for d in 0:(D-1)]
    end
    name_priority = reduce(vcat,boundary_names)
    fine_node_groups = gk.physical_nodes(fine_mesh;merge_dims=true,disjoint=true,name_priority)

    # Recover boundary info
    d_to_local_dface_to_fine_nodes = Vector{Vector{Vector{Int}}}(undef,D+1)
    fine_node_mask = fill(true,n_fine_nodes)
    for d in 0:(D-1)
        n_local_dfaces = gk.num_faces(gk.boundary(coarse_refcell),d)
        local_dface_to_fine_nodes = Vector{Vector{Int}}(undef,n_local_dfaces)
        for local_dface in 1:n_local_dfaces
            fine_nodes = fine_node_groups[boundary_names[d+1][local_dface]]
            local_dface_to_fine_nodes[local_dface] = fine_nodes
            fine_node_mask[fine_nodes] .= false
        end
        d_to_local_dface_to_fine_nodes[d+1] = local_dface_to_fine_nodes
    end
    d_to_local_dface_to_fine_nodes[D+1] = [findall(fine_node_mask)]

    ## Ensure consistent mapping of periodic nodes on opposite faces 
    # Get periodicity information
    fine_node_to_master_node = collect(1:n_fine_nodes)
    d_to_local_dface_to_permutation = Vector{Vector{Vector{Int}}}(undef, D+1)
    periodic_node_to_fine_node, periodic_node_to_master_node = gk.periodic_nodes(
        fine_mesh)
    fine_node_to_master_node[periodic_node_to_fine_node] = periodic_node_to_master_node

    # Fixing vertices indirection for unit cell
    # Assumes one level of indirection for periodicity (e.g., vertex -> master -> master)
    # TODO: d-1 levels of indirection
    for _ in 1:(D-1)
        for fine_node in 1:n_fine_nodes
            master_node = fine_node_to_master_node[fine_node]
            if fine_node == master_node
                continue  
            end
            master_master_node = fine_node_to_master_node[master_node]
            fine_node_to_master_node[fine_node] = master_master_node 
        end 
    end

    fine_node_to_permuted_node = zeros(Int, n_fine_nodes)
    fine_node_to_face_node = zeros(Int, n_fine_nodes)
    d_to_local_dface_to_opposite_dface = gk.opposite_faces(gk.geometry(coarse_refcell))
    d_to_local_dface_to_is_master = Vector{Vector{Bool}}(undef, D+1)
    for d in 0:(D-1)
        local_dface_to_is_master = similar(d_to_local_dface_to_opposite_dface[d+1], Bool) 
        n_local_dfaces = gk.num_faces(gk.boundary(coarse_refcell),d)
        for local_dface in 1:n_local_dfaces 
            local_dface_to_is_master[local_dface] = true
        end
     
        for local_dface in 1:n_local_dfaces 
            if  local_dface_to_is_master[local_dface] == false
                continue
            end
        
            opposite_local_dface = d_to_local_dface_to_opposite_dface[d+1][local_dface]
            local_dface_to_is_master[opposite_local_dface] = false
        end

        d_to_local_dface_to_is_master[d+1] = local_dface_to_is_master
    end 
    
    for d in 0:(D-1)
        n_local_dfaces = gk.num_faces(gk.boundary(coarse_refcell),d)
        local_dface_to_permutation = Vector{Vector{Int}}(undef, n_local_dfaces)
        local_dface_to_fine_nodes = d_to_local_dface_to_fine_nodes[d+1]
        for local_dface_1 in 1:n_local_dfaces

            fine_nodes_1 = local_dface_to_fine_nodes[local_dface_1]
            master_nodes_1 = fine_node_to_master_node[fine_nodes_1]
           
            # Handle nonperiodic case with identity permutation 
            # assumes consistent numbering of node ids on opposite faces... 
            # Also handles case where there are periodic nodes but the master nodes 
            # will use the identity permutation
            # TODO: in 3d geometry, since surfaces represent the periodic copies,
            # the opposite edges can end up being a map from slaves to slave, and therefore
            # this condition below fails because there is no obvious master (as in the 2d case) 
            if length(periodic_node_to_fine_node) == 0 || d_to_local_dface_to_is_master[d+1][local_dface_1] 

                permutation = collect(1:length(local_dface_to_fine_nodes[local_dface_1]))
                local_dface_to_permutation[local_dface_1] = permutation
                fine_node_to_permuted_node[fine_nodes_1] = permutation
                fine_node_to_face_node[fine_nodes_1] = 1:length(fine_nodes_1)
                continue 
            end

            # Handle periodic case: use a local reference cell and its opposite 
            # face to get the corresponding fine nodes and the master of those fine nodes.
            # Then ensure that a given local d-face has the same ordering as the opposite
            # face....
            local_dface_2 = d_to_local_dface_to_opposite_dface[d+1][local_dface_1]
            fine_nodes_2 = local_dface_to_fine_nodes[local_dface_2]
            master_nodes_2 = fine_node_to_master_node[fine_nodes_2]
            permutation = indexin(master_nodes_2, master_nodes_1)
            local_dface_to_permutation[local_dface_1] = permutation
            fine_node_to_permuted_node[fine_nodes_1] = permutation
            fine_node_to_face_node[fine_nodes_1] = 1:length(fine_nodes_1)
        end
        d_to_local_dface_to_permutation[d+1] = local_dface_to_permutation
    end
    d_to_local_dface_to_permutation[D+1] = [
        collect(1:length(d_to_local_dface_to_fine_nodes[D+1][1]))]

    
    node_ids = collect(1:length(fine_node_to_master_node))
    vtk_grid(
        outpath,
        gk.vtk_args(fine_mesh, D)...) do vtk
        gk.vtk_physical_faces!(vtk, fine_mesh, D)
        gk.vtk_physical_nodes!(vtk, fine_mesh, D)
        vtk["periodic_master_id"] = fine_node_to_master_node
        vtk["node_id"] = node_ids
        vtk["fine_node_to_permuted_node"] = fine_node_to_permuted_node 
        vtk["fine_node_to_face_node"] = fine_node_to_face_node
    end 
end 

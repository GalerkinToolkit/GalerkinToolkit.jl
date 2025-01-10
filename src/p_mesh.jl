
struct PMesh{A,B,C,D} <: GT.AbstractType
    mesh_partition::A
    node_partition::B
    face_partition::C
    partition_strategy::D
end
partition_strategy(a) = a.partition_strategy
PartitionedArrays.partition(m::PMesh) = m.mesh_partition
face_partition(a::PMesh,d) = a.face_partition[d+1]
node_partition(a::PMesh) = a.node_partition
function index_partition(a::PMesh)
    function setup(nodes,faces...)
        PMeshLocalIds(nodes,faces)
    end
    map(setup,a.node_partition,a.face_partition...)
end

function num_nodes(mesh::PMesh)
    length(PRange(mesh.node_partition))
end

function num_faces(mesh::PMesh)
    D = num_dims(mesh)
    [ num_faces(mesh,d) for d in 0:D]
end

function num_faces(mesh::PMesh,d)
    length(PRange(mesh.face_partition[d+1]))
end

function num_dims(mesh::PMesh)
    length(mesh.face_partition) - 1
end

function num_ambient_dims(mesh::PMesh)
    PartitionedArrays.getany(map(num_ambient_dims,mesh.mesh_partition))
end

function physical_names(pmesh::PMesh,d)
     map(pmesh.mesh_partition) do mesh
        GT.physical_names(mesh,d)
    end |> PartitionedArrays.getany
end

function physical_faces(pmesh::PMesh,d)
    names = physical_names(pmesh,d)
    map(collect(names)) do name
        name => map(pmesh.mesh_partition) do mesh
            GT.physical_faces(mesh,d)[name]
        end
    end |> Dict
end

function node_coordinates(pmesh::PMesh)
    data = map(GT.node_coordinates,pmesh.mesh_partition)
    PVector(data,pmesh.node_partition)
end

function face_nodes(mesh::PMesh)
    D = num_dims(mesh)
    [ face_nodes(mesh,d) for d in 0:D]
end

function face_nodes(pmesh::PMesh,d)
    data = map(mesh->GT.face_nodes(mesh,d),pmesh.mesh_partition)
    PVector(data,pmesh.face_partition[d+1])
end

function face_reference_id(mesh::PMesh)
    D = num_dims(mesh)
    [ face_reference_id(mesh,d) for d in 0:D]
end

function face_reference_id(pmesh::PMesh,d)
    data = map(mesh->GT.face_reference_id(mesh,d),pmesh.mesh_partition)
    PVector(data,pmesh.face_partition[d+1])
end

function reference_spaces(mesh::PMesh)
    D = num_dims(mesh)
    [ reference_spaces(mesh,d) for d in 0:D]
end

function reference_spaces(pmesh::PMesh,d)
    map(mesh->GT.reference_spaces(mesh,d),pmesh.mesh_partition)
end

function periodic_nodes(pmesh::PMesh)
    map(GT.periodic_nodes,pmesh.mesh_partition)
end

function outwards_normals(pmesh::PMesh)
    data = map(GT.outwards_normals,pmesh.mesh_partition)
    if eltype(data) <: Nothing
        nothing
    else
        data
    end
end

#function visualization_mesh(pmesh::PMesh,d,domface_to_face;kwargs...)
#    map(pmesh.mesh_partition,domface_to_face) do mesh, faces
#        visualization_mesh(mesh,d,faces;kwargs...)
#    end
#end

"""
    label_boundary_faces!(mesh::PMesh;physical_name="boundary")

Update `mesh` inplace by using partition ownership of faces to label only the boundary of 
meshes in a parallel mesh where the boundary is defined as a face owned by a given partition 
and not incident with any other faces on another partition (i.e., not on the interface). 
"""
function label_boundary_faces!(mesh::PMesh;physical_name="boundary")
    D = num_dims(mesh)
    d = D - 1
    face_parts = face_partition(mesh, d)
    v = pfill(1,face_parts)
    assemble!(v) |> wait
    map(partition(mesh),partition(v),face_parts) do mymesh, myv, myfaces
        topo = topology(mymesh)
        face_to_cells = face_incidence(topo,d,D)
        local_to_owner_face = local_to_owner(myfaces)
        part = part_id(myfaces)
        let face = 0
            faces = findall(face_to_cells) do cells
                face += 1
                myv[face] == 1 && local_to_owner_face[face] == part && length(cells) == 1
            end
            physical_faces(mymesh,d)[physical_name] = faces
        end
    end
    mesh
end

struct PMeshLocalIds{A,B} <: GT.AbstractType
    node_indices::A
    face_indices::B
end
node_indices(a::PMeshLocalIds) = a.node_indices
face_indices(a::PMeshLocalIds,d) = a.face_indices[d+1]

function partition_mesh(mesh,np;
    partition_strategy=GT.partition_strategy(),
    parts = LinearIndices((np,)),
    graph = mesh_graph(mesh;partition_strategy),
    graph_partition = Metis.partition(graph,np),
    renumber = true,
    )
    graph_nodes = partition_strategy.graph_nodes
    graph_edges = partition_strategy.graph_edges
    if graph_nodes === :nodes
        partition_mesh_nodes(graph_partition,parts,mesh,graph,partition_strategy,renumber)
    elseif graph_nodes === :cells
        partition_mesh_cells(graph_partition,parts,mesh,graph,partition_strategy,renumber)
    else
        error("Case not implemented")
    end
end

function scatter_mesh(pmeshes_on_main;source=MAIN)
    snd = map_main(pmeshes_on_main;main=source) do pmesh
        map(tuple,pmesh.mesh_partition,pmesh.node_partition,map(tuple,pmesh.face_partition...))
    end
    rcv = scatter(snd;source)
    mesh_partition, node_partition, face_partition_array = rcv |> tuple_of_arrays
    face_partition = face_partition_array |> tuple_of_arrays
    np = length(snd)
    snd2 = map_main(pmeshes_on_main;main=source) do pmesh
        fill(partition_strategy(pmesh),np)
    end
    rcv2 = scatter(snd2)
    partition_strategy_mesh = PartitionedArrays.getany(rcv2)
    PMesh(mesh_partition,node_partition,face_partition,partition_strategy_mesh)
end

function partition_mesh_nodes(node_to_color,parts,mesh,graph,partition_strategy,renumber)
    ghost_layers = partition_strategy.ghost_layers
    face_nodes_mesh = face_nodes(mesh)
    nnodes = num_nodes(mesh)
    function setup(part)
        onode_to_node = findall(color->color==part,node_to_color)
        node_to_mask = fill(false,nnodes)
        if ghost_layers == 0
            node_to_mask[onode_to_node] .= true
        elseif ghost_layers == 1
            for node in onode_to_node
                pini = graph.colptr[node]
                pend = graph.colptr[node+1]-1
                for p in pini:pend
                    node2 = graph.rowval[p]
                    node_to_mask[node2] = true
                end
            end
        else
            error("case not implemented")
        end
        lnode_to_node = findall(node_to_mask)
        lnode_to_color = node_to_color[lnode_to_node]
        onode_to_lnode = findall(color->color==part,lnode_to_color)
        hnode_to_lnode = findall(color->color!=part,lnode_to_color)
        onode_to_node = lnode_to_node[onode_to_lnode]
        hnode_to_node = lnode_to_node[hnode_to_lnode]
        hnode_to_color = lnode_to_color[hnode_to_lnode]
        own = OwnIndices(nnodes,part,onode_to_node)
        ghost = GhostIndices(nnodes,hnode_to_node,hnode_to_color)
        local_nodes = OwnAndGhostIndices(own,ghost,node_to_color)
        local_faces = map(face_nodes_mesh) do face_to_nodes
            nfaces = length(face_to_nodes)
            face_to_mask = fill(false,nfaces)
            for face in 1:nfaces
                nodes = face_to_nodes[face]
                if all(node->node_to_mask[node],nodes)
                    face_to_mask[face] = true
                end
            end
            lface_to_face = findall(face_to_mask)
            nlfaces = length(lface_to_face)
            lface_to_color = zeros(Int32,nlfaces)
            for (lface,face) in enumerate(lface_to_face)
                nodes = face_to_nodes[face]
                color = maximum(node->node_to_color[node],nodes)
                lface_to_color[lface] = color
            end
            oface_to_lface = findall(color->color==part,lface_to_color)
            hface_to_lface = findall(color->color!=part,lface_to_color)
            oface_to_face = lface_to_face[oface_to_lface]
            hface_to_face = lface_to_face[hface_to_lface]
            hface_to_color = lface_to_color[hface_to_lface]
            own = OwnIndices(nfaces,part,oface_to_face)
            ghost = GhostIndices(nfaces,hface_to_face,hface_to_color)
            OwnAndGhostIndices(own,ghost)
        end
        lface_to_face_mesh = map(local_to_global,local_faces)
        lnode_to_node_mesh = local_to_global(local_nodes)
        lmesh = restrict_mesh(mesh,lnode_to_node_mesh,lface_to_face_mesh)
        lmesh, local_nodes, Tuple(local_faces)
    end
    mesh_partition, node_partition, face_partition_array = map(setup,parts) |> tuple_of_arrays
    face_partition = face_partition_array |> tuple_of_arrays
    if renumber
        node_partition = renumber_partition(node_partition)
        face_partition = map(renumber_partition,face_partition)
    end
    # TODO here we have the opportunity to provide the parts rcv
    assembly_graph(node_partition)
    pmesh = PMesh(mesh_partition,node_partition,face_partition,partition_strategy)
    pmesh
end

function partition_mesh_cells(cell_to_color,parts,mesh,graph,partition_strategy,renumber)
    ghost_layers = partition_strategy.ghost_layers
    D = num_dims(mesh)
    ncells = num_faces(mesh,D)
    topo = topology(mesh)
    function setup(part)
        ocell_to_cell = findall(color->color==part,cell_to_color)
        cell_to_mask = fill(false,ncells)
        if ghost_layers == 0
            cell_to_mask[ocell_to_cell] .= true
        elseif ghost_layers == 1
            for cell in ocell_to_cell
                pini = graph.colptr[cell]
                pend = graph.colptr[cell+1]-1
                for p in pini:pend
                    cell2 = graph.rowval[p]
                    cell_to_mask[cell2] = true
                end
            end
        else
            error("case not implemented")
        end
        lcell_to_cell = findall(cell_to_mask)
        lcell_to_color = cell_to_color[lcell_to_cell]
        ocell_to_lcell = findall(color->color==part,lcell_to_color)
        hcell_to_lcell = findall(color->color!=part,lcell_to_color)
        ocell_to_cell = lcell_to_cell[ocell_to_lcell]
        hcell_to_cell = lcell_to_cell[hcell_to_lcell]
        hcell_to_color = lcell_to_color[hcell_to_lcell]
        own = OwnIndices(ncells,part,ocell_to_cell)
        ghost = GhostIndices(ncells,hcell_to_cell,hcell_to_color)
        local_cells = OwnAndGhostIndices(own,ghost)
        local_faces = map(0:(D-1)) do d
            cell_to_faces = face_incidence(topo,D,d)
            nfaces = num_faces(mesh,d)
            face_to_mask = fill(false,nfaces)
            for cell in lcell_to_cell
                faces = cell_to_faces[cell]
                face_to_mask[faces] .= true
            end
            lface_to_face = findall(face_to_mask)
            nlfaces = length(lface_to_face)
            lface_to_color = zeros(Int32,nlfaces)
            face_to_cells = face_incidence(topo,d,D)
            for (lface,face) in enumerate(lface_to_face)
                cells = face_to_cells[face]
                color = maximum(cell->cell_to_color[cell],cells)
                lface_to_color[lface] = color
            end
            oface_to_lface = findall(color->color==part,lface_to_color)
            hface_to_lface = findall(color->color!=part,lface_to_color)
            oface_to_face = lface_to_face[oface_to_lface]
            hface_to_face = lface_to_face[hface_to_lface]
            hface_to_color = lface_to_color[hface_to_lface]
            own = OwnIndices(nfaces,part,oface_to_face)
            ghost = GhostIndices(nfaces,hface_to_face,hface_to_color)
            OwnAndGhostIndices(own,ghost)
        end
        push!(local_faces,local_cells)
        nnodes = num_nodes(mesh)
        cell_to_nodes = face_nodes(mesh,D)
        node_to_mask = fill(false,nnodes)
        for cell in lcell_to_cell
            nodes = cell_to_nodes[cell]
            node_to_mask[nodes] .= true
        end
        lnode_to_node = findall(node_to_mask)
        nlnodes = length(lnode_to_node)
        lnode_to_color = zeros(Int32,nlnodes)
        node_to_cells = generate_face_coboundary(cell_to_nodes,nnodes)
        for (lnode,node) in enumerate(lnode_to_node)
            cells = node_to_cells[node]
            color = maximum(cell->cell_to_color[cell],cells)
            lnode_to_color[lnode] = color
        end
        onode_to_lnode = findall(color->color==part,lnode_to_color)
        hnode_to_lnode = findall(color->color!=part,lnode_to_color)
        onode_to_node = lnode_to_node[onode_to_lnode]
        hnode_to_node = lnode_to_node[hnode_to_lnode]
        hnode_to_color = lnode_to_color[hnode_to_lnode]
        own = OwnIndices(nnodes,part,onode_to_node)
        ghost = GhostIndices(nnodes,hnode_to_node,hnode_to_color)
        local_nodes = OwnAndGhostIndices(own,ghost)
        lnode_to_node_mesh = local_to_global(local_nodes)
        lface_to_face_mesh = map(local_to_global,local_faces)
        lmesh = restrict_mesh(mesh,lnode_to_node_mesh,lface_to_face_mesh)
        lmesh, local_nodes, Tuple(local_faces)
    end
    mesh_partition, node_partition, face_partition_array = map(setup,parts) |> tuple_of_arrays
    face_partition = face_partition_array |> tuple_of_arrays
    if renumber
        node_partition = renumber_partition(node_partition)
        face_partition = map(renumber_partition,face_partition)
    end
    # TODO here we have the opportunity to provide the parts rcv
    assembly_graph(node_partition)
    map(assembly_graph,face_partition)
    pmesh = PMesh(mesh_partition,node_partition,face_partition,partition_strategy)
    pmesh
end

struct PartitionStrategy{A,B} <: GT.AbstractType
    graph_nodes::Symbol
    graph_edges::Symbol
    graph_nodes_dim::A
    graph_edges_dim::B
    ghost_layers::Int
end

function partition_strategy(;
    graph_nodes=:cells,
    graph_edges=:nodes,
    ghost_layers=1,
    graph_nodes_dim=nothing,
    graph_edges_dim=nothing)

    @assert graph_nodes in (:cells,:nodes,:faces)
    @assert graph_edges in (:cells,:nodes,:faces)

    if graph_nodes ∉ (:cells,:nodes) && graph_nodes_dim === nothing
        error("graph_nodes_dim needs to be defined")
    end

    if graph_edges ∉ (:cells,:nodes) && graph_edges_dim === nothing
        error("graph_edges_dim needs to be defined")
    end

    PartitionStrategy(
                      graph_nodes,
                      graph_edges,
                      graph_nodes_dim,
                      graph_edges_dim,
                      ghost_layers)

end

function mesh_graph(mesh::AbstractMesh;
    partition_strategy=GT.partition_strategy())
    graph_nodes = partition_strategy.graph_nodes
    graph_edges = partition_strategy.graph_edges
    graph_nodes_dim = partition_strategy.graph_nodes_dim
    graph_edges_dim = partition_strategy.graph_edges_dim
    function barrier(nnodes,d_to_cell_to_nodes)
        ndata = 0
        for cell_to_nodes in d_to_cell_to_nodes
            ncells = length(cell_to_nodes)
            for cell in 1:ncells
                nodes = cell_to_nodes[cell]
                nlnodes = length(nodes)
                ndata += nlnodes*nlnodes
            end
        end
        I = zeros(Int32,ndata)
        J = zeros(Int32,ndata)
        p = 0
        for cell_to_nodes in d_to_cell_to_nodes
            ncells = length(cell_to_nodes)
            for cell in 1:ncells
                nodes = cell_to_nodes[cell]
                nlnodes = length(nodes)
                for j in 1:nlnodes
                    for i in 1:nlnodes
                        p += 1
                        I[p] = nodes[i]
                        J[p] = nodes[j]
                    end
                end
            end
        end
        V = ones(Int8,ndata)
        g = sparse(I,J,V,nnodes,nnodes)
        fill!(g.nzval,Int8(1))
        g
    end

    if graph_nodes === :nodes
        nnodes = num_nodes(mesh)
        if graph_edges === :cells
            face_nodes_mesh = face_nodes(mesh,num_dims(mesh))
            return barrier(nnodes,[face_nodes_mesh])
        elseif graph_edges === :faces && graph_edges_dim === :all
            face_nodes_mesh = face_nodes(mesh)
            return barrier(nnodes,face_nodes_mesh)
        elseif graph_edges === :faces
            face_nodes_mesh = face_nodes(mesh,graph_edges_dim)
            return barrier(nnodes,[face_nodes_mesh])
        else
            error("case not implemented")
        end
    elseif graph_nodes === :cells
        D = num_dims(mesh)
        ndfaces = num_faces(mesh,D)
        if graph_edges === :nodes
            dface_to_nodes = face_nodes(mesh,D)
            nnodes = num_nodes(mesh)
            node_to_dfaces = generate_face_coboundary(dface_to_nodes,nnodes)
            return barrier(ndfaces,[node_to_dfaces])
        elseif graph_edges === :faces && graph_edges_dim !== :all
            topo = topology(mesh)
            node_to_dfaces = face_incidence(topo,graph_edges_dim,D)
            return barrier(ndfaces,[node_to_dfaces])
        else
            error("case not implemented")
        end
    else
        error("case not implemented")
    end
end





#function two_level_mesh(coarse_mesh,fine_mesh;boundary_names=nothing)
#
#    D = num_dims(fine_mesh)
#    n_fine_cells = num_faces(fine_mesh,D)
#    n_fine_nodes = num_nodes(fine_mesh)
#    fine_refcell = first(reference_spaces(fine_mesh,D))
#    coarse_refcell = first(reference_spaces(coarse_mesh,D))
#    n_coarse_cells = num_faces(coarse_mesh,D)
#    fine_cell_local_node_to_fine_node = face_nodes(fine_mesh,D)
#    coarse_cell_lnode_to_coarse_node = face_nodes(coarse_mesh,D)
#    coarse_node_to_x = node_coordinates(coarse_mesh)
#    topo = topology(coarse_mesh)
#
#    # The user can provide custom names for physical groups on the boundary
#    # but we give some default value
#    if boundary_names === nothing
#        boundary_names = [
#            [ "$d-face-$face" for face in 1:num_faces(boundary(coarse_refcell),d)] 
#            for d in 0:(D-1)]
#    end
#    name_priority = reduce(vcat,boundary_names)
#    fine_node_groups = physical_nodes(fine_mesh;merge_dims=true,disjoint=true,name_priority)
#
#    # Recover boundary info
#    d_to_local_dface_to_fine_nodes = Vector{Vector{Vector{Int}}}(undef,D+1)
#    fine_node_mask = fill(true,n_fine_nodes)
#    for d in 0:(D-1)
#        n_local_dfaces = num_faces(boundary(coarse_refcell),d)
#        local_dface_to_fine_nodes = Vector{Vector{Int}}(undef,n_local_dfaces)
#        for local_dface in 1:n_local_dfaces
#            fine_nodes = fine_node_groups[boundary_names[d+1][local_dface]]
#            local_dface_to_fine_nodes[local_dface] = fine_nodes
#            fine_node_mask[fine_nodes] .= false
#        end
#        d_to_local_dface_to_fine_nodes[d+1] = local_dface_to_fine_nodes
#    end
#    d_to_local_dface_to_fine_nodes[D+1] = [findall(fine_node_mask)]
#    fine_node_to_d = zeros(Int,n_fine_nodes)
#    fine_node_to_local_dface = zeros(Int,n_fine_nodes)
#    for d in 0:D
#        n_local_dfaces = length(d_to_local_dface_to_fine_nodes[d+1])
#        for local_dface in 1:n_local_dfaces
#            fine_nodes = d_to_local_dface_to_fine_nodes[d+1][local_dface]
#            fine_node_to_d[fine_nodes] .= d
#            fine_node_to_local_dface[fine_nodes] .= local_dface
#        end
#    end
#
#    ## Ensure consistent mapping of periodic nodes on opposite faces 
#    # Get periodicity information
#    fine_node_to_master_node = collect(1:n_fine_nodes)
#    d_to_local_dface_to_permutation = Vector{Vector{Vector{Int}}}(undef, D+1)
#    periodic_node_to_fine_node, periodic_node_to_master_node = periodic_nodes(
#        fine_mesh)
#    fine_node_to_master_node[periodic_node_to_fine_node] = periodic_node_to_master_node
#
#    # Fixing vertices indirection for unit cell
#    # Assumes one level of indirection for periodicity (e.g., vertex -> master -> master)
#    # TODO: d-1 levels of indirection
#    for _ in 1:(D-1)
#        for fine_node in 1:n_fine_nodes
#            master_node = fine_node_to_master_node[fine_node]
#            if fine_node == master_node
#                continue  
#            end
#            master_master_node = fine_node_to_master_node[master_node]
#            fine_node_to_master_node[fine_node] = master_master_node 
#        end 
#    end
#
#    # Sort edges according master node (only needed in 3D)
#    if D==3
#        let d=1
#            n_local_dfaces = length(d_to_local_dface_to_fine_nodes[d+1])
#            for local_dface in 1:n_local_dfaces
#                fine_nodes = d_to_local_dface_to_fine_nodes[d+1][local_dface]
#                master_nodes = fine_node_to_master_node[fine_nodes]
#                d_to_local_dface_to_fine_nodes[d+1][local_dface] = fine_nodes[sortperm(master_nodes)]
#            end
#        end
#    end
#
#    fine_node_to_permuted_node = zeros(Int, n_fine_nodes)
#    fine_node_to_face_node = zeros(Int, n_fine_nodes)
#    d_to_local_dface_to_opposite_dface = opposite_faces(geometry(coarse_refcell))
#    d_to_local_dface_to_is_master = Vector{Vector{Bool}}(undef, D+1)
#    for d in 0:(D-1)
#        local_dface_to_is_master = similar(d_to_local_dface_to_opposite_dface[d+1], Bool) 
#        n_local_dfaces = num_faces(boundary(coarse_refcell),d)
#        for local_dface in 1:n_local_dfaces 
#            local_dface_to_is_master[local_dface] = true
#        end
#        for local_dface in 1:n_local_dfaces 
#            if  local_dface_to_is_master[local_dface] == false
#                continue
#            end
#            opposite_local_dface = d_to_local_dface_to_opposite_dface[d+1][local_dface]
#            local_dface_to_is_master[opposite_local_dface] = false
#        end
#        d_to_local_dface_to_is_master[d+1] = local_dface_to_is_master
#    end 
#    
#    for d in 0:(D-1)
#        n_local_dfaces = num_faces(boundary(coarse_refcell),d)
#        local_dface_to_permutation = Vector{Vector{Int}}(undef, n_local_dfaces)
#        local_dface_to_fine_nodes = d_to_local_dface_to_fine_nodes[d+1]
#        for local_dface_1 in 1:n_local_dfaces
#
#            fine_nodes_1 = local_dface_to_fine_nodes[local_dface_1]
#            master_nodes_1 = fine_node_to_master_node[fine_nodes_1]
#
#            # NB! The permutation thing only works for D-1 objects!
#            # Handle nonperiodic case with identity permutation 
#            # assumes consistent numbering of node ids on opposite faces... 
#            # Also handles case where there are periodic nodes but the master nodes 
#            # will use the identity permutation
#            # TODO: in 3d geometry, since surfaces represent the periodic copies,
#            # the opposite edges can end up being a map from slaves to slave, and therefore
#            # this condition below fails because there is no obvious master (as in the 2d case) 
#            if d !=D-1 || length(periodic_node_to_fine_node) == 0 || d_to_local_dface_to_is_master[d+1][local_dface_1] 
#
#                permutation = collect(1:length(local_dface_to_fine_nodes[local_dface_1]))
#                local_dface_to_permutation[local_dface_1] = permutation
#                fine_node_to_permuted_node[fine_nodes_1] = permutation
#                fine_node_to_face_node[fine_nodes_1] = 1:length(fine_nodes_1)
#                continue 
#            end
#
#            # Handle periodic case: use a local reference cell and its opposite 
#            # face to get the corresponding fine nodes and the master of those fine nodes.
#            # Then ensure that a given local d-face has the same ordering as the opposite
#            # face....
#            local_dface_2 = d_to_local_dface_to_opposite_dface[d+1][local_dface_1]
#            fine_nodes_2 = local_dface_to_fine_nodes[local_dface_2]
#            master_nodes_2 = fine_node_to_master_node[fine_nodes_2]
#            permutation = indexin(master_nodes_2, master_nodes_1)
#            local_dface_to_permutation[local_dface_1] = permutation
#            fine_node_to_permuted_node[fine_nodes_1] = permutation
#            fine_node_to_face_node[fine_nodes_1] = 1:length(fine_nodes_1)
#        end
#        d_to_local_dface_to_permutation[d+1] = local_dface_to_permutation
#    end
#    d_to_local_dface_to_permutation[D+1] = [
#        collect(1:length(d_to_local_dface_to_fine_nodes[D+1][1]))]
#
#    # Setup the map of fine nodes to physical coordinates via finite element interpolation
#    fine_node_to_x = node_coordinates(fine_mesh)
#    A = tabulator(coarse_refcell)(value,fine_node_to_x)
#    coarse_cell_fine_node_to_x = Vector{Vector{SVector{D,Float64}}}(undef,n_coarse_cells)
#    for coarse_cell in 1:n_coarse_cells
#        lnode_to_coarse_node = coarse_cell_lnode_to_coarse_node[coarse_cell]
#        lnode_to_x = coarse_node_to_x[lnode_to_coarse_node]
#        coarse_cell_fine_node_to_x[coarse_cell] = A*lnode_to_x
#    end
#
#    # Glue fine node ids
#    # TODO ordering in the physical group
#    # start with a 2x2 unit cell
#    final_node = 0
#    d_coarse_dface_to_offset = Vector{Vector{Int}}(undef,D+1)
#    for d in 0:D
#        local_dface_to_fine_nodes = d_to_local_dface_to_fine_nodes[d+1]
#        coarse_dface_to_coarse_cells = face_incidence(topo,d,D)
#        coarse_cell_to_coarse_dfaces = face_incidence(topo,D,d)
#        n_coarse_dfaces = num_faces(coarse_mesh,d)
#        coarse_dface_to_offset = zeros(Int,n_coarse_dfaces)
#        for coarse_dface in 1:n_coarse_dfaces
#            coarse_cells = coarse_dface_to_coarse_cells[coarse_dface]
#            for coarse_cell in coarse_cells
#                coarse_dfaces = coarse_cell_to_coarse_dfaces[coarse_cell]
#                local_dface = findfirst(i->coarse_dface==i,coarse_dfaces)
#                fine_nodes = local_dface_to_fine_nodes[local_dface]
#                if coarse_dface_to_offset[coarse_dface] == 0
#                    coarse_dface_to_offset[coarse_dface] = final_node
#                    final_node += length(fine_nodes)
#                end
#            end    
#        end
#        d_coarse_dface_to_offset[d+1] = coarse_dface_to_offset
#    end
#    n_final_nodes = final_node # renaming of counter variable for clarity 
#
#    # Initialize map of coarse cell to final node mediated by fine node
#    coarse_cell_fine_node_to_final_node = Vector{Vector{Int}}(undef,n_coarse_cells)
#    for coarse_cell in 1:n_coarse_cells
#        coarse_cell_fine_node_to_final_node[coarse_cell] = zeros(Int,n_fine_nodes)
#    end
#
#    # Map each dimension to final nodes mediated by coarse d-dimensional face
#    n_coarse_nodes = num_nodes(coarse_mesh)
#    d_to_coarse_dface_to_final_nodes = Vector{Vector{Vector{Int}}}(undef,D+1)
#    for d in 0:D
#        local_dface_to_fine_nodes = d_to_local_dface_to_fine_nodes[d+1]
#        coarse_dface_to_coarse_cells = face_incidence(topo,d,D)
#        coarse_cell_to_coarse_dfaces = face_incidence(topo,D,d)
#        n_coarse_dfaces = num_faces(coarse_mesh,d)
#        coarse_dface_to_final_nodes = Vector{Vector{Int}}(undef,n_coarse_dfaces)
#        for coarse_dface in 1:n_coarse_dfaces
#            offset = d_coarse_dface_to_offset[d+1][coarse_dface]
#            coarse_cells = coarse_dface_to_coarse_cells[coarse_dface]
#            for coarse_cell in coarse_cells
#                coarse_dfaces = coarse_cell_to_coarse_dfaces[coarse_cell]
#                local_dface = findfirst(i->coarse_dface==i,coarse_dfaces)
#                permutation = d_to_local_dface_to_permutation[d+1][local_dface]
#                fine_nodes = local_dface_to_fine_nodes[local_dface][permutation]
#                final_nodes =  offset .+ (1:length(fine_nodes)) 
#                coarse_cell_fine_node_to_final_node[coarse_cell][fine_nodes] = final_nodes
#                coarse_dface_to_final_nodes[coarse_dface] = final_nodes
#            end    
#        end
#        d_to_coarse_dface_to_final_nodes[d+1] = coarse_dface_to_final_nodes
#    end
#
#    # Apply the coordinate transformation to the final nodes 
#    final_node_to_x = zeros(SVector{D,Float64},n_final_nodes)
#    for coarse_cell in 1:n_coarse_cells
#        fine_node_to_final_node = coarse_cell_fine_node_to_final_node[coarse_cell]  
#        fine_node_to_x = coarse_cell_fine_node_to_x[coarse_cell] 
#        final_node_to_x[fine_node_to_final_node] = fine_node_to_x
#    end
#
#    # Map final cells to final node mediated by local (finite element) reference entities
#    n_final_cells = n_coarse_cells*n_fine_cells
#    final_cell_local_node_to_final_node = Vector{Vector{Int}}(undef,n_final_cells)
#
#    final_cell = 0
#    final_cell_to_coarse_cell = Vector{Int}(undef, n_final_cells)
#    for coarse_cell in 1:n_coarse_cells
#        for fine_cell in 1:n_fine_cells
#            local_node_to_fine_node = fine_cell_local_node_to_fine_node[fine_cell]
#            local_node_to_final_node = coarse_cell_fine_node_to_final_node[
#                coarse_cell][local_node_to_fine_node]
#            final_cell += 1
#            final_cell_local_node_to_final_node[final_cell] = local_node_to_final_node
#            final_cell_to_coarse_cell[final_cell] = coarse_cell
#        end
#    end
#    final_cell_to_refid = fill(1,n_final_cells)
#    refid_to_fine_refcell = [fine_refcell]
#
#    chain = chain_from_arrays(
#                final_node_to_x,
#                JaggedArray(final_cell_local_node_to_final_node),
#                final_cell_to_refid,
#                refid_to_fine_refcell)
#
#    # TODO we could avoid this call to complexify by using
#    # the fine faces (just as we did with the fine nodes)
#    # TODO maybe we don't need to complexify and only find the faces
#    # needed for the physical groups
#    final_mesh = mesh_from_chain(chain) |> complexify
#
#    # Refine physical groups
#    final_node_to_mask = fill(false,n_final_nodes)
#    for d in 0:D
#        final_physical_dfaces = physical_faces(final_mesh,d)
#        coarse_phsycial_dfaces = physical_faces(coarse_mesh,d)
#        final_dfaces_to_final_nodes = face_nodes(final_mesh,d)
#        for (name,coarse_dfaces_in_group) in coarse_phsycial_dfaces
#            fill!(final_node_to_mask,false)
#            for n in 0:d
#                coarse_dface_to_coarse_nfaces = face_incidence(topo,d,n)
#                coarse_cell_to_coarse_nfaces = face_incidence(topo,D,n)
#                coarse_nface_to_coarse_cells = face_incidence(topo,n,D)
#                local_nface_to_fine_nodes = d_to_local_dface_to_fine_nodes[n+1]
#                for coarse_dface in coarse_dfaces_in_group
#                    for coarse_nface  in coarse_dface_to_coarse_nfaces[coarse_dface]
#                        coarse_cells = coarse_nface_to_coarse_cells[coarse_nface]
#                        coarse_cell = first(coarse_cells)
#                        coarse_nfaces = coarse_cell_to_coarse_nfaces[coarse_cell]
#                        local_nface = findfirst(i->coarse_nface==i,coarse_nfaces)
#                        fine_nodes = local_nface_to_fine_nodes[local_nface]
#                        final_nodes = coarse_cell_fine_node_to_final_node[coarse_cell][fine_nodes]
#                        final_node_to_mask[final_nodes] .=true
#                    end
#                end
#                final_faces_in_group = findall(final_nodes->all(final_node->final_node_to_mask[final_node],final_nodes),final_dfaces_to_final_nodes)
#                final_physical_dfaces[name] = final_faces_in_group
#            end
#        end
#    end
#
#    glue = (;
#        d_to_coarse_dface_to_final_nodes,
#        coarse_cell_fine_node_to_final_node,
#        d_to_local_dface_to_fine_nodes,
#        final_cell_to_coarse_cell)
#
#    final_mesh, glue
#end
#
#function two_level_mesh(coarse_mesh::PMesh,fine_mesh;kwargs...)
#    # TODO for the moment we assume a cell-based partition without ghosts
#    # TODO: NEED TO ASSERT CELL BASED PARTITION???
#    D = num_dims(fine_mesh)
#
#    function setup_local_meshes(my_coarse_mesh)
#        two_level_mesh(my_coarse_mesh,fine_mesh; kwargs...)
#    end
#    mesh_partition, glue = map(setup_local_meshes, partition(coarse_mesh)) |> tuple_of_arrays
#
#    # mark owernship of nodes using a local final mesh and local glue
#    function mark_nodes(final_mesh, local_glue, coarse_indices)
#        d_to_coarse_dface_to_final_nodes = local_glue.d_to_coarse_dface_to_final_nodes
#        my_final_node_to_owner = fill(0,num_nodes(final_mesh))
#        for d in 0:D
#            coarse_dface_to_final_nodes = d_to_coarse_dface_to_final_nodes[d+1]
#            coarse_dfaces = face_indices(coarse_indices,d)
#            coarse_dface_to_owner = local_to_owner(coarse_dfaces)
#            n_coarse_dfaces = length(coarse_dface_to_owner)
#            for coarse_dface in 1:n_coarse_dfaces
#                owner = coarse_dface_to_owner[coarse_dface]
#                final_nodes = coarse_dface_to_final_nodes[coarse_dface]
#                my_final_node_to_owner[final_nodes] .= owner
#            end
#        end
#        my_final_node_to_owner
#    end
#    index_partition_coarse_mesh = index_partition(coarse_mesh) # node/face ixs per part
#    final_node_to_owner = map(mark_nodes, mesh_partition, glue, index_partition_coarse_mesh) 
#    parts = linear_indices(final_node_to_owner)
#
#    # count the number of owned final mesh nodes per partition
#    n_own_final_nodes = map(
#        (owners,part)->count(owner->owner==part,owners),final_node_to_owner,parts)
#    n_final_nodes = sum(n_own_final_nodes) 
#
#    # owned indices of final mesh nodes for each partition
#    own_node_partition = variable_partition(n_own_final_nodes,n_final_nodes)
#
#    # assign a non-zero global id to final mesh nodes if they are owned by a partition
#    final_node_to_gid = map(v->zeros(Int,length(v)),final_node_to_owner)
#    map(final_node_to_gid,
#        final_node_to_owner,
#        own_node_partition,
#        parts) do gids, owners, own_nodes, part
#        own = 0
#        for i in 1:length(gids)
#            owner = owners[i]
#            if owner == part
#                own += 1
#                gids[i] = own_nodes[own]
#            end
#        end
#    end
#
#    # for each dimension of the coarse d-face, handle ownership and ensure consistent data 
#    for d in 0:D
#        # For a given d-dimensional (i.e., point, edge, surface, vol) coarse face,
#        # get a map from coarse face to final node, iterate through the coarse faces,
#        # and get a global id corresponding to that final node
#        function get_coarse_dface_to_gids(coarse_dfaces, local_glue, gids)
#            coarse_dface_to_final_nodes = local_glue.d_to_coarse_dface_to_final_nodes[d+1]
#            n_coarse_dfaces = local_length(coarse_dfaces)
#            coarse_dface_to_gids = Vector{Vector{Int}}(undef,n_coarse_dfaces)
#            for coarse_dface in 1:n_coarse_dfaces
#                final_nodes = coarse_dface_to_final_nodes[coarse_dface]
#                coarse_dface_to_gids[coarse_dface] = gids[final_nodes]
#            end
#            JaggedArray(coarse_dface_to_gids)
#        end
#        coarse_dface_partition = face_partition(coarse_mesh, d)
#        coarse_dface_to_gids_data = map( # this is a map over partition specific data
#            get_coarse_dface_to_gids, coarse_dface_partition, glue, final_node_to_gid)
#        coarse_dface_to_gids = PVector(coarse_dface_to_gids_data, coarse_dface_partition)
#        consistent!(coarse_dface_to_gids) |> wait # owners give their values to ghosts?
#
#        # For each d-dimensional coarse face, update the non-partitioned gids 
#        # with the partitioned and consistent gids
#        function update_gids!(coarse_dfaces, coarse_dface_to_gids, local_glue, gids)
#            coarse_dface_to_final_nodes = local_glue.d_to_coarse_dface_to_final_nodes[d+1]
#            n_coarse_dfaces = local_length(coarse_dfaces) 
#            for coarse_dface in 1:n_coarse_dfaces
#                final_nodes = coarse_dface_to_final_nodes[coarse_dface]
#                gids[final_nodes] = coarse_dface_to_gids[coarse_dface]
#            end
#        end
#        map(update_gids!, 
#            coarse_dface_partition, 
#            coarse_dface_to_gids_data, 
#            glue, 
#            final_node_to_gid)
#    end
#
#    # For each part, use the final node (i.e., local id) to global id map
#    # and final node to owner map to construct arbitrary indices for partitioned 
#    # mesh nodes 
#    function finalize_node_partition(part, gids, owners)
#        LocalIndices(n_final_nodes, part, gids, owners)
#    end
#    node_partition = map(
#        finalize_node_partition, parts, final_node_to_gid, final_node_to_owner)
#
#    # TODO: variable_partition has ghosts and periodic positional args... how to handle this?
#    n_own_cells = map(final_mesh -> num_faces(final_mesh, D), mesh_partition) 
#    n_cells = sum(n_own_cells)
#    cell_partition = variable_partition(n_own_cells, n_cells) 
#
#   
#    # TODO: 0 and 1 dimensional faces have ghost cells, so the use of
#    # variable_partition here is just a placeholder and is not accurate
#    function dummy_face_partition(d)
#        n_own_dfaces = map(mesh -> num_faces(mesh, d), mesh_partition)
#        n_dfaces = sum(n_own_dfaces)
#        return variable_partition(n_own_dfaces, n_dfaces)
#    end
#    _face_partition = ntuple(
#        i-> i == (D+1) ? cell_partition : dummy_face_partition(i-1),
#        D+1)
#
#    final_glue = nothing # TODO: placeholder for parallel glue
#    partition_strategy = coarse_mesh.partition_strategy # TODO: is this right? 
#    final_mesh = PMesh(
#        mesh_partition, node_partition, _face_partition, partition_strategy)
#    return final_mesh, final_glue
#end

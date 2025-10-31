
"""
    cartesian_mesh(domain,cells_per_dir)

Create a multi-dimensional Cartesian mesh. The dimension of the mesh
is defined by the length of `cells_per_dir`.
The number of cells in direction `i` is given by `cells_per_dir[i]`. The extends
of the Cartesian mesh are given in `domain`. The range in direction `i` covered
by the mesh is given by `domain[2*i-1,2*i]`.

# Keyword arguments
- `boundary=true` [optional]: Include faces on the boundary and generate face groups identifying which faces are  on which face of  bounding box of the mesh. The groups are named `\$d-face-\$i` for the face `i` of dimension `d` of the bounding box.
- `complexify=true` [optional]: Generate all low dimensional faces so that the mesh is  a cell complex.
- `simplexify=false` [optional]: Generate a mesh of simplex faces instead of hyper-cubes.
"""
function cartesian_mesh(
    domain,cells_per_dir;
    boundary=Val(true),
    complexify=Val(true),
    simplexify=Val(false),
    periodic=nothing,
    )
    mesh0 = if val_parameter(boundary)
        if val_parameter(simplexify)
            structured_simplex_mesh_with_boundary(domain,cells_per_dir)
        else
            cartesian_mesh_with_boundary(domain,cells_per_dir)
        end
    else
        if val_parameter(simplexify)
            chain = structured_simplex_chain(domain,cells_per_dir)
        else
            chain = cartesian_chain(domain,cells_per_dir)
        end
        GT.mesh(chain)
    end
    if periodic !== nothing
        if val_parameter(periodic) == true
            periodic_nodes = collect(Int32,1:num_nodes(mesh0))
        elseif val_parameter(periodic) == false
            periodic_nodes = 1:num_nodes(mesh0)
        else
            periodic_nodes = cartesian_periodic_nodes(domain,cells_per_dir,periodic)
        end
        mesh = replace_periodic_nodes(mesh0,periodic_nodes)
    else
        mesh = mesh0
    end
    if val_parameter(complexify)
        GT.complexify(mesh)
    else
        mesh
    end
end

#function cartesian_pmesh(domain,cells_per_dir,parts,parts_per_dir;
#    boundary=true,
#    complexify=true,
#    simplexify=false,
#    partition_strategy = GT.partition_strategy()
#    )
#
#    mesh = cartesian_mesh(domain,cells_per_dir;boundary,complexify,simplexify)
#
#    graph_nodes = partition_strategy.graph_nodes
#    graph_edges = partition_strategy.graph_edges
#    ghost_layers = partition_strategy.ghost_layers
#    @assert graph_nodes === :cells
#    @assert graph_edges === :nodes
#    np = prod(parts_per_dir)
#    graph = mesh_graph(mesh;partition_strategy)
#    parts_seq = LinearIndices((np,))
#    graph_partition = zeros(Int,prod(cells_per_dir))
#    for ids in uniform_partition(LinearIndices((np,)),parts_per_dir,cells_per_dir)
#        graph_partition[local_to_global(ids)] = local_to_owner(ids)
#    end
#    partition_mesh(mesh,np;partition_strategy,parts,graph_partition)
#end

function bounding_box_from_domain(domain)
    l = length(domain)
    D = div(l,2)
    pmin = SVector(ntuple(d->domain[2*(d-1)+1],Val(D)))
    pmax = SVector(ntuple(d->domain[2*d],Val(D)))
    (pmin,pmax)
end

function domain_from_bounding_box(box)
    l = sum(length,box)
    ntuple(Val(l)) do i
        vector = mod(i-1,2)+1
        component = div(i-1,2)+1
        box[vector][component]
    end
end

function cartesian_mesh_with_boundary(domain,cells_per_dir)
    if any(i->i!=1,cells_per_dir) && any(i->i<2,cells_per_dir)
        error("At least 2 cells in any direction (or 1 cell in all directions)")
    end
    function barrier(
      D,
      cell_to_nodes,
      nnodes,
      d_to_ldface_to_lnodes)

      node_to_n = zeros(Int32,nnodes)
      for nodes in cell_to_nodes
        for node in nodes
          node_to_n[node] += Int32(1)
        end
      end
      J = typeof(JaggedArray(Vector{Int32}[]))
      face_to_nodes = Vector{J}(undef,D)
      groups = [ Dict{String,Vector{Int32}}() for d in 0:(D-1) ]
      ngroups = 0
      for d in 0:(D-1)
        nmax = 2^d
        ldface_to_lnodes = d_to_ldface_to_lnodes[d+1]
        ndfaces = 0
        for nodes in cell_to_nodes
          for (ldface,lnodes) in enumerate(ldface_to_lnodes)
            isboundary = true
            for lnode in lnodes
              node = nodes[lnode]
              if node_to_n[node] > nmax
                isboundary = false
                break
              end
            end
            if isboundary
              ndfaces += 1
            end
          end
        end
        ptrs = zeros(Int32,ndfaces+1)
        for dface in 1:ndfaces
          ptrs[dface+1] += Int32(nmax)
        end
        length_to_ptrs!(ptrs)
        ndata = ptrs[end]-1
        data = zeros(Int32,ndata)
        dface_to_physical_group = zeros(Int32,ndfaces)
        ndfaces = 0
        for nodes in cell_to_nodes
          for (ldface,lnodes) in enumerate(ldface_to_lnodes)
            isboundary = true
            for lnode in lnodes
              node = nodes[lnode]
              if node_to_n[node] > nmax
                isboundary = false
                break
              end
            end
            if isboundary
              ndfaces += 1
              group = ngroups + ldface
              dface_to_physical_group[ndfaces] = group
              p = ptrs[ndfaces]-Int32(1)
              for (i,lnode) in enumerate(lnodes)
                node = nodes[lnode]
                data[p+i] = node
              end
            end
          end
        end
        nldfaces = length(ldface_to_lnodes)
        face_to_nodes[d+1] = JaggedArray(data,ptrs)
        for ldface in 1:nldfaces
          group = ngroups + ldface
          group_name = "$(d)-face-$ldface"
          faces_in_physical_group = findall(g->g==group,dface_to_physical_group)
          groups[d+1][group_name] = faces_in_physical_group
        end
        ngroups += nldfaces
        if d == (D-1)
            groups[d+1]["boundary"] = 1:length(dface_to_physical_group)
        end
      end # d
      ngroups += 1
      groups, face_to_nodes
    end # barrier
    chain = cartesian_chain(domain,cells_per_dir)
    interior_mesh = GT.mesh(chain)
    D = num_dims(interior_mesh)
    cell_to_nodes = face_nodes(interior_mesh,D)
    reference_cells = reference_spaces(interior_mesh,D)
    node_coords = node_coordinates(interior_mesh)
    ref_cell = first(reference_cells)
    refid_to_refface = reference_spaces(remove_interior(complexify(ref_cell)))
    nnodes = num_nodes(interior_mesh)
    d_to_ldface_to_lnodes = [face_nodes(remove_interior(complexify(ref_cell)),d) for d in 0:(D-1)]
    groups, face_to_nodes = barrier(
      D,
      cell_to_nodes,
      nnodes,
      d_to_ldface_to_lnodes)
    face_to_refid = [ ones(Int8,length(face_to_nodes[d+1]))  for d in 0:(D-1)]
    mesh_face_nodes = push(face_to_nodes,face_nodes(interior_mesh,D))
    mesh_face_reference_id = push(face_to_refid,face_reference_id(interior_mesh,D))
    mesh_reference_spaces = push(refid_to_refface,reference_spaces(interior_mesh,D))
    mesh_groups = push(groups,group_faces(interior_mesh,D))
    GT.mesh(
     node_coordinates = node_coords,
     face_nodes = mesh_face_nodes,
     face_reference_id = mesh_face_reference_id,
     reference_spaces = mesh_reference_spaces,
     geometry_names = [ [ "$d-face-$face" for face in 1:num_faces(GT.mesh(GT.domain(ref_cell)),d) ] for d in 0:D],
     group_faces=mesh_groups,
    )
end

function cartesian_chain(domain,cells_per_dir)
    box = bounding_box_from_domain(domain)
    D = length(cells_per_dir)
    nodes_per_dir = cells_per_dir .+ 1
    pmin = first(box)
    pmax = last(box)
    extent_per_dir = pmax .- pmin
    h_per_dir = SVector(extent_per_dir ./ cells_per_dir)
    nnodes = prod(nodes_per_dir)
    ncells = prod(cells_per_dir)
    nlnodes = 2^D
    cell_nodes_ptrs = fill(Int32(nlnodes),ncells+1)
    cell_nodes_ptrs[1] = 0
    length_to_ptrs!(cell_nodes_ptrs)
    cell_nodes_data = zeros(Int32,ncells*nlnodes)
    cell_nodes = JaggedArray(cell_nodes_data,cell_nodes_ptrs)
    cell_cis = CartesianIndices(cells_per_dir)
    cell_lis = LinearIndices(cells_per_dir)
    node_cis = CartesianIndices(nodes_per_dir)
    node_lis = LinearIndices(nodes_per_dir)
    lnode_cis = CartesianIndices(ntuple(i->0:1,Val(D)))
    for (cell_li,cell_ci) in enumerate(cell_cis)
        nodes = cell_nodes[cell_li]
        for (lnode_li,lnode_ci) in enumerate(lnode_cis)
            node_ci = CartesianIndex(Tuple(cell_ci) .+ Tuple(lnode_ci))
            node_li = node_lis[node_ci]
            nodes[lnode_li] = node_li
        end
    end
    node_coords = zeros(SVector{D,Float64},nnodes)
    for (node_li,node_ci) in enumerate(node_cis)
        anchor = SVector(Tuple(node_ci) .- 1)
        x = pmin .+ h_per_dir .* anchor
        node_coords[node_li] = x
    end
    cell_reference_id = fill(Int32(1),ncells)
    cell_geometry = unit_n_cube(Val(D))
    order = 1
    ref_cell = lagrange_space(cell_geometry,order)
    reference_cells = [ref_cell]
    interior_cells = collect(Int32,1:length(cell_nodes))
    groups = Dict(["interior"=>interior_cells,"$D-face-1"=>interior_cells])
    chain = GT.chain(
        node_coordinates = node_coords,
        face_nodes = cell_nodes,
        face_reference_id = cell_reference_id,
        reference_spaces = reference_cells,
        group_faces=groups,
       )
    chain
end

function structured_simplex_chain(domain,cells_per_dir)
    box = bounding_box_from_domain(domain)
    D = length(cells_per_dir)
    nodes_per_dir = cells_per_dir .+ 1
    pmin = first(box)
    pmax = last(box)
    extent_per_dir = pmax .- pmin
    h_per_dir = SVector(extent_per_dir ./ cells_per_dir)
    nnodes = prod(nodes_per_dir)
    nlnodes = 2^D
    cell_geometry = unit_n_cube(Val(D))
    ref_simplex_mesh = simplexify(cell_geometry)
    rscell_to_lnodes = face_nodes(ref_simplex_mesh,D)
    nrscells = length(rscell_to_lnodes)
    nslnodes = length(first(rscell_to_lnodes))
    ncells = prod(cells_per_dir)*nrscells
    cell_nodes_ptrs = fill(Int32(nslnodes),ncells+1)
    cell_nodes_ptrs[1] = 0
    length_to_ptrs!(cell_nodes_ptrs)
    ndata = cell_nodes_ptrs[end]-1
    cell_nodes_data = zeros(Int32,ndata)
    cell_nodes = JaggedArray(cell_nodes_data,cell_nodes_ptrs)
    cell_cis = CartesianIndices(cells_per_dir)
    cell_lis = LinearIndices(cells_per_dir)
    node_cis = CartesianIndices(nodes_per_dir)
    node_lis = LinearIndices(nodes_per_dir)
    lnode_cis = CartesianIndices(ntuple(i->0:1,Val(D)))
    clnodes = zeros(Int,nlnodes)
    scell = 0
    for (cell_li,cell_ci) in enumerate(cell_cis)
        for (lnode_li,lnode_ci) in enumerate(lnode_cis)
            node_ci = CartesianIndex(Tuple(cell_ci) .+ Tuple(lnode_ci))
            node_li = node_lis[node_ci]
            clnodes[lnode_li] = node_li
        end
        for srcell in 1:nrscells
            scell += 1
            nodes = cell_nodes[scell]
            lnodes = rscell_to_lnodes[srcell]
            for (i,lnode) in enumerate(lnodes)
                nodes[i] = clnodes[lnode]
            end
        end
    end
    node_coords = zeros(SVector{D,Float64},nnodes)
    for (node_li,node_ci) in enumerate(node_cis)
        anchor = SVector(Tuple(node_ci) .- 1)
        x = pmin .+ h_per_dir .* anchor
        node_coords[node_li] = x
    end
    cell_reference_id = fill(Int32(1),ncells)
    order = 1
    reference_cells = reference_spaces(ref_simplex_mesh,D)
    interior_cells = collect(Int32,1:length(cell_nodes))
    groups = Dict(["interior"=>interior_cells,"$D-face-1"=>interior_cells])
    chain = GT.chain(
        node_coordinates = node_coords,
        face_nodes = cell_nodes,
        face_reference_id = cell_reference_id,
        reference_spaces = reference_cells;
        group_faces=groups,
       )
    chain
end

function structured_simplex_mesh_with_boundary(domain,cells_per_dir)
    if any(i->i!=1,cells_per_dir) && any(i->i<2,cells_per_dir)
        error("At least 2 cells in any direction (or 1 cell in all directions)")
    end
    function barrier(
      D,
      cell_to_nodes,
      nnodes,
      d_to_ldface_to_lnodes,
      d_to_ldface_to_sldface_to_lnodes)

        node_to_n = zeros(Int32,nnodes)
        for nodes in cell_to_nodes
            for node in nodes
                node_to_n[node] += Int32(1)
            end
        end
        J = typeof(JaggedArray(Vector{Int32}[]))
        face_to_nodes = Vector{J}(undef,D)
        groups = [ Dict{String,Vector{Int32}}() for d in 0:(D-1) ]
        ngroups = 0
        for d in 0:(D-1)
            nmax = 2^d
            ldface_to_lnodes = d_to_ldface_to_lnodes[d+1]
            ldface_to_sldface_to_lnodes = d_to_ldface_to_sldface_to_lnodes[d+1]
            ndfaces = 0
            for nodes in cell_to_nodes
                for (ldface,lnodes) in enumerate(ldface_to_lnodes)
                    isboundary = true
                    for lnode in lnodes
                        node = nodes[lnode]
                        if node_to_n[node] > nmax
                            isboundary = false
                            break
                        end
                    end
                    if isboundary
                        ndfaces += length(ldface_to_sldface_to_lnodes[ldface])
                    end
                end
            end
            nslnodes = length(ldface_to_sldface_to_lnodes[begin][begin])
            ptrs = zeros(Int32,ndfaces+1)
            for dface in 1:ndfaces
                ptrs[dface+1] += Int32(nslnodes)
            end
            length_to_ptrs!(ptrs)
            ndata = ptrs[end]-1
            data = zeros(Int32,ndata)
            dface_to_physical_group = zeros(Int32,ndfaces)
            ndfaces = 0
            for nodes in cell_to_nodes
                for (ldface,lnodes) in enumerate(ldface_to_lnodes)
                    isboundary = true
                    for lnode in lnodes
                        node = nodes[lnode]
                        if node_to_n[node] > nmax
                            isboundary = false
                            break
                        end
                    end
                    if isboundary
                        group = ngroups + ldface
                        sldface_to_lnodes = ldface_to_sldface_to_lnodes[ldface]
                        nsldfaces = length(sldface_to_lnodes)
                        for sldface in 1:nsldfaces
                            ndfaces += 1
                            dface_to_physical_group[ndfaces] = group
                            p = ptrs[ndfaces]-Int32(1)
                            mylnodes = sldface_to_lnodes[sldface]
                            for (i,lnode) in enumerate(mylnodes)
                                node = nodes[lnode]
                                data[p+i] = node
                            end
                        end
                    end
                end
            end
            nldfaces = length(ldface_to_lnodes)
            face_to_nodes[d+1] = JaggedArray(data,ptrs)
            for ldface in 1:nldfaces
                group = ngroups + ldface
                group_name = "$(d)-face-$ldface"
                faces_in_physical_group = findall(g->g==group,dface_to_physical_group)
                groups[d+1][group_name] = faces_in_physical_group
            end
            ngroups += nldfaces
            if d == (D-1)
                #groups[d+1]["boundary"] = 1:length(dface_to_physical_group)
            end
        end # d
        ngroups += 1
        groups, face_to_nodes
    end # barrier
    chain = cartesian_chain(domain,cells_per_dir)
    interior_mesh = GT.mesh(chain)
    D = num_dims(interior_mesh)
    cell_to_nodes = face_nodes(interior_mesh,D)
    reference_cells = reference_spaces(interior_mesh,D)
    ref_cell = first(reference_cells)
    #refid_to_refface = reference_spaces(boundary(ref_cell))
    nnodes = num_nodes(interior_mesh)

    cell_geometry = unit_n_cube(Val(D))
    ref_simplex_mesh = simplexify(cell_geometry)
    d_to_ldface_to_sldface_to_lnodes = [
      [ face_nodes(ref_simplex_mesh,d)[group_faces(ref_simplex_mesh,d)["$d-face-$ldface"]]
      for ldface in 1:num_faces(complexify(ref_cell),d) ] for d in 0:(D-1)]
    d_to_ldface_to_lnodes = [face_nodes(complexify(ref_cell),d) for d in 0:(D-1)]
    groups, face_to_nodes = barrier(
      D,
      cell_to_nodes,
      nnodes,
      d_to_ldface_to_lnodes,
      d_to_ldface_to_sldface_to_lnodes
     )
    simplex_chain = structured_simplex_chain(domain,cells_per_dir)
    node_coords = node_coordinates(simplex_chain)
    face_to_refid = [ ones(Int8,length(face_to_nodes[d+1]))  for d in 0:(D-1)]
    mesh_face_nodes = push(face_to_nodes,face_nodes(simplex_chain,D))
    mesh_face_reference_id = push(face_to_refid,face_reference_id(simplex_chain,D))
    mesh_reference_faces = reference_spaces(ref_simplex_mesh)
    mesh_groups = push(groups,group_faces(simplex_chain,D))
    GT.mesh(
     node_coordinates = node_coords,
     face_nodes = mesh_face_nodes,
     face_reference_id = mesh_face_reference_id,
     reference_spaces = mesh_reference_faces,
     geometry_names = [ [ "$d-face-$face" for face in 1:num_faces(GT.mesh(cell_geometry),d) ] for d in 0:D],
     group_faces=mesh_groups,
    )
end

function cartesian_periodic_nodes(domain,cells,periodic)
    nodes = cells .+ 1
    nnodes = prod(nodes)
    node_owner = zeros(Int32,nnodes)
    node_owner .= 1:nnodes
    lis = LinearIndices(nodes)
    D = length(cells)
    for (d,flag) in enumerate(periodic)
        if ! flag
            continue
        end
        t = ntuple(Val(D)) do n
            if n != d
                Int(nodes[n])
            else
                1
            end
        end
        cis = CartesianIndices(t)
        for ci in cis
            owner = lis[ci]
            t2 = ntuple(Val(D)) do n
                if n != d
                    Int(Tuple(ci)[n])
                else
                    Int(nodes[d])
                end
            end
            ci2 = CartesianIndex(t2)
            node = lis[ci2]
            node_owner[node] = owner
        end
    end
    fold_periodic_nodes!(node_owner)
    node_owner
end


function moebius_strip(cells;complexify=Val(true))
    domain = (-1.,1.,0.0,2*pi)
    mesh0 = cartesian_mesh(domain,cells;complexify=Val(false),boundary=Val(false))
    #https://simple.wikipedia.org/wiki/M%C3%B6bius_strip
    function f(x0)
        r,θ = x0
        θ2 = 0.5*θ
        r2 = 0.5*0.5*r
        t = (1+r2*cos(θ2))
        x = t*cos(θ)
        y = t*sin(θ)
        z = r2*sin(θ2)
        SVector(x,y,z)
    end
    x0 = GT.node_coordinates(mesh0)
    nnodes0 = num_nodes(mesh0)
    coupling = collect(Int32,1:nnodes0)
    nodesx = cells[1] + 1
    s =  (nnodes0+1) .- (1:nodesx)
    m =  1:nodesx
    coupling[s] = m
    D = num_dims(mesh0)
    face_nodes = JaggedArray(GT.face_nodes(mesh0,D))
    face_nodes.data .= (i->coupling[i]).(face_nodes.data)
    x = f.(x0)
    x = x[1:(nnodes0-nodesx)]
    mesh1 = replace_node_coordinates(mesh0,x)
    if val_parameter(complexify)
        mesh2 = GT.complexify(mesh1)
    else
        mesh2 = mesh1
    end
end


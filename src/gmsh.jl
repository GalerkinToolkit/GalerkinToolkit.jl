const INVALID = 0

function default_gmsh_options()
    [
     "General.Terminal"=>1,
     "Mesh.SaveAll"=>1,
     "Mesh.MedImportGroupsOfNodes"=>1
    ]
end

function with_gmsh(f;options=default_gmsh_options())
    gmsh.initialize()
    for (k,v) in options
        gmsh.option.setNumber(k,v)
    end
    try
        return f()
    finally
        gmsh.finalize()
    end
end

function gmsh_mesh(file;renumber=true,kwargs...)
    @assert ispath(file) "File not found: $(file)"
    with_gmsh(;kwargs...) do
        gmsh.open(file)
        renumber && gmsh.model.mesh.renumberNodes()
        renumber && gmsh.model.mesh.renumberElements()
        mesh_from_gmsh_module()
    end
end

function mesh_from_gmsh_module()
    entities = gmsh.model.getEntities()
    nodeTags, coord, parametricCoord = gmsh.model.mesh.getNodes()

    # find num_dims
    ddim = -1
    for e in entities
        ddim = max(ddim,e[1])
    end
    if ddim == -1
        error("No entities in the msh file.")
    end
    D = ddim

    # find embedded_dimension
    dtouched = [false,false,false]
    for node in nodeTags
        if !(coord[(node-1)*3+1] + 1 ≈ 1)
            dtouched[1] = true
        end
        if !(coord[(node-1)*3+2] + 1 ≈ 1)
            dtouched[2] = true
        end
        if !(coord[(node-1)*3+3] + 1 ≈ 1)
            dtouched[3] = true
        end
    end
    if dtouched[3]
        adim = 3
    elseif dtouched[2]
        adim = 2
    elseif dtouched[1]
        adim = 1
    else
        adim = 0
    end

    # Setup node coords
    nmin = minimum(nodeTags)
    nmax = maximum(nodeTags)
    nnodes = length(nodeTags)
    if !(nmax == nnodes && nmin == 1)
        error("Only consecutive node tags allowed.")
    end
    my_node_to_coords = zeros(SVector{adim,Float64},nnodes)
    m = zero(MVector{adim,Float64})
    for node in nodeTags
        for j in 1:adim
            k = (node-1)*3 + j
            xj = coord[k]
            m[j] = xj
        end
        my_node_to_coords[node] = m
    end

    # Setup face nodes
    offsets = zeros(Int32,D+1)
    my_face_nodes = Vector{JaggedArray{Int32,Int32}}(undef,D+1)
    for d in 0:D
        elemTypes, elemTags, nodeTags = gmsh.model.mesh.getElements(d)
        ndfaces = 0
        for t in 1:length(elemTypes)
            ndfaces += length(elemTags[t])
        end
        if ndfaces != 0
            nmin::Int = minimum( minimum, elemTags )
            nmax::Int = maximum( maximum, elemTags )
            if !( (nmax-nmin+1) == ndfaces)
                error("Only consecutive elem tags allowed.")
            end
            offsets[d+1] = nmin-1
        end
        ptrs = zeros(Int32,ndfaces+1)
        dface = 0
        for t in 1:length(elemTypes)
            elementName, dim, order, numNodes, nodeCoord =
            gmsh.model.mesh.getElementProperties(elemTypes[t])
            for e in 1:length(elemTags[t])
                dface += 1
                ptrs[dface+1] = numNodes
            end
        end
        length_to_ptrs!(ptrs)
        ndata = ptrs[end]-1
        data = zeros(Int32,ndata)
        dface = 1
        for t in 1:length(elemTypes)
            p = ptrs[dface]-Int32(1)
            for (i,node) in enumerate(nodeTags[t])
                data[p+i] = node
            end
            dface += length(elemTags[t])
        end
        my_face_nodes[d+1] = JaggedArray(data,ptrs)
    end

    # Setup face_reference_id
    my_face_reference_id = Vector{Vector{Int32}}(undef,D+1)
    for d in 0:D
        elemTypes, elemTags, nodeTags = gmsh.model.mesh.getElements(d)
        ndfaces = length(my_face_nodes[d+1])
        dface_to_refid = zeros(Int8,ndfaces)
        refid = 0
        dface = 0
        for t in 1:length(elemTypes)
            refid += 1
            for e in 1:length(elemTags[t])
                dface += 1
                dface_to_refid[dface] = refid
            end
        end
        my_face_reference_id[d+1] = dface_to_refid
    end

    # Setup reference faces
    my_reference_faces = ()
    for d in D:-1:0
        elemTypes, elemTags, nodeTags = gmsh.model.mesh.getElements(d)
        refdfaces = ()
        for t in 1:length(elemTypes)
            refface = reference_face_from_gmsh_eltype(elemTypes[t])
            refdfaces = (refdfaces...,refface)
        end
        if refdfaces == ()
            refdfaces = reference_faces(first(first(my_reference_faces)),d)
        end
        my_reference_faces = (refdfaces,my_reference_faces...)
    end

    # Setup periodic nodes
    node_to_main_node = fill(Int32(INVALID),nnodes)
    for (dim,tag) in entities
        tagMaster, nodeTags, nodeTagsMaster, = gmsh.model.mesh.getPeriodicNodes(dim,tag)
        for i in 1:length(nodeTags)
            node = nodeTags[i]
            main_node = nodeTagsMaster[i]
            node_to_main_node[node] = main_node
        end
    end
    my_periodic_nodes = collect(Int32,findall(i->i!=INVALID,node_to_main_node))
    my_periodic_to_master = node_to_main_node[my_periodic_nodes]
    my_periodic_to_coeff = ones(length(my_periodic_to_master))

    # Setup physical groups
    T_group = typeof(physical_group(Int32[],""))
    my_groups = [ Dict{Int,T_group}() for d in 0:D]
    for d in 0:D
        offset = Int32(offsets[d+1])
        dimTags = gmsh.model.getPhysicalGroups(d)
        for (dim,tag) in dimTags
            @boundscheck @assert dim == d
            g_entities = gmsh.model.getEntitiesForPhysicalGroup(dim,tag)
            ndfaces_in_physical_group = 0
            for entity in g_entities
                elemTypes, elemTags, nodeTags = gmsh.model.mesh.getElements(dim,entity)
                for t in 1:length(elemTypes)
                    ndfaces_in_physical_group += length(elemTags[t])
                end
            end
            dfaces_in_physical_group = zeros(Int32,ndfaces_in_physical_group)
            ndfaces_in_physical_group = 0
            for entity in g_entities
                elemTypes, elemTags, nodeTags = gmsh.model.mesh.getElements(dim,entity)
                for t in 1:length(elemTypes)
                    for etag in elemTags[t]
                        ndfaces_in_physical_group += 1
                        dfaces_in_physical_group[ndfaces_in_physical_group] = Int32(etag)-offset
                    end
                end
            end
            groupname = gmsh.model.getPhysicalName(dim,tag)
            my_group = physical_group(dfaces_in_physical_group,groupname)
            my_groups[d+1][Int(tag)] = my_group
        end
    end

    mesh = new_mesh(my_node_to_coords,my_face_nodes,my_face_reference_id,my_reference_faces)
    if length(my_periodic_nodes) != 0
        constraints = periodic_node_constraints(my_periodic_nodes,my_periodic_to_master,my_periodic_to_master)
        mesh2 = set_periodic_node_constraints(mesh,constraints)
    else
        mesh2 = mesh
    end
    mesh3 = set_phyisical_groups(mesh2,my_groups)
    mesh3
end

function reference_face_from_gmsh_eltype(eltype)
    if eltype == 1
        Meshes.Segment(Meshes.Point(0),Meshes.Point(1))
    elseif eltype == 2
        Meshes.Triangle(Meshes.Point.([(0,0),(1,0),(0,1)]))
    elseif eltype == 3
        Meshes.Quadrangle(Meshes.Point.([(0,0),(1,0),(1,1),(0,1)]))
    elseif eltype == 4
        Meshes.Tetrahedron(Meshes.Point.([(0,0,0),(1,0,0),(0,1,0),(0,0,1)]))
    elseif eltype == 5
        Meshes.Hexahedron(Meshes.Point.([(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)]))
    elseif eltype == 15
        Meshes.Point(SVector{0,Float64}())
    elseif eltype == 8
        linear = reference_face_from_gmsh_eltype(1)
        order = 2
        GmshHighOrderSimplex(linear,order)
    elseif eltype == 9
        linear = reference_face_from_gmsh_eltype(2)
        order = 2
        GmshHighOrderSimplex(linear,order)
    else
        en, = gmsh.model.mesh.getElementProperties(eltype)
        error("Unsupported element type. elemType: $eltype ($en)")
    end
end

struct GmshHighOrderSimplex{A}
    linear_mesh::A
    order::Int
end

linear_mesh(a::GmshHighOrderSimplex) = a.linear_mesh
dimension(a::GmshHighOrderSimplex) = dimension(linear_mesh(a))
is_simplex(a::GmshHighOrderSimplex) = true

function node_coordinates(a::GmshHighOrderSimplex)
    linear = linear_mesh(a)
    x = node_coordinates(linear)
    D = dimension(a)
    if D == 1
        if a.order == 2
            [1.0*x[1],1.0*x[2],0.5*(x[1]+x[2])]
        else
            error("Order not yet implemented")
        end
    elseif D == 2
        if a.order == 2
            [1.0*x[1],1.0*x[2],1.0*x[3],0.5*(x[1]+x[2]),0.5*(x[2]+x[3]),0.5*(x[3]+x[1])]
        else
            error("Order not yet implemented")
        end
    else
        error("Dimension not yet implemented")
    end
end

face_reference_id(a::GmshHighOrderSimplex,d) = face_reference_id(linear_mesh(a),d)

function reference_faces(a::GmshHighOrderSimplex,::Val{d}) where d
    order = a.order
    map(reference_faces(linear_mesh(a),Val(d))) do l
        GmshHighOrderSimplex(l,order)
    end
end

function face_nodes(a::GmshHighOrderSimplex,d)
    D = dimension(a)
    if D == 1
        if a.order == 2
            d==0 && return JaggedArray{Int32,Int32}(Vector{Int32}[[1],[2]])
            d==1 && return JaggedArray{Int32,Int32}(Vector{Int32}[[1,2,3]])
        else
            error("Order not yet implemented")
        end
    elseif D == 2
        if a.order == 2
            d==0 && return JaggedArray{Int32,Int32}(Vector{Int32}[[1],[2],[3]])
            d==1 && return JaggedArray{Int32,Int32}(Vector{Int32}[[1,2,4],[2,3,5],[3,1,6]])
            d==2 && return JaggedArray{Int32,Int32}(Vector{Int32}[[1,2,3,4,5,6]])
        else
            error("Order not yet implemented")
        end
    else
        error("Dimension not yet implemented")
    end
end

function vtk_mesh_cell(a::GmshHighOrderSimplex)
    nodes -> WriteVTK.MeshCell(vtk_cell_type(a),nodes)
end

function vtk_cell_type(a::GmshHighOrderSimplex)
    D = dimension(a)
    if D == 1
        if a.order == 2
            WriteVTK.VTKCellTypes.VTK_QUADRATIC_EDGE
        else
            error("Order not yet implemented")
        end
    elseif D == 2
        if a.order == 2
            WriteVTK.VTKCellTypes.VTK_QUADRATIC_TRIANGLE
        else
            error("Order not yet implemented")
        end
    else
        error("Dimension not yet implemented")
    end
end



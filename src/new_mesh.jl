

num_dims(m::AbstractMesh) = length(reference_spaces(m))-1
num_ambient_dims(m::AbstractMesh) = length(eltype(node_coordinates(m)))
options(m::AbstractMesh) = options(first(last(reference_spaces(m))))
num_faces(m::AbstractMesh,d) = length(face_reference_id(m,d))
num_nodes(fe::AbstractMesh) = length(node_coordinates(fe))
num_faces(mesh::AbstractMesh) = map(length,face_reference_id(mesh))

function label_faces_in_dim!(m::AbstractMesh,d;physical_name="__$d-FACES__")
    groups = physical_faces(m,d)
    if haskey(groups,physical_name)
        return physical_name
    end
    Ti = int_type(options(m))
    faces = collect(Ti,1:num_faces(m,d))
    groups[physical_name] = faces
    physical_name
end

function label_interior_faces!(mesh::AbstractMesh;physical_name="__INTERIOR_FACES__")
    D = num_dims(mesh)
    d = D-1
    groups = physical_faces(mesh,d)
    if haskey(groups,physical_name)
        return physical_name
    end
    topo = topology(mesh)
    face_to_cells = face_incidence(topo,d,D)
    faces = findall(cells->length(cells)==2,face_to_cells)
    groups[physical_name] = faces
    physical_name
end

function label_boundary_faces!(mesh::AbstractMesh;physical_name="__BOUNDARY_FACES__")
    D = num_dims(mesh)
    d = D-1
    groups = physical_faces(mesh,d)
    if haskey(groups,physical_name)
        return physical_name
    end
    topo = topology(mesh)
    face_to_cells = face_incidence(topo,d,D)
    faces = findall(cells->length(cells)==1,face_to_cells)
    groups[physical_name] = faces
    physical_name
end

function domain(mesh::AbstractMesh,d;
    mesh_id = objectid(mesh),
    face_around=nothing,
    is_reference_domain=Val(false)
    )
    physical_name = label_faces_in_dim!(mesh,d)
    mesh_domain(;
        mesh,
        mesh_id,
        num_dims=Val(val_parameter(d)),
        physical_names=[physical_name],
        is_reference_domain)
end

function interior(mesh::AbstractMesh;
    mesh_id = objectid(mesh),
    physical_names=[label_faces_in_dim!(mesh,num_dims(mesh))],
    is_reference_domain=Val(false)
    )
    d = num_dims(mesh)
    mesh_domain(;
        mesh,
        mesh_id,
        physical_names,
        num_dims=Val(val_parameter(d)),
        is_reference_domain)
end

function skeleton(mesh::AbstractMesh;
    mesh_id = objectid(mesh),
    physical_names=[label_interior_faces!(mesh)],
    is_reference_domain=Val(false)
    )
    d = num_dims(mesh) - 1
    mesh_domain(;
        mesh,
        mesh_id,
        physical_names,
        num_dims=Val(val_parameter(d)),
        is_reference_domain)
end

function boundary(mesh::AbstractMesh;
    mesh_id = objectid(mesh),
    physical_names=[label_boundary_faces!(mesh)],
    is_reference_domain=Val(false)
    )
    d = num_dims(mesh) - 1
    mesh_domain(;
        mesh,
        mesh_id,
        physical_names,
        num_dims=Val(val_parameter(d)),
        is_reference_domain)
end

function reference_domains(a::AbstractMesh,d)
    map(domain,reference_spaces(a,d))
end

function reference_domains(a::AbstractMesh)
    D = num_dims(a)
    ntuple(d->reference_domains(a,d-1),Val(D+1))
end

function remove_interior(mesh::AbstractMesh)
    GT.mesh(;
         node_coordinates = node_coordinates(mesh),
         face_nodes = face_nodes(mesh)[1:end-1],
         face_reference_id = face_reference_id(mesh)[1:end-1],
         reference_spaces = reference_spaces(mesh)[1:end-1],
         physical_faces = physical_faces(mesh)[1:end-1],
         outward_normals = outward_normals(mesh),
         is_cell_complex = Val(is_cell_complex(mesh)),
         periodic_nodes = periodic_nodes(mesh)
        )
end

function simplexify(mesh::AbstractMesh;glue=Val(false))

    # TODO add attributes to mesh to figure out if we really
    # need to simplexify and/or complexify
    #mesh, = complexify(mesh)

    D = num_dims(mesh)
    refid_to_refcell = reference_spaces(mesh,D)
    #@assert length(refid_to_refcell) == 1
    refcell = first(refid_to_refcell)
    #@assert order(refcell) == 1
    ltmesh = simplexify(domain(refcell))

    ltcell_to_lnodes = face_nodes(ltmesh,D)
    cell_to_nodes = face_nodes(mesh,D)

    ncells = length(cell_to_nodes)
    nltcells = length(ltcell_to_lnodes)
    ntcells = ncells * nltcells
    tcell_to_nodes_ptrs = zeros(Int32,ntcells+1)
    tcell_to_cell = zeros(Int32,ntcells)

    tcell = 0
    for cell in 1:ncells
        for lnodes in ltcell_to_lnodes
            tcell +=1
            tcell_to_nodes_ptrs[tcell+1] = length(lnodes)
            tcell_to_cell[tcell] = cell
        end
    end
    length_to_ptrs!(tcell_to_nodes_ptrs)
    ndata = tcell_to_nodes_ptrs[end]-1
    T = eltype(eltype(cell_to_nodes))
    tcell_to_nodes_data = zeros(T,ndata)
    k = 1
    for cell in 1:ncells
        nodes = cell_to_nodes[cell]
        for lnodes in ltcell_to_lnodes
            for lnode in lnodes
                node = nodes[lnode]
                tcell_to_nodes_data[k] = node
                k += 1
            end
        end
    end
    tcell_to_nodes = JaggedArray(tcell_to_nodes_data,tcell_to_nodes_ptrs)
    reftcell = first(reference_spaces(ltmesh,D))
    tcell_to_refid = fill(Int8(1),ntcells)
    refid_to_reftcell = [reftcell]
    node_to_coords = node_coordinates(mesh)

    tchain = chain(;
        node_coordinates = node_to_coords,
        face_nodes = tcell_to_nodes,
        face_reference_id = tcell_to_refid,
        reference_spaces = refid_to_reftcell,)
    tmesh = complexify(GT.mesh(tchain))
    topo = topology(mesh)
    ttopo = topology(tmesh)
    ltopo = topology(ltmesh)
    d_to_tface_to_face = map(0:(D-1)) do d
        cell_to_faces = face_incidence(topo,D,d)
        tcell_to_tfaces = JaggedArray(face_incidence(ttopo,D,d))
        lface_to_lnodes = face_nodes(complexify(refcell),d)
        ltface_to_ltnodes = face_nodes(complexify(reftcell),d)
        if length(tcell_to_tfaces.data) != 0
            ntfaces = maximum(tcell_to_tfaces.data)
        else
            ntfaces = 0
        end
        tface_to_face = simplexify_generate_tface_to_face(
                                                          cell_to_faces,
                                                          tcell_to_tfaces,
                                                          ltcell_to_lnodes,
                                                          ltface_to_ltnodes,
                                                          lface_to_lnodes,
                                                          ntfaces)

    end
    push!(d_to_tface_to_face,tcell_to_cell)
    for d in 0:D
        groups = physical_faces(mesh,d)
        tgroups = physical_faces(tmesh,d)
        nfaces = num_faces(mesh,d)
        ntfaces = num_faces(tmesh,d)
        face_to_mask = fill(false,nfaces)
        tface_to_mask = fill(false,ntfaces)
        tface_to_face = d_to_tface_to_face[d+1]
        for (k,v) in groups
            fill!(face_to_mask,false)
            fill!(tface_to_mask,false)
            face_to_mask[v] .= true
            for tface in 1:ntfaces
                face = tface_to_face[tface]
                if face != 0
                    tface_to_mask[tface] = face_to_mask[face]
                end
            end
            tgroups[k] = findall(tface_to_mask)
        end
    end
    if val_parameter(glue)
        tmesh, d_to_tface_to_face
    else
        tmesh
    end
end

function simplexify_generate_tface_to_face(
        cell_to_faces,
        tcell_to_tfaces,
        ltcell_to_lnodes,
        ltface_to_ltnodes,
        lface_to_lnodes,
        ntfaces)

    nltcells = length(ltcell_to_lnodes)
    nltfaces = length(ltface_to_ltnodes)
    nlfaces = length(lface_to_lnodes)
    ltcell_ltface_to_lface = [ zeros(Int32,nltfaces) for ltcell in 1:nltcells  ]
    for ltcell in 1:nltcells
        ltnode_to_lnode = ltcell_to_lnodes[ltcell]
        for ltface in 1:nltfaces
            ltnodes = ltface_to_ltnodes[ltface]
            for lface in 1:nlfaces
                lnodes = lface_to_lnodes[lface]
                allin = true
                for ltnode in ltnodes
                    lnode = ltnode_to_lnode[ltnode]
                    if !(lnode in lnodes)
                        allin = false
                        break
                    end
                end
                if allin
                    ltcell_ltface_to_lface[ltcell][ltface] = lface
                    break
                end
            end
        end
    end
    ltcell_to_lfaces = ltcell_ltface_to_lface
    tface_to_face = zeros(Int32,ntfaces)
    ncells = length(cell_to_faces)
    nltcells = length(ltcell_to_lfaces)
    tcell = 1
    for cell in 1:ncells
        faces = cell_to_faces[cell]
        for ltcell in 1:nltcells
            tfaces = tcell_to_tfaces[tcell]
            nltfaces = length(tfaces)
            ltface_to_lface = ltcell_to_lfaces[ltcell]
            for ltface in 1:nltfaces
                tface = tfaces[ltface]
                lface = ltface_to_lface[ltface]
                if lface != 0
                    face = faces[lface]
                    tface_to_face[tface] = face
                end
            end
            tcell += 1
        end
    end
    tface_to_face
end

function restrict(mesh::AbstractMesh,args...)
    restrict_mesh(mesh,args...)
end

function restrict_mesh(mesh,lnode_to_node,lface_to_face_mesh)
    nnodes = num_nodes(mesh)
    node_to_lnode = zeros(Int32,nnodes)
    node_to_lnode[lnode_to_node] = 1:length(lnode_to_node)
    lnode_to_coords = node_coordinates(mesh)[lnode_to_node]
    lface_to_lnodes_mesh = map(lface_to_face_mesh,face_nodes(mesh)) do lface_to_face,face_to_nodes
        lface_to_nodes = view(face_to_nodes,lface_to_face)
        lface_to_lnodes = JaggedArray(lface_to_nodes)
        f = node->node_to_lnode[node]
        lface_to_lnodes.data .= f.(lface_to_lnodes.data)
        lface_to_lnodes
    end
    lface_to_refid_mesh = map((a,b)->b[a],lface_to_face_mesh,face_reference_id(mesh))
    D = num_dims(mesh)
    lgroups_mesh = map(lface_to_face_mesh,num_faces(mesh),physical_faces(mesh)) do lface_to_face, nfaces, groups
        lgroups = Dict{String,Vector{Int32}}()
        face_to_lface = zeros(Int32,nfaces)
        face_to_lface[lface_to_face] = 1:length(lface_to_face)
        for (k,faces) in groups
            lgroups[k] = filter(i->i!=0,face_to_lface[faces])
        end
        lgroups
    end
    pnode_to_node,pnode_to_master = periodic_nodes(mesh)
    plnode_to_lnode = filter(i->i!=0,node_to_lnode[pnode_to_node])
    plnode_to_lmaster = filter(i->i!=0,node_to_lnode[pnode_to_master])
    if outward_normals(mesh) !== nothing
        lnormals = outward_normals(mesh)[lface_to_face_mesh[end]]
    else
        lnormals = nothing
    end

    lmesh = GT.mesh(;
        node_coordinates = lnode_to_coords,
        face_nodes = lface_to_lnodes_mesh,
        face_reference_id = lface_to_refid_mesh,
        reference_spaces = reference_spaces(mesh),
        physical_faces = lgroups_mesh,
        periodic_nodes = (plnode_to_lnode=>plnode_to_lmaster),
        outward_normals = lnormals
        )

    lmesh
end

struct Mesh{A} <: AbstractMesh
    contents::A
end

function mesh(;
        node_coordinates,
        face_nodes,
        face_reference_id,
        reference_spaces,
        periodic_nodes = default_periodic_nodes(reference_spaces),
        physical_faces = default_physical_faces(reference_spaces),
        outward_normals = nothing,
        is_cell_complex = Val(false),
        workspace = nothing,
    )
    contents = (;
                node_coordinates,
                face_nodes,
                face_reference_id,
                reference_spaces,
                periodic_nodes,
                physical_faces,
                outward_normals,
                is_cell_complex,
                workspace,
               )
    mesh = Mesh(contents)
end

function replace_workspace(mesh::Mesh,workspace)
    contents = (;
                node_coordinates=node_coordinates(mesh),
                face_nodes=face_nodes(mesh),
                face_reference_id=face_reference_id(mesh),
                reference_spaces=reference_spaces(mesh),
                periodic_nodes=periodic_nodes(mesh),
                physical_faces=physical_faces(mesh),
                outward_normals=outward_normals(mesh),
                is_cell_complex=Val(is_cell_complex(mesh)),
                workspace,
               )
    Mesh(contents)
end

node_coordinates(m::Mesh) = m.contents.node_coordinates
face_nodes(m::Mesh) = m.contents.face_nodes
face_nodes(m::Mesh,d) = m.contents.face_nodes[d+1]
face_reference_id(m::Mesh) = m.contents.face_reference_id
face_reference_id(m::Mesh,d) = m.contents.face_reference_id[d+1]
reference_spaces(m::Mesh) = m.contents.reference_spaces
reference_spaces(m::Mesh,d) = m.contents.reference_spaces[d+1]
physical_faces(m::Mesh) = m.contents.physical_faces
physical_faces(m::Mesh,d) = m.contents.physical_faces[d+1]
periodic_nodes(m::Mesh) = m.contents.periodic_nodes
is_cell_complex(m::Mesh) = val_parameter(m.contents.is_cell_complex)
workspace(m::Mesh) = m.contents.workspace

function default_physical_faces(reference_spaces)
    [ Dict{String,Vector{int_type(options(first(last(reference_spaces))))}}() for _ in 1:length(reference_spaces) ]
end

function default_periodic_nodes(reference_spaces)
    Ti = int_type(options(first(last(reference_spaces))))
    Ti[] => Ti[]
end

function outward_normals(m::Mesh)
    m.contents.outward_normals
end

"""
"""
function physical_names(mesh,d)
    groups = physical_faces(mesh,d)
    Set(keys(groups))
end

function physical_names(mesh;merge_dims=Val(false))
    D = num_dims(mesh)
    d_to_names = [ physical_names(mesh,d) for d in 0:D]
    if val_parameter(merge_dims) == false
        return d_to_names
    end
    reduce(union,d_to_names)
end

num_dims(m::AbstractChain) = num_dims(domain(first(reference_spaces(m))))
num_ambient_dims(m::AbstractChain) = length(eltype(node_coordinates(m)))
options(m::AbstractChain) = options(first(reference_spaces(m)))
num_faces(m::AbstractChain) = length(face_reference_id(m))

"""
"""
function mesh(chain::AbstractChain)
    D = num_dims(chain)
    cell_nodes = face_nodes(chain)
    cell_reference_id = face_reference_id(chain)
    reference_cells = reference_spaces(chain)
    node_coords = node_coordinates(chain)
    face_to_nodes = Vector{typeof(cell_nodes)}(undef,D+1)
    face_to_refid = Vector{typeof(cell_reference_id)}(undef,D+1)
    for d in 0:D-1
        face_to_nodes[d+1] = Vector{Int}[]
        face_to_refid[d+1] = Int[]
    end
    face_to_nodes[end] = cell_nodes
    face_to_refid[end] = cell_reference_id
    ref_cell = first(reference_cells)
    ref_faces = reference_spaces(complexify(ref_cell))
    refid_to_refface = push(ref_faces[1:end-1],reference_cells)
    cell_groups = physical_faces(chain)
    groups = [ typeof(cell_groups)() for d in 0:D]
    groups[end] = cell_groups
    pnodes = periodic_nodes(chain)
    onormals = outward_normals(chain)
    GT.mesh(;
      node_coordinates = node_coords,
      face_nodes = face_to_nodes,
      face_reference_id = face_to_refid,
      reference_spaces = refid_to_refface,
      periodic_nodes = pnodes,
      physical_faces = groups,
      outward_normals = onormals)
end

struct Chain{A} <: AbstractChain
    contents::A
end

function chain(;
        node_coordinates,
        face_nodes,
        face_reference_id,
        reference_spaces,
        periodic_nodes = default_periodic_nodes((reference_spaces,)),
        physical_faces = default_physical_faces((reference_spaces,))[end],
        outward_normals = nothing,
    )
    contents = (;
                node_coordinates,
                face_nodes,
                face_reference_id,
                reference_spaces,
                periodic_nodes,
                physical_faces,
                outward_normals,
               )
    Chain(contents)
end

node_coordinates(m::Chain) = m.contents.node_coordinates
face_nodes(m::Chain) = m.contents.face_nodes
face_reference_id(m::Chain) = m.contents.face_reference_id
reference_spaces(m::Chain) = m.contents.reference_spaces
physical_faces(m::Chain) = m.contents.physical_faces
periodic_nodes(m::Chain) = m.contents.periodic_nodes
outward_normals(m::Chain) = m.contents.outward_normals

function chain(mesh::AbstractMesh,D=Val(num_dims(mesh)))
    d = val_parameter(D)
    chain(;
          node_coordinates=node_coordinates(mesh),
          face_nodes=face_nodes(mesh,d),
          face_reference_id=face_reference_id(mesh,d),
          reference_spaces=reference_spaces(mesh,d),
          periodic_nodes=periodic_nodes(mesh),
          physical_faces=physical_faces(mesh,d),
          outward_normals=outward_normals(mesh),
         )
end



num_dims(m::AbstractMesh) = length(reference_spaces(m))-1
num_ambient_dims(m::AbstractMesh) = length(eltype(node_coordinates(m)))
options(m::AbstractMesh) = options(first(last(reference_spaces(m))))
num_faces(m::AbstractMesh,d) = length(face_reference_id(m,d))
num_nodes(fe::AbstractMesh) = length(node_coordinates(fe))
num_faces(mesh::AbstractMesh) = map(length,face_reference_id(mesh))
face_nodes(m::AbstractMesh,d) = face_nodes(m)[d+1]
face_reference_id(m::AbstractMesh,d) = face_reference_id(m)[d+1]
reference_spaces(m::AbstractMesh,d) = reference_spaces(m)[d+1]
group_faces(m::AbstractMesh,d) = group_faces(m)[d+1]
geometry_names(m::AbstractMesh,d) = geometry_names(m)[d+1]

function geometries(m::AbstractMesh,d)
    map(name->domain(m,d;group_names=[name]),geometry_names(m,d))
end

function face_offset(a)
    D = num_dims(a)
    offsets = zeros(Int,D+1)
    for d in 1:D
        offsets[d+1] = offsets[d] + num_faces(a,d-1)
    end
    offsets
end
function face_offset(a,d)
    face_offset(a)[d+1]
end
function face_range(a,d)
    o = face_offset(a,d)
    o .+ (1:num_faces(a,d))
end
function face_range(a)
    D = num_dims(a)
    map(d->face_range(a,d),0:D)
end
#function face_dim(a,d)
#    n = num_faces(a,d)
#    fill(d,n)
#end
#function face_dim(a)
#    D = num_dims(a)
#    reduce(vcat,map(d->face_dim(a,d),0:D))
#end

function label_faces_in_dim!(m::AbstractMesh,d;physical_name="__$d-FACES__")
    groups = group_faces(m,d)
    if haskey(groups,physical_name)
        return physical_name
    end
    Ti = int_type(options(m))
    faces = collect(Ti,1:num_faces(m,d))
    groups[physical_name] = faces
    physical_name
end

function label_faces_in_dim!(m::AbstractPMesh,d;physical_name="__$d-FACES__")
    p_mesh = partition(m)
    foreach(p_mesh) do mesh
        label_faces_in_dim!(mesh,d;physical_name)
    end
    physical_name
end

function label_interior_faces!(mesh::AbstractMesh;physical_name="__INTERIOR_FACES__")
    D = num_dims(mesh)
    d = D-1
    groups = group_faces(mesh,d)
    if haskey(groups,physical_name)
        return physical_name
    end
    topo = topology(mesh)
    face_to_cells = face_incidence(topo,d,D)
    faces = findall(cells->length(cells)==2,face_to_cells)
    groups[physical_name] = faces
    physical_name
end

function label_interior_faces!(pmesh::AbstractPMesh;physical_name="__INTERIOR_FACES__")
    D = num_dims(pmesh)
    d = D - 1
    vals = map(partition(pmesh)) do mesh
        topo = topology(mesh)
        nfaces = num_faces(topo,d)
        face_to_v = zeros(Int32,nfaces)
        cell_to_faces = face_incidence(topo,D,d)
        ncells = num_faces(topo,D)
        cell_ids = face_local_indices(mesh,D)
        cell_to_owner = local_to_owner(cell_ids)
        part = part_id(cell_ids)
        for cell in 1:ncells
            owner = cell_to_owner[cell]
            if part != owner
                continue
            end
            faces = cell_to_faces[cell]
            for face in faces
                face_to_v[face] += 1
            end
        end
        face_to_v
    end
    ids = face_partition(pmesh, d)
    v = PVector(vals,ids)
    assemble!(v) |> wait
    map(partition(pmesh),vals) do mesh, face_to_v
       group_faces(mesh,d)[physical_name] = findall(i->i==2,face_to_v)
    end
    physical_name
end

function label_boundary_faces!(mesh::AbstractMesh;physical_name="__BOUNDARY_FACES__")
    D = num_dims(mesh)
    d = D-1
    groups = group_faces(mesh,d)
    if haskey(groups,physical_name)
        return physical_name
    end
    topo = topology(mesh)
    face_to_cells = face_incidence(topo,d,D)
    faces = findall(cells->length(cells)==1,face_to_cells)
    groups[physical_name] = faces
    physical_name
end

function label_boundary_faces!(domain::AbstractDomain;physical_name="__BOUNDARY_$(objectid(domain))__")
    mesh = GT.mesh(domain)
    D = num_dims(mesh)
    d = D-1
    groups = group_faces(mesh,d)
    if haskey(groups,physical_name)
        return physical_name
    end
    topo = topology(mesh)
    cell_to_faces = face_incidence(topo,D,d)
    nfaces = num_faces(mesh,d)
    face_count = zeros(Int32,nfaces)
    for faces in cell_to_faces
        for face in faces
            face_count[face] += 1
        end
    end
    faces = findall(count->count==1,face_count)
    groups[physical_name] = faces
    physical_name
end

function label_boundary_faces!(pmesh::AbstractPMesh;physical_name="__BOUNDARY_FACES__")
    D = num_dims(pmesh)
    d = D - 1
    vals = map(partition(pmesh)) do mesh
        topo = topology(mesh)
        nfaces = num_faces(topo,d)
        face_to_v = zeros(Int32,nfaces)
        cell_to_faces = face_incidence(topo,D,d)
        ncells = num_faces(topo,D)
        cell_ids = face_local_indices(mesh,D)
        cell_to_owner = local_to_owner(cell_ids)
        part = part_id(cell_ids)
        for cell in 1:ncells
            owner = cell_to_owner[cell]
            if part != owner
                continue
            end
            faces = cell_to_faces[cell]
            for face in faces
                face_to_v[face] += 1
            end
        end
        face_to_v
    end
    ids = face_partition(pmesh, d)
    v = PVector(vals,ids)
    assemble!(v) |> wait
    map(partition(pmesh),vals) do mesh, face_to_v
       group_faces(mesh,d)[physical_name] = findall(i->i==1,face_to_v)
    end
    physical_name
end

"""
"""
function domain(mesh::AbstractMesh,d;
    mesh_id = objectid(mesh),
    face_around=nothing,
    is_reference_domain=Val(false),
    group_names=[label_faces_in_dim!(mesh,val_parameter(d))],
    )
    mesh_domain(;
        mesh,
        mesh_id,
        num_dims=Val(val_parameter(d)),
        group_names,
        is_reference_domain)
end

function interior(mesh::AbstractMesh;
    mesh_id = objectid(mesh),
    group_names=[label_faces_in_dim!(mesh,num_dims(mesh))],
    is_reference_domain=Val(false)
    )
    d = num_dims(mesh)
    mesh_domain(;
        mesh,
        mesh_id,
        group_names,
        num_dims=Val(val_parameter(d)),
        is_reference_domain)
end

function skeleton(mesh::AbstractMesh;
    mesh_id = objectid(mesh),
    group_names=[label_interior_faces!(mesh)],
    is_reference_domain=Val(false)
    )
    d = num_dims(mesh) - 1
    mesh_domain(;
        mesh,
        mesh_id,
        group_names,
        num_dims=Val(val_parameter(d)),
        is_reference_domain)
end

function boundary(mesh::AbstractMesh;
    mesh_id = objectid(mesh),
    group_names=[label_boundary_faces!(mesh)],
    is_reference_domain=Val(false),
    face_around = 1,
    )
    d = num_dims(mesh) - 1
    mesh_domain(;
        mesh,
        mesh_id,
        group_names,
        face_around,
        num_dims=Val(val_parameter(d)),
        is_reference_domain)
end

function boundary(domain::AbstractDomain;
    mesh_id = objectid(GT.mesh(domain)),
    group_names=[label_boundary_faces!(domain)],
    is_reference_domain=Val(false),
    face_around = 1,
    )
    mesh = GT.mesh(domain)
    d = num_dims(mesh) - 1
    mesh_domain(;
        mesh,
        mesh_id,
        group_names,
        face_around,
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
         group_faces = group_faces(mesh)[1:end-1],
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
        groups = group_faces(mesh,d)
        tgroups = group_faces(tmesh,d)
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

function restrict(mesh::AbstractMesh,args...;kwargs...)
    restrict_mesh(mesh,args...;kwargs...)
end

function restrict_mesh(mesh,lnode_to_node,lface_to_face_mesh;kwargs...)
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
    lgroups_mesh = map(lface_to_face_mesh,num_faces(mesh),group_faces(mesh)) do lface_to_face, nfaces, groups
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
        group_faces = lgroups_mesh,
        periodic_nodes = (plnode_to_lnode=>plnode_to_lmaster),
        outward_normals = lnormals,
        kwargs...
        )

    lmesh
end

struct Mesh{A} <: AbstractMesh
    contents::A
end

"""
    create_mesh(;kwargs...)

Build an arbitrary mesh object.

See also [`cartesian_mesh`](@ref) and [`mesh_from_msh`](@ref).

# Level

Intermediate

# Keyword arguments

- `node_coordinates`: The vector containing the coordinates of all mesh nodes. `node_coordinates[i]` is the coordinate vector for global node number `i`.
- `face_nodes`: A highly-nested vector containing the node ids for each face in the mesh. `node_coordinates[n]` with `n=face_nodes[d+1][i][k]` is the global node coordinate for local node number `k` in face `i` of dimension `d`. The object `face_nodes[d+1]` is a long vector of small vectors of integers. It is often represented using a `JaggedArray` object that uses continuous linear memory for performance.
- `reference_spaces`: A nested tuple containing the reference spaces for faces. `reference_spaces[d+1][i]` is the reference space number `i` of dimension `d`.
- `face_reference_id` [optional]: A nested vector containing which reference space is assigned to each face. `reference_sapces[d+1][r]` with `r=face_reference_id[d+1][i]` is the reference space associated with face number `i` of dimension `d`. By default, all faces are assigned to the first reference space in its dimension.
- `group_faces` [optional]: A vector of dictionaries containing groups labeled groups of faces. `group_faces[d+1][label]` is a vector of integers containing the ids  of the faces labeled as `label` in dimension `d`. These labels might overlap. By default, no faces groups are created.
- `is_cell_complex=Val(false)` [optional]: `Val(true)` if the input data represents a cell complex, `Val(false)` otherwise.
- `outward_normals=nothing` [optinal]: Vector containing the normal vectors for the faces of maximum dimension of the mesh. This is relevant for meshes of dimension `d` embedded in `d+1` dimensions as there is no way to tell which should be the orientation of the normals from the other quantities defining the mesh.  `outward_normals[f]` gives the normal vector of face number `f` of dimension `d=length(face_nodes)-1`.
"""
function create_mesh end


function create_mesh(;kwargs...)
    mesh(;kwargs...)
end

function default_face_reference_id(face_nodes,reference_spaces)
    D = length(face_nodes)-1
    [ fill(Int8(1),length(face_nodes[i+1]))   for i in 0:D ] 
end

function mesh(;
        node_coordinates,
        face_nodes,
        face_reference_id = default_face_reference_id(face_nodes,reference_spaces),
        reference_spaces,
        periodic_nodes = default_periodic_nodes(reference_spaces),
        group_faces = default_group_faces(reference_spaces),
        outward_normals = nothing,
        geometry_names = [ String[] for d in 1:length(face_reference_id)],
        is_cell_complex = Val(false),
        node_local_indices = PartitionedArrays.block_with_constant_size(1,(1,),(length(node_coordinates),)),
        face_local_indices = [ PartitionedArrays.block_with_constant_size(1,(1,),(length(face_reference_id[d]),)) for d in 1:length(face_reference_id)],
        workspace = nothing,
    )
    contents = (;
                node_coordinates,
                face_nodes,
                face_reference_id,
                reference_spaces,
                periodic_nodes,
                group_faces,
                outward_normals,
                geometry_names,
                is_cell_complex,
                node_local_indices,
                face_local_indices,
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
                geometry_names=geometry_names(mesh),
                group_faces=group_faces(mesh),
                outward_normals=outward_normals(mesh),
                is_cell_complex=Val(is_cell_complex(mesh)),
                node_local_indices=node_local_indices(mesh),
                face_local_indices=face_local_indices(mesh),
                workspace,
               )
    Mesh(contents)
end

function replace_node_coordinates(mesh::Mesh,node_coordinates)
    contents = (;
                node_coordinates,
                face_nodes=face_nodes(mesh),
                face_reference_id=face_reference_id(mesh),
                reference_spaces=reference_spaces(mesh),
                periodic_nodes=periodic_nodes(mesh),
                group_faces=group_faces(mesh),
                geometry_names=geometry_names(mesh),
                outward_normals=outward_normals(mesh),
                is_cell_complex=Val(is_cell_complex(mesh)),
                node_local_indices=node_local_indices(mesh),
                face_local_indices=face_local_indices(mesh),
                workspace=workspace(mesh),
               )
    Mesh(contents)
end

node_coordinates(m::Mesh) = m.contents.node_coordinates
face_nodes(m::Mesh) = m.contents.face_nodes
face_nodes(m::Mesh,d) = m.contents.face_nodes[d+1]
face_reference_id(m::Mesh) = m.contents.face_reference_id
face_reference_id(m::Mesh,d) = m.contents.face_reference_id[d+1]
reference_spaces(m::Mesh) = m.contents.reference_spaces
reference_spaces(m::Mesh,d) = m.contents.reference_spaces[val_parameter(d)+1]
group_faces(m::Mesh) = m.contents.group_faces
group_faces(m::Mesh,d) = m.contents.group_faces[d+1]
periodic_nodes(m::Mesh) = m.contents.periodic_nodes
is_cell_complex(m::Mesh) = val_parameter(m.contents.is_cell_complex)
workspace(m::Mesh) = m.contents.workspace

node_local_indices(m::Mesh) = m.contents.node_local_indices
face_local_indices(m::Mesh) = m.contents.face_local_indices
face_local_indices(m::Mesh,d) = m.contents.face_local_indices[d+1]

function default_group_faces(reference_spaces)
    [ Dict{String,Vector{int_type(options(first(last(reference_spaces))))}}() for _ in 1:length(reference_spaces) ]
end

function default_periodic_nodes(reference_spaces)
    Ti = int_type(options(first(last(reference_spaces))))
    Ti[] => Ti[]
end

function outward_normals(m::Mesh)
    m.contents.outward_normals
end

function group_names(mesh,d)
    groups = group_faces(mesh,d)
    Set(keys(groups))
end

function group_names(mesh;merge_dims=Val(false))
    D = num_dims(mesh)
    d_to_names = [ group_names(mesh,d) for d in 0:D]
    if val_parameter(merge_dims) == false
        return d_to_names
    end
    reduce(union,d_to_names)
end

function geometry_names(a::Mesh)
    a.contents.geometry_names
end

struct Chain{A} <: AbstractMesh
    contents::A
end

"""
"""
function chain(;
        node_coordinates,
        face_nodes,
        face_reference_id,
        reference_spaces,
        periodic_nodes = default_periodic_nodes((reference_spaces,)),
        group_faces = default_group_faces((reference_spaces,))[end],
        outward_normals = nothing,
    )
    contents = (;
                node_coordinates,
                face_nodes,
                face_reference_id,
                reference_spaces,
                periodic_nodes,
                group_faces,
                outward_normals,
               )
    Chain(contents)
end

node_coordinates(m::Chain) = m.contents.node_coordinates

function face_nodes(m::Chain)
    D = num_dims(m)
    T = typeof(m.contents.face_nodes)
    d_face_nodes = Vector{T}(undef,D+1)
    for d in 0:D
        d_face_nodes[d+1] = Vector{Int}[]
    end
    d_face_nodes[end] = m.contents.face_nodes
    d_face_nodes
end

function face_reference_id(m::Chain)
    D = num_dims(m)
    T = typeof(m.contents.face_reference_id)
    face_to_refid = Vector{T}(undef,D+1)
    for d in 0:D
        face_to_refid[d+1] = Int[]
    end
    face_to_refid[end] = m.contents.face_reference_id
    face_to_refid
end

function reference_spaces(m::Chain)
    reference_cells = m.contents.reference_spaces
    ref_cell = first(reference_cells)
    ref_faces = reference_spaces(complexify(ref_cell))
    refid_to_refface = push(ref_faces[1:end-1],reference_cells)
end

function group_faces(m::Chain)
    D = num_dims(m)
    cell_groups = m.contents.group_faces
    groups = [ typeof(cell_groups)() for d in 0:D]
    groups[end] = cell_groups
    groups
end

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
          group_faces=group_faces(mesh,d),
          outward_normals=outward_normals(mesh),
         )
end

"""
"""
function mesh(chain::Chain)
    #D = num_dims(chain)
    #cell_nodes = face_nodes(chain)
    #cell_reference_id = face_reference_id(chain)
    #reference_cells = reference_spaces(chain)
    #node_coords = node_coordinates(chain)
    #face_to_nodes = Vector{typeof(cell_nodes)}(undef,D+1)
    #face_to_refid = Vector{typeof(cell_reference_id)}(undef,D+1)
    #for d in 0:D-1
    #    face_to_nodes[d+1] = Vector{Int}[]
    #    face_to_refid[d+1] = Int[]
    #end
    #face_to_nodes[end] = cell_nodes
    #face_to_refid[end] = cell_reference_id
    #ref_cell = first(reference_cells)
    #ref_faces = reference_spaces(complexify(ref_cell))
    #refid_to_refface = push(ref_faces[1:end-1],reference_cells)
    #cell_groups = group_faces(chain)
    #groups = [ typeof(cell_groups)() for d in 0:D]
    #groups[end] = cell_groups
    #pnodes = periodic_nodes(chain)
    #onormals = outward_normals(chain)
    #GT.mesh(;
    #  node_coordinates = node_coords,
    #  face_nodes = face_to_nodes,
    #  face_reference_id = face_to_refid,
    #  reference_spaces = refid_to_refface,
    #  periodic_nodes = pnodes,
    #  group_faces = groups,
    #  outward_normals = onormals)
    GT.mesh(;
            node_coordinates = node_coordinates(chain),
            face_nodes = face_nodes(chain),
            face_reference_id = face_reference_id(chain),
            reference_spaces = reference_spaces(chain),
            periodic_nodes = periodic_nodes(chain),
            group_faces = group_faces(chain),
            outward_normals = outward_normals(chain))
end

function mesh_space(mesh::AbstractMesh,D)
    vD = Val(val_parameter(D))
    MeshSpace(mesh,vD)
end

struct MeshSpace{A,B} <: AbstractSpace{A}
    mesh::A
    num_dims::Val{B}
end

function num_dims(space::MeshSpace)
    val_parameter(space.num_dims)
end

function reference_spaces(space::MeshSpace)
    reference_spaces(space.mesh,space.num_dims)
end

function face_reference_id(space::MeshSpace)
    D = num_dims(space)
    face_reference_id(space.mesh,D)
end

function domain(space::MeshSpace)
    domain(space.mesh,space.num_dims)
end

num_free_dofs(space::MeshSpace) = num_nondes(space.mesh)
num_dirichlet_dofs(space::MeshSpace) = 0

function face_dofs(space::MeshSpace)
    D = num_dims(space)
    face_dofs(space.mesh,D)
end








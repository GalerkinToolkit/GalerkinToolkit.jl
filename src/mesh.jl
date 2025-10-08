

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

"""
"""
function group_faces_in_dim! end

function group_faces_in_dim!(m::AbstractMesh,d;group_name="__$d-FACES__")
    groups = group_faces(m,d)
    if haskey(groups,group_name)
        return group_name
    end
    if is_partitioned(m)
        faces = map(partition(GT.face_ids(m,d))) do ids
            Ti = int_type(options(m))
            collect(Ti,own_to_global(ids))
        end |> Partitioned
    else
        Ti = int_type(options(m))
        faces = collect(Ti,1:num_faces(m,d))
    end
    groups[group_name] = faces
    group_name
end

#function group_faces_in_dim!(m::AbstractPMesh,d;group_name="__$d-FACES__")
#    p_mesh = partition(m)
#    foreach(p_mesh) do mesh
#        group_faces_in_dim!(mesh,d;group_name)
#    end
#    group_name
#end

"""
"""
function group_interior_faces! end

function group_interior_faces!(mesh::AbstractMesh;group_name="__INTERIOR_FACES__")
    D = num_dims(mesh)
    d = D-1
    groups = group_faces(mesh,d)
    if haskey(groups,group_name)
        return group_name
    end
    topo = topology(mesh)
    face_to_cells = face_incidence(topo,d,D)
    faces = findall(cells->length(cells)==2,face_to_cells)
    groups[group_name] = faces
    group_name
end

#function group_interior_faces!(pmesh::AbstractPMesh;group_name="__INTERIOR_FACES__")
#    D = num_dims(pmesh)
#    d = D - 1
#    vals = map(partition(pmesh)) do mesh
#        topo = topology(mesh)
#        nfaces = num_faces(topo,d)
#        face_to_v = zeros(Int32,nfaces)
#        cell_to_faces = face_incidence(topo,D,d)
#        ncells = num_faces(topo,D)
#        cell_ids = face_local_indices(mesh,D)
#        cell_to_owner = local_to_owner(cell_ids)
#        part = part_id(cell_ids)
#        for cell in 1:ncells
#            owner = cell_to_owner[cell]
#            if part != owner
#                continue
#            end
#            faces = cell_to_faces[cell]
#            for face in faces
#                face_to_v[face] += 1
#            end
#        end
#        face_to_v
#    end
#    ids = face_partition(pmesh, d)
#    v = PVector(vals,ids)
#    assemble!(v) |> wait
#    map(partition(pmesh),vals) do mesh, face_to_v
#       group_faces(mesh,d)[group_name] = findall(i->i==2,face_to_v)
#    end
#    group_name
#end

"""
"""
function group_boundary_faces! end

function group_boundary_faces!(mesh::AbstractMesh;group_name="__BOUNDARY_FACES__")
    D = num_dims(mesh)
    d = D-1
    groups = group_faces(mesh,d)
    if haskey(groups,group_name)
        return group_name
    end
    topo = topology(mesh)
    face_to_cells = face_incidence(topo,d,D)
    faces = findall(cells->length(cells)==1,face_to_cells)
    groups[group_name] = faces
    group_name
end

function group_boundary_faces!(domain::AbstractDomain;group_name="__BOUNDARY_$(objectid(domain))__")
    mesh = GT.mesh(domain)
    D = num_dims(domain)
    d = D-1
    groups = group_faces(mesh,d)
    if haskey(groups,group_name)
        return group_name
    end
    topo = topology(mesh)
    cell_to_faces = face_incidence(topo,D,d)
    nfaces = num_faces(mesh,d)
    face_count = zeros(Int32,nfaces)
    for cell in GT.faces(domain)
        faces = cell_to_faces[cell]
        for face in faces
            face_count[face] += 1
        end
    end
    faces = findall(count->count==1,face_count)
    groups[group_name] = faces
    group_name
end

#function group_boundary_faces!(pmesh::AbstractPMesh;group_name="__BOUNDARY_FACES__")
#    D = num_dims(pmesh)
#    d = D - 1
#    vals = map(partition(pmesh)) do mesh
#        topo = topology(mesh)
#        nfaces = num_faces(topo,d)
#        face_to_v = zeros(Int32,nfaces)
#        cell_to_faces = face_incidence(topo,D,d)
#        ncells = num_faces(topo,D)
#        cell_ids = face_local_indices(mesh,D)
#        cell_to_owner = local_to_owner(cell_ids)
#        part = part_id(cell_ids)
#        for cell in 1:ncells
#            owner = cell_to_owner[cell]
#            if part != owner
#                continue
#            end
#            faces = cell_to_faces[cell]
#            for face in faces
#                face_to_v[face] += 1
#            end
#        end
#        face_to_v
#    end
#    ids = face_partition(pmesh, d)
#    v = PVector(vals,ids)
#    assemble!(v) |> wait
#    map(partition(pmesh),vals) do mesh, face_to_v
#       group_faces(mesh,d)[group_name] = findall(i->i==1,face_to_v)
#    end
#    group_name
#end

"""
"""
function domain(mesh::AbstractMesh,d;
    mesh_id = objectid(mesh),
    faces_around=nothing,
    is_reference_domain=Val(false),
    group_names=[group_faces_in_dim!(mesh,val_parameter(d))],
    )
    mesh_domain(mesh;
        mesh_id,
        num_dims=Val(val_parameter(d)),
        group_names,
        is_reference_domain)
end

function interior(mesh::AbstractMesh;
    mesh_id = objectid(mesh),
    group_names=[group_faces_in_dim!(mesh,num_dims(mesh))],
    is_reference_domain=Val(false)
    )
    d = num_dims(mesh)
    mesh_domain(mesh;
        mesh_id,
        group_names,
        num_dims=Val(val_parameter(d)),
        is_reference_domain)
end

function skeleton(mesh::AbstractMesh;
    mesh_id = objectid(mesh),
    group_names=[group_interior_faces!(mesh)],
    is_reference_domain=Val(false)
    )
    d = num_dims(mesh) - 1
    mesh_domain(mesh;
        mesh_id,
        group_names,
        num_dims=Val(val_parameter(d)),
        is_reference_domain)
end

function boundary(mesh::AbstractMesh;
    mesh_id = objectid(mesh),
    group_names=[group_boundary_faces!(mesh)],
    is_reference_domain=Val(false),
    faces_around = nothing,
    )
    d = num_dims(mesh) - 1
    domain = mesh_domain(mesh;
        mesh_id,
        group_names,
        faces_around,
        num_dims=Val(val_parameter(d)),
        is_reference_domain)
    if faces_around === nothing
        nfaces = num_faces(domain)
        new_faces_around = FillArrays.Fill(1,nfaces)
        domain2 = replace_faces_around(domain,new_faces_around)
    else
        domain2 = domain
    end
end

function boundary(domain::AbstractDomain;
    group_names=[group_boundary_faces!(domain)],
    kwargs...)
    mesh = GT.mesh(domain)
    nfaces = num_faces(domain)
    domain2 = boundary(mesh;group_names,kwargs...)
    faces_around = boundary_faces_around(domain,domain2)
    replace_faces_around(domain2,faces_around)
end

function boundary_faces_around(domain_D,domain_d)
    mesh = GT.mesh(domain_D)
    D = num_dims(domain_D)
    d = num_dims(domain_d)
    topo = topology(mesh)
    dface_Dfaces = GT.face_incidence(topo,d,D)
    Dface_mask = fill(false,num_faces(mesh,D))
    Dface_mask[GT.faces(domain_D)] .= true
    map(GT.faces(domain_d)) do dface
        Dfaces = dface_Dfaces[dface]
        i = -1
        for (j,Dface) in enumerate(Dfaces)
            mask = Dface_mask[Dface]
            if mask
                i = j
                break
            end
        end
        @boundscheck @assert i != -1
        Int32(i)
    end
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
         normals = normals(mesh),
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
    node_owner = periodic_nodes(mesh)
    if isa(node_owner,AbstractRange)
        lnode_lowner = 1:nlnodes
    else
        lnode_lowner = node_to_lnode[node_owner[lnode_to_node]]
    end
    if normals(mesh) !== nothing
        lnormals = normals(mesh)[lface_to_face_mesh[end]]
    else
        lnormals = nothing
    end

    lmesh = GT.mesh(;
        node_coordinates = lnode_to_coords,
        face_nodes = lface_to_lnodes_mesh,
        face_reference_id = lface_to_refid_mesh,
        reference_spaces = reference_spaces(mesh),
        group_faces = lgroups_mesh,
        periodic_nodes = lnode_lowner,
        normals = lnormals,
        kwargs...
        )

    lmesh
end

struct Mesh{A,B,C,D,E,F,G,H,I,J,K,W} <: AbstractMesh
    node_coordinates::A
    face_nodes::B
    face_reference_id::C
    reference_spaces::D
    periodic_nodes::E
    group_faces::F
    normals::G
    geometry_names::H
    is_cell_complex::I
    node_local_indices::J
    face_local_indices::K
    workspace::W
end

"""
    create_mesh(;kwargs...)

Build an arbitrary mesh object.

See also [`cartesian_mesh`](@ref), [`mesh_from_msh`](@ref), and [`mesh_from_gmsh`](@ref).

# Level

Intermediate

# Keyword arguments

- `node_coordinates`: The vector containing the coordinates of all mesh nodes. `node_coordinates[i]` is the coordinate vector for *global* node number `i`.
- `face_nodes`: A highly-nested vector containing the node ids for each face in the mesh. `node_coordinates[n]` with `n=face_nodes[d+1][i][k]` is the global node coordinate vector for *local* node number `k` in face `i` of dimension `d`. The object `face_nodes[d+1]` is a long vector of small vectors of integers. It is often represented using a `JaggedArray` object that uses continuous linear memory for performance.
- `reference_spaces`: A nested tuple containing the reference spaces for faces. `reference_spaces[d+1][i]` is the reference space number `i` in dimension `d`. Reference interpolation spaces are defined with functions like [`lagrange_space`](@ref).
- `face_reference_id` [optional]: A nested vector containing which reference space is assigned to each face. `reference_sapces[d+1][r]` with `r=face_reference_id[d+1][i]` is the reference space associated with face number `i` of dimension `d`. By default, all faces are assigned to the first reference space in its dimension.
- `group_faces` [optional]: A vector of dictionaries containing labeled groups of faces. `group_faces[d+1][group_name]` is a vector of integers containing the ids  of the faces of dimension `d` in the group named `group_name`. These groups might overlap. By default, no faces groups are created.
- `is_cell_complex=Val(false)` [optional]: `Val(true)` if the input data represents a cell complex, `Val(false)` otherwise.
- `normals=nothing` [optinal]: Vector containing the normal vectors for the faces of maximum dimension of the mesh. This is relevant for meshes of dimension `d` embedded in `d+1` dimensions as there is no way to tell which should be the orientation of the normals from the other quantities in the mesh.  `normals[f]` gives the normal vector of face number `f` of dimension `d=length(face_nodes)-1`.
"""
function create_mesh end


function create_mesh(;kwargs...)
    mesh(;kwargs...)
end

function default_face_reference_id(face_nodes)
    D = length(face_nodes)-1
    [ fill(Int8(1),length(face_nodes[i+1]))   for i in 0:D ] 
end

function mesh(;
        node_coordinates,
        face_nodes,
        face_reference_id = default_face_reference_id(face_nodes),
        reference_spaces,
        periodic_nodes = default_periodic_nodes(node_coordinates),
        group_faces = default_group_faces(reference_spaces),
        normals = nothing,
        geometry_names = [ String[] for d in 1:length(face_reference_id)],
        is_cell_complex = Val(false),
        node_local_indices = PartitionedArrays.block_with_constant_size(1,(1,),(length(node_coordinates),)),
        face_local_indices = [ PartitionedArrays.block_with_constant_size(1,(1,),(length(face_reference_id[d]),)) for d in 1:length(face_reference_id)],
        workspace = nothing,
    )
    mesh = Mesh(
                node_coordinates,
                face_nodes,
                face_reference_id,
                reference_spaces,
                periodic_nodes,
                group_faces,
                normals,
                geometry_names,
                is_cell_complex,
                node_local_indices,
                face_local_indices,
                workspace,
               )
end

function replace_workspace(mesh::Mesh,workspace)
    Mesh(
         mesh.node_coordinates,
         mesh.face_nodes,
         mesh.face_reference_id,
         mesh.reference_spaces,
         mesh.periodic_nodes,
         mesh.group_faces,
         mesh.normals,
         mesh.geometry_names,
         mesh.is_cell_complex,
         mesh.node_local_indices,
         mesh.face_local_indices,
         workspace,
        )
end

function mesh_workspace(;topology)
    MeshWorkspace(topology)
end

struct MeshWorkspace{A} <: AbstractType
    topology::A
end

function replace_node_coordinates(mesh::Mesh,node_coordinates)
    Mesh(
         node_coordinates,
         mesh.face_nodes,
         mesh.face_reference_id,
         mesh.reference_spaces,
         mesh.periodic_nodes,
         mesh.group_faces,
         mesh.normals,
         mesh.geometry_names,
         mesh.is_cell_complex,
         mesh.node_local_indices,
         mesh.face_local_indices,
         mesh.workspace,
        )
end

function replace_periodic_nodes(mesh::Mesh,periodic_nodes)
    @assert mesh.workspace === nothing
    Mesh(
         mesh.node_coordinates,
         mesh.face_nodes,
         mesh.face_reference_id,
         mesh.reference_spaces,
         periodic_nodes,
         mesh.group_faces,
         mesh.normals,
         mesh.geometry_names,
         mesh.is_cell_complex,
         mesh.node_local_indices,
         mesh.face_local_indices,
         mesh.workspace,
        )
end

node_coordinates(m::Mesh) = m.node_coordinates
face_nodes(m::Mesh) = m.face_nodes
face_nodes(m::Mesh,d) = m.face_nodes[val_parameter(d)+1]
face_reference_id(m::Mesh) = m.face_reference_id
face_reference_id(m::Mesh,d) = m.face_reference_id[val_parameter(d)+1]
reference_spaces(m::Mesh) = m.reference_spaces
reference_spaces(m::Mesh,d) = m.reference_spaces[val_parameter(d)+1]
group_faces(m::Mesh) = m.group_faces
group_faces(m::Mesh,d) = m.group_faces[val_parameter(d)+1]
periodic_nodes(m::Mesh) = m.periodic_nodes
is_cell_complex(m::Mesh) = val_parameter(m.is_cell_complex)
workspace(m::Mesh) = m.workspace

node_local_indices(m::Mesh) = m.node_local_indices
face_local_indices(m::Mesh) = m.face_local_indices
face_local_indices(m::Mesh,d) = m.face_local_indices[d+1]

function default_group_faces(reference_spaces)
    [ Dict{String,Vector{int_type(options(first(last(reference_spaces))))}}() for _ in 1:length(reference_spaces) ]
end

function default_periodic_nodes(node_coordinates)
    nnodes = length(node_coordinates)
    1:nnodes
end

function normals(m::Mesh)
    m.normals
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
    a.geometry_names
end

struct Chain{A} <: AbstractMesh
    contents::A
end

"""
    create_chain(;kwargs...)

Build an arbitrary mesh object, containing all faces of the same dimension. This function
is similar to [`create_mesh`](@ref) but it only receives face arrays of one dimension.

See also [`create_mesh`](@ref).

# Level

Intermediate

# Keyword arguments

- `node_coordinates`: Like for [`create_mesh`](@ref).
- `face_nodes`: A nested vector containing the node ids for each face in the mesh. `node_coordinates[n]` with `n=face_nodes[i][k]` is the global node coordinate vector for *local* node number `k` in face `i`. 
- `reference_spaces`: A tuple containing the reference spaces for faces. `reference_spaces[i]` is the reference space number `i`.
- `face_reference_id` [optional]: A vector containing which reference space is assigned to each face. `reference_sapces[r]` with `r=face_reference_id[i]` is the reference space associated with face number `i`. By default, all faces are assigned to the first reference space in its dimension.
- `group_faces` [optional]: A Dictionary containing labeled groups of faces. `group_faces[group_name]` is a vector of integers containing the ids  of the faces in the group named `group_name`. These groups might overlap. By default, no faces groups are created.
- `normals=nothing` [optinal]: Like for [`create_mesh`](@ref).
"""
function create_chain(;kwargs...)
    mesh(chain(;kwargs...))
end

function chain(;
        node_coordinates,
        face_nodes,
        face_reference_id = default_face_reference_id([face_nodes])[1],
        reference_spaces,
        periodic_nodes = default_periodic_nodes(node_coordinates),
        group_faces = default_group_faces((reference_spaces,))[end],
        normals = nothing,
    )
    contents = (;
                node_coordinates,
                face_nodes,
                face_reference_id,
                reference_spaces,
                periodic_nodes,
                group_faces,
                normals,
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

normals(m::Chain) = m.contents.normals

function chain(mesh::AbstractMesh,D=Val(num_dims(mesh)))
    d = val_parameter(D)
    chain(;
          node_coordinates=node_coordinates(mesh),
          face_nodes=face_nodes(mesh,d),
          face_reference_id=face_reference_id(mesh,d),
          reference_spaces=reference_spaces(mesh,d),
          periodic_nodes=periodic_nodes(mesh),
          group_faces=group_faces(mesh,d),
          normals=normals(mesh),
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
    #onormals = normals(chain)
    #GT.mesh(;
    #  node_coordinates = node_coords,
    #  face_nodes = face_to_nodes,
    #  face_reference_id = face_to_refid,
    #  reference_spaces = refid_to_refface,
    #  periodic_nodes = pnodes,
    #  group_faces = groups,
    #  normals = onormals)
    GT.mesh(;
            node_coordinates = node_coordinates(chain),
            face_nodes = face_nodes(chain),
            face_reference_id = face_reference_id(chain),
            reference_spaces = reference_spaces(chain),
            periodic_nodes = periodic_nodes(chain),
            group_faces = group_faces(chain),
            normals = normals(chain))
end

function mesh_space(mesh::AbstractMesh,D)
    num_dims = Val(val_parameter(D))
    domain = GT.mesh_domain(mesh;num_dims)
    MeshSpace(num_dims,domain)
end

struct MeshSpace{A,B} <: AbstractSpace
    num_dims::Val{A}
    domain::B
end

function mesh(space::MeshSpace)
    mesh(space.domain)
end

function num_dims(space::MeshSpace)
    val_parameter(space.num_dims)
end

function reference_spaces(space::MeshSpace)
    reference_spaces(mesh(space),space.num_dims)
end

function face_reference_id(space::MeshSpace)
    D = num_dims(space)
    face_reference_id(mesh(space),D)
end

function domain(space::MeshSpace)
    space.domain
end

num_free_dofs(space::MeshSpace) = num_nondes(mesh(space))
num_dirichlet_dofs(space::MeshSpace) = 0

function face_dofs(space::MeshSpace)
    D = num_dims(space)
    face_nodes(mesh(space),D)
end

# PartitionedRelated

#TODO eventually rename to faces
# faces to mesh_faces
# and inverse_faces to domain_faces
function face_ids(mesh::AbstractMesh,D)
    axes(face_reference_id(mesh,D),1)
end

function nodes(mesh::AbstractMesh)
    axes(node_coordinates(mesh),1)
end

@inline function is_partitioned(mesh::AbstractMesh)
    isa(node_coordinates(mesh),PVector)
end

# Constructors

"""
    with_mesh_partitioner(mesher[,partitioner];[parts])

Generate a mesh calling `mesher()` partition it, and distribute it
over the part ids in `parts`.

### Arguments

* Function `mesher()` should have no arguments and returns a sequential mesh object. This function is called only on one process.
* `partitioner ` [optional]: A function that takes a graph encoded as a sparse matrix, and returns a vector containing the part id of each node in the graph. Defaults to `Metis.partition`.

### Keyword arguments

* `parts` [optional]: A vector containing the part indices `1:P` where `P` is the number of parts in the data distribution. By default, `P` is the number of MPI ranks and `1:P` is distributed one item per rank.
"""
function with_mesh_partitioner(mesher,
    partitioner = Metis.partition;
    parts=default_parts()
    )
    np = length(parts)
    mesh_main = map_main(parts) do part
        mesh = mesher()
        g = node_graph(mesh)
        node_part = partitioner(g,np)
        parts_seq = LinearIndices((np,))
        node_partition = PA.partition_from_color(parts_seq,node_part)
        partitioned_mesh(mesh,node_partition)
    end
    scatter_mesh(mesh_main)
end

function partitioned_mesh(mesh,node_partition)
    parts = PA.linear_indices(node_partition)
    node_coordinates = partitioned_vector(GT.node_coordinates(mesh),node_partition)
    node_part = PA.getany(map(PA.global_to_owner,node_partition))
    face_part = map(GT.face_nodes(mesh)) do dface_nodes
        map(dface_nodes) do nodes
            maximum(node->node_part[node],nodes)
        end
    end
    face_partition = map(face_part) do dface_part
        PA.partition_from_color(parts,dface_part)
    end
    face_nodes = partitioned_faces(GT.face_nodes(mesh),face_partition)
    face_reference_id = partitioned_faces(GT.face_reference_id(mesh),face_partition)
    group_faces = partitioned_groups(GT.group_faces(mesh),face_part,parts)
    reference_spaces = GT.reference_spaces(mesh)
    is_cell_complex = GT.is_cell_complex(mesh)
    if GT.workspace(mesh) !== nothing
        topology = partitioned_topology(GT.topology(mesh),face_partition)
        workspace = mesh_workspace(;topology)
    else
        workspace = nothing
    end
    create_mesh(;
        node_coordinates,
        face_nodes,
        face_reference_id,
        reference_spaces,
        group_faces,
        is_cell_complex,
        workspace
       )
end

function partitioned_vector(x,row_partition)
    pvector(row_partition) do ids
        convert(typeof(x),x[ids])
    end
end

function partitioned_faces(d_dface_x,d_dface_partition)
    map(d_dface_x,d_dface_partition) do dface_x,dface_partition
        partitioned_vector(dface_x,dface_partition)
    end
end

function partitioned_groups(group_faces,face_part,parts)
    map(group_faces,face_part) do group_dfaces, dface_part
        map(collect(keys(group_dfaces))) do group
            f_dface = group_dfaces[group]
            dfaces = map(parts) do part
                f_part = dface_part[f_dface]
                fs = findall(p->p==part,f_part)
                f_dface[fs]
            end
            group => Partitioned(dfaces)
        end |> Dict
    end
end

function partitioned_topology(topo,face_partition)
    d_dface_partition = face_partition
    D = num_dims(topo)
    ijs = CartesianIndices((D+1,D+1))
    face_incidence = map(ijs) do ij
        i,j = Tuple(ij)
        iface_jfaces = GT.face_incidence(topo)[i,j]
        iface_partition = d_dface_partition[i]
        partitioned_vector(iface_jfaces,iface_partition)
    end
    face_permutation_ids = map(ijs) do ij
        i,j = Tuple(ij)
        if i>=j
            iface_perms = GT.face_permutation_ids(topo)[i,j]
            iface_partition = d_dface_partition[i]
            partitioned_vector(iface_perms,iface_partition)
        else
            nothing
        end
    end
    face_reference_id = partitioned_faces(GT.face_reference_id(topo),face_partition)
    reference_topologies = GT.reference_topologies(topo)
    mesh_topology(;
                  face_incidence,
                  face_reference_id,
                  face_permutation_ids,
                  reference_topologies,)
end

function node_graph(mesh)
    d_to_cell_to_nodes = GT.face_nodes(mesh)
    nnodes = GT.num_nodes(mesh)
    node_graph_impl(nnodes,d_to_cell_to_nodes)
end

function node_graph(mesh,d)
    d_to_cell_to_nodes = [GT.face_nodes(mesh,d)]
    nnodes = GT.num_nodes(mesh)
    node_graph_impl(nnodes,d_to_cell_to_nodes)
end

function node_graph_impl(nnodes,d_to_cell_to_nodes)
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

function scatter_mesh(mesh_main)
    node_coordinates = scatter_node_coordinates(mesh_main)
    face_partition = scatter_face_partition(mesh_main)
    face_nodes = scatter_faces(GT.face_nodes,mesh_main,face_partition)
    face_reference_id = scatter_faces(GT.face_reference_id,mesh_main,face_partition)
    group_faces = scatter_group_faces(mesh_main)
    reference_spaces = PA.getany(multicast(map_main(GT.reference_spaces,mesh_main)))
    is_cell_complex = PA.getany(multicast(map_main(GT.is_cell_complex,mesh_main)))
    workspace = scatter_mesh_workspace(mesh_main,face_partition)
    create_mesh(;
        node_coordinates,
        face_nodes,
        face_reference_id,
        reference_spaces,
        group_faces,
        is_cell_complex,
        workspace
       )
end

function scatter_mesh_workspace(mesh_main,face_partition)
    compute_workspace = PA.getany(multicast(map_main(mesh->GT.workspace(mesh)!==nothing,mesh_main)))
    if compute_workspace
        topology = scatter_topology(mesh_main,face_partition)
        workspace = mesh_workspace(;topology)
    else
        workspace = nothing
    end
end

function scatter_node_coordinates(mesh_main)
    x = PA.getany(multicast(map_main(mesh->similar(GT.node_coordinates(mesh),0),mesh_main)))
    parts = linear_indices(mesh_main)
    node_part_main = map_main(mesh_main) do mesh
        PA.getany(map(PA.global_to_owner,partition(GT.nodes(mesh))))
    end
    node_part = PA.getany(multicast(node_part_main))
    node_partition = PA.partition_from_color(parts,node_part)
    otherwise = mesh -> [x]
    node_coordinates_partition = map_main(mesh_main;otherwise) do mesh
        partition(GT.node_coordinates(mesh))
    end |> scatter
    @assert eltype(node_coordinates_partition) != Any
    PVector(node_coordinates_partition,node_partition)
end

function scatter_face_partition(mesh_main)
    parts = linear_indices(mesh_main)
    D = PA.getany(multicast(map_main(num_dims,mesh_main)))
    dims = ntuple(d->d-1,Val(D+1))
    map(dims) do d
        dface_part_main = map_main(mesh_main) do mesh
            dfaces = GT.face_ids(mesh,d)
            PA.getany(map(PA.global_to_owner,partition(dfaces)))
        end
        dface_part = PA.getany(multicast(dface_part_main))
        PA.partition_from_color(parts,dface_part)
    end
end

function scatter_faces(f,mesh_main,face_partition)
    D = PA.getany(multicast(map_main(num_dims,mesh_main)))
    dims = ntuple(identity,Val(D+1))
    map(dims) do d
        a = PA.getany(multicast(map_main(mesh->PA.getany(partition(f(mesh)[d])),mesh_main)))
        otherwise = mesh -> [a]
        part_dface_nodes_main = map_main(mesh_main;otherwise) do mesh
            partition(f(mesh)[d])
        end
        part_dface_nodes = scatter(part_dface_nodes_main)
        part_dface_nodes_2 = map(v->convert(typeof(a),v),part_dface_nodes)
        PVector(part_dface_nodes_2,face_partition[d])
    end
end

function scatter_group_faces(mesh_main)
    D = PA.getany(multicast(map_main(num_dims,mesh_main)))
    dims = ntuple(i->i-1,Val(D+1))
    d_groups = map(dims) do d
        groups_main = map_main(mesh_main) do mesh
            group_dfaces = GT.group_faces(mesh,d)
            collect(keys(group_dfaces))
        end
        PA.getany(multicast(groups_main))
    end
    map(dims,d_groups) do d,groups
        map(groups) do group
            otherwise = mesh -> [Int32[]]
            dfaces_main = map_main(mesh_main;otherwise) do mesh
                group_dfaces = GT.group_faces(mesh,d)
                dfaces = group_dfaces[group]
                partition(dfaces)
            end
            part_dfaces = scatter(dfaces_main)
            group => Partitioned(part_dfaces)
        end |> Dict
    end
end

function scatter_topology(mesh_main,face_partition)
    topo_main = map_main(topology,mesh_main)
    face_reference_id = scatter_faces(GT.face_reference_id,topo_main,face_partition)
    reference_topologies = PA.getany(multicast(map_main(GT.reference_topologies,topo_main)))
    face_incidence = scatter_face_incidence(topo_main,face_partition)
    face_permutation_ids = scatter_face_permutation_ids(topo_main,face_partition)
    mesh_topology(;
                  face_incidence,
                  face_reference_id,
                  face_permutation_ids,
                  reference_topologies,)
end

function scatter_face_incidence(topo_main,face_partition)
    D = PA.getany(multicast(map_main(num_dims,topo_main)))
    ijs = CartesianIndices((D+1,D+1))
    map(ijs) do ij
        i,j = Tuple(ij)
        a = PA.getany(multicast(map_main(topo->PA.getany(partition(GT.face_incidence(topo)[i,j])),topo_main)))
        otherwise = mesh -> [a]
        part_iface_jfaces_main = map_main(topo_main;otherwise) do topo
            iface_jfaces = GT.face_incidence(topo)[i,j]
            partition(iface_jfaces)
        end
        part_iface_jfaces = scatter(part_iface_jfaces_main)
        part_iface_jfaces_2 = map(v->convert(typeof(a),v),part_iface_jfaces)
        PVector(part_iface_jfaces_2,face_partition[i])
    end
end

function scatter_face_permutation_ids(topo_main,face_partition)
    D = PA.getany(multicast(map_main(num_dims,topo_main)))
    ijs = CartesianIndices((D+1,D+1))
    map(ijs) do ij
        i,j = Tuple(ij)
        if i>=j
            a = PA.getany(multicast(map_main(topo->PA.getany(partition(GT.face_permutation_ids(topo)[i,j])),topo_main)))
            otherwise = mesh -> [a]
            part_iface_jfaces_main = map_main(topo_main;otherwise) do topo
                iface_jfaces = GT.face_permutation_ids(topo)[i,j]
                partition(iface_jfaces)
            end
            part_iface_jfaces = scatter(part_iface_jfaces_main)
            part_iface_jfaces_2 = map(v->convert(typeof(a),v),part_iface_jfaces)
            PVector(part_iface_jfaces_2,face_partition[i])
        else
            nothing
        end
    end
end

function PartitionedArrays.centralize(mesh::AbstractMesh)
    if ! is_partitioned(mesh)
        return mesh
    end
    node_coordinates = collect(GT.node_coordinates(mesh))
    face_nodes = map(collect,GT.face_nodes(mesh))
    face_reference_id = map(collect,GT.face_reference_id(mesh))
    group_faces = map(GT.group_faces(mesh)) do group_dfaces
        map(collect(keys(group_dfaces))) do group
            group => collect(group_dfaces[group])
        end |> Dict
    end
    reference_spaces = GT.reference_spaces(mesh)
    is_cell_complex = GT.is_cell_complex(mesh)
    workspace = GT.centralize(GT.workspace(mesh))
    create_mesh(;
        node_coordinates,
        face_nodes,
        face_reference_id,
        reference_spaces,
        group_faces,
        is_cell_complex,
        workspace
       )
end

function PartitionedArrays.centralize(w::MeshWorkspace)
    MeshWorkspace(centralize(w.topology))
end

function PartitionedArrays.centralize(topo::AbstractTopology)
    face_reference_id = map(collect,GT.face_reference_id(topo))
    reference_topologies = GT.reference_topologies(topo)
    face_incidence = map(collect,GT.face_incidence(topo))
    face_permutation_ids = map(GT.face_permutation_ids(topo)) do v
        if v !== nothing
            collect(v)
        else
            nothing
        end
    end
    mesh_topology(;
                  face_incidence,
                  face_reference_id,
                  face_permutation_ids,
                  reference_topologies,)
end

struct Partitioned{A} <: AbstractType
    partition::A
end

function PartitionedArrays.partition(a::Partitioned)
    a.partition
end

function Base.collect(a::Partitioned)
    reduce(vcat,collect(partition(a)))
end

# TODO move to PartitionedArrays

function foreach_part(f,args...)
    foreach(f,map(partition,args)...)
end

function map_parts(f,args...)
    map(f,map(partition,args)...)
end

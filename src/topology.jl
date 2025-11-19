
topology(a::AbstractTopology) = a
num_dims(t::AbstractTopology) = length(reference_topologies(t))-1
num_faces(t::AbstractTopology,d) = length(face_reference_id(t,d))

@auto_hash_equals struct FaceTopology{A,B} <: AbstractFaceTopology
    boundary::A
    vertex_permutations::B
end

function face_topology(;boundary,vertex_permutations)
    face_topology(boundary,vertex_permutations)
end

function face_topology(boundary,vertex_permutations)
    FaceTopology(boundary,vertex_permutations)
end

boundary(a::FaceTopology) = a.boundary
vertex_permutations(a::FaceTopology) = a.vertex_permutations

function num_dims(a::AbstractFaceTopology)
    if boundary(a) === nothing
        0
    else
        Int(num_dims(boundary(a))) + 1
    end
end

function num_faces(a::AbstractFaceTopology,d)
    m = num_dims(a)
    if d == m
        1
    else
        Int(num_faces(boundary(a),d))
    end
end

function face_incidence(a::AbstractFaceTopology)
    m = num_dims(a)
    if m==0
        return fill([[1]],1,1)
    end
    fi_boundary = face_incidence(boundary(a))
    T = typeof(fi_boundary[1,1])
    fi = Matrix{T}(undef,m+1,m+1)
    for j in 0:(m-1)
        for i in 0:(m-1)
            fi[i+1,j+1] = fi_boundary[i+1,j+1]
        end
        fi[j+1,end] = fill([1],num_faces(a,j))
        fi[end,j+1] = map(k->[k],1:num_faces(a,j))
    end
    fi
end

function face_incidence(a::AbstractFaceTopology,d,D)
    face_incidence(a)[d+1,D+1]
end

function face_reference_id(a::AbstractFaceTopology,d)
    m = num_dims(a)
    if m == d
        Int[1]
    else
        collect(Int,face_reference_id(boundary(a),d))
    end
end

function reference_topologies(a::AbstractFaceTopology,d)
    m = num_dims(a)
    if m == d
        (a,)
    else
        reference_topologies(boundary(a),d)
    end
end

function face_permutation_ids(a::AbstractFaceTopology)
    m = num_dims(a)
    if m==0
        return fill([[[1]]],1,1)
    end
    fi_boundary = face_permutation_ids(boundary(a))
    T = typeof(fi[1,1])
    fi = Matrix{T}(undef,m+1,m+1)
    for j in 0:(m-1)
        for i in 0:(m-1)
            fi[i+1,j+1] = fi_boundary[i+1,j+1]
        end
        fi[end,j] = fill([1],num_faces(a,j))
    end
    fi
end

function face_permutation_ids(a::AbstractFaceTopology,d,D)
    face_permutation_ids(a)[d+1,D+1]
end

@auto_hash_equals struct MeshTopologyWorkspace{A} <: AbstractType
    complexify_glue::A
end

complexify_glue(a::MeshTopologyWorkspace) = a.complexify_glue

#function replace_complexify_glue(a::MeshTopologyWorkspace,complexify_glue)
#    MeshTopologyWorkspace(complexify_glue)
#end

@auto_hash_equals struct MeshTopology{A,B,C,D,E,F,G} <: AbstractTopology
    face_incidence::A
    face_reference_id::B
    face_permutation_ids::C
    reference_topologies::D
    periodic_faces::E
    periodic_faces_permutation_id::F
    workspace::G
end

workspace(topo::MeshTopology) = topo.workspace
complexify_glue(topo::MeshTopology) = complexify_glue(workspace(topo))

function replace_periodic_nodes(workspace::MeshTopologyWorkspace,periodic_nodes)
    complexify_glue = replace_periodic_nodes(workspace.complexify_glue,periodic_nodes)
    MeshTopologyWorkspace(complexify_glue)
end

function replace_periodic_nodes(topo::MeshTopology,periodic_nodes)
    D = num_dims(topo)
    Ti = eltype(eltype(eltype(topo.face_incidence)))
    periodic_faces = Vector{Vector{Ti}}(undef,D+1)
    periodic_faces_permutation_id = Vector{Vector{Ti}}(undef,D+1)
    workspace = replace_periodic_nodes(GT.workspace(topo),periodic_nodes)
    MeshTopology(
                 topo.face_incidence,
                 topo.face_reference_id,
                 topo.face_permutation_ids,
                 topo.reference_topologies,
                 periodic_faces,
                 periodic_faces_permutation_id,
                 workspace,
                )
end


function mesh_topology(;
        face_incidence,
        face_reference_id,
        face_permutation_ids,
        reference_topologies,
        periodic_faces = map(a->1:length(a),face_reference_id),
        periodic_faces_permutation_id = map(a->ones(Int32,length(a)),face_reference_id),
        workspace = nothing
    )

    MeshTopology(
                 face_incidence,
                 face_reference_id,
                 face_permutation_ids,
                 reference_topologies,
                 periodic_faces,
                 periodic_faces_permutation_id,
                 workspace,
                )
end

#function replace_workspace(topo::MeshTopology,workspace)
#    MeshTopology(
#                 topo.face_incidence,
#                 topo.face_reference_id,
#                 topo.face_permutation_ids,
#                 topo.reference_topologies,
#                 topo.periodic_faces,
#                 topo.periodic_faces_permutation_id,
#                 workspace)
#end

#function replace_complexify_glue(topo::MeshTopology,complexify_glue)
#    ws = replace_complexify_glue(workspace(topo),complexify_glue)
#    replace_workspace(topo,ws)
#end

function num_faces(topo::MeshTopology)
    D = num_dims(topo)
    map(d->num_faces(topo,d),ntuple(i->i-1,Val(D+1)))
end

function num_faces(topo::MeshTopology,vd)
    d = val_parameter(vd)
    if GT.workspace(topo) === nothing
        return length(GT.face_reference_id(topo,d))
    end
    glue = GT.complexify_glue(GT.workspace(topo))
    D = num_dims(glue)
    if glue.num_faces[d+1] < 0
        if d == 0
            create_vertices!(glue)
        elseif d == D
            glue.num_faces[D+1] = num_faces(glue.parent_mesh,D)
        else
            create_face_boundary!(topo,d)
        end
    end
    glue.num_faces[d+1]
end

function face_incidence(a::MeshTopology)
    D = num_dims(a)
    # This is just to make sure we return
    # an array fully filled
    for m in 0:D
        for n in 0:D
            face_incidence(a,m,n)
        end
    end
    a.face_incidence
end

function face_incidence(topo::MeshTopology,vd,vD)
    d = val_parameter(vd)
    D = val_parameter(vD)
    if ! isassigned(topo.face_incidence,d+1,D+1)
        if D == d
            fill_face_interior_mesh_topology!(topo,d)
        elseif d == num_dims(topo) && D == 0
            glue = GT.complexify_glue(GT.workspace(topo))
            topo.face_incidence[d+1,D+1] = GT.parent_face_vertices(glue,d)
        elseif d == 0 && D == num_dims(topo)
            glue = GT.complexify_glue(GT.workspace(topo))
            topo.face_incidence[d+1,D+1] = GT.vertex_parent_faces(glue,D)
        elseif D == 0
            create_face_vertices!(topo,d)
        elseif d == D+1
            create_face_boundary!(topo,D)
        elseif D < d
            fill_face_boundary_mesh_topology!(topo,d,D)
        else
            GT.face_incidence(topo,D,d)
            #fill_face_boundary_mesh_topology!(topo,D,d)
            fill_face_coboundary_mesh_topology!(topo,D,d)
        end
    end
    topo.face_incidence[d+1,D+1]
end

function face_reference_id(a::MeshTopology)
    D = num_dims(a)
    # This is just to make sure we return
    # an array fully filled
    for d in 0:D
        face_reference_id(a,d)
    end
    a.face_reference_id
end

function face_reference_id(topo::MeshTopology,vd)
    d = val_parameter(vd)
    if ! isassigned(topo.face_reference_id, d+1)
        create_face_reference_id!(topo,d)
    end
    topo.face_reference_id[d+1]
end


function face_permutation_ids(a::MeshTopology)
    D = num_dims(a)
    # This is just to make sure we return
    # an array fully filled
    for d in 0:D
        for n in 0:d
            face_permutation_ids(a,d,n)
        end
    end
    a.face_permutation_ids
end
function face_permutation_ids(topo::MeshTopology,vd,vD)
    d = val_parameter(vd)
    D = val_parameter(vD)
    if ! isassigned(topo.face_permutation_ids,d+1,D+1)
        fill_face_permutation_ids!(topo,d,D)
    end
    topo.face_permutation_ids[d+1,D+1]
end
reference_topologies(a::MeshTopology) = a.reference_topologies
reference_topologies(a::MeshTopology,d) = a.reference_topologies[val_parameter(d)+1]

function periodic_faces(topo::MeshTopology,d)
    if eltype(topo.periodic_faces) <: AbstractRange
        if first(topo.periodic_faces[d+1]) < 0
            fill_periodic_faces!(topo,d)
        end
    else
        if ! isassigned(topo.periodic_faces,d+1) || first(topo.periodic_faces[d+1]) == -1
            fill_periodic_faces!(topo,d)
        end
    end
    topo.periodic_faces[d+1]
end

function periodic_faces_permutation_id(topo::MeshTopology,d)
    if eltype(topo.periodic_faces) <: AbstractRange
        if topo.periodic_faces_permutation_id[d+1].value < 0
            fill_periodic_faces_permutation_id!(topo,d)
        end
    else
        if ! isassigned(topo.periodic_faces_permutation_id,d+1) || first(topo.periodic_faces_permutation_id[d+1]) == -1
            fill_periodic_faces_permutation_id!(topo,d)
        end
    end
    topo.periodic_faces_permutation_id[d+1]
end

function face_local_faces(topo,d,D)
    dface_to_Dface_around_to_Dface = GT.face_incidence(topo,d,D) |> JaggedArray
    Dface_to_ldface_to_dface = GT.face_incidence(topo,D,d)
    data = copy(dface_to_Dface_around_to_Dface.data)
    fill!(data,0)
    ptrs = dface_to_Dface_around_to_Dface.ptrs
    dface_to_Dface_around_to_ldface = JaggedArray(data,ptrs)
    ndfaces = num_faces(topo,d)
    for dface in 1:ndfaces
        Dface_around_to_Dface = dface_to_Dface_around_to_Dface[dface]
        for (Dface_around, Dface) in enumerate(Dface_around_to_Dface)
            ldface_to_dface = Dface_to_ldface_to_dface[Dface]
            for (ldface2,dface2) in enumerate(ldface_to_dface)
                if dface == dface2
                    dface_to_Dface_around_to_ldface[dface][Dface_around] = ldface2
                    break
                end
            end
        end
    end
    dface_to_Dface_around_to_ldface
end

function face_incidence_ext(topo,d,D)
    dface_to_Dface_around_to_Dface = GT.face_incidence(topo,d,D) |> JaggedArray
    Dface_to_ldface_to_dface = GT.face_incidence(topo,D,d)
    data = copy(dface_to_Dface_around_to_Dface.data)
    fill!(data,0)
    ptrs = dface_to_Dface_around_to_Dface.ptrs
    dface_to_Dface_around_to_ldface = JaggedArray(data,ptrs)
    ndfaces = num_faces(topo,d)
    for dface in 1:ndfaces
        Dface_around_to_Dface = dface_to_Dface_around_to_Dface[dface]
        for (Dface_around, Dface) in enumerate(Dface_around_to_Dface)
            ldface_to_dface = Dface_to_ldface_to_dface[Dface]
            for (ldface2,dface2) in enumerate(ldface_to_dface)
                if dface == dface2
                    dface_to_Dface_around_to_ldface[dface][Dface_around] = ldface2
                    break
                end
            end
        end
    end
    dface_to_Dface_around_to_Dface,dface_to_Dface_around_to_ldface
end

#function topology(mesh::Mesh)
#    @assert workspace(mesh) !== nothing
#    mesh_workspace = GT.workspace(mesh)
#    complexify_glue = GT.complexify_glue(mesh_workspace)
#    topology = GT.topology(mesh_workspace)
#    topology_workspace = GT.workspace(topology)
#    replace_complexify_glue(topology,complexify_glue)
#end

"""
"""
function topology(mesh::AbstractMesh)
    if workspace(mesh) !== nothing
        return GT.topology(workspace(mesh))
    end
    # Assumes that the input is a face complex
    # We end up here when computing face topologies
    # The result is fully initialized in this case.
    Ti = Int32
    T = JaggedArray{Ti,Ti}
    D = num_dims(mesh)
    my_face_incidence = Matrix{T}(undef,D+1,D+1)
    my_face_reference_id  = [ face_reference_id(mesh,d) for d in 0:D ]
    dims = ntuple(d->d-1,Val(D+1))
    my_reference_faces = map(d->map(GT.topology,reference_domains(mesh,d)),dims)
    my_face_permutation_ids = Matrix{T}(undef,D+1,D+1)
    periodic_nodes = GT.periodic_nodes(mesh)
    if periodic_nodes isa AbstractRange
        periodic_faces = map(t->(-1:0),0:D)
        periodic_faces_permutation_id = map(t->Fill(-Ti(1),1),0:D)
    else
        periodic_faces = Vector{Vector{Ti}}(undef,D+1)
        periodic_faces_permutation_id = Vector{Vector{Ti}}(undef,D+1)
    end
    topo = mesh_topology(;
        face_incidence = my_face_incidence,
        face_reference_id = my_face_reference_id,
        face_permutation_ids = my_face_permutation_ids,
        reference_topologies = my_reference_faces,
        periodic_faces,
        periodic_faces_permutation_id,
       )
    for d in 0:D
        fill_face_interior_mesh_topology!(topo,d)
    end
    for d in 1:D
        fill_face_vertices_mesh_topology!(topo,mesh,d)
        fill_face_coboundary_mesh_topology!(topo,d,0)
    end
    for d in 1:(D-1)
        for n in (D-d):-1:1
            m = n+d
            fill_face_boundary_mesh_topology!(topo,m,n)
            fill_face_coboundary_mesh_topology!(topo,m,n)
        end
    end
    for d in 0:D
        for n in 0:d
            fill_face_permutation_ids!(topo,d,n)
        end
        # This loop is to make sure that we 
        # fully initialize all entries in face_permutation_ids
        for n in (d+1):D
            topo.face_permutation_ids[d+1,n+1] = JaggedArray([Int32[]])
        end
    end
    #fill_periodic_vertices!(topo,mesh)
    #fill_periodic_vertices_permutation_id!(topo)
    for d in 0:D
        fill_periodic_faces!(topo,d)
        fill_periodic_faces_permutation_id!(topo,d)
    end
    topo
end

function fill_face_interior_mesh_topology!(mesh_topology,d)
    n = num_faces(mesh_topology,d)
    ptrs = collect(Int32,1:(n+1))
    data = collect(Int32,1:n)
    mesh_topology.face_incidence[d+1,d+1] = JaggedArray(data,ptrs)
end

function generate_face_coboundary(nface_to_mfaces,nmfaces)
    ptrs = zeros(Int32,nmfaces+1)
    nnfaces = length(nface_to_mfaces)
    for nface in 1:nnfaces
        mfaces = nface_to_mfaces[nface]
        for mface in mfaces
            ptrs[mface+1] += Int32(1)
        end
    end
    length_to_ptrs!(ptrs)
    ndata = ptrs[end]-1
    data = zeros(Int32,ndata)
    for nface in 1:nnfaces
        mfaces = nface_to_mfaces[nface]
        for mface in mfaces
            p = ptrs[mface]
            data[p] = nface
            ptrs[mface] += Int32(1)
        end
    end
    rewind_ptrs!(ptrs)
    mface_to_nfaces = JaggedArray(data,ptrs)
    mface_to_nfaces
end

function fill_face_coboundary_mesh_topology!(topology,n,m)
    nmfaces = num_faces(topology,m)
    nface_to_mfaces = face_incidence(topology,n,m)
    topology.face_incidence[m+1,n+1] = generate_face_coboundary(nface_to_mfaces,nmfaces)
end

function fill_face_vertices_mesh_topology!(topo,mesh,d)
    function barrier(nnodes,vertex_to_nodes,dface_to_nodes,dface_to_refid,refid_to_lvertex_to_lnodes)
        vertex_to_node = JaggedArray(vertex_to_nodes).data
        #if vertex_to_node == 1:nnodes
        #    return dface_to_nodes
        #end
        node_to_vertex = zeros(Int32,nnodes)
        nvertices = length(vertex_to_nodes)
        node_to_vertex[vertex_to_node] = 1:nvertices
        ndfaces = length(dface_to_nodes)
        dface_to_vertices_ptrs = zeros(Int32,ndfaces+1)
        for dface in 1:ndfaces
            refid = dface_to_refid[dface]
            nlvertices = length(refid_to_lvertex_to_lnodes[refid])
            dface_to_vertices_ptrs[dface+1] = nlvertices
        end
        length_to_ptrs!(dface_to_vertices_ptrs)
        ndata = dface_to_vertices_ptrs[end]-1
        dface_to_vertices_data = zeros(Int32,ndata)
        for dface in 1:ndfaces
            refid = dface_to_refid[dface]
            lvertex_to_lnodes = refid_to_lvertex_to_lnodes[refid]
            nlvertices = length(lvertex_to_lnodes)
            lnode_to_node = dface_to_nodes[dface]
            offset = dface_to_vertices_ptrs[dface]-1
            for lvertex in 1:nlvertices
                lnode = first(lvertex_to_lnodes[lvertex])
                node = lnode_to_node[lnode]
                vertex = node_to_vertex[node]
                dface_to_vertices_data[offset+lvertex] = vertex
            end
        end
        dface_to_vertices = JaggedArray(dface_to_vertices_data,dface_to_vertices_ptrs)
    end
    nnodes = num_nodes(mesh)
    vertex_to_nodes = face_nodes(mesh,0)
    dface_to_nodes = face_nodes(mesh,d)
    dface_to_refid = face_reference_id(mesh,d)
    refid_refface = reference_spaces(mesh,d)
    refid_to_lvertex_to_lnodes = map(refface->face_nodes(remove_interior(complexify(refface)),0),refid_refface)
    topo.face_incidence[d+1,0+1] = barrier(nnodes,vertex_to_nodes,dface_to_nodes,dface_to_refid,refid_to_lvertex_to_lnodes)
end

function fill_face_boundary_mesh_topology!(topo,D,d)
    function barrier(
            Dface_to_vertices,
            vertex_to_Dfaces,
            dface_to_vertices,
            vertex_to_dfaces,
            Dface_to_refid,
            Drefid_to_ldface_to_lvertices)

        # Count
        ndfaces = length(dface_to_vertices)
        nDfaces = length(Dface_to_vertices)
        # Allocate output
        ptrs = zeros(Int32,nDfaces+1)
        for Dface in 1:nDfaces
            Drefid = Dface_to_refid[Dface]
            ldface_to_lvertices = Drefid_to_ldface_to_lvertices[Drefid]
            ptrs[Dface+1] = length(ldface_to_lvertices)
        end
        length_to_ptrs!(ptrs)
        ndata = ptrs[end]-1
        data = fill(Int32(INVALID_ID),ndata)
        Dface_to_dfaces = JaggedArray(data,ptrs)
        for Dface in 1:nDfaces
            Drefid = Dface_to_refid[Dface]
            ldface_to_lvertices = Drefid_to_ldface_to_lvertices[Drefid]
            lvertex_to_vertex = Dface_to_vertices[Dface]
            ldface_to_dface = Dface_to_dfaces[Dface]
            for (ldface,lvertices) in enumerate(ldface_to_lvertices)
                # Find the global d-face for this local d-face
                dface2 = Int32(INVALID_ID)
                #vertices = view(lvertex_to_vertex,lvertices)
                for (i,lvertex) in enumerate(lvertices)
                    vertex = lvertex_to_vertex[lvertex]
                    dfaces = vertex_to_dfaces[vertex]
                    for dface1 in dfaces
                        vertices1 = dface_to_vertices[dface1]
                        if same_valid_ids(lvertex_to_vertex,vertices1,lvertices,1:length(vertices1))
                            dface2 = dface1
                            break
                        end
                    end
                    if dface2 != Int32(INVALID_ID)
                        break
                    end
                end
                @boundscheck begin
                    msg = """


                    Error in: topology_from_mesh

                    The given mesh is provably not a face complex.
                    """
                    @assert dface2 != Int32(INVALID_ID) msg
                end
                ldface_to_dface[ldface] = dface2
            end # (ldface,lvertices)
        end # Dface
        Dface_to_dfaces
    end
    Dface_to_vertices = face_incidence(topo,D,0)
    vertex_to_Dfaces = face_incidence(topo,0,D)
    dface_to_vertices = face_incidence(topo,d,0)
    vertex_to_dfaces = face_incidence(topo,0,d)
    Dface_to_refid = face_reference_id(topo,D)
    refid_refface = reference_topologies(topo,D)
    Drefid_to_ldface_to_lvertices = map(refface->face_incidence(boundary(refface),d,0),refid_refface)
    Dface_to_dfaces = barrier(
            Dface_to_vertices,
            vertex_to_Dfaces,
            dface_to_vertices,
            vertex_to_dfaces,
            Dface_to_refid,
            Drefid_to_ldface_to_lvertices)
    topo.face_incidence[D+1,d+1] = Dface_to_dfaces
end

function fill_face_permutation_ids!(top,D,d)
    function barrier!(
            cell_to_lface_to_pindex,
            cell_to_lface_to_face,
            cell_to_cvertex_to_vertex,
            cell_to_ctype,
            ctype_to_lface_to_cvertices,
            face_to_fvertex_to_vertex,
            face_to_ftype,
            ftype_to_pindex_to_cfvertex_to_fvertex)

        ncells = length(cell_to_lface_to_face)
        for cell in 1:ncells
            ctype = cell_to_ctype[cell]
            lface_to_cvertices = ctype_to_lface_to_cvertices[ctype]
            a = cell_to_lface_to_face.ptrs[cell]-1
            c = cell_to_cvertex_to_vertex.ptrs[cell]-1
            for (lface,cfvertex_to_cvertex) in enumerate(lface_to_cvertices)
                face = cell_to_lface_to_face.data[a+lface]
                ftype = face_to_ftype[face]
                b = face_to_fvertex_to_vertex.ptrs[face]-1
                pindex_to_cfvertex_to_fvertex = ftype_to_pindex_to_cfvertex_to_fvertex[ftype]
                pindexfound = false
                for (pindex, cfvertex_to_fvertex) in enumerate(pindex_to_cfvertex_to_fvertex)
                    found = true
                    for (cfvertex,fvertex) in enumerate(cfvertex_to_fvertex)
                        vertex1 = face_to_fvertex_to_vertex.data[b+fvertex]
                        cvertex = cfvertex_to_cvertex[cfvertex]
                        vertex2 = cell_to_cvertex_to_vertex.data[c+cvertex]
                        if vertex1 != vertex2
                            found = false
                            break
                        end
                    end
                    if found
                        cell_to_lface_to_pindex.data[a+lface] = pindex
                        pindexfound = true
                        break
                    end
                end
                @assert pindexfound "Valid pindex not found"
            end
        end
    end
    @assert D >= d
    cell_to_lface_to_face = JaggedArray(face_incidence(top,D,d))
    data = similar(cell_to_lface_to_face.data,Int8)
    ptrs = cell_to_lface_to_face.ptrs
    cell_to_lface_to_pindex = JaggedArray(data,ptrs)
    if d == D || d == 0
        fill!(cell_to_lface_to_pindex.data,Int8(1))
        top.face_permutation_ids[D+1,d+1] = cell_to_lface_to_pindex
        return top
    end
    face_to_fvertex_to_vertex = JaggedArray(face_incidence(top,d,0))
    face_to_ftype = face_reference_id(top,d)
    ref_dfaces = reference_topologies(top,d)
    ftype_to_pindex_to_cfvertex_to_fvertex = map(vertex_permutations,ref_dfaces)
    cell_to_cvertex_to_vertex = JaggedArray(face_incidence(top,D,0))
    cell_to_ctype = face_reference_id(top,D)
    ref_Dfaces = reference_topologies(top,D)
    ctype_to_lface_to_cvertices = map(p->face_incidence(boundary(p),d,0),ref_Dfaces)
    barrier!(
             cell_to_lface_to_pindex,
             cell_to_lface_to_face,
             cell_to_cvertex_to_vertex,
             cell_to_ctype,
             ctype_to_lface_to_cvertices,
             face_to_fvertex_to_vertex,
             face_to_ftype,
             ftype_to_pindex_to_cfvertex_to_fvertex)
    top.face_permutation_ids[D+1,d+1] = cell_to_lface_to_pindex
    top
end

function intersection!(a,b,na,nb)
  function findeq!(i,a,b,nb)
    for j in 1:nb
      if a[i] == b[j]
        return
      end
    end
    a[i] = INVALID_ID
    return
  end
  for i in 1:na
    if a[i] == INVALID_ID
      continue
    end
    findeq!(i,a,b,nb)
  end
end

function find_eq(v,b,idsb)
    for i in idsb
        vs = b[i]
        if v == vs
            return true
        end
    end
    return false
end

function is_subset(a,b,idsa,idsb)
    for i in 1:length(idsa)
        v = a[idsa[i]]
        if v == INVALID_ID
            continue
        end
        if !find_eq(v,b,idsb)
            return false
        end
    end
    return true
end

function same_valid_ids(a,b,idsa,idsb)
    if !is_subset(a,b,idsa,idsb)
        return false
    end
    if !is_subset(b,a,idsb,idsa)
        return false
    end
    return true
end

#function fill_periodic_vertices!(topo,mesh)
#    d = 0
#    node_owner = periodic_nodes(mesh)
#    nnodes = num_nodes(mesh)
#    nvertices = num_faces(mesh,d)
#    vertex_nodes = JaggedArray(face_nodes(mesh,d))
#    vertex_node = vertex_nodes.data
#    node_vertex = zeros(Int32,nnodes)
#    node_vertex[vertex_node] = 1:nvertices
#    vertex_owner = zeros(Int32,nvertices)
#    vertex_owner .= 1:nvertices
#    for node in 1:nnodes
#        vertex = node_vertex[node]
#        if vertex == 0
#            continue
#        end
#        node2 = node_owner[node]
#        vertex2 = node_vertex[node2]
#        @boundscheck @assert vertex2 != 0
#        vertex_owner[vertex] = vertex2
#    end
#    topo.periodic_faces[d+1] = vertex_owner
#end
#
#function fill_periodic_vertices_permutation_id!(topo)
#    d = 0
#    nvertices = num_faces(topo,d)
#    topo.periodic_faces_permutation_id[d+1] = ones(Int32,nvertices)
#    topo
#end

function fill_periodic_faces!(topo,d)
    if eltype(topo.periodic_faces) <: AbstractRange
        topo.periodic_faces[d+1] = 1:num_faces(topo,d)
        return
    end
    if !isassigned(topo.periodic_faces,d+1)
        topo.periodic_faces[d+1] = zeros(Int32,num_faces(topo,d))
    end
    periodic_dfaces = topo.periodic_faces[d+1]
    if d == 0
        glue = complexify_glue(topo)
        parent_mesh = glue.parent_mesh
        periodic_nodes = GT.periodic_nodes(parent_mesh)
        vertex_node = GT.vertex_node(glue)
        node_vertex = GT.node_vertex(glue)
        periodic_dfaces[:] = node_vertex[periodic_nodes[vertex_node]]
        return
    end
    nfaces = num_faces(topo,d)
    periodic_dfaces .=  1:nfaces #collect(Int32,1:nfaces)
    face_vertices = face_incidence(topo,d,0)
    vertex_faces = face_incidence(topo,0,d)
    vertex_owner = periodic_faces(topo,0)
    maxfaces = maximum(faces->length(faces),vertex_faces)
    faces1 = fill(Int32(INVALID_ID),maxfaces)
    faces2 = fill(Int32(INVALID_ID),maxfaces)
    nfaces1 = 0
    nfaces2 = 0
    for face in 1:nfaces
        vertices = face_vertices[face]
        if all(vertex->vertex_owner[vertex]!=vertex,vertices)
            fill!(faces1,Int32(INVALID_ID))
            fill!(faces2,Int32(INVALID_ID))
            for (i,v) in enumerate(vertices)
                vertex = vertex_owner[v]
                faces = vertex_faces[vertex]
                if i == 1
                    copyto!(faces1,faces)
                    nfaces1 = length(faces)
                else
                    copyto!(faces2,faces)
                    nfaces2 = length(faces)
                    intersection!(faces1,faces2,nfaces1,nfaces2)
                end
            end
            for face1 in faces1
                if face1 != INVALID_ID
                    periodic_dfaces[face] = face1
                    break
                end
            end
        end
    end
    nothing
end

function fill_periodic_faces_permutation_id!(topo,d)
    if eltype(topo.periodic_faces) <: AbstractRange
        topo.periodic_faces_permutation_id[d+1] = Fill(1,num_faces(topo,d))
        return
    end
    if !isassigned(topo.periodic_faces_permutation_id,d+1)
        topo.periodic_faces_permutation_id[d+1] = ones(Int32,num_faces(topo,d))
    end
    periodic_dfaces_permutation_id = topo.periodic_faces_permutation_id[d+1]
    if d == 0
        #topo.periodic_faces_permutation_id[d+1] = fill(1,num_faces(topo,d))
        fill!(periodic_dfaces_permutation_id,1)
        return
    end
    rid_pindex_lvert1_lvert2 = map(vertex_permutations,reference_topologies(topo,d))
    face_rid = face_reference_id(topo,d)
    face_owner = periodic_faces(topo,d)
    face_vertices = face_incidence(topo,d,0)
    nfaces = length(face_owner)
    periodic_dfaces_permutation_id = ones(Int32,nfaces)
    periodic_vertices = GT.periodic_faces(topo,0)
    for face in 1:nfaces
        owner = face_owner[face]
        if face == owner
            continue
        end
        rid = face_rid[face]
        pindex_lvert1_lvert2 = rid_pindex_lvert1_lvert2[rid]
        vertices1 = face_vertices[face]
        vertices2 = face_vertices[owner]
        pindexfound = false
        for pindex in 1:length(pindex_lvert1_lvert2)
            lvert1_lvert2 = pindex_lvert1_lvert2[pindex]
            found = true
            for (lvert1,lvert2) in enumerate(lvert1_lvert2)
                vert1 = periodic_vertices[vertices1[lvert1]]
                vert2 = vertices2[lvert2]
                if vert1 != vert2
                    found = false
                    break
                end
            end
            if found
                periodic_dfaces_permutation_id[face] = pindex
                pindexfound = true
                break
            end
        end
        @assert pindexfound "Valid pindex not found"
    end
    #topo.periodic_faces_permutation_id[d+1]  = face_perm_id
    nothing
end

# mesh contains a workspace with topopolgy (with no workspace), glue and parent_mesh
# calling topology(mesh) returns a topology with glue and parent_mesh in the workspace
# In this way, we do not need to compute anything in advance, everything can be lazy
# We also need to remove the function mesh_topology(mesh) that computes it from a mesh
# or keep it, but in this case, mesh and parent_mesh will be the same in topology
# workspace

function complexify(parent_mesh::AbstractMesh;glue=Val(false))
    Ti = int_type(options(parent_mesh))
    T = JaggedArray{Ti,Ti}
    D = num_dims(parent_mesh)

    reference_spaces = unique_reference_spaces(parent_mesh)

    # Allocate complexify_glue
    vertex_node_ref = Base.RefValue{Vector{Ti}}()
    node_vertex_ref = Base.RefValue{Vector{Ti}}()
    parent_face_vertices = Vector{T}(undef,D+1)
    vertex_parent_faces = Vector{T}(undef,D+1)
    num_faces = fill(-1,D+1)
    parent_face_face = Vector{Vector{Ti}}(undef,D+1)
    complexify_glue = ComplexifyGlue(
                                      num_faces,
                                      parent_mesh,
                                      vertex_node_ref,
                                      node_vertex_ref,
                                      parent_face_vertices,
                                      vertex_parent_faces,
                                      parent_face_face,
                                     )

    # Allocate mesh topology
    face_incidence = Matrix{T}(undef,D+1,D+1)
    reference_topologies = map(spaces->Tuple(unique(map(space->GT.topology(GT.domain(space)),spaces))),reference_spaces)
    face_permutation_ids = Matrix{T}(undef,D+1,D+1)
    face_reference_id = Vector{Vector{Ti}}(undef,D+1)
    periodic_nodes = GT.periodic_nodes(parent_mesh)
    if periodic_nodes isa AbstractRange
        periodic_faces = map(t->(-1:0),0:D)
        periodic_faces_permutation_id = map(t->Fill(-Ti(1),1),0:D)
    else
        periodic_faces = Vector{Vector{Ti}}(undef,D+1)
        periodic_faces_permutation_id = Vector{Vector{Ti}}(undef,D+1)
    end
    topology = GT.mesh_topology(;
                                face_incidence,
                                face_reference_id,
                                face_permutation_ids,
                                reference_topologies,
                                periodic_faces,
                                periodic_faces_permutation_id,
                                workspace = MeshTopologyWorkspace(complexify_glue)
                               )


    # Allocate mesh
    node_coordinates = GT.node_coordinates(parent_mesh)
    face_nodes = Vector{T}(undef,D+1)
    normals = GT.normals(parent_mesh)
    group_faces = Vector{Dict{String,Vector{Ti}}}(undef,D+1)
    mesh = GT.create_mesh(;
                          node_coordinates,
                          face_nodes,
                          face_reference_id,
                          reference_spaces,
                          normals,
                          periodic_nodes,
                          group_faces,
                          is_face_complex = Val(true),
                          workspace = MeshWorkspace(topology)
                         )

    if val_parameter(glue)
        mesh, GT.parent_face_face(topology)
    else
        mesh
    end
end

struct ComplexifyGlue{A,B,C,D,E,F,G} <: AbstractType
    num_faces::A
    parent_mesh::B
    vertex_node_ref::C
    node_vertex_ref::D
    parent_face_vertices::E
    vertex_parent_faces::F
    parent_face_face::G
end

function replace_periodic_nodes(glue::ComplexifyGlue,periodic_nodes)
    parent_mesh = replace_periodic_nodes(glue.parent_mesh,periodic_nodes)
    ComplexifyGlue(
                   glue.num_faces,
                   parent_mesh,
                   glue.vertex_node_ref,
                   glue.node_vertex_ref,
                   glue.parent_face_vertices,
                   glue.vertex_parent_faces,
                   glue.parent_face_face,
                  )
end

#function complexify_glue(workspace::MeshWorkspace)
#    workspace.complexify_glue
#end
#
#function complexify_glue(mesh::AbstractMesh)
#    complexify_glue(workspace(mesh))
#end

function node_vertex(a::ComplexifyGlue)
    (;node_vertex_ref) = a
    if ! isassigned(node_vertex_ref)
        create_vertices!(a)
    end
    node_vertex_ref[]
end

function vertex_node(a::ComplexifyGlue)
    (;vertex_node_ref) = a
    if ! isassigned(vertex_node_ref)
        create_vertices!(a)
    end
    vertex_node_ref[]
end

function parent_face_vertices(a::ComplexifyGlue,d)
    (;parent_face_vertices,parent_mesh) = a
    if ! isassigned(parent_face_vertices,d+1)
        node_vertex = GT.node_vertex(a)
        parent_face_vertices[d+1] = fill_face_vertices(parent_mesh,d,node_vertex)
    end
    parent_face_vertices[d+1]
end

function vertex_parent_faces(a::ComplexifyGlue,d)
    (;vertex_parent_faces,) = a
    if ! isassigned(vertex_parent_faces,d+1)
        n_vertices = length(GT.vertex_node(a))
        parent_dface_vertices = GT.parent_face_vertices(a,d)
        vertex_parent_faces[d+1] = generate_face_coboundary(parent_dface_vertices,n_vertices)
    end
    vertex_parent_faces[d+1]
end

function num_dims(a::ComplexifyGlue)
    num_dims(a.parent_mesh)
end

function parent_face_face(topo::MeshTopology)
    D = num_dims(topo)
    for d in 0:D
        parent_face_face(topo,d)
    end
    complexify_glue(topo).parent_face_face
end

function parent_face_face(topo::MeshTopology,vd)
    a = complexify_glue(workspace(topo))
    d = val_parameter(vd)
    D = num_dims(a)
    if ! isassigned(a.parent_face_face,d+1)
        if d == 0
            create_vertices!(a)
        elseif d == D
            nDfaces = num_faces(a.parent_mesh,D)
            a.parent_face_face[D+1] = collect(1:nDfaces)
        else
            create_face_boundary!(topo,d)
        end
    end
    a.parent_face_face[d+1]
end

function create_vertices!(a::ComplexifyGlue)
    # TODO this can be optimized for linear meshes
    function barrier!(node_vertex,face_to_nodes,face_to_refid,refid_to_lvertex_to_lnodes)
        valid_id = -one(eltype(node_vertex))
        for face in eachindex(face_to_nodes)
            nodes = face_to_nodes[face]
            refid = face_to_refid[face]
            lvertex_to_lnodes = refid_to_lvertex_to_lnodes[refid]
            for lnodes in lvertex_to_lnodes
                lnode = first(lnodes)
                node = nodes[lnode]
                node_vertex[node] = valid_id
            end
        end
    end
    (;node_vertex_ref,vertex_node_ref,parent_mesh,num_faces,parent_face_face) = a
    Ti = Int32
    nnodes = num_nodes(parent_mesh)
    node_vertex = zeros(Ti,nnodes)
    D = num_dims(parent_mesh)
    for d in 0:D
        face_to_nodes = face_nodes(parent_mesh,d)
        if length(face_to_nodes) == 0
            continue
        end
        face_to_refid = face_reference_id(parent_mesh,d)
        refid_to_lvertex_to_lnodes = map(reference_spaces(parent_mesh,d)) do a
            if num_dims(domain(a)) != 0
                face_nodes(complexify(a),0)
            else
                [interior_nodes(a)]
            end
        end
        barrier!(node_vertex,face_to_nodes,face_to_refid,refid_to_lvertex_to_lnodes)
    end
    #valid_id = -one(eltype(node_vertex))
    vertex = Ti(0)
    # Loop needed to get the already existing vertices
    for nodes in GT.face_nodes(parent_mesh,0)
        node = nodes[1]
        if node_vertex[node] == -1
            vertex += Ti(1)
            node_vertex[node] = vertex
        end
    end
    for node in eachindex(node_vertex)
        if node_vertex[node] == -1
            vertex += Ti(1)
            node_vertex[node] = vertex
        end
    end
    nvertices = vertex
    vertex_node = zeros(Ti,nvertices)
    for (node,vertex) in enumerate(node_vertex)
        if vertex != 0
            vertex_node[vertex] = node
        end
    end
    node_vertex_ref[] = node_vertex
    vertex_node_ref[] = vertex_node
    num_faces[0+1] = nvertices
    parent_face_face[0+1] = collect(1:GT.num_faces(parent_mesh,0))
    nothing
end

function create_face_boundary!(topo::MeshTopology,d)
    glue = complexify_glue(workspace(topo))
    (;num_faces,parent_face_face) = glue
    n = val_parameter(d)+1
    if val_parameter(d) == 0
        num_faces[0+1] = nvertices
        create_face_vertices!(topo,1)
    end
    nface_vertices = GT.face_incidence(topo,n,0)
    vertex_nfaces = GT.face_incidence(topo,0,n)
    parent_dface_vertices = GT.parent_face_vertices(glue,d)
    vertex_parent_dfaces = GT.vertex_parent_faces(glue,d)
    nface_nrefid = GT.face_reference_id(topo,n)
    nrefid_refntopo = GT.reference_topologies(topo,Val(n))
    nrefid_ldface_lvertices = map(a->face_incidence(a,d,0),nrefid_refntopo)
    nface_dfaces, n_dfaces, parent_dface_dface = generate_face_boundary(
                                                                        nface_vertices,
                                                                        vertex_nfaces,
                                                                        parent_dface_vertices,
                                                                        vertex_parent_dfaces,
                                                                        nface_nrefid,
                                                                        nrefid_ldface_lvertices)
    topo.face_incidence[n+1,d+1] = nface_dfaces
    num_faces[d+1] = n_dfaces
    parent_face_face[d+1] = parent_dface_dface
    nothing
end

function create_face_vertices!(topo,d)
    glue = complexify_glue(workspace(topo))
    n = val_parameter(d)+1
    if val_parameter(d) == 0
        fill_face_interior_mesh_topology!(topo,0)
        return
    end
    D = num_dims(topo)
    if val_parameter(d) == D 
        topo.face_incidence[D+1,0+1] = GT.parent_face_vertices(glue,D)
        return
    end
    parent_dface_dface = GT.parent_face_face(topo,d)
    parent_dface_vertices = GT.parent_face_vertices(glue,d)
    n_dfaces = GT.num_faces(topo,d)
    nface_vertices = face_incidence(topo,n,0)
    nface_dfaces = face_incidence(topo,n,d)
    nrefid_refntopo = GT.reference_topologies(topo,Val(n))
    nrefid_ldface_lvertices = map(a->face_incidence(a,d,0),nrefid_refntopo)
    nface_nrefid = GT.face_reference_id(topo,n)
    dface_vertices = generate_face_vertices(
                                            n_dfaces,
                                            parent_dface_dface,
                                            parent_dface_vertices,
                                            nface_vertices,
                                            nface_dfaces,
                                            nface_nrefid,
                                            nrefid_ldface_lvertices)
    topo.face_incidence[d+1,0+1] = dface_vertices
    nothing
end

function create_face_nodes!(mesh,d)
    @assert d <= num_dims(mesh)
    topo = GT.topology(GT.workspace(mesh))
    glue = complexify_glue(GT.workspace(topo))
    (;parent_mesh,) = glue
    n = val_parameter(d)+1
    if val_parameter(d) == 0
        vertex_node = GT.vertex_node(glue)
        nvertices = length(vertex_node)
        ptrs = collect(Int32,1:(nvertices+1))
        mesh.face_nodes[0+1] = JaggedArray(vertex_node,ptrs)
        return
    end
    D = num_dims(topo)
    if val_parameter(d) == D 
        mesh.face_nodes[D+1] = GT.face_nodes(parent_mesh,D)
        return
    end
    parent_dface_dface = GT.parent_face_face(topo,d)
    parent_dface_nodes = GT.face_nodes(parent_mesh,d)
    n_dfaces = GT.num_faces(topo,d)
    nface_vertices = face_nodes(mesh,n)
    nface_dfaces = face_incidence(topo,n,d)
    nrefid_refnspace = GT.reference_spaces(mesh,Val(n))
    nrefid_ldface_lnodes = map(a->face_nodes(complexify(a),d),nrefid_refnspace)
    nface_nrefid = GT.face_reference_id(mesh,n)
    dface_nodes = generate_face_vertices(
                                            n_dfaces,
                                            parent_dface_dface,
                                            parent_dface_nodes,
                                            nface_vertices,
                                            nface_dfaces,
                                            nface_nrefid,
                                            nrefid_ldface_lnodes)
    mesh.face_nodes[d+1] = dface_nodes
    nothing
end

function create_face_reference_id!(mesh::Mesh,d)
    D = num_dims(mesh)
    topo = GT.topology(mesh)
    glue = GT.complexify_glue(topo)
    (;parent_mesh) = glue
    if D == val_parameter(d)
        mesh.face_reference_id[D+1] = GT.face_reference_id(parent_mesh,D)
        return
    end
    drefid_refdspace = GT.reference_spaces(mesh,d)
    if length(drefid_refdspace) == 1
        ndfaces = num_faces(topo,d)
        mesh.face_reference_id[d+1] = fill(Int32(1),ndfaces)
    else
        n = val_parameter(d) + 1
        nrefid_refnspace = GT.reference_spaces(mesh,Val(n))
        drefid_refdspace = GT.reference_spaces(mesh,d)
        nface_nrefid = GT.face_reference_id(mesh,Val(n))
        nface_dfaces = GT.face_incidence(topo,n,val_parameter(d))
        nrefid_ldrefid_refdspace = map(refnspace->GT.reference_spaces(complexify(refnspace),d),nrefid_refnspace)
        nredid_ldrefid_drefid = map(ldrefid_refdspace->indexin(ldrefid_refdspace,[drefid_refdspace...]),nrefid_ldrefid_refdspace)
        nrefid_ldface_ldrefid = map(refnspace->GT.face_reference_id(complexify(refnspace),d),nrefid_refnspace)
        ndfaces = num_faces(topo,d)
        nnfaces = num_faces(topo,n)
        dface_drefid = zeros(Int32,ndfaces)
        for nface in 1:nnfaces
            ldface_dface = nface_dfaces[nface]
            nrefid = nface_nrefid[nface]
            ldrefid_drefid = nredid_ldrefid_drefid[nrefid]
            ldface_ldrefid = nrefid_ldface_ldrefid[nrefid]
            nldfaces = length(ldface_dface)
            for ldface in 1:nldfaces
                ldrefid = ldface_ldrefid[ldface]
                drefid = ldrefid_drefid[ldrefid]
                dface = ldface_dface[ldface]
                dface_drefid[dface] = drefid
            end
        end
        mesh.face_reference_id[d+1] = dface_drefid
    end
    nothing
end

function create_face_reference_id!(topo::MeshTopology,d)
    D = num_dims(topo)
    glue = GT.complexify_glue(topo)
    (;parent_mesh) = glue
    if D == val_parameter(d)
        Dface_rid = GT.face_reference_id(parent_mesh,D)
        Drefid_refDtopo = GT.reference_topologies(topo,Val(D))
        if length(Drefid_refDtopo) == 1
            nDfaces = length(Dface_rid)
            Dface_Drefid = fill(Int32(1),nDfaces)
        else
            rid_rtopo = map(space->GT.topology(GT.domain(space)),GT.reference_spaces(parent_mesh,Val(D)))
            rid_Drefid = indexin(rid_rtopo,[Drefid_refDtopo...])
            Dface_Drefid = rid_Drefid[Dface_rid]
        end
        topo.face_reference_id[D+1] = Dface_Drefid
        return
    end
    drefid_refdtopo = GT.reference_topologies(topo,d)
    if length(drefid_refdtopo) == 1
        ndfaces = num_faces(topo,d)
        topo.face_reference_id[d+1] = fill(Int32(1),ndfaces)
    else
        n = val_parameter(d) + 1
        nrefid_refntopo = GT.reference_topologies(topo,Val(n))
        drefid_refdtopo = GT.reference_topologies(topo,d)
        nface_nrefid = GT.face_reference_id(topo,Val(n))
        nface_dfaces = GT.face_incidence(topo,n,val_parameter(d))
        nrefid_ldrefid_refdtopo = map(refntopo->GT.reference_topologies(refntopo,d),nrefid_refntopo)
        nredid_ldrefid_drefid = map(ldrefid_refdtopo->indexin(ldrefid_refdtopo,[drefid_refdtopo...]),nrefid_ldrefid_refdtopo)
        nrefid_ldface_ldrefid = map(refntopo->GT.face_reference_id(refntopo,d),nrefid_refntopo)
        ndfaces = num_faces(topo,d)
        nnfaces = num_faces(topo,n)
        dface_drefid = zeros(Int32,ndfaces)
        for nface in 1:nnfaces
            ldface_dface = nface_dfaces[nface]
            nrefid = nface_nrefid[nface]
            ldrefid_drefid = nredid_ldrefid_drefid[nrefid]
            ldface_ldrefid = nrefid_ldface_ldrefid[nrefid]
            nldfaces = length(ldface_dface)
            for ldface in 1:nldfaces
                ldrefid = ldface_ldrefid[ldface]
                drefid = ldrefid_drefid[ldrefid]
                dface = ldface_dface[ldface]
                dface_drefid[dface] = drefid
            end
        end
        topo.face_reference_id[d+1] = dface_drefid
    end
    nothing
end

function unique_reference_spaces(mesh)
    D = num_dims(mesh)
    ntuple(d->unique_reference_spaces(mesh,d-1),Val(D+1))
end

function unique_reference_spaces(mesh,d)
    if d == num_dims(mesh)
        return GT.reference_spaces(mesh,d)
    end
    drefid_to_ref_dface = GT.reference_spaces(mesh,d)
    n = d + 1
    nrefid_to_ref_nface = GT.reference_spaces(mesh,n)
    nrefid_to_drefrefid_to_ref_dface = map(a->reference_spaces(complexify(a),d),nrefid_to_ref_nface)
    i_to_ref_dface = collect(Any,drefid_to_ref_dface)
    drefid_to_i = collect(1:length(drefid_to_ref_dface))
    #i = length(i_to_ref_dface)
    #Ti = Int32
    #nrefid_to_drefrefid_to_i = map(a->zeros(Ti,length(a)),nrefid_to_drefrefid_to_ref_dface)
    for (nrefid,drefrefid_to_ref_dface) in enumerate(nrefid_to_drefrefid_to_ref_dface)
        for (drefrefid, ref_dface) in enumerate(drefrefid_to_ref_dface)
            push!(i_to_ref_dface,ref_dface)
            #i += 1
            #nrefid_to_drefrefid_to_i[nrefid][drefrefid] = i
        end
    end
    u_to_ref_dface = unique(i_to_ref_dface)
    Tuple(u_to_ref_dface)
end

function create_group_faces!(mesh,d)
    topo = topology(mesh)
    glue = complexify_glue(topo)
    (;parent_mesh) = glue
    parent_groups = GT.group_faces(parent_mesh,d)
    parent_dface_dface = GT.parent_face_face(topo,d)
    group_faces = Dict{String,Vector{Int32}}()
    for (group_name,parent_group_faces) in parent_groups
        group_faces[group_name] = parent_dface_dface[parent_group_faces]
    end
    mesh.group_faces[d+1] = group_faces
    nothing
end


#function complexify(mesh::AbstractMesh;glue=Val(false))
#    Ti = int_type(options(mesh))
#    T = JaggedArray{Ti,Ti}
#    D = num_dims(mesh)
#    oldface_to_newvertices = Vector{T}(undef,D+1)
#    newvertex_to_oldfaces = Vector{T}(undef,D+1)
#    newface_incidence = Matrix{T}(undef,D+1,D+1)
#    nnewfaces = zeros(Int,D+1)
#    newface_refid = Vector{Vector{Ti}}(undef,D+1)
#    newreffaces = Vector{Any}(undef,D+1)
#    newface_nodes = Vector{T}(undef,D+1)
#    old_to_new = Vector{Vector{Ti}}(undef,D+1)
#    node_to_newvertex, n_new_vertices = find_node_to_vertex(mesh) # Optimizable for linear meshes
#    for d in 0:D
#        oldface_to_newvertices[d+1] = fill_face_vertices(mesh,d,node_to_newvertex) # Optimizable for linear meshes
#        newvertex_to_oldfaces[d+1] = generate_face_coboundary(oldface_to_newvertices[d+1],n_new_vertices) # Optimizable for linear meshes
#    end
#    newface_incidence[D+1,0+1] = oldface_to_newvertices[D+1]
#    newface_incidence[0+1,D+1] = newvertex_to_oldfaces[D+1]
#    nnewfaces[D+1] = length(oldface_to_newvertices[D+1])
#    newface_refid[D+1] = face_reference_id(mesh,D)
#    newreffaces[D+1] = reference_spaces(mesh,D)
#    newface_nodes[D+1] = face_nodes(mesh,D)
#    old_to_new[D+1] = collect(Ti,1:length(newface_nodes[D+1]))
#    # TODO optimize for d==0
#    # Going through 0 is a bug since vertex ids are already created.
#    for d in (D-1):-1:1
#        n = d+1
#        new_nface_to_new_vertices = newface_incidence[n+1,0+1]
#        new_vertex_to_new_nfaces = newface_incidence[0+1,n+1]
#        old_dface_to_new_vertices = oldface_to_newvertices[d+1]
#        new_vertex_to_old_dfaces = newvertex_to_oldfaces[d+1]
#        new_nface_to_nrefid = newface_refid[n+1]
#        old_dface_to_drefid = face_reference_id(mesh,d)
#        drefid_to_ref_dface = reference_spaces(mesh,d)
#        old_dface_to_nodes = face_nodes(mesh,d)
#        new_nface_to_nodes = newface_nodes[n+1]
#        nrefid_to_ldface_to_lvertices = map(a->face_incidence(topology(GT.mesh(domain(a))),d,0),newreffaces[n+1])
#        nrefid_to_ldface_to_lnodes = map(a->face_nodes(complexify(a),d),newreffaces[n+1])
#        nrefid_to_ldface_to_drefrefid = map(a->face_reference_id(complexify(a),d),newreffaces[n+1])
#        nrefid_to_drefrefid_to_ref_dface = map(a->reference_spaces(complexify(a),d),newreffaces[n+1])
#        new_nface_to_new_dfaces, n_new_dfaces, old_dface_to_new_dface = generate_face_boundary(
#            new_nface_to_new_vertices,
#            new_vertex_to_new_nfaces,
#            old_dface_to_new_vertices,
#            new_vertex_to_old_dfaces,
#            new_nface_to_nrefid,
#            nrefid_to_ldface_to_lvertices)
#        new_dface_to_new_vertices = generate_face_vertices(
#            n_new_dfaces,
#            old_dface_to_new_dface,
#            old_dface_to_new_vertices,
#            new_nface_to_new_vertices,
#            new_nface_to_new_dfaces,
#            new_nface_to_nrefid,
#            nrefid_to_ldface_to_lvertices)
#        new_vertex_to_new_dfaces = generate_face_coboundary(new_dface_to_new_vertices,n_new_vertices)
#        new_dface_to_nodes = generate_face_vertices(
#            n_new_dfaces,
#            old_dface_to_new_dface,
#            old_dface_to_nodes,
#            new_nface_to_nodes,
#            new_nface_to_new_dfaces,
#            new_nface_to_nrefid,
#            nrefid_to_ldface_to_lnodes)
#        new_dface_to_new_drefid, new_refid_to_ref_dface = generate_reference_spaces(
#            n_new_dfaces,
#            old_dface_to_new_dface,
#            old_dface_to_drefid,
#            drefid_to_ref_dface,
#            new_nface_to_new_dfaces,
#            new_nface_to_nrefid,
#            nrefid_to_ldface_to_drefrefid,
#            nrefid_to_drefrefid_to_ref_dface)
#        if d > 0
#            newface_incidence[n+1,d+1] = new_nface_to_new_dfaces
#            newface_incidence[d+1,0+1] = new_dface_to_new_vertices
#            newface_incidence[0+1,d+1] = new_vertex_to_new_dfaces
#        end
#        newface_refid[d+1] = new_dface_to_new_drefid
#        newreffaces[d+1] = new_refid_to_ref_dface
#        nnewfaces[d+1] = n_new_dfaces
#        newface_nodes[d+1] = new_dface_to_nodes
#        old_to_new[d+1] = old_dface_to_new_dface
#    end
#
#    newvertex_to_node = zeros(Int32,n_new_vertices)
#    [node_to_newvertex]
#
#
#    node_to_coords = node_coordinates(mesh)
#    old_group_faces = GT.group_faces(mesh)
#    new_group_faces = [ Dict{String,Vector{Int32}}() for d in 0:D] # TODO hardcoded
#    for d in 0:D
#        old_groups = old_group_faces[d+1]
#        for (group_name,old_group_faces) in old_groups
#            old_to_new_d = old_to_new[d+1]
#            new_group_faces[d+1][group_name] = ld_to_new_d[old_group_faces]
#        end
#    end
#    new_mesh = GT.mesh(;
#            node_coordinates = node_to_coords,
#            face_nodes = newface_nodes,
#            face_reference_id = newface_refid,
#            reference_spaces = Tuple(newreffaces),
#            group_faces = new_group_faces,
#            periodic_nodes = periodic_nodes(mesh),
#            normals = normals(mesh),
#            geometry_names = geometry_names(mesh),
#            is_face_complex = Val(true),
#           )
#    mtopology = mesh_topology(;
#                          face_incidence = newface_incidence, #m_face_incidence,
#                          face_reference_id = newface_refid,
#                          face_permutation_ids = Matrix{JaggedArray{Int32,Int32}}(undef,D+1,D+1),
#                          reference_topologies = map(reffaces->map(refface->GT.topology(GT.domain(refface)),reffaces),Tuple(newreffaces)),
#                          periodic_faces = Vector{Vector{Int32}}(undef,D+1),
#                          periodic_faces_permutation_id = Vector{Vector{Int32}}(undef,D+1)
#                         )
#    fill_periodic_vertices!(mtopology,new_mesh)
#    fill_periodic_vertices_permutation_id!(mtopology)
#    workspace = mesh_workspace(;topology=mtopology)
#    new_mesh2 = replace_workspace(new_mesh,workspace)
#    if val_parameter(glue)
#        new_mesh2, old_to_new
#    else
#        new_mesh2
#    end
#end

const INVALID_ID = 0

function generate_face_vertices(
    n_new_dfaces,
    old_dface_to_new_dface,
    old_dface_to_new_vertices,
    new_nface_to_new_vertices,
    new_nface_to_new_dfaces,
    new_nface_to_nrefid,
    nrefid_to_ldface_to_lvertices
    )

    Ti = eltype(eltype(old_dface_to_new_vertices))
    new_dface_to_touched = fill(false,n_new_dfaces)
    new_dface_to_new_vertices_ptrs = zeros(Ti,n_new_dfaces+1)
    n_old_dfaces = length(old_dface_to_new_dface)
    for old_dface in 1:n_old_dfaces
        new_dface = old_dface_to_new_dface[old_dface]
        new_vertices = old_dface_to_new_vertices[old_dface]
        new_dface_to_new_vertices_ptrs[new_dface+1] = length(new_vertices)
        new_dface_to_touched[new_dface] = true
    end
    n_new_nfaces = length(new_nface_to_new_vertices)
    for new_nface in 1:n_new_nfaces
        nrefid = new_nface_to_nrefid[new_nface]
        ldface_to_lvertices = nrefid_to_ldface_to_lvertices[nrefid]
        ldface_to_new_dface = new_nface_to_new_dfaces[new_nface]
        n_ldfaces = length(ldface_to_new_dface)
        for ldface in 1:n_ldfaces
            new_dface = ldface_to_new_dface[ldface]
            if new_dface_to_touched[new_dface]
                continue
            end
            lvertices = ldface_to_lvertices[ldface]
            new_dface_to_new_vertices_ptrs[new_dface+1] = length(lvertices)
            new_dface_to_touched[new_dface] = true
        end

    end
    length_to_ptrs!(new_dface_to_new_vertices_ptrs)
    ndata = new_dface_to_new_vertices_ptrs[end]-1
    new_dface_to_new_vertices_data = zeros(Ti,ndata)
    new_dface_to_new_vertices = JaggedArray(new_dface_to_new_vertices_data,new_dface_to_new_vertices_ptrs)
    fill!(new_dface_to_touched,false)
    for old_dface in 1:n_old_dfaces
        new_dface = old_dface_to_new_dface[old_dface]
        new_vertices_in = old_dface_to_new_vertices[old_dface]
        new_vertices_out = new_dface_to_new_vertices[new_dface]
        for i in 1:length(new_vertices_in)
            new_vertices_out[i] = new_vertices_in[i]
        end
        new_dface_to_touched[new_dface] = true
    end
    for new_nface in 1:n_new_nfaces
        nrefid = new_nface_to_nrefid[new_nface]
        ldface_to_lvertices = nrefid_to_ldface_to_lvertices[nrefid]
        ldface_to_new_dface = new_nface_to_new_dfaces[new_nface]
        n_ldfaces = length(ldface_to_new_dface)
        new_vertices_in = new_nface_to_new_vertices[new_nface]
        for ldface in 1:n_ldfaces
            new_dface = ldface_to_new_dface[ldface]
            if new_dface_to_touched[new_dface]
                continue
            end
            new_vertices_out = new_dface_to_new_vertices[new_dface]
            lvertices = ldface_to_lvertices[ldface]
            for i in 1:length(lvertices)
                new_vertices_out[i] = new_vertices_in[lvertices[i]]
            end
            new_dface_to_touched[new_dface] = true
        end

    end
    new_dface_to_new_vertices
end

#function generate_reference_spaces(
#        n_new_dfaces,
#        old_dface_to_new_dface,
#        old_dface_to_drefid,
#        drefid_to_ref_dface,
#        new_nface_to_new_dfaces,
#        new_nface_to_nrefid,
#        nrefid_to_ldface_to_drefrefid,
#        nrefid_to_drefrefid_to_ref_dface)
#
#    i_to_ref_dface = collect(Any,drefid_to_ref_dface)
#    drefid_to_i = collect(1:length(drefid_to_ref_dface))
#    i = length(i_to_ref_dface)
#    Ti = Int32
#    nrefid_to_drefrefid_to_i = map(a->zeros(Ti,length(a)),nrefid_to_drefrefid_to_ref_dface)
#    for (nrefid,drefrefid_to_ref_dface) in enumerate(nrefid_to_drefrefid_to_ref_dface)
#        for (drefrefid, ref_dface) in enumerate(drefrefid_to_ref_dface)
#            push!(i_to_ref_dface,ref_dface)
#            i += 1
#            nrefid_to_drefrefid_to_i[nrefid][drefrefid] = i
#        end
#    end
#    u_to_ref_dface = unique(i_to_ref_dface)
#    i_to_u = indexin(i_to_ref_dface,u_to_ref_dface)
#    new_dface_to_u = zeros(Ti,n_new_dfaces)
#    new_dface_to_touched = fill(false,n_new_dfaces)
#    for (old_dface,new_dface) in enumerate(old_dface_to_new_dface)
#        drefid = old_dface_to_drefid[old_dface]
#        i = drefid_to_i[drefid]
#        u = i_to_u[i]
#        new_dface_to_u[new_dface] = u
#        new_dface_to_touched[new_dface] = true
#    end
#    for (new_nface,nrefid) in enumerate(new_nface_to_nrefid)
#        ldface_to_drefrefid = nrefid_to_ldface_to_drefrefid[nrefid]
#        ldface_to_new_dface = new_nface_to_new_dfaces[new_nface]
#        drefrefid_to_i = nrefid_to_drefrefid_to_i[nrefid]
#        for (ldface,new_dface) in enumerate(ldface_to_new_dface)
#            new_dface = ldface_to_new_dface[ldface]
#            if new_dface_to_touched[new_dface]
#                continue
#            end
#            drefrefid = ldface_to_drefrefid[ldface]
#            i = drefrefid_to_i[drefrefid]
#            u = i_to_u[i]
#            new_dface_to_u[new_dface] = u
#            new_dface_to_touched[new_dface] = true
#        end
#    end
#    new_dface_to_u, Tuple(u_to_ref_dface)
#end

function generate_face_boundary(
    Dface_to_vertices,
    vertex_to_Dfaces,
    dface_to_vertices,
    vertex_to_dfaces,
    Dface_to_refid,
    Drefid_to_ldface_to_lvertices)

    # Count
    ndfaces = length(dface_to_vertices)
    nDfaces = length(Dface_to_vertices)
    nvertices = length(vertex_to_Dfaces)
    maxldfaces = 0
    for ldface_to_lvertices in Drefid_to_ldface_to_lvertices
        maxldfaces = max(maxldfaces,length(ldface_to_lvertices))
    end
    maxDfaces = 0
    for vertex in 1:length(vertex_to_Dfaces)
        Dfaces = vertex_to_Dfaces[vertex]
        maxDfaces = max(maxDfaces,length(Dfaces))
    end
    # Allocate output
    ptrs = zeros(Int32,nDfaces+1)
    for Dface in 1:nDfaces
        Drefid = Dface_to_refid[Dface]
        ldface_to_lvertices = Drefid_to_ldface_to_lvertices[Drefid]
        ptrs[Dface+1] = length(ldface_to_lvertices)
    end
    length_to_ptrs!(ptrs)
    ndata = ptrs[end]-1
    data = fill(Int32(INVALID_ID),ndata)
    Dface_to_dfaces = GenericJaggedArray(data,ptrs)
    # Main loop
    Dfaces1 = fill(Int32(INVALID_ID),maxDfaces)
    Dfaces2 = fill(Int32(INVALID_ID),maxDfaces)
    ldfaces1 = fill(Int32(INVALID_ID),maxDfaces)
    nDfaces1 = 0
    nDfaces2 = 0
    newdface = Int32(ndfaces)
    old_to_new = collect(Int32,1:ndfaces)
    for Dface in 1:nDfaces
        Drefid = Dface_to_refid[Dface]
        ldface_to_lvertices = Drefid_to_ldface_to_lvertices[Drefid]
        lvertex_to_vertex = Dface_to_vertices[Dface]
        ldface_to_dface = Dface_to_dfaces[Dface]
        for (ldface,lvertices) in enumerate(ldface_to_lvertices)
            # Do nothing if this local face has already been processed by
            # a neighbor
            if ldface_to_dface[ldface] != Int32(INVALID_ID)
                continue
            end
            # Find if there is already a global d-face for this local d-face
            # if yes, then use the global id of this d-face
            # if not, create a new one
            dface2 = Int32(INVALID_ID)
            fill!(Dfaces1,Int32(INVALID_ID))
            fill!(Dfaces2,Int32(INVALID_ID))
            #vertices = view(lvertex_to_vertex,lvertices)
            for (i,lvertex) in enumerate(lvertices)
                vertex = lvertex_to_vertex[lvertex]
                dfaces = vertex_to_dfaces[vertex]
                for dface1 in dfaces
                    vertices1 = dface_to_vertices[dface1]
                    if same_valid_ids(lvertex_to_vertex,vertices1,lvertices,1:length(vertices1))
                        dface2 = dface1
                        break
                    end
                end
                if dface2 != Int32(INVALID_ID)
                    break
                end
            end
            if dface2 == Int32(INVALID_ID)
                newdface += Int32(1)
                dface2 = newdface
            end
            # Find all D-faces around this local d-face
            for (i,lvertex) in enumerate(lvertices)
                vertex = lvertex_to_vertex[lvertex]
                Dfaces = vertex_to_Dfaces[vertex]
                if i == 1
                    copyto!(Dfaces1,Dfaces)
                    nDfaces1 = length(Dfaces)
                else
                    copyto!(Dfaces2,Dfaces)
                    nDfaces2 = length(Dfaces)
                    intersection!(Dfaces1,Dfaces2,nDfaces1,nDfaces2)
                end
            end
            # Find their correspondent local d-face and set the d-face
            for Dface1 in Dfaces1
                if Dface1 != INVALID_ID
                    Drefid1 = Dface_to_refid[Dface1]
                    lvertex1_to_vertex1 = Dface_to_vertices[Dface1]
                    ldface1_to_lvertices1 = Drefid_to_ldface_to_lvertices[Drefid1]
                    ldface2 = Int32(INVALID_ID)
                    for (ldface1,lvertices1) in enumerate(ldface1_to_lvertices1)
                        #vertices1 = view(lvertex1_to_vertex1,lvertices1)
                        if same_valid_ids(lvertex_to_vertex,lvertex1_to_vertex1,lvertices,lvertices1)
                            ldface2 = ldface1
                            break
                        end
                    end
                    @boundscheck @assert ldface2 != INVALID_ID # TODO: issue gmsh quad eles
                    Dface_to_dfaces[Dface1][ldface2] = dface2
                end
            end
        end # (ldface,lvertices)
    end # Dface
    Dface_to_dfaces, newdface, old_to_new
end

function fill_face_vertices(mesh,d,node_to_vertex)
    function barrier(node_to_vertex,face_to_nodes,face_to_refid,refid_to_lvertex_to_lnodes)
        Ti = eltype(node_to_vertex)
        nfaces = length(face_to_nodes)
        face_to_vertices_ptrs = zeros(Ti,nfaces+1)
        for face in 1:nfaces
            nodes = face_to_nodes[face]
            refid = face_to_refid[face]
            lvertex_to_lnodes = refid_to_lvertex_to_lnodes[refid]
            nlvertices = length(lvertex_to_lnodes)
            face_to_vertices_ptrs[face+1] = nlvertices
        end
        length_to_ptrs!(face_to_vertices_ptrs)
        ndata = face_to_vertices_ptrs[end]-1
        face_to_vertices_data = zeros(Ti,ndata)
        face_to_vertices = JaggedArray(face_to_vertices_data,face_to_vertices_ptrs)
        for face in 1:nfaces
            vertices = face_to_vertices[face]
            nodes = face_to_nodes[face]
            refid = face_to_refid[face]
            lvertex_to_lnodes = refid_to_lvertex_to_lnodes[refid]
            for (lvertex,lnodes) in enumerate(lvertex_to_lnodes)
                lnode = first(lnodes)
                vertex = node_to_vertex[nodes[lnode]]
                @boundscheck @assert vertex != INVALID_ID
                vertices[lvertex] = vertex
            end
        end
        face_to_vertices
    end
    face_to_nodes = face_nodes(mesh,d)
    face_to_refid = face_reference_id(mesh,d)
    refid_to_lvertex_to_lnodes = map(reference_spaces(mesh,d)) do a
        if num_dims(domain(a)) != 0
            face_nodes(complexify(a),0)
        else
            [interior_nodes(a)]
        end
    end
    barrier(node_to_vertex,face_to_nodes,face_to_refid,refid_to_lvertex_to_lnodes)
end

function find_node_to_vertex(mesh)
    function barrier!(node_to_vertex,face_to_nodes,face_to_refid,refid_to_lvertex_to_lnodes)
        valid_id = -one(eltype(node_to_vertex))
        for face in eachindex(face_to_nodes)
            nodes = face_to_nodes[face]
            refid = face_to_refid[face]
            lvertex_to_lnodes = refid_to_lvertex_to_lnodes[refid]
            for lnodes in lvertex_to_lnodes
                lnode = first(lnodes)
                node = nodes[lnode]
                node_to_vertex[node] = valid_id
            end
        end
    end
    Ti = Int32
    nnodes = num_nodes(mesh)
    node_to_vertex = zeros(Ti,nnodes)
    fill!(node_to_vertex,Ti(0))
    D = num_dims(mesh)
    for d in 0:D
        face_to_nodes = face_nodes(mesh,d)
        if length(face_to_nodes) == 0
            continue
        end
        face_to_refid = face_reference_id(mesh,d)
        refid_to_lvertex_to_lnodes = map(reference_spaces(mesh,d)) do a
            if num_dims(domain(a)) != 0
                face_nodes(complexify(a),0)
            else
                [interior_nodes(a)]
            end
        end
        barrier!(node_to_vertex,face_to_nodes,face_to_refid,refid_to_lvertex_to_lnodes)
    end
    valid_id = -one(eltype(node_to_vertex))
    vertex = Ti(0)
    # Loop needed to get the already existing vertices
    for nodes in GT.face_nodes(mesh,0)
        node = nodes[1]
        if node_to_vertex[node] == -1
            vertex += Ti(1)
            node_to_vertex[node] = vertex
        end
    end
    for node in eachindex(node_to_vertex)
        if node_to_vertex[node] == -1
            vertex += Ti(1)
            node_to_vertex[node] = vertex
        end
    end
    node_to_vertex, vertex
end

# The following is a bit ugly but is necessary to dramatically decrease
# the type complexity of meshes and thus of many other types in the code.

function face_topology(boundary::Nothing,vertex_permutations)
    VertexTopology(vertex_permutations)
end

@auto_hash_equals struct VertexTopology{Ti} <: AbstractFaceTopology
    vertex_permutations::Vector{Vector{Ti}}
end
boundary(a::VertexTopology) = nothing
vertex_permutations(a::VertexTopology) = a.vertex_permutations

function face_topology(boundary::MeshTopology,vertex_permutations)
    D = num_dims(boundary) + 1
    if D == 1
        EdgeTopology(boundary,vertex_permutations)
    elseif D == 2
        SurfaceTopology(boundary,vertex_permutations)
    elseif D == 3
        VolumeTopology(boundary,vertex_permutations)
    else
        FaceTopology(boundary,vertex_permutations)
    end
end

@auto_hash_equals struct EdgeTopology{Ti,Tr} <: AbstractFaceTopology
    boundary::MeshTopology{Matrix{JaggedArray{Ti,Ti}},Vector{Vector{Tr}},Matrix{JaggedArray{Ti,Ti}},Tuple{Tuple{VertexTopology{Ti}}}}
    vertex_permutations::Vector{Vector{Ti}}
end
boundary(a::EdgeTopology) = a.boundary
vertex_permutations(a::EdgeTopology) = a.vertex_permutations

@auto_hash_equals struct SurfaceTopology{Ti,Tr} <: AbstractFaceTopology
    boundary::MeshTopology{Matrix{JaggedArray{Ti,Ti}},Vector{Vector{Tr}},Matrix{JaggedArray{Ti,Ti}}, Tuple{Tuple{VertexTopology{Ti}},Tuple{EdgeTopology{Ti,Tr}}}}
    vertex_permutations::Vector{Vector{Ti}}
end
boundary(a::SurfaceTopology) = a.boundary
vertex_permutations(a::SurfaceTopology) = a.vertex_permutations

@auto_hash_equals struct VolumeTopology{Ti,Tr} <: AbstractFaceTopology
    boundary::MeshTopology{Matrix{JaggedArray{Ti,Ti}},Vector{Vector{Tr}},Matrix{JaggedArray{Ti,Ti}}, Tuple{Tuple{VertexTopology{Ti}},Tuple{EdgeTopology{Ti,Tr}},Tuple{SurfaceTopology{Ti,Tr}}}}
    vertex_permutations::Vector{Vector{Ti}}
end
boundary(a::VolumeTopology) = a.boundary
vertex_permutations(a::VolumeTopology) = a.vertex_permutations




topology(a::AbstractTopology) = a
num_dims(t::AbstractTopology) = length(reference_topologies(t))-1
num_faces(t::AbstractTopology,d) = length(face_reference_id(t,d))

struct FaceTopology{A} <: AbstractFaceTopology
    contents::A
end

function face_topology(;boundary,vertex_permutations)
    contents = (;boundary,vertex_permutations)
    FaceTopology(contents)
end

boundary(a::FaceTopology) = a.contents.boundary
vertex_permutations(a::FaceTopology) = a.contents.vertex_permutations

function num_dims(a::FaceTopology)
    if boundary(a) === nothing
        0
    else
        Int(num_dims(boundary(a))) + 1
    end
end

function num_faces(a::FaceTopology,d)
    m = num_dims(a)
    if d == m
        1
    else
        Int(num_faces(boundary(a),d))
    end
end

function face_incidence(a::FaceTopology)
    m = num_dims(a)
    if m==0
        return fill([[1]],1,1)
    end
    fi_boundary = face_incidence(boundary(a))
    T = typeof(fi[1,1])
    fi = Matrix{T}(undef,m+1,m+1)
    for j in 0:(m-1)
        for i in 0:(m-1)
            fi[i+1,j+1] = fi_boundary[i+1,j+1]
        end
        fi[j,end] = fill([1],num_faces(a,j))
        fi[end,j] = collect(1:num_faces(a,j))
    end
    fi
end

function face_incidence(a::FaceTopology,d,D)
    face_incidence(a)[d+1,D+1]
end

function face_reference_id(a::FaceTopology,d)
    m = num_dims(a)
    if m == d
        Int[1]
    else
        collect(Int,face_reference_id(boundary(a),d))
    end
end

function reference_topologies(a::FaceTopology,d)
    m = num_dims(a)
    if m == d
        (a,)
    else
        reference_topologies(boundary(a),d)
    end
end

function face_permutation_ids(a::FaceTopology)
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

function face_permutation_ids(a::FaceTopology,d,D)
    face_permutation_ids(a)[d+1,D+1]
end

struct MeshTopology{A} <: AbstractTopology
    contents::A
end

function mesh_topology(;
        face_incidence,
        face_reference_id,
        face_permutation_ids,
        reference_topologies,)

    contents = (;
                face_incidence,
                face_reference_id,
                face_permutation_ids,
                reference_topologies,)

    MeshTopology(contents)
end

face_incidence(a::MeshTopology) = a.contents.face_incidence
face_incidence(a::MeshTopology,d,D) = a.contents.face_incidence[d+1,D+1]
face_reference_id(a::MeshTopology) = a.contents.face_reference_id
face_reference_id(a::MeshTopology,d) = a.contents.face_reference_id[d+1]
face_permutation_ids(a::MeshTopology) = a.contents.face_permutation_ids
face_permutation_ids(a::MeshTopology,d,D) = a.contents.face_permutation_ids[d+1,D+1]
reference_topologies(a::MeshTopology) = a.contents.reference_topologies
reference_topologies(a::MeshTopology,d) = a.contents.reference_topologies[d+1]

"""
"""
function topology(mesh::AbstractMesh)
    if workspace(mesh) !== nothing
        return workspace(mesh).topology
    end
    # Assumes that the input is a cell complex
    T = JaggedArray{Int32,Int32}
    D = num_dims(mesh)
    my_face_incidence = Matrix{T}(undef,D+1,D+1)
    my_face_reference_id  = [ face_reference_id(mesh,d) for d in 0:D ]
    dims = ntuple(d->d-1,Val(D+1))
    my_reference_faces = map(d->map(topology,reference_domains(mesh,d)),dims)
    my_face_permutation_ids = Matrix{T}(undef,D+1,D+1)
    topo = mesh_topology(;
        face_incidence = my_face_incidence,
        face_reference_id = my_face_reference_id,
        face_permutation_ids = my_face_permutation_ids,
        reference_topologies = my_reference_faces,
       )
    for d in 0:D
        fill_face_interior_mesh_topology!(topo,mesh,d)
    end
    for d in 1:D
        fill_face_vertices_mesh_topology!(topo,mesh,d)
        fill_face_coboundary_mesh_topology!(topo,mesh,d,0)
    end
    for d in 1:(D-1)
        for n in (D-d):-1:1
            m = n+d
            fill_face_boundary_mesh_topology!(topo,mesh,m,n)
            fill_face_coboundary_mesh_topology!(topo,mesh,m,n)
        end
    end
    for d in 0:D
        for n in 0:d
            fill_face_permutation_ids!(topo,d,n)
        end
    end
    topo
end

function fill_face_interior_mesh_topology!(mesh_topology,mesh,d)
    n = num_faces(mesh,d)
    ptrs = collect(Int32,1:(n+1))
    data = collect(Int32,1:n)
    face_incidence(mesh_topology)[d+1,d+1] = JaggedArray(data,ptrs)
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

function fill_face_coboundary_mesh_topology!(topology,mesh,n,m)
    nmfaces = num_faces(mesh,m)
    nface_to_mfaces = face_incidence(topology,n,m)
    face_incidence(topology)[m+1,n+1] = generate_face_coboundary(nface_to_mfaces,nmfaces)
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
    face_incidence(topo)[d+1,0+1] = barrier(nnodes,vertex_to_nodes,dface_to_nodes,dface_to_refid,refid_to_lvertex_to_lnodes)
end

function fill_face_boundary_mesh_topology!(topo,mesh,D,d)
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
                vertices = view(lvertex_to_vertex,lvertices)
                for (i,lvertex) in enumerate(lvertices)
                    vertex = lvertex_to_vertex[lvertex]
                    dfaces = vertex_to_dfaces[vertex]
                    for dface1 in dfaces
                        vertices1 = dface_to_vertices[dface1]
                        if same_valid_ids(vertices,vertices1)
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

                    The given mesh is provably not a cell complex.
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
    face_incidence(topo)[D+1,d+1] = Dface_to_dfaces
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
        face_permutation_ids(top)[D+1,d+1] = cell_to_lface_to_pindex
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
    face_permutation_ids(top)[D+1,d+1] = cell_to_lface_to_pindex
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

function find_eq(v,b)
    for vs in b
        if v == vs
            return true
        end
    end
    return false
end

function is_subset(a,b)
    for i in 1:length(a)
        v = a[i]
        if v == INVALID_ID
            continue
        end
        if !find_eq(v,b)
            return false
        end
    end
    return true
end

function same_valid_ids(a,b)
    if !is_subset(a,b)
        return false
    end
    if !is_subset(b,a)
        return false
    end
    return true
end

function complexify(mesh::AbstractMesh;glue=Val(false))
    Ti = int_type(options(mesh))
    T = JaggedArray{Ti,Ti}
    D = num_dims(mesh)
    oldface_to_newvertices = Vector{T}(undef,D+1)
    newvertex_to_oldfaces = Vector{T}(undef,D+1)
    newface_incidence = Matrix{T}(undef,D+1,D+1)
    nnewfaces = zeros(Int,D+1)
    newface_refid = Vector{Vector{Ti}}(undef,D+1)
    newreffaces = Vector{Any}(undef,D+1)
    newface_nodes = Vector{T}(undef,D+1)
    old_to_new = Vector{Vector{Ti}}(undef,D+1)
    node_to_newvertex, n_new_vertices = find_node_to_vertex(mesh) # Optimizable for linear meshes
    for d in 0:D
        oldface_to_newvertices[d+1] = fill_face_vertices(mesh,d,node_to_newvertex) # Optimizable for linear meshes
        newvertex_to_oldfaces[d+1] = generate_face_coboundary(oldface_to_newvertices[d+1],n_new_vertices) # Optimizable for linear meshes
    end
    newface_incidence[D+1,0+1] = oldface_to_newvertices[D+1]
    newface_incidence[0+1,D+1] = newvertex_to_oldfaces[D+1]
    nnewfaces[D+1] = length(oldface_to_newvertices[D+1])
    newface_refid[D+1] = face_reference_id(mesh,D)
    newreffaces[D+1] = reference_spaces(mesh,D)
    newface_nodes[D+1] = face_nodes(mesh,D)
    old_to_new[D+1] = collect(Ti,1:length(newface_nodes[D+1]))
    # TODO optimize for d==0
    for d in (D-1):-1:0
        n = d+1
        new_nface_to_new_vertices = newface_incidence[n+1,0+1]
        new_vertex_to_new_nfaces = newface_incidence[0+1,n+1]
        old_dface_to_new_vertices = oldface_to_newvertices[d+1]
        new_vertex_to_old_dfaces = newvertex_to_oldfaces[d+1]
        new_nface_to_nrefid = newface_refid[n+1]
        old_dface_to_drefid = face_reference_id(mesh,d)
        drefid_to_ref_dface = reference_spaces(mesh,d)
        old_dface_to_nodes = face_nodes(mesh,d)
        new_nface_to_nodes = newface_nodes[n+1]
        nrefid_to_ldface_to_lvertices = map(a->face_incidence(topology(GT.mesh(domain(a))),d,0),newreffaces[n+1])
        nrefid_to_ldface_to_lnodes = map(a->face_nodes(complexify(a),d),newreffaces[n+1])
        nrefid_to_ldface_to_drefrefid = map(a->face_reference_id(complexify(a),d),newreffaces[n+1])
        nrefid_to_drefrefid_to_ref_dface = map(a->reference_spaces(complexify(a),d),newreffaces[n+1])
        new_nface_to_new_dfaces, n_new_dfaces, old_dface_to_new_dface = generate_face_boundary(
            new_nface_to_new_vertices,
            new_vertex_to_new_nfaces,
            old_dface_to_new_vertices,
            new_vertex_to_old_dfaces,
            new_nface_to_nrefid,
            nrefid_to_ldface_to_lvertices)
        new_dface_to_new_vertices = generate_face_vertices(
            n_new_dfaces,
            old_dface_to_new_dface,
            old_dface_to_new_vertices,
            new_nface_to_new_vertices,
            new_nface_to_new_dfaces,
            new_nface_to_nrefid,
            nrefid_to_ldface_to_lvertices)
        new_vertex_to_new_dfaces = generate_face_coboundary(new_dface_to_new_vertices,n_new_vertices)
        new_dface_to_nodes = generate_face_vertices(
            n_new_dfaces,
            old_dface_to_new_dface,
            old_dface_to_nodes,
            new_nface_to_nodes,
            new_nface_to_new_dfaces,
            new_nface_to_nrefid,
            nrefid_to_ldface_to_lnodes)
        new_dface_to_new_drefid, new_refid_to_ref_dface = generate_reference_spaces(
            n_new_dfaces,
            old_dface_to_new_dface,
            old_dface_to_drefid,
            drefid_to_ref_dface,
            new_nface_to_new_dfaces,
            new_nface_to_nrefid,
            nrefid_to_ldface_to_drefrefid,
            nrefid_to_drefrefid_to_ref_dface)
        newface_incidence[n+1,d+1] = new_nface_to_new_dfaces
        newface_incidence[d+1,0+1] = new_dface_to_new_vertices
        newface_incidence[0+1,d+1] = new_vertex_to_new_dfaces
        newface_refid[d+1] = new_dface_to_new_drefid
        newreffaces[d+1] = new_refid_to_ref_dface
        nnewfaces[d+1] = n_new_dfaces
        newface_nodes[d+1] = new_dface_to_nodes
        old_to_new[d+1] = old_dface_to_new_dface
    end
    node_to_coords = node_coordinates(mesh)
    old_physical_faces = physical_faces(mesh)
    new_physical_faces = [ Dict{String,Vector{Int32}}() for d in 0:D] # TODO hardcoded
    for d in 0:D
        old_groups = old_physical_faces[d+1]
        for (group_name,old_group_faces) in old_groups
            new_group_faces = similar(old_group_faces)
            new_group_faces .= old_to_new[d+1][old_group_faces]
            new_physical_faces[d+1][group_name] = new_group_faces
        end
    end
    new_mesh = GT.mesh(;
            node_coordinates = node_to_coords,
            face_nodes = newface_nodes,
            face_reference_id = newface_refid,
            reference_spaces = Tuple(newreffaces),
            physical_faces = new_physical_faces,
            periodic_nodes = periodic_nodes(mesh),
            outward_normals = outward_normals(mesh),
            geometry_names = geometry_names(mesh),
            is_cell_complex = Val(true),
           )
    mtopology = GT.topology(new_mesh)
    workspace = (;topology=mtopology)
    new_mesh2 = replace_workspace(new_mesh,workspace)
    if val_parameter(glue)
        new_mesh2, old_to_new
    else
        new_mesh2
    end
end

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

function generate_reference_spaces(
        n_new_dfaces,
        old_dface_to_new_dface,
        old_dface_to_drefid,
        drefid_to_ref_dface,
        new_nface_to_new_dfaces,
        new_nface_to_nrefid,
        nrefid_to_ldface_to_drefrefid,
        nrefid_to_drefrefid_to_ref_dface)

    i_to_ref_dface = collect(Any,drefid_to_ref_dface)
    drefid_to_i = collect(1:length(drefid_to_ref_dface))
    i = length(i_to_ref_dface)
    Ti = Int32
    nrefid_to_drefrefid_to_i = map(a->zeros(Ti,length(a)),nrefid_to_drefrefid_to_ref_dface)
    for (nrefid,drefrefid_to_ref_dface) in enumerate(nrefid_to_drefrefid_to_ref_dface)
        for (drefrefid, ref_dface) in enumerate(drefrefid_to_ref_dface)
            push!(i_to_ref_dface,ref_dface)
            i += 1
            nrefid_to_drefrefid_to_i[nrefid][drefrefid] = i
        end
    end
    u_to_ref_dface = unique(i_to_ref_dface)
    i_to_u = indexin(i_to_ref_dface,u_to_ref_dface)
    new_dface_to_u = zeros(Ti,n_new_dfaces)
    new_dface_to_touched = fill(false,n_new_dfaces)
    for (old_dface,new_dface) in enumerate(old_dface_to_new_dface)
        drefid = old_dface_to_drefid[old_dface]
        i = drefid_to_i[drefid]
        u = i_to_u[i]
        new_dface_to_u[new_dface] = u
        new_dface_to_touched[new_dface] = true
    end
    for (new_nface,nrefid) in enumerate(new_nface_to_nrefid)
        ldface_to_drefrefid = nrefid_to_ldface_to_drefrefid[nrefid]
        ldface_to_new_dface = new_nface_to_new_dfaces[new_nface]
        drefrefid_to_i = nrefid_to_drefrefid_to_i[nrefid]
        for (ldface,new_dface) in enumerate(ldface_to_new_dface)
            new_dface = ldface_to_new_dface[ldface]
            if new_dface_to_touched[new_dface]
                continue
            end
            drefrefid = ldface_to_drefrefid[ldface]
            i = drefrefid_to_i[drefrefid]
            u = i_to_u[i]
            new_dface_to_u[new_dface] = u
            new_dface_to_touched[new_dface] = true
        end
    end
    new_dface_to_u, Tuple(u_to_ref_dface)
end

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
            vertices = view(lvertex_to_vertex,lvertices)
            for (i,lvertex) in enumerate(lvertices)
                vertex = lvertex_to_vertex[lvertex]
                dfaces = vertex_to_dfaces[vertex]
                for dface1 in dfaces
                    vertices1 = dface_to_vertices[dface1]
                    if same_valid_ids(vertices,vertices1)
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
                        vertices1 = view(lvertex1_to_vertex1,lvertices1)
                        if same_valid_ids(vertices,vertices1)
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
        valid_id = one(eltype(node_to_vertex))
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
    fill!(node_to_vertex,Ti(INVALID_ID))
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
    vertex = Ti(0)
    for node in eachindex(node_to_vertex)
        if node_to_vertex[node] != Ti(INVALID_ID)
            vertex += Ti(1)
            node_to_vertex[node] = vertex
        end
    end
    node_to_vertex, vertex
end



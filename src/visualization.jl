

function visualization_mesh(mesh::AbstractMesh,dim,ids=num_faces(mesh,dim);order=nothing,refinement=nothing)
    function barrier(
            refid_to_tabulation,
            refid_to_scell_to_snodes,
            refid_to_scell_to_srefid,
            refid_to_srefid_to_oid,
            refid_to_srefid_to_vrefface,
            refid_to_snode_to_coords,
            node_to_coords,
            cell_to_nodes,
            cell_to_refid,
            ::Val{Dn}) where Dn

        ncells = length(cell_to_refid)
        nvnodes = 0
        nvcells = 0
        for cell in 1:ncells
            refid = cell_to_refid[cell]
            nvnodes += length(refid_to_snode_to_coords[refid])
            nvcells += length(refid_to_scell_to_srefid[refid])

        end
        nrefids = length(refid_to_srefid_to_oid)
        i_to_oid = reduce(vcat,refid_to_srefid_to_oid)
        i_to_vrefface = reduce(vcat,refid_to_srefid_to_vrefface)
        refid_to_srefid_to_i = Vector{Vector{Int}}(undef,nrefids)
        i = 0
        for refid in 1:nrefids
            srefid_to_oid = refid_to_srefid_to_oid[refid]
            nsrefids = length(srefid_to_oid)
            srefid_to_i = zeros(Int,nsrefids)
            for srefid in 1:nsrefids
                i += 1
                srefid_to_i[srefid] = i
            end
            refid_to_srefid_to_i[refid] = srefid_to_i
        end
        vrefid_to_oid = unique(i_to_oid)
        i_to_vrefid = indexin(i_to_oid,vrefid_to_oid)
        vrefid_to_i = indexin(vrefid_to_oid,i_to_oid)
        vrefid_to_vrefface = i_to_vrefface[vrefid_to_i]
        Tx = SVector{Dn,Float64}
        vnode_to_coords = zeros(Tx,nvnodes)
        vcell_to_vnodes_ptrs = zeros(Int32,nvcells+1)
        vcell_to_vrefid = zeros(Int32,nvcells)
        vcell_to_cell = zeros(Int32,nvcells)
        cell_to_vnodes = fill(0:1,ncells)
        vcell = 0
        for cell in 1:ncells
            refid = cell_to_refid[cell]
            scell_to_snodes = refid_to_scell_to_snodes[refid]
            nscells = length(scell_to_snodes)
            scell_to_srefid = refid_to_scell_to_srefid[refid]
            srefid_to_i = refid_to_srefid_to_i[refid]
            for scell in 1:nscells
                srefid = scell_to_srefid[scell]
                i = srefid_to_i[srefid]
                vrefid = i_to_vrefid[i]
                snodes = scell_to_snodes[scell]
                vcell += 1
                vcell_to_vnodes_ptrs[vcell+1] = length(snodes)
                vcell_to_vrefid[vcell] = vrefid
                vcell_to_cell[vcell] = cell
            end
        end
        length_to_ptrs!(vcell_to_vnodes_ptrs)
        ndata = vcell_to_vnodes_ptrs[end]-1
        vcell_to_vnodes_data = zeros(Int32,ndata)
        vcell = 0
        vnode = 0
        vnode_prev = 1
        for cell in 1:ncells
            refid = cell_to_refid[cell]
            scell_to_snodes = refid_to_scell_to_snodes[refid]
            nscells = length(scell_to_snodes)
            for scell in 1:nscells
                snodes = scell_to_snodes[scell]
                vcell += 1
                p = vcell_to_vnodes_ptrs[vcell]
                for (i,snode) in enumerate(snodes)
                    vcell_to_vnodes_data[p-1+i] = snode + vnode
                end
            end
            tabulation = refid_to_tabulation[refid]
            nsnodes = size(tabulation,1)
            nodes = cell_to_nodes[cell]
            for snode in 1:nsnodes
                y = zero(Tx)
                for (i,node) in enumerate(nodes)
                    coeff = tabulation[snode,i]
                    x = node_to_coords[node]
                    y += coeff*x
                end
                vnode += 1
                vnode_to_coords[vnode] = y
            end
            cell_to_vnodes[cell] = vnode_prev:vnode
            vnode_prev = vnode + 1
        end
        vcell_to_vnodes = JaggedArray(vcell_to_vnodes_data,vcell_to_vnodes_ptrs)
        vchain = GT.chain(;
                        node_coordinates = vnode_to_coords,
                        face_nodes = vcell_to_vnodes,
                        face_reference_id = vcell_to_vrefid,
                        reference_spaces = vrefid_to_vrefface)
        vmesh = GT.mesh(vchain)
        vglue = (;parent_face=vcell_to_cell,
                 reference_coordinates=refid_to_snode_to_coords,
                 face_fine_nodes = cell_to_vnodes,
                 num_dims=Val(dim))
        vmesh, vglue
    end # barrier
    refid_to_refface = reference_spaces(mesh,dim)
    refid_to_refmesh = map(refid_to_refface) do ref_face
        if order === nothing && refinement === nothing
            # Use the given cells as visualization cells
            complexify(ref_face)
        elseif order !== nothing && refinement === nothing
            # Use cells of given order as visualization cells
            geo = GT.domain(ref_face)
            ref_face_ho = lagrange_space(geo,order)
            complexify(ref_face_ho)
        elseif order === nothing && refinement !== nothing
            # Use linear sub-cells with $refinement per direction per direction
            geom = GT.domain(ref_face)
            refine_reference_geometry(geom,refinement)
        else
            error("order and refinement kw-arguments can not be given at the same time")
        end
    end
    refid_to_tabulation = map(refid_to_refface,refid_to_refmesh) do refface,refmesh
        x = node_coordinates(refmesh)
        tabulator(refface)(value,x)
    end
    refid_to_scell_to_snodes = map(refmesh->face_nodes(refmesh,dim),refid_to_refmesh)
    refid_to_scell_to_srefid = map(refmesh->face_reference_id(refmesh,dim),refid_to_refmesh)
    refid_to_srefid_to_oid = map(refmesh->collect(map(objectid,reference_spaces(refmesh,dim))),refid_to_refmesh)
    refid_to_srefid_to_vrefface = map(refmesh->reference_spaces(refmesh,dim),refid_to_refmesh)
    refid_to_snode_to_coords = map(node_coordinates,refid_to_refmesh)
    node_to_coords = node_coordinates(mesh)
    cell_to_nodes = view(face_nodes(mesh,dim),ids)
    cell_to_refid = view(face_reference_id(mesh,dim),ids)
    Dn = num_ambient_dims(mesh)
    barrier(
            refid_to_tabulation,
            refid_to_scell_to_snodes,
            refid_to_scell_to_srefid,
            refid_to_srefid_to_oid,
            refid_to_srefid_to_vrefface,
            refid_to_snode_to_coords,
            node_to_coords,
            cell_to_nodes,
            cell_to_refid,
            Val(Dn))
end

function refine_reference_geometry(geo,resolution)
    function refine_n_cube_aligned(geo,n)
        box = bounding_box(geo)
        domain = domain_from_bounding_box(box)
        D = num_dims(geo)
        cells = ntuple(i->n,Val(D))
        cartesian_mesh(domain,cells)
    end
    function refine_unit_triangle(geo,n)
        # Copyed + adapted from Gridap
        tri_num(n) = n*(n+1)÷2
        v(n,i,j) = tri_num(n) - tri_num(n-i+1) + j
        D = 2
        quad_to_tris = ((1,2,3),(2,4,3))
        quad = CartesianIndices( (0:1,0:1) )
        Tp = SVector{2,Float64}
        n_verts = tri_num(n+1)
        n_cells = tri_num(n)+tri_num(n-1)
        n_verts_x_cell = 3
        X = zeros(Tp,n_verts)
        T = [ zeros(Int,n_verts_x_cell) for i in 1:n_cells ]
        for i in 1:n+1
          for j in 1:n+1-i+1
            vert = v(n+1,i,j)
            X[vert] = SVector((i-1)/n,(j-1)/n)
          end
        end
        for i in 1:n
          for j in 1:n-(i-1)
            verts = ntuple( lv-> v(n+1, (i,j).+quad[lv].I ...), Val{2^D}() )
            cell = v(n,i,j)
            T[cell] .= map(i->verts[i],quad_to_tris[1])
            if (i-1)+(j-1) < n-1
              cell = tri_num(n) + v(n-1,i,j)
              T[cell] .= map(i->verts[i],quad_to_tris[2])
            end
          end
        end
        refface = lagrange_mesh_face(geo,1)
        chain = chain_from_arrays(
                       X,
                       T,
                       fill(1,length(T)),
                       [refface]
                      )
        mesh_from_chain(chain)
    end
    function refine_unit_tet(geo,n)
        # Copyed + adapted from Gridap
        tri_num(n) = n*(n+1)÷2
        tet_num(n) = n*(n+1)*(n+2)÷6
        v(n,i,j) = tri_num(n) - tri_num(n-i+1) + j
        v(n,i,j,k) = tet_num(n) - tet_num(n-i+1) + v(n-i+1,j,k)
        D = 3
        cube_to_tets = ((1,2,3,5),(2,4,3,6),(3,5,7,6),(2,3,5,6),(3,4,7,6),(4,6,7,8))
        cube = CartesianIndices( (0:1,0:1,0:1) )
        n_core_tets = length(cube_to_tets)-2
        Tp = SVector{3,Float64}
        n_verts = tet_num(n+1)
        n_cells = tet_num(n)+n_core_tets*tet_num(n-1)+tet_num(n-2)
        n_verts_x_cell = 4
        X = zeros(Tp,n_verts)
        T = [ zeros(Int,n_verts_x_cell) for i in 1:n_cells ]
        for i in 1:n+1
          for j in 1:n+1-(i-1)
            for k in 1:n+1-(i-1)-(j-1)
              vert = v(n+1,i,j,k)
              X[vert] = SVector((i-1)/n,(j-1)/n,(k-1)/n)
            end
          end
        end
        for i in 1:n
          for j in 1:n-(i-1)
            for k in 1:n-(i-1)-(j-1)
              verts = ntuple( lv-> v(n+1, (i,j,k).+cube[lv].I ...), Val{2^D}() )
              cell = v(n,i,j,k)
              T[cell] .= map(i->verts[i],cube_to_tets[1])
              if (i-1)+(j-1)+(k-1) < n-1
                cell = tet_num(n) + (v(n-1,i,j,k)-1)*n_core_tets
                for t in 1:n_core_tets
                  T[cell+t] .= map(i->verts[i],cube_to_tets[t+1])
                end
              end
              if (i-1)+(j-1)+(k-1) < n-2
                cell = tet_num(n) + n_core_tets*tet_num(n-1) + v(n-2,i,j,k)
                T[cell] .= map(i->verts[i],cube_to_tets[end])
              end
            end
          end
        end
        refface = lagrange_mesh_face(geo,1)
        chain = chain_from_arrays(
                       X,
                       T,
                       fill(1,length(T)),
                       [refface],
                      )
        mesh_from_chain(chain)
    end
    if is_n_cube(geo) && is_axis_aligned(geo)
        refine_n_cube_aligned(geo,resolution)
    elseif is_unit_simplex(geo) && num_dims(geo) == 2
        refine_unit_triangle(geo,resolution)
    elseif is_unit_simplex(geo) && num_dims(geo) == 3
        refine_unit_tet(geo,resolution)
    else
        error("Case not implemented (yet)")
    end
end





struct Plot{A,B,C,D} <: AbstractType
    mesh::A
    face_data::B
    node_data::C
    cache::D
    function Plot(mesh,face_data,node_data,cache=nothing)
        A = typeof(mesh)
        B = typeof(face_data)
        C = typeof(node_data)
        D = typeof(cache)
        new{A,B,C,D}(mesh,face_data,node_data,cache)
    end
end

function visualization_mesh(plt::Plot;glue=Val(false))
    vis_mesh = plt.mesh
    vis_glue = plt.cache.glue
    if val_parameter(glue)
        vis_mesh, vis_glue
    else
        vis_mesh
    end
end

function replace_mesh(plt::Plot,mesh)
    Plot(mesh,plt.face_data,plt.node_data,plt.cache)
end

struct PPlot{A,B} <: AbstractType
    partition::A
    cache::B
    function PPlot(p,cache=nothing)
        A = typeof(p)
        B = typeof(cache)
        new{A,B}(p,cache)
    end
end

function num_ambient_dims(plt::Plot)
    num_ambient_dims(plt.mesh)
end

function num_ambient_dims(plt::PPlot)
    PartitionedArrays.getany(map(num_ambient_dims,plt.partition))
end

function PartitionedArrays.partition(plt::PPlot)
    plt.partition
end

function simplexify(plt::PPlot)
    PPlot(map(simplexify,plt.partition))
end

function face_data(plt::Plot,d;merge_dims=Val(false))
    if val_parameter(merge_dims)
        mesh = plt.mesh
        dict = face_data(plt;merge_dims)
        offsets = face_offset(mesh)
        offset = offsets[d+1]
        faces = 1:num_faces(mesh,d)
        for (k,v) in dict
            dict[k] = v[faces.+offset]
        end
        dict
    else
        plt.face_data[d+1]
    end
end

function face_data(plt::Plot;merge_dims=Val(false))
    if val_parameter(merge_dims)
        mesh = plt.mesh
        nfaces = sum(num_faces(mesh))
        offsets = face_offset(mesh)
        D = num_dims(mesh)
        dict = Dict{String,Any}()
        for d in 0:D
            offset = offsets[d+1]
            for group in plt.face_data[d+1]
                name,data = group
                if !haskey(dict,name)
                    face_mask = zeros(eltype(data),nfaces)
                    dict[name] = face_mask
                end
                face_mask = dict[name]
                faces = 1:num_faces(mesh,d)
                face_mask[faces.+offset] .= data
            end
        end
        dict
    else
        plt.face_data
    end
end

node_data(plt::Plot) = plt.node_data

function face_color(plt::Plot,name,d)
    mesh = plt.mesh
    allcolors = face_data(plt;merge_dims=true)[name]
    offset = face_offset(mesh,d)
    ndfaces = num_faces(mesh,d)
    dfaces = 1:ndfaces
    dface_to_color = allcolors[dfaces .+ offset]
    node_to_color = similar(dface_to_color,num_nodes(mesh))
    dface_to_nodes = face_nodes(mesh,d)
    for face in 1:ndfaces
        nodes = dface_to_nodes[face]
        node_to_color[nodes] .= dface_to_color[face]
    end
    color = node_to_color
    color
end

struct FaceData
    name::String
end

struct NodeData
    name::String
end

function face_colorrange(plt::Plot,name)
    mesh = plt.mesh
    allcolors = face_data(plt;merge_dims=true)[name]
    minc = minimum(allcolors)
    maxc = maximum(allcolors)
    colorrange = (minc,maxc)
end

function node_color(plt::Plot,name)
    node_data(plt)[name]
end

function node_colorrange(plt::Plot,name)
    mesh = plt.mesh
    allcolors = node_data(plt)[name]
    minc = minimum(allcolors)
    maxc = maximum(allcolors)
    colorrange = (minc,maxc)
end

function complexify(plt::Plot)
    mesh = plt.mesh
    tmesh,old_to_new = complexify(mesh;glue=Val(true))
    D = num_dims(mesh)
    fdata = [Dict{String,Any}() for d in 0:D]
    for d in 0:D
        groups = plt.face_data[d+1]
        tgroups = fdata[d+1]
        nfaces = num_faces(mesh,d)
        ntfaces = num_faces(tmesh,d)
        face_to_tface = old_to_new[d+1]
        tface_to_face = zeros(Int32,ntfaces)
        tface_to_face[face_to_tface] = 1:nfaces
        for (k,v) in groups
            tface_to_mask = similar(v,ntfaces)
            fill!(tface_to_mask,zero(eltype(tface_to_mask)))
            for tface in 1:ntfaces
                face = tface_to_face[tface]
                if face != 0
                    tface_to_mask[tface] = v[face]
                end
            end
            tgroups[k] = tface_to_mask
        end
    end
    Plot(tmesh,fdata,plt.node_data)
end

function simplexify(plt::Plot)
    plt = complexify(plt)
    mesh = plt.mesh
    tmesh, d_to_tface_to_face  = simplexify(mesh;glue=Val(true))
    D = num_dims(tmesh)
    fdata = [Dict{String,Any}() for d in 0:D]
    for d in 0:D
        groups = plt.face_data[d+1]
        tgroups = fdata[d+1]
        nfaces = num_faces(mesh,d)
        ntfaces = num_faces(tmesh,d)
        tface_to_face = d_to_tface_to_face[d+1]
        for (k,v) in groups
            tface_to_mask = similar(v,ntfaces)
            fill!(tface_to_mask,zero(eltype(tface_to_mask)))
            for tface in 1:ntfaces
                face = tface_to_face[tface]
                if face != 0
                    tface_to_mask[tface] = v[face]
                end
            end
            tgroups[k] = tface_to_mask
        end
    end
    Plot(tmesh,fdata,plt.node_data)
end

function plot(mesh::AbstractMesh)
    Plot(mesh,face_data(mesh),node_data(mesh))
end

function face_data(mesh::AbstractMesh,d)
    ndfaces = num_faces(mesh,d)
    dict = Dict{String,Vector{Int32}}()
    for group in physical_faces(mesh,d)
        name,faces = group
        face_mask = zeros(Int32,ndfaces)
        face_mask[faces] .= 1
        dict[name] = face_mask
    end
    #dict["__LOCAL_DFACE__"] = collect(1:ndfaces)
    dict
end

function node_data(mesh::AbstractMesh)
    nnodes = num_nodes(mesh)
    dict = Dict{String,Vector{Int32}}()
    #for group in physical_nodes(mesh;merge_dims=Val(true))
    #    name,nodes = group
    #    node_mask = zeros(Int32,nnodes)
    #    node_mask[nodes] .= 1
    #    dict[name] = node_mask
    #end
    #dict["__LOCAL_NODE__"] = collect(1:nnodes)
    dict
end

function face_data(mesh::AbstractMesh)
    D = num_dims(mesh)
    map(d->face_data(mesh,d),0:D)
end

function face_data(mesh::AbstractMesh,ids::PMeshLocalIds,d)
    facedata = face_data(mesh,d)
    faceids = face_indices(ids,d)
    part = part_id(faceids)
    nfaces = local_length(faceids)
    facedata["__OWNER__"] = local_to_owner(faceids)
    facedata["__IS_LOCAL__"] = local_to_owner(faceids) .== part
    facedata["__PART__"] = fill(part,nfaces)
    #facedata["__GLOBAL_DFACE__"] = local_to_global(faceids)
    facedata
end

function face_data(mesh::AbstractMesh,ids::PMeshLocalIds)
    D = num_dims(mesh)
    map(d->face_data(mesh,ids,d),0:D)
end

function node_data(mesh::AbstractMesh,ids::PMeshLocalIds)
    nodedata = node_data(mesh)
    nodeids = node_indices(ids)
    part = part_id(nodeids)
    nfaces = local_length(nodeids)
    nodedata["__OWNER__"] = local_to_owner(nodeids)
    nodedata["__IS_LOCAL__"] = local_to_owner(nodeids) .== part
    nodedata["__PART__"] = fill(part,nfaces)
    #nodedata["__GLOBAL_NODE__"] = local_to_global(nodeids)
    nodedata
end

function plot(pmesh::PMesh)
    plts = map(partition(pmesh),index_partition(pmesh)) do mesh,ids
        Plot(mesh,face_data(mesh,ids),node_data(mesh,ids))
    end
    PPlot(plts)
end

function restrict_to_dim(mesh::AbstractMesh,d)
    chain = GT.chain(;
    node_coordinates = node_coordinates(mesh),
    face_nodes = face_nodes(mesh,d),
    face_reference_id = face_reference_id(mesh,d),
    reference_spaces = reference_spaces(mesh,d),
    periodic_nodes = periodic_nodes(mesh),
    physical_faces = physical_faces(mesh,d))
    GT.mesh(chain)
end

function restrict_to_dim(plt::Plot,d)
    mesh = restrict_to_dim(plt.mesh,d)
    dfacedata = face_data(plt,d;merge_dims=true)
    fd = map(0:d) do i
        if i == d
            dfacedata
        else
            typeof(dfacedata)()
        end
    end
    pltd = Plot(mesh,fd,node_data(plt))
    pltd
end

function restrict(plt::Plot,newnodes,newfaces)
    mesh = restrict(plt.mesh,newnodes,newfaces)
    nodedata = node_data(plt)[newnodes]
    facedata = map(face_data(plt),newfaces) do data,nfs
        dict = copy(data)
        for (k,v) in dict
            dict[k] = dict[k][nfs]
        end
        dict
    end
    Plot(mesh,facedata,nodedata)
end

function shrink(plt::Plot;scale=0.75)
    mesh = plt.mesh
    D = num_dims(mesh)
    nnewnodes = 0
    for d in 0:D
        ndfaces = num_faces(mesh,d)
        dface_to_nodes = face_nodes(mesh,d)
        for dface in 1:ndfaces
            nodes = dface_to_nodes[dface]
            nnewnodes += length(nodes)
        end
    end
    node_to_x = node_coordinates(mesh)
    newnode_to_x = similar(node_to_x,nnewnodes)
    face_newnodes = map(deepcopy,face_nodes(mesh)) # copy not properly implemented for JaggedArray?
    newnode = 0
    newnode_to_node = zeros(Int32,nnewnodes)
    for d in 0:D
        ndfaces = num_faces(mesh,d)
        dface_to_nodes = face_nodes(mesh,d)
        dface_to_newnodes = face_newnodes[d+1]
        for dface in 1:ndfaces
            nodes = dface_to_nodes[dface]
            xm = sum(node->node_to_x[node],nodes) / length(nodes)
            newnodes = dface_to_newnodes[dface]
            for (i,node) in enumerate(nodes)
                newnode += 1
                x = node_to_x[node]
                newnodes[i] = newnode
                newnode_to_x[newnode] = scale*(x-xm) + xm
                newnode_to_node[newnode] = node
            end
        end
    end
    new_mesh = GT.mesh(;
            node_coordinates = newnode_to_x,
            face_nodes = face_newnodes,
            face_reference_id = face_reference_id(mesh),
            reference_spaces = reference_spaces(mesh),
            physical_faces = physical_faces(mesh), # TODO propagate also periodic nodes??
            outward_normals = outward_normals(mesh)
           )
    newnode_data = Dict{String,Any}()
    for (k,node_to_val) in node_data(plt)
        newnode_data[k] = node_to_val[newnode_to_node]
    end
    Plot(new_mesh,plt.face_data,newnode_data)
end

function shrink(pplt::PPlot;scale=0.75)
    plts = map(partition(pplt)) do plt
        shrink(plt;scale)
    end
    PPlot(plts)
end

function warp_by_vector(plt::Plot,vec::NodeData;scale=1)
    mesh = plt.mesh
    node_to_x = node_coordinates(mesh)
    node_to_vec = plt.node_data[vec.name]
    node_to_x = node_to_x .+ scale .* node_to_vec
    mesh2 = replace_node_coordinates(mesh,node_to_x)
    plt2 = replace_mesh(plt,mesh2)
    plt2
end

function plot(domain::AbstractDomain;kwargs...)
    mesh = GT.mesh(domain)
    d = GT.num_dims(domain)
    domface_to_face = GT.faces(domain)
    vismesh = GT.visualization_mesh(mesh,d,domface_to_face;kwargs...)
    node_data = Dict{String,Any}()
    #quad = PlotQuadrature(mesh,domain,vismesh)
    #foreach(pairs(fields)) do (sym,field)
    #    node_data[string(sym)] = plot_field(quad,field)
    #end
    vmesh,glue = vismesh
    cache = (;glue,domain)
    Plot(vmesh,face_data(vmesh),node_data,cache)
end

function plot(domain::AbstractMeshDomain{<:PMesh};kwargs...)
    mesh = GT.mesh(domain)
    plts = map(partition(domain)) do mydom
        plt = plot(mydom;kwargs...)
    end
    cache = (;domain)
    PPlot(plts,cache)
end

function plot!(field,plt::Plot;label)
    plot!(plt,field;label)
end

function plot!(plt::Plot,field;label)
    q = GT.coordinates(plt)
    f_q = field(q)
    term = GT.term(f_q)
    plot_impl!(plt,term,label)
end

function plot!(plt::PPlot,field;label)
    q = GT.coordinates(plt)
    f_q = field(q)
    term = GT.term(f_q)
    map(partition(plt),term) do myplt, myterm
        plot_impl!(myplt,myterm,label)
    end
    plt
end

function plot_impl!(plt,term,label)
    vmesh  = plt.mesh
    vglue = plt.cache.glue
    nnodes = GT.num_nodes(vmesh)
    index = GT.generate_index(domain(plt))
    d = target_dim(index)
    face_to_nodes = GT.get_symbol!(index,vglue.face_fine_nodes,"face_to_nodes")
    t = term(index)
    T = typeof(prototype(t))
    data = zeros(T,nnodes)
    expr1 = expression(t) |> simplify
    face = face_index(index,d)
    point = point_index(index)
    deps = (face,point)
    v = GT.topological_sort(expr1,deps)
    expr = quote
        (data,state) -> begin
            $(unpack_index_storage(index,:state))
            $(v[1])
            for $face in 1:length($face_to_nodes)
                nodes = $face_to_nodes[$face]
                $(v[2])
                for $point in 1:length(nodes)
                    data[nodes[$point]] = $(v[3])
                end
            end
        end
    end
    filldata! = eval(expr)
    invokelatest(filldata!,data,GT.index_storage(index))
    plt.node_data[label] = data
    plt
end

domain(plt::Union{Plot,PPlot}) = plt.cache.domain

function coordinates(plt::Union{Plot,PPlot})
    domain = plt |> GT.domain
    q = GT.reference_coordinates(plt)
    if is_reference_domain(domain)
        return q
    end
    d = num_dims(domain(plt))
    phi = physical_map(mesh(domain(plt)),d)
    phi(q)
end

#function coordinates(plt::Union{Plot,PPlot},::ReferenceDomain)
#    GT.reference_coordinates(plt)
#end
#
#function coordinates(plt::Union{Plot,PPlot},::PhysicalDomain)
#    domain_phys = plt |> GT.domain
#    d = num_dims(domain(plt))
#    phi = physical_map(mesh(domain(plt)),d)
#    q = GT.reference_coordinates(plt)
#    phi(q)
#    #domain_ref = domain_phys |> reference_domain
#    #phi = GT.domain_map(domain_ref,domain_phys)
#    #q = GT.reference_coordinates(plt)
#    #phi(q)
#end

function reference_coordinates(plt::Plot)
    domain = GT.reference_domain(plt.cache.domain)
    d = GT.num_dims(domain)
    domface_to_face = GT.faces(domain)
    mesh = GT.mesh(domain)
    vmesh, vglue = GT.visualization_mesh(plt,glue=true)
    refid_to_snode_to_coords = vglue.reference_coordinates
    d = GT.num_dims(vmesh)
    face_reference_id = GT.face_reference_id(mesh,d)
    GT.point_quantity(refid_to_snode_to_coords,domain;reference=true,face_reference_id)
    #prototype = first(first(refid_to_snode_to_coords))
    #GT.quantity(prototype,domain) do index
    #    domface = index.face
    #    point = index.point
    #    dict = index.dict
    #    domface_to_face_sym = get!(dict,domface_to_face,gensym("domface_to_face"))
    #    face_to_refid_sym = get!(dict,face_to_refid,gensym("face_to_refid"))
    #    refid_to_snode_to_coords_sym = get!(dict,refid_to_snode_to_coords,gensym("refid_to_snode_to_coords"))
    #    @term begin
    #        face = $domface_to_face_sym[$domface]
    #        coords = reference_value(
    #            $refid_to_snode_to_coords_sym,
    #            $face_to_refid_sym,
    #            $face)
    #        coords[$point]
    #    end
    #end
end

function reference_coordinates(plt::PPlot)
    q = map(GT.reference_coordinates,partition(plt))
    term = map(GT.term,q)
    prototype = map(GT.prototype,q) |> PartitionedArrays.getany
    GT.quantity(term,prototype,plt.cache.domain)
end

# VTK

struct VTKPlot{A,B} <: AbstractType
    plot::A
    vtk::B
end

struct PVD{A} <: AbstractType
    pvd::A
end

struct PPVD{A} <: AbstractType
    partition::A
end

function vtk_points end
function vtk_points! end
function vtk_cells end
function vtk_cells! end
function vtk_args end
function vtk_args! end
function vtk_physical_faces end
function vtk_physical_faces! end
#function vtk_physical_nodes end
#function vtk_physical_nodes! end
function vtk_mesh_cell end
function vtk_mesh_cell! end
function translate_vtk_data end
function translate_vtk_data! end
function vtk_close_impl end
function vtk_close_impl! end

# Makie prototype functions to be defined inside Makie's extension module; see ext/GalerkinToolkitMakieExt.jl

function makieplot end
function makieplot! end
function makie0d end
function makie0d! end
function makie1d end
function makie1d! end
function makie2d end
function makie2d! end
function makie2d1d end
function makie2d1d! end
function makie3d end
function makie3d! end
function makie3d1d end
function makie3d1d! end

# handler for displaying a hint in case the user tries to call functions defined in the extension module
# https://github.com/JuliaLang/julia/blob/b9d9b69165493f6fc03870d975be05c67f14a30b/base/errorshow.jl#L1039
#
#function __init__()
#    @static if !isdefined(Base, :get_extension)
#        @require Makie="ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" begin
#            include("../ext/GalerkinToolkitMakieExt.jl")
#        end
#    end
#    Base.Experimental.register_error_hint(MethodError) do io, exc, argtypes, kwargs
#        if exc.f in [makieplot, makieplot!, makie0d, makie0d!, makie1d, makie1d!, makie2d, makie2d!, makie2d1d, makie2d1d!, makie3d, makie3d!, makie3d1d, makie3d1d!]
#            if isempty(methods(exc.f))
#                print(io, "\n$(exc.f) has no methods, yet. Makie has to be loaded for the plotting extension to be activated. Run `using Makie`, `using CairoMakie`, `using GLMakie` or any other package that also loads Makie.")
#            end
#        end
#    end
#end

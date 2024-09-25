
struct PlotNew{A,B,C} <: AbstractType
    mesh::A
    face_data::B
    node_data::C
end

face_data(plt::PlotNew,d) = plt.face_data[d+1]

function face_data(plt::PlotNew;merge_dims=Val(false))
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

node_data(plt::PlotNew) = plt.node_data

function plot(mesh::AbstractMesh)
    PlotNew(mesh,face_data(mesh),node_data(mesh))
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
    dict
end

function node_data(mesh::AbstractMesh)
    nnodes = num_nodes(mesh)
    dict = Dict{String,Vector{Int32}}()
    for group in physical_nodes(mesh;merge_dims=Val(true))
        name,nodes = group
        node_mask = zeros(Int32,nnodes)
        node_mask[nodes] .= 1
        dict[name] = node_mask
    end
    dict
end

function face_data(mesh::AbstractMesh)
    D = num_dims(mesh)
    map(d->face_data(mesh,d),0:D)
end

function restrict_to_dim(plt::PlotNew,d)
    mesh = restrict_to_dim(plt.mesh,d)
    dfacedata = face_data(plt,d)
    fd = map(0:d) do i
        if i == d
            dfacedata
        else
            typeof(dfacedata)()
        end
    end
    pltd = PlotNew(mesh,fd,node_data(plt))
    pltd
end

function restrict(plt::PlotNew,newnodes,newfaces)
    mesh = restrict(plt.mesh,newnodes,newfaces)
    nodedata = node_data(plt)[newnodes]
    facedata = map(face_data(plt),newfaces) do data,nfs
        dict = copy(data)
        for (k,v) in dict
            dict[k] = dict[k][nfs]
        end
        dict
    end
    PlotNew(mesh,facedata,nodedata)
end

function shrink(plt::PlotNew;scale=0.75)
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
    face_newnodes = map(copy,face_nodes(mesh))
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
    new_mesh = GT.mesh_from_arrays(
            newnode_to_x,
            face_newnodes,
            face_reference_id(mesh),
            reference_faces(mesh);
            physical_faces = physical_faces(mesh), # TODO propagate also periodic nodes??
            outwards_normals = outwards_normals(mesh)
           )
    newnode_data = Dict{String,Any}()
    for (k,node_to_val) in node_data(plt)
        newnode_data[k] = node_to_val[newnode_to_node]
    end
    PlotNew(new_mesh,plt.face_data,newnode_data)
end

# VTK

function save_vtk(f,filename,plt::PlotNew)
    function translate(v)
        v
    end
    function translate(v::AbstractVector{<:SVector{2}})
        z = zero(eltype(eltype(v)))
        map(vi->SVector((vi...,z)),v)
    end
    mesh = plt.mesh
    D = num_dims(mesh)
    vtk_grid(filename,GT.vtk_args(mesh)...) do vtk
        for (k,v) in node_data(plt)
            vtk[k,WriteVTK.VTKPointData()] = translate(v)
        end
        for (k,v) in face_data(plt;merge_dims=true)
            vtk[k,WriteVTK.VTKCellData()] = translate(v)
        end
        f(vtk)
    end
end

function save_vtk(filename,plt::PlotNew)
    save_vtk(identity,filename,plt)
end

function save_vtk(f,filename,mesh::AbstractMesh)
    plt = plot(mesh)
    save_vtk(f,filename,plt)
end

function save_vtk(filename,mesh::AbstractMesh)
    save_vtk(identity,filename,mesh)
end

# Makie


Makie.@recipe(Makie3d) do scene
    dt = Makie.default_theme(scene, Makie.Mesh)
    dt
end

Makie.preferred_axis_type(plot::Makie3d) = Makie.LScene

function Makie.plot!(sc::Makie3d{<:Tuple{<:PlotNew}})
    plt = sc[1]
    plt2 = Makie.lift(makie_volumes_impl,plt)
    valid_attributes = Makie.shared_attributes(sc, Makie2d)
    makie2d!(sc,valid_attributes,plt2)
end

function makie_volumes_impl(plt::PlotNew)
    D=3
    d=2
    mesh, = complexify(restrict_to_dim(plt.mesh,D))
    topo = topology(mesh)
    face_to_cells = face_incidence(topo,d,D)
    face_isboundary = map(cells->length(cells)==1,face_to_cells)
    mesh2 = restrict_to_dim(mesh,d)
    newnodes = 1:num_nodes(mesh)
    newfaces = [ Int[] for _ in 0:d ]
    newfaces[end] = findall(face_isboundary)
    mesh3 = restrict(mesh,newnodes,newfaces)
    face_to_cell = map(first,face_to_cells)
    newface_to_cell = face_to_cell[newfaces[end]]
    celldata = copy(face_data(plt,D))
    for (k,v) in celldata
        celldata[k] = v[newface_to_cell]
    end
    fd = map(0:d) do i
        if i == d
            celldata
        else
            typeof(celldata)()
        end
    end
    plt3 = PlotNew(mesh3,fd,node_data(plt))
    plt3
end

Makie.@recipe(Makie3d1d) do scene
    dt = Makie.default_theme(scene, Makie2d1d)
    dt
end

Makie.preferred_axis_type(plot::Makie3d1d) = Makie.LScene

function Makie.plot!(sc::Makie3d1d{<:Tuple{<:PlotNew}})
    plt = sc[1]
    plt2 = Makie.lift(makie_volumes_impl,plt)
    valid_attributes = Makie.shared_attributes(sc, Makie2d1d)
    makie2d1d!(sc,valid_attributes,plt2)
end

Makie.@recipe(Makie2d) do scene
    dt = Makie.default_theme(scene, Makie.Mesh)
    dt
end

function Makie.plot!(sc::Makie2d{<:Tuple{<:PlotNew}})
    plt = sc[1]
    args = Makie.lift(makie_faces_impl,plt)
    vert = Makie.lift(a->a.vert,args)
    conn = Makie.lift(a->a.conn,args)
    valid_attributes = Makie.shared_attributes(sc, Makie.Mesh)
    Makie.mesh!(sc,valid_attributes,vert,conn)
end

Makie.preferred_axis_type(plot::Makie2d) = Makie.LScene

function makie_faces_impl(plt)
    mesh = plt.mesh
    D = num_ambient_dims(mesh)
    nnodes = num_nodes(mesh)
    vert = zeros(Float64,nnodes,D)
    node_to_x = node_coordinates(mesh)
    for (node,x) in enumerate(node_to_x)
        vert[node,:] = x
    end
    d = 2
    nfaces = num_faces(mesh,d)
    face_to_nodes = face_nodes(mesh,d)
    conn = zeros(Int32,nfaces,3)
    for (face,nodes) in enumerate(face_to_nodes)
        conn[face,:] = nodes
    end
    (;vert,conn)
end

Makie.@recipe(Makie2d1d) do scene
    dt = Makie.default_theme(scene, Makie1d)
    dt
end

Makie.preferred_axis_type(plot::Makie2d1d) = Makie.LScene

function Makie.plot!(sc::Makie2d1d{<:Tuple{<:PlotNew}})
    plt = sc[1]
    plt2 = Makie.lift(makie_face_edges_impl,plt)
    valid_attributes = Makie.shared_attributes(sc, Makie1d)
    makie1d!(sc,valid_attributes,plt2)
end

function makie_face_edges_impl(plt)
    # TODO maybe already complexified
    D=2
    d=1
    mesh2, = complexify(restrict_to_dim(plt.mesh,D))
    topo = topology(mesh2)
    edge_to_faces = face_incidence(topo,d,D)
    facedata = copy(face_data(plt,D))
    for (k,v) in facedata
        facedata[k] = map(faces -> sum(v[faces])/length(faces) ,edge_to_faces)
    end
    mesh3 = restrict_to_dim(mesh2,d)
    fd = map(0:d) do i
        if i == d
            facedata
        else
            typeof(facedata)()
        end
    end
    plt3 = PlotNew(mesh3,fd,node_data(plt))
    plt3
end

Makie.@recipe(Makie1d) do scene
    dt = Makie.default_theme(scene, Makie.LineSegments)
    dt
end

Makie.preferred_axis_type(plot::Makie1d) = Makie.LScene

function Makie.plot!(sc::Makie1d{<:Tuple{<:PlotNew}})
    plt = sc[1]
    x = Makie.lift(makie_edges_impl,plt)
    valid_attributes = Makie.shared_attributes(sc, Makie.LineSegments)
    Makie.linesegments!(sc,valid_attributes,x)
end

function makie_edges_impl(plt)
    mesh = plt.mesh
    d = 1
    nedges = num_faces(mesh,d)
    node_to_x = node_coordinates(mesh)
    edge_to_nodes = face_nodes(mesh,d)
    T = eltype(eltype(node_to_x))
    S = num_ambient_dims(mesh)
    V = eltype(node_to_x)
    p = Vector{V}(undef,2*nedges)
    k = 0
    for edge in 1:nedges
        nodes = edge_to_nodes[edge]
        for node in nodes
            k += 1
            p[k] = node_to_x[node]
        end
    end
    p
end

Makie.@recipe(Makie0d) do scene
    dt = Makie.default_theme(scene, Makie.Scatter)
    dt
end

Makie.preferred_axis_type(plot::Makie0d) = Makie.LScene

function Makie.plot!(sc::Makie0d{<:Tuple{<:PlotNew}})
    plt = sc[1]
    x = Makie.lift(makie_vertices_impl,plt)
    valid_attributes = Makie.shared_attributes(sc, Makie.Scatter)
    Makie.scatter!(sc,valid_attributes,x)
end

function makie_vertices_impl(plt)
    mesh = plt.mesh
    d = 0
    nedges = num_faces(mesh,d)
    node_to_x = node_coordinates(mesh)
    edge_to_nodes = face_nodes(mesh,d)
    T = eltype(eltype(node_to_x))
    S = num_ambient_dims(mesh)
    V = eltype(node_to_x)
    p = Vector{V}(undef,nedges)
    k = 0
    for edge in 1:nedges
        nodes = edge_to_nodes[edge]
        for node in nodes
            k += 1
            p[k] = node_to_x[node]
        end
    end
    p
end



#function makie_mesh_args(plt)
#    mesh = plt.mesh
#    # assumes that mesh is a "visualizable" mesh
#    # i.e., all 2d objects are triangles.
#    # Only 2d objects (triangles) are handled
#    D = num_ambient_dims(mesh)
#    nnodes = num_nodes(mesh)
#    vert = zeros(Float64,nnodes,D)
#    node_to_x = node_coordinates(mesh)
#    for (node,x) in enumerate(node_to_x)
#        vert[node,:] = x
#    end
#    d = 2
#    nfaces = num_faces(mesh,d)
#    face_to_nodes = face_nodes(mesh,d)
#    conn = zeros(Int32,nfaces,3)
#    for (face,nodes) in enumerate(face_to_nodes)
#        conn[face,:] = nodes
#    end
#    (vert,conn)
#end
#
#function makie_linesegments_args(plt)
#    mesh = plt.mesh
#    d = 1
#    nedges = num_faces(mesh,d)
#    node_to_x = node_coordinates(mesh)
#    edge_to_nodes = face_nodes(mesh,d)
#    T = eltype(eltype(node_to_x))
#    S = num_ambient_dims(mesh)
#    V = eltype(node_to_x)
#    p = Vector{V}(undef,2*nedges)
#    k = 0
#    for edge in 1:nedges
#        nodes = edge_to_nodes[edge]
#        for node in nodes
#            k += 1
#            p[k] = node_to_x[node]
#        end
#    end
#    (p,)
#end
#
#function makie_scatter_args(plt)
#    mesh = plt.mesh
#    d = 0
#    nedges = num_faces(mesh,d)
#    node_to_x = node_coordinates(mesh)
#    edge_to_nodes = face_nodes(mesh,d)
#    T = eltype(eltype(node_to_x))
#    S = num_ambient_dims(mesh)
#    V = eltype(node_to_x)
#    p = Vector{V}(undef,nedges)
#    k = 0
#    for edge in 1:nedges
#        nodes = edge_to_nodes[edge]
#        for node in nodes
#            k += 1
#            p[k] = node_to_x[node]
#        end
#    end
#    (p,)
#end
#
#struct FaceColor
#    name::String
#end
#
#struct NodeColor
#    name::String
#end
#
#function makie_mesh_color(plt,color)
#    (;color)
#end
#
#function makie_linesegments_color(plt,color)
#    (;color)
#end
#
#function makie_scatter_color(plt,color)
#    (;color)
#end
#
#function makie_mesh_color(plt,::Nothing)
#    (;)
#end
#
#function makie_linesegments_color(plt,::Nothing)
#    (;)
#end
#
#function makie_scatter_color(plt,::Nothing)
#    (;)
#end
#
#function makie_mesh_color(plt,fc::FaceColor)
#    mesh = plt.mesh
#    allcolors = face_data(plt;merge_dims=true)[fc.name]
#    minc = minimum(allcolors)
#    maxc = maximum(allcolors)
#    colorrange = (minc,maxc)
#    d = 2
#    offset = face_offset(mesh,d)
#    npoints = num_faces(mesh,d)
#    points = 1:npoints
#    dface_to_color = allcolors[points .+ offset]
#    node_to_color = similar(dface_to_color,num_nodes(mesh))
#    dface_to_nodes = face_nodes(mesh,d)
#    for face in 1:npoints
#        nodes = dface_to_nodes[face]
#        node_to_color[nodes] .= dface_to_color[face]
#    end
#    color = node_to_color
#    (;color,colorrange)
#end
#
#function makie_linesegments_color(plt,fc::FaceColor)
#    mesh = plt.mesh
#    allcolors = face_data(plt;merge_dims=true)[fc.name]
#    minc = minimum(allcolors)
#    maxc = maximum(allcolors)
#    colorrange = (minc,maxc)
#    d = 1
#    offset = face_offset(mesh,d)
#    npoints = num_faces(mesh,d)
#    points = 1:npoints
#    dface_to_color = allcolors[points .+ offset]
#    color = similar(dface_to_color,2*npoints)
#    color[1:2:(end-1)] = dface_to_color
#    color[2:2:end] = dface_to_color
#    (;color,colorrange)
#end
#
#function makie_scatter_color(plt,fc::FaceColor)
#    mesh = plt.mesh
#    allcolors = face_data(plt;merge_dims=true)[fc.name]
#    minc = minimum(allcolors)
#    maxc = maximum(allcolors)
#    colorrange = (minc,maxc)
#    d = 0
#    offset = face_offset(mesh,d)
#    npoints = num_faces(mesh,d)
#    points = 1:npoints
#    color = allcolors[points .+ offset]
#    (;color,colorrange)
#end
#
#function render_with_makie!(ax,plt::PlotNew;color=nothing)
#    plt = shrink(plt;scale=1)
#    args1 = makie_mesh_args(plt)
#    args2 = makie_linesegments_args(plt)
#    args3 = makie_scatter_args(plt)
#    color1 = makie_mesh_color(plt,color)
#    color2 = makie_linesegments_color(plt,color)
#    color3 = makie_scatter_color(plt,color)
#    Makie.mesh!(ax,args1...;color1...)
#    Makie.linesegments!(ax,args2...;color2...)
#    Makie.scatter!(ax,args3...;color3...)
#    ax
#end
#
#Makie.@recipe(MakiePlot) do scene
#    dt = Makie.default_theme(scene, Makie.Mesh)
#    Makie.Attributes(
#        color=dt[:color],
#
#      )
#end
#
#function makie_mesh_color_2(plt,fc)
#    fc
#end
#
#function makie_mesh_color_2(plt,fc::FaceColor)
#    mesh = plt.mesh
#    allcolors = face_data(plt;merge_dims=true)[fc.name]
#    minc = minimum(allcolors)
#    maxc = maximum(allcolors)
#    colorrange = (minc,maxc)
#    d = 2
#    offset = face_offset(mesh,d)
#    npoints = num_faces(mesh,d)
#    points = 1:npoints
#    dface_to_color = allcolors[points .+ offset]
#    node_to_color = similar(dface_to_color,num_nodes(mesh))
#    dface_to_nodes = face_nodes(mesh,d)
#    for face in 1:npoints
#        nodes = dface_to_nodes[face]
#        node_to_color[nodes] .= dface_to_color[face]
#    end
#    color = node_to_color
#    color
#end
#
#function Makie.plot!(sc::MakiePlot{<:Tuple{<:PlotNew}})
#    plt = sc[1]
#    args1 = Makie.lift(makie_mesh_args,plt)
#    vert = Makie.lift(first,args1)
#    conn = Makie.lift(last,args1)
#    valid_attributes = Makie.shared_attributes(sc, Makie.Mesh)
#    valid_attributes[:color] = Makie.lift(makie_mesh_color_2,plt,valid_attributes[:color])
#    Makie.mesh!(sc,valid_attributes,vert,conn)
#end




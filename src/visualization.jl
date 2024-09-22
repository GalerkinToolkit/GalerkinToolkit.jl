
struct PlotNew{A,B,C,D}
    mesh::A
    face_data::B
    node_data::C
    params::D
end

face_data(a::PlotNew,d) = a.face_data[d+1]

function face_data(a::PlotNew;merge_dims=true)
    # TODO if all dims are empty except one
    # do not create new data
    if merge_dims
        mesh = a.mesh
        nfaces = sum(num_faces(mesh))
        offsets = face_offset(mesh)
        D = num_dims(mesh)
        dict = Dict{String,Any}()
        for d in 0:D
            offset = offsets[d+1]
            for group in a.face_data[d+1]
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
        a.face_data
    end
end

node_data(a::PlotNew) = a.node_data

function plot(mesh;params...)
    D = num_dims(mesh)
    face_data = [ Dict{String,Any}() for _ in 0:D ]
    node_data = Dict{String,Any}()
    PlotNew(mesh,face_data,node_data,params)
end

function plot(mesh,face_data,node_data;params...)
    PlotNew(mesh,face_data,node_data,params)
end

function mesh_of_dim(mesh,d)
    error("todo")
end

function plot(mesh,d;params...)
    dmesh = mesh_of_dim(mesh,d)
    plot(dmesh;params...)
end

function plot!(plt::PlotNew,::typeof(physical_faces),d)
    mesh = plt.mesh
    ndfaces = num_faces(mesh,d)
    dict = face_data(plt,d)
    for group in physical_faces(mesh,d)
        name,faces = group
        face_mask = zeros(Int,ndfaces)
        face_mask[faces] .= 1
        dict[name] = face_mask
    end
    plt
end

function plot!(plt::PlotNew,::typeof(physical_faces))
    mesh = plt.mesh
    D = num_dims(mesh)
    for d in 0:D
        plot!(plt,physical_faces,d)
    end
    plt
end

function shrink(plt;coeff=0.75)
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
                newnode_to_x[newnode] = coeff*(x-xm) + xm
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
    plot(new_mesh,plt.face_data,newnode_data)
end


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
    plot!(plt,physical_faces)
    save_vtk(f,filename,plt)
end

function save_vtk(filename,mesh::AbstractMesh)
    save_vtk(identity,filename,mesh)
end

function makie_mesh_args(plt)
    mesh = plt.mesh
    # assumes that mesh is a "visualizable" mesh
    # i.e., all 2d objects are triangles.
    # Only 2d objects (triangles) are handled
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
    (vert,conn)
end

function makie_linesegments_args(plt)
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
    (p,)
end

function makie_scatter_args(plt)
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
    (p,)
end

struct FaceColor
    name::String
end

struct NodeColor
    name::String
end

function makie_mesh_color(plt,color)
    (;color)
end

function makie_linesegments_color(plt,color)
    (;color)
end

function makie_scatter_color(plt,color)
    (;color)
end

function makie_mesh_color(plt,::Nothing)
    (;)
end

function makie_linesegments_color(plt,::Nothing)
    (;)
end

function makie_scatter_color(plt,::Nothing)
    (;)
end

function makie_mesh_color(plt,fc::FaceColor)
    mesh = plt.mesh
    allcolors = face_data(plt;merge_dims=true)[fc.name]
    minc = minimum(allcolors)
    maxc = maximum(allcolors)
    colorrange = (minc,maxc)
    d = 2
    offset = face_offset(mesh,d)
    npoints = num_faces(mesh,d)
    points = 1:npoints
    dface_to_color = allcolors[points .+ offset]
    node_to_color = similar(dface_to_color,num_nodes(mesh))
    dface_to_nodes = face_nodes(mesh,d)
    for face in 1:npoints
        nodes = dface_to_nodes[face]
        node_to_color[nodes] .= dface_to_color[face]
    end
    color = node_to_color
    (;color,colorrange)
end

function makie_linesegments_color(plt,fc::FaceColor)
    mesh = plt.mesh
    allcolors = face_data(plt;merge_dims=true)[fc.name]
    minc = minimum(allcolors)
    maxc = maximum(allcolors)
    colorrange = (minc,maxc)
    d = 1
    offset = face_offset(mesh,d)
    npoints = num_faces(mesh,d)
    points = 1:npoints
    dface_to_color = allcolors[points .+ offset]
    color = similar(dface_to_color,2*npoints)
    color[1:2:(end-1)] = dface_to_color
    color[2:2:end] = dface_to_color
    (;color,colorrange)
end

function makie_scatter_color(plt,fc::FaceColor)
    mesh = plt.mesh
    allcolors = face_data(plt;merge_dims=true)[fc.name]
    minc = minimum(allcolors)
    maxc = maximum(allcolors)
    colorrange = (minc,maxc)
    d = 0
    offset = face_offset(mesh,d)
    npoints = num_faces(mesh,d)
    points = 1:npoints
    color = allcolors[points .+ offset]
    (;color,colorrange)
end

function render_with_makie!(ax,plt::PlotNew;color=nothing)
    plt = shrink(plt;coeff=1)
    args1 = makie_mesh_args(plt)
    args2 = makie_linesegments_args(plt)
    args3 = makie_scatter_args(plt)
    color1 = makie_mesh_color(plt,color)
    color2 = makie_linesegments_color(plt,color)
    color3 = makie_scatter_color(plt,color)
    MakieCore.mesh!(ax,args1...;color1...)
    MakieCore.linesegments!(ax,args2...;color2...)
    MakieCore.scatter!(ax,args3...;color3...)
    ax
end





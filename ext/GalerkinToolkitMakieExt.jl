module GalerkinToolkitMakieExt

#TODOS
#shadows
#linewidth
#arrows
#pplot,pmesh,pdomain

import GalerkinToolkit as GT
import GalerkinToolkit: makieplot, makieplot!
import GalerkinToolkit: makie0d, makie0d!, makie1d, makie1d!
import GalerkinToolkit: makie2d, makie2d!, makie2d1d, makie2d1d!
import GalerkinToolkit: makie3d, makie3d!, makie3d1d, makie3d1d!
using PartitionedArrays
using Makie

Makie.@recipe Makie0d begin
    Makie.documented_attributes(Makie.Scatter)...
end

function Makie.preferred_axis_type(plot::Makie0d)
    d = GT.num_ambient_dims(plot[1][])
    d == 3 ? Makie.LScene : Makie.Axis
end

function Makie.plot!(p::Makie0d{<:Tuple{<:GT.Plot}})
    valid_attributes = Makie.shared_attributes(p, Makie.Scatter)
    attrs_in = [:converted_1,:color,:colorrange]
    attrs_out = :newcolorrange
    map!(setup_colorrange_impl,p.attributes,attrs_in,attrs_out)
    attrs_in = [:converted_1,:color]
    attrs_out = [:points,:newcolor]
    map!(makie_vertices_impl,p.attributes,attrs_in,attrs_out)
    color = p.newcolor
    colorrange = p.newcolorrange
    points = p.points
    Makie.scatter!(p,valid_attributes,points;color,colorrange)
end

Makie.@recipe Makie1d begin
    Makie.documented_attributes(Makie.LineSegments)...
end

function Makie.preferred_axis_type(plot::Makie1d)
    d = GT.num_ambient_dims(plot[1][])
    d == 3 ? Makie.LScene : Makie.Axis
end

function Makie.plot!(p::Makie1d{<:Tuple{<:GT.Plot}})
    valid_attributes = Makie.shared_attributes(p, Makie.LineSegments)
    attrs_in = [:converted_1,:color,:colorrange]
    attrs_out = :newcolorrange
    map!(setup_colorrange_impl,p.attributes,attrs_in,attrs_out)
    attrs_in = [:converted_1,:color]
    attrs_out = [:points,:newcolor]
    map!(makie_edges_impl,p.attributes,attrs_in,attrs_out)
    color = p.newcolor
    colorrange = p.newcolorrange
    points = p.points
    Makie.linesegments!(p,valid_attributes,points;color,colorrange)
end

Makie.@recipe Makie2d begin
    Makie.documented_attributes(Makie.Mesh)...
end

function Makie.preferred_axis_type(plot::Makie2d)
    d = GT.num_ambient_dims(plot[1][])
    d == 3 ? Makie.LScene : Makie.Axis
end

function Makie.plot!(p::Makie2d{<:Tuple{<:GT.Plot}})
    valid_attributes = Makie.shared_attributes(p, Makie.Mesh)
    attrs_in = [:converted_1,:color,:colorrange]
    attrs_out = :newcolorrange
    map!(setup_colorrange_impl,p.attributes,attrs_in,attrs_out)
    attrs_in = [:converted_1,:color]
    attrs_out = [:vert,:conn,:newcolor]
    map!(makie_faces_impl,p.attributes,attrs_in,attrs_out)
    color = p.newcolor
    colorrange = p.newcolorrange
    Makie.mesh!(p,valid_attributes,p.vert,p.conn;color,colorrange)
    p
end

Makie.@recipe Makie3d begin
    Makie.documented_attributes(Makie2d)...
end

Makie.preferred_axis_type(plot::Makie3d) = Makie.LScene

function Makie.plot!(p::Makie3d{<:Tuple{<:GT.Plot}})
    valid_attributes = Makie.shared_attributes(p, Makie2d)
    attrs_in = [:converted_1,]
    attrs_out = :plt
    map!(makie_volumes_impl,p.attributes,attrs_in,attrs_out)
    plt = p.plt
    makie2d!(p,valid_attributes,plt)
end

Makie.@recipe Makie2d1d begin
    Makie.documented_attributes(Makie1d)...
end

function Makie.preferred_axis_type(plot::Makie2d1d)
    d = GT.num_ambient_dims(plot[1][])
    d == 3 ? Makie.LScene : Makie.Axis
end

function Makie.plot!(p::Makie2d1d{<:Tuple{<:GT.Plot}})
    valid_attributes = Makie.shared_attributes(p, Makie1d)
    attrs_in = [:converted_1,]
    attrs_out = :plt
    map!(makie_face_edges_impl,p.attributes,attrs_in,attrs_out)
    plt = p.plt
    makie1d!(p,valid_attributes,plt)
end

Makie.@recipe Makie3d1d begin
    Makie.documented_attributes(Makie2d1d)...
end

Makie.preferred_axis_type(plot::Makie3d1d) = Makie.LScene

function Makie.plot!(p::Makie3d1d{<:Tuple{<:GT.Plot}})
    attrs_in = [:converted_1,]
    attrs_out = :plt
    map!(makie_volumes_impl,p.attributes,attrs_in,attrs_out)
    plt = p.plt
    valid_attributes = Makie.shared_attributes(p, Makie2d)
    makie2d1d!(p,valid_attributes,plt)
end

Makie.@recipe MakiePlot begin
    Makie.documented_attributes(Makie.Mesh)...
    dim = Makie.Automatic()
    shrink = false
    warp_by_vector = nothing
    warp_by_scalar = nothing
    warp_scale = 1
    strokecolor = nothing
    strokewidth = nothing
    color      = :lightblue
    colormap   = :bluesreds
    colorrange = Makie.Automatic()
    shading    = Makie.NoShading
    cycle      = nothing
    refinement = nothing
end

function Makie.preferred_axis_type(plot::MakiePlot)
    plt = plot[1][]
    d = GT.num_ambient_dims(plt)
    d == 3 ? Makie.LScene : Makie.Axis
end

function Makie.plot!(p::MakiePlot{<:Tuple{<:GT.Plot}})
    attrs_in = [:converted_1,:warp_by_vector,:warp_scale]
    attrs_out = :plt_warp_vector
    map!(p.attributes,attrs_in,attrs_out) do plt,vector,scale
        GT.warp_by_vector(plt,vector;scale)
    end
    attrs_in = [:plt_warp_vector,:warp_by_scalar,:warp_scale]
    attrs_out = :plt_warp_scalar
    map!(p.attributes,attrs_in,attrs_out) do plt,scalar,scale
        GT.warp_by_scalar(plt,scalar;scale)
    end
    attrs_in = [:plt_warp_scalar,:shrink]
    attrs_out = :plt_shrink
    map!(p.attributes,attrs_in,attrs_out) do plt,scale
        GT.shrink(plt;scale)
    end
    plt = p.plt_shrink
    plt_now = plt[]
    dim_now = p.dim[]
    d_now = GT.num_dims(plt_now.mesh)
    dim_now = dim_now == Makie.Automatic() ? d_now : dim_now
    if 0 in dim_now 
        if p.color[] !== nothing
            valid_attributes = Makie.shared_attributes(p,Makie0d)
            makie0d!(p,valid_attributes,plt)
        end
    end
    if 1 in dim_now
        if p.color[] !== nothing
            valid_attributes = Makie.shared_attributes(p,Makie1d)
            makie1d!(p,valid_attributes,plt)
        end
    end
    if 2 in dim_now
        if p.color[] !== nothing
            valid_attributes = Makie.shared_attributes(p,Makie2d)
            makie2d!(p,valid_attributes,plt)
        end
        if p.strokecolor[] !== nothing
            color = p.strokecolor
            linewidth = p.strokewidth
            valid_attributes = Makie.shared_attributes(p,Makie2d1d)
            makie2d1d!(p,valid_attributes,plt;color)#,linewidth)
        end
    end
    if 3 in dim_now
        if p.color[] !== nothing
            valid_attributes = Makie.shared_attributes(p,Makie3d)
            makie3d!(p,valid_attributes,plt)
        end
        if p.strokecolor[] !== nothing
            color = p.strokecolor
            linewidth = p.strokewidth
            valid_attributes = Makie.shared_attributes(p,Makie3d1d)
            makie3d1d!(p,valid_attributes,plt;color)#,linewidth)
        end
    end
    p
end

Makie.plottype(::GT.Plot) = MakiePlot
Makie.plottype(::GT.PPlot) = MakiePlot
Makie.plottype(::GT.AbstractMesh) = MakiePlot
Makie.plottype(::GT.PMesh) = MakiePlot

function Makie.convert_arguments(::Type{<:MakiePlot},mesh::GT.AbstractMesh)
    plt = GT.plot(mesh)
    (plt,)
end

function Makie.convert_arguments(::Type{<:MakiePlot},mesh::GT.PMesh)
    plt = GT.plot(mesh)
    (plt,)
end

Makie.plottype(::GT.AbstractDomain) = MakiePlot

function Makie.plot!(p::MakiePlot{<:Tuple{<:GT.AbstractDomain}})
    attrs_in = [:converted_1,:refinement,:color,:warp_by_vector,:warp_by_scalar]
    attrs_out = [:plt,:newcolor,:new_warp_by_vector,:new_warp_by_scalar]
    map!(p.attributes,attrs_in,attrs_out) do domain,refinement,color,vec,scal
        plt = GT.plot(domain;refinement)
        if isa(color,GT.AbstractQuantity) || isa(color,Function)
            label = string(gensym())
            GT.plot!(plt,color;label)
            newcolor = GT.node_color(plt,label)
        else
            newcolor = color
        end
        if isa(vec,GT.AbstractQuantity) || isa(vec,Function)
            label = string(gensym())
            GT.plot!(plt,vec;label)
            newvec = GT.NodeData(label)
        else
            newvec = vec
        end
        if isa(scal,GT.AbstractQuantity) || isa(scal,Function)
            label = string(gensym())
            GT.plot!(plt,scal;label)
            newscal = GT.NodeData(label)
        else
            newscal = scal
        end
        (plt,newcolor,newvec,newscal)
    end
    color = p.newcolor
    warp_by_vector = p.new_warp_by_vector
    warp_by_scalar = p.new_warp_by_scalar
    valid_attributes = Makie.shared_attributes(p,MakiePlot)
    plot!(p,valid_attributes,p.plt;color,warp_by_vector,warp_by_scalar)
    p
end

# TODO not sure about this
# what if u is not scalar-valued?
#Makie.plottype(::GT.AbstractQuantity) = MakiePlot
#
#function Makie.plot!(sc::MakiePlot{<:Tuple{<:GT.AbstractQuantity}})
#    u = sc[1]
#    valid_attributes = Makie.shared_attributes(sc, MakiePlot)
#    dom = Makie.lift(domain,u)
#    valid_attributes[:color] = u
#    makieplot!(sc,valid_attributes,dom)
#end

function vector_of_observables(a)
    #TODO remove trick for DebugArray
    function start(v)
        first(v)
    end
    function rest(v)
        v[2:end]
    end
    if length(a[]) == 1
        b = Makie.lift(first,a)
        return [b,]
    else
        b = Makie.lift(start,a)
        c = Makie.lift(rest,a)
        return [b,vector_of_observables(c)...]
    end
end

function makie_volumes_impl(plt::GT.Plot;simplexify=Val(false))
    @assert GT.num_dims(plt.mesh) == 3
    D=3
    d=2
    mesh = GT.complexify(GT.restrict_to_dim(plt.mesh,D))
    topo = GT.topology(mesh)
    face_to_cells = GT.face_incidence(topo,d,D)
    face_isboundary = map(cells->length(cells)==1,face_to_cells)
    mesh2 = GT.restrict_to_dim(mesh,d)
    newnodes = 1:GT.num_nodes(mesh2)
    newfaces = [ Int[] for _ in 0:d ]
    newfaces[end] = findall(face_isboundary)
    mesh3 = GT.restrict(mesh2,newnodes,newfaces)
    face_to_cell = map(first,face_to_cells)
    newface_to_cell = face_to_cell[newfaces[end]]
    celldata = copy(GT.face_data(plt,D;merge_dims=true))
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
    plt3 = GT.Plot(mesh3,fd,GT.node_data(plt))
    if GT.val_parameter(simplexify)
        GT.simplexify(plt3)
    else
        plt3
    end
end

function setup_colorrange_impl(plt,color,colorrange)
    if colorrange != Makie.Automatic()
        return colorrange
    end
    if isa(color,GT.FaceData)
        colorrange = GT.face_colorrange(plt,color.name)
    end
    if isa(color,GT.NodeData)
        colorrange = GT.node_colorrange(plt,color.name)
    end
    colorrange
end

function setup_colors_impl(plt,color,d)
    if isa(color,GT.FaceData)
        if d == 2
            plt = GT.shrink(plt;scale=1)
        end
        color2 = GT.face_color(plt,color.name,d)
    end
    if isa(color,GT.NodeData)
        color2 = GT.node_color(plt,color.name)
    end
    if isa(color,GT.NodeData) || isa(color,GT.FaceData)
        face_to_nodes = GT.face_nodes(plt.mesh,d)
        nfaces = length(face_to_nodes)
        if d == 0
            color = similar(color2,nfaces)
            for (face,nodes) in enumerate(face_to_nodes)
                node = first(nodes)
                color[face] = color2[node]
            end
        elseif d == 1
            color = similar(color2,2*nfaces)
            k = 0
            for nodes in face_to_nodes
                for node in nodes
                    k+=1
                    color[k] = color2[node]
                end
            end
        elseif d == 2
            color = color2
        else
            error()
        end
    end
    plt,color
end

function makie_faces_impl(plt,color)
    plt = GT.simplexify(plt)
    d = 2
    #plt = shrink(plt,scale=0.995)
    plt,color = setup_colors_impl(plt,color,d)
    mesh = plt.mesh
    D = GT.num_ambient_dims(mesh)
    nnodes = GT.num_nodes(mesh)
    vert = zeros(Float64,nnodes,D)
    node_to_x = GT.node_coordinates(mesh)
    for (node,x) in enumerate(node_to_x)
        vert[node,:] = x
    end
    nfaces = GT.num_faces(mesh,d)
    face_to_nodes = GT.face_nodes(mesh,d)
    conn = zeros(Int32,nfaces,3)
    for (face,nodes) in enumerate(face_to_nodes)
        conn[face,:] = nodes
    end
    if nfaces == 0
        vert = ones(Float64,3,D)
        conn = zeros(Int32,1,3)
        conn[:] = 1:3
        color =:pink
    end
    (vert,conn,color)
end

function makie_face_edges_impl(plt)
    # TODO maybe already complexified
    D=2
    d=1
    mesh2 = GT.complexify(GT.restrict_to_dim(plt.mesh,D))
    topo = GT.topology(mesh2)
    edge_to_faces = GT.face_incidence(topo,d,D)
    K = Any
    facedata = Dict{String,K}()
    for (k,v) in GT.face_data(plt,D;merge_dims=true)
        facedata[k] = map(faces -> sum(v[faces])/length(faces) ,edge_to_faces)
    end
    mesh3 = GT.restrict_to_dim(mesh2,d)
    fd = map(0:d) do i
        if i == d
            facedata
        else
            typeof(facedata)()
        end
    end
    plt3 = GT.Plot(mesh3,fd,GT.node_data(plt))
    plt3
end

function makie_edges_impl(plt,color)
    d = 1
    plt,color = setup_colors_impl(plt,color,d)
    mesh = plt.mesh
    nedges = GT.num_faces(mesh,d)
    node_to_x = GT.node_coordinates(mesh)
    edge_to_nodes = GT.face_nodes(mesh,d)
    T = eltype(eltype(node_to_x))
    S = GT.num_ambient_dims(mesh)
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
    (p,color)
end

function makie_vertices_impl(plt,color)
    d = 0
    plt,color = setup_colors_impl(plt,color,d)
    mesh = plt.mesh
    nedges = GT.num_faces(mesh,d)
    node_to_x = GT.node_coordinates(mesh)
    edge_to_nodes = GT.face_nodes(mesh,d)
    T = eltype(eltype(node_to_x))
    S = GT.num_ambient_dims(mesh)
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
    (p,color)
end

# Vector-valued

function Makie.plot!(p::Arrows{<:Tuple{<:GT.Plot,GT.NodeData}})
    valid_attributes = Makie.shared_attributes(p,Arrows)
    attrs_in = [:converted_1,:converted_2,:color]
    attrs_out = [:coords,:vecs,:newcolor]
    map!(makie_arrows_impl,p.attributes,attrs_in,attrs_out)
    color = p.newcolor
    Makie.arrows!(p,valid_attributes,p.coords,p.vecs;color)
end

function makie_arrows_impl(plt,vecs::GT.NodeData,color)
    D = GT.num_ambient_dims(plt.mesh)
    x = GT.node_coordinates(plt.mesh)
    nnodes = length(x)
    if D == 2
        node_to_coord = [ Makie.Point2f(xi) for xi in x]
    elseif D == 3
        node_to_coord = [ Makie.Point3f(xi) for xi in x]
    else
        error("not implemented")
    end
    node_to_vec = plt.node_data[vecs.name]
    if isa(color,GT.NodeData)
        color2 = plt.node_data[color.name]
    else
        color2 = color
    end
    (node_to_coord,node_to_vec,color2)
end

function Makie.plot!(sc::Arrows{<:Tuple{<:GT.AbstractField}})
    valid_attributes = Makie.shared_attributes(p, Arrows)
    attrs_in = [:converted_1,:color]
    attrs_out = [:plt,:vecs,:newcolor]
    args = map!(p.attributes,attrs_in,attrs_out) do q,color
        dom = GT.domain(q)
        plt = GT.plot(dom)
        label_q = string(gensym())
        GT.plot!(plt,q;label=label_q)
        if isa(color,GT.AbstractQuantity) || isa(color,Function)
            label_color = string(gensym())
            GT.plot!(plt,color;label=label_color)
            color = GT.node_color(plt,label_color)
        end
        vecs = GT.NodeData(label_q)
        (plt,vecs,color)
    end
    plt = p.plt
    vecs = p.vecs
    color = p.newcolor
    arrows!(p,valid_attributes,plt,vecs;color)
end

end # modules

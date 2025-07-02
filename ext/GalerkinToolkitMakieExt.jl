module GalerkinToolkitMakieExt
using LinearAlgebra
using StaticArrays

#TODOS
#shadows
#physical groups
#linewidth
#strokecolor also for 2d vertices
# 3d vertices and 3d edges
#arrows
#pplot,pmesh,pdomain

#FaceData -> FaceColor
#NodeData -> NodeColor


import GalerkinToolkit as GT
import GalerkinToolkit: makie0d, makie0d!
import GalerkinToolkit: makie1d, makie1d!
import GalerkinToolkit: makie2d, makie2d!
import GalerkinToolkit: makie3d, makie3d!
using PartitionedArrays
using Makie


#function Makie.convert_arguments(::Type{<:Mesh},plt::GT.Plot,dim::Integer)
#    #color = :blue
#    #colorrange = Makie.Automatic
#    #vert, conn, color, colorrange = setup_makie2d(plt,dim,color,colorrange)
#    #(vert,conn)
#end

function Makie.plot!(p::Makie.Mesh{<:Tuple{<:GT.Plot}})
    yyy
#    attrs_in = [:converted_1]
#    attrs_out = :dim
#    map!(GT.num_dims,p.attributes,attrs_in,attrs_out)
#    valid_attributes = Makie.shared_attributes(p, Makie.Mesh)
#    Makie.mesh!(p,valid_attributes,p[1],p.dim)
end
#
function Makie.plot!(p::Mesh{<:Tuple{<:GT.Plot,<:Integer}})
    yyy
    attrs_in = [:converted_1,:converted_2,:color,:colorrange]
    attrs_out = [:vert,:conn,:newcolor,:newcolorrange]
    map!(setup_makie2d,p.attributes,attrs_in,attrs_out)
    valid_attributes = Makie.shared_attributes(p, Makie.Mesh)
    color = p.newcolor
    colorrange = p.newcolorrange
    Makie.mesh!(p,valid_attributes,p.vert,p.conn;color,colorrange)
end




function shared_attributes_mixin()
    Makie.@DocumentedAttributes begin
        dim = Makie.Automatic()
        shrink = nothing
        warp_by_vector = nothing
        warp_by_scalar = nothing
        warp_scale = 1
        refinement = nothing
    end
end

Makie.@recipe Makie0d begin
    Makie.documented_attributes(Makie.Scatter)...
    shared_attributes_mixin()...
end

function Makie.preferred_axis_type(p::Makie0d)
    plt = p[1][]
    plot_preferred_axis_type(plt)
end

function Makie.plot!(p::Makie0d{<:Tuple{<:GT.Plot}})
    attrs_in = [:converted_1,:warp_by_vector,:warp_by_scalar,:warp_scale,:shrink]
    attrs_out = [:plt]
    map!(setup_plt_changes,p.attributes,attrs_in,attrs_out)
    attrs_in = [:plt,:dim,:color,:colorrange]
    attrs_out = [:points,:newcolor,:newcolorrange]
    map!(setup_makie0d,p.attributes,attrs_in,attrs_out)
    valid_attributes = Makie.shared_attributes(p, Makie.Scatter)
    points = p.points
    color = p.newcolor
    colorrange = p.newcolorrange
    Makie.scatter!(p,valid_attributes,points;color,colorrange)
end

function setup_makie0d(plt_in,dim_in,color,colorrange)
    d = 0
    D = dim_in == Makie.Automatic() ? GT.num_dims(plt_in.mesh) : dim_in
    if D == 3
        plt_in = GT.skin(plt_in)
        dim_in = 2
    else
        plt_in = plt_in
    end
    plt = restrict_to_dim_for_makie(plt_in,dim_in,d)
    colorrange = setup_colorrange_impl(plt,color,colorrange)
    plt, color = setup_colors_impl(plt,color,d)
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
    (p,color,colorrange)
end

Makie.@recipe Makie1d begin
    Makie.documented_attributes(Makie.LineSegments)...
    shared_attributes_mixin()...
end

function Makie.preferred_axis_type(p::Makie1d)
    plt = p[1][]
    plot_preferred_axis_type(plt)
end

function Makie.plot!(p::Makie1d{<:Tuple{<:GT.Plot}})
    attrs_in = [:converted_1,:warp_by_vector,:warp_by_scalar,:warp_scale,:shrink]
    attrs_out = [:plt]
    map!(setup_plt_changes,p.attributes,attrs_in,attrs_out)
    attrs_in = [:plt,:dim,:color,:colorrange]
    attrs_out = [:points,:newcolor,:newcolorrange]
    map!(setup_makie1d,p.attributes,attrs_in,attrs_out)
    valid_attributes = Makie.shared_attributes(p, Makie.LineSegments)
    points = p.points
    color = p.newcolor
    colorrange = p.newcolorrange
    Makie.linesegments!(p,valid_attributes,points;color,colorrange)
end

function setup_makie1d(plt_in,dim_in,color,colorrange)
    d = 1
    D = dim_in == Makie.Automatic() ? GT.num_dims(plt_in.mesh) : dim_in
    if D == 3
        plt_in = GT.skin(plt_in)
        dim_in = 2
    else
        plt_in = plt_in
    end
    plt = restrict_to_dim_for_makie(plt_in,dim_in,d)
    colorrange = setup_colorrange_impl(plt,color,colorrange)
    plt, color = setup_colors_impl(plt,color,d)
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
    (p,color,colorrange)
end

Makie.@recipe Makie2d begin
    Makie.documented_attributes(Makie.Mesh)...
    shared_attributes_mixin()...
end

function Makie.preferred_axis_type(p::Makie2d)
    plt = p[1][]
    plot_preferred_axis_type(plt)
end

function Makie.plot!(p::Makie2d{<:Tuple{<:GT.Plot}})
    attrs_in = [:converted_1,:warp_by_vector,:warp_by_scalar,:warp_scale,:shrink]
    attrs_out = [:plt]
    map!(setup_plt_changes,p.attributes,attrs_in,attrs_out)
    attrs_in = [:plt,:dim,:color,:colorrange]
    attrs_out = [:vert,:conn,:newcolor,:newcolorrange]
    map!(setup_makie2d,p.attributes,attrs_in,attrs_out)
    valid_attributes = Makie.shared_attributes(p, Makie.Mesh)
    color = p.newcolor
    colorrange = p.newcolorrange
    Makie.mesh!(p,valid_attributes,p.vert,p.conn;color,colorrange)
end

function setup_makie2d(plt,dim_in,color_in,colorrange_in)
    d = 2
    D = dim_in == Makie.Automatic() ? GT.num_dims(plt.mesh) : dim_in
    if D == 3
        plt = GT.skin(plt)
    else
        plt = plt
    end
    plt = restrict_to_dim_for_makie(plt,d,d)
    plt = GT.shrink(plt;scale=1)
    plt = GT.simplexify(plt)
    colorrange = setup_colorrange_impl(plt,color_in,colorrange_in)
    plt, color = setup_colors_impl(plt,color_in,d)
    vert,conn = makie_faces_mesh(plt)
    if GT.num_ambient_dims(plt.mesh) == 3
        makie_faces_mesh_orient!(conn,plt)
    end
    (vert,conn,color,colorrange)
end

function makie_faces_mesh(plt)
    d = 2
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
    end
    (vert,conn)
end

function makie_faces_mesh_orient!(conn,plt)
    d=2
    if ! haskey(GT.face_data(plt,d),GT.PLOT_NORMALS_KEY)
        return
    end
    face_n = GT.face_data(plt,d)[GT.PLOT_NORMALS_KEY]
    mesh = plt.mesh
    node_x = GT.node_coordinates(mesh)
    face_nodes = GT.face_nodes(mesh,2)
    nfaces = size(conn,1)
    for face in 1:nfaces
        nodes = face_nodes[face]
        x1 = node_x[nodes[1]]
        x2 = node_x[nodes[2]]
        x3 = node_x[nodes[3]]
        v1 = x2-x1
        v2 = x3-x1
        n1 = cross(v1,v2)
        n2 = face_n[face]
        if dot(n1,n2) < 0
            c3 = conn[face,3]
            conn[face,3] = conn[face,2]
            conn[face,2] = c3
        end
    end
end

function plot_preferred_axis_type(plt)
    d = GT.num_ambient_dims(plt)
    d == 3 ? Makie.Axis3 : Makie.Axis
end

function restrict_to_dim_for_makie(plt,dim_in,dim_plot)
    D = dim_in == Makie.Automatic() ? GT.num_dims(plt.mesh) : dim_in
    d = dim_plot
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
    nnodes = GT.num_nodes(mesh3)
    node_count = zeros(Int32,nnodes)
    dface_nodes = GT.face_nodes(mesh3,d)
    for nodes in dface_nodes
        for node in nodes
            node_count[node] += 1
        end
    end
    newnodes = findall(count->count!=0,node_count)
    newfaces = [ Int[] for _ in 0:d ]
    newfaces[end] = collect(1:GT.num_faces(mesh3,d))
    GT.restrict(plt3,newnodes,newfaces)
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
    if GT.num_faces(plt.mesh,d) == 0
        return (plt,:pink)
    end
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

function setup_plt_changes(plt,vector,scalar,scale,shrink)
    if vector !== nothing
        plt = GT.warp_by_vector(plt,vector;scale)
    end
    if scalar !== nothing
        plt = GT.warp_by_scalar(plt,scalar;scale)
    end
    if shrink !== nothing
        plt = GT.shrink(plt;scale=shrink)
    end
    (plt,)
end


#Makie.@recipe Makie0d begin
#    Makie.documented_attributes(Makie.Scatter)...
#end
#
#function Makie.preferred_axis_type(plot::Makie0d)
#    d = GT.num_ambient_dims(plot[1][])
#    d == 3 ? Makie.Axis3 : Makie.Axis
#end
#
#function Makie.plot!(p::Makie0d{<:Tuple{<:GT.Plot}})
#    valid_attributes = Makie.shared_attributes(p, Makie.Scatter)
#    attrs_in = [:converted_1,:color,:colorrange]
#    attrs_out = :newcolorrange
#    map!(setup_colorrange_impl,p.attributes,attrs_in,attrs_out)
#    attrs_in = [:converted_1,:color]
#    attrs_out = [:points,:newcolor]
#    map!(makie_vertices_impl,p.attributes,attrs_in,attrs_out)
#    color = p.newcolor
#    colorrange = p.newcolorrange
#    points = p.points
#    Makie.scatter!(p,valid_attributes,points;color,colorrange)
#end
#
#Makie.@recipe Makie1d begin
#    Makie.documented_attributes(Makie.LineSegments)...
#end
#
#function Makie.preferred_axis_type(plot::Makie1d)
#    d = GT.num_ambient_dims(plot[1][])
#    d == 3 ? Makie.Axis3 : Makie.Axis
#end
#
#function Makie.plot!(p::Makie1d{<:Tuple{<:GT.Plot}})
#    valid_attributes = Makie.shared_attributes(p, Makie.LineSegments)
#    attrs_in = [:converted_1,:color,:colorrange]
#    attrs_out = :newcolorrange
#    map!(setup_colorrange_impl,p.attributes,attrs_in,attrs_out)
#    attrs_in = [:converted_1,:color]
#    attrs_out = [:points,:newcolor]
#    map!(makie_edges_impl,p.attributes,attrs_in,attrs_out)
#    color = p.newcolor
#    colorrange = p.newcolorrange
#    points = p.points
#    Makie.linesegments!(p,valid_attributes,points;color,colorrange)
#end
#
#Makie.@recipe Makie2d begin
#    Makie.documented_attributes(Makie.Mesh)...
#end
#
#function Makie.preferred_axis_type(plot::Makie2d)
#    d = GT.num_ambient_dims(plot[1][])
#    d == 3 ? Makie.Axis3 : Makie.Axis
#end
#
#function Makie.plot!(p::Makie2d{<:Tuple{<:GT.Plot}})
#    valid_attributes = Makie.shared_attributes(p, Makie.Mesh)
#    attrs_in = [:converted_1,:color,:colorrange]
#    attrs_out = :newcolorrange
#    map!(setup_colorrange_impl,p.attributes,attrs_in,attrs_out)
#    attrs_in = [:converted_1,:color]
#    attrs_out = [:vert,:conn,:newcolor]
#    map!(makie_faces_impl,p.attributes,attrs_in,attrs_out)
#    color = p.newcolor
#    colorrange = p.newcolorrange
#    Makie.mesh!(p,valid_attributes,p.vert,p.conn;color,colorrange)
#    p
#end
#
#Makie.@recipe Makie3d begin
#    Makie.documented_attributes(Makie2d)...
#end
#
#Makie.preferred_axis_type(plot::Makie3d) = Makie.Axis3
#
#function Makie.plot!(p::Makie3d{<:Tuple{<:GT.Plot}})
#    valid_attributes = Makie.shared_attributes(p, Makie2d)
#    attrs_in = [:converted_1,]
#    attrs_out = :plt
#    map!(makie_volumes_impl,p.attributes,attrs_in,attrs_out)
#    plt = p.plt
#    makie2d!(p,valid_attributes,plt)
#end
#
#Makie.@recipe Makie2d1d begin
#    Makie.documented_attributes(Makie1d)...
#end
#
#function Makie.preferred_axis_type(plot::Makie2d1d)
#    d = GT.num_ambient_dims(plot[1][])
#    d == 3 ? Makie.Axis3 : Makie.Axis
#end
#
#function Makie.plot!(p::Makie2d1d{<:Tuple{<:GT.Plot}})
#    valid_attributes = Makie.shared_attributes(p, Makie1d)
#    attrs_in = [:converted_1,]
#    attrs_out = :plt
#    map!(makie_face_edges_impl,p.attributes,attrs_in,attrs_out)
#    plt = p.plt
#    makie1d!(p,valid_attributes,plt)
#end
#
#Makie.@recipe Makie3d1d begin
#    Makie.documented_attributes(Makie2d1d)...
#end
#
#Makie.preferred_axis_type(plot::Makie3d1d) = Makie.Axis3
#
#function Makie.plot!(p::Makie3d1d{<:Tuple{<:GT.Plot}})
#    attrs_in = [:converted_1,]
#    attrs_out = :plt
#    map!(makie_volumes_impl,p.attributes,attrs_in,attrs_out)
#    plt = p.plt
#    valid_attributes = Makie.shared_attributes(p, Makie2d1d)
#    makie2d1d!(p,valid_attributes,plt)
#end
#
#Makie.@recipe MakiePlot begin
#    Makie.documented_attributes(Makie2d)...
#    dim = Makie.Automatic()
#    shrink = false
#    warp_by_vector = nothing
#    warp_by_scalar = nothing
#    warp_scale = 1
#    refinement = nothing
#    strokecolor = nothing
#    strokewidth = 1.5
#    transparency = false
#    overdraw = false
#    fxaa = false
#end
#
#function Makie.preferred_axis_type(plot::MakiePlot)
#    plt = plot[1][]
#    d = GT.num_ambient_dims(plt)
#    d == 3 ? Makie.Axis3 : Makie.Axis
#end
#
#function Makie.plot!(p::MakiePlot{<:Tuple{<:GT.Plot}})
#    attrs_in = [:converted_1,:warp_by_vector,:warp_scale]
#    attrs_out = :plt_warp_vector
#    map!(p.attributes,attrs_in,attrs_out) do plt,vector,scale
#        GT.warp_by_vector(plt,vector;scale)
#    end
#    attrs_in = [:plt_warp_vector,:warp_by_scalar,:warp_scale]
#    attrs_out = :plt_warp_scalar
#    map!(p.attributes,attrs_in,attrs_out) do plt,scalar,scale
#        GT.warp_by_scalar(plt,scalar;scale)
#    end
#    attrs_in = [:plt_warp_scalar,:shrink]
#    attrs_out = :plt_shrink
#    map!(p.attributes,attrs_in,attrs_out) do plt,scale
#        GT.shrink(plt;scale)
#    end
#    plt = p.plt_shrink
#    plt_now = plt[]
#    dim_now = p.dim[]
#    d_now = GT.num_dims(plt_now.mesh)
#    dim_now = dim_now == Makie.Automatic() ? d_now : dim_now
#    if 0 in dim_now 
#        if p.color[] !== nothing
#            valid_attributes = Makie.shared_attributes(p,Makie0d)
#            makie0d!(p,valid_attributes,plt)
#        end
#    end
#    if 1 in dim_now
#        if p.color[] !== nothing
#            valid_attributes = Makie.shared_attributes(p,Makie1d)
#            makie1d!(p,valid_attributes,plt)
#        end
#    end
#    if 2 in dim_now
#        if p.color[] !== nothing
#            valid_attributes = Makie.shared_attributes(p,Makie2d)
#            makie2d!(p,valid_attributes,plt)
#        end
#        if p.strokecolor[] !== nothing
#            color = p.strokecolor
#            linewidth = p.strokewidth
#            valid_attributes = Makie.shared_attributes(p,Makie2d1d)
#            makie2d1d!(p,valid_attributes,plt;color,linewidth)
#        end
#    end
#    if 3 in dim_now
#        if p.color[] !== nothing
#            valid_attributes = Makie.shared_attributes(p,Makie3d)
#            makie3d!(p,valid_attributes,plt)
#        end
#        if p.strokecolor[] !== nothing
#            color = p.strokecolor
#            linewidth = p.strokewidth
#            valid_attributes = Makie.shared_attributes(p,Makie3d1d)
#            makie3d1d!(p,valid_attributes,plt;color,linewidth)
#        end
#    end
#    p
#end
#
#Makie.plottype(::GT.Plot) = MakiePlot
#Makie.plottype(::GT.PPlot) = MakiePlot
#Makie.plottype(::GT.AbstractMesh) = MakiePlot
#Makie.plottype(::GT.PMesh) = MakiePlot
#
#function Makie.convert_arguments(::Type{<:MakiePlot},mesh::GT.AbstractMesh)
#    plt = GT.plot(mesh)
#    (plt,)
#end
#
#function Makie.convert_arguments(::Type{<:MakiePlot},mesh::GT.PMesh)
#    plt = GT.plot(mesh)
#    (plt,)
#end
#
#Makie.plottype(::GT.AbstractDomain) = MakiePlot
#
#function Makie.plot!(p::MakiePlot{<:Tuple{<:GT.AbstractDomain}})
#    attrs_in = [:converted_1,:refinement,:color,:warp_by_vector,:warp_by_scalar]
#    attrs_out = [:plt,:newcolor,:new_warp_by_vector,:new_warp_by_scalar]
#    map!(p.attributes,attrs_in,attrs_out) do domain,refinement,color,vec,scal
#        plt = plot_for_makie(domain,refinement)
#        if isa(color,GT.AbstractQuantity) || isa(color,Function)
#            label = string(gensym())
#            GT.plot!(plt,color;label)
#            newcolor = GT.node_color(plt,label)
#        else
#            newcolor = color
#        end
#        if isa(vec,GT.AbstractQuantity) || isa(vec,Function)
#            label = string(gensym())
#            GT.plot!(plt,vec;label)
#            newvec = GT.NodeData(label)
#        else
#            newvec = vec
#        end
#        if isa(scal,GT.AbstractQuantity) || isa(scal,Function)
#            label = string(gensym())
#            GT.plot!(plt,scal;label)
#            newscal = GT.NodeData(label)
#        else
#            newscal = scal
#        end
#        (plt,newcolor,newvec,newscal)
#    end
#    color = p.newcolor
#    warp_by_vector = p.new_warp_by_vector
#    warp_by_scalar = p.new_warp_by_scalar
#    valid_attributes = Makie.shared_attributes(p,MakiePlot)
#    plot!(p,valid_attributes,p.plt;color,warp_by_vector,warp_by_scalar)
#    p
#end
#
#function plot_for_makie(domain::GT.AbstractDomain,refinement)
#    if GT.num_dims(domain) == 3
#        domain2 = GT.boundary(domain)
#    else
#        domain2 = domain
#    end
#    plt = GT.plot(domain2;refinement)
#    if GT.num_dims(domain) == 3
#        face_n = collect_face_normals(domain2)
#        vface_face = plt.cache.glue.parent_face
#        vface_n = face_n[vface_face]
#        GT.face_data(plt,2)[GT.PLOT_NORMALS_KEY] = vface_n
#    end
#    plt
#end
#
#function collect_face_normals(Γ::GT.AbstractDomain)
#    dΓ = GT.measure(Γ,0)
#    dface_point_n = GT.unit_normal_accessor(dΓ)
#    Tn = typeof(GT.prototype(dface_point_n))
#    ndfaces = GT.num_faces(Γ)
#    dface_to_n = zeros(Tn,ndfaces)
#    for dface in 1:ndfaces
#        dface_to_n[dface] = dface_point_n(dface,1)(1)
#    end
#    dface_to_n
#end
#
## TODO not sure about this
## what if u is not scalar-valued?
##Makie.plottype(::GT.AbstractQuantity) = MakiePlot
##
##function Makie.plot!(sc::MakiePlot{<:Tuple{<:GT.AbstractQuantity}})
##    u = sc[1]
##    valid_attributes = Makie.shared_attributes(sc, MakiePlot)
##    dom = Makie.lift(domain,u)
##    valid_attributes[:color] = u
##    makieplot!(sc,valid_attributes,dom)
##end
#
##function vector_of_observables(a)
##    #TODO remove trick for DebugArray
##    function start(v)
##        first(v)
##    end
##    function rest(v)
##        v[2:end]
##    end
##    if length(a[]) == 1
##        b = Makie.lift(first,a)
##        return [b,]
##    else
##        b = Makie.lift(start,a)
##        c = Makie.lift(rest,a)
##        return [b,vector_of_observables(c)...]
##    end
##end
#
#function makie_volumes_impl(plt::GT.Plot)
#    GT.skin(plt)
#end
#
#function setup_colorrange_impl(plt,color,colorrange)
#    if colorrange != Makie.Automatic()
#        return colorrange
#    end
#    if isa(color,GT.FaceData)
#        colorrange = GT.face_colorrange(plt,color.name)
#    end
#    if isa(color,GT.NodeData)
#        colorrange = GT.node_colorrange(plt,color.name)
#    end
#    colorrange
#end
#
#function setup_colors_impl(plt,color,d)
#    if GT.num_faces(plt.mesh,d) == 0
#        return (plt,:pink)
#    end
#    if isa(color,GT.FaceData)
#        if d == 2
#            plt = GT.shrink(plt;scale=1)
#        end
#        color2 = GT.face_color(plt,color.name,d)
#    end
#    if isa(color,GT.NodeData)
#        color2 = GT.node_color(plt,color.name)
#    end
#    if isa(color,GT.NodeData) || isa(color,GT.FaceData)
#        face_to_nodes = GT.face_nodes(plt.mesh,d)
#        nfaces = length(face_to_nodes)
#        if d == 0
#            color = similar(color2,nfaces)
#            for (face,nodes) in enumerate(face_to_nodes)
#                node = first(nodes)
#                color[face] = color2[node]
#            end
#        elseif d == 1
#            color = similar(color2,2*nfaces)
#            k = 0
#            for nodes in face_to_nodes
#                for node in nodes
#                    k+=1
#                    color[k] = color2[node]
#                end
#            end
#        elseif d == 2
#            color = color2
#        else
#            error()
#        end
#    end
#    plt,color
#end
#
#function makie_faces_impl(plt,color)
#    D = GT.num_ambient_dims(plt.mesh)
#    # Makie seems to not like vertices not touched by any element
#    # for the shading.
#    plt = GT.restrict_to_dim(plt,2)
#    plt = GT.shrink(plt;scale=1)
#    plt = GT.simplexify(plt)
#    d = 2
#    plt,color = setup_colors_impl(plt,color,d)
#    vert,conn = makie_faces_mesh(plt)
#    if D == 3
#        makie_faces_mesh_orient!(conn,plt)
#    end
#    (vert,conn,color)
#end
#
#function makie_faces_mesh(plt)
#    d = 2
#    mesh = plt.mesh
#    D = GT.num_ambient_dims(mesh)
#    nnodes = GT.num_nodes(mesh)
#    vert = zeros(Float64,nnodes,D)
#    node_to_x = GT.node_coordinates(mesh)
#    for (node,x) in enumerate(node_to_x)
#        vert[node,:] = x
#    end
#    nfaces = GT.num_faces(mesh,d)
#    face_to_nodes = GT.face_nodes(mesh,d)
#    conn = zeros(Int32,nfaces,3)
#    for (face,nodes) in enumerate(face_to_nodes)
#        conn[face,:] = nodes
#    end
#    if nfaces == 0
#        vert = ones(Float64,3,D)
#        conn = zeros(Int32,1,3)
#        conn[:] = 1:3
#    end
#    (vert,conn)
#end
#
#function makie_faces_mesh_orient!(conn,plt)
#    d=2
#    mesh = plt.mesh
#    node_x = GT.node_coordinates(mesh)
#    face_nodes = GT.face_nodes(mesh,2)
#    nfaces = size(conn,1)
#    face_n = GT.face_data(plt,d)[GT.PLOT_NORMALS_KEY]
#    for face in 1:nfaces
#        nodes = face_nodes[face]
#        x1 = node_x[nodes[1]]
#        x2 = node_x[nodes[2]]
#        x3 = node_x[nodes[3]]
#        v1 = x2-x1
#        v2 = x3-x1
#        n1 = cross(v1,v2)
#        n2 = face_n[face]
#        if dot(n1,n2) < 0
#            c3 = conn[face,3]
#            conn[face,3] = conn[face,2]
#            conn[face,2] = c3
#        end
#    end
#end
#
#function makie_face_edges_impl(plt)
#    # TODO maybe already complexified
#    D=2
#    d=1
#    mesh2 = GT.complexify(GT.restrict_to_dim(plt.mesh,D))
#    topo = GT.topology(mesh2)
#    edge_to_faces = GT.face_incidence(topo,d,D)
#    K = Any
#    facedata = Dict{String,K}()
#    for (k,v) in GT.face_data(plt,D;merge_dims=true)
#        facedata[k] = map(faces -> sum(v[faces])/length(faces) ,edge_to_faces)
#    end
#    mesh3 = GT.restrict_to_dim(mesh2,d)
#    fd = map(0:d) do i
#        if i == d
#            facedata
#        else
#            typeof(facedata)()
#        end
#    end
#    plt3 = GT.Plot(mesh3,fd,GT.node_data(plt))
#    plt3
#end
#
#function makie_edges_impl(plt,color)
#    d = 1
#    plt,color = setup_colors_impl(plt,color,d)
#    mesh = plt.mesh
#    nedges = GT.num_faces(mesh,d)
#    node_to_x = GT.node_coordinates(mesh)
#    edge_to_nodes = GT.face_nodes(mesh,d)
#    T = eltype(eltype(node_to_x))
#    S = GT.num_ambient_dims(mesh)
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
#    (p,color)
#end
#
#function makie_vertices_impl(plt,color)
#    d = 0
#    plt,color = setup_colors_impl(plt,color,d)
#    mesh = plt.mesh
#    nedges = GT.num_faces(mesh,d)
#    node_to_x = GT.node_coordinates(mesh)
#    edge_to_nodes = GT.face_nodes(mesh,d)
#    T = eltype(eltype(node_to_x))
#    S = GT.num_ambient_dims(mesh)
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
#    (p,color)
#end
#
## Vector-valued
#
#function Makie.plot!(p::Arrows{<:Tuple{<:GT.Plot,GT.NodeData}})
#    valid_attributes = Makie.shared_attributes(p,Arrows)
#    attrs_in = [:converted_1,:converted_2,:color]
#    attrs_out = [:coords,:vecs,:newcolor]
#    map!(makie_arrows_impl,p.attributes,attrs_in,attrs_out)
#    color = p.newcolor
#    Makie.arrows!(p,valid_attributes,p.coords,p.vecs;color)
#end
#
#function makie_arrows_impl(plt,vecs::GT.NodeData,color)
#    D = GT.num_ambient_dims(plt.mesh)
#    x = GT.node_coordinates(plt.mesh)
#    nnodes = length(x)
#    if D == 2
#        node_to_coord = [ Makie.Point2f(xi) for xi in x]
#    elseif D == 3
#        node_to_coord = [ Makie.Point3f(xi) for xi in x]
#    else
#        error("not implemented")
#    end
#    node_to_vec = plt.node_data[vecs.name]
#    if isa(color,GT.NodeData)
#        color2 = plt.node_data[color.name]
#    else
#        color2 = color
#    end
#    (node_to_coord,node_to_vec,color2)
#end
#
#function Makie.plot!(sc::Arrows{<:Tuple{<:GT.AbstractField}})
#    valid_attributes = Makie.shared_attributes(p, Arrows)
#    attrs_in = [:converted_1,:color]
#    attrs_out = [:plt,:vecs,:newcolor]
#    args = map!(p.attributes,attrs_in,attrs_out) do q,color
#        dom = GT.domain(q)
#        plt = GT.plot(dom)
#        label_q = string(gensym())
#        GT.plot!(plt,q;label=label_q)
#        if isa(color,GT.AbstractQuantity) || isa(color,Function)
#            label_color = string(gensym())
#            GT.plot!(plt,color;label=label_color)
#            color = GT.node_color(plt,label_color)
#        end
#        vecs = GT.NodeData(label_q)
#        (plt,vecs,color)
#    end
#    plt = p.plt
#    vecs = p.vecs
#    color = p.newcolor
#    arrows!(p,valid_attributes,plt,vecs;color)
#end

end # modules

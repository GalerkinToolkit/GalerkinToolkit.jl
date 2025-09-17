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

#FaceColor -> FaceColor
#NodeColor -> NodeColor


import GalerkinToolkit as GT
import GalerkinToolkit: makie_vertices, makie_vertices!
import GalerkinToolkit: makie_edges, makie_edges!
import GalerkinToolkit: makie_surfaces, makie_surfaces!
import GalerkinToolkit: makie_arrows2d, makie_arrows2d!
import GalerkinToolkit: makie_arrows3d, makie_arrows3d!
using PartitionedArrays
using Makie

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

Makie.@recipe Makie_vertices begin
    Makie.documented_attributes(Makie.Scatter)...
    shared_attributes_mixin()...
end

function Makie.preferred_axis_type(p::Makie_vertices)
    plt = p[1][]
    plot_preferred_axis_type(plt)
end

function Makie.plot!(p::Makie_vertices{<:Tuple{<:GT.Plot}})
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

Makie.@recipe Makie_edges begin
    Makie.documented_attributes(Makie.LineSegments)...
    shared_attributes_mixin()...
end

function Makie.preferred_axis_type(p::Makie_edges)
    plt = p[1][]
    plot_preferred_axis_type(plt)
end

function Makie.plot!(p::Makie_edges{<:Tuple{<:GT.Plot}})
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

Makie.@recipe Makie_surfaces begin
    Makie.documented_attributes(Makie.Mesh)...
    shared_attributes_mixin()...
end

function Makie.preferred_axis_type(p::Makie_surfaces)
    plt = p[1][]
    plot_preferred_axis_type(plt)
end

function Makie.plot!(p::Makie_surfaces{<:Tuple{<:GT.Plot}})
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
    vert = zeros(SVector{D,Float64},nnodes)
    node_to_x = GT.node_coordinates(mesh)
    for (node,x) in enumerate(node_to_x)
        vert[node] = x
    end
    nfaces = GT.num_faces(mesh,d)
    face_to_nodes = GT.face_nodes(mesh,d)
    conn = zeros(Int32,nfaces,3)
    for (face,nodes) in enumerate(face_to_nodes)
        conn[face,:] = nodes
    end
    if nfaces == 0
        vert = zeros(SVector{D,Float64},3)
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

Makie.@recipe Makie_arrows2d begin
    Makie.documented_attributes(Makie.Arrows2D)...
    shared_attributes_mixin()...
end

function Makie.preferred_axis_type(p::Makie_arrows2d)
    plt = p[1][]
    plot_preferred_axis_type(plt)
end

function Makie.plot!(p::Makie_arrows2d{<:Tuple{<:GT.Plot,<:GT.NodeColor}})
    attrs_in = [:converted_1,:warp_by_vector,:warp_by_scalar,:warp_scale,:shrink]
    attrs_out = [:plt]
    map!(setup_plt_changes,p.attributes,attrs_in,attrs_out)
    valid_attributes = Makie.shared_attributes(p,Makie.Arrows2D)
    attrs_in = [:plt,:converted_2,:color]
    attrs_out = [:coords,:vecs,:newcolor]
    map!(makie_arrows_impl,p.attributes,attrs_in,attrs_out)
    color = p.newcolor
    Makie.arrows2d!(p,valid_attributes,p.coords,p.vecs;color)
end

Makie.@recipe Makie_arrows3d begin
    Makie.documented_attributes(Makie.Arrows3D)...
    shared_attributes_mixin()...
end

function Makie.preferred_axis_type(p::Makie_arrows3d)
    plt = p[1][]
    plot_preferred_axis_type(plt)
end

function Makie.plot!(p::Makie_arrows3d{<:Tuple{<:GT.Plot,<:GT.NodeColor}})
    attrs_in = [:converted_1,:warp_by_vector,:warp_by_scalar,:warp_scale,:shrink]
    attrs_out = [:plt]
    map!(setup_plt_changes,p.attributes,attrs_in,attrs_out)
    valid_attributes = Makie.shared_attributes(p,Makie.Arrows2D)
    attrs_in = [:plt,:converted_2,:color]
    attrs_out = [:coords,:vecs,:newcolor]
    map!(makie_arrows_impl,p.attributes,attrs_in,attrs_out)
    color = p.newcolor
    Makie.arrows3d!(p,valid_attributes,p.coords,p.vecs;color)
end

function makie_arrows_impl(plt,vecs::GT.NodeColor,color)
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
    if isa(color,GT.NodeColor)
        color2 = plt.node_data[color.name]
    else
        color2 = color
    end
    (node_to_coord,node_to_vec,color2)
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
    if isa(color,GT.FaceColor)
        colorrange = GT.face_colorrange(plt,color.name)
    end
    if isa(color,GT.NodeColor)
        colorrange = GT.node_colorrange(plt,color.name)
    end
    colorrange
end

function setup_colors_impl(plt,color,d)
    if GT.num_faces(plt.mesh,d) == 0
        return (plt,:pink)
    end
    if isa(color,GT.FaceColor)
        if d == 2
            plt = GT.shrink(plt;scale=1)
        end
        color2 = GT.face_color(plt,color.name,d)
    end
    if isa(color,GT.NodeColor)
        color2 = GT.node_color(plt,color.name)
    end
    if isa(color,GT.NodeColor) || isa(color,GT.FaceColor)
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
    if GT.is_partitioned(plt)
        plt = GT.centralize(plt)
    end
    (plt,)
end

# Conversions

function Makie.convert_arguments(::Type{<:Makie_surfaces},mesh::GT.AbstractMesh)
    plt = GT.plot(mesh)
    (plt,)
end

function Makie.convert_arguments(::Type{<:Makie_edges},mesh::GT.AbstractMesh)
    plt = GT.plot(mesh)
    (plt,)
end

function Makie.convert_arguments(::Type{<:Makie_vertices},mesh::GT.AbstractMesh)
    plt = GT.plot(mesh)
    (plt,)
end

#function Makie.convert_arguments(::Type{<:Makie_surfaces},mesh::GT.AbstractPMesh)
#    plt = PartitionedArrays.centralize(GT.plot(mesh))
#    (plt,)
#end
#
#function Makie.convert_arguments(::Type{<:Makie_edges},mesh::GT.AbstractPMesh)
#    plt = PartitionedArrays.centralize(GT.plot(mesh))
#    (plt,)
#end
#
#function Makie.convert_arguments(::Type{<:Makie_vertices},mesh::GT.AbstractPMesh)
#    plt = PartitionedArrays.centralize(GT.plot(mesh))
#    (plt,)
#end

function Makie.convert_single_argument(::Type{<:Makie.AbstractPlot},pplot::GT.PPlot)
    plt = PartitionedArrays.centralize(pplot)
    plt
end

function Makie.plot!(p::Makie_surfaces{<:Tuple{<:GT.AbstractDomain}})
    attrs_in = [:converted_1,:refinement,:color,:warp_by_vector,:warp_by_scalar]
    attrs_out = [:plt,:newcolor,:new_warp_by_vector,:new_warp_by_scalar]
    map!(setup_makie_domain,p.attributes,attrs_in,attrs_out)
    color = p.newcolor
    warp_by_vector = p.new_warp_by_vector
    warp_by_scalar = p.new_warp_by_scalar
    valid_attributes = Makie.shared_attributes(p,Makie_surfaces)
    makie_surfaces!(p,valid_attributes,p.plt;color,warp_by_vector,warp_by_scalar)
    p
end

function Makie.plot!(p::Makie_edges{<:Tuple{<:GT.AbstractDomain}})
    attrs_in = [:converted_1,:refinement,:color,:warp_by_vector,:warp_by_scalar]
    attrs_out = [:plt,:newcolor,:new_warp_by_vector,:new_warp_by_scalar]
    map!(setup_makie_domain,p.attributes,attrs_in,attrs_out)
    color = p.newcolor
    warp_by_vector = p.new_warp_by_vector
    warp_by_scalar = p.new_warp_by_scalar
    valid_attributes = Makie.shared_attributes(p,Makie_edges)
    makie_edges!(p,valid_attributes,p.plt;color,warp_by_vector,warp_by_scalar)
    p
end

function Makie.plot!(p::Makie_vertices{<:Tuple{<:GT.AbstractDomain}})
    attrs_in = [:converted_1,:refinement,:color,:warp_by_vector,:warp_by_scalar]
    attrs_out = [:plt,:newcolor,:new_warp_by_vector,:new_warp_by_scalar]
    map!(setup_makie_domain,p.attributes,attrs_in,attrs_out)
    color = p.newcolor
    warp_by_vector = p.new_warp_by_vector
    warp_by_scalar = p.new_warp_by_scalar
    valid_attributes = Makie.shared_attributes(p,Makie_vertices)
    makie_vertices!(p,valid_attributes,p.plt;color,warp_by_vector,warp_by_scalar)
    p
end

function setup_makie_domain(domain,refinement,color,vec,scal)
    plt = plot_for_makie(domain,refinement)
    if isa(color,GT.AbstractQuantity) || isa(color,Function)
        label = string(gensym())
        GT.plot!(plt,color;label)
        newcolor = GT.NodeColor(label)#GT.node_color(plt,label)
    else
        newcolor = color
    end
    if isa(vec,GT.AbstractQuantity) || isa(vec,Function)
        label = string(gensym())
        GT.plot!(plt,vec;label)
        newvec = GT.NodeColor(label)
    else
        newvec = vec
    end
    if isa(scal,GT.AbstractQuantity) || isa(scal,Function)
        label = string(gensym())
        GT.plot!(plt,scal;label)
        newscal = GT.NodeColor(label)
    else
        newscal = scal
    end
    (plt,newcolor,newvec,newscal)
end

function plot_for_makie(domain::GT.AbstractDomain,refinement)
    return GT.plot(domain;refinement) 
    # TODO not sure about this one
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
    #if GT.num_dims(domain) == 3
    #    domain2 = GT.boundary(domain)
    #else
    #    domain2 = domain
    #end
    #plt = GT.plot(domain2;refinement)
    #if GT.num_dims(domain) == 3
    #    face_n = collect_face_normals(domain2)
    #    vface_face = plt.cache.glue.parent_face
    #    vface_n = face_n[vface_face]
    #    GT.face_data(plt,2)[GT.PLOT_NORMALS_KEY] = vface_n
    #end
    #plt
end

function Makie.plot!(p::Makie_arrows2d{<:Tuple{<:GT.AbstractDomain,<:GT.AbstractField}})
    attrs_in = [:converted_1,:converted_2,:color,:refinement,:warp_by_vector,:warp_by_scalar]
    attrs_out = [:plt,:vecs,:newcolor,:new_warp_by_vector,:new_warp_by_scalar]
    args = map!(setup_domain_arrows,p.attributes,attrs_in,attrs_out)
    plt = p.plt
    vecs = p.vecs
    color = p.newcolor
    warp_by_vector = p.new_warp_by_vector
    warp_by_scalar = p.new_warp_by_scalar
    valid_attributes = Makie.shared_attributes(p, Makie_arrows2d)
    makie_arrows2d!(p,valid_attributes,plt,vecs;color,warp_by_vector,warp_by_scalar)
end

function Makie.plot!(p::Makie_arrows3d{<:Tuple{<:GT.AbstractDomain,<:GT.AbstractField}})
    attrs_in = [:converted_1,:converted_2,:color,:refinement,:warp_by_vector,:warp_by_scalar]
    attrs_out = [:plt,:vecs,:newcolor,:new_warp_by_vector,:new_warp_by_scalar]
    args = map!(setup_domain_arrows,p.attributes,attrs_in,attrs_out)
    plt = p.plt
    vecs = p.vecs
    color = p.newcolor
    warp_by_vector = p.new_warp_by_vector
    warp_by_scalar = p.new_warp_by_scalar
    valid_attributes = Makie.shared_attributes(p, Makie_arrows3d)
    makie_arrows3d!(p,valid_attributes,plt,vecs;color,warp_by_vector,warp_by_scalar)
end

function setup_domain_arrows(dom,q,color,refinement,vec,scal)
    plt = plot_for_makie(dom,refinement)
    label_q = string(gensym())
    GT.plot!(plt,q;label=label_q)
    if isa(color,GT.AbstractQuantity) || isa(color,Function)
        label_color = string(gensym())
        GT.plot!(plt,color;label=label_color)
        color = GT.node_color(plt,label_color)
    end
    if isa(vec,GT.AbstractQuantity) || isa(vec,Function)
        label = string(gensym())
        GT.plot!(plt,vec;label)
        newvec = GT.NodeColor(label)
    else
        newvec = vec
    end
    if isa(scal,GT.AbstractQuantity) || isa(scal,Function)
        label = string(gensym())
        GT.plot!(plt,scal;label)
        newscal = GT.NodeColor(label)
    else
        newscal = scal
    end
    varrow = GT.NodeColor(label_q)
    (plt,varrow,color,newvec,newscal)
end

end # modules

module GalerkinToolkitMakieExt

using GalerkinToolkit
import GalerkinToolkit: makieplot, makieplot!
import GalerkinToolkit: makie0d, makie0d!, makie1d, makie1d!
import GalerkinToolkit: makie2d, makie2d!, makie2d1d, makie2d1d!
import GalerkinToolkit: makie3d, makie3d!, makie3d1d, makie3d1d!
using PartitionedArrays
using Makie

Makie.@recipe(MakiePlot) do scene
    t1 = Makie.Theme(
        dim = Makie.Automatic(),
        shrink = false,
        strokecolor = nothing,
        color = :lightblue,
        colormap = :bluesreds,
        shading = Makie.NoShading,
        cycle = nothing,
    )
    t2 = Makie.default_theme(scene, Makie.Mesh)
    merge(t1, t2)
end

Makie.preferred_axis_type(plot::MakiePlot) = Makie.LScene

function Makie.plot!(sc::MakiePlot{<:Tuple{<:GalerkinToolkit.Plot}})
    plt = sc[1]
    # TODO these are not reactive
    if sc[:shrink][] != false
        scale = sc[:shrink]
        plt = Makie.lift((a, b) -> GalerkinToolkit.shrink(a; scale = b), plt, scale)
    end
    dim = sc[:dim][]
    d = GalerkinToolkit.num_dims(plt[].mesh)
    strokecolor = sc[:strokecolor][]
    if dim == Makie.Automatic()
        dim = d
    end
    cmp = y -> y in dim
    if cmp(0)
        valid_attributes = Makie.shared_attributes(sc, Makie0d)
        makie0d!(sc, valid_attributes, plt)
    end
    if cmp(1)
        valid_attributes = Makie.shared_attributes(sc, Makie1d)
        makie1d!(sc, valid_attributes, plt)
    end
    if cmp(2)
        valid_attributes = Makie.shared_attributes(sc, Makie2d)
        makie2d!(sc, valid_attributes, plt)
        if strokecolor !== nothing
            valid_attributes = Makie.shared_attributes(sc, Makie2d1d)
            valid_attributes[:color] = sc[:strokecolor]
            makie2d1d!(sc, valid_attributes, plt)
        end
    end
    if cmp(3)
        valid_attributes = Makie.shared_attributes(sc, Makie3d)
        makie3d!(sc, valid_attributes, plt)
        if strokecolor !== nothing
            valid_attributes = Makie.shared_attributes(sc, Makie3d1d)
            valid_attributes[:color] = sc[:strokecolor]
            makie3d1d!(sc, valid_attributes, plt)
        end
    end
end

function Makie.plot!(sc::MakiePlot{<:Tuple{<:GalerkinToolkit.PPlot}})
    valid_attributes = Makie.shared_attributes(sc, MakiePlot)
    pplt = sc[1]
    plts = makie_pplot_setup(pplt)
    if plts !== nothing
        foreach(plts) do plt
            makieplot!(sc, valid_attributes, plt)
        end
    end
end

Makie.plottype(::GalerkinToolkit.Plot) = MakiePlot
Makie.plottype(::GalerkinToolkit.PPlot) = MakiePlot

Makie.plottype(::GalerkinToolkit.AbstractMesh) = MakiePlot

function Makie.convert_arguments(::Type{<:MakiePlot}, mesh::GalerkinToolkit.AbstractMesh)
    plt = GalerkinToolkit.plot(mesh)
    (plt,)
end

Makie.plottype(::GalerkinToolkit.AbstractDomain) = MakiePlot

function Makie.plot!(sc::MakiePlot{<:Tuple{<:GalerkinToolkit.AbstractDomain}})
    dom = sc[1]
    valid_attributes = Makie.shared_attributes(sc, MakiePlot)
    color = valid_attributes[:color]
    args = Makie.lift(dom, color) do dom, color
        if isa(color, GalerkinToolkit.AbstractQuantity)
            label = string(gensym())
            plt = GalerkinToolkit.plot(dom)
            GalerkinToolkit.plot!(plt, color; label)
            color = GalerkinToolkit.NodeData(label)
        else
            plt = GalerkinToolkit.plot(dom)
        end
        (; plt, color)
    end
    plt = Makie.lift(i -> i.plt, args)
    color = Makie.lift(i -> i.color, args)
    valid_attributes[:color] = color
    makieplot!(sc, valid_attributes, plt)
end

# TODO not sure about this
# what if u is not scalar-valued?
#Makie.plottype(::GalerkinToolkit.AbstractQuantity) = MakiePlot
#
#function Makie.plot!(sc::MakiePlot{<:Tuple{<:GalerkinToolkit.AbstractQuantity}})
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
        b = Makie.lift(first, a)
        return [b]
    else
        b = Makie.lift(start, a)
        c = Makie.lift(rest, a)
        return [b, vector_of_observables(c)...]
    end
end

function makie_pplot_setup(pplt)
    obs = Makie.Observable{Any}(nothing)
    function update!(pplt)
        map_main(PartitionedArrays.gather(pplt.partition)) do myplts
            obs[] = myplts
        end
    end
    Makie.on(update!, pplt)
    update!(pplt[])
    if obs[] !== nothing
        plts = vector_of_observables(obs)
    else
        plts = nothing
    end
    plts
end

Makie.@recipe(Makie3d) do scene
    dt = Makie.default_theme(scene, Makie.Mesh)
    dt
end

Makie.preferred_axis_type(plot::Makie3d) = Makie.LScene

function Makie.plot!(sc::Makie3d{<:Tuple{<:GalerkinToolkit.Plot}})
    plt = sc[1]
    valid_attributes = Makie.shared_attributes(sc, Makie2d)
    colorrange = valid_attributes[:colorrange]
    color = valid_attributes[:color]
    colorrange = Makie.lift(setup_colorrange_impl, plt, color, colorrange)
    valid_attributes[:colorrange] = colorrange
    plt2 = Makie.lift(makie_volumes_impl, plt)
    makie2d!(sc, valid_attributes, plt2)
end

function Makie.plot!(sc::Makie3d{<:Tuple{<:GalerkinToolkit.PPlot}})
    valid_attributes = Makie.shared_attributes(sc, Makie3d)
    pplt = sc[1]
    plts = makie_pplot_setup(pplt)
    if plts !== nothing
        foreach(plts) do plt
            makie3d!(sc, valid_attributes, plt)
        end
    end
end

function makie_volumes_impl(plt::GalerkinToolkit.Plot)
    @assert GalerkinToolkit.num_dims(plt.mesh) == 3
    D = 3
    d = 2
    mesh, = GalerkinToolkit.complexify(GalerkinToolkit.restrict_to_dim(plt.mesh, D))
    topo = GalerkinToolkit.topology(mesh)
    face_to_cells = GalerkinToolkit.face_incidence(topo, d, D)
    face_isboundary = map(cells -> length(cells) == 1, face_to_cells)
    mesh2 = GalerkinToolkit.restrict_to_dim(mesh, d)
    newnodes = 1:GalerkinToolkit.num_nodes(mesh2)
    newfaces = [Int[] for _ = 0:d]
    newfaces[end] = findall(face_isboundary)
    mesh3 = GalerkinToolkit.restrict(mesh2, newnodes, newfaces)
    face_to_cell = map(first, face_to_cells)
    newface_to_cell = face_to_cell[newfaces[end]]
    celldata = copy(GalerkinToolkit.face_data(plt, D; merge_dims = true))
    for (k, v) in celldata
        celldata[k] = v[newface_to_cell]
    end
    fd = map(0:d) do i
        if i == d
            celldata
        else
            typeof(celldata)()
        end
    end
    plt3 = GalerkinToolkit.Plot(mesh3, fd, GalerkinToolkit.node_data(plt))
    plt3
end

Makie.@recipe(Makie3d1d) do scene
    dt = Makie.default_theme(scene, Makie2d1d)
    dt
end

Makie.preferred_axis_type(plot::Makie3d1d) = Makie.LScene

function Makie.plot!(sc::Makie3d1d{<:Tuple{<:GalerkinToolkit.Plot}})
    plt = sc[1]
    valid_attributes = Makie.shared_attributes(sc, Makie2d1d)
    colorrange = valid_attributes[:colorrange]
    color = valid_attributes[:color]
    colorrange = Makie.lift(setup_colorrange_impl, plt, color, colorrange)
    valid_attributes[:colorrange] = colorrange
    plt2 = Makie.lift(makie_volumes_impl, plt)
    makie2d1d!(sc, valid_attributes, plt2)
end

function Makie.plot!(sc::Makie3d1d{<:Tuple{<:GalerkinToolkit.PPlot}})
    valid_attributes = Makie.shared_attributes(sc, Makie3d1d)
    pplt = sc[1]
    plts = makie_pplot_setup(pplt)
    if plts !== nothing
        foreach(plts) do plt
            makie3d1d!(sc, valid_attributes, plt)
        end
    end
end

Makie.@recipe(Makie2d) do scene
    dt = Makie.default_theme(scene, Makie.Mesh)
    dt
end

Makie.preferred_axis_type(plot::Makie2d) = Makie.LScene

function Makie.plot!(sc::Makie2d{<:Tuple{<:GalerkinToolkit.Plot}})
    plt = sc[1]
    valid_attributes = Makie.shared_attributes(sc, Makie.Mesh)
    color = valid_attributes[:color]
    colorrange = valid_attributes[:colorrange]
    colorrange = Makie.lift(setup_colorrange_impl, plt, color, colorrange)
    args = Makie.lift(makie_faces_impl, plt, color)
    vert = Makie.lift(a -> a.vert, args)
    conn = Makie.lift(a -> a.conn, args)
    color = Makie.lift(a -> a.color, args)
    valid_attributes[:color] = color
    valid_attributes[:colorrange] = colorrange
    Makie.mesh!(sc, valid_attributes, vert, conn)
end

function Makie.plot!(sc::Makie2d{<:Tuple{<:GalerkinToolkit.PPlot}})
    valid_attributes = Makie.shared_attributes(sc, Makie2d)
    pplt = sc[1]
    plts = makie_pplot_setup(pplt)
    if plts !== nothing
        foreach(plts) do plt
            makie2d!(sc, valid_attributes, plt)
        end
    end
end

function setup_colorrange_impl(plt, color, colorrange)
    if colorrange != Makie.Automatic()
        return colorrange
    end
    if isa(color, GalerkinToolkit.FaceData)
        colorrange = GalerkinToolkit.face_colorrange(plt, color.name)
    end
    if isa(color, GalerkinToolkit.NodeData)
        colorrange = GalerkinToolkit.node_colorrange(plt, color.name)
    end
    colorrange
end

function setup_colors_impl(plt, color, d)
    if isa(color, GalerkinToolkit.FaceData)
        if d == 2
            plt = GalerkinToolkit.shrink(plt; scale = 1)
        end
        color2 = GalerkinToolkit.face_color(plt, color.name, d)
    end
    if isa(color, GalerkinToolkit.NodeData)
        color2 = GalerkinToolkit.node_color(plt, color.name)
    end
    if isa(color, GalerkinToolkit.NodeData) || isa(color, GalerkinToolkit.FaceData)
        face_to_nodes = GalerkinToolkit.face_nodes(plt.mesh, d)
        nfaces = length(face_to_nodes)
        if d == 0
            color = similar(color2, nfaces)
            for (face, nodes) in enumerate(face_to_nodes)
                node = first(nodes)
                color[face] = color2[node]
            end
        elseif d == 1
            color = similar(color2, 2 * nfaces)
            k = 0
            for nodes in face_to_nodes
                for node in nodes
                    k += 1
                    color[k] = color2[node]
                end
            end
        elseif d == 2
            color = color2
        else
            error()
        end
    end
    plt, color
end

function makie_faces_impl(plt, color)
    d = 2
    #plt = shrink(plt,scale=0.995)
    plt, color = setup_colors_impl(plt, color, d)
    mesh = plt.mesh
    D = GalerkinToolkit.num_ambient_dims(mesh)
    nnodes = GalerkinToolkit.num_nodes(mesh)
    vert = zeros(Float64, nnodes, D)
    node_to_x = GalerkinToolkit.node_coordinates(mesh)
    for (node, x) in enumerate(node_to_x)
        vert[node, :] = x
    end
    nfaces = GalerkinToolkit.num_faces(mesh, d)
    face_to_nodes = GalerkinToolkit.face_nodes(mesh, d)
    conn = zeros(Int32, nfaces, 3)
    for (face, nodes) in enumerate(face_to_nodes)
        conn[face, :] = nodes
    end
    if nfaces == 0
        vert = ones(Float64, 3, D)
        conn = zeros(Int32, 1, 3)
        conn[:] = 1:3
        color = :pink
    end
    (; vert, conn, color)
end

Makie.@recipe(Makie2d1d) do scene
    dt = Makie.default_theme(scene, Makie1d)
    dt
end

Makie.preferred_axis_type(plot::Makie2d1d) = Makie.LScene

function Makie.plot!(sc::Makie2d1d{<:Tuple{<:GalerkinToolkit.Plot}})
    plt = sc[1]
    valid_attributes = Makie.shared_attributes(sc, Makie1d)
    colorrange = valid_attributes[:colorrange]
    color = valid_attributes[:color]
    colorrange = Makie.lift(setup_colorrange_impl, plt, color, colorrange)
    valid_attributes[:colorrange] = colorrange
    plt2 = Makie.lift(makie_face_edges_impl, plt)
    makie1d!(sc, valid_attributes, plt2)
end

function Makie.plot!(sc::Makie2d1d{<:Tuple{<:GalerkinToolkit.PPlot}})
    valid_attributes = Makie.shared_attributes(sc, Makie2d1d)
    pplt = sc[1]
    plts = makie_pplot_setup(pplt)
    if plts !== nothing
        foreach(plts) do plt
            makie2d1d!(sc, valid_attributes, plt)
        end
    end
end

function makie_face_edges_impl(plt)
    # TODO maybe already complexified
    D = 2
    d = 1
    mesh2, = GalerkinToolkit.complexify(GalerkinToolkit.restrict_to_dim(plt.mesh, D))
    topo = GalerkinToolkit.topology(mesh2)
    edge_to_faces = GalerkinToolkit.face_incidence(topo, d, D)
    K = Any
    facedata = Dict{String,K}()
    for (k, v) in GalerkinToolkit.face_data(plt, D; merge_dims = true)
        facedata[k] = map(faces -> sum(v[faces]) / length(faces), edge_to_faces)
    end
    mesh3 = GalerkinToolkit.restrict_to_dim(mesh2, d)
    fd = map(0:d) do i
        if i == d
            facedata
        else
            typeof(facedata)()
        end
    end
    plt3 = GalerkinToolkit.Plot(mesh3, fd, GalerkinToolkit.node_data(plt))
    plt3
end

Makie.@recipe(Makie1d) do scene
    dt = Makie.default_theme(scene, Makie.LineSegments)
    dt
end

Makie.preferred_axis_type(plot::Makie1d) = Makie.LScene

function Makie.plot!(sc::Makie1d{<:Tuple{<:GalerkinToolkit.Plot}})
    plt = sc[1]
    valid_attributes = Makie.shared_attributes(sc, Makie.LineSegments)
    color = valid_attributes[:color]
    colorrange = valid_attributes[:colorrange]
    colorrange = Makie.lift(setup_colorrange_impl, plt, color, colorrange)
    args = Makie.lift(makie_edges_impl, plt, color)
    p = Makie.lift(a -> a.p, args)
    color = Makie.lift(a -> a.color, args)
    valid_attributes[:color] = color
    valid_attributes[:colorrange] = colorrange
    Makie.linesegments!(sc, valid_attributes, p)
end

function Makie.plot!(sc::Makie1d{<:Tuple{<:GalerkinToolkit.PPlot}})
    valid_attributes = Makie.shared_attributes(sc, Makie1d)
    pplt = sc[1]
    plts = makie_pplot_setup(pplt)
    if plts !== nothing
        foreach(plts) do plt
            makie1d!(sc, valid_attributes, plt)
        end
    end
end

function makie_edges_impl(plt, color)
    d = 1
    plt, color = setup_colors_impl(plt, color, d)
    mesh = plt.mesh
    nedges = GalerkinToolkit.num_faces(mesh, d)
    node_to_x = GalerkinToolkit.node_coordinates(mesh)
    edge_to_nodes = GalerkinToolkit.face_nodes(mesh, d)
    T = eltype(eltype(node_to_x))
    S = GalerkinToolkit.num_ambient_dims(mesh)
    V = eltype(node_to_x)
    p = Vector{V}(undef, 2 * nedges)
    k = 0
    for edge = 1:nedges
        nodes = edge_to_nodes[edge]
        for node in nodes
            k += 1
            p[k] = node_to_x[node]
        end
    end
    (; p, color)
end

Makie.@recipe(Makie0d) do scene
    dt = Makie.default_theme(scene, Makie.Scatter)
    dt
end

Makie.preferred_axis_type(plot::Makie0d) = Makie.LScene

function Makie.plot!(sc::Makie0d{<:Tuple{<:GalerkinToolkit.Plot}})
    plt = sc[1]
    valid_attributes = Makie.shared_attributes(sc, Makie.Scatter)
    color = valid_attributes[:color]
    colorrange = valid_attributes[:colorrange]
    colorrange = Makie.lift(setup_colorrange_impl, plt, color, colorrange)
    args = Makie.lift(makie_vertices_impl, plt, color)
    p = Makie.lift(a -> a.p, args)
    color = Makie.lift(a -> a.color, args)
    valid_attributes[:color] = color
    valid_attributes[:colorrange] = colorrange
    Makie.scatter!(sc, valid_attributes, p)
end

function Makie.plot!(sc::Makie0d{<:Tuple{<:GalerkinToolkit.PPlot}})
    valid_attributes = Makie.shared_attributes(sc, Makie0d)
    pplt = sc[1]
    plts = makie_pplot_setup(pplt)
    if plts !== nothing
        foreach(plts) do plt
            makie0d!(sc, valid_attributes, plt)
        end
    end
end

function makie_vertices_impl(plt, color)
    d = 0
    plt, color = setup_colors_impl(plt, color, d)
    mesh = plt.mesh
    nedges = GalerkinToolkit.num_faces(mesh, d)
    node_to_x = GalerkinToolkit.node_coordinates(mesh)
    edge_to_nodes = GalerkinToolkit.face_nodes(mesh, d)
    T = eltype(eltype(node_to_x))
    S = GalerkinToolkit.num_ambient_dims(mesh)
    V = eltype(node_to_x)
    p = Vector{V}(undef, nedges)
    k = 0
    for edge = 1:nedges
        nodes = edge_to_nodes[edge]
        for node in nodes
            k += 1
            p[k] = node_to_x[node]
        end
    end
    (; p, color)
end

end # modules

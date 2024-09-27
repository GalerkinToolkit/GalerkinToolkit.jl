
struct PlotNew{A,B,C,D} <: AbstractType
    mesh::A
    face_data::B
    node_data::C
    cache::D
    function PlotNew(mesh,face_data,node_data,cache=nothing)
        A = typeof(mesh)
        B = typeof(face_data)
        C = typeof(node_data)
        D = typeof(cache)
        new{A,B,C,D}(mesh,face_data,node_data,cache)
    end
end

struct PPlot{A} <: AbstractType
    partition::A
end

function PartitionedArrays.partition(plt::PPlot)
    plt.partition
end

function face_data(plt::PlotNew,d;merge_dims=Val(false))
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

function face_color(plt::PlotNew,name,d)
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

function face_colorrange(plt::PlotNew,name)
    mesh = plt.mesh
    allcolors = face_data(plt;merge_dims=true)[name]
    minc = minimum(allcolors)
    maxc = maximum(allcolors)
    colorrange = (minc,maxc)
end

function node_color(plt::PlotNew,name)
    node_data(plt)[name]
end

function node_colorrange(plt::PlotNew,name)
    mesh = plt.mesh
    allcolors = node_data(plt)[name]
    minc = minimum(allcolors)
    maxc = maximum(allcolors)
    colorrange = (minc,maxc)
end

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
    #dict["__LOCAL_DFACE__"] = collect(1:ndfaces)
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
        PlotNew(mesh,face_data(mesh,ids),node_data(mesh,ids))
    end
    PPlot(plts)
end

function restrict_to_dim(mesh::AbstractMesh,d)
    chain = chain_from_arrays(
    node_coordinates(mesh),
    face_nodes(mesh,d),
    face_reference_id(mesh,d),
    reference_faces(mesh,d);
    periodic_nodes = periodic_nodes(mesh),
    physical_faces = physical_faces(mesh,d))
    mesh_from_chain(chain)
end

function restrict_to_dim(plt::PlotNew,d)
    mesh = restrict_to_dim(plt.mesh,d)
    dfacedata = face_data(plt,d;merge_dims=true)
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

function shrink(pplt::PPlot;scale=0.75)
    plts = map(partition(pplt)) do plt
        shrink(plt;scale)
    end
    PPlot(plts)
end

function plotnew(domain::AbstractDomain;kwargs...)
    mesh = GT.mesh(domain)
    d = GT.face_dim(domain)
    domface_to_face = GT.faces(domain)
    vismesh = GT.visualization_mesh(mesh,d,domface_to_face;kwargs...)
    node_data = Dict{String,Any}()
    #quad = PlotQuadrature(mesh,domain,vismesh)
    #foreach(pairs(fields)) do (sym,field)
    #    node_data[string(sym)] = plot_field(quad,field)
    #end
    vmesh,glue = vismesh
    cache = (;glue,domain)
    PlotNew(vmesh,face_data(vmesh),node_data,cache)
end

function plot(domain::AbstractDomain{<:PMesh};kwargs...)
    mesh = GT.mesh(domain)
    plts = map(partition(domain)) do mydom
        plt = plotnew(mydom;kwargs...)
    end
    PPlot(plts)
end

function plot!(field,plt::PlotNew;label)
    plot!(plt,field;label)
end

function plot!(plt::PlotNew,field;label)
    q = GT.coordinates(plt)
    f_q = field(q)
    term = GT.term(f_q)
    T = typeof(GT.prototype(f_q))
    plot_impl!(plt,term,label,T)
end

function plot_impl!(plt,term,label,::Type{T}) where T
    vmesh  = plt.mesh
    vglue = plt.cache.glue
    nnodes = GT.num_nodes(vmesh)
    data = zeros(T,nnodes)
    face_to_nodes = vglue.face_fine_nodes
    face = :face
    point = :point
    index = GT.index(;face,point)
    dict = index.dict
    dict[face_to_nodes] = :face_to_nodes
    deps = (face,point)
    expr1 = term(index) |> simplify
    v = GT.topological_sort(expr1,deps)
    expr = quote
        function filldata!(data,state)
            $(unpack_storage(dict,:state))
            $(v[1])
            for $face in 1:length(face_to_nodes)
                nodes = face_to_nodes[$face]
                $(v[2])
                for $point in 1:length(nodes)
                    data[nodes[$point]] = $(v[3])
                end
            end
        end
    end
    filldata! = eval(expr)
    invokelatest(filldata!,data,GT.storage(index))
    plt.node_data[label] = data
    plt
end

function plot!(plt::PPlot,field;label)
    q = GT.coordinates(plt)
    f_q = field(q)
    term = GT.term(f_q)
    T = typeof(GT.prototype(f_q))
    map(partition(plt),term) do myplt, myterm
        plot_impl!(myplt,myterm,label,T)
    end
    plt
end

domain(plt::PlotNew) = plt.cache.domain
visualization_mesh(plt::PlotNew) = (plt.mesh, plt.cache.glue)

function coordinates(plt::PlotNew)
    domain = plt |> GT.domain
    GT.coordinates(plt,domain)
end

function coordinates(plt::PlotNew,::ReferenceDomain)
    GT.reference_coordinates(plt)
end

function coordinates(plt::PlotNew,::PhysicalDomain)
    domain_phys = plt |> GT.domain
    domain_ref = domain_phys |> reference_domain
    phi = GT.domain_map(domain_ref,domain_phys)
    q = GT.reference_coordinates(plt)
    phi(q)
end

function reference_coordinates(plt::PlotNew)
    domain = GT.reference_domain(plt.domain)
    d = GT.face_dim(domain)
    domface_to_face = GT.faces(domain)
    mesh = GT.mesh(domain)
    vmesh, vglue = GT.visualization_mesh(plt)
    refid_to_snode_to_coords = vglue.reference_coordinates
    d = GT.num_dims(vmesh)
    face_to_refid = GT.face_reference_id(mesh,d)
    prototype = first(first(refid_to_snode_to_coords))
    GT.quantity(prototype,domain) do index
        domface = index.face
        point = index.point
        dict = index.dict
        domface_to_face_sym = get!(dict,domface_to_face,gensym("domface_to_face"))
        face_to_refid_sym = get!(dict,face_to_refid,gensym("face_to_refid"))
        refid_to_snode_to_coords_sym = get!(dict,refid_to_snode_to_coords,gensym("refid_to_snode_to_coords"))
        @term begin
            face = $domface_to_face_sym[$domface]
            coords = reference_value(
                $refid_to_snode_to_coords_sym,
                $face_to_refid_sym,
                $face)
            coords[$point]
        end
    end
end

function reference_coordinates(plt::PPlot)
    q = map(GT.reference_coordinates,partition(plt))
    term = map(GT.term,q)
    prototype = map(GT.prototype,q) |> PartitionedArrays.getany
    GT.quantity(term,prototype,plt.domain)
end

# VTK

function translate_vtk_data(v)
    v
end
function translate_vtk_data(v::AbstractVector{<:SVector{2}})
    z = zero(eltype(eltype(v)))
    map(vi->SVector((vi...,z)),v)
end

function WriteVTK.vtk_grid(filename,plt::PlotNew;kwargs...)
    mesh = plt.mesh
    D = num_dims(mesh)
    vtk = vtk_grid(filename,GT.vtk_args(mesh)...;kwargs...)
    for (k,v) in node_data(plt)
        vtk[k,WriteVTK.VTKPointData()] = translate_vtk_data(v)
    end
    for (k,v) in face_data(plt;merge_dims=true)
        vtk[k,WriteVTK.VTKCellData()] = translate_vtk_data(v)
    end
    cache = (;vtk,plt.cache...)
    PlotNew(plt.mesh,plt.face_data,plt.node_data,cache)
end

function WriteVTK.vtk_grid(f,filename,plt::Union{PlotNew,PPlot};kwargs...)
    vtk = vtk_grid(filename,plt)
    try
        f(vtk)
    finally
        files = close(vtk)
    end
    files
end

function WriteVTK.vtk_grid(filename,pplt::PPlot;kwargs...)
    plts = partition(pplt)
    parts = linear_indices(plts)
    nparts = length(parts)
    plts  = map(plts,parts) do plt,part
        mesh = plt.mesh
        vtk = pvtk_grid(filename,GT.vtk_args(mesh)...;part,nparts;kwargs...) do vtk
            for (k,v) in node_data(plt)
                vtk[k,WriteVTK.VTKPointData()] = translate_vtk_data(v)
            end
            for (k,v) in face_data(plt;merge_dims=true)
                vtk[k,WriteVTK.VTKCellData()] = translate_vtk_data(v)
            end
        end
        cache = (;vtk,plt.cache...)
        PlotNew(plt.mesh,plt.face_data,plt.node_data,cache)
    end
    PPlot(plts)
end

function WriteVTK.close(plt::PlotNew)
    WriteVTK.close(pplot.cache.vtk)
end

function WriteVTK.close(plt::PPlot)
    map(WriteVTK.close,pplot.partition)
end

struct PVD{A} <: AbstractType
    pvd::A
end

struct PPVD{A} <: AbstractType
    partition::A
end

function Base.setindex!(a::PVD,time,plt::PlotNew)
    a.pvd[time] = plt.cache.vtk
end

function Base.setindex!(a::PPVD,time,plt::PPlot)
    map_main(a.partition,plt.partition) do a,plt
        a[time] = plt
    end
end

function WriteVTK.paraview_collection(filename,pplt::PlotNew;kwargs...)
    WriteVTK.paraview_collection(filename;kwargs...)
end

function WriteVTK.paraview_collection(filename,pplt::PPlot;kwargs...)
    pvds = map_main(partition(pplt)) do _
        WriteVTK.paraview_collection(filename;kwargs...)
    end
    PPVD(pvds)
end

function WriteVTK.close(ppvd::PPVD)
    map_main(WriteVTK.close,ppvd.partition)
end

function WriteVTK.paraview_collection(f,filename,plt::Union{PlotNew,PPlot};kwargs...)
    vtk = WriteVTK.paraview_collection(filename,plt;kwargs...)
    try
        f(vtk)
    finally
        files = close(vtk)
    end
    files
end

function WriteVTK.vtk_grid(filename,mesh::Union{AbstractMesh,PMesh};plot_params=(;),vtk_grid_params=(;))
    plt = plotnew(mesh;plot_params...)
    WriteVTK.vtk_grid(filename,plt;vtk_grid_params...)
end

function WriteVTK.vtk_grid(filename,mesh::AbstractDomain;plot_params=(;),vtk_grid_params=(;))
    plt = plotnew(mesh;plot_params...)
    WriteVTK.vtk_grid(filename,plt;vtk_grid_params...)
end

# Makie

Makie.@recipe(MakiePlot) do scene
    t1 = Makie.Theme(
        dim = nothing,
        shrink = false,
        edgecolor = nothing,
        color      = :lightblue,
        colormap   = :bluesreds,
        #shading    = Makie.NoShading,
        cycle      = nothing,
       )
    t2 = Makie.default_theme(scene, Makie.Mesh)
    merge(t1,t2)
end

Makie.preferred_axis_type(plot::MakiePlot) = Makie.LScene

function Makie.plot!(sc::MakiePlot{<:Tuple{<:PlotNew}})
    plt = sc[1]
    # TODO these are not reactive
    if sc[:shrink][] != false
        scale = sc[:shrink]
        plt = Makie.lift((a,b)->shrink(a;scale=b),plt,scale)
    end
    dim = sc[:dim][]
    d = num_dims(plt[].mesh)
    edgecolor = sc[:edgecolor][]
    if dim === nothing
        cmp = (x,y) -> x >= y
    else
        cmp = (x,y) -> y==dim
    end
    if cmp(d,0)
        valid_attributes = Makie.shared_attributes(sc, Makie0d)
        makie0d!(sc,valid_attributes,plt)
    end
    if cmp(d,1)
        valid_attributes = Makie.shared_attributes(sc, Makie1d)
        makie1d!(sc,valid_attributes,plt)
    end
    if cmp(d,2)
        valid_attributes = Makie.shared_attributes(sc, Makie2d)
        makie2d!(sc,valid_attributes,plt)
        if edgecolor !== nothing
            valid_attributes = Makie.shared_attributes(sc, Makie2d1d)
            valid_attributes[:color] = sc[:edgecolor]
            makie2d1d!(sc,valid_attributes,plt)
        end
    end
    if cmp(d,3)
        valid_attributes = Makie.shared_attributes(sc, Makie3d)
        makie3d!(sc,valid_attributes,plt)
        if edgecolor !== nothing
            valid_attributes = Makie.shared_attributes(sc, Makie3d1d)
            valid_attributes[:color] = sc[:edgecolor]
            makie3d1d!(sc,valid_attributes,plt)
        end
    end
end

Makie.plottype(::PlotNew) = MakiePlot

Makie.plottype(::AbstractMesh) = MakiePlot

function Makie.convert_arguments(::Type{<:MakiePlot},mesh::AbstractMesh)
    plt = plot(mesh)
    (plt,)
end

Makie.plottype(::AbstractDomain) = MakiePlot

function Makie.plot!(sc::MakiePlot{<:Tuple{<:AbstractDomain}})
    dom = sc[1]
    valid_attributes = Makie.shared_attributes(sc, MakiePlot)
    color = valid_attributes[:color]
    args = Makie.lift(dom,color) do dom,color
        if isa(color,AbstractQuantity)
            sym = gensym()
            plt = plotnew(dom;fields=Dict(sym=>color))
            color = NodeData(string(sym))
        else
            plt = plotnew(dom)
        end
        (;plt,color)
    end
    plt = Makie.lift(i->i.plt,args)
    color = Makie.lift(i->i.color,args)
    valid_attributes[:color] = color
    makieplot!(sc,valid_attributes,plt)
end

# TODO not sure about this
# what if u is not scalar-valued?
#Makie.plottype(::AbstractQuantity) = MakiePlot
#
#function Makie.plot!(sc::MakiePlot{<:Tuple{<:AbstractQuantity}})
#    u = sc[1]
#    valid_attributes = Makie.shared_attributes(sc, MakiePlot)
#    dom = Makie.lift(domain,u)
#    valid_attributes[:color] = u
#    makieplot!(sc,valid_attributes,dom)
#end

## TODO not sure about this
# function vector_of_observables(a)
#    #TODO remove trick for DebugArray
#    function start(v)
#        first(v)
#    end
#    function start(v::DebugArray)
#        first(v.items)
#    end
#    function rest(v)
#        v[2:end]
#    end
#    function rest(v::DebugArray)
#        v.items[2:end]
#    end
#    if length(a[]) == 1
#        b = Makie.lift(first,a)
#        return [b,]
#    else
#        b = Makie.lift(start,a)
#        c = Makie.lift(rest,a)
#        return [b,vector_of_observables(c)...]
#    end
#end

Makie.@recipe(Makie3d) do scene
    dt = Makie.default_theme(scene, Makie.Mesh)
    dt
end

Makie.preferred_axis_type(plot::Makie3d) = Makie.LScene

function Makie.plot!(sc::Makie3d{<:Tuple{<:PlotNew}})
    plt = sc[1]
    valid_attributes = Makie.shared_attributes(sc, Makie2d)
    colorrange = valid_attributes[:colorrange]
    color = valid_attributes[:color]
    colorrange = Makie.lift(setup_colorrange_impl,plt,color,colorrange)
    valid_attributes[:colorrange] = colorrange
    plt2 = Makie.lift(makie_volumes_impl,plt)
    makie2d!(sc,valid_attributes,plt2)
end

# TODO not sure about this
#function Makie.plot!(sc::Makie3d{<:Tuple{<:PPlot}})
#    valid_attributes = Makie.shared_attributes(sc,Makie3d)
#    pplt = sc[1]
#    plts = vector_of_observables(Makie.lift(partition,pplt))
#    foreach(plts) do plt
#        makie3d!(sc,valid_attributes,plt)
#    end
#end

function makie_volumes_impl(plt::PlotNew)
    @assert num_dims(plt.mesh) == 3
    D=3
    d=2
    mesh, = complexify(restrict_to_dim(plt.mesh,D))
    topo = topology(mesh)
    face_to_cells = face_incidence(topo,d,D)
    face_isboundary = map(cells->length(cells)==1,face_to_cells)
    mesh2 = restrict_to_dim(mesh,d)
    newnodes = 1:num_nodes(mesh2)
    newfaces = [ Int[] for _ in 0:d ]
    newfaces[end] = findall(face_isboundary)
    mesh3 = restrict(mesh2,newnodes,newfaces)
    face_to_cell = map(first,face_to_cells)
    newface_to_cell = face_to_cell[newfaces[end]]
    celldata = copy(face_data(plt,D;merge_dims=true))
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
    valid_attributes = Makie.shared_attributes(sc, Makie2d1d)
    colorrange = valid_attributes[:colorrange]
    color = valid_attributes[:color]
    colorrange = Makie.lift(setup_colorrange_impl,plt,color,colorrange)
    valid_attributes[:colorrange] = colorrange
    plt2 = Makie.lift(makie_volumes_impl,plt)
    makie2d1d!(sc,valid_attributes,plt2)
end

Makie.@recipe(Makie2d) do scene
    dt = Makie.default_theme(scene, Makie.Mesh)
    dt
end

Makie.preferred_axis_type(plot::Makie2d) = Makie.LScene

function Makie.plot!(sc::Makie2d{<:Tuple{<:PlotNew}})
    plt = sc[1]
    valid_attributes = Makie.shared_attributes(sc, Makie.Mesh)
    color = valid_attributes[:color]
    colorrange = valid_attributes[:colorrange]
    colorrange = Makie.lift(setup_colorrange_impl,plt,color,colorrange)
    args = Makie.lift(makie_faces_impl,plt,color)
    vert = Makie.lift(a->a.vert,args)
    conn = Makie.lift(a->a.conn,args)
    color = Makie.lift(a->a.color,args)
    valid_attributes[:color] = color
    valid_attributes[:colorrange] = colorrange
    Makie.mesh!(sc,valid_attributes,vert,conn)
end

# TODO Not sure about this
#function Makie.plot!(sc::Makie2d{<:Tuple{<:PPlot}})
#    valid_attributes = Makie.shared_attributes(sc,Makie2d)
#    pplt = sc[1]
#    plts = vector_of_observables(Makie.lift(partition,pplt))
#    foreach(plts) do plt
#        makie2d!(sc,valid_attributes,plt)
#    end
#end

function setup_colorrange_impl(plt,color,colorrange)
    if colorrange != Makie.Automatic()
        return colorrange
    end
    if isa(color,FaceData)
        colorrange = face_colorrange(plt,color.name)
    end
    if isa(color,NodeData)
        colorrange = node_colorrange(plt,color.name)
    end
    colorrange
end

function setup_colors_impl(plt,color,d)
    if isa(color,FaceData)
        if d == 2
            plt = shrink(plt;scale=1)
        end
        color2 = face_color(plt,color.name,d)
    end
    if isa(color,NodeData)
        color2 = node_color(plt,color.name)
    end
    if isa(color,NodeData) || isa(color,FaceData)
        face_to_nodes = face_nodes(plt.mesh,d)
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
    d = 2
    #plt = shrink(plt,scale=0.995)
    plt,color = setup_colors_impl(plt,color,d)
    mesh = plt.mesh
    D = num_ambient_dims(mesh)
    nnodes = num_nodes(mesh)
    vert = zeros(Float64,nnodes,D)
    node_to_x = node_coordinates(mesh)
    for (node,x) in enumerate(node_to_x)
        vert[node,:] = x
    end
    nfaces = num_faces(mesh,d)
    face_to_nodes = face_nodes(mesh,d)
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
    (;vert,conn,color)
end

Makie.@recipe(Makie2d1d) do scene
    dt = Makie.default_theme(scene, Makie1d)
    dt
end

Makie.preferred_axis_type(plot::Makie2d1d) = Makie.LScene

function Makie.plot!(sc::Makie2d1d{<:Tuple{<:PlotNew}})
    plt = sc[1]
    valid_attributes = Makie.shared_attributes(sc, Makie1d)
    colorrange = valid_attributes[:colorrange]
    color = valid_attributes[:color]
    colorrange = Makie.lift(setup_colorrange_impl,plt,color,colorrange)
    valid_attributes[:colorrange] = colorrange
    plt2 = Makie.lift(makie_face_edges_impl,plt)
    makie1d!(sc,valid_attributes,plt2)
end

function makie_face_edges_impl(plt)
    # TODO maybe already complexified
    D=2
    d=1
    mesh2, = complexify(restrict_to_dim(plt.mesh,D))
    topo = topology(mesh2)
    edge_to_faces = face_incidence(topo,d,D)
    K = Any
    facedata = Dict{String,K}()
    for (k,v) in face_data(plt,D;merge_dims=true)
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
    valid_attributes = Makie.shared_attributes(sc, Makie.LineSegments)
    color = valid_attributes[:color]
    colorrange = valid_attributes[:colorrange]
    colorrange = Makie.lift(setup_colorrange_impl,plt,color,colorrange)
    args = Makie.lift(makie_edges_impl,plt,color)
    p = Makie.lift(a->a.p,args)
    color = Makie.lift(a->a.color,args)
    valid_attributes[:color] = color
    valid_attributes[:colorrange] = colorrange
    Makie.linesegments!(sc,valid_attributes,p)
end

function makie_edges_impl(plt,color)
    d = 1
    plt,color = setup_colors_impl(plt,color,d)
    mesh = plt.mesh
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
    (;p,color)
end

Makie.@recipe(Makie0d) do scene
    dt = Makie.default_theme(scene, Makie.Scatter)
    dt
end

Makie.preferred_axis_type(plot::Makie0d) = Makie.LScene

function Makie.plot!(sc::Makie0d{<:Tuple{<:PlotNew}})
    plt = sc[1]
    valid_attributes = Makie.shared_attributes(sc, Makie.Scatter)
    color = valid_attributes[:color]
    colorrange = valid_attributes[:colorrange]
    colorrange = Makie.lift(setup_colorrange_impl,plt,color,colorrange)
    args = Makie.lift(makie_vertices_impl,plt,color)
    p = Makie.lift(a->a.p,args)
    color = Makie.lift(a->a.color,args)
    valid_attributes[:color] = color
    valid_attributes[:colorrange] = colorrange
    Makie.scatter!(sc,valid_attributes,p)
end

function makie_vertices_impl(plt,color)
    d = 0
    plt,color = setup_colors_impl(plt,color,d)
    mesh = plt.mesh
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
    (;p,color)
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
#struct FaceData
#    name::String
#end
#
#struct NodeData
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
#function makie_mesh_color(plt,fc::FaceData)
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
#function makie_linesegments_color(plt,fc::FaceData)
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
#function makie_scatter_color(plt,fc::FaceData)
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
#function makie_mesh_color_2(plt,fc::FaceData)
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




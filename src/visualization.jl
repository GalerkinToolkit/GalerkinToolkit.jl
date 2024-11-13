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
        Plot(mesh,face_data(mesh,ids),node_data(mesh,ids))
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
    Plot(vmesh,face_data(vmesh),node_data,cache)
end

function plot(domain::AbstractDomain{<:PMesh};kwargs...)
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
    T = typeof(GT.prototype(f_q))
    plot_impl!(plt,term,label,T)
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

domain(plt::Union{Plot,PPlot}) = plt.cache.domain
visualization_mesh(plt::Plot) = (plt.mesh, plt.cache.glue)

function coordinates(plt::Union{Plot,PPlot})
    domain = plt |> GT.domain
    GT.coordinates(plt,domain)
end

function coordinates(plt::Union{Plot,PPlot},::ReferenceDomain)
    GT.reference_coordinates(plt)
end

function coordinates(plt::Union{Plot,PPlot},::PhysicalDomain)
    domain_phys = plt |> GT.domain
    domain_ref = domain_phys |> reference_domain
    phi = GT.domain_map(domain_ref,domain_phys)
    q = GT.reference_coordinates(plt)
    phi(q)
end

function reference_coordinates(plt::Plot)
    domain = GT.reference_domain(plt.cache.domain)
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

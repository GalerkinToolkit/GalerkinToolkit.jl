module GalerkinToolkitWriteVTKExt

using GalerkinToolkit
import GalerkinToolkit: plot, plot!, translate_vtk_data, translate_vtk_data!, vtk_close_impl, vtk_close_impl!
import GalerkinToolkit: vtk_points, vtk_points!, vtk_cells, vtk_cells!, vtk_args, vtk_args!
import GalerkinToolkit: vtk_group_faces, vtk_group_faces! #, vtk_physical_nodes, vtk_physical_nodes!
import GalerkinToolkit: vtk_mesh_cell, vtk_mesh_cell!
using PartitionedArrays
using StaticArrays
using WriteVTK

# Visualization

function translate_vtk_data(v)
    v
end

function translate_vtk_data(v::AbstractVector{<:SVector{2}})
    z = zero(eltype(eltype(v)))
    map(vi->SVector((vi...,z)),v)
end

function plot!(field,plt::GalerkinToolkit.VTKPlot{<:GalerkinToolkit.Plot};label)
    plot!(plt.plot,field;label)
    v = plt.plot.node_data[label]
    plt.vtk[label,WriteVTK.VTKPointData()] = translate_vtk_data(v)
    plt
end

function plot!(field,plt::GalerkinToolkit.VTKPlot{<:GalerkinToolkit.PPlot};label)
    plot!(plt.plot,field;label)
    foreach(plt.plot.partition,plt.vtk) do myplt, myvtk
        v = myplt.node_data[label]
        myvtk[label,WriteVTK.VTKPointData()] = translate_vtk_data(v)
    end
    plt
end

function plot!(plt::GalerkinToolkit.VTKPlot,field;label)
    plot!(field,plt;label)
end

function WriteVTK.vtk_grid(filename,plt::GalerkinToolkit.Plot;kwargs...)
    mesh = plt.mesh
    D = GalerkinToolkit.num_dims(mesh)
    vtk = vtk_grid(filename,GalerkinToolkit.vtk_args(mesh)...;kwargs...)
    for (k,v) in GalerkinToolkit.node_data(plt)
        vtk[k,WriteVTK.VTKPointData()] = translate_vtk_data(v)
    end
    for (k,v) in GalerkinToolkit.face_data(plt;merge_dims=true)
        vtk[k,WriteVTK.VTKCellData()] = translate_vtk_data(v)
    end
    GalerkinToolkit.VTKPlot(plt,vtk)
end

function WriteVTK.vtk_grid(f::Function,filename,plt::Union{GalerkinToolkit.Plot,GalerkinToolkit.PPlot};kwargs...)
    vtk = vtk_grid(filename,plt)
    files = nothing
    try
        f(vtk)
    finally
        files = close(vtk)
    end
    files
end

function WriteVTK.vtk_grid(filename,pplt::GalerkinToolkit.PPlot;kwargs...)
    plts = partition(pplt)
    parts = linear_indices(plts)
    nparts = length(parts)
    vtks  = map(plts,parts) do plt,part
        mesh = plt.mesh
        vtk = pvtk_grid(filename,GalerkinToolkit.vtk_args(mesh)...;part,nparts,kwargs...)
        for (k,v) in GalerkinToolkit.node_data(plt)
            vtk[k,WriteVTK.VTKPointData()] = translate_vtk_data(v)
        end
        for (k,v) in GalerkinToolkit.face_data(plt;merge_dims=true)
            vtk[k,WriteVTK.VTKCellData()] = translate_vtk_data(v)
        end
        vtk
    end
    plts = GalerkinToolkit.PPlot(plts,pplt.cache)
    GalerkinToolkit.VTKPlot(plts,vtks)
end

function WriteVTK.pvtk_grid(filename::AbstractString,pplt::GalerkinToolkit.Plot;kwargs...)
    plts = GalerkinToolkit.subdomains(pplt)
    parts = linear_indices(plts)
    nparts = length(parts)
    vtks  = map(plts,parts) do plt,part
        mesh = plt.mesh
        vtk = pvtk_grid(filename,GalerkinToolkit.vtk_args(mesh)...;part,nparts,kwargs...)
        for (k,v) in GalerkinToolkit.node_data(plt)
            vtk[k,WriteVTK.VTKPointData()] = translate_vtk_data(v)
        end
        for (k,v) in GalerkinToolkit.face_data(plt;merge_dims=true)
            vtk[k,WriteVTK.VTKCellData()] = translate_vtk_data(v)
        end
        vtk
    end
    GalerkinToolkit.PVTKPlot(vtks)
end

function WriteVTK.pvtk_grid(f::Function,filename,plt::GalerkinToolkit.Plot;kwargs...)
    vtk = pvtk_grid(filename,plt)
    files = nothing
    try
        f(vtk)
    finally
        files = close(vtk)
    end
    files
end

function WriteVTK.close(plt::GalerkinToolkit.VTKPlot)
    vtk_close_impl(plt.plot,plt.vtk)
end

function vtk_close_impl(plt::GalerkinToolkit.Plot,vtk)
    WriteVTK.close(vtk)
end

function vtk_close_impl(plt::GalerkinToolkit.PPlot,vtks)
    map(WriteVTK.close,vtks)
end

function WriteVTK.close(plt::GalerkinToolkit.PVTKPlot)
    map(WriteVTK.close,plt.vtks)
end

function WriteVTK.vtk_grid(filename,mesh::Union{GalerkinToolkit.AbstractMesh,GalerkinToolkit.PMesh};plot_params=(;),vtk_grid_params...)
    plt = plot(mesh;plot_params...)
    WriteVTK.vtk_grid(filename,plt;vtk_grid_params...)
end

function WriteVTK.vtk_grid(filename,dom::GalerkinToolkit.AbstractDomain;plot_params=(;),vtk_grid_params...)
    plt = plot(dom;plot_params...)
    WriteVTK.vtk_grid(filename,plt;vtk_grid_params...)
end

function WriteVTK.vtk_grid(f::Function,filename,mesh::Union{GalerkinToolkit.AbstractMesh,GalerkinToolkit.PMesh};plot_params=(;),vtk_grid_params...)
    plt = plot(mesh;plot_params...)
    WriteVTK.vtk_grid(f,filename,plt;vtk_grid_params...)
end

function WriteVTK.vtk_grid(f::Function,filename,dom::GalerkinToolkit.AbstractDomain;plot_params=(;),vtk_grid_params...)
    plt = plot(dom;plot_params...)
    WriteVTK.vtk_grid(f,filename,plt;vtk_grid_params...)
end

function WriteVTK.pvtk_grid(filename,mesh::GalerkinToolkit.AbstractMesh;plot_params=(;),vtk_grid_params...)
    plt = plot(mesh;plot_params...)
    WriteVTK.pvtk_grid(filename,plt;vtk_grid_params...)
end

function WriteVTK.pvtk_grid(f::Function,filename,mesh::GalerkinToolkit.AbstractMesh;plot_params=(;),vtk_grid_params...)
    plt = plot(mesh;plot_params...)
    WriteVTK.pvtk_grid(f,filename,plt;vtk_grid_params...)
end

function WriteVTK.close(pvd::GalerkinToolkit.PVD)
    WriteVTK.close(pvd.pvd)
end

function WriteVTK.close(ppvd::GalerkinToolkit.PPVD)
    map_main(WriteVTK.close,ppvd.partition)
end

function Base.setindex!(a::GalerkinToolkit.PVD,plt::GalerkinToolkit.VTKPlot,time)
    a.pvd[time] = plt.vtk
end

function Base.setindex!(a::GalerkinToolkit.PVD,plt::GalerkinToolkit.PVTKPlot,time)
    map_main(plt.vtks) do vtk
        a.pvd[time] = vtk
    end
end

function Base.setindex!(a::GalerkinToolkit.PPVD,plt::GalerkinToolkit.VTKPlot,time)
    map_main(a.partition,plt.vtk) do a,vtk
        a[time] = vtk
    end
end

function WriteVTK.paraview_collection(filename,pplt::GalerkinToolkit.Plot;kwargs...)
    WriteVTK.paraview_collection(filename;kwargs...) |> GalerkinToolkit.PVD
end

function WriteVTK.paraview_collection(filename,pplt::GalerkinToolkit.PPlot;kwargs...)
    pvds = map_main(partition(pplt)) do _
        WriteVTK.paraview_collection(filename;kwargs...)
    end
    GalerkinToolkit.PPVD(pvds)
end

function WriteVTK.paraview_collection(f::Function,filename,plt::Union{GalerkinToolkit.Plot,GalerkinToolkit.PPlot};kwargs...)
    vtk = WriteVTK.paraview_collection(filename,plt;kwargs...)
    files = nothing
    try
        f(vtk)
    finally
        files = close(vtk)
    end
    files
end

function WriteVTK.paraview_collection(filename,mesh::Union{GalerkinToolkit.AbstractMesh,GalerkinToolkit.PMesh};plot_params=(;),vtk_grid_params...)
    plt = GalerkinToolkit.plot(mesh;plot_params...)
    WriteVTK.paraview_collection(filename,plt;vtk_grid_params...)
end

function WriteVTK.paraview_collection(filename,dom::GalerkinToolkit.AbstractDomain;plot_params=(;),vtk_grid_params...)
    plt = GalerkinToolkit.plot(dom;plot_params...)
    WriteVTK.paraview_collection(filename,plt;vtk_grid_params...)
end

function WriteVTK.paraview_collection(f::Function,filename,mesh::Union{GalerkinToolkit.AbstractMesh,GalerkinToolkit.PMesh};plot_params=(;),vtk_grid_params...)
    plt = GalerkinToolkit.plot(mesh;plot_params...)
    WriteVTK.paraview_collection(f,filename,plt;vtk_grid_params...)
end

function WriteVTK.paraview_collection(f::Function,filename,dom::GalerkinToolkit.AbstractDomain;plot_params=(;),vtk_grid_params...)
    plt = GalerkinToolkit.plot(dom;plot_params...)
    WriteVTK.paraview_collection(f,filename,plt;vtk_grid_params...)
end

# Mesh

function vtk_points(mesh)
    function barrirer(coords)
        nnodes = length(coords)
        points = zeros(3,nnodes)
        for node in 1:nnodes
            coord = coords[node]
            for i in 1:length(coord)
                points[i,node] = coord[i]
            end
        end
        points
    end
    coords = GalerkinToolkit.node_coordinates(mesh)
    barrirer(coords)
end

function vtk_cells(mesh,d)
    function barrirer(face_to_refid,face_to_nodes,refid_mesh_cell)
        cells = map(face_to_refid,face_to_nodes) do refid, nodes
            mesh_cell = refid_mesh_cell[refid]
            if mesh_cell === nothing
                msg = """
                Not enough information to visualize this mesh via vtk:
                vtk_mesh_cell returns nothing for the reference face in position $refid in dimension $d.
                """
                error(msg)
            end
            mesh_cell(nodes)
        end
        cells
    end
    face_to_nodes = GalerkinToolkit.face_nodes(mesh,d)
    face_to_refid = GalerkinToolkit.face_reference_id(mesh,d)
    refid_refface = GalerkinToolkit.reference_spaces(mesh,d)
    refid_mesh_cell = map(vtk_mesh_cell,refid_refface)
    barrirer(face_to_refid,face_to_nodes,refid_mesh_cell)
end

"""
    args = vtk_args(mesh[,d])

Return the arguments `args` to be passed in final position
to functions like `WriteVTK.vtk_grid`.
"""
function vtk_args(mesh,d)
    points = vtk_points(mesh)
    cells = vtk_cells(mesh,d)
    points, cells
end

function vtk_args(mesh)
    points = vtk_points(mesh)
    D = GalerkinToolkit.num_dims(mesh)
    allcells = [vtk_cells(mesh,d) for d in 0:D if GalerkinToolkit.num_faces(mesh,d) != 0]
    cells = reduce(vcat,allcells)
    points, cells
end

"""
"""
function vtk_group_faces!(vtk,mesh,d;physical_faces=GalerkinToolkit.physical_faces(mesh,d))
    ndfaces = GalerkinToolkit.num_faces(mesh,d)
    for group in physical_faces
        name,faces = group
        face_mask = zeros(Int,ndfaces)
        face_mask[faces] .= 1
        vtk[name,WriteVTK.VTKCellData()] = face_mask
    end
    vtk
end

function vtk_group_faces!(vtk,mesh;physical_faces=GalerkinToolkit.physical_faces(mesh))
    nfaces = sum(GalerkinToolkit.num_faces(mesh))
    offsets = GalerkinToolkit.face_offset(mesh)
    D = GalerkinToolkit.num_dims(mesh)
    data = Dict{String,Vector{Int}}()
    for d in 0:D
        for group in physical_faces[d+1]
            name, = group
            if !haskey(data,name)
                face_mask = zeros(Int,nfaces)
                data[name] = face_mask
            end
        end
    end
    for d in 0:D
        for group in physical_faces[d+1]
            offset = offsets[d+1]
            name,faces = group
            face_mask = data[name]
            face_mask[faces.+offset] .= 1
        end
    end
    for (name,face_mask) in data
        vtk[name,WriteVTK.VTKCellData()] = face_mask
    end
    vtk
end

#function vtk_physical_nodes!(vtk,mesh,d;physical_nodes=GalerkinToolkit.physical_nodes(mesh,d))
#    nnodes = GalerkinToolkit.num_nodes(mesh)
#    for group in physical_nodes
#        name,nodes = group
#        nodes_mask = zeros(Int,nnodes)
#        nodes_mask[nodes] .= 1
#        vtk[name,WriteVTK.VTKPointData()] = nodes_mask
#    end
#    vtk
#end
#
#function vtk_physical_nodes!(vtk,mesh;physical_nodes=GalerkinToolkit.physical_nodes(mesh))
#    D = GalerkinToolkit.num_dims(mesh)
#    for d in 0:D
#        vtk_physical_nodes!(vtk,mesh,d,physical_nodes=physical_nodes[d+1])
#    end
#    vtk
#end

"""
"""
function vtk_mesh_cell(ref_face)
    geom = GalerkinToolkit.domain(ref_face)
    d = GalerkinToolkit.num_dims(geom)
    nnodes = GalerkinToolkit.num_nodes(ref_face)
    lib_to_user = GalerkinToolkit.lib_to_user_nodes(ref_face)
    if d == 0 && nnodes == 1
        cell_type = WriteVTK.VTKCellTypes.VTK_VERTEX
        vtk_to_lib = [1]
    elseif d == 1 && (GalerkinToolkit.is_simplex(geom) || GalerkinToolkit.is_n_cube(geom)) && nnodes == 2
        cell_type = WriteVTK.VTKCellTypes.VTK_LINE
        vtk_to_lib = [1,2]
    elseif d == 1 && (GalerkinToolkit.is_simplex(geom) || GalerkinToolkit.is_n_cube(geom)) && nnodes == 3
        cell_type = WriteVTK.VTKCellTypes.VTK_QUADRATIC_EDGE
        vtk_to_lib = [1,3,2]
    elseif d == 2 && GalerkinToolkit.is_n_cube(geom) && nnodes == 4
        cell_type = WriteVTK.VTKCellTypes.VTK_QUAD
        vtk_to_lib = [1,2,4,3]
    elseif d == 2 && GalerkinToolkit.is_simplex(geom) && nnodes == 3
        cell_type = WriteVTK.VTKCellTypes.VTK_TRIANGLE
        vtk_to_lib = [1,2,3]
    elseif d == 2 && GalerkinToolkit.is_simplex(geom) && nnodes == 6
        cell_type = WriteVTK.VTKCellTypes.VTK_QUADRATIC_TRIANGLE
        vtk_to_lib = [1,3,6,2,5,4]
    elseif d == 3 && GalerkinToolkit.is_n_cube(geom) && nnodes == 8
        cell_type = WriteVTK.VTKCellTypes.VTK_HEXAHEDRON
        vtk_to_lib = [1,2,4,3,5,6,8,7]
    elseif d == 3 && GalerkinToolkit.is_simplex(geom) && nnodes == 4
        cell_type = WriteVTK.VTKCellTypes.VTK_TETRA
        vtk_to_lib = [1,2,3,4]
    elseif d == 3 && GalerkinToolkit.is_simplex(geom) && nnodes == 10
        cell_type = WriteVTK.VTKCellTypes.VTK_QUADRATIC_TETRA
        vtk_to_lib = [1,3,6,10,2,5,4,7,8,9]
    else
        return nothing
    end
    nodes -> WriteVTK.MeshCell(cell_type,(nodes[lib_to_user])[vtk_to_lib])
end

end # module

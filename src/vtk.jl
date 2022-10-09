
"""
    f = vtk_mesh_cell(elem)
Return a function `f(nodes)` that given the node ids `nodes` returns
the `WriteVTK.MeshCell` corresponding to the element `elem`.
"""
function vtk_mesh_cell end

"""
    ct = vtk_cell_type(elem)
Return the constant in `WriteVTK.MeshCellTypes` corresponding to `elem`.
"""
function vtk_cell_type end

"""
    args = vtk_args(mesh,d)

Return the arguments `args` to be passed in final position
to functions like `WriteVTK.vtk_grid`.
"""
function vtk_args(mesh,d)
    function _barrirer(
            coords,
            face_to_refid,
            face_to_nodes,
            refid_mesh_cell)
        cells = map(face_to_refid,face_to_nodes) do refid, nodes
            mesh_cell = refid_mesh_cell[refid]
            mesh_cell(nodes)
        end
        nnodes = length(coords)
        points = zeros(3,nnodes)
        for node in 1:nnodes
            coord = coords[node]
            for i in 1:length(coord)
                points[i,node] = coord[i]
            end
        end
        points, cells
    end
    coords = node_coordinates(mesh)
    face_to_nodes = face_nodes(mesh,d)
    face_to_refid = face_ref_id(mesh,d)
    refid_refface = ref_faces(mesh,d)
    refid_mesh_cell = map(vtk_mesh_cell,refid_refface)
    _barrirer(
              coords,
              face_to_refid,
              face_to_nodes,
              refid_mesh_cell)
end

"""
    physical_groups!(vtk::WriteVTK.DatasetFile,geo,dim[,groups])

Add physical groups for visualization in vtk format. The `vtk` object
needs to be associated with `geo` in dimension `dim`. If the last argument is
not given, it will be computed as `physical_groups(geo)`.

# Examples
    using GalerkinToolkit
    using Meshes
    using WriteVTK
    t = Triangle((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))
    d = 1
    groups = physical_groups(t)
    group = PhysicalGroup([1,3],d,"axes")
    push!(groups,group)
    vtk_grid("triangle_edges",vtk_args(t,d)...) do vtk
        physical_groups!(vtk,t,d,groups)
    end
"""
function physical_groups!(
    vtk::WriteVTK.DatasetFile,geo,dim,groups=physical_groups(geo))
    ndfaces = num_faces(geo,dim)
    for (i,group) in enumerate(groups)
        if group.dim == dim
            face_mask = zeros(Int,ndfaces)
            face_mask[group.faces] .= i
            vtk[group.name,WriteVTK.VTKCellData()] = face_mask
        end
    end
    vtk
end



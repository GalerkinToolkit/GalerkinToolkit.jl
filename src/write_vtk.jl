
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
    perm = vtk_node_permutation(elem)
    nodeids_vtk_order = nodeids_elem_order[perm]
"""
function vtk_cell_node_permutation end

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
    coords = node_coordinates(mesh)
    barrirer(coords)
end

function vtk_cells(mesh,d)
    function barrirer(face_to_refid,face_to_nodes,refid_mesh_cell)
        cells = map(face_to_refid,face_to_nodes) do refid, nodes
            mesh_cell = refid_mesh_cell[refid]
            mesh_cell(nodes)
        end
        cells
    end
    face_to_nodes = face_nodes(mesh,d)
    face_to_refid = face_reference_id(mesh,d)
    refid_refface = reference_faces(mesh,d)
    refid_mesh_cell = map(vtk_mesh_cell,refid_refface)
    barrirer(face_to_refid,face_to_nodes,refid_mesh_cell)
end

"""
    args = vtk_args(mesh,d)

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
    D = dimension(mesh)
    allcells = [vtk_cells(mesh,d) for d in 0:D]
    cells = reduce(vcat,allcells)
    points, cells
end

function vtk_physical_groups!(vtk,geo,d;physical_groups=physical_groups(geo,d))
    ndfaces = num_faces(geo,d)
    for group in physical_groups
        name,faces = group
        face_mask = zeros(Int,ndfaces)
        face_mask[faces] .= 1
        vtk[name,WriteVTK.VTKCellData()] = face_mask
    end
    vtk
end

function vtk_physical_groups!(vtk,geo;physical_groups=physical_groups(geo))
    nfaces = sum(num_faces(geo))
    offsets = face_offsets(geo)
    D = dimension(geo)
    data = Dict{String,Vector{Int}}()
    for d in 0:D
        for group in physical_groups[d+1]
            name, = group
            if !haskey(data,name)
                face_mask = zeros(Int,nfaces)
                data[name] = face_mask
            end
        end
    end
    for d in 0:D
        for group in physical_groups[d+1]
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


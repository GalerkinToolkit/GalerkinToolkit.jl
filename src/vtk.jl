
function vtk_mesh_cell end

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

function physical_groups!(vtk::WriteVTK.DatasetFile,mesh,d)
  groups = physical_groups(mesh)
  ndfaces = num_faces(mesh,d)
  ids = physical_group_ids(groups,d)
  for id in ids
    name = physical_group_name(groups,d,id)
    faces = physical_group_faces(groups,d,id)
    face_mask = zeros(Int,ndfaces)
    face_mask[faces] .= id
    vtk[name,WriteVTK.VTKCellData()] = face_mask
  end
  vtk
end




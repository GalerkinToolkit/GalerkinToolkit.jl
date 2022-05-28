
function vtk_cell_type end
function vtk_node_perm end

default_vtk_node_perm(a) = 1:num_nodes(a)
vtk_node_perm(a) = default_vtk_node_perm(a)

function vtk_args(mesh,d)
  function _barrirer(
    coords,
    face_to_refid,
    face_to_nodes,
    refid_vtktype,
    refid_vtkperm)
    cells = map(face_to_refid,face_to_nodes) do refid, nodes
      vtktype = refid_vtktype[refid]
      vtkperm = refid_vtkperm[refid]
      WriteVTK.MeshCell(vtktype,nodes[vtkperm])
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
  refid_vtktype = map(vtk_cell_type,refid_refface)
  refid_vtkperm = map(vtk_node_perm,refid_refface)
  _barrirer(
    coords,
    face_to_refid,
    face_to_nodes,
    refid_vtktype,
    refid_vtkperm)
end


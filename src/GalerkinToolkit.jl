module GalerkinToolkit

using StaticArrays

export prefix!
export rewind!
export JaggedArray
export JArray
include("jagged_array.jl")

export domain_dim
export ambient_dim
export face_ref_id
export ref_faces
export face_incidence
export face_nodes
export face_own_nodes
export node_coordinates
export periodic_nodes
export hanging_nodes
export is_hanging
export is_periodic
export num_faces
export num_nodes
export is_simplex
export is_hypercube
export vertex_node
export node_vertex
export add_group!
export group_faces!
export group_faces
export group_name
export group_names
export group_id
export group_ids
export GroupCollection
export physical_groups
export void
export SimpleFEMesh
export fe_mesh
include("fe_mesh.jl")

#import Meshes
#include("meshes.jl")

end

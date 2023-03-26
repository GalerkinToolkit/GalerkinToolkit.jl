module GalerkinToolkit

using Meshes
using WriteVTK
using PartitionedArrays
using StaticArrays
using Gmsh
using P4est
using MPI
using LinearAlgebra
using FillArrays

export JaggedArray
export dimension
export embedded_dimension
export node_coordinates
export reference_faces
export face_nodes
export face_reference_id
export vtk_mesh_cell
export vtk_cell_type
export vtk_args
export vtk_physical_groups!
export has_periodic_nodes
export has_hanging_nodes
export has_physical_groups
export is_simplex
export is_hypercube
export num_nodes
export num_faces
export GenericMesh
export periodic_node_constraints
export set_phyisical_groups
export set_periodic_node_constraints
export set_haning_node_constraints
export hanging_node_constraints
export physical_group
export physical_groups
export faces_in_group
export group_name
export periodic_nodes
export periodic_to_master
export periodic_to_coeff
export hanging_nodes
export hanging_to_masters
export hanging_to_coeffs
export gmsh_mesh
export forest_from_mesh
export anchor
export level
export refine!
export balance!
export coarsen!
export partition!
export node_coordinates!
export CONNECT_CORNER
export CONNECT_FACE
export CONNECT_FULL
export find_ghost_leafs
export mesh_from_forest
export TwoPartPartition
export partition_from_mask
export dof_glue_from_forest
export permutation
export classify_mesh_nodes
export classify_mesh_nodes!
export classify_mesh_faces
export classify_mesh_faces!

include("mesh_interface.jl")
include("dof_glue.jl")
include("write_vtk.jl")
include("meshes.jl")
include("gmsh.jl")
include("p4est.jl")

end # module

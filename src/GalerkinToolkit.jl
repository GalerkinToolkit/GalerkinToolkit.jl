module GalerkinToolkit

using StaticArrays
import Meshes
using Gmsh: gmsh
import WriteVTK
import P4est
import CBinding
import MPI

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
export face_vertices
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
export group_dim
export group_dims
export group_id
export group_ids
export GroupCollection
export physical_groups
export physical_groups!
export group_nodes
export classify_nodes
export VOID
export SimpleFEMesh
export fe_mesh
export mesh_face_vertices
include("fe_mesh.jl")

export polytopal_complex
include("polytopal_complex.jl")

include("meshes.jl")

export msh_file
export with_gmsh
include("gmsh.jl")

export vtk_args
include("vtk.jl")

export p4est_mesh_refiner
export destroy!
export refine!
export coarsen!
export balance!
include("p4est.jl")

export call
export evaluate
export sample
export inverse_map
export scale
export linear_combination
export Monomial
export NodalValue
export AffineMap
export Operator
export FunctionSpace
export basis_functions
export Q_space
export P_space
export P̃_space
export S̃_space
export cartesian_product
export direct_sum
include("functions.jl")

export MathTuple
export math_tuple
include("math_tuple.jl")

end

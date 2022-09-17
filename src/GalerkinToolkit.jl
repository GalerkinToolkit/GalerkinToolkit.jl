module GalerkinToolkit

using StaticArrays
import Meshes
using Gmsh: gmsh
import WriteVTK
import P4est
import CBinding
import MPI
import BlockArrays

export prefix!
export rewind!
export jagged_array
export GenericJaggedArray
export JaggedArray
include("jagged_array.jl")

export INVALID_ID
export num_dims
export num_ambient_dims
export num_faces
export num_edges
export num_vertices
export face_nodes
export face_nodes!
export face_own_nodes
export face_own_nodes!
export face_faces
export face_vertices
export ref_faces
export ref_faces!
export face_ref_id
export face_ref_id!
export ref_face_faces
export ref_face_nodes
export ref_face_own_nodes
export num_nodes
export node_coordinates
export periodic_nodes
export periodic_nodes!
export has_periodic_nodes
export num_periodic_nodes
export PeriodicNode
export hanging_nodes
export hanging_nodes!
export has_hanging_nodes
export num_periodic_nodes
export HangingNode
export is_simplex
export is_hypercube
export physical_groups
export physical_groups!
export PhysicalGroup
export vertex_node
export node_vertex
include("mesh_interface.jl")



#export FEMesh
#export fe_mesh
#include("fe_mesh.jl")


#export vertex_node
#export node_vertex
#export polytopal_complex
#include("polytopal_complex.jl")
#
#include("meshes.jl")
#
#export msh_file
#export with_gmsh
#include("gmsh.jl")
#
#export vtk_args
#include("vtk.jl")
#
#export p4est_mesh_refiner
#export destroy!
#export refine!
#export coarsen!
#export balance!
#include("p4est.jl")
#
#export call
#export evaluate
#export sample
#export inverse_map
#export scale
#export linear_combination
#export Monomial
#export NodalValue
#export AffineMap
#export Operator
#export cartesian_product
#export vector_valued_basis
#export direct_sum
#include("functions.jl")
#
#export q_monomial_basis
#export p_monomial_basis
#export p̃_monomial_basis
#export s̃_basis
#export q_equispaced_nodes
#export p_equispaced_nodes
#include("polynomial_bases.jl")
#
#export MathTuple
#export math_tuple
#include("math_tuple.jl")

end

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
export add_physical_group!
export physical_group_faces!
export physical_group_faces
export physical_group_name
export physical_group_names
export physical_group_dim
export physical_group_dims
export physical_group_id
export physical_group_ids
export PhysicalGroupCollection
export physical_groups
export physical_groups!
export physical_group_nodes
export classify_nodes
export FEMesh
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
export cartesian_product
export vector_valued_basis
export direct_sum
include("functions.jl")

export q_monomial_basis
export p_monomial_basis
export p̃_monomial_basis
export s̃_basis
export q_equispaced_nodes
export p_equispaced_nodes
include("polynomial_bases.jl")

export MathTuple
export math_tuple
include("math_tuple.jl")

end

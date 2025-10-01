#
# ```@meta
# CurrentModule=GalerkinToolkit
# ```
#
# # Mesh partitioning
#
# One can partition mesh objects and distribute them over different process for computations
# on large-scale parallel systems. [PartitionedArrays.jl](https://github.com/PartitionedArrays/PartitionedArrays.jl)
# is used under-the-hood to distribute the vectors defining the mesh object. This allows one to use several backends
# including MPI (by default), and the debug-backend of PartitionedArrays for developing and debugging
# purposes.
#
# ## Data partition model
#
# A partitioned mesh is made of the same vectors as a sequential mesh, except for the fact that these
# vectors are of type `PVector`, the type representing a distributed vector in PartitionedArrays.
# Thus, one can work with a partitioned mesh as with a sequential mesh, as long as
# you don't index the `PVector` objects directly (which is intentionally forbidden by PartitionedArrays
# for performance reasons). The usage of partitioned meshes is transparent to the user in most of
# the cases, as one rarely works with the underlying data directly. Functions that you
# call on a partitioned mesh and access the underlying mesh data, are dispatched to parallel implementations that carefully access the partitioned data. In particular, types that build on top 
# of a partitioned mesh also distribute data using `PVector` objects. The same usage remarks apply
# to them.
#
# The mesh data in a `mesh` object is partitioned as follows:
#
# * `GT.node_coordinates(mesh)` is a `PVector` distributed over nodes containing node coordinates.
# * `GT.face_nodes(mesh,d)` is a `PVector` distributed over `d`-faces of containing *global* node ids.
# * `GT.reference_spaces(mesh,d)` is typically a small tuple which is replicated across all processes.
# * `GT.face_reference_id(mesh,d)` is a `PVector` distributed over `d`-faces  containing *global* reference indices.
# * `GT.group_faces(mesh,d)` is a `Dict{String,PVector}` object whose keys are replicated on all processes. `GT.group_faces(mesh,d)[key]` is a `PVector` distributed over the faces in this group containing *global* face ids.
# * Other vectors like `normals(mesh)` are distributed following the same principle.
#
# All these `PVector` objects only include own indices and point to global indices.  This allows us to distribute node and face data in an arbitrary way.
#
# ### Creating partitioned meshes
#
# A partitioned mesh can be created with function [`with_mesh_partitioner`](@ref).
#
# ```@docs; canonical=false
# with_mesh_partitioner
# ```
# ### Example 
#
# Create a 3D mesh from a `msh` file partitioned into 5 parts.
# Visualize the mesh using colors according to part owner.
#

module mesh_partitioning_1 # hide
import GalerkinToolkit as GT
import GLMakie
import Makie
import FileIO # hide

parts = 1:5
mesh = GT.with_mesh_partitioner(;parts) do 
    GT.mesh_from_msh(msh_file)
end

#Visualize it
fig = Makie.Figure()
elevation = 0.24π
azimuth = -0.55π
aspect = :data
ax = Makie.Axis3(fig[1,1];aspect,elevation,azimuth)
Makie.hidespines!(ax)
Makie.hidedecorations!(ax)
color = GT.FaceColor("__OWNER__")
GT.makie_surfaces!(mesh;color)
FileIO.save(joinpath(@__DIR__,"fig_meshes_12.png"),Makie.current_figure()) # hide
end # hide

# ![](fig_meshes_12.png)

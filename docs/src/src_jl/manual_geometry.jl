#
# ```@meta
# CurrentModule=GalerkinToolkit
# ```
#
# # Geometric foundations
#
# We use the following dependencies in code snippets in this page.

import GalerkinToolkit as GT
import GLMakie
import Makie
import FileIO # hide

#
#
# ## Computational Domains
#
#  The type [`AbstractDomain`](@ref) is the parent of all types representing a computational domain $\Omega\subset\mathbb{R}^d$, where $d$ is the number of space dimensions. Domains provide all information needed to perform two key operations: numerical integration, and definition of interpolation spaces.
#
# ### Single-face domains
#
# The simplest domains are defined as a single primitive geometrical shape. One can build such domains with the functions [`unit_n_cube`](@ref) and [`unit_simplex`](@ref). Both function take `N` or `Val(N)`, being `N` an integer equal to the number of spatial dimensions of the object we want to create.

# For example:
cube = GT.unit_n_cube(3)

# and

triangle = GT.unit_simplex(Val(2))

#
# !!! note
#     Many functions in GalerkinToolkit accept `N` or `Val(N)`, when expecting `N` be the number of spatial dimensions. The latter case is needed for type-stability, and the former is provided for convenience.
#
#
# Functions [`unit_n_cube`](@ref) and [`unit_simplex`](@ref) build unit hypercubes and unit simplices of any arbitrary dimension. These are the two primitive domains needed in GalerkinToolkit. The rest of the domains are obtained by transforming them.

# ### Visualizing domains
#
# Domains are easily visualized using Makie. For instance, we can visualize the cube generated above as folows.

axis = (;aspect=:data)
GT.makie_surface(cube;axis)
#GT.makie_lines!(cube;color=:black) #TODO
FileIO.save(joinpath(@__DIR__,"fig_geom_1.png"),Makie.current_figure()) # hide
nothing # hide
#
# ![](fig_geom_1.png)

# And the triangle as follows.
#
axis = (;aspect=Makie.DataAspect())
#GT.makie_surface(triangle) #TODO
GT.makie_lines(triangle;axis,color=:black)
FileIO.save(joinpath(@__DIR__,"fig_geom_2.png"),Makie.current_figure()) # hide
nothing # hide
#
# ![](fig_geom_2.png)
#
# The full details about visualization are given in Section Visualization.
#
#
# ### Multi-face domains
# 
# More complex domains are defined as the union of (a subset of) the faces in a computational mesh. Each one of the faces in a mesh
# is in turn defined as the mapping of a single-face domain (a unit hypercupe or unit simplex in practice). To build a domain from a mesh one needs to define which faces from within the mesh to use. This is done by specifying the dimension of the faces (0, 1, 2, or 3) plus a vector of strings containing the names of the face groups (called *physical faces*) to use. These face groups are included in the mesh object and can be defined as explained in section [Computational meshes](@ref). This action is performed with function [`domain`](@ref). 
#
# For instance, the following build a computational domain taking all faces of dimension 3 in the computational mesh.

assets_dir = normpath(joinpath(@__DIR__,"..","..","..","assets"))
msh_file = joinpath(assets_dir,"model.msh")
mesh = GT.mesh_from_msh(msh_file)
Ω = GT.domain(mesh,Val(3))

# The domain can be visualized as

axis = (;aspect=:data)
GT.makie_surface(Ω;axis)
GT.makie_lines!(Ω;color=:black)
FileIO.save(joinpath(@__DIR__,"fig_geom_3.png"),Makie.current_figure()) # hide
nothing # hide

# ![](fig_geom_3.png)

# In this other example, we take the faces of dimension 2 including in the face groups named `"circle"`, `"triangle"`, and `"square"`.
physical_names = ["circle", "triangle", "square"]
Γ = GT.domain(mesh,Val(2);physical_names)

#
# We also visualize this domain:

GT.makie_surface(Γ;axis)
GT.makie_lines!(Γ;color=:black)
FileIO.save(joinpath(@__DIR__,"fig_geom_4.png"),Makie.current_figure()) # hide
nothing # hide

# ![](fig_geom_4.png)

#
# ### Using computational domains
#
#
# ## Computational meshes
#
#
#

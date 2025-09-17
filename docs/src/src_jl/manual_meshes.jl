#
# ```@meta
# CurrentModule=GalerkinToolkit
# ```
#
# # Meshes
#
# A mesh object in GalerkinToolkit contains all geometrical information needed in a finite element (FE) computation.
# A mesh is a set of polygons (or polytopes in general) which we refer to as *faces* or $d$*-faces*, where $d$ is the face parametric dimension.
# We call *vertices*, *edges*, *surfaces*, and *volumes* to faces of 0, 1, 2, and 3 dimensions respectively.  
# Meshes also
# include additional metadata, including *face groups* used to identify particular faces in the mesh, e.g., to impose boundary conditions.
#
# It is worth noting that GalerkinToolkit is not a mesh generation library. The mesh implementation is designed to provide the rich geometrical information needed in FE methods, rather than mesh generation. Meshes are often generated with external tools and then transformed into GalerkinToolkit objects with helper functions such as [`mesh_from_gmsh`](@ref).
#
#
# ### Code dependencies
#
# We use the following dependencies in the code snippets in this page.

import GalerkinToolkit as GT
import GLMakie
import Makie
import StaticArrays
import FileIO # hide

# In GalerkinToolkit we load dependencies from the Julia standard library with `using` statements, and from other packages with `import` statements.
# The latter forces to qualify functions with the package name, which explicitly reveals their origin. We do not qualify functions from the standard library since they are well known by Julia programmers.

#
# ## Mesh specification
#
#
#
# All types implementing meshes are subtypes of [`AbstractMesh`](@ref).
# Important features of a mesh include:
#
# - A mesh can potentially contain faces with different number of dimensions. I.e., the same mesh object can include vertices, edges, surfaces, and volumes. The number of dimensions of a mesh is the maximum number of dimension of its faces.
# - The number of dimensions of a mesh can be smaller or equal to the number of *ambient* dimension. The latter is the number of components in a node coordinate vector.
# - A mesh might or might not represent a cell complex. However, many algorithms require to work with a cell complex.
# - *Physical* faces are defined using reference interpolation spaces and node coordinates. A physical face $F=\varphi(\hat F)$ is defined by transforming a *reference* face $\hat F$ with a mapping $\varphi: \hat F \rightarrow \mathbb{R}^D$, where D is the number of ambient dimensions.  The mapping is defined as $\varphi(\hat x) = \sum_i \hat s_i(\hat x) x_{(F,i)}$. Function $\hat s_i: \hat F \rightarrow \mathbb{R}$ is the scalar basis function number $i$ in the reference (interpolation) space of $F$. The vector $x_{(F,i)}$ contains the coordinates of the local node $i$ in face $F$.
#
#
# ## Creating a mesh
# 
#
# Arbitrary mesh objects are defined from low-level quantities with function [`create_mesh`](@ref).
#
# ```@docs; canonical=false
# create_mesh
# ```
#
# ### Example
#
# In the following example, we generate and visualize a mesh of three first order triangles. Only faces of dimension 2 are present in this example. The arrays for vertices and edges are empty.
#

#Node coordinates
T = StaticArrays.SVector{2,Float64}
node_coordinates = T[(0,0),(1,0),(0,1),(1,1),(2,0)]

#Face nodes
face_nodes_0 = Vector{Int}[]
face_nodes_1 = Vector{Int}[]
face_nodes_2 = [[1,2,3],[2,3,4],[2,4,5]]
face_nodes = [
    face_nodes_0,
    face_nodes_1,
    face_nodes_2]

#Reference spaces
reference_spaces_0 = ()
reference_spaces_1 = ()
order = 1
triangle = GT.unit_simplex(Val(2))
triangle3 = GT.lagrange_space(triangle,order)
reference_spaces_2 = (triangle3,)
reference_spaces = (
    reference_spaces_0,
    reference_spaces_1,
    reference_spaces_2)

#Create mesh
mesh = GT.create_mesh(;
    node_coordinates,
    face_nodes,
    reference_spaces)

#Visualize
axis = (;aspect=Makie.DataAspect())
shading = Makie.NoShading
GT.makie_surfaces(mesh;axis,shading)
GT.makie_edges!(mesh;color=:black)
FileIO.save(joinpath(@__DIR__,"fig_meshes_1.png"),Makie.current_figure()) # hide
nothing # hide

# ![](fig_meshes_1.png)
#
# ### Example
#
# In this other slightly more complex example, we define a mesh including
# faces of different dimensions: surfaces, edges and vertices. To be able to see
# all faces in the visualization, we need to "shrink" them. Otherwise, the surfaces would hide the edges and vertices.

#Face nodes
face_nodes_0 = [[1],[3]]
face_nodes_1 = [[1,2],[2,5],[5,4]]
face_nodes_2 = [[1,2,3],[2,3,4],[2,4,5]]
face_nodes = [
    face_nodes_0,
    face_nodes_1,
    face_nodes_2]

#Reference spaces
vertex = GT.unit_simplex(Val(0))
vertex1 = GT.lagrange_space(vertex,order)
segment = GT.unit_simplex(Val(1))
segment2 = GT.lagrange_space(segment,order)
reference_spaces_0 = (vertex1,)
reference_spaces_1 = (segment2,)
reference_spaces_2 = (triangle3,)
reference_spaces = (
    reference_spaces_0,
    reference_spaces_1,
    reference_spaces_2)

#Create mesh
mesh = GT.create_mesh(;
    node_coordinates,
    face_nodes,
    reference_spaces)

#Visualize
axis = (;aspect=Makie.DataAspect())
shrink = 0.8
GT.makie_surfaces(mesh;axis,shading,shrink)
GT.makie_edges!(mesh;dim=1,shrink)
GT.makie_vertices!(mesh;dim=0)
FileIO.save(joinpath(@__DIR__,"fig_meshes_2.png"),Makie.current_figure()) # hide
nothing # hide
#
# ![](fig_meshes_2.png)
#
# ## Creating a chain
#
# Using function [`create_mesh`](@ref) might be tedious if all faces are of the same dimension. In this case, we can use the simpler constructor [`create_chain`](@ref). It works like 
# [`create_mesh`](@ref), but we pass data only for one face dimension. The resulting object
# is still a mesh object whose type is a subtype of [`AbstractMesh`](@ref).
#
# ```@docs; canonical=false
# create_chain
# ```
# ### Example
#
# We create the mesh of the first example, but using [`create_chain`](@ref).
#
#Face nodes
face_nodes = [[1,2,3],[2,3,4],[2,4,5]]

#Reference spaces
reference_spaces = (triangle3,)

#Create mesh
chain = GT.create_chain(;
    node_coordinates,
    face_nodes,
    reference_spaces)

#Visualize
GT.makie_surfaces(chain;axis,shading)
GT.makie_edges!(chain;color=:black)
FileIO.save(joinpath(@__DIR__,"fig_meshes_1a.png"),Makie.current_figure()) # hide
nothing # hide

# ![](fig_meshes_1a.png)
#
# ## Accessing mesh data
#
# The underlying data defining a mesh can be accessed with functions listed in the docstring of the type [`AbstractMesh`](@ref).
#
# ```@docs; canonical=false
# AbstractMesh
# ```
#
#
# Several of these queries return information for faces of all dimensions present in the mesh. For instance, `GT.num_faces(mesh)` returns a vector containing the number of faces in each dimension. It is often possible to restrict these queries to a given dimension. The call `GT.num_faces(mesh,d)` returns an integer with the number of faces of dimension `d`. When information for all dimensions is returned in a vector, position `d+1` in the vector contains the information for dimension `d`. Remember that the lowest dimension is `d=0` for the vertices. In particular, `GT.num_faces(mesh,d) == GT.num_faces(mesh)[d+1]`.
#
# ### Example
#
#  Let us get some information from the mesh object we created in the previous example. First, let us get the number of dimensions, number of nodes and number of faces

GT.num_dims(mesh)

#
#

GT.num_nodes(mesh)

#
#

GT.num_faces(mesh)

# To get the number of faces in a dimension 2, we do

GT.num_faces(mesh,2)

# Similarly we can get the node coordinates and face nodes for faces of dimension 2.

GT.node_coordinates(mesh)

#
#

GT.face_nodes(mesh,2)

# To get the face node ids for all dimensions, we skip the second argument:

GT.face_nodes(mesh)

# ### Reference spaces
#
# Information about the interpolation spaces used to define the mesh physical geometry is available via the function [`reference_spaces`](@ref) and [`face_reference_id`](@ref). Given a mesh object `mesh`,
# one accesses the reference space of face number `f` in dimension `d` calling `ref_space = GT.reference_spaces(mesh,d)[r]` with `r=GT.face_reference_id(mesh,d)`. A mesh has typically a small number of reference spaces, and vector `GT.face_reference_id(mesh,d)` indicates which reference space is assigned to each mesh face. The object `ref_space` is of a type that implements the [`AbstractSpace`](@ref) interface. All details about interpolation spaces are given in the section about [Interpolation](@ref). In particular, we can get some basic information from this space like the number of degrees of freedom (DOFs), and the shape functions.
#
# ### Example
#
# We manually build the physical map that transforms a reference face into a physical face. We consider face number `f=2` in dimension `d=1` from the mesh previously defined. We need to get the reference shape functions and the face node coordinates for the given face.

#Get data
d = 1
f = 2
r = GT.face_reference_id(mesh,d)[f]
ref_space = GT.reference_spaces(mesh,d)[r]
lnode_s = GT.shape_functions(ref_space)
lnode_node = GT.face_nodes(mesh,d)[f]
nlnodes = length(lnode_node)
node_x = GT.node_coordinates(mesh)

#Define the map
φ = y -> begin
    sum(1:nlnodes) do lnode
        s = lnode_s[lnode]
        x = node_x[lnode_node[lnode]]
        s(y)*x
    end
end

#Apply the map to a given 1-d point
y = StaticArrays.SVector(0.5,)
φ(y)

# The returned value should be the mid point of the edge number 2.

# Geting information from a mesh might be tedious for the many array indirection present.
# To fix this, the library provides iterators to visit the faces of the mesh.
# These functions are fully explained in Section Iterators. We rewrite this example using an iterator object.

#Face iterator
d = 1
mesh_faces = GT.each_face(mesh,d)

#Restrict iterator at current face
face = 2
mesh_face = mesh_faces[face]
lnode_s = GT.shape_functions(mesh_face)
lnode_x = GT.node_coordinates(mesh_face)

#Define the map
φ = y -> begin
    sum(1:nlnodes) do lnode
        s = lnode_s[lnode]
        x = lnode_x[lnode]
        s(y)*x
    end
end

#Apply the map to a given 1-d point
y = StaticArrays.SVector(0.5,)
φ(y)

# ## Cell complexes
#
# The mesh object described above is general enough to describe cell or [polyhedral complexes](https://en.wikipedia.org/wiki/Polyhedral_complex), but it is not guaranteed that the mesh is indeed a cell complex. For instance, none of the meshes in the two last examples is a cell complex. The first one has no vertices nor edges. The second one has only few vertices and edges, but many are missing.
# One can complete a given mesh with all low-dimensional faces needed to be a cell complex with function [`complexify`](@ref). Calling [`is_cell_complex`](@ref) on the returned mesh, will give `true`.
#
# For function `complexify` to work, neighboring faces should share node ids. I.e., duplicated nodes are not allowed in the input mesh. Otherwise, duplicated faces might be generated, or the functions might not work at all.
#
# ### Example
#
# Let us complete the mesh we generated in the last example into a cell complex.

#Convert
mesh2 = GT.complexify(mesh)
@assert GT.is_cell_complex(mesh2)

#Visualize
GT.makie_surfaces(mesh2;axis,shading,shrink)
GT.makie_edges!(mesh2;dim=1,shrink)
GT.makie_vertices!(mesh2;dim=0)
FileIO.save(joinpath(@__DIR__,"fig_meshes_3.png"),Makie.current_figure()) # hide
nothing # hide

# ![](fig_meshes_3.png)
#
# Note that the mesh contains now all low-dimensional faces.


# ## Mesh topology
#
#
# When a mesh is a cell complex, there are well-defined face incidence relationships.
# All face incidence relations are stored in an object called mesh *topology*. A mesh topology is represented by the type [`AbstractTopology`](@ref).
#
# ```@docs; canonical=false
# AbstractTopology
# ```
#
# ### Incidence relations
#
# A mesh topology is obtained with
# function [`topology`](@ref) called on a given mesh object: `topo = GT.topology(mesh)`. The mesh needs to be a cell complex for this to work.  Then, one uses function [`face_incidence`](@ref) on the topology object to get the incidence relations. 
# `GT.face_incidence(topo,D,d)` is a long vector of small vectors of integers, often implemented with a `JaggedArray`.
# - For `d<D`, `GT.face_incidence(topo,D,d)[F]` is a vector of integers containing the ids of the faces of dimension `d` on the boundary of face number `F` of dimension `D`.
# - For `d<D`, `GT.face_incidence(topo,d,D)[f]` is a vector of integers containing the ids of the faces of dimension `D` around the face number `f` of dimension `d`.
# - For `d==D`, `GT.face_incidence(topo,D,D)[F] == [F]`.
#
# ### Example
#
# Let us get some of the incidence relations for the cell complex we generated above.

topo = GT.topology(mesh2)
surface_to_edges = GT.face_incidence(topo,2,1)

# This output is read as follows. Surface 1 has edges 1, 4, and 5 on its boundary; surface 2 has edges 5, 6, 7 on its boundary; etc.

edge_to_surfaces = GT.face_incidence(topo,1,2)

# According to this output, edge 1 touches surface 1, edge 2 touches surface 3, etc. We can also see that there are two interior edges touching two surfaces. Edge 5 touches surfaces 1 and 2, and edge 6 touches surfaces 2 and 3.
#
# ### Reference topologies
#
# The reference faces in a mesh are also cell complexes. For instance, a reference volume has surfaces, edges, and vertices on its boundary. The incidence relation between the faces in a reference face are obtained using  a *reference topology*. Reference topologies are accessed with function [`reference_topologies`](@ref) given a topology object.
#
# ```@docs;canonical=false
# reference_topologies
# ```
# The rationale behind accessing reference topologies is analogous to accessing reference spaces in a mesh. Note, however, that the face reference ids in  a mesh `mesh`, `GT.face_reference_id(mesh,d)` might be different from the ones in the corresponding mesh topology `topo`, `GT.face_reference_id(topo,d)`, since different interpolation spaces can be defined on the same reference face (e.g., in p-adaptive methods).
#
# A reference topology object behaves like any other topology object and can be queried with the methods from the [`AbstractTopology`](@ref) interface.
#
# ### Example
#
# We get the reference topology of the first 2-face of the previously generated topology, which corresponds to a reference triangle. Then, we show the incidence relation between edges and vertices, i.e., for each edge which are the vertices on its boundary.
#

#Get the reference topology
f = 1
r = GT.face_reference_id(topo,2)[f]
ref_topo = GT.reference_topologies(topo,2)[r]

#See the edge to vertex relations
edge_to_vertices = GT.face_incidence(ref_topo,1,0)

# The first edge goes from vertex 1 to 2, the second edge from vertex 1 to 3, and the third edge from vertex 2 to 3.


#
# ### Permutations
#
# If a mesh is a cell complex, each face in the boundary of a face is also explicitly contained in the mesh. However, the incidence relations of these two faces are the same but might be stored in different order.
#
# Let us consider a topology  object `topo` and two integers `d<D`. We get the incident relations `D_d = GT.face_incidence(topo,D,d)`, `d_0 = GT.face_incidence(topo,d,0)` and `D_0 = GT.face_incidence(topo,D,0)`. Consider face number `v` in dimension `D`. The `d`-faces on its boundary are given in vector `D_d[v]`. Consider the integer in `l` position in this list, namely `s=D_d[v][l]`. `s` is the id of a face of dimension `d`. The 0-faces on the boundary of `v` are  in `D_0[v]` and on the boundary of `s` are in `d_0[s]`.  Now, consider the reference topology of face `v`, namely `ref_topo_v`. We get this incidence relation from the reference topology: `ref_d_0 = GT.face_incidence(ref_topo_v,d,0)`. This gives us an alternative way of obtaining the 0-faces of `s`, namely `D_0[v][ref_d_0[l]]`. That is, we can take the id `s` and compute directly `d_0[s]`, or we can go to the neighbor face `v` and compute `D_0[v][ref_d_0[l]]` using the local id `l` corresponding to `s` in `v`. The vectors `d_0[s]` and  `D_0[v][ref_d_0[l]]` contain the same vertex ids, but not in the same order!
#
# To fix this issue, we provide the permutation `P` that transforms one vector into the other, namely `d_0[s][P] == D_0[v][ref_d_0[l]]`. For `d==1`, the permutation vector `P` is either `[1,2]` or `[2,1]` since an edge has two vertices. In general, the possible permutations are enumerated and stored in the reference topology associated with face `s`, namely `ref_topo_s`. They are accessed with function [`vertex_permutations`](@ref) in this way: `k_P = GT.vertex_permutations(ref_topo_s)`. This is a vector of vectors containing all permutations. To get the permutation `P` from this list, we use function [`face_permutation_ids`](@ref). First we get an index into the list of permutations with `k=GT.face_permutation_ids(topo,D,0)[v][l]` and using the index `k`, we get the permutation from the list `P = k_P[k]`.
#
# This information is needed in many situations, including the generation of high-order interpolation spaces and integration of jump and average terms on interior faces in discontinuous Galerkin methods.
#
# ## Face groups
#
# Face groups allow us to select specific faces in a mesh for different modeling purposes: impose boundary conditions, define different equations in different parts of the mesh etc. A face group is a vector of integers containing the ids of the faces in this group plus a string containing a name for this group. This groups are stored using a dictionary that maps strings to vectors (group names to group definitions) in a per dimension basis (one dictionary per dimension). The vector contains faces of the same dimension, but it is possible define groups containing faces of different dimensions by splitting them in a vector per dimension.
# Face groups can overlap and can be added after the mesh object is created.
#
# Face groups are accessed and added using the function [`group_faces`](@ref).
#
# ```@docs; canonical=false
# group_faces
# ```
#
# ### Common face groups
#
# GalerkinToolkit provides a number of functions that generate commonly used face groups.
# - [`group_boundary_faces!`](@ref)
# - [`group_interior_faces!`](@ref)
# - [`group_faces_in_dim!`](@ref)
#
# These functions do what the name suggests (see the docstrings for further details). The first one is often used to impose boundary
# conditions and the second one in discontinuous Galerkin methods to define interior penalty terms. They are often called under the hood when calling functions like [`boundary`](@ref) and [`skeleton`](@ref).
#
#
#
#
# ### Example
#
# Let us add some face groups to the last mesh we created. We add a group for the boundary edges, for the interior edges,
# and a group for the two first surfaces.
#
#

GT.group_boundary_faces!(mesh2;group_name="boundary")
GT.group_interior_faces!(mesh2;group_name="interior")
GT.group_faces(mesh2,2)["foo"] = [1,2]
nothing # hide

# We can also visualize the faces with colors telling if a face belongs to a group or not. We can visualize the mesh faces
# labeled as "boundary" in orange color, and the rest in blue:


color = GT.FaceColor("boundary")
blue = Makie.wong_colors()[1]
orange = Makie.wong_colors()[2]
colormap = [blue,orange]
GT.makie_surfaces(mesh2;axis,shrink,shading,color,colormap)
GT.makie_edges!(mesh2;dim=1,shrink,color,colormap)
GT.makie_vertices!(mesh2;dim=0,color,colormap)
FileIO.save(joinpath(@__DIR__,"fig_meshes_4.png"),Makie.current_figure()) # hide
nothing # hide

# ![](fig_meshes_4.png)

# Idem, but now visualizing the group "foo".

color = GT.FaceColor("foo")
GT.makie_surfaces(mesh2;axis,shrink,shading,color,colormap)
GT.makie_edges!(mesh2;dim=1,shrink,color,colormap)
GT.makie_vertices!(mesh2;dim=0,color,colormap)
FileIO.save(joinpath(@__DIR__,"fig_meshes_5.png"),Makie.current_figure()) # hide
nothing # hide
#
# ![](fig_meshes_5.png)
#
# ## Gmsh meshes
#
# Meshes generated with gmsh can be transformed into GalerkinToolkit mesh objects using functions [`mesh_from_gmsh`](@ref) and [`mesh_from_msh`](@ref).
# The *physical groups* defined within gmsh will be transformed into face groups in the GalerkinToolkit mesh, which is useful to impose boundary conditions.
# We use [Gmsh.jl](https://github.com/JuliaFEM/Gmsh.jl) under the hood as the wrapper to the Julia API of gmsh.
#
# ### A mesh from a msh file
#
# Function [`mesh_from_msh`](@ref) reads and crates a mesh object from a `.msh` file (the default format used by gmsh to export meshes).
#
# ```@docs; canonical=false
# mesh_from_msh
# ```
#
# ### Example
# 
# We create a mesh from a `.msh` file and visualize it.
# In this case, we only visualize the 2-faces in the mesh.
# We color them according to the face group named `"sides"`.
# Faces in the group are visualized in orange and other faces in blue.
# This face group is only defined for 2-faces. If you visualize the
# 3-faces (as by default), you would not see this face group.
#

#Read the mesh
assets_dir = normpath(joinpath(@__DIR__,"..","..","..","assets"))
msh_file = joinpath(assets_dir,"model.msh")
mesh = GT.mesh_from_msh(msh_file)

#Visualize it
fig = Makie.Figure()
ax = Makie.Axis3(fig[1,1];aspect=:data)
Makie.hidespines!(ax)
Makie.hidedecorations!(ax)
color = GT.FaceColor("sides")
GT.makie_surfaces!(mesh;dim=2,color,colormap)
GT.makie_edges!(mesh;color=:black)
FileIO.save(joinpath(@__DIR__,"fig_meshes_6.png"),Makie.current_figure()) # hide
nothing # hide

# ![](fig_meshes_6.png)

# ### Meshes from the gmsh API
#
# It is also possible to generate meshes in Julia code using the gmsh API and then convert them to GalerkinToolkit objects.
# This is done with two functions: See also [`with_gmsh`](@ref) and [`mesh_from_gmsh`](@ref).
#
# ```@docs; canonical=false
# with_gmsh
# mesh_from_gmsh
# ```
#
# ### Example
#
# We generate a simple 2d mesh with the gmsh Julia API.
#

#Generate mesh with GMSH Julia API
mesh = GT.with_gmsh() do gmsh
    mesh_size=0.04
    T=2
    N=100
    R = 0.15
    dim = 2
    gmsh.option.set_number("General.Verbosity", 2)
    rect_tag = gmsh.model.occ.add_rectangle(0,0,0,1,1)
    circle_tag = gmsh.model.occ.add_circle(0.5,0.5,0,R)
    circle_curve_tag = gmsh.model.occ.add_curve_loop([circle_tag])
    circle_surf_tag = gmsh.model.occ.add_plane_surface([circle_curve_tag])
    gmsh.model.occ.cut([(dim,rect_tag)],[(dim,circle_surf_tag)]);
    gmsh.model.occ.synchronize()
    domain_tags = [1]
    outer_tags = [6,7,8,9]
    inner_tags = [5]
    gmsh.model.model.add_physical_group(dim,domain_tags,-1,"domain")
    gmsh.model.model.add_physical_group(dim-1,outer_tags,-1,"outer")
    gmsh.model.model.add_physical_group(dim-1,inner_tags,-1,"inner")
    gmsh.option.set_number("Mesh.MeshSizeMax",mesh_size)
    gmsh.model.mesh.generate(dim)
    #Transform it to a mesh object
    GT.mesh_from_gmsh(gmsh)
end

#Visualize
fig = Makie.Figure()
ax,sc = GT.makie_surfaces(fig[1,1],mesh;axis,shading)
GT.makie_edges!(mesh;color=:black)
Makie.hidespines!(ax)
Makie.hidedecorations!(ax)
FileIO.save(joinpath(@__DIR__,"fig_meshes_7.png"),Makie.current_figure()) # hide
nothing # hide

# ![](fig_meshes_7.png)
#
#
# ## Cartesian meshes
#
# GalerkinToolkit comes with a built-in Cartesian mesh generator implemented in function [`cartesian_mesh`](@ref).
#
# ```@docs
# cartesian_mesh
# ```
#
# ### Example
#
# Generate a 3D Cartesian mesh  of the box $(0,1)\times(-1,1)\times(3,4)$ with
# 10, 20 , and 10 cells in each axis.

#Create mesh
domain = (0,1,-1,1,3,4)
cells = (10,20,10)
mesh = GT.cartesian_mesh(domain,cells)

#Visualize it
GT.makie_surfaces(mesh;axis=(;aspect=:data))
GT.makie_edges!(mesh;color=:black)
FileIO.save(joinpath(@__DIR__,"fig_meshes_8.png"),Makie.current_figure()) # hide
nothing # hide
#
# ![](fig_meshes_8.png)
#
# ### Example
#
# Create a mesh with the same geometry as before, but using simplex cells instead. We also manually add a face group named `"foo"` containing all volumes to the left of the plain $y=0$. To build this group, we check if the cell mid point is on the left of the plain. This example also shows how to use accessors to compute
# the midpoint of each volume in the mesh.

#Create mesh
domain = (0,1,-1,1,3,4)
cells = (10,20,10)
mesh = GT.cartesian_mesh(domain,cells;simplexify=true)

#Find faces in new group
mesh_faces = GT.each_face(mesh)
new_group_faces = findall(mesh_faces) do mesh_face
    lnode_x = GT.node_coordinates(mesh_face)
    xm = sum(lnode_x) / length(lnode_x)
    xm[2] < 0
end

#Add new group to mesh
d = GT.num_dims(mesh)
GT.group_faces(mesh,d)["foo"] = new_group_faces

#Visualize
color = GT.FaceColor("foo")
GT.makie_surfaces(mesh;color,colormap,axis=(;aspect=:data))
GT.makie_edges!(mesh;color=:black)
FileIO.save(joinpath(@__DIR__,"fig_meshes_9.png"),Makie.current_figure()) # hide
nothing # hide

# ![](fig_meshes_9.png)
#
# ### Example
#
# Create a coarse cartesian mesh of the unit square with 4 and 4 cells in each axis. Visualize faces in all dimensions, shrinking them to avoid overlaps.
#

#Generate mesh
domain = (0,1,0,1)
cells = (4,4)
mesh = GT.cartesian_mesh(domain,cells)

#Visualize it
GT.makie_surfaces(mesh;axis,shading,shrink)
GT.makie_edges!(mesh;dim=1,shrink)
GT.makie_vertices!(mesh;dim=0)
FileIO.save(joinpath(@__DIR__,"fig_meshes_10.png"),Makie.current_figure()) # hide
nothing # hide
#
# ![](fig_meshes_10.png)
#
# Now, the same as before, but only generate low-dimensional faces on the boundary.


#Generate mesh
domain = (0,1,0,1)
cells = (4,4)
mesh = GT.cartesian_mesh(domain,cells;complexify=false)

#Visualize it
GT.makie_surfaces(mesh;axis,shading,shrink)
GT.makie_edges!(mesh;dim=1,shrink)
GT.makie_vertices!(mesh;dim=0)
FileIO.save(joinpath(@__DIR__,"fig_meshes_11.png"),Makie.current_figure()) # hide
nothing # hide
#
# ![](fig_meshes_11.png)

# # Mesh topology
# !!! note
#     TODOs
#     - The API for $\text{mesh}(\hat F)$ can be improved and for $\text{mesh}(F)$ is not available. The latter is not needed.
#     - Rename `cell_complex` by `polyhedral_complex` ? 
#
# Geometry alone is not enough in many Finite Element methods.
# One also needs an efficient way of describing the incidence relations between the faces.
# Informally, two faces are incident if one "touches" the other.
# Incidence relations endows a mesh with a graph-like structure that we refer
# to as the *mesh topology*. This structure is completely independent from the values of the node
# coordinates and only depends on topological information such as node and face ids.
# The mesh topology is often needed to "glue" faces together, 
# e.g., to build (high-order) continuous FE spaces,
# or to integrate at face interfaces.
#
# ## Topology objects
#
# Given a mesh object `M::AbstractMesh`, one extracts its topological information with
# function `T=topology(M)`. The returned object `T` is an instance of a type `<:AbstractTopology`
# that encodes the incidence relations of the faces in mesh `M`.
# In particular, it encodes which faces in the mesh have non-empty intersections, and also
# additional data used to define a common parametrization of the intersection for all incident faces.
#
# ## Face incidence
#
# Two faces $F_1$ and $F_2$ are *incident* if their closures have non-empty intersection,
# $\bar F_2 \cap \bar F_2 \neq \empty$.
# The incidence relations among all faces in a mesh are obtained from 
# a topology object `T::AbstractTopology`  with function `face_incidence(T,m,n)`.
# The return value is a large vector of small vectors of integers, often represented with a `PartitionedArrays.JaggedArray` object. 
# The meaning of this vector depends on the values `m` and `n`.
#
# - For `m<n`, `nfaces = face_incidence(T,d1,d2)[mface]` is a vector containing the ids of the faces of dimension `n` that are incident to the face of dimension `m` with id `mface`. We call the faces with id `nfaces` the *faces around* the face with id `mface`.
#
# ## Face-local faces
#
# Let $\hat F$ a reference face of dimension $d$ and $\text{mesh}(\hat F)$ the set containing 
# $\hat F$ and all faces of dimension lower than $d$ on the boundary of $\hat F$.
# For a reference edge, $\text{mesh}(\hat F)$ contains the edge and the two end vertices.
# For a reference square, it contains the square, the four edges defining the boundary, and the four vertices
# at the intersection of the edges. For a reference cube, it contains the cube, six surfaces, twelve
# edges, and eight vertices (see the next figure).  We call $\text{mesh}(\hat F)$ the
# *face-local* faces (or simply the local faces) of $\hat F$.
# We define  the local faces of a physical face $F$, namely $\text{mesh}(F)$, by mapping
# the local faces of its reference face $\hat F$ with the physical map $\phi:\hat F \rightarrow F$
# that transforms $\hat F$ into $F$.
# The sets $\text{mesh}(\hat F)$ and $\text{mesh}(F)$ are meshes as defined in the previous page, and thus they can be represented using a mesh instance with type `<:AbstractMesh`.
# In the API, the mesh $\text{mesh}(\hat F)$ for the reference face $\hat F$ is obtained from a physical face `F::AbstractMeshFace` as `mesh(reference_space(F))`.
#
#
# ## Cell Complex
#
#
#
#
#
# 
# 
#
# Two faces $F_1$ and $F_2$ are *incident* if their closures have non-empty intersection,
# $\bar F_2 \cap \bar F_2 \neq \empty$. Face incidence relations
#
#
# Given a faces $F$ in a mesh $M$, our goal is to find which other faces in the mesh $M$ have a non-empty intersection with the local faces of $F$. In other words,
# For each local face $f\in\text{mesh}(F)$ of a face $F\in M$ in a mesh $M$, our goal is to find another face $F\prime \in M$ such that $\bar f \cap \bar F\prime\neq \empty$
#
# we look for faces $F\prime\in M$ such that
#
#
# ### Cell complex
#
# A mesh $M$ is a *cell complex* (or [polyhedral complex](https://en.wikipedia.org/wiki/Polyhedral_complex))
#
# Let `F,F\prime\in M` two faces in a mesh $M$ and $\hat F$, $\hat F\prime$  their respective
# reference faces.
# 

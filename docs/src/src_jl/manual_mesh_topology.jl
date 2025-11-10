# ```@meta
# CurrentModule=GalerkinToolkit
# ```
#
# # Mesh topology
# !!! note
#     TODOs
#     - The API for $\text{mesh}(\hat F)$ can be improved and for $\text{mesh}(F)$ is not available. The latter is not needed.
#     - Rename `cell_complex` by `polyhedral_complex` ?  or face complex?
#
#
# The mesh topology describes how the faces in a mesh are connected, i.e., how faces
# "touch" each other. We formalize the meaning of "touching" by
# introducing the concept of *face incidence*.
# This information is needed to "glue" faces together, 
# e.g., to build (high-order) continuous finite element (FE) spaces,
# or to integrate at face interfaces.
# In the API, one extracts the face connections 
# of a given mesh object `M::AbstractMesh` with
# function `T=topology(M)`.
# The returned object `T` is an instance of a type
# that implements the [`AbstractTopology`](@ref) interface.
# We detail here the main methods in this interface.
#
#
# ## Face incidence
#
#
# In section [Mesh geometry](@ref), we defined the local faces $\hat M$ of a reference $\hat F$.
# We extend this definition to a physical face $F$ as follows.
# The set of local faces $\tilde F$ of a physical face $F$  is
# the set  containing the image of all
# faces in $\hat M$ via the physical map $\phi$ of $F$,
# namely $\tilde F := \{ \phi(f) \text{ for } f\in\hat M \}$.
#
# We say that faces $F_1\in M$ and $F_2\in M$ in a mesh $M$ are *incident*
# (or that one is *incident to* the other)
# if they share a common local face, namely
# $\tilde F_1 \cap \tilde F_2 \neq \empty$.
# Note that this requirement is stronger than having a non-empty
# interface $\bar F_1 \cap \bar F_2 \neq \empty$ since the interface needs to coincide
# with a local face both of $F_1$ and $F_2$.
#
# When two faces
# $F_1$ and $F_2$ are incident, the common face 
# $\tilde F_1 \cap \tilde F_2$ is identified with two different ids, namely
# $i_1$ and $i_2$, one according to the face ids for $\tilde F_1$
# and another according to the ids for $\tilde F_2$. We say that face $F_1$
# is incident to the local $d$-face of $F_2$ with id $i_2$,
# where $d$ is the dimension
# of $\tilde F_1 \cap \tilde F_2$.
# Conversely, $F_2$ is incident to the local
# $d$-face of $F_1$ with id $i_1$.
#
# Given a mesh `M::AbstractMesh` and its corresponding topology object
# `T=topology(M)`, one extract the incidence relations between faces faces in mesh
# `M` with function `face_incidence(T,m,n)`.
# The return value is a large vector of small vectors of integers, often represented with a `PartitionedArrays.JaggedArray` object.
# The item `nfaces = face_incidence(T,m,n)[mface]` is a vector containing the
# ids of the faces of dimension `n` that are incident to the face  of dimension `m` 
# and
# id `mface::Integer`.
# The meaning of this vector depends on the values `m` and `n`.
#
#
# - For `m>n`, the faces in `nfaces` are on the boundary of face `mface`. For an integer `i`, `nfaces[i]` is the face incident to the local face `i` of face  `mface`. This encodes which faces are on the boundary of another face, and also which is the local face of `mface` being touched.
# - For `m<n`, we call the faces in `nfaces` the *faces around* `mface`. For an integer `i`, we call `nfaces[i]` the `i`-th face around `mface`. Faces in `nfaces` are also often called the *co-boundary* of `mface`. E.g., for `m=3` and `n=1`, this information describes several volumes `nfaces` sharing a common edge `mface`.
# - For `m=n`, `nfaces` are the faces that overlap with `mface`. We always assume that faces are non-overlapping. Thus, `nfaces == [mface]`, i.e., a face overlaps with itself.
#
#
# ## Face complex
#
# Our implementation of function `face_incidence`, and many other algorithms
# in the library, assume that the underlying mesh is a *face complex*
# (or [polyhedral complex](https://en.wikipedia.org/wiki/Polyhedral_complex)).
# We say that a mesh $M$ is face complex if it has the following properties:
#
# * The mesh $M$ is *conforming*, and
# * the local faces of any face are also faces of the mesh, namely $\tilde F \subset M$ for all $F\in M$.
#
# A conforming mesh is such that if two faces have a non-empty interface,
# then the two faces are incident, namely  $\bar F_1 \cap \bar F_2\neq\empty$ implies $\tilde F_1 \cap \tilde F_2\neq\empty$ for any $F_1,F_2\in M$.
#
# If needed, one can transform a given mesh `M::AbstractMesh`  into another mesh
# that is a face complex with function `complexify(M)`. This function adds the missing
# faces in `M` to make it a Polytopal Complex. Note that this function assumes
# that the input mesh is already conforming. Note that function `mesh_from_gmsh`,
# `mesh_from_msh`, and `cartesian_mesh` return a mesh that is a face complex (assumingthat the mesh generator was able to create a conforming mesh).
#
#
# 

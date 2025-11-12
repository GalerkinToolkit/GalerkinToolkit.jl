# ```@meta
# CurrentModule=GalerkinToolkit
# ```
#
# # Mesh topology
# !!! note
#     TODOs
#     - The API for $\text{mesh}(\hat F)$ can be improved and for $\text{mesh}(F)$ is not available. The latter is not needed.
#     - Rename `cell_complex` by `polyhedral_complex` ?  or face complex?
#     - Define $\hat x_n$ and $g^F_n$ in section Mesh geometry
#
#
# The mesh topology describes how the faces in a mesh are connected.
# This information is needed to "glue" faces together, 
# e.g., to build (high-order) conforming finite element (FE) spaces,
# or to integrate at face interfaces.
# In the API, one extracts the topology information 
#  of a given mesh `M::AbstractMesh` with
#  function `T=topology(M)`.
#  The returned object `T` is an instance of a type
#  that implements the [`AbstractTopology`](@ref) interface.
#  We detail here the main methods in this interface.
#
# ## Gluing faces
#
# By "glue" the interface $Γ_{1,2}$ of two faces $F_1\in M$ and $F_2\in M$.
#  
#  1. Find the local faces $f_i\in L(R(M,F_i))$ such that $\Gamma_{1,2}=\phi(M,F_i)(f_i)$ for $i=1,2$.
#
#  1. Find two maps $\phi^?_1$ and $\phi^?_2$ such that $\phi(M,F_1)(\phi^?_1(x))=\phi(M,F_2)(\phi^?_2(x))$  for all $x\in R(L(R(M,F_i)),f_i)$ for $i=1,2$.
#
#  This is needed to integrate on the interface quantities functions defined on the faces $F_i$ and to build finite element spaces
#  that are conforming on the interface $\Gamma_{1,2}$. How to achive this from the solution on problem X is explained in the corresponding sections.
#
#  Assumptions:
#
#  1. To be able to find the local faces $f_i$ the mesh needs to be conforming: if $\bar F_1\cap \bar F_2$ then exist the local faces $f_i\in L(R(M,F_i))$ such that $\Gamma_{1,2}=\phi(M,F_i)(f_i)$ for $i=1,2$.
#
#  1. For convenience, we assume that the interface $\Gamma_{1,2}\in M$ is also represented explicitly as a face in mesh $M$.
#
#  A face that fulfills 1 and 2 for all pairs of faces is called a *face complex*
# (or [polyhedral complex](https://en.wikipedia.org/wiki/Polyhedral_complex)).
#
# If needed one can use function complexify (assuming that the input mesh is already conforming).
#
# ## Face incidence
#
# The topology is a graph $T$, where the $M$ are the vertices of the graph. The edges of the graph as follows.
#
# The edges $(F_1,F_2)$ and $(F_2,F_1)$ exists if and only if
#
# - the two faces are the same, $F_1 = F_2$, or
# - the two faces are of different dimensions and have a non-empty interface, $\Gamma(F_1,F_2) \neq \empty$.
#
#  If two faces are connected by an edge in the graph $T$, we say that the faces are *incident*, or
#  *adjacent*.
#
#  In the API, the object `T::AbstractTopology` that represents the graph $T$ is obtained
#  with `T=topology(M)` from a mesh object `M::AbstractMesh`.
#  The graph `T` is encoded using [adjacency lists](https://en.wikipedia.org/wiki/Adjacency_list)
#  organized according to face dimensions. We have as many lists as pairs $(m,n)$. The adjacency list for a given pair
#  $(m,n)$ contains the edges $(F_m,F_n)$ where $F_m$ and $F_n$ are faces of dimensions $m$ and $n$ respectively.
#  In the API, the adjacency list for the pair `(m,n)` is returned by function `face_incidence(T,m,n)`.
#  It is a long vector of small vectors of integers encoded via a JaggedArray.
#
#  Let $F_m$ be the `m`-face in $M$ with id `mface`. The vector `nfaces=face_incidence(T,m,n)[mface]`
#  contains the ids of the `n`-faces $F_n\in M$ that are incident to face $F_m$ in the graph T.
#
#
#  For `m>n`, the face
#  ids in `nfaces` are sorted according to the id of the local faces of $F_m$. Let $f_n$ be the local
#  `n`-face in $L(R(M,F_m))$ with id `i`. Then, `nfaces[i]` is the id of the `n`-face $F_n\in M$ on the boundary of
#  $F_m$ adjacent to local face $f_n$.
#  For `m<n`, the ids in `nfaces`
#  are arbitrarily sorted. In this case, we call the `n`-face with id `nfaces[i]` the `i`-th *face around* $F_m$.
#   For `m=n`, the vector `nfaces` is simply equal to `[mface]` (a vector with a single integer).
#
#  Form this information we can "glue" two faces as follows: Starting from the interface , we get the faces around, then for each
#  face around we can find the index of the local face around.
#
#  ## Vertex permutations
#
#  The map $\phi^?_i$ is defined as $\phi(L(R(M,F_i)),f_i)\circ \varphi(L(R(L(R(M,F_i)),f_i)),P)$ where $P$ is a permutataion
#  matrix.
#
# ## Reference topologies
#
# ## Definition
#
# We define the topology $T$ of a mesh $M$ as a graph, where the vertices of the graph
# are the faces in $M$ and the edges of the graph are defined as follows.
# Two faces $F_1\in M$ and $F_2\in M$ are connected by an edge in the graph $T$
# if $F_1=F_2$ or the faces are of different dimensions and exist a local face $f_1 \in L(F_1)$ and a local face
# $f_2\in L(F_2)$ with equivalent global node ids, namely $N(F_1,f_1)\equiv N(F_2,f_2)$.
# Moreover, to vectors of node ids $N_1$ and $N_2$ are equivalent
#  $N_1\equiv N_2$ if they contain the same ids up to a permutation,
#  i.e., 
#  exists a permutation matrix $P$ such that $N_1 = P N_2$.
#  If two faces are connected by an edge in the graph $T$, we say that the faces are *incident*,
#  *adjacent*, or *neighbors*.
#
#  ## Face incidence
#
#  In the API, the object `T::AbstractTopology` that represents the graph $T$ is obtained
#  with `T=topology(M)` from a mesh object `M::AbstractMesh`.
#  The graph `T` is encoded using [adjacency lists](https://en.wikipedia.org/wiki/Adjacency_list)
#  organized according to face dimensions.
#  The adjacency list encoding the connection of faces of dimension `n` with faces
#  of dimension `m`  is returned by function `face_incidence(T,m,n)`.
#
# Let $F_m$ be the `m`-face in $M$ with id `mface`. The vector `nfaces=face_incidence(T,m,n)[mface]`
# contains the ids of the `n`-faces $F_n\in M$ that are incident to face $F_m$. For `m>n`, the face
# ids in `nfaces` are sorted according to the id of the local faces of $F_m$. Let $f_n$ be the local
# `n`-face in $L(F_m)$ with id `i`. Then, `nfaces[i]` is the id of the `n`-face $F_n\in M$ connected
# to $F_m$ via the local face $f_n$, namely $N(F_n)\equiv N(F_m,f_n)$. For `m<n`, the ids in `nfaces`
# are arbitrarily sorted. In this case, we call the `n`-face with id `nfaces[i]` the `i`-th *face around* $F_m$.
# For `m=n`, the vector `nfaces` is simply equal to `[mface]` (a vector with a single integer), assuming that
# all faces in the mesh have unique node ids.
#
#
# ### Notation
#
# * notation: $M$ is a mesh
# * notation: $F\in M$ is a face in mesh $M$
# * notation: A face is uniquely identified by the mesh it belongs, its face id, and dimension.
# * notation: $N(F)$ is the vector of node ids for a face $F\in M$
# * notation: $N(F,f)$ is the vector of node ids of face $F\in M$ restricted to the local face $f\in\hat M$.
# * notation: Note that $N(F)=N(F,\hat F)$, where $\hat F$ is the reference face of $F$.
# * notation: $\hat M(F)$ are the (reference) local faces of $F$.
# * notation: Two vector of node ids $N_1$ and $N_2$ are equivalent $N_1\equiv N_2$ if they contain the same ids up to a permutation. I.e., exists a permutation matrix $P$ such that $N_1 = P N_2$.
# 
#  Note that this requirement is stronger than having a non-empty
#  interface $\bar F_1 \cap \bar F_2 \neq \empty$ since the interface needs to coincide
#  with a local face both of $F_1$ and $F_2$ with matching global ids.
#
#  If two faces $F_1$ and $F_2$ are incident, then they share a local face,
#  meaning $L(F_1) \cap L(F_2) \neq \empty$. Note however that the converse is not true,
#  two faces might share a local face, but the node ids of the local face might be different seen
#  from $F_1$ and $F_2$, namely $L(F_1) \cap L(F_2) \neq \empty$ does not necessarily imply $N(F_1,f)\equiv N(F_2,f)$
#  with $f\in L(F_1) \cap L(F_2)$. Examples in which the implication is not true are meshes
#  with duplicated node coordinates or with faces of non-matching interpolation order.
#
#
# `nfaces = face_incidence(T,m,n)[mface]`
#
#
# $F_m$ is the `m`-face in $M$ with id `mface`.
# The vector `nfaces` contains the ids of the `n`-faces $F_n\in M$ that are incident to face $F_m$.
#
#  $f_n$ is the `n`-face in $\hat M$ with id `i`.
# `nfaces[i]` is the id of the `n`-face $F_n\in M$ such that $N(F_n)\equiv N(F_m,f_n)$
#
# Two faces are incident if they touch via a local face of 
#
# Two faces $F_1$ and $F_2$ are incident if exist a local face $f_1 \in \hat M(F_1)$ and a local face
# $f_2\in \hat M(F_1)$ such that $N(F_1,f_1)\equiv N(F_2,f_2)$.
#
# Geometrical version
#
# Two faces $F_1$ and $F_2$ are incident if they share a common local face, namely $L(F_1)\cap L(F_2)\neq \empty$.
#
#
# `nfaces[i]` is the id of the `n`-face $F_n\in M$ such that $F_n=f_n$.
#
# ## Local face node ids
#
#  Let us consider a physical face $F\in M$ in a mesh $M$ together with its reference
#  space $\hat V$, reference face $\hat F$, and reference local faces $\hat M$
#  In addition, let $\hat x_n$ be the coordinate of node number $n$ in the reference space $\hat V$
#  and $g^F_n$ the corresponding global node id in face $F$.
#  See section [Mesh geometry](@ref) for the definition of these objects.
#
#  With this notation, we introduce $\hat N(f)$ the set of node ids of space $\hat V$ at a local face $f\in \hat M$,
#  as the set containing ids of the nodes on $\bar f$,
#  namely $\hat N(f):=\{ n \text{ such that } \hat x_n\in\bar f\}$.
#
#  $id^{-1}(\texttt{mface},\texttt{m})$ `mface` `m` $\mathcal{N}$
#
#
# $\hat N(f)$
#
# $\tilde F := \{ N^F(f) \text{ for }f\in\hat M\}$
#
#
# `nfaces` are the ids 
# The vector `nfaces` contains the ids of the `n`-faces $F_n\in M$ that are incident to face $F_m$.
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

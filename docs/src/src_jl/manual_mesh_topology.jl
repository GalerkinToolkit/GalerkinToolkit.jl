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
#     - The code is actually working with the inverse of the permutation matrix $P$.
#
#
# The mesh topology describes how the faces in a mesh are connected.
# This information is needed to "glue" faces together, 
# e.g., to build (high-order) conforming finite element (FE) spaces,
# or to integrate at face interfaces.
#
# ## Local faces
#
# Let $F\in M$ be a physical face in mesh $M$, where $\hat F:=R(F)$ is its reference
# face. Let $\hat C := C(\hat F)$ be the chain used to define the geometry of
# the reference face $\hat F$ and let $d$ be the dimension of $\hat F$. See
# Section [Mesh geometry](@ref) for the definition of these objects.
# Let us also define the interface
# $\Gamma(\mathcal{F})$ of a set of faces $\mathcal{F}$
# as the open set such that its closure is the intersection of the closures
# of the faces in $\mathcal{F}$, namely 
# $\overline{\Gamma(\mathcal{F})} = \cap_{G\in \mathcal{F}} \bar G$.
# If the interface is not empty, $\Gamma(\mathcal{F})\neq\empty$, we call the faces in
# $\mathcal{F}$ the *faces around* $\Gamma(\mathcal{F})$.
#
# With these notations, we define the set of *local faces* of the reference face $\hat F$
# as the mesh $L(\hat F)$ containing:
# 1. the reference face $\hat F$,
# 2. the faces in the chain $\hat C$,
# 3. and the interfaces $\Gamma(\{A,B\})$ of any pair of faces $A,B\in\hat C$.
#
# The set $L(\hat F)$ contains all faces in the boundary of the polytope $\hat F$. E.g.,
# for a reference edge, it contains the edge and the two end vertices.
# For a reference square, it contains the square, four edges, and the four vertices
# at the intersection of the edges. For a reference cube, it contains the cube, six surfaces, twelve edges, and eight vertices (see the next figure).
#
# In the code, one can access the mesh of local faces
# from a reference face `Fref::AbstractDomain` with `Mref = mesh(Fref)`. The mesh
# `Mref::AbstractMesh` is like any other mesh used in the code and, e.g., it can be visualized as
# any other mesh shown in next Figure.
# 
# ![](fig_meshes_defs_5.png)
#
# **Figure**: Visualization of the local faces of a reference cube, square, and segment.
# Faces are shrunk in the visualization to illustrate that they are open sets.
#
module Meshes05 # hide
#Code used to generate the figure
import GalerkinToolkit as GT
import GLMakie as Makie
import FileIO # hide
domain = (0,2,0,1,0,2)
cells = (4,2,4)
F3 = GT.unit_n_cube(Val(3))
F2 = GT.unit_n_cube(Val(2))
F1 = GT.unit_n_cube(Val(1))
M3 = GT.mesh(F3)
M2 = GT.mesh(F2)
M1 = GT.mesh(F1)
fig = Makie.Figure()
aspect = :data
ax = Makie.Axis3(fig[1,1];aspect)
Makie.hidespines!(ax)
Makie.hidedecorations!(ax)
shrink = 0.8
GT.makie_surfaces!(ax,M3;dim=2,shrink)
GT.makie_surfaces!(ax,M3;dim=3,shrink)
GT.makie_edges!(ax,M3;dim=1,shrink)
GT.makie_vertices!(ax,M3;dim=0,shrink)
aspect = Makie.DataAspect()
axis = (;aspect)
ax, = GT.makie_surfaces(fig[1,2],M2;dim=2,shrink,axis)
Makie.hidespines!(ax)
Makie.hidedecorations!(ax)
GT.makie_edges!(M2;dim=1,shrink)
GT.makie_vertices!(M2;dim=0,shrink)
ax = Makie.Axis(fig[1,3])
Makie.hidespines!(ax)
Makie.hidedecorations!(ax)
GT.makie_edges!(ax,M1;dim=1,shrink)
GT.makie_vertices!(ax,M1;dim=0,shrink)
FileIO.save(joinpath(@__DIR__,"fig_meshes_defs_5.png"),fig) # hide
end # hide
nothing # hide

# ## Conforming meshes
#
# A mesh $M$ is *conforming* if for any pair of faces $F_1,F_2\in M$
# with non-empty interface $\Gamma(\{F_1,F_2\})\neq\empty$ exists a local
# face $f_1\in L(R(F_1))$ and a local face  $f_2\in L(R(F_2))$ such that
#
# * their images via the physical map coincide with the interface, $\Gamma(\{F_1,F_2\})=\phi^{F_1}(f_1)=\phi^{F_2}(f_2)$, and
# * they share global node ids, $N(F_1,f_1)=PN(F_2,f_2)$.
#
# If these local faces exist, we call $f_1$ (resp. $f_2$) the local
# face of $F_1$ (resp $F_2$) *equivalent* to $\Gamma(\{F_1,F_2\})$.
# In this definition, $P$ is a permutation matrix and $N(F,f)$ is a vector containing the node ids of $F$ restricted to the local face $f\in L(R(F))$, namely $[N(F,f)]_r=[N(F)]_l$ with $l=[N(f)]_r$. Both local faces have the same global ids but these ids are allowed to be order differently in each local face. This is why we include the permutation matrix $P$.
#
# ## Face complexes
#
# A mesh $M$ is a *face complex* (or [polyhedral complex](https://en.wikipedia.org/wiki/Polyhedral_complex)) if
# * it is conforming, and
# * it contains the images of all local faces, namely $\phi^F(f)\in M$ for all $f\in L(R(F))$ and all $F\in M$.
#
# Given a mesh $M$ that is conforming, but not a face complex, it is always possible to add additional faces to make it a face complex. In the API, this is done
# with function `M2 = complexify(M)` from a mesh object `M::AbstractMesh`. The result
# is `M2::AbstractMesh`. In addition function `is_face_complex(M)` checks if
# `M` is a face complex. It returns true form meshes created with
# function `complexify`.
# If the chain $\hat C$ used to define the local faces $L(\hat F)$ above is conforming,
# the set of local faces $L(\hat F)$ is a face complex. In fact, we
# create the local faces in the library by applying `complexify` to the chain $\hat C$.
#
#
# ## Gluing faces
#
# Let us consider a face complex $M$ and a set of faces $\mathcal{F}\subset M$ of dimension $n$ with a non-empty interface $\Gamma(\mathcal{F})\neq\empty$, being a face of dimension $m$.
# We refer to *gluing* faces $\mathcal{F}$ at the their interface $\Gamma(\mathcal{F})$ as solving the following two-step problem for each face around $F\in\mathcal{F}$:
#
#
# 1. Find the local face $f\in L(R(F))$ and the permutation matrix $P$ such that $N(\Gamma(\mathcal{F}))=PN(F,f)$.
# 2. Find a map $\varphi$ such that $\phi^F(\phi^f(\varphi(x)))=\phi^{\Gamma(\mathcal{F})}(x)$ for all $x\in R(\Gamma(\mathcal{F}))$.
#
# For this problem to have a solution, the mesh $M$ needs to be conforming. This is what allows
# us to find the local face $f$, the permutation matrix $P$, and the map $\varphi$. In addition, mesh $M$ needs to be a face complex so that the interface $\Gamma(\mathcal{F}))\in M$ is also
# a face in the mesh. In this case, it makes sense to talk about the node ids of the interface $N(\Gamma(\mathcal{F})))$, its physical map $\phi^{\Gamma(\mathcal{F})}$, and its reference face $R(\Gamma(\mathcal{F}))$.
#
# By solving this problem, we are building a parametrization of the interface
# seen from each face around. This is the key ingredient needed to integrate quantities
# defined on the faces around at the interface (see section TODO),
# and to build interpolations that are conforming at the interface (see section TODO).
#
# Once the step 1. of the problem is solved and we have the local face $f$
# and the permutation matrix $P$, 
# the map $\varphi$ solution of step 2. is readily
# computed as follows
#
# $\varphi(x) = \sum_{n=1}^{\hat N} [P\hat X]_n [\hat S]_n(x),$
#
# where $\hat V$ is the reference space of the local face $f$, $\hat X$ is the vector
# of node coordinates of $\hat V$, $\hat S$ is the vector of shape functions of $\hat V$, and
# $\hat N$ is the length of $\hat X$ and $\hat S$.
# In the API, this map can be built from a mesh object `M::AbstractMesh`, the id `F::Integer` and dimension `m` of face $F$, the local face id `f::Integer` and dimension `n` of $f$, and the permutation vector `P::Vector{<:Integer}` encoding the permutation matrix $P$ as follows:
#
# ```julia
# rm = face_reference_id(M,m)[F]
# Fref = reference_domain(M,m)[rm]
# L = mesh(Fref)
# rn = face_reference_id(L,n)[f]
# V = reference_spaces(L,n)[rn]
# X = node_coordinates(V)
# S = shape_functions(V)
# N = length(X)
# φ = x -> sum(n->X[P[n]]*S[n],1:N)
# ```
#
# Thus, the only remaining problem 
# is finding the local face id of $f$
# and the permutation vector `P` representing the matrix $P$.
# In the code, this information is encoded for all possible interfaces
# in the object returned by `T=topology(M)` for a mesh
# object `M::AbstractMesh`.
#  The returned object `T` is an instance of a type
#  that implements the [`AbstractTopology`](@ref) interface.
#  We detail now the main methods in this interface.
#
# ## Face incidence
#
# An object `T::AbstractTopology` represents the *topology* of a given face complex
# $M$. We define the topology of $M$ 
# as the graph $T$, where the vertices of $T$ are the faces in $M$.
# The edges $(F_1,F_2)$ and $(F_2,F_1)$ connecting two faces $F_1,F_2\in M$ exists in the graph $T$ if
# - the two faces are the same, $F_1 = F_2$, or
# - the two faces are of different dimensions and have a non-empty interface, $\Gamma(\{F_1,F_2\}) \neq \empty$.
#
#  If two faces are connected by an edge in the graph $T$, we say that the faces are *incident*, or  *adjacent*. Note that this graph is symmetric.
#
#  In the API,
#  The graph is encoded using [adjacency lists](https://en.wikipedia.org/wiki/Adjacency_list)
#  organized according to face dimensions.
#  The adjacency list containing the edges starting at faces of dimension `m` and ending at a faces of dimension `n` is returned by function
#  `face_incidence(T,m,n)`. One recovers all edges in the graph `T` by calling this function for all possible pairs `(m,n)` with `m in 0:D` and `n in 0:D`, being `D=num_dims(M)`.
#  The result of `face_incidence(T,m,n)` is a long vector of small vectors of integers encoded via a `JaggedArray` object, since the inner vectors often have different lengths.
#
# The vector `Fns=face_incidence(T,m,n)[Fm]` contains the ids of the `n`-faces $F_n\in M$ that are incident to `m`-face $F_m$ with id `Fm`.
# For `m>n`, the face ids in `Fns` are sorted according to the id of the local faces of $F_m$,
# namely `Fn=Fns[fn]` is the id of the `n`-face $F_n\in M$ equivalent to the local
# face  $f_n\in L(R(F_m))$ with id `fn`.
# For `m<n`, the ids in `Fns`
# are arbitrarily sorted.
# For `m==n`, the vector `Fns` is equal to `[Fm]` (a vector with a single integer).
#
# From the topology object `T`, we can obtain part of the solution of the face-gluing problem above as follows.  We start with the id `Fm::Integer` of the face of dimension `m` that
# stands for the interface $F_m :=\Gamma(\mathcal{F})\in M$. The ids of the `n`-faces in $\mathcal{F}$
# are obtained as `Fns=face_incidence(T,m,n)[Fm]`. For a given id `Fn in Fns`, say that $F_n$ is the face in $\mathcal{F}$ with id `Fn`. The problem is to 
# find `fn` the id of the local face $f_n$ of $F_n$ equivalent to the interface $F_m$.
# This is done
# by getting the ids `Fms = face_incidence(T,n,m)[Fn]` of the `m`-faces on the boundary of $F_n$ and look for in which position is the id of interface $F_m$, namely 
# `fn=find(F->F==Fm,Fms)`.
#
# ## Node permutations
#
# Finally, the permutation `P` of the face-gluing problem is encoded as follows.
# Given the id `Fm::Integer` of an `m`-face in a mesh `M::AbstractMesh`, and the
# id `fn::Integer` of local `n`-face of `F`, all possible permutations for this local
# face are enumerated in the reference space of `fn`. Let `V` the reference space
# associated with local face `fn` in face `F`, all possible node permutations 
# are given in vector `node_permutations(V)`. From all these, the permutation `P`
# corresponding to the local face `fn` is obtained as `P=node_permutations(V)[k])`,
# where `k` is an index obtained from the topology object `k=face_permutation_id(T,m,n)[Fm][fn]`.
#

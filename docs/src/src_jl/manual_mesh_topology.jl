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
#     - A way of getting the local ids directly?
#
#
# The mesh topology provides the information needed to "glue" faces
# in a *face complex*. A face complex is a mesh with some
# additional properties (see definition below) that allow us 
# to formalize the meaning of "gluing" faces.
# This information is needed in many finite element (FE) methods, 
# e.g., to build (high-order) conforming FE spaces,
# or to integrate at face interfaces.
#
#
#
# ## Local faces
#
# In order to define what a face complex is, we need first to introduce
# the concept of *local faces*. To this end,
# let $F\in M$ be a physical face in mesh $M$, where $\hat\Omega:=\hat\Omega(F)$ is its reference
# domain. Let us also consider the chain $\hat C:=C(\hat \Omega)$ used to define the geometry of
# the reference domain $\hat \Omega$ and let $d$ be the dimension of $\hat \Omega$. See
# Section [Mesh geometry](@ref) for the definition of these objects.
# Let us also define the interface
# $\Gamma(\mathcal{A})$ of a set of faces $\mathcal{A}$
# as the open set such that its closure is the intersection of the closures
# of the faces in $\mathcal{A}$, namely 
# $\overline{\Gamma(\mathcal{A})} = \cap_{A\in \mathcal{A}} \bar A$.
#
# With these notations, we define the set of *local faces* of the reference domain $\hat \Omega$
# as the mesh $L(\hat \Omega)$ containing:
# 1. the reference domain $\hat \Omega$,
# 2. the faces in the chain $\hat C$,
# 3. and the interfaces $\Gamma(\{F_1,F_2\})$ of any pair of faces $F_1,F_2\in\hat C$.
#
# The set $L(\hat \Omega)$ contains all faces in the boundary of the polytope $\hat \Omega$. E.g.,
# for a reference edge, it contains the edge and the two end vertices.
# For a reference square, it contains the square, four edges, and the four vertices
# at the intersection of the edges. For a reference cube, it contains the cube, six surfaces, twelve edges, and eight vertices (see the next figure).
#
# In the code, one can access the mesh of local faces
# from a reference domain `Ωref::AbstractDomain` with `Mref = mesh(Ωref)`. The mesh
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
# One of the key properties of a face complex is that it is *conforming*.
# A mesh $M$ is conforming if for any pair of faces $F_1,F_2\in M$
# with non-empty interface $\Gamma_{1,2}:=\Gamma(\{F_1,F_2\})\neq\empty$ exists a local
# face $f_1\in L(\hat\Omega(F_1))$ and a local face  $f_2\in L(\hat\Omega(F_2))$ such that
#
# * their images via the physical map coincide with the interface, $\Gamma_{1,2}=\phi^{F_1}(f_1)=\phi^{F_2}(f_2)$, and
# * they share global node ids, $n(F_1,f_1)=Pn(F_2,f_2)$.
#
# In this definition, $P$ is a permutation matrix and $n(F,f)$ is a vector containing the node
# ids of $F$ restricted to the local face $f\in L(\hat\Omega(F))$, namely $[n(F,f)]_r:=[n(F)]_l$ with $l=[n(f)]_r$.
# Both local faces have the same global ids but these ids are allowed to be order differently in each local face.
# This is why we include the permutation matrix $P$.
#
# If these local faces exist, we say that face $F_1$ is *conforming* to face $F_2$ at the interface
# $\Gamma_{1,2}$, *via* the local face $f_1$ (idem for $F_2$). We also say that
# the local faces $f_1$ and $f_2$, and the interface $\Gamma_{1,2}$ are *equivalent*.
#
# ## Face complexes
#
# With these notations we can introduce the concept of face complex as follows.
# A mesh $M$ is a *face complex* (or a [polyhedral complex](https://en.wikipedia.org/wiki/Polyhedral_complex)) if
# * it is conforming, and
# * it contains the images of all local faces, namely $\phi^F(f)\in M$ for all $f\in L(\hat\Omega(F))$ and all $F\in M$.
#
# Given a mesh $M$ that is conforming, but not a face complex, it is always possible to add additional faces to make it a face complex. In the API, this is done
# with function `M2 = complexify(M)` from a mesh object `M::AbstractMesh`. The result
# is `M2::AbstractMesh`. In addition function `is_face_complex(M)` checks if
# `M` is a face complex. It returns true form meshes created with
# function `complexify`.
# If the chain $\hat C$ used to define the local faces $L(\hat \Omega)$ above is conforming,
# the set of local faces $L(\hat \Omega)$ is a face complex. In fact, we
# create the local faces in the library by applying `complexify` to the chain $\hat C$.
#
#
# ## Gluing faces
#
# Let us consider a face complex $M$, a face $F\in M$,
# and the set of faces $\mathcal{A}(F)$ such that their interface $\Gamma(\mathcal{A}(F))=F$  coincides with face
# $F$. If these face exist, we call $\mathcal{A}(F)$ the *faces around* $F$.
# We refer to *gluing* faces $\mathcal{A}(F)$ at the their interface $F$ as solving the following two-step
# problem for each face around $A\in\mathcal{A}(F)$:
#
# 1. Find the local face $a\in L(\hat\Omega(A))$ and the permutation matrix $P$ such that $n(F)=P n(A,a)$.
# 2. Build a map $\varphi$ such that $\phi^A(\phi^a(\varphi(x)))=\phi^{F}(x)$ for all $x\in \hat\Omega(F)$.
#
# For this problem to have a solution, the mesh $M$ needs to be conforming. This is what allows
# us to find the local face $a$, the permutation matrix $P$, and the map $\varphi$.
# In addition, mesh $M$ needs to be a face complex so that the interface $\Gamma(\mathcal{A}(F))\in M$ is also
# a face in the mesh, i.e., face $F$.
# In this case, it makes sense to talk about the node ids of the interface $n(F)$,
# its physical map $\phi^{F}$, and its reference domain $\hat\Omega(F)$.
#
# By solving this problem, we are building a common parametrization of the interface
# seen from each face around. This is the key ingredient needed to compute integrals
# at the interface of quantities defined on the faces around.
# (see section TODO),
# and to build interpolations defined on the faces around that are conforming at the interface (see section TODO).
#
# Once the step 1. of the problem is solved and we have the local face $a$
# and the permutation matrix $P$, 
# the map $\varphi$ solution of step 2. is readily
# computed as
#
# $\varphi(\xi) := \sum_{n=1}^{\hat N} [P\hat x]_n [\hat s]_n(\xi),$
#
# where $\hat x:=x(\hat V(a))$, $\hat s:=s(\hat V(a))$, and
# $\hat N$ is the length of $\hat x$ and $\hat s$.
# Thus, the only remaining  part of the problem 
# is step 1.
# In the code, the solution of this step is encoded for all possible interfaces
# in the object returned by `T=topology(M)` for a mesh
# object `M::AbstractMesh`.
#  The returned object `T` is an instance of a type
#  that implements the [`AbstractTopology`](@ref) interface.
#  We detail now the main methods in this interface.
#
# ## Face incidence
#
# Given a face $F\in M$ in a face complex $M$, we detail how to get
# the set of faces around $\mathcal{A}(F)\subset M$, and  the local face $a\in L(\hat\Omega(A))$
# that is equivalent to $F$ for each face around $A\in\mathcal{A}(F)$.
# This information is encoded in a graph $T$ called the mesh *topology*.
#
# We define the topology of $M$ 
# as the graph $T$, where the vertices of $T$ are the faces in $M$ and its edges are refined as follows:
# The edges $(F_1,F_2)$ and $(F_2,F_1)$ connecting two faces $F_1,F_2\in M$ exists in the graph $T$ if
# - the two faces are the same, $F_1 = F_2$, or
# - the two faces are of different dimensions and have a non-empty interface, $d(F_1)\neq d(F_2)$ and $\Gamma(\{F_1,F_2\}) \neq \empty$.
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
# The vector `As=face_incidence(T,m,n)[F]` has the following interpretation.
# For `m<n`, the ids in `As` are the ids of the faces around $\mathcal{A}(F)$, which are arbitrarily sorted in vector `As`.
# For `m>n`, the face ids in `As` of the n-faces on the boundary of $F$. The ids in `As`
# are sorted according to the id of the local faces of $F$,
# namely `A=As[f]` means that 
# the local
# `n`-face  $f$ of $F$  with id `f` is equivalent to the 
# `n`-face $A\in M$ with id `A`.
# For `m==n`, the vector `As` is equal to `[F]` (a vector with a single integer).
#
# In summary, given a topology object `T::AbstractTopology`, we can get the following information for the `m`-face with id `F::Integer`.
# We get the ids of the `n`-faces around `F` as `As=face_incidence(T,m,n)[F]`. For each `A in As`, we get the id `a` of the local
# face  of `A` equivalent to `F` as follows. We get all `m`-faces at the boundary of `A` with `Bs=face_incidence(T,n,m)[A]`, and we
# find in which position in `Bs` the id `F` is located, namely `a=find(B->B==F,Bs)`.
#
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

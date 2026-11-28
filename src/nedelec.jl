# Nédélec spaces of the first kind (H(curl)-conforming) of any order on simplices.
#
# Reference element of order k (k = 0 is the lowest order), on the reference
# triangle or tetrahedron:
#   NED_k = P_k^D ⊕ S_{k+1},  S_{k+1} = {p ∈ homogeneous P_{k+1}^D : p⋅x = 0}.
# (The primal basis is written in powers of x - c, c the centroid: the same space,
# better conditioned at higher order.)
# Degrees of freedom, per entity, with the entity's vertices in ascending order. Each
# is a mean over the entity's parameter domain against a weight from a basis that is
# orthonormal for that mean:
#   edge (a,b), t = x_b - x_a:        ∫₀¹ (v⋅t) ℓⱼ(s) ds, ℓⱼ(s) = √(2j+1) Lⱼ(2s-1) (Legendre), j = 0..k
#   face (a,b,c) (3D, k ≥ 1):         2 ∫ (v⋅tᵢ) qᵣ(u,w) dÂ over the unit triangle, x = x_a + u t₁ + w t₂,
#                                     tᵢ ∈ {t₁,t₂} = {x_b-x_a, x_c-x_a}, qᵣ of degree ≤ k-1
#   interior (2D k ≥ 1, 3D k ≥ 2):    |K̂|⁻¹ ∫ v_d qᵣ dx̂, d = 1..D, qᵣ of degree ≤ k-1 (2D) or ≤ k-2 (3D)
# The qᵣ come from Gram–Schmidt, in order of degree, of products of shifted Legendre
# polynomials (see nedelec_orthonormal_values). The weights span the same spaces as
# monomial weights, so the element is the same; only the basis changes. With monomial
# face and interior weights and plain integrals, cond(K+M) (K the curl-curl and M the
# mass matrix) on the reference element at k = 1 to 3 was 4.3 to about 2400 times that
# of Gridap's Nédélec element; with these weights it is 0.45 to 2.6 times.
# Face weights depend only on (u,w), i.e. on the face's vertex order, so every cell
# sharing a face sees the same moments.
#
# Mesh space and orientation: each cell is given a "sorted frame", the reference
# map whose i-th reference vertex goes to the cell vertex with the i-th smallest
# global node ID. In that frame every shared edge and face is seen with the same
# vertex order from all cells that contain it, so the edge and face moments define
# global DOFs with no signs or permutations. With φ the cell's geometry map and
# R the affine map of the reference simplex such that φ∘R is the sorted frame,
# the cell's shape function in local slot s is
#   φ_s(x) = J⁻ᵀ ψᴾ_s(ξ),  ψᴾ_s(ξ) = A_R⁻ᵀ φ̂_{perm(s)}(R⁻¹ξ),  ξ = φ⁻¹(x),
# where perm(s) is the reference DOF of the sorted-frame entity that R maps to the
# slot's entity. NED_k is invariant under such maps, so ψᴾ_s = Σⱼ Cᴾ[s,j] φ̂ⱼ with
# a matrix Cᴾ precomputed for each of the (D+1)! vertex permutations P. At order
# 0, Cᴾ is diagonal with entries ±1 (the classical edge signs). Jacobians transform
# as J⁻ᵀ ∇̂ψᴾ J⁻¹, which is exact for the affine (linear simplex) cells accepted.
# On a face of dimension D-1 (traces), the same formulas hold with ξ the face
# quadrature point mapped into the selected incident cell and J that cell's D x D
# Jacobian; Cᴾ is the cell's matrix. The reference values are tabulated per local
# face and face permutation.

"""
    nedelec_space(domain::AbstractFaceDomain,order)
    nedelec_space(domain::AbstractMeshDomain,order;dirichlet_boundary=nothing)

Return the Nédélec space of the first kind (H(curl)-conforming) of the given `order` on `domain`.

`order=0` is the lowest-order (Whitney) element with one DOF per edge; order `k` contains the polynomials of degree
`k` and has edge, face and interior moments. Implemented on triangles and tetrahedra. The mesh version requires linear
simplex geometry. Its shape functions use the covariant Piola map, and each cell is oriented by the global node IDs
of its vertices (see [`nedelec_cell_orientation`](@ref)). A `dirichlet_boundary` fixes the tangential trace: the DOFs
of the edges and faces in the closure of the given boundary faces become Dirichlet DOFs.

Shape functions support `value` and `ForwardDiff.jacobian` (from which the curl follows). Form arguments and
discrete fields can be integrated on cells and on faces of dimension `D-1` (boundary and skeleton domains). On a
face, the full vector of the selected incident cell is evaluated: the cell is given by the `faces_around` of the
domain (on a boundary domain made of interior faces, `boundary(mesh;...)` takes the first incident cell unless
`faces_around` is given) or by the side index on a skeleton. Lower-dimensional domains raise an `ArgumentError`.

# Level

Basic
"""
function nedelec_space end

struct CovariantPiolaMap end

function nedelec_check_order(order)
    (order isa Integer && order >= 0) || throw(ArgumentError(
        "The order of a Nédélec space must be a non-negative integer, got order=$order"))
    nothing
end

function nedelec_check_geometry(geometry)
    is_simplex(geometry) || throw(ArgumentError("Nédélec spaces are only implemented on simplices"))
    num_dims(geometry) in (2,3) || throw(ArgumentError("Nédélec spaces are implemented in 2 and 3 dimensions"))
    nothing
end

function nedelec_space(domain::AbstractFaceDomain,order)
    nedelec_check_order(order)
    nedelec_check_geometry(domain)
    workspace = nedelec_setup(domain,Int(order))
    NedelecFaceSpace(domain,order,workspace)
end

struct NedelecFaceSpace{A,B,C} <: AbstractFaceSpace
    domain::A
    order::B
    workspace::C
end

domain(fe::NedelecFaceSpace) = fe.domain
order(fe::NedelecFaceSpace) = fe.order
num_dofs(fe::NedelecFaceSpace) = length(fe.workspace.dual_basis)
primal_basis(fe::NedelecFaceSpace) = fe.workspace.primal_basis
dual_basis(fe::NedelecFaceSpace) = fe.workspace.dual_basis
shape_functions(fe::NedelecFaceSpace) = fe.workspace.shape_functions

function tabulator(fe::NedelecFaceSpace)
    primal_t = permutedims(fe.workspace.primal_basis)
    B = fe.workspace.change_of_basis
    (f,x) -> broadcast(f,primal_t,x)*B
end

# Exponent tuples of the monomials of total degree ≤ k (or == k) in D variables.
function nedelec_exponents(D,k;homogeneous=false)
    k < 0 && return NTuple{D,Int}[]
    all = [Tuple(i) for i in CartesianIndices(ntuple(_->0:k,D))][:]
    filter(e -> homogeneous ? sum(e) == k : sum(e) <= k, all)
end

@inline nedelec_monomial(e,x) = prod(ntuple(i->x[i]^e[i],length(e)))

function nedelec_legendre(j,t)
    p0,p1 = one(t),t
    j == 0 && return p0
    for n in 1:j-1
        p0,p1 = p1,((2n+1)*t*p1-n*p0)/(n+1)
    end
    p1
end

# Values at the quadrature points xs (weights ws) of a basis of the polynomials of
# degree ≤ n in dim variables, orthonormal for the mean over the quadrature's domain:
# Σ_q ws[q] p_r(xs[q]) p_s(xs[q]) / Σ_q ws[q] = δ_rs. Gram–Schmidt, in order of degree,
# of products of Legendre polynomials in 2x_i - 1, done by a QR factorization of the
# weighted values (R with a positive diagonal is the Cholesky factor of the Gram matrix,
# without forming that matrix). The quadrature must integrate degree 2n exactly.
function nedelec_orthonormal_values(dim,n,xs,ws)
    exponents = sort(nedelec_exponents(dim,n),by=sum)
    P = [prod(ntuple(i->nedelec_legendre(e[i],2x[i]-1),dim)) for x in xs, e in exponents]
    sqrtW = Diagonal(sqrt.(collect(ws) ./ sum(ws)))
    R = qr(sqrtW*P).R
    R = Diagonal(sign.(diag(R)))*R
    P/UpperTriangular(R)
end

function nedelec_unit_vector(::Val{D},d,m) where D
    SVector{D}(ntuple(i -> i == d ? m : zero(m),Val(D)))
end

# Written in y = x - c: P_k and S_{k+1}(y) span NED_k for any shift c, because
# S_{k+1}(x - c) - S_{k+1}(x) ⊂ P_k.
function nedelec_primal_basis(::Val{D},k,c) where D
    polynomials = Function[]
    for d in 1:D, e in nedelec_exponents(D,k)
        push!(polynomials,x->nedelec_unit_vector(Val(D),d,nedelec_monomial(e,x-c)))
    end
    rotations = Function[]
    if D == 2
        # S_{k+1} = {q(y)(-y₂,y₁) : q homogeneous of degree k}.
        for e in nedelec_exponents(2,k;homogeneous=true)
            push!(rotations,x->(y = x-c; nedelec_monomial(e,y)*SVector(-y[2],y[1])))
        end
    else
        # S_{k+1} is spanned by x × (m e_d), m homogeneous of degree k; keep a
        # linearly independent subset (integer coefficients over the monomials).
        hom = nedelec_exponents(3,k+1;homogeneous=true)
        position = Dict((c,e)=>i for (i,(c,e)) in enumerate(Iterators.product(1:3,hom)))
        shift(e,j) = ntuple(i->e[i]+(i==j),3)
        selected = zeros(Float64,0,length(position))
        for e in nedelec_exponents(3,k;homogeneous=true), d in 1:3
            coeffs = zeros(Float64,length(position))
            # x × (m e_d) = (x2 v3 - x3 v2, x3 v1 - x1 v3, x1 v2 - x2 v1) with v = m e_d.
            terms = d == 1 ? ((2,3,1),(3,2,-1)) : d == 2 ? ((1,3,-1),(3,1,1)) : ((1,2,1),(2,1,-1))
            for (component,variable,sign) in terms
                coeffs[position[(component,shift(e,variable))]] += sign
            end
            candidate = vcat(selected,coeffs')
            if rank(candidate) > size(selected,1)
                selected = candidate
                push!(rotations,x->(y = x-c; cross(y,nedelec_unit_vector(Val(3),d,nedelec_monomial(e,y)))))
            end
        end
    end
    @assert length(rotations) == (D == 2 ? k+1 : (k+1)*(k+3))
    vcat(polynomials,rotations)
end

function nedelec_setup(domain,k)
    D = num_dims(domain)
    cell_mesh = mesh(domain)
    node_x = node_coordinates(cell_mesh)
    @assert length(node_x) == D+1 # the reference simplex mesh has its vertices as nodes
    # Entities by dimension, each with its vertices in ascending order.
    entity_nodes = [map(nodes->sort(Int.(collect(nodes))),face_nodes(cell_mesh,d)) for d in 0:D]
    dual_basis = Function[]
    dof_entity = Tuple{Int,Int}[] # dof -> (dimension, entity)
    # Each moment is Σ_q ω[q] f(x_q): ω folds the quadrature weight, the division by the
    # entity's parametric measure and the weight polynomial (see the top of this file).
    # Edge moments.
    s,w = gausslegendre(k+3)
    ss = (s .+ 1) ./ 2
    edge_ω = [w[q]/2*sqrt(2j+1)*nedelec_legendre(j,s[q]) for q in eachindex(s), j in 0:k]
    for (edge,(a,b)) in enumerate(entity_nodes[2]), j in 0:k
        xa,t = node_x[a],node_x[b]-node_x[a]
        ω = edge_ω[:,j+1]
        push!(dual_basis,v -> sum(q -> ω[q]*dot(t,v(xa+ss[q]*t)),1:length(ss)))
        push!(dof_entity,(1,edge))
    end
    # Face moments (3D).
    if D == 3 && k >= 1
        face_quadrature = quadrature(unit_simplex(Val(2)),2k+3)
        uw,wf = coordinates(face_quadrature),weights(face_quadrature)
        qf = nedelec_orthonormal_values(2,k-1,uw,wf)
        face_ω = [wf[q]/sum(wf)*qf[q,r] for q in eachindex(wf), r in axes(qf,2)]
        for (face,(a,b,c)) in enumerate(entity_nodes[3]), r in axes(face_ω,2), i in 1:2
            xa = node_x[a]
            t1,t2 = node_x[b]-xa,node_x[c]-xa
            t = i == 1 ? t1 : t2
            ω = face_ω[:,r]
            push!(dual_basis,v -> sum(q -> ω[q]*dot(t,v(xa+uw[q][1]*t1+uw[q][2]*t2)),1:length(wf)))
            push!(dof_entity,(2,face))
        end
    end
    # Interior moments.
    interior_degree = D == 2 ? k-1 : k-2
    if interior_degree >= 0
        cell_quadrature = quadrature(domain,2k+3)
        xs,wc = coordinates(cell_quadrature),weights(cell_quadrature)
        qc = nedelec_orthonormal_values(D,interior_degree,xs,wc)
        cell_ω = [wc[q]/sum(wc)*qc[q,r] for q in eachindex(wc), r in axes(qc,2)]
        for r in axes(cell_ω,2), d in 1:D
            ω = cell_ω[:,r]
            push!(dual_basis,v -> sum(q -> ω[q]*v(xs[q])[d],1:length(wc)))
            push!(dof_entity,(D,1))
        end
    end
    centroid = sum(node_x)/length(node_x)
    primal_basis = nedelec_primal_basis(Val(D),k,centroid)
    ndofs = length(dual_basis)
    expected = D == 2 ? (k+1)*(k+3) : div((k+1)*(k+3)*(k+4),2)
    @assert length(primal_basis) == ndofs == expected
    # Shape functions dual to the moments.
    A = [dual_basis[i](primal_basis[j]) for i in 1:ndofs, j in 1:ndofs]
    B = inv(A)
    shape_functions = map(1:ndofs) do i
        x -> sum(j -> B[j,i]*primal_basis[j](x),1:ndofs)
    end
    entity_dofs = [[Int32[] for _ in entity_nodes[d+1]] for d in 0:D]
    for (dof,(d,entity)) in enumerate(dof_entity)
        push!(entity_dofs[d+1][entity],dof)
    end
    reference_topology = topology(domain)
    workspace = (;entity_nodes,entity_dofs,dof_entity,primal_basis,dual_basis,shape_functions,
        change_of_basis=B,reference_topology)
    orientation = nedelec_orientation_matrices(D,node_x,workspace)
    (;workspace...,orientation...)
end

# The matrices Cᴾ for every permutation P of the reference vertices (see the top of this file).
function nedelec_orientation_matrices(D,node_x,ws)
    (;entity_nodes,entity_dofs,dof_entity,dual_basis,primal_basis,change_of_basis) = ws
    ndofs = length(dual_basis)
    perms = [collect(p) for p in permutations(1:D+1)]
    entity_index = [Dict(nodes=>entity for (entity,nodes) in enumerate(entity_nodes[d+1])) for d in 0:D]
    own_index = zeros(Int,ndofs)
    for d in 0:D, dofs in entity_dofs[d+1], (i,dof) in enumerate(dofs)
        own_index[dof] = i
    end
    matrices = zeros(Float64,ndofs,ndofs,length(perms))
    vertices = [SVector{D,Float64}(x) for x in node_x]
    for (P,perm) in enumerate(perms)
        inverse = invperm(perm)
        b = vertices[perm[1]]
        A = hcat(ntuple(i->vertices[perm[i+1]]-b,Val(D))...)
        Ainv = inv(SMatrix{D,D}(A))
        # G[j,i] = σ̂ⱼ(A⁻ᵀ pᵢ∘R⁻¹); the push-forward of φ̂_src = Σᵢ B[i,src] pᵢ has moments G*B[:,src].
        G = [dual_basis[j](ξ -> transpose(Ainv)*primal_basis[i](Ainv*(ξ-b))) for j in 1:ndofs, i in 1:ndofs]
        GB = G*change_of_basis
        for slot in 1:ndofs
            d,entity = dof_entity[slot]
            sorted_nodes = sort([inverse[v] for v in entity_nodes[d+1][entity]])
            source = entity_dofs[d+1][entity_index[d+1][sorted_nodes]][own_index[slot]]
            matrices[slot,:,P] = GB[:,source]
        end
    end
    (;orientation_permutations=perms,orientation_matrices=matrices)
end

function nedelec_num_face_permutations(fe::NedelecFaceSpace,d)
    topo = fe.workspace.reference_topology
    face_to_rid = face_reference_id(topo,d)
    rid_to_topo = reference_topologies(topo,d)
    map(face_to_rid) do rid
        length(vertex_permutations(rid_to_topo[rid]))
    end
end

face_own_dofs(fe::NedelecFaceSpace,d) = fe.workspace.entity_dofs[d+1]

function face_dofs(fe::NedelecFaceSpace,d)
    # Closure: the DOFs of every entity whose vertices belong to the d-face.
    (;entity_nodes,dof_entity) = fe.workspace
    map(entity_nodes[d+1]) do nodes
        Int32[dof for (dof,(dd,entity)) in enumerate(dof_entity) if dd <= d && all(in(nodes),entity_nodes[dd+1][entity])]
    end
end

function face_own_dof_permutations(fe::NedelecFaceSpace,d)
    # Each cell is evaluated in its sorted frame, where a shared entity's own DOFs are
    # ordered the same way from every cell, so no permutation is needed.
    own_dofs = face_own_dofs(fe,d)
    nperms = nedelec_num_face_permutations(fe,d)
    map(own_dofs,nperms) do dofs,n
        [collect(Int32,1:length(dofs)) for _ in 1:n]
    end
end

# Oriented reference shape function ψᴾ_slot = Σⱼ Cᴾ[slot,j] φ̂ⱼ.
function nedelec_oriented_shape_function(fe::NedelecFaceSpace,slot,P)
    C = fe.workspace.orientation_matrices
    φ̂ = fe.workspace.shape_functions
    ξ -> sum(j -> C[slot,j,P]*φ̂[j](ξ),1:length(φ̂))
end

function nedelec_space(domain::AbstractMeshDomain,order;dirichlet_boundary=nothing)
    nedelec_check_order(order)
    mesh = GT.mesh(domain)
    D = num_dims(domain)
    num_ambient_dims(mesh) == D || throw(ArgumentError("Nédélec spaces require full-dimensional geometry"))
    ctype_to_refface = reference_spaces(mesh,D)
    for refface in ctype_to_refface
        geometry = GT.domain(refface)
        nedelec_check_geometry(geometry)
        num_dofs(refface) == D+1 || throw(ArgumentError(
            "Nédélec spaces require linear (affine) simplex geometry"))
    end
    ctype_to_reffe = map(refface->nedelec_space(GT.domain(refface),order),ctype_to_refface)
    orientation = nedelec_cell_orientation(mesh,D,ctype_to_reffe)
    workspace = nothing
    NedelecMeshSpace(domain,order,dirichlet_boundary,ctype_to_reffe,orientation,workspace) |> setup_space
end

struct NedelecMeshSpace{A,B,C,E,S,W} <: AbstractSpace
    domain::A
    order::B
    dirichlet_boundary::C
    reference_spaces::E
    cell_orientation::S # Vector{Int8}: [mesh_cell] -> index of the vertex permutation of its sorted frame
    workspace::W
end

domain(a::NedelecMeshSpace) = a.domain
order(a::NedelecMeshSpace) = a.order
continuous(a::NedelecMeshSpace) = true
dirichlet_boundary(a::NedelecMeshSpace) = a.dirichlet_boundary
reference_spaces(a::NedelecMeshSpace) = a.reference_spaces
workspace(a::NedelecMeshSpace) = a.workspace

function face_reference_id(a::NedelecMeshSpace)
    D = num_dims(a.domain)
    face_reference_id(GT.mesh(a.domain),D)
end

function replace_workspace(a::NedelecMeshSpace,workspace)
    NedelecMeshSpace(a.domain,a.order,a.dirichlet_boundary,a.reference_spaces,a.cell_orientation,workspace)
end

"""
    nedelec_cell_orientation(space)

Return the vector `orientation[mesh_cell]` with, for each cell, the index of the permutation of the reference vertices
that sorts the cell's vertices by global node ID. The cell's shape functions are evaluated in that sorted frame, which
makes every shared edge and face oriented consistently (for the lowest order, edges point from the lower to the higher
global node ID).

# Level

Advanced
"""
nedelec_cell_orientation(a::NedelecMeshSpace) = a.cell_orientation

function nedelec_cell_orientation(mesh,D,ctype_to_reffe)
    Dface_to_nodes = face_nodes(mesh,D)
    Dface_to_ctype = face_reference_id(mesh,D)
    lookup = map(ctype_to_reffe) do fe
        Dict(perm=>Int8(P) for (P,perm) in enumerate(fe.workspace.orientation_permutations))
    end
    map(1:length(Dface_to_nodes)) do Dface
        lookup[Dface_to_ctype[Dface]][sortperm(collect(Dface_to_nodes[Dface]))]
    end
end

# Code generation for the compiler (serial, threaded and CUDA execution).

shape_function_map(::NedelecMeshSpace) = CovariantPiolaMap()

# The generated code below evaluates on cells only. Traces (faces of dimension
# D-1) take the accessor path (nedelec_shape_function_accessor_physical), which
# tabulates the incident cell's shape functions at the face points.
function nedelec_check_integration_domain(V,Q)
    if num_dims(GT.domain(Q)) != num_dims(GT.domain(V))
        error("Nédélec cell code generation reached with a quadrature of dimension $(num_dims(GT.domain(Q)))")
    end
    nothing
end

function check_integration_domain(space::NedelecMeshSpace,domain::AbstractDomain)
    D = num_dims(GT.domain(space))
    d = num_dims(domain)
    if d < D-1
        throw(ArgumentError(
            "Nédélec shape functions can be integrated on cells and on faces of dimension D-1 only (D=$D, got $d)"))
    end
    nothing
end

"""
    nedelec_reference_table(f,V,Q)

Tabulate `f` (`value` or `ForwardDiff.jacobian`) of the oriented reference shape functions of the Nédélec space `V`
at the reference points of the quadrature `Q`. The result `table[slot,point,P]` holds one table per vertex
permutation `P` of the reference cell (see [`nedelec_cell_orientation`](@ref)).

# Level

Advanced
"""
function nedelec_reference_table(f,V::NedelecMeshSpace,Q::AbstractQuadrature)
    mesh = GT.mesh(GT.domain(Q))
    d = num_dims(GT.domain(Q))
    rids = face_reference_id(mesh,d)
    rid = first(rids)
    @assert all(==(rid),rids) # one cell type, as in space_integration_data
    points = coordinates(reference_quadratures(Q)[rid])
    fe = reference_spaces(V)[rid]
    T0 = tabulator(fe)(f,points) # [point,dof] in the reference frame
    C = fe.workspace.orientation_matrices
    ndofs,npoints,nperms = size(C,1),length(points),size(C,3)
    table = Array{eltype(T0),3}(undef,ndofs,npoints,nperms)
    for P in 1:nperms, q in 1:npoints, slot in 1:ndofs
        table[slot,q,P] = sum(j -> C[slot,j,P]*T0[q,j],1:ndofs)
    end
    table
end

function mapped_shape_function_expr(::CovariantPiolaMap,func,mesh,V,Q,mesh_expr,V_expr,Q_expr,f_expr,q_expr,i_expr)
    nedelec_check_integration_domain(V,Q)
    orientation_expr = :(Int(nedelec_cell_orientation($V_expr)[$f_expr]))
    is_tabulated = func === value || func === ForwardDiff.jacobian
    if !is_tabulated
        func === ForwardDiff.gradient && throw(ArgumentError(
            "Nédélec shape functions are vector-valued; use ForwardDiff.jacobian instead of ForwardDiff.gradient"))
        d = num_dims(GT.domain(V))
        rid = first(face_reference_id(Q))
        xi0 = zero(eltype(coordinates(reference_quadratures(Q)[rid])))
        shape_expr = :(nedelec_oriented_shape_function(reference_spaces($V_expr)[$rid],$i_expr,$orientation_expr))
        geometry_functions_expr = :(shape_functions(reference_spaces($mesh_expr, Val($d))[$rid]))
        mesh_data_expr = :(mesh_integration_data($mesh_expr, $Q_expr))
        node_x_expr = :(node_coordinates($mesh_data_expr))
        nodes_expr = :(face_nodes($mesh_data_expr)[$f_expr])
        geometry_space_expr = :(mesh_space($mesh_expr, Val($d)))
        geometry_data_expr = :(space_integration_data(value, $geometry_space_expr, $Q_expr))
        geometry_values_expr = :(reference_tabulator($geometry_data_expr))
        return :(evaluate_physical_covariant_shape($func, $shape_expr, $geometry_functions_expr,
            $node_x_expr, $nodes_expr, $geometry_values_expr, $q_expr, $xi0))
    end
    v_expr = :(nedelec_reference_table($func, $V_expr, $Q_expr)[$i_expr, $q_expr, $orientation_expr])
    J_expr = jacobian_expr(mesh, Q, mesh_expr, Q_expr, f_expr, q_expr)
    invJ_expr = :(adjugate($J_expr) * inv(det($J_expr)))
    if func === value
        :(transpose($invJ_expr) * $v_expr)
    else
        :(transpose($invJ_expr) * $v_expr * $invJ_expr)
    end
end

function mapped_shape_function_expr(::CovariantPiolaMap,V,mesh_expr,V_expr,f_expr,i_expr)
    rid_expr = :(face_reference_id($V_expr)[$f_expr])
    orientation_expr = :(Int(nedelec_cell_orientation($V_expr)[$f_expr]))
    shape_expr = :(nedelec_oriented_shape_function(reference_spaces($V_expr)[$rid_expr],$i_expr,$orientation_expr))
    d = num_dims(GT.domain(V))
    geometry_rid_expr = :(face_reference_id($mesh_expr, $d)[$f_expr])
    geometry_space_expr = :(reference_spaces($mesh_expr, Val($d))[$geometry_rid_expr])
    functions_expr = :(shape_functions($geometry_space_expr))
    node_x_expr = :(node_coordinates($mesh_expr))
    nodes_expr = :(face_nodes($mesh_expr, $d)[$f_expr])
    xi0_expr = :(zero(first(node_coordinates($geometry_space_expr))))
    :(physical_covariant_shape_function($shape_expr, $functions_expr, $node_x_expr, $nodes_expr, $xi0_expr))
end

function physical_covariant_shape_function(s_ref, geometry_functions, node_x, nodes, xi0)
    function phi(xi)
        physical_map_value(xi, geometry_functions, node_x, nodes)
    end
    affine_covariant_shape_function(s_ref, phi, xi0)
end

# Nédélec spaces only accept affine cells, so the geometry map is inverted in
# closed form, ξ = ξ0 + J⁻¹(x - φ(ξ0)). This is exact and, unlike a Newton
# inverse, can be differentiated with ForwardDiff.
function affine_covariant_shape_function(s_ref, phi, xi0)
    J = ForwardDiff.jacobian(phi, xi0)
    x0 = phi(xi0)
    x -> begin
        xi = xi0 + J \ (x - x0)
        transpose(J) \ s_ref(xi)
    end
end

function evaluate_physical_covariant_shape(func, s_ref, geometry_functions, node_x,
        nodes, reference_values, point, xi0)
    s_phys = physical_covariant_shape_function(s_ref, geometry_functions, node_x, nodes, xi0)
    x_q = physical_point_value(node_x, nodes, reference_values, point)
    func(s_phys, x_q)
end

# Accessors (discrete fields used as data, sampling).

# Per cell (mesh face ID), per local DOF: (dof_v,invJ) -> transform(Σⱼ Cᴾ[dof,j] dof_v(j), invJ),
# where dof_v gives the reference-frame values of all DOFs at the current point.
function nedelec_shape_function_accessor_modifier(space::NedelecMeshSpace,transform)
    orientation = nedelec_cell_orientation(space)
    Dface_to_ctype = face_reference_id(space)
    ctype_to_C = map(fe->fe.workspace.orientation_matrices,reference_spaces(space))
    C0 = first(ctype_to_C)
    prototype = (dof_v,invJ) -> transform(C0[1,1,1]*dof_v(1),invJ)
    function face_dof_modifier(face,face_around=nothing)
        C = ctype_to_C[Dface_to_ctype[face]]
        P = Int(orientation[face])
        ndofs = size(C,1)
        function dof_modifier(dof)
            (dof_v,invJ) -> transform(sum(j -> C[dof,j,P]*dof_v(j),1:ndofs),invJ)
        end
    end
    accessor(face_dof_modifier,prototype)
end

function shape_function_accessor_modifier(f::typeof(value),space::NedelecMeshSpace)
    nedelec_shape_function_accessor_modifier(space,(v,invJ)->transpose(invJ)*v)
end

function shape_function_accessor_modifier(f::typeof(ForwardDiff.jacobian),space::NedelecMeshSpace)
    nedelec_shape_function_accessor_modifier(space,(v,invJ)->transpose(invJ)*v*invJ)
end

function shape_function_accessor_modifier(f::typeof(ForwardDiff.gradient),space::NedelecMeshSpace)
    throw(ArgumentError(
        "Nédélec shape functions are vector-valued; use ForwardDiff.jacobian instead of ForwardDiff.gradient"))
end

# The generic tabulated accessor passes no Jacobian for `value` and maps each DOF on
# its own; the Nédélec map needs the Jacobian and combines the DOFs of a cell.
for T in (:value,:(ForwardDiff.jacobian))
    for B in (false,true)
        @eval function shape_function_accessor_physical(f::typeof($T),space::NedelecMeshSpace,
                measure::AbstractQuadrature,disable_tabulation::Val{$B})
            nedelec_shape_function_accessor_physical(f,space,measure)
        end
    end
end

# On a face of dimension D-1 the reference values are those of the incident cell
# (selected by the domain's faces_around, or by the side on a skeleton) at the
# face points mapped into that cell, and J is the cell's D x D Jacobian there.
function nedelec_shape_function_accessor_physical(f,space,measure)
    check_integration_domain(space,GT.domain(measure))
    face_point_dof_v = shape_function_accessor_reference(f,space,measure)
    dface_to_modif = shape_function_accessor_modifier(f,space,GT.domain(measure))
    D = num_dims(GT.domain(space))
    face_point_Dphi = jacobian_accessor(measure,Val(D))
    J0 = GT.prototype(face_point_Dphi)
    v0 = GT.prototype(face_point_dof_v)
    prototype = GT.prototype(dface_to_modif)(j->v0,inv(J0))
    function face_point_dof_s(face,face_around=nothing)
        point_dof_v = face_point_dof_v(face,face_around)
        dof_modif = dface_to_modif(face,face_around)
        point_Dphi = face_point_Dphi(face,face_around)
        function point_dof_s(point,J=nothing)
            J2 = J === nothing ? point_Dphi(point) : J
            invJ = inv(J2)
            dof_v = point_dof_v(point)
            function dof_s(dof)
                dof_modif(dof)(dof_v,invJ)
            end
        end
    end
    accessor(face_point_dof_s,prototype)
end

# Untabulated physical shape functions: x -> J⁻ᵀ ψᴾ(ξ) with ξ = φ⁻¹(x) (affine cells).
function shape_function_accessor_physical(f,space::NedelecMeshSpace)
    domain = GT.domain(space)
    D = num_dims(domain)
    mesh = GT.mesh(domain)
    face_phi = physical_map_accessor(value,mesh,Val(D))
    orientation = nedelec_cell_orientation(space)
    Dface_to_ctype = face_reference_id(space)
    ctype_to_reffe = reference_spaces(space)
    x0 = zero(SVector{D,Float64})
    s0 = first(shape_functions(first(ctype_to_reffe)))
    prototype = q->f(x->s0(x),q)
    function face_dof_s_phys(face,face_around=nothing)
        phi = face_phi(face)
        fe = ctype_to_reffe[Dface_to_ctype[face]]
        P = Int(orientation[face])
        function dof_s_phys(dof)
            s_phys = affine_covariant_shape_function(nedelec_oriented_shape_function(fe,dof,P),phi,x0)
            x->f(s_phys,x)
        end
    end
    accessor(face_dof_s_phys,prototype)
end

# Space and field accessors (space_accessor, field_accessor). The generic loop maps
# each DOF on its own (map_shape_function); here all DOFs of the cell are combined.
function shape_functions(f,a::SpaceAccessor{AtInterior,<:NedelecMeshSpace})
    (;space,mesh_accessor,reference_space_accessor) = a
    dof_sref = GT.shape_functions(f,reference_space_accessor)
    nedelec_map_shape_functions!(workspace(f,a),f,space,mesh_accessor,dof_sref)
end

function shape_functions(f,a::SpaceAccessor{AtSkeleton,<:NedelecMeshSpace})
    (;space,mesh_accessor,reference_space_accessor) = a
    face_around = mesh_accessor.space_accessor.location.face_around
    dof_sref = GT.shape_functions(f,reference_space_accessor)
    nedelec_map_shape_functions!(workspace(f,a)[face_around],f,space,mesh_accessor,dof_sref)
end

function nedelec_map_shape_functions!(dof_sphys,f,space,mesh_accessor,dof_sref)
    cell = mesh_accessor.space_accessor.location.Dface
    fe = reference_spaces(space)[face_reference_id(space)[cell]]
    C = fe.workspace.orientation_matrices
    P = Int(nedelec_cell_orientation(space)[cell])
    J = coordinate(ForwardDiff.jacobian,mesh_accessor)
    ndofs = length(dof_sref)
    for dof in 1:ndofs
        s = sum(j -> C[dof,j,P]*dof_sref[j],1:ndofs)
        dof_sphys[dof] = nedelec_covariant_map(f,J,s)
    end
    view(dof_sphys,1:ndofs)
end

nedelec_covariant_map(::typeof(GT.value),J,s) = transpose(J)\s
nedelec_covariant_map(::typeof(ForwardDiff.jacobian),J,s) = (transpose(J)\s)/J

# The generic accessors call map_shape_function with dof = AnyIndex() only to size
# their workspaces; a concrete DOF cannot be mapped on its own.
for T in (:(GT.value),:(ForwardDiff.jacobian))
    @eval function map_shape_function(f::typeof($T),space::NedelecMeshSpace,dof,mesh_accessor,sref)
        dof isa AnyIndex || throw(ArgumentError(
            "Nédélec shape functions of a cell are mapped together, not one DOF at a time"))
        nedelec_covariant_map(f,coordinate(ForwardDiff.jacobian,mesh_accessor),sref)
    end
end

function map_shape_function(f::typeof(ForwardDiff.gradient),space::NedelecMeshSpace,dof,mesh_accessor,sref)
    throw(ArgumentError(
        "Nédélec shape functions are vector-valued; use ForwardDiff.jacobian instead of ForwardDiff.gradient"))
end

# Interpolation: in each cell's sorted frame, the DOF values are the reference
# moments of the covariant pull-back û(ξ) = J_sᵀ u(φ_s(ξ)). Shared edges and faces
# get the same values from every cell.

function interpolate_impl!(
    f::AnalyticalField,
    u::DiscreteField,
    space::NedelecMeshSpace,
    free_or_diri::FreeOrDirichlet;location=nothing)
    is_piecewise(f) && throw(ArgumentError("Piecewise interpolation into Nédélec spaces is not implemented"))
    fun = f.definition
    free_vals = GT.free_values(u)
    diri_vals = GT.dirichlet_values(u)
    domain = GT.domain(space)
    @assert is_physical_domain(domain)
    D = num_dims(domain)
    mesh = GT.mesh(domain)
    node_x = node_coordinates(mesh)
    Dface_to_nodes = face_nodes(mesh,D)
    Dface_to_dofs = face_dofs(space)
    Dface_to_ctype = face_reference_id(space)
    ctype_to_reffe = reference_spaces(space)
    orientation = nedelec_cell_orientation(space)
    dirichlet_dof_location = GT.dirichlet_dof_location(space)
    for Dface in faces(domain)
        fe = ctype_to_reffe[Dface_to_ctype[Dface]]
        (;dual_basis,entity_nodes,entity_dofs,dof_entity,orientation_permutations) = fe.workspace
        perm = orientation_permutations[orientation[Dface]]
        nodes = Dface_to_nodes[Dface]
        # Sorted frame: φ_s(ξ) = x₁ + J_s ξ, with the vertices in ascending global node ID.
        xs = [node_x[nodes[perm[i]]] for i in 1:D+1]
        Js = hcat(ntuple(i->xs[i+1]-xs[1],Val(D))...)
        û = ξ -> transpose(Js)*fun(xs[1]+Js*ξ)
        inverse = invperm(perm)
        dofs = Dface_to_dofs[Dface]
        for (slot,(d,entity)) in enumerate(dof_entity)
            dof = dofs[slot]
            if dof > 0
                free_or_diri == DIRICHLET && continue
            else
                free_or_diri == FREE && continue
                # A piece of a PiecewiseField (location given) writes only the DOFs of its own
                # boundary piece; an ordinary field (location === nothing) writes all of them.
                location === nothing || dirichlet_dof_location[-dof] == location || continue
            end
            # The sorted-frame DOF that corresponds to this slot (as in nedelec_orientation_matrices).
            sorted_nodes = sort([inverse[v] for v in entity_nodes[d+1][entity]])
            sorted_entity = findfirst(==(sorted_nodes),entity_nodes[d+1])
            own = findfirst(==(slot),entity_dofs[d+1][entity])
            v = dual_basis[entity_dofs[d+1][sorted_entity][own]](û)
            if dof > 0
                free_vals[dof] = v
            else
                diri_vals[-dof] = v
            end
        end
    end
    u
end

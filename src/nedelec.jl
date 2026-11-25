# Nédélec spaces of the first kind (H(curl)-conforming), lowest order, on simplices.
#
# Reference element: one DOF per edge, the tangential moment
#   σ_e(v) = ∫₀¹ v(x_a + s (x_b - x_a)) ⋅ (x_b - x_a) ds,
# where (a,b) are the edge's vertices in the reference edge's node order.
# The shape functions dual to these moments are the Whitney edge functions.
#
# Mesh space: each global edge is oriented from its lower to its higher global
# node ID. A cell's shape function for local edge e is
#   φ_e = s_e J⁻ᵀ φ̂_e,
# with J the Jacobian of the geometry map and s_e = ±1 the agreement between the
# local and global edge orientations. Jacobians transform as s_e J⁻ᵀ ∇̂φ̂_e J⁻¹,
# which is exact for the affine (linear simplex) cells accepted here.

"""
    nedelec_space(domain::AbstractFaceDomain,order)
    nedelec_space(domain::AbstractMeshDomain,order;dirichlet_boundary=nothing)

Return the Nédélec space of the first kind (H(curl)-conforming) of the given `order` on `domain`.

Only `order=0` (one DOF per edge, the lowest-order or Whitney element) on triangles and tetrahedra is implemented.
The mesh version requires linear simplex geometry. Its shape functions use the covariant Piola map and the edge
orientation from the lower to the higher global node ID. A `dirichlet_boundary` fixes the tangential trace: the DOFs of
the edges in the closure of the given boundary faces become Dirichlet DOFs.

Shape functions support `value` and `ForwardDiff.jacobian` (from which the curl follows). Trace (boundary or
skeleton) integration of form arguments is not available on this branch for any space.

# Level

Basic
"""
function nedelec_space end

struct CovariantPiolaMap end

function nedelec_check_order(order)
    order == 0 || throw(ArgumentError("Only lowest-order Nédélec spaces (order=0) are implemented, got order=$order"))
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
    workspace = nedelec_setup(domain)
    NedelecFaceSpace(domain,order,workspace)
end

struct NedelecFaceSpace{A,B,C} <: AbstractFaceSpace
    domain::A
    order::B
    workspace::C
end

domain(fe::NedelecFaceSpace) = fe.domain
order(fe::NedelecFaceSpace) = fe.order
num_dofs(fe::NedelecFaceSpace) = length(fe.workspace.edge_nodes)
primal_basis(fe::NedelecFaceSpace) = fe.workspace.primal_basis
dual_basis(fe::NedelecFaceSpace) = fe.workspace.dual_basis

function nedelec_primal_basis(::Val{3})
    # NED0 = P0³ ⊕ {b×x}: six functions.
    e = (SVector(1.0,0.0,0.0),SVector(0.0,1.0,0.0),SVector(0.0,0.0,1.0))
    constants = [x -> e[d] + zero(eltype(x))*x for d in 1:3]
    rotations = [x -> cross(e[d],x) for d in 1:3]
    vcat(constants,rotations)
end

function nedelec_primal_basis(::Val{2})
    # NED0 = P0² ⊕ {b(-y,x)}: three functions.
    [x -> SVector(one(eltype(x)),zero(eltype(x))),
     x -> SVector(zero(eltype(x)),one(eltype(x))),
     x -> SVector(-x[2],x[1])]
end

function nedelec_edge_quadrature()
    # Two Gauss-Legendre points on [0,1]: exact for the linear tangential traces of NED0.
    s,w = gausslegendre(2)
    (s .+ 1) ./ 2, w ./ 2
end

function nedelec_setup(domain)
    D = num_dims(domain)
    cell_mesh = mesh(domain)
    node_x = node_coordinates(cell_mesh)
    @assert length(node_x) == D+1 # the reference simplex mesh has its vertices as nodes
    edge_nodes = map(nodes->(Int(nodes[1]),Int(nodes[2])),face_nodes(cell_mesh,1))
    ss,ws = nedelec_edge_quadrature()
    dual_basis = map(edge_nodes) do (a,b)
        xa,xb = node_x[a],node_x[b]
        t = xb - xa
        v -> sum(1:length(ss)) do k
            ws[k]*dot(v(xa + ss[k]*t),t)
        end
    end
    primal_basis = nedelec_primal_basis(Val(D))
    @assert length(primal_basis) == length(dual_basis)
    reference_topology = topology(domain)
    (;edge_nodes,primal_basis,dual_basis,reference_topology)
end

function nedelec_num_face_permutations(fe::NedelecFaceSpace,d)
    topo = fe.workspace.reference_topology
    face_to_rid = face_reference_id(topo,d)
    rid_to_topo = reference_topologies(topo,d)
    map(face_to_rid) do rid
        length(vertex_permutations(rid_to_topo[rid]))
    end
end

function face_own_dofs(fe::NedelecFaceSpace,d)
    nfaces = num_faces(mesh(fe.domain),d)
    if d == 1
        [Int32[edge] for edge in 1:nfaces]
    else
        [Int32[] for _ in 1:nfaces]
    end
end

function face_dofs(fe::NedelecFaceSpace,d)
    D = num_dims(fe.domain)
    edge_nodes = fe.workspace.edge_nodes
    nedges = length(edge_nodes)
    if d == 0
        [Int32[] for _ in 1:num_faces(mesh(fe.domain),0)]
    elseif d == 1
        face_own_dofs(fe,d)
    elseif d == D
        [collect(Int32,1:nedges)]
    else
        # Closure: the edges whose two nodes belong to the face.
        map(face_nodes(mesh(fe.domain),d)) do nodes
            Int32[edge for edge in 1:nedges if all(in(nodes),edge_nodes[edge])]
        end
    end
end

function face_own_dof_permutations(fe::NedelecFaceSpace,d)
    # One DOF per edge: reversing an edge changes the sign of its shape function,
    # which the mesh space handles with the edge signs, not with a permutation.
    own_dofs = face_own_dofs(fe,d)
    nperms = nedelec_num_face_permutations(fe,d)
    map(own_dofs,nperms) do dofs,n
        [collect(Int32,1:length(dofs)) for _ in 1:n]
    end
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
    dof_signs = nedelec_dof_signs(mesh,D,ctype_to_reffe)
    workspace = nothing
    NedelecMeshSpace(domain,order,dirichlet_boundary,ctype_to_reffe,dof_signs,workspace) |> setup_space
end

struct NedelecMeshSpace{A,B,C,E,S,W} <: AbstractSpace
    domain::A
    order::B
    dirichlet_boundary::C
    reference_spaces::E
    dof_signs::S # Matrix{Int8}: [local_dof,mesh_cell] -> ±1
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
    NedelecMeshSpace(a.domain,a.order,a.dirichlet_boundary,a.reference_spaces,a.dof_signs,workspace)
end

"""
    nedelec_dof_signs(space)

Return the matrix `signs[local_dof,mesh_cell]` with the factor ±1 that relates the local edge orientation of
`mesh_cell` to the global one (from the lower to the higher global node ID).

# Level

Advanced
"""
nedelec_dof_signs(a::NedelecMeshSpace) = a.dof_signs

function nedelec_dof_signs(mesh,D,ctype_to_reffe)
    Dface_to_nodes = face_nodes(mesh,D)
    Dface_to_ctype = face_reference_id(mesh,D)
    nDfaces = length(Dface_to_nodes)
    max_ndofs = maximum(num_dofs,ctype_to_reffe)
    signs = ones(Int8,max_ndofs,nDfaces)
    for Dface in 1:nDfaces
        fe = ctype_to_reffe[Dface_to_ctype[Dface]]
        nodes = Dface_to_nodes[Dface]
        for (edge,(a,b)) in enumerate(fe.workspace.edge_nodes)
            signs[edge,Dface] = nodes[a] < nodes[b] ? Int8(1) : Int8(-1)
        end
    end
    signs
end

# Code generation for the compiler (serial, threaded and CUDA execution).

shape_function_map(::NedelecMeshSpace) = CovariantPiolaMap()

function nedelec_check_integration_domain(V,Q)
    if num_dims(GT.domain(Q)) != num_dims(GT.domain(V))
        throw(ArgumentError("Trace integration of Nédélec form arguments is not implemented"))
    end
    nothing
end

function mapped_shape_function_expr(::CovariantPiolaMap,func,mesh,V,Q,mesh_expr,V_expr,Q_expr,f_expr,q_expr,i_expr)
    nedelec_check_integration_domain(V,Q)
    sign_expr = :(nedelec_dof_signs($V_expr)[$i_expr,$f_expr])
    is_tabulated = func === value || func === ForwardDiff.jacobian
    if !is_tabulated
        func === ForwardDiff.gradient && throw(ArgumentError(
            "Nédélec shape functions are vector-valued; use ForwardDiff.jacobian instead of ForwardDiff.gradient"))
        d = num_dims(GT.domain(V))
        rid = first(face_reference_id(Q))
        xi0 = zero(eltype(coordinates(reference_quadratures(Q)[rid])))
        shape_expr = :(shape_functions(reference_spaces($V_expr)[$rid])[$i_expr])
        geometry_functions_expr = :(shape_functions(reference_spaces($mesh_expr, Val($d))[$rid]))
        mesh_data_expr = :(mesh_integration_data($mesh_expr, $Q_expr))
        node_x_expr = :(node_coordinates($mesh_data_expr))
        nodes_expr = :(face_nodes($mesh_data_expr)[$f_expr])
        geometry_space_expr = :(mesh_space($mesh_expr, Val($d)))
        geometry_data_expr = :(space_integration_data(value, $geometry_space_expr, $Q_expr))
        geometry_values_expr = :(reference_tabulator($geometry_data_expr))
        return :(evaluate_physical_covariant_shape($func, $shape_expr, $sign_expr, $geometry_functions_expr,
            $node_x_expr, $nodes_expr, $geometry_values_expr, $q_expr, $xi0))
    end
    space_data_expr = :(space_integration_data($func, $V_expr, $Q_expr))
    v_expr = :(reference_tabulator($space_data_expr)[$i_expr, $q_expr])
    J_expr = jacobian_expr(mesh, Q, mesh_expr, Q_expr, f_expr, q_expr)
    invJ_expr = :(adjugate($J_expr) * inv(det($J_expr)))
    if func === value
        :($sign_expr * (transpose($invJ_expr) * $v_expr))
    else
        :($sign_expr * (transpose($invJ_expr) * $v_expr * $invJ_expr))
    end
end

function mapped_shape_function_expr(::CovariantPiolaMap,V,mesh_expr,V_expr,f_expr,i_expr)
    rid_expr = :(face_reference_id($V_expr)[$f_expr])
    shape_expr = :(shape_functions(reference_spaces($V_expr)[$rid_expr])[$i_expr])
    sign_expr = :(nedelec_dof_signs($V_expr)[$i_expr,$f_expr])
    d = num_dims(GT.domain(V))
    geometry_rid_expr = :(face_reference_id($mesh_expr, $d)[$f_expr])
    geometry_space_expr = :(reference_spaces($mesh_expr, Val($d))[$geometry_rid_expr])
    functions_expr = :(shape_functions($geometry_space_expr))
    node_x_expr = :(node_coordinates($mesh_expr))
    nodes_expr = :(face_nodes($mesh_expr, $d)[$f_expr])
    xi0_expr = :(zero(first(node_coordinates($geometry_space_expr))))
    :(physical_covariant_shape_function($shape_expr, $sign_expr, $functions_expr, $node_x_expr, $nodes_expr, $xi0_expr))
end

function physical_covariant_shape_function(s_ref, sign, geometry_functions, node_x, nodes, xi0)
    function phi(xi)
        physical_map_value(xi, geometry_functions, node_x, nodes)
    end
    inverse_phi = inv_map(phi, xi0)
    x -> begin
        xi = inverse_phi(x)
        J = ForwardDiff.jacobian(phi, xi)
        sign * (transpose(J) \ s_ref(xi))
    end
end

function evaluate_physical_covariant_shape(func, s_ref, sign, geometry_functions, node_x,
        nodes, reference_values, point, xi0)
    s_phys = physical_covariant_shape_function(s_ref, sign, geometry_functions, node_x, nodes, xi0)
    x_q = physical_point_value(node_x, nodes, reference_values, point)
    func(s_phys, x_q)
end

# Accessors (discrete fields used as data, sampling).

function nedelec_shape_function_accessor_modifier(space::NedelecMeshSpace,transform)
    signs = nedelec_dof_signs(space)
    prototype = (v,invJ) -> transform(v,invJ)
    function face_dof_modifier(face,face_around=nothing)
        function dof_modifier(dof)
            sign = signs[dof,face]
            (v,invJ) -> sign*transform(v,invJ)
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

# The generic tabulated accessor passes no Jacobian for `value`; the covariant map needs it.
for T in (:value,:(ForwardDiff.jacobian))
    for B in (false,true)
        @eval function shape_function_accessor_physical(f::typeof($T),space::NedelecMeshSpace,
                measure::AbstractQuadrature,disable_tabulation::Val{$B})
            nedelec_shape_function_accessor_physical(f,space,measure)
        end
    end
end

function nedelec_shape_function_accessor_physical(f,space,measure)
    face_point_dof_v = shape_function_accessor_reference(f,space,measure)
    dface_to_modif = shape_function_accessor_modifier(f,space,GT.domain(measure))
    D = num_dims(GT.domain(space))
    face_point_Dphi = jacobian_accessor(measure,Val(D))
    J0 = GT.prototype(face_point_Dphi)
    prototype = GT.prototype(dface_to_modif)(GT.prototype(face_point_dof_v),inv(J0))
    function face_point_dof_s(face,face_around=nothing)
        point_dof_v = face_point_dof_v(face,face_around)
        dof_modif = dface_to_modif(face,face_around)
        point_Dphi = face_point_Dphi(face,face_around)
        function point_dof_s(point,J=nothing)
            J2 = J === nothing ? point_Dphi(point) : J
            invJ = inv(J2)
            dof_v = point_dof_v(point)
            function dof_s(dof)
                dof_modif(dof)(dof_v(dof),invJ)
            end
        end
    end
    accessor(face_point_dof_s,prototype)
end

# Interpolation: each global DOF is the tangential moment along its globally oriented edge.

function interpolate_impl!(
    f::AnalyticalField,
    u::DiscreteField,
    space::NedelecMeshSpace,
    free_or_diri::FreeOrDirichlet;location=1)
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
    s,w = gausslegendre(4)
    ss,ws = (s .+ 1) ./ 2, w ./ 2
    for Dface in faces(domain)
        fe = ctype_to_reffe[Dface_to_ctype[Dface]]
        nodes = Dface_to_nodes[Dface]
        dofs = Dface_to_dofs[Dface]
        for (edge,(a,b)) in enumerate(fe.workspace.edge_nodes)
            dof = dofs[edge]
            if dof > 0
                free_or_diri == DIRICHLET && continue
            else
                free_or_diri == FREE && continue
            end
            na,nb = minmax(nodes[a],nodes[b])
            xa,xb = node_x[na],node_x[nb]
            t = xb - xa
            v = sum(1:length(ss)) do k
                ux = fun(xa + ss[k]*t)
                ws[k]*sum(map(*,ux,t))
            end
            if dof > 0
                free_vals[dof] = v
            else
                diri_vals[-dof] = v
            end
        end
    end
    u
end

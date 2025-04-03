
function accessor(a,b)
    Accessor(a,b)
end

struct Accessor{A,B} <: AbstractAccessor
    definition::A
    prototype::B
end

function prototype(a::Accessor)
    a.prototype
end

function (a::Accessor)(face,face_around=nothing)
    a.definition(face,face_around)
end

# Untabulated
function shape_function_accessor(f,space::AbstractSpace)
    if is_reference_domain(domain(space))
        shape_function_accessor_reference(f,space)
    else
        shape_function_accessor_physical(f,space)
    end
end

# Tabulated
function shape_function_accessor(f,space::AbstractSpace,measure::AbstractQuadrature)
    if is_reference_domain(domain(space))
        shape_function_accessor_reference(f,space,measure)
    else
        shape_function_accessor_physical(f,space,measure)
    end
end

# Untabulated
function shape_function_accessor_reference(f,space::AbstractSpace)
    rid_to_dof_s = map(shape_functions,reference_spaces(space))
    face_to_rid = face_reference_id(space)
    function face_dof_g(face,face_around=nothing)
        rid = face_to_rid[face]
        dof_s = rid_to_dof_s[rid]
        function dof_g(dof)
            s = dof_s[dof]
            function g(x)
                f(s,x)
            end
        end
    end
    prototype = first(first(rid_to_dof_s))
    accessor(face_dof_g,prototype)
end

# Tabulated
function shape_function_accessor_reference(f,space::AbstractSpace,measure::AbstractQuadrature)
    domain_measure = domain(measure)
    domain_space = domain(space)
    d = num_dims(domain_measure)
    D = num_dims(domain_space)
    face_around = GT.face_around(domain_measure)
    if d == D
        shape_function_accessor_reference_interior(f,space,measure)
    elseif d+1==D && face_around !== nothing
        shape_function_accessor_reference_boundary(f,space,measure)
    else
        shape_function_accessor_reference_skeleton(f,space,measure)
    end
end

function shape_function_accessor_reference_interior(f,space::AbstractSpace,measure::AbstractQuadrature)
    mesh = GT.mesh(measure)
    dom = GT.domain(measure)
    d = num_dims(dom)
    rid_to_point_to_x = map(coordinates,reference_quadratures(measure))
    # NB the TODOs below can be solved by introducing an extra nesting level
    # TODO assumes same reference elems for integration and for interpolation
    rid_to_tab = map(rid_to_point_to_x,reference_spaces(space)) do point_to_x, refface
        tabulator(refface)(f,point_to_x)
    end
    prototype = first(first(rid_to_tab))
    face_to_rid = face_reference_id(mesh,d)
    sface_to_face = faces(dom)
    function face_point_dof_s(sface,face_around=nothing)
        face = sface_to_face[sface]
        rid = face_to_rid[face]
        tab = rid_to_tab[rid]
        function point_dof_s(point,J=nothing)
            function dof_s(dof)
                tab[point,dof]
            end
        end
    end
    accessor(face_point_dof_s,prototype)
end

function shape_function_accessor_reference_skeleton(f,space::AbstractSpace,measure::AbstractQuadrature)
    d = num_dims(domain(measure))
    D = num_dims(domain(space))
    mesh = GT.mesh(space)
    topo = topology(mesh)
    drid_to_refdface = reference_spaces(mesh,Val(d))
    Drid_to_refDface = reference_spaces(mesh,Val(D))
    Drid_to_reffe = reference_spaces(space)
    drid_to_point_to_x = map(coordinates,reference_quadratures(measure))
    dface_to_Dfaces, dface_to_ldfaces = GT.face_incidence_ext(topo,d,D)
    Dface_to_ldface_to_perm = GT.face_permutation_ids(topo,D,d)
    face_to_dface = faces(domain(measure))
    Dface_to_Drid = face_reference_id(mesh,D)
    dface_to_drid = face_reference_id(mesh,d)
    # NB the TODOs below can be solved by introducing two extra nesting levels
    # TODO this assumes the same reffes for mesh and quadrature
    drid_Drid_ldface_perm_to_tab = map(drid_to_point_to_x,drid_to_refdface) do point_to_x,refdface
        # TODO this assumes the same reffes for mesh and interpolation
        map(Drid_to_reffe,Drid_to_refDface) do reffe,refDface
            ldface_perm_varphi = reference_map(refdface,refDface)
            map(ldface_perm_varphi) do perm_varphi
                map(perm_varphi) do varphi
                    point_to_q = varphi.(point_to_x)
                    tabulator(reffe)(f,point_to_q)
                end
            end
        end
    end
    prototype = first(first(first(first(first(drid_Drid_ldface_perm_to_tab)))))
    function face_point_dof_s(face,face_around)
        dface = face_to_dface[face]
        Dfaces = dface_to_Dfaces[dface]
        ldfaces = dface_to_ldfaces[dface]
        Dface = Dfaces[face_around]
        ldface = ldfaces[face_around]
        drid = dface_to_drid[dface]
        Drid = Dface_to_Drid[Dface]
        perm = Dface_to_ldface_to_perm[Dface][ldface]
        tab = drid_Drid_ldface_perm_to_tab[drid][Drid][ldface][perm]
        function point_dof_s(point)
            function dof_s(dof)
                tab[point,dof]
            end
        end
    end
    accessor(face_point_dof_s,prototype)
end

function shape_function_accessor_reference_boundary(f,space::AbstractSpace,measure::AbstractQuadrature)
    face_around = GT.face_around(domain(measure))
    face_point_dof_s = shape_function_accessor_reference_skeleton(f,space,measure)
    prototype = GT.prototype(face_point_dof_s)
    function face_point_dof_b(face,dummy=nothing)
        face_point_dof_s(face,face_around)
    end
    accessor(face_point_dof_b,prototype)
end

function nodes_accessor(mesh::AbstractMesh,vD,domain::AbstractDomain)
    D = val_parameter(vD)
    d = num_dims(domain)
    face_around = GT.face_around(domain)
    if d == D
        nodes_accessor_interior(mesh,vD,domain)
    elseif d+1==D && face_around !== nothing
        nodes_accessor_boundary(mesh,vD,domain)
    else
        nodes_accessor_skeleton(mesh,vD,domain)
    end
end

function nodes_accessor_interior(mesh::AbstractMesh,vD,domain::AbstractDomain)
    D = val_parameter(vD)
    Dface_nodes = face_nodes(mesh,D)
    face_to_Dface = faces(domain)
    function face_to_nodes(face,face_around=nothing)
        Dface = face_to_Dface[face]
        Dface_nodes[Dface]
    end
    prototype = Int32[1]
    accessor(face_to_nodes,prototype)
end

function nodes_accessor_skeleton(mesh::AbstractMesh,vD,domain::AbstractDomain)
    D = val_parameter(vD)
    d = num_dims(domain)
    Dface_nodes = face_nodes(mesh,D)
    face_to_dface = faces(domain)
    topo = topology(mesh)
    dface_to_Dfaces = face_incidence(topo,d,D)
    function face_to_nodes(face,face_around)
        dface = face_to_dface[face]
        Dfaces = dface_to_Dfaces[dface]
        Dface = Dfaces[face_around]
        Dface_nodes[Dface]
    end
    prototype = Int32[1]
    accessor(face_to_nodes,prototype)
end

function nodes_accessor_boundary(mesh::AbstractMesh,vD,domain::AbstractDomain)
    face_to_nodes = nodes_accessor_skeleton(mesh,vD,domain)
    face_around = GT.face_around(domain)
    prototype = GT.prototype(face_to_nodes)
    function face_to_nodes_2(face,dummy=nothing)
        face_to_nodes(face,face_around)
    end
    accessor(face_to_nodes_2,prototype)
end

# Untabulated
function physical_map_accessor(f,mesh::AbstractMesh,vD)
    D = val_parameter(vD)
    space = mesh_space(mesh,vD)
    face_lnode_s = shape_function_accessor_reference(value,space)
    node_to_x = node_coordinates(mesh)
    Dface_nodes = face_nodes(mesh,D)
    z = zero(eltype(node_to_x))
    prototype = q->f(y->GT.prototype(face_lnode_s)(y)*z,q)
    function face_phi(Dface,face_around=nothing)
        lnode_to_s = face_lnode_s(Dface,face_around)
        lnode_to_node = Dface_nodes[Dface]
        function phi(q)
            nlnodes = length(lnode_to_node)
            sum(1:nlnodes) do lnode
                node = lnode_to_node[lnode]
                x = node_to_x[node]
                s = lnode_to_s(lnode)
                s(q)*x
            end
        end
        q->f(phi,q)
    end
    accessor(face_phi,prototype)
end

# Tabulated

function physical_map_accessor(f,measure::AbstractQuadrature,vD)
    error("Case not implemented and not needed unless you want higher order derivatives of the physical map.")
end

function physical_map_accessor(f::typeof(value),measure::AbstractQuadrature,vD)
    mesh = GT.mesh(domain(measure))
    space = mesh_space(mesh,vD)
    face_point_lnode_s = shape_function_accessor_reference(f,space,measure)
    face_to_nodes = nodes_accessor(mesh,vD,domain(measure))
    node_to_x = node_coordinates(mesh)
    z = zero(eltype(node_to_x))
    prototype = z*GT.prototype(face_point_lnode_s)
    function face_point_phi(face,face_around=nothing)
        point_lnode_s = face_point_lnode_s(face,face_around)
        lnode_to_node = face_to_nodes(face,face_around)
        function point_s(point)
            lnode_s = point_lnode_s(point)
            nlnodes = length(lnode_to_node)
            sum(1:nlnodes) do lnode
                node = lnode_to_node[lnode]
                x = node_to_x[node]
                s = lnode_s(lnode)
                x*s
            end
        end
    end
    accessor(face_point_phi,prototype)
end

function physical_map_accessor(f::typeof(ForwardDiff.jacobian),measure::AbstractQuadrature,vD)
    mesh = GT.mesh(domain(measure))
    space = mesh_space(mesh,vD)
    face_point_lnode_s = shape_function_accessor_reference(ForwardDiff.gradient,space,measure)
    face_to_nodes = nodes_accessor(mesh,vD,domain(measure))
    node_to_x = node_coordinates(mesh)
    z = zero(eltype(node_to_x))
    prototype = outer(z,GT.prototype(face_point_lnode_s))
    function face_point_phi(face,face_around=nothing)
        point_lnode_s = face_point_lnode_s(face,face_around)
        lnode_to_node = face_to_nodes(face,face_around)
        function point_s(point)
            lnode_s = point_lnode_s(point)
            nlnodes = length(lnode_to_node)
            sum(1:nlnodes) do lnode
                node = lnode_to_node[lnode]
                x = node_to_x[node]
                s = lnode_s(lnode)
                outer(x,s)
            end
        end
    end
    accessor(face_point_phi,prototype)
end

function coordinate_accessor(measure::AbstractQuadrature)
    vD=Val(num_dims(domain(measure)))
    physical_map_accessor(value,measure,vD)
end

function jacobian_accessor(measure::AbstractQuadrature,vD=Val(num_dims(domain(measure))))
    physical_map_accessor(ForwardDiff.jacobian,measure,vD)
end

function num_points_accessor(measure::AbstractQuadrature)
    rid_to_point_to_w = map(weights,reference_quadratures(measure))
    rid_to_n = map(length,rid_to_point_to_w)
    dface_to_rid = face_reference_id(measure)
    face_to_dface = faces(domain(measure))
    prototype = 1
    function face_npoints(face,face_around=nothing)
        dface = face_to_dface[face]
        rid = dface_to_rid[dface]
        rid_to_n[rid]
    end
    accessor(face_npoints,prototype)
end

function weight_accessor(measure::AbstractQuadrature)
    if is_reference_domain(domain(measure))
        weight_accessor_reference(measure)
    else
        weight_accessor_physical(measure)
    end
end

function weight_accessor_reference(measure::AbstractQuadrature)
    face_point_J = jacobian_accessor(measure)
    rid_to_point_to_w = map(weights,reference_quadratures(measure))
    dface_to_rid = face_reference_id(measure)
    face_to_dface = faces(domain(measure))
    prototype = first(first(rid_to_point_to_w))
    function face_point_v(face,face_around=nothing)
        dface = face_to_dface[face]
        rid = dface_to_rid[dface]
        point_to_w = rid_to_point_to_w[rid]
        function point_v(point)
            point_to_w[point]
        end
    end
    accessor(face_point_v,prototype)
end

function weight_accessor_physical(measure::AbstractQuadrature)
    face_point_v = weight_accessor_reference(measure)
    face_point_J = jacobian_accessor(measure)
    prototype = change_of_measure(GT.prototype(face_point_J))*GT.prototype(face_point_v)
    function face_point_w(face,face_around=nothing)
        point_v = face_point_v(face)
        point_J = face_point_J(face)
        function point_w(point,J = point_J(point))
            v = point_v(point)
            change_of_measure(J)*v
        end
        return point_w
    end
    accessor(face_point_w,prototype)
end

# push reference shape functions to physical ones

function shape_function_accessor_modifier(f::typeof(value),space::AbstractSpace)
    function modifier(v,J)
        v
    end
    function face_dof_modifier(face,face_around=nothing)
        function dof_modifier(dof)
            return modifier
        end
    end
    accessor(face_dof_modifier,modifier)
end

function shape_function_accessor_modifier(f::typeof(ForwardDiff.gradient),space::AbstractSpace)
    function modifier(v,J)
        transpose(J)\v
    end
    function face_dof_modifier(face,face_around=nothing)
        function dof_modifier(dof)
            return modifier
        end
    end
    accessor(face_dof_modifier,modifier)
end

function shape_function_accessor_modifier(f::typeof(ForwardDiff.jacobian),space::AbstractSpace)
    function modifier(v,J)
        v/J
    end
    function face_dof_modifier(face,face_around=nothing)
        function dof_modifier(dof)
            return modifier
        end
    end
    accessor(face_dof_modifier,modifier)
end


function shape_function_accessor_modifier(f,space::AbstractSpace,domain::AbstractDomain)
    D = num_dims(space)
    d = num_dims(domain)
    face_around = GT.face_around(domain)
    if d == D
        shape_function_accessor_modifier_interior(f,space,domain)
    elseif d+1==D && face_around !== nothing
        shape_function_accessor_modifier_boundary(f,space,domain)
    else
        shape_function_accessor_modifier_skeleton(f,space,domain)
    end
end

function shape_function_accessor_modifier_interior(f,space::AbstractSpace,domain::AbstractDomain)
    Dface_to_modif = shape_function_accessor_modifier(f,space)
    face_to_Dface = faces(domain)
    function face_to_modif(face,face_around=nothing)
        Dface = face_to_Dface[face]
        Dface_to_modif(Dface)
    end
    accessor(face_to_modif,GT.prototype(Dface_to_modif))
end

function shape_function_accessor_modifier_skeleton(f,space::AbstractSpace,domain::AbstractDomain)
    D = num_dims(GT.domain(space))
    d = num_dims(domain)
    face_to_dface = faces(domain)
    topo = topology(GT.mesh(space))
    dface_to_Dfaces = face_incidence(topo,d,D)
    Dface_to_modif = shape_function_accessor_modifier(f,space)
    function face_to_modif(face,face_around)
        dface = face_to_dface[face]
        Dfaces = dface_to_Dfaces[dface]
        Dface = Dfaces[face_around]
        Dface_to_modif(Dface)
    end
    accessor(face_to_modif,GT.prototype(Dface_to_modif))
end

function shape_function_accessor_modifier_boundary(f,space::AbstractSpace,domain::AbstractDomain)
    face_modif = shape_function_accessor_modifier_skeleton(f,space,domain)
    face_around = GT.face_around(domain)
    function face_to_modif(face,dummy=nothing)
        face_modif(face,face_around)
    end
    accessor(face_to_modif,GT.prototype(face_modif))
end

# Untabulated
function shape_function_accessor_physical(f,space::AbstractSpace)
    domain = GT.domain(space)
    D = num_dims(domain)
    mesh = GT.mesh(domain)
    face_dof_s_ref = shape_function_accessor_reference(value,space)
    face_dof_modif = shape_function_accessor_modifier(value,space)
    face_phi = physical_map_accessor(value,mesh,Val(D))
    face_Dphi = physical_map_accessor(ForwardDiff.jacobian,mesh,Val(D))
    x0 = zero(SVector{D,Float64})
    prototype = q->f(x->GT.prototype(face_dof_modif)(GT.prototype(face_dof_s_ref)(x),GT.prototype(face_Dphi)(x)),q)
    function face_dof_s_phys(face,face_around=nothing)
        phi = face_phi(face)
        Dphi = face_Dphi(face)
        invphi = inv_map(phi,x0)
        dof_s_ref = face_dof_s_ref(face)
        dof_modif = face_dof_modif(face)
        function dof_s_phys(dof)
            modif = dof_modif(dof)
            s_ref = dof_s_ref(dof)
            function s_phys(x)
                v = s_ref(x)
                J = Dphi(x)
                modif(v,J)
            end
            x->f(s_phys,x)
        end
    end
    accessor(face_dof_s_phys,prototype)
end

# This gets tabulated only for particular instances of f.
function shape_function_accessor_physical(f,space::AbstractSpace,measure::AbstractQuadrature)
    if num_dims(domain(space)) != num_dims(domain(measure))
        error("case not implemented")
    end
    Dface_dof_x_s = shape_function_accessor_physical(f,space)
    face_point_x = coordinate_accessor(measure)
    face_dface = faces(domain(measure))
    function face_point_dof_s(face,face_around=nothing)
        point_x = face_point_x(face)
        Dface = face_to_dface[face]
        dof_x_s = Dface_dof_x_s(Dface)
        function point_dof_s(point,J=nothing)
            x = point_x(point)
            function dof_s(dof)
                dof_x_s(dof)(x)
            end
        end
    end
    prototype = GT.prototype(Dface_dof_x_s)(GT.prototype(face_point_x))
    accessor(face_point_dof_s,prototype)
end

# Tabulated
for T in (:value,:(ForwardDiff.gradient),:(ForwardDiff.jacobian))
    @eval begin
        function shape_function_accessor_physical(f::typeof($T),space::AbstractSpace,measure::AbstractQuadrature)
            face_point_dof_v = shape_function_accessor_reference(f,space,measure)
            face_ndofs = num_dofs_accessor(space,GT.domain(measure))
            dface_to_modif = shape_function_accessor_modifier(f,space,domain(measure))
            D = num_dims(domain(space))
            face_point_Dphi = jacobian_accessor(measure,Val(D))
            prototype = GT.prototype(dface_to_modif)(GT.prototype(face_point_dof_v),GT.prototype(face_point_Dphi))
            P = typeof(prototype)
            max_n_faces_around = 2
            face_around_dof_s = fill(zeros(P,max_num_reference_dofs(space)),max_n_faces_around)
            function face_point_dof_s(face,face_around=nothing)
                point_dof_v = face_point_dof_v(face,face_around)
                dof_modif = dface_to_modif(face,face_around)
                ndofs = face_ndofs(face,face_around)
                point_Dphi = face_point_Dphi(face,face_around)
                dof_s = if face_around === nothing
                    face_around_dof_s[1]
                else
                    face_around_dof_s[face_around]
                end
                function point_J_dof_s(point,J)
                    dof_v = point_dof_v(point)
                    for dof in 1:ndofs # TODO: do the jacobian only once for simplices
                        v = dof_v(dof)
                        modif = dof_modif(dof)
                        dof_s[dof] = modif(v,J)
                    end
                    function dof_f(dof)
                        dof_s[dof]
                    end
                end
                function point_dof_s(point,J = nothing)
                    J2 = J === nothing ? point_Dphi(point) : J
                    point_J_dof_s(point,J2)
                end
                return point_dof_s
            end
            accessor(face_point_dof_s,prototype)
        end # function
    end # @eval
end # for

function dofs_accessor(space::AbstractSpace,domain::AbstractDomain)
    D = num_dims(space)
    d = num_dims(domain)
    face_around = GT.face_around(domain)
    if d == D
        dofs_accessor_interior(space,domain)
    elseif d+1==D && face_around !== nothing
        dofs_accessor_boundary(space,domain)
    else
        dofs_accessor_skeleton(space,domain)
    end
end

function dofs_accessor_interior(space::AbstractSpace,dom::AbstractDomain)
    face_to_Dface = faces(dom)
    Dface_to_dofs_ = face_dofs(space)
    prototype = Int32[1]
    function face_to_dofs(face,face_around=nothing)
        Dface = face_to_Dface[face]
        dofs = Dface_to_dofs_[Dface]
        dofs
    end
    prototype = Int32[1]
    accessor(face_to_dofs,prototype)
end

function dofs_accessor_skeleton(space::AbstractSpace,domain::AbstractDomain)
    D = num_dims(GT.domain(space))
    d = num_dims(domain)
    Dface_dofs = face_dofs(space)
    face_to_dface = faces(domain)
    topo = topology(GT.mesh(space))
    dface_to_Dfaces = face_incidence(topo,d,D)
    function face_to_dofs(face,face_around)
        dface = face_to_dface[face]
        Dfaces = dface_to_Dfaces[dface]
        Dface = Dfaces[face_around]
        Dface_dofs[Dface]
    end
    prototype = Int32[1]
    accessor(face_to_dofs,prototype)
end

function dofs_accessor_boundary(space::AbstractSpace,domain::AbstractDomain)
    face_to_dofs = dofs_accessor_skeleton(space,domain)
    face_around = GT.face_around(domain)
    prototype = GT.prototype(face_to_dofs)
    function face_to_dofs_2(face,dummy=nothing)
        face_to_dofs(face,face_around)
    end
    accessor(face_to_dofs_2,prototype)
end

function num_dofs_accessor(space::AbstractSpace,domain::AbstractDomain)
    D = num_dims(space)
    d = num_dims(domain)
    face_around = GT.face_around(domain)
    if d == D
        num_dofs_accessor_interior(space,domain)
    elseif d+1==D && face_around !== nothing
        num_dofs_accessor_boundary(space,domain)
    else
        num_dofs_accessor_skeleton(space,domain)
    end
end

function num_dofs_accessor_interior(space::AbstractSpace,domain::AbstractDomain)
    rid_to_n = map(num_dofs,reference_spaces(space))
    Dface_to_rid = face_reference_id(space)
    face_to_Dface = faces(domain)
    function face_to_n(face,face_around=nothing)
        Dface = face_to_Dface[face]
        rid = Dface_to_rid[Dface]
        rid_to_n[rid]
    end
    accessor(face_to_n,1)
end

function num_dofs_accessor_skeleton(space::AbstractSpace,domain::AbstractDomain)
    mesh = GT.mesh(space)
    rid_to_n = map(num_dofs,reference_spaces(space))
    Dface_to_rid = face_reference_id(space)
    face_to_dface = faces(domain)
    topo = topology(mesh)
    d = num_dims(domain)
    D = num_dims(GT.domain(space))
    dface_to_Dfaces = face_incidence(topo,d,D)
    function face_to_n(face,face_around)
        dface = face_to_dface[face]
        Dfaces = dface_to_Dfaces[dface]
        Dface = Dfaces[face_around]
        rid = Dface_to_rid[Dface]
        rid_to_n[rid]
    end
    accessor(face_to_n,1)
end

function num_dofs_accessor_boundary(space::AbstractSpace,domain::AbstractDomain)
    a = num_dofs_accessor_skeleton(space,domain)
    face_around = GT.face_around(domain)
    function face_to_n(face,dummy=nothing)
        a(face,face_around)
    end
    accessor(face_to_n,GT.prototype(a))
end

struct DiscreteFieldAccessor{A,B} <: AbstractAccessor
    update::A
    accessor::B
end

prototype(f::DiscreteFieldAccessor) = prototype(f.accessor)

function (f::DiscreteFieldAccessor)(face,face_around=nothing)
    f.accessor(face,face_around)
end

function update(f::DiscreteFieldAccessor;discrete_field)
    uh = discrete_field
    accessor = f.update(uh)
    DiscreteFieldAccessor(f.update,accessor)
end

function discrete_field_accessor(f,uh::DiscreteField,measure::AbstractQuadrature)
    face_to_dofs = dofs_accessor(GT.space(uh),GT.domain(measure))
    face_to_point_to_ldof_to_s = shape_function_accessor(f,GT.space(uh),measure)
    function field_to_accessor(uh)
        space = GT.space(uh)
        free_values = GT.free_values(uh)
        dirichlet_values = GT.dirichlet_values(uh)
        prototype = zero(eltype(free_values))*GT.prototype(face_to_point_to_ldof_to_s)
        function face_point_u(face,face_around=nothing)
            ldof_to_dof = face_to_dofs(face,face_around)
            point_to_ldof_to_s = face_to_point_to_ldof_to_s(face,face_around)
            function point_u(point,J=nothing)
                ldof_to_s = point_to_ldof_to_s(point,J)
                nldofs = length(ldof_to_dof)
                sum(1:nldofs) do ldof
                    dof = ldof_to_dof[ldof]
                    s = ldof_to_s(ldof)
                    if dof > 0
                        v = free_values[dof]
                    else
                        v = dirichlet_values[-dof]
                    end
                    v*s
                end
            end
        end
        GT.accessor(face_point_u,prototype)
    end
    accessor = field_to_accessor(uh)
    DiscreteFieldAccessor(field_to_accessor,accessor)
end

## TODO not needed anymore
#struct DirichletAccessor{A,B} <: AbstractAccessor
#    update::A
#    accessor::B
#end
#
#prototype(f::DirichletAccessor) = prototype(f.accessor)
#
#function (f::DirichletAccessor)(face)
#    f.accessor(face)
#end
#
#function update(f::DirichletAccessor;discrete_field)
#    uh = discrete_field
#    accessor = f.update(uh)
#    DirichletAccessor(f.update,accessor)
#end
#
#function dirichlet_accessor(uh::DiscreteField,domain::AbstractDomain)
#    face_to_dofs = dofs_accessor(GT.space(uh),domain)
#    function field_to_accessor(uh)
#        space = GT.space(uh)
#        dirichlet_values = GT.dirichlet_values(uh)
#        prototype = nothing
#        function face_dirichlet!(face,face_around=nothing)
#            ldof_to_dof = face_to_dofs(face,face_around)
#            nldofs = length(ldof_to_dof)
#            function dirichlet!(A,b)
#                m = size(A,1)
#                z = zero(eltype(b))
#                for i in 1:m
#                    bi = z
#                    for j in 1:nldofs
#                        dof = ldof_to_dof[j]
#                        if dof < 0
#                            uj = dirichlet_values[-dof]
#                            bi += A[i,j]*uj
#                        end
#                    end
#                    b[i] -= bi
#                end
#                nothing
#            end
#        end
#        GT.accessor(face_dirichlet!,prototype)
#    end
#    accessor = field_to_accessor(uh)
#    DirichletAccessor(field_to_accessor,accessor)
#end

function unit_normal_accessor(measure::AbstractQuadrature)
    domain = GT.domain(measure)
    mesh = GT.mesh(domain)
    D = num_dims(mesh)
    d = num_dims(domain)
    if D==d && is_physical_domain(domain)
        unit_normal_accessor_physical_interior(measure)
    elseif D==d+1 && is_physical_domain(domain)
        unit_normal_accessor_physical(measure)
    elseif D==d+1 && is_reference_domain(domain)
        unit_normal_accessor_reference(measure)
    else
        error()
    end
end

function unit_normal_accessor_physical_interior(measure::AbstractQuadrature)
    error("not implemented")
end

function unit_normal_accessor_physical(measure::AbstractQuadrature)
    domain = GT.domain(measure)
    d = num_dims(domain)
    D = d+1
    face_n_ref = unit_normal_accessor_reference(measure)
    face_point_J = jacobian_accessor(measure,Val(D))
    function face_point_n(face,face_around=nothing)
        point_J = face_point_J(face,face_around)
        n_ref = face_n_ref(face,face_around)
        function n_phys(point, J=nothing)
            J2 = J === nothing ? point_J(point) : J
            map_unit_normal(J2,n_ref)
        end
        return n_phys
    end
    accessor(face_point_n,GT.prototype(face_n_ref))
end

function unit_normal_accessor_reference(measure::AbstractQuadrature)
    D = num_dims(mesh(measure))
    d = num_dims(domain(measure))
    @assert D == d+1
    if GT.face_around(domain(measure)) === nothing
        unit_normal_accessor_reference_skeleton(measure)
    else
        unit_normal_accessor_reference_boundary(measure)
    end
end

function unit_normal_accessor_reference_skeleton(measure::AbstractQuadrature)
    domain = GT.domain(measure)
    d = num_dims(domain)
    D = d+1
    mesh = GT.mesh(domain)
    topo = topology(mesh)
    dface_to_Dfaces, dface_to_ldfaces = GT.face_incidence_ext(topo,d,D)
    Dface_to_Drid = face_reference_id(mesh,D)
    face_to_dface = faces(domain)
    Drid_to_ldface_to_n = map(GT.reference_spaces(mesh,Val(D))) do refface
        boundary = refface |> GT.domain |> GT.mesh
        boundary |> GT.outward_normals # TODO also rename?
    end
    prototype = first(first(Drid_to_ldface_to_n))
    function face_n(face,face_around)
        dface = face_to_dface[face]
        Dfaces = dface_to_Dfaces[dface]
        ldfaces = dface_to_ldfaces[dface]
        Dface = Dfaces[face_around]
        ldface = ldfaces[face_around]
        Drid = Dface_to_Drid[Dface]
        n = Drid_to_ldface_to_n[Drid][ldface]
    end
    accessor(face_n,prototype)
end

function unit_normal_accessor_reference_boundary(measure::AbstractQuadrature)
    face_n = unit_normal_accessor_reference_skeleton(measure)
    face_around = GT.face_around(domain(measure))
    function face_n_2(face,dummy=nothing)
        face_n(face,face_around)
    end
    accessor(face_n_2,GT.prototype(face_n))
end

function face_diameter(domain::AbstractDomain)
    d = num_dims(domain)
    dinv = 1/d
    measure = GT.measure(domain,1)
    face_point_w = weight_accessor(measure)
    face_npoints = num_points_accessor(measure)
    z = zero(prototype(face_point_w))
    nfaces = num_faces(domain)
    diams = fill(z,nfaces)
    for face in 1:nfaces
        point_w = face_point_w(face)
        npoints = face_npoints(face)
        s = z
        for point in 1:npoints
            w = point_w(point)
            s += w
        end
        diams[face] =  s^dinv
    end
    diams
end

#function form_argument_accessor(f,space::AbstractSpace)
#    shape_function_accessor(f,space)
#end
#
# TODO: type instability when using accessor for customized functions
function form_argument_accessor(f,space::AbstractSpace,measure::AbstractQuadrature,field=1)
    face_point_dof_s = shape_function_accessor(f,space,measure)
    prototype = GT.prototype(face_point_dof_s)
    the_field = field
    function face_point_dof_a(face,face_around=nothing)
        the_face_around = face_around
        point_dof_s = face_point_dof_s(face,face_around)
        function point_dof_a(point, J=nothing)
            dof_s = point_dof_s(point, J)
            function dof_a(dof,field=1,face_around=nothing)
                s = dof_s(dof)
                mask = face_around == the_face_around && field == the_field
                if mask
                    s
                else
                    zero(s)
                end
            end
        end
    end
    accessor(face_point_dof_a,prototype)
end

function num_faces_around_accesor(domain_space,domain)
    d = num_dims(domain)
    D = num_dims(domain_space)
    face_around = GT.face_around(domain)
    if d == D
        num_faces_around_accesor_interior(domain_space,space)
    elseif d+1==D && face_around !== nothing
        num_faces_around_accesor_interior(domain_space,space)
    else
        num_faces_around_accesor_skeleton(domain_space,space)
    end
end

function num_faces_around_accesor_interior(space_domain,domain)
    function n_faces_around(face,face_around=nothing)
        1
    end
    accessor(n_faces_around,1)
end

function num_faces_around_accesor_skeleton(space_domain,domain)
    function n_faces_around(face,face_around=nothing)
        2
    end
    accessor(n_faces_around,1)
end

# remove it as we use quadrature directly
# function num_points_accessor(measure::Measure)
#     num_points_accessor(quadrature(measure))
# end

# function coordinate_accessor(measure::Measure)
#     coordinate_accessor(quadrature(measure))
# end

# function jacobian_accessor(measure::Measure,args...)
#     jacobian_accessor(quadrature(measure),args...)
# end

# function weight_accessor(measure::Measure)
#     weight_accessor(quadrature(measure))
# end

# function discrete_field_accessor(f,uh::DiscreteField,measure::Measure)
#     discrete_field_accessor(f,uh,quadrature(measure))
# end

# function shape_function_accessor(f,space::AbstractSpace,measure::Measure)
#     shape_function_accessor(f,space,quadrature(measure))
# end

# function form_argument_accessor(f,space::AbstractSpace,measure::Measure)
#     form_argument_accessor(f,space,quadrature(measure))
# end

# function physical_map_accessor(f,measure::Measure,vD)
#     physical_map_accessor(f,quadrature(measure),vD)
# end

# function unit_normal_accessor(measure::Measure)
#     unit_normal_accessor(quadrature(measure))
# end


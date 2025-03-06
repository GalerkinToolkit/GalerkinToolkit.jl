
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

function (a::Accessor)(face)
    a.definition(face)
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
        Drid = Dface_to_Drid[dface]
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
    face_point_dof_s = shape_function_accessor_physical_skeleton(f,space,measure)
    prototype = GT.prototype(face_point_dof_s)
    function face_point_dof_b(face)
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
        Dfaces = face_to_Dfaces[face]
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
    space = mesh_space(mesh,vD)
    face_lnode_s = shape_function_accessor_reference(value,space)
    node_to_x = node_coordinates(mesh)
    Dface_nodes = face_nodes(mesh,D)
    z = zero(eltype(node_to_x))
    prototype = q->f(y->GT.prototype(face_lnode_s)(y)*z,q)
    function face_phi(Dface)
        lnode_to_s = face_lnode_s(face)
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
                s = lnode_to_s(lnode)
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
                s = lnode_to_s(lnode)
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
        rid = Dface_to_rid[dface]
        rid_to__n[rid]
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
        rid = Dface_to_rid[dface]
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
    prototype = change_of_measure(prototype(face_point_J))*prototype(face_point_v)
    function face_point_w(face,face_around=nothing)
        point_v = face_point_v(face)
        point_J = face_point_J(face)
        function point_w(point)
            v = point_v(point)
            J = point_J(point)
            change_of_measure(J)*v
        end
        function point_w(point,J)
            v = point_v(point)
            change_of_measure(J)*v
        end
        return point_w
    end
    accessor(face_point_w,prototype)
end

# push reference shape functions to physical ones

function shape_function_accessor_modifier(f::ForwardDiff.value,space::AbstractSpace)
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

function shape_function_accessor_modifier(f::ForwardDiff.gradient,space::AbstractSpace)
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

function shape_function_accessor_modifier(f::ForwardDiff.jacobian,space::AbstractSpace)
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

# Untabulated
function shape_function_accessor_physical(f,space::AbstractSpace)
    domain = GT.domain(space)
    D = num_dims(domain)
    mesh = GT.mesh(domain)
    face_dof_s_ref = shape_function_accessor_reference(value,space)
    face_dof_modif = shape_function_accessor_modifier(value,space)
    face_phi = physical_map_accessor(mesh,Val(D))
    face_Dphi = physical_map_accessor(ForwardDiff.jacobian,mesh,Val(D))
    x0 = zero(SVector{D,Float64})
    prototype = q->f(x->GT.prototype(face_dof_modif)(GT.prototype(face_dof_s_ref)(x),GT.prototype(face_Dphi)(x)),q)
    function face_dof_s_phys(face)
        phi = face_phi(face)
        Dphi = face_Dphi(face)
        invphi = inv_map(phi,x0)
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
    Dface_dof_x_s = shape_function_accessor_physical(f,space)
    face_point_x = coordinate_accessor(measure)
    face_Dface = faces(domain(measure))
    function face_point_dof_s(face)
        point_x = face_point_x(face)
        Dface = face_to_Dface[face]
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
for T in (value,ForwardDiff.gradient,ForwardDiff.jacobian)
    @eval begin
        function shape_function_accessor_physical(f::$T,space::AbstractSpace,measure::AbstractQuadrature)
            face_point_dof_v = shape_function_accessor_reference(f,space,measure)
            face_ndofs = num_dofs_accessor(space,measure)
            Dface_dof_modif = shape_function_accessor_modifier(f,space)
            D = num_dims(domain(space))
            face_point_Dphi = jacobian_accessor(measure,Val(D))
            face_to_Dface = faces(domain(measure))
            Dface_to_rid = face_reference_id(space)
            D = num_dims(domain(space))
            prototype = GT.prototype(Dface_dof_modif)(GT.prototype(face_point_dof_v),GT.prototype(face_point_Dphi))
            rid_to_dof_s = map(reference_spaces(space)) do reffe
                P = typeof(prototype)
                zeros(P,num_dofs(reffe))
            end
            function face_point_dof_s(face,face_around=nothing)
                Dface = face_to_Dface[face]
                point_dof_v = face_point_dof_v(face,face_around)
                dof_modif = Dface_dof_modif(Dface,face_around)
                rid = Dface_to_rid[Dface]
                dof_s = rid_to_dof_s[rid]
                ndofs = face_ndofs(face,face_around)
                function point_dof_s(point,J=nothing)
                    dof_v = point_dof_v(point)
                    ndofs = length(dof_s)
                    for dof in 1:ndofs
                        v = dof_v(dof)
                        modif = dof_modif(dof)
                        dof_s[dof] = modif(v,J)
                    end
                    function dof_f(dof)
                        dof_s[dof]
                    end
                end
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
        face = face_to_Dface[face]
        dofs = Dface_to_dofs_[Dface]
        dofs
    end
end

function dofs_accessor_skeleton(space::AbstractSpace,domain::AbstractDomain)
    D = num_dims(GT.domain(space))
    d = num_dims(domain)
    Dface_dofs = face_dofs(space)
    face_to_dface = faces(domain)
    topo = topology(mesh)
    dface_to_Dfaces = face_incidence(topo,d,D)
    function face_to_dofs(face,face_around)
        Dfaces = face_to_Dfaces[face]
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

function num_dofs_accessor_interior(space::AbstractSpace,dom::AbstractDomain)
    rid_to_n = map(num_dofs,reference_spaces(space))
    Dface_to_rid = face_reference_id(space)
    face_to_Dface = faces(dom)
    function face_to_n(face,face_around=nothing)
        Dface = face_to_Dface[face]
        rid = Dface_to_rid[Dface]
        rid_to_n[rid]
    end
    accessor(face_to_n,1)
end

function num_dofs_accessor_skeleton(space::AbstractSpace,dom::AbstractDomain)
    mesh = GT.mesh(space)
    rid_to_n = map(num_dofs,reference_spaces(space))
    Dface_to_rid = face_reference_id(space)
    face_to_dface = faces(dom)
    topo = topology(mesh)
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

function (f::DiscreteFieldAccessor)(face)
    f.accessor(face)
end

function update(f::DiscreteFieldAccessor;discrete_field)
    uh = discrete_field
    accessor = f.update(uh)
    DiscreteFieldAccessor(f.update,accessor)
end

function discrete_field_accessor(f,uh::DiscreteField,measure::AbstractQuadrature)
    face_to_dofs = dofs_accessor(GT.space(uh),measure)
    face_to_point_to_ldof_to_s = shape_function_accessor(f,GT.space(uh),measure)
    function field_to_accessor(uh)
        space = GT.space(uh)
        free_values = GT.free_values(uh)
        dirichlet_values = GT.dirichlet_values(uh)
        prototype = zero(eltype(free_values))*GT.prototype(face_to_point_to_ldof_to_s)
        function face_point_u(face,face_around=nothing)
            ldof_to_dof = face_to_dofs(face,face_around)
            point_to_ldof_to_s = face_to_point_to_ldof_to_s(face,face_around)
            function point_u(point)
                ldof_to_s point_to_ldof_to_s(point)
                nldofs = length(ldof_to_dof)
                sum(1:nldofs) do ldof
                    s = ldof_to_s(ldof)
                    if dof > 0
                        v = f.free_vals[dof]
                    else
                        v = f.diri_vals[-dof]
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


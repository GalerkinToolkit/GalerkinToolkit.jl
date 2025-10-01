
function each_face(a::NewAbstractAccessor)
    EachFace(a)
end

function each_face(mesh::AbstractMesh,args...;kwargs...)
    mesh_acc = mesh_accessor(mesh,args...;kwargs...)
    each_face(mesh_acc)
end

function each_face(quadrature::AbstractQuadrature,args...;kwargs...)
    acc = quadrature_accessor(quadrature,args...;kwargs...)
    each_face(acc)
end

function each_face(space::AbstractSpace,args...;kwargs...)
    acc = space_accessor(space,args...;kwargs...)
    each_face(acc)
end

function each_face(uh::AbstractField,args...;kwargs...)
    acc = field_accessor(uh,args...;kwargs...)
    each_face(acc)
end

struct EachFace{A} <: AbstractType
    accessor::A
end

function tabulate(f,iter::EachFace)
    accessor = tabulate(f,iter.accessor)
    EachFace(accessor)
end

function compute(f,iter::EachFace)
    accessor = compute(f,iter.accessor)
    EachFace(accessor)
end

Base.length(iter::EachFace) = num_faces(iter.accessor)
Base.isdone(iter::EachFace,face) = face > length(iter)
Base.getindex(iter::EachFace,face) = at_face(iter.accessor,face)
Base.keys(iter::EachFace) = LinearIndices((length(iter),))

function Base.iterate(iter::EachFace,face=1)
    if Base.isdone(iter,face)
        nothing
    else
        accessor = iter[face]
        (accessor,face+1)
    end
end

function each_face_around(a::NewAbstractAccessor)
    EachFaceAround(a)
end

struct EachFaceAround{A} <: AbstractType
    accessor::A
end

Base.length(iter::EachFaceAround) = num_faces_around(iter.accessor)
Base.isdone(iter::EachFaceAround,face_around) = face_around > length(iter)
Base.getindex(iter::EachFaceAround,face_around) = at_face_around(iter.accessor,face_around)
Base.keys(iter::EachFaceAround) = LinearIndices((length(iter),))

function Base.iterate(iter::EachFaceAround,face_around=1)
    if Base.isdone(iter,face_around)
        nothing
    else
        accessor = iter[face_around]
        (accessor,face_around+1)
    end
end

function each_point(a::NewAbstractAccessor)
    EachPoint(a)
end

struct EachPoint{A} <: AbstractType
    accessor::A
end

Base.length(iter::EachPoint) = num_points(iter.accessor)
Base.isdone(iter::EachPoint,point) = point > length(iter)
Base.getindex(iter::EachPoint,point) = at_point(iter.accessor,point)
Base.keys(iter::EachPoint) = LinearIndices((length(iter),))

function Base.iterate(iter::EachPoint,point=1)
    if Base.isdone(iter,point)
        nothing
    else
        accessor = iter[point]
        (accessor,point+1)
    end
end

function reference_space_accessor(space::AbstractSpace,quadrature::AbstractQuadrature;compute=(),tabulate=())
    domain = GT.domain(quadrature)
    mesh = GT.mesh(space)
    d = Val(num_dims(domain))
    D = Val(num_dims(GT.domain(space)))
    if val_parameter(D) == val_parameter(d)
        loop_case = AtInterior()
        dface_ldfaces = nothing
    else
        loop_case = AtSkeleton()
        topo = topology(mesh)
        dface_ldfaces = face_local_faces(topo,d,D)
    end
    values = nothing
    gradients = nothing
    jacobians = nothing
    face = -1
    dface = -1
    Dface = -1
    face_around = -1
    ldface = -1
    point = -1
    location = ReferenceSpaceAccessorLocation(
                                              face,
                                              dface,
                                              Dface,
                                              face_around,
                                              ldface,
                                              point,
                                             )
    workspace = ReferenceSpaceAccessorWorkspace(
                                                d,
                                                D,
                                                dface_ldfaces,
                                                values,
                                                gradients,
                                                jacobians
                                               )
    a = ReferenceSpaceAccessor(
                               loop_case,
                               space,
                               quadrature,
                               location,
                               workspace
                              )
    setup_accessor(a,tabulate,compute)
end

function setup_accessor(a0,tabulate,compute)
    if GT.value in tabulate
        a1 = GT.tabulate(GT.value,a0)
    else
        a1 = a0
    end
    if ForwardDiff.gradient in tabulate
        a2 = GT.tabulate(ForwardDiff.gradient,a1)
    else
        a2 = a1
    end
    if ForwardDiff.jacobian in tabulate
        a3 = GT.tabulate(ForwardDiff.jacobian,a2)
    else
        a3 = a2
    end
    if unit_normal in compute
        a4 = GT.compute(unit_normal,a3)
    else
        a4 = a3
    end
    if coordinate in compute
        a5 = GT.compute(coordinate,a4)
    else
        a5 = a4
    end
    if ForwardDiff.jacobian in compute
        a6 = GT.compute(ForwardDiff.jacobian,a5)
    else
        a6 = a5
    end
    a6
end

struct AtInterior end
struct AtSkeleton end

struct ReferenceSpaceAccessor{A,B,C,D,E} <: NewAbstractAccessor
    loop_case::A
    space::B
    quadrature::C
    location::D
    workspace::E
end

struct ReferenceSpaceAccessorLocation <: AbstractType
    face::Int
    dface::Int
    Dface::Int
    face_around::Int
    ldface::Int
    point::Int
end

struct ReferenceSpaceAccessorWorkspace{A,B,C,X,E,F} <: AbstractType
    d::Val{A}
    D::Val{B}
    dface_ldfaces::C
    values::X
    gradients::E
    jacobians::F
end

function replace_face_dface(a::ReferenceSpaceAccessor,face,dface)
    location = replace_face_dface(a.location,face,dface)
    replace_location(a,location)
end

function replace_face_dface_Dface(a::ReferenceSpaceAccessor,face,dface,Dface)
    location = replace_face_dface_Dface(a.location,face,dface,Dface)
    replace_location(a,location)
end

function replace_Dface_face_around_ldface(a::ReferenceSpaceAccessor,Dface,face_around,ldface)
    location = replace_Dface_face_around_ldface(a.location,Dface,face_around,ldface)
    replace_location(a,location)
end

function replace_point(a::ReferenceSpaceAccessor,point)
    location = replace_point(a.location,point)
    replace_location(a,location)
end

function replace_values(a::ReferenceSpaceAccessor,values)
    workspace = replace_values(a.workspace,values)
    replace_workspace(a,workspace)
end

function replace_gradients(a::ReferenceSpaceAccessor,gradients)
    workspace = replace_gradients(a.workspace,gradients)
    replace_workspace(a,workspace)
end

function replace_jacobians(a::ReferenceSpaceAccessor,jacobians)
    workspace = replace_jacobians(a.workspace,jacobians)
    replace_workspace(a,workspace)
end

function replace_location(a::ReferenceSpaceAccessor,location)
    ReferenceSpaceAccessor(
                           a.loop_case,
                           a.space,
                           a.quadrature,
                           location,
                           a.workspace
                          )
end

function replace_workspace(a::ReferenceSpaceAccessor,workspace)
    ReferenceSpaceAccessor(
                           a.loop_case,
                           a.space,
                           a.quadrature,
                           a.location,
                           workspace
                          )
end

function replace_face_dface(a::ReferenceSpaceAccessorLocation,face,dface)
    ReferenceSpaceAccessorLocation(
                                   face,
                                   dface,
                                   a.Dface,
                                   a.face_around,
                                   a.ldface,
                                   a.point,
                                  )
end

function replace_face_dface_Dface(a::ReferenceSpaceAccessorLocation,face,dface,Dface)
    ReferenceSpaceAccessorLocation(
                                   face,
                                   dface,
                                   Dface,
                                   a.face_around,
                                   a.ldface,
                                   a.point,
                                  )
end

function replace_Dface_face_around_ldface(a::ReferenceSpaceAccessorLocation,Dface,face_around,ldface)
    ReferenceSpaceAccessorLocation(
                                   a.face,
                                   a.dface,
                                   Dface,
                                   face_around,
                                   ldface,
                                   a.point,
                                  )
end

function replace_point(a::ReferenceSpaceAccessorLocation,point)
    ReferenceSpaceAccessorLocation(
                                   a.face,
                                   a.dface,
                                   a.Dface,
                                   a.face_around,
                                   a.ldface,
                                   point,
                                  )
end

function replace_values(a::ReferenceSpaceAccessorWorkspace,values)
    ReferenceSpaceAccessorWorkspace(
                                    a.d,
                                    a.D,
                                    a.dface_ldfaces,
                                    values,
                                    a.gradients,
                                    a.jacobians
                                   )
end

function replace_gradients(a::ReferenceSpaceAccessorWorkspace,gradients)
    ReferenceSpaceAccessorWorkspace(
                                    a.d,
                                    a.D,
                                    a.dface_ldfaces,
                                    a.values,
                                    gradients,
                                    a.jacobians
                                   )
end

function replace_jacobians(a::ReferenceSpaceAccessorWorkspace,jacobians)
    ReferenceSpaceAccessorWorkspace(
                                    a.d,
                                    a.D,
                                    a.dface_ldfaces,
                                    a.values,
                                    a.gradients,
                                    jacobians
                                   )
end

function num_faces(a::ReferenceSpaceAccessor)
    (;quadrature) = a
    domain = GT.domain(quadrature)
    num_faces(domain)
end

function at_face(a::ReferenceSpaceAccessor{AtInterior},face)
    (;quadrature) = a
    domain = GT.domain(quadrature)
    face_dface = GT.faces(domain)
    dface = face_dface[face]
    Dface = dface
    replace_face_dface_Dface(a,face,dface,Dface)
end

function at_face(a::ReferenceSpaceAccessor{AtSkeleton},face)
    (;space,quadrature) = a
    (;d,D) = a.workspace
    mesh = GT.mesh(space)
    domain = GT.domain(quadrature)
    face_dface = GT.faces(domain)
    dface = face_dface[face]
    a2 = replace_face_dface(a,face,dface)
    faces_around = GT.faces_around(domain)
    if faces_around !== nothing
        face_around = faces_around[face]
        a3 = at_face_around(a2,face_around)
    else
        a3 = a2
    end
    a3
end

struct AnyIndex end

function at_any_index(a::ReferenceSpaceAccessor{AtInterior})
    # TODO the asserts can be removed by playing with AnyIndex
    @assert num_faces(a) > 0
    a2 = at_face(a,1)
    @assert num_points(a2) > 0
    a3 = at_point(a2,1)
end

function at_any_index(a::ReferenceSpaceAccessor{AtSkeleton})
    # TODO the asserts can be removed by playing with AnyIndex
    @assert num_faces(a) > 0
    a2 = at_face(a,1)
    @assert num_faces_around(a2) > 0
    a3 = at_face_around(a2,1)
    @assert num_points(a3) > 0
    a4 = at_point(a3,1)
end

function num_faces_around(a::ReferenceSpaceAccessor)
    (;space) = a
    (;d,D) = a.workspace
    (;dface) = a.location
    mesh = GT.mesh(space)
    topo = topology(mesh)
    dface_Dfaces = face_incidence(topo,d,D)
    Dfaces = dface_Dfaces[dface]
    length(Dfaces)
end

function at_face_around(a::ReferenceSpaceAccessor,face_around_0)
    (;space,quadrature) = a
    (;d,D,dface_ldfaces) = a.workspace
    (;face,dface) = a.location
    domain = GT.domain(quadrature)
    mesh = GT.mesh(space)
    topo = topology(mesh)
    faces_around_permutation = GT.faces_around_permutation(domain)
    if faces_around_permutation !== nothing
        face_around = faces_around_permutation[face][face_around_0]
    else
        face_around = face_around_0
    end
    dface_Dfaces = face_incidence(topo,d,D)
    Dface = dface_Dfaces[dface][face_around]
    ldface = dface_ldfaces[dface][face_around]
    replace_Dface_face_around_ldface(a,Dface,face_around,ldface)
end

function dofs(a::ReferenceSpaceAccessor)
    (;space) = a
    (;Dface) = a.location
    Dface_dofs = face_dofs(space)
    Dface_dofs[Dface]
end

function num_dofs(a::ReferenceSpaceAccessor)
    length(nodes(a))
end

function tabulate(f,a::ReferenceSpaceAccessor{AtInterior})
    (;quadrature,space) = a
    (;D) = a.workspace
    rid_point_x = map(coordinates,reference_quadratures(quadrature))
    # NB the TODOs below can be solved by introducing an extra nesting level
    # TODO assumes same reference elems for integration and for interpolation
    rid_to_tab = map(rid_point_x,reference_spaces(space)) do point_to_x, refface
        collect(permutedims(tabulator(refface)(f,point_to_x))) # TODO fix this globally
    end
    replace_tabulators(f,a,rid_to_tab)
end

function tabulate(f,a::ReferenceSpaceAccessor{AtSkeleton})
    (;quadrature,space) = a
    (;D,d) = a.workspace
    mesh = GT.mesh(space)
    topo = topology(mesh)
    drid_refdface = reference_spaces(mesh,d)
    Drid_refDface = reference_spaces(mesh,D)
    Drid_reffe = reference_spaces(space)
    drid_point_x = map(coordinates,reference_quadratures(quadrature))
    # NB the TODOs below can be solved by introducing two extra nesting levels
    # TODO this assumes the same reffes for mesh and quadrature
    drid_Drid_ldface_perm_tab = map(drid_point_x,drid_refdface) do point_to_x,refdface
        # TODO this assumes the same reffes for mesh and interpolation
        map(Drid_reffe,Drid_refDface) do reffe,refDface
            ldface_perm_varphi = reference_map(refdface,refDface)
            map(ldface_perm_varphi) do perm_varphi
                map(perm_varphi) do varphi
                    point_to_q = varphi.(point_to_x)
                    collect(permutedims(tabulator(reffe)(f,point_to_q)))
                end
            end
        end
    end
    replace_tabulators(f,a,drid_Drid_ldface_perm_tab)
end

function replace_tabulators(::typeof(GT.value),a::ReferenceSpaceAccessor,values)
    replace_values(a,values)
end

function replace_tabulators(::typeof(ForwardDiff.gradient),a::ReferenceSpaceAccessor,gradients)
    replace_gradients(a,gradients)
end

function replace_tabulators(::typeof(ForwardDiff.jacobian),a::ReferenceSpaceAccessor,jacobians)
    replace_jacobians(a,jacobians)
end

function tabulator(f,a::ReferenceSpaceAccessor{AtInterior})
    (;space) = a
    (;D) = a.workspace
    (;Dface) = a.location
    Drid = face_reference_id(space)[Dface]
    Drid_tab = tabulators(f,a)
    tab = Drid_tab[Drid]
end

function tabulator(f,a::ReferenceSpaceAccessor{AtSkeleton})
    (;space,quadrature) = a
    (;Dface,dface,ldface) = a.location
    (;d,D) = a.workspace
    mesh = GT.mesh(space)
    dface_drid = face_reference_id(mesh,d)
    Dface_Drid = face_reference_id(mesh,D)
    #Dface_Drid = face_reference_id(space) # TODO
    Drid = Dface_Drid[Dface]
    drid = dface_drid[dface]
    topo = topology(mesh)
    Dface_ldface_perm = GT.face_permutation_ids(topo,D,d)
    perm = Dface_ldface_perm[Dface][ldface]
    drid_Drid_ldface_perm_tab = tabulators(f,a)
    tab = drid_Drid_ldface_perm_tab[drid][Drid][ldface][perm]
end

function tabulators(::typeof(GT.value),a::ReferenceSpaceAccessor)
    (;values) = a.workspace
    values
end

function tabulators(::typeof(ForwardDiff.gradient),a::ReferenceSpaceAccessor)
    (;gradients) = a.workspace
    gradients
end

function tabulators(::typeof(ForwardDiff.jacobian),a::ReferenceSpaceAccessor)
    (;jacobians) = a.workspace
    jacobians
end

function shape_functions(a::ReferenceSpaceAccessor)
    (;space) = a
    (;Dface) = a.location
    Dface_Drid = face_reference_id(space)
    Drid = Dface_Drid[Dface]
    Drid_rspace = reference_spaces(space)
    rspace = Drid_rspace[Drid]
    shape_functions(rspace)
end

function num_points(a::ReferenceSpaceAccessor)
    num_points(quadrature(a))
end

function quadrature(a::ReferenceSpaceAccessor)
    (;quadrature) = a
    (;dface) = a.location
    dface_drid = GT.face_reference_id(quadrature)
    drid_refqua = reference_quadratures(quadrature)
    drid = dface_drid[dface]
    drid_refqua[drid]
end

function at_point(a::ReferenceSpaceAccessor,point)
    replace_point(a,point)
end

function shape_functions(f,a::ReferenceSpaceAccessor)
    (;point) = a.location
    tab = tabulator(f,a)
    view(tab,:,point)
end

# Mesh accessor

function mesh_accessor(mesh::AbstractMesh,
    D=Val(num_dims(mesh)),
    domain=GT.domain(mesh,Val(val_parameter(D)));kwargs...)
    degree = 0
    quadrature = GT.quadrature(domain,degree)
    a = mesh_accessor(mesh,Val(val_parameter(D)),quadrature)
end

function mesh_accessor(mesh::AbstractMesh,D,quadrature::AbstractQuadrature;
        tabulate=(), compute=())
    space = mesh_space(mesh,D)
    space_accessor = reference_space_accessor(space,quadrature)
    loop_case = space_accessor.loop_case
    J = nothing
    reference_normals = nothing
    workspace = MeshAccessorWorkspace(J,reference_normals)
    a = MeshAccessor(loop_case,space_accessor,workspace)
    setup_accessor(a,tabulate,compute)
end

struct MeshAccessor{A,B,C} <: NewAbstractAccessor
    loop_case::A
    space_accessor::B
    workspace::C
end

struct MeshAccessorWorkspace{A,B}
    jacobian::A
    reference_normals::B
end

function replace_space_accessor(a::MeshAccessor,space_accessor)
    MeshAccessor(
                 a.loop_case,
                 space_accessor,
                 a.workspace)
end

function replace_workspace(a::MeshAccessor,workspace)
    MeshAccessor(
                 a.loop_case,
                 a.space_accessor,
                 workspace)
end

function replace_jacobian(a::MeshAccessor,jacobian)
    workspace = replace_jacobian(a.workspace,jacobian)
    replace_workspace(a,workspace)
end

function replace_jacobian(a::MeshAccessorWorkspace,jacobian)
    MeshAccessorWorkspace(
                          jacobian,
                          a.reference_normals
                         )
end

function replace_reference_normals(a::MeshAccessor,reference_normals)
    workspace = replace_reference_normals(a.workspace,reference_normals)
    replace_workspace(a,workspace)
end

function replace_reference_normals(a::MeshAccessorWorkspace,reference_normals)
    MeshAccessorWorkspace(
                          a.jacobian,
                          reference_normals
                         )
end

function num_faces(a::MeshAccessor)
    (;space_accessor) = a
    num_faces(space_accessor)
end

function at_face(a::MeshAccessor,face)
    space_accessor = at_face(a.space_accessor,face)
    replace_space_accessor(a,space_accessor)
end

function at_any_index(a::MeshAccessor{AtInterior})
    # TODO the asserts can be removed by playing with AnyIndex
    @assert num_faces(a) > 0
    a2 = at_face(a,1)
    @assert num_points(a2) > 0
    a3 = at_point(a2,1)
end

function at_any_index(a::MeshAccessor{AtSkeleton})
    # TODO the asserts can be removed by playing with AnyIndex
    @assert num_faces(a) > 0
    a2 = at_face(a,1)
    @assert num_faces_around(a2) > 0
    a3 = at_face_around(a2,1)
    @assert num_points(a3) > 0
    a4 = at_point(a3,1)
end

function num_faces_around(a::MeshAccessor)
    (;space_accessor) = a
    num_faces_around(space_accessor)
end

function at_face_around(a::MeshAccessor,face_around)
    space_accessor = at_face_around(a.space_accessor,face_around)
    replace_space_accessor(a,space_accessor)
end

function nodes(a::MeshAccessor)
    (;space_accessor) = a
    dofs(space_accessor)
end

function num_nodes(a::MeshAccessor)
    length(nodes(a))
end

function node_coordinates(a::MeshAccessor)
    (;space_accessor) = a
    (;space) = space_accessor
    mesh = GT.mesh(space) 
    node_x = node_coordinates(mesh)
    view(node_x,nodes(a))
end

function tabulate(f,a::MeshAccessor)
    space_accessor = tabulate(f,a.space_accessor)
    replace_space_accessor(a,space_accessor)
end

function tabulate(f::typeof(ForwardDiff.gradient),a::MeshAccessor)
    space_accessor = tabulate(f,a.space_accessor)
    space_accessor_2 = at_any_index(space_accessor)
    mesh = GT.mesh(space_accessor.space)
    x = zero(eltype(node_coordinates(mesh)))
    g = zero(eltype(shape_functions(f,space_accessor_2)))
    J = zero(outer(x,g))
    a2 = replace_space_accessor(a,space_accessor)
    replace_jacobian(a2,J)
end

function unit_normal end

function compute(::typeof(unit_normal),a::MeshAccessor{AtSkeleton})
    (;space_accessor) = a
    (;space) = space_accessor
    (;D,d) = space_accessor.workspace
    mesh = GT.mesh(space)
    @assert val_parameter(d) + 1 == val_parameter(D)
    Drid_ldface_nref = map(GT.reference_spaces(mesh,D)) do refface
        GT.normals(GT.mesh(GT.domain(refface)))# TODO rename normals?
    end
    replace_reference_normals(a,Drid_ldface_nref)
end

function coordinate end

function compute(::typeof(GT.coordinate),a::MeshAccessor)
    tabulate(GT.value,a)
end

function compute(::typeof(ForwardDiff.jacobian),a::MeshAccessor)
    tabulate(ForwardDiff.gradient,a)
end

function shape_functions(a::MeshAccessor)
    (;space_accessor) = a
    shape_functions(space_accessor)
end

function diameter(a::MeshAccessor)
    nodes = GT.nodes(a)
    mesh = GT.mesh(a.space_accessor.space)
    node_x = GT.node_coordinates(mesh)
    nnodes = length(nodes)
    diam = zero(eltype(eltype(node_x)))
    for i in 1:nnodes
        xi = node_x[nodes[i]]
        for j in 1:nnodes
            xj = node_x[nodes[j]]
            diam = max(diam,norm(xi-xj))
        end
    end
    diam
end

function num_points(a::MeshAccessor)
    (;space_accessor) = a
    num_points(space_accessor)
end

function quadrature(a::MeshAccessor)
    (;space_accessor) = a
    quadrature(space_accessor)
end

function at_point(a::MeshAccessor,point)
    space_accessor = at_point(a.space_accessor,point)
    a2 = replace_space_accessor(a,space_accessor)
    compute_jacobian(a2)
end

# Do not remove the @noinline
# it seems to be performance relevant
@noinline function sum_jacobian(node_x,i_node,i_s,n,J0)
    J = zero(J0)
    i = 0
    while i < n
        i += 1
        J += outer(node_x[i_node[i]],i_s[i])
    end
    J
    ##TODO sum leads to much faster than hand-written loop, but why?
    #the answer was the noinline above. But why?
    #J = sum(i->outer(node_x[i_node[i]],i_s[i]),1:n;init=zero(J0))
end

 function compute_jacobian(a::MeshAccessor)
    J0 = a.workspace.jacobian
    if J0 === nothing
        return a
    end
    space_accessor = a.space_accessor
    point = space_accessor.location.point
    i_s = shape_functions(ForwardDiff.gradient,space_accessor)
    i_node = nodes(a)
    mesh = GT.mesh(space_accessor.space)
    node_x = node_coordinates(mesh)
    n = length(i_s)
    J = sum_jacobian(node_x,i_node,i_s,n,J0)
    replace_jacobian(a,J)
end

function shape_functions(f,a::MeshAccessor)
    (;space_accessor) = a
    shape_functions(f,space_accessor)
end

function coordinate(a::MeshAccessor)
    coordinate(GT.value,a)
end

function ForwardDiff.jacobian(a::MeshAccessor)
    coordinate(ForwardDiff.jacobian,a)
end

function coordinate(::typeof(value),a::MeshAccessor)
    s = shape_functions(value,a)
    x = node_coordinates(a)
    n = length(s)
    sum(i->x[i]*s[i],1:n)
end

function coordinate(::typeof(ForwardDiff.jacobian),a::MeshAccessor)
    J = a.workspace.jacobian
    @assert J !== nothing
    J
    #s = shape_functions(ForwardDiff.gradient,a)
    #x = node_coordinates(a)
    #n = num_nodes(a)
    #sum(i->outer(x[i],s[i]),1:n)
end

function weight(a::MeshAccessor)
    (;space_accessor) = a
    (;point) = space_accessor.location
    qua = quadrature(a)
    w = weights(qua)[point]
    J = coordinate(ForwardDiff.jacobian,a)
    change_of_measure(J)*w
end

function unit_normal(a::MeshAccessor)
    (;space_accessor) = a
    (;reference_normals,) = a.workspace
    (;space) = space_accessor
    (;D) = space_accessor.workspace
    (;Dface,ldface) = space_accessor.location
    Drid_ldface_nref = reference_normals
    mesh = GT.mesh(space)
    Dface_Drid = face_reference_id(mesh,D)
    Drid = Dface_Drid[Dface]
    nref = Drid_ldface_nref[Drid][ldface]
    J = coordinate(ForwardDiff.jacobian,a)
    nphys = map_unit_normal(J,nref)
    nphys
end

function map_unit_normal(J,n)
    Jt = transpose(J)
    pinvJt = transpose(inv(Jt*J)*Jt)
    v = pinvJt*n
    m = sqrt(v⋅v)
    if m < eps()
        return zero(v)
    else
        return v/m
    end
end

function quadrature_accessor(quadrature::AbstractQuadrature)
    domain = GT.domain(quadrature)
    mesh = GT.mesh(domain)
    d = Val(num_dims(domain))
    acc1 = mesh_accessor(mesh,d,quadrature)
    acc2 = tabulate(GT.value,acc1)
    acc3 = tabulate(ForwardDiff.gradient,acc2)
    acc3
end

struct SpaceAccessor{A,B,C,D,E} <: NewAbstractAccessor
    loop_case::A
    space::B
    mesh_accessor::C
    reference_space_accessor::D
    workspace::E
end

struct SpaceAccessorWorkspace{A,B,C} <: AbstractType
    values::A
    gradients::B
    jacobians::C
end

function replace_mesh_accessor(a::SpaceAccessor,mesh_accessor)
    SpaceAccessor(
                  a.loop_case,
                  a.space,
                  mesh_accessor,
                  a.reference_space_accessor,
                  a.workspace,
                 )
end

function replace_reference_space_accessor(a::SpaceAccessor,reference_space_accessor)
    SpaceAccessor(
                  a.loop_case,
                  a.space,
                  a.mesh_accessor,
                  reference_space_accessor,
                  a.workspace,
                 )
end

function replace_workspace(a::SpaceAccessor,workspace)
    SpaceAccessor(
                  a.loop_case,
                  a.space,
                  a.mesh_accessor,
                  a.reference_space_accessor,
                  workspace,
                 )
end

function replace_values(a::SpaceAccessor,values)
    workspace = replace_values(a.workspace,values)
    replace_workspace(a,workspace)
end

function replace_gradients(a::SpaceAccessor,gradients)
    workspace = replace_gradients(a.workspace,gradients)
    replace_workspace(a,workspace)
end

function replace_jacobians(a::SpaceAccessor,jacobians)
    workspace = replace_jacobians(a.workspace,jacobians)
    replace_workspace(a,workspace)
end

function replace_values(a::SpaceAccessorWorkspace,values)
    SpaceAccessorWorkspace(
                           values,
                           a.gradients,
                           a.jacobians,
                          )
end

function replace_gradients(a::SpaceAccessorWorkspace,gradients)
    SpaceAccessorWorkspace(
                           a.values,
                           gradients,
                           a.jacobians,
                          )
end

function replace_jacobians(a::SpaceAccessorWorkspace,jacobians)
    SpaceAccessorWorkspace(
                           a.values,
                           a.gradients,
                           jacobians,
                          )
end

function space_accessor(space::AbstractSpace,domain::AbstractDomain=GT.domain(space);kwargs...)
    degree = 0
    quadrature = GT.quadrature(domain,degree)
    space_accessor(space,quadrature;kwargs...)
end

function space_accessor(space::AbstractSpace,quadrature::AbstractQuadrature;kwargs...)
    D = num_dims(GT.domain(space))
    mesh = GT.mesh(space)
    mesh_accessor = GT.mesh_accessor(mesh,Val(D),quadrature)
    mesh_accessor_2 = tabulate(ForwardDiff.gradient,mesh_accessor)
    space_accessor(space,mesh_accessor_2;kwargs...)
end

function space_accessor(space::AbstractSpace,faces::EachFace;kwargs...)
    space_accessor(space,faces.accessor;kwargs...)
end

function space_accessor(space::AbstractSpace,mesh_accessor::MeshAccessor;
        tabulate=(), compute = ())
    quadrature = mesh_accessor.space_accessor.quadrature
    reference_space_accessor = GT.reference_space_accessor(space,quadrature)
    loop_case = reference_space_accessor.loop_case
    values = nothing
    gradients = nothing
    jacobians = nothing
    workspace = SpaceAccessorWorkspace(
                                       values,
                                       gradients,
                                       jacobians
                                      )
    a = SpaceAccessor(
                      loop_case,
                      space,
                      mesh_accessor,
                      reference_space_accessor,
                      workspace)
    setup_accessor(a,tabulate,compute)
end

function tabulate(f,a::SpaceAccessor{AtInterior})
    (;space,mesh_accessor) = a
    reference_space_accessor = tabulate(f,a.reference_space_accessor)
    reference_space_accessor_2 = at_any_index(reference_space_accessor)
    mesh_accessor_2 = at_any_index(mesh_accessor)
    dof_sref = shape_functions(f,reference_space_accessor_2)
    sref = zero(eltype(dof_sref))
    dof = AnyIndex()
    sphys = map_shape_function(f,space,dof,mesh_accessor_2,sref)
    ndofs = max_num_reference_dofs(space)
    dof_sphys = zeros(typeof(sphys),ndofs)
    a2 = replace_reference_space_accessor(a,reference_space_accessor)
    replace_workspace(f,a2,dof_sphys)
end

function tabulate(f,a::SpaceAccessor{AtSkeleton})
    (;space,mesh_accessor) = a
    reference_space_accessor = tabulate(f,a.reference_space_accessor)
    reference_space_accessor_2 = at_any_index(reference_space_accessor)
    mesh_accessor_2 = at_any_index(mesh_accessor)
    dof_sref = shape_functions(f,reference_space_accessor_2)
    sref = zero(eltype(dof_sref))
    dof = AnyIndex()
    sphys = map_shape_function(f,space,dof,mesh_accessor_2,sref)
    ndofs = max_num_reference_dofs(space)
    max_num_faces_around = 2 # TODO
    nfa = max_num_faces_around
    face_around_dof_sphys = [ zeros(typeof(sphys),ndofs) for _ in 1:nfa]
    a2 = replace_reference_space_accessor(a,reference_space_accessor)
    replace_workspace(f,a2,face_around_dof_sphys)
end

function compute(::typeof(coordinate),a::SpaceAccessor)
    mesh_accessor = tabulate(GT.value,a.mesh_accessor)
    replace_mesh_accessor(a,mesh_accessor)
end

function compute(::typeof(unit_normal),a::SpaceAccessor)
    mesh_accessor = compute(unit_normal,a.mesh_accessor)
    replace_mesh_accessor(a,mesh_accessor)
end

function replace_workspace(::typeof(GT.value),a::SpaceAccessor,values)
    replace_values(a,values)
end

function replace_workspace(::typeof(ForwardDiff.gradient),a::SpaceAccessor,gradients)
    replace_gradients(a,gradients)
end

function replace_workspace(::typeof(ForwardDiff.jacobian),a::SpaceAccessor,jacobians)
    replace_jacobians(a,jacobians)
end

function workspace(::typeof(GT.value),a::SpaceAccessor)
    a.workspace.values
end

function workspace(::typeof(ForwardDiff.gradient),a::SpaceAccessor)
    a.workspace.gradients
end

function workspace(::typeof(ForwardDiff.jacobian),a::SpaceAccessor)
    a.workspace.jacobians
end

function num_faces(a::SpaceAccessor)
    num_faces(a.mesh_accessor)
end

function at_face(a::SpaceAccessor,face)
    mesh_accessor = at_face(a.mesh_accessor,face)
    reference_space_accessor = at_face(a.reference_space_accessor,face)
    a2 = replace_mesh_accessor(a,mesh_accessor)
    replace_reference_space_accessor(a2,reference_space_accessor)
end

function num_faces_around(a::SpaceAccessor)
    num_faces_around(a.mesh_accessor)
end

function at_face_around(a::SpaceAccessor,face)
    mesh_accessor = at_face_around(a.mesh_accessor,face)
    reference_space_accessor = at_face_around(a.reference_space_accessor,face)
    a2 = replace_mesh_accessor(a,mesh_accessor)
    replace_reference_space_accessor(a2,reference_space_accessor)
end

function dofs(a::SpaceAccessor)
    dofs(a.reference_space_accessor)
end

function num_dofs(a::SpaceAccessor)
    length(dofs(a))
end

function num_points(a::SpaceAccessor)
    num_points(a.mesh_accessor)
end

function unit_normal(a::SpaceAccessor)
    unit_normal(a.mesh_accessor)
end

function at_point(a::SpaceAccessor,point)
    mesh_accessor = at_point(a.mesh_accessor,point)
    reference_space_accessor = at_point(a.reference_space_accessor,point)
    a2 = replace_mesh_accessor(a,mesh_accessor)
    replace_reference_space_accessor(a2,reference_space_accessor)
end

function at_point(a::SpaceAccessor,mesh_accessor::MeshAccessor)
    location = mesh_accessor.space_accessor.location
    reference_space_accessor = replace_location(a.reference_space_accessor,location)
    a2 = replace_mesh_accessor(a,mesh_accessor)
    replace_reference_space_accessor(a2,reference_space_accessor)
end


@inline function shape_functions(f,a::SpaceAccessor{AtInterior})
    (;space,mesh_accessor,reference_space_accessor) = a
    @show reference_space_accessor.location.point
    dof_sref = GT.shape_functions(f,reference_space_accessor)
    dof_sphys = workspace(f,a)
    ndofs = length(dof_sref)
    dof = 0
    while dof < ndofs
        dof += 1
        sref = dof_sref[dof]
        sphys = map_shape_function(f,space,dof,mesh_accessor,sref)
        dof_sphys[dof] = sphys
    end
    dof_sphys
    view(dof_sphys,1:ndofs)
end

function shape_functions(f,a::SpaceAccessor{AtSkeleton})
    (;space,mesh_accessor,reference_space_accessor) = a
    face_around = mesh_accessor.space_accessor.location.face_around
    dof_sref = GT.shape_functions(f,reference_space_accessor)
    dof_sphys = workspace(f,a)[face_around]
    ndofs = length(dof_sref)
    for dof in 1:ndofs
        sref = dof_sref[dof]
        sphys = map_shape_function(f,space,dof,mesh_accessor,sref)
        dof_sphys[dof] = sphys
    end
    dof_sphys
    view(dof_sphys,1:ndofs)
end

function map_shape_function(::typeof(GT.value),space,dof,mesh_accessor,sref)
    sref
end

@inline function map_shape_function(::typeof(ForwardDiff.gradient),space,dof,mesh_accessor,sref)
    J = jacobian(mesh_accessor)
    sphys = transpose(J)\sref
end

function map_shape_function(::typeof(ForwardDiff.jacobian),space,dof,mesh_accessor,sref)
    J = coordinate(ForwardDiff.jacobian,mesh_accessor)
    sphys = sref/J
end

function weight(a::SpaceAccessor)
    (;mesh_accessor) = a
    weight(mesh_accessor)
end

function coordinate(a::SpaceAccessor)
    (;mesh_accessor) = a
    coordinate(mesh_accessor)
end

struct NewDiscreteFieldAccessor{A,B,C,D} <: NewAbstractAccessor
    loop_case::A
    field::B
    space_accessor::C
    workspace::D
end

function replace_workspace(a::NewDiscreteFieldAccessor,workspace)
    NewDiscreteFieldAccessor(
                             a.loop_case,
                             a.field,
                             a.space_accessor,
                             workspace,
                            )
end

function replace_space_accessor(a::NewDiscreteFieldAccessor,space_accessor)
    NewDiscreteFieldAccessor(
                             a.loop_case,
                             a.field,
                             space_accessor,
                             a.workspace,
                            )
end

function field_accessor(field::DiscreteField,args...;kwargs...)
    space = GT.space(field)
    space_accessor = GT.space_accessor(space,args...;kwargs...)
    loop_case = space_accessor.loop_case
    workspace = nothing
    a = NewDiscreteFieldAccessor(loop_case,field,space_accessor,workspace)
    setup_workspace(a)
end

function setup_workspace(a::NewDiscreteFieldAccessor{AtInterior})
    (;field) = a
    space = GT.space(field)
    n = max_num_reference_dofs(space)
    T = eltype(free_values(field))
    dof_value = zeros(T,n)
    workspace = dof_value
    replace_workspace(a,workspace)
end

function setup_workspace(a::NewDiscreteFieldAccessor{AtSkeleton})
    (;field) = a
    space = GT.space(field)
    n = max_num_reference_dofs(space)
    max_num_faces_around = 2 # TODO
    m = max_num_faces_around
    T = eltype(free_values(space))
    face_around_dof_value = [zeros(T,n) for _ in 1:m]
    workspace = face_around_dof_value
    replace_workspace(a,workspace)
end

function num_faces(a::NewDiscreteFieldAccessor)
    (;space_accessor) = a
    num_faces(space_accessor)
end

function at_face(a::NewDiscreteFieldAccessor,face)
    space_accessor = at_face(a.space_accessor,face)
    replace_space_accessor(a,space_accessor)
end

function at_any_index(a::NewDiscreteFieldAccessor)
    space_accessor = at_any_index(a.space_accessor)
    replace_space_accessor(a,space_accessor)
end

function num_faces_around(a::NewDiscreteFieldAccessor)
    (;space_accessor) = a
    num_faces_around(space_accessor)
end

function at_face_around(a::NewDiscreteFieldAccessor,face_around)
    space_accessor = at_face_around(a.space_accessor,face_around)
    replace_space_accessor(a,space_accessor)
end

function dofs(a::NewDiscreteFieldAccessor)
    (;space_accessor) = a
    dofs(space_accessor)
end

function num_dofs(a::NewDiscreteFieldAccessor)
    length(dofs(a))
end

function values(a::NewDiscreteFieldAccessor{AtInterior})
    (;space_accessor,field,workspace) = a
    free_dof_value = free_values(field)
    diri_dof_value = dirichlet_values(field)
    idof_value = workspace
    dofs = GT.dofs(a)
    ndofs = length(dofs)
    for idof in 1:ndofs
        dof = dofs[idof]
        if dof > 0
            free_dof = dof
            value = free_dof_value[dof]
        else
            diri_dof = -dof
            value = diri_dof_value[diri_dof]
        end
        idof_value[idof] = value
    end
    view(idof_value,1:ndofs)
end

function values(a::NewDiscreteFieldAccessor{AtSkeleton})
    (;space_accessor,field,workspace) = a
    free_dof_value = free_values(field)
    diri_dof_value = dirichlet_values(field)
    face_around = GT.face_around(a)
    idof_value = workspace[face_around]
    dofs = GT.dofs(a)
    ndofs = length(dofs)
    for idof in 1:ndofs
        dof = dofs[idof]
        if dof > 0
            free_dof = dof
            value = free_dof_value[dof]
        else
            diri_dof = -dof
            value = diri_dof_value[diri_dof]
        end
        idof_value[idof] = value
    end
    view(idof_value,1:ndofs)
end

function tabulate(f,a::NewDiscreteFieldAccessor)
    space_accessor = tabulate(f,aspace_accessor)
    replace_space_accessor(a,space_accessor)
end

function compute(f,a::NewDiscreteFieldAccessor)
    space_accessor = compute(f,a.space_accessor)
    replace_space_accessor(a,space_accessor)
end

function num_points(a::NewDiscreteFieldAccessor)
    num_points(a.space_accessor)
end

function at_point(a::NewDiscreteFieldAccessor,point)
    space_accessor = at_point(a.space_accessor,point)
    replace_space_accessor(a,space_accessor)
end

function shape_functions(f,a::NewDiscreteFieldAccessor)
    (;space_accessor) = a
    shape_functions(f,space_accessor)
end

function field(f,a::NewDiscreteFieldAccessor)
    s = shape_functions(f,a)
    x = values(a)
    n = num_dofs(a)
    sum(i->x[i]*s[i],1:n)
end

function weight(a::NewDiscreteFieldAccessor)
    (;space_accessor) = a
    weight(space_accessor)
end

function coordinate(a::NewDiscreteFieldAccessor)
    (;space_accessor) = a
    coordinate(space_accessor)
end

struct AnalyticalFieldAccessor{L,A,B} <: NewAbstractAccessor
    loop_case::L
    field::A
    mesh_accessor::B
end

function field_accessor(field::AnalyticalField,quadrature;kwargs...)
    domain = GT.domain(field)
    mesh = GT.mesh(domain)
    D = num_dims(domain)
    mesh_accessor = GT.mesh_accessor(mesh,Val(D),quadrature;kwargs...)
    loop_case = mesh_accessor.loop_case
    AnalyticalFieldAccessor(loop_case,field,mesh_accessor)
end

function num_faces(a::AnalyticalFieldAccessor)
    (;mesh_accessor) = a
    num_faces(mesh_accessor)
end

function at_face(a::AnalyticalFieldAccessor,face)
    mesh_accessor = at_face(a.mesh_accessor,face)
    replace_mesh_accessor(a,mesh_accessor)
end

function at_any_index(a::AnalyticalFieldAccessor)
    mesh_accessor = at_any_index(a.mesh_accessor)
    replace_mesh_accessor(a,mesh_accessor)
end

function num_faces_around(a::AnalyticalFieldAccessor)
    (;mesh_accessor) = a
    num_faces_around(mesh_accessor)
end

function at_face_around(a::AnalyticalFieldAccessor,face_around)
    mesh_accessor = at_face_around(a.mesh_accessor,face_around)
    replace_mesh_accessor(a,mesh_accessor)
end

function at_point(a::AnalyticalFieldAccessor,point)
    mesh_accessor = at_point(a.mesh_accessor,point)
    replace_mesh_accessor(a,mesh_accessor)
end

function field(a::AnalyticalFieldAccessor)
    location = a.field.location
    def = f.field.definition
    if location == nothing
        def
    else
        group_names = GT.group_names(a.field.domain)
        Dface = a.mesh_accessor.space_accessor.location
        loc = location[Dface]
        name = group_names
        y->def(y,name)
    end
end

function field(f,a::AnalyticalFieldAccessor)
    x = coordinate(a.mesh_accessor)
    g = field(a)
    f(g,x)
end

function weight(a::AnalyticalFieldAccessor)
    (;mesh_accessor) = a
    weight(mesh_accessor)
end

function coordinate(a::AnalyticalFieldAccessor)
    (;mesh_accessor) = a
    coordinate(mesh_accessor)
end

#### Old stuff
###
###function accessor(a,b)
###    Accessor(a,b)
###end
###
###struct Accessor{A,B} <: AbstractAccessor
###    definition::A
###    prototype::B
###end
###
###function prototype(a::Accessor)
###    a.prototype
###end
###
###function (a::Accessor)(face,face_around=nothing)
###    a.definition(face,face_around)
###end
###
###"""
###    reference_space_accessor(g,mesh,d)
###
###Return an accessor object `face_g` giving access to a quantity
###in the reference space of dimension `d` in the mesh object `mesh`.
###Calling `face_g(face)` on a face id `face` is equivalent to calling
###`g(ref_space)` with `ref_space = GT.reference_spaces(mesh,d)[r]` and
###`r=GT.face_reference_id(mesh,d)[d]`.
###
#### Level
###
###Intermediate
###
###"""
###function reference_space_accessor end
###
###function reference_space_accessor(mesh,d)
###    reference_space_accessor(identity,mesh,d)
###end
###
###function reference_space_accessor(f,mesh,d)
###    rid_to_space = reference_spaces(mesh,d)
###    rid_to_v = map(f,rid_to_space)
###    face_to_rid = GT.face_reference_id(mesh,d)
###    function face_to_v(face,face_around)
###        rid = face_to_rid[face]
###        v = rid_to_v[rid]
###    end
###    prototype = first(rid_to_v)
###    accessor(face_to_v,prototype)
###end
###
###"""
###    reference_topology_accessor(g,topo,d)
###
###Return an accessor object `face_g` giving access to a quantity
###in the reference topology of dimension `d` in the object `topo`.
###Calling `face_g(face)` on a face id `face` is equivalent to calling
###`g(ref_topo)` with `ref_topo = GT.reference_topologies(topo,d)[r]` and
###`r=GT.face_reference_id(topo,d)[d]`.
###
#### Level
###
###Intermediate
###
###"""
###function reference_topology_accessor end
###
###function reference_topology_accessor(mesh,d)
###    reference_topology_accessor(identity,mesh,d)
###end
###
###function reference_topology_accessor(f,mesh,d)
###    rid_to_topology = reference_topologies(mesh,d)
###    rid_to_v = map(f,rid_to_topology)
###    face_to_rid = GT.face_reference_id(mesh,d)
###    function face_to_v(face,face_around)
###        rid = face_to_rid[face]
###        v = rid_to_v[rid]
###    end
###    prototype = first(rid_to_v)
###    accessor(face_to_v,prototype)
###end
###
###"""
###    node_coordinate_accessor(mesh,d)
###
###Return an accessor object `face_lnode_x` that gives access to the coordinate vector
###`x` of the local node `lnode` of face `f` in dimension `d` in the object `mesh` as 
###`x = face_lnode_x(face)(lnode)`.
###
###See also [`reference_space_accessor`](@ref).
###
#### Level
###
###Intermediate
###"""
###function node_coordinate_accessor(mesh::AbstractMesh,d)
###    node_x = GT.node_coordinates(mesh)
###    face_nodes = GT.face_nodes(mesh,d)
###    function face_lnode_x(face,face_around)
###        nodes = face_nodes[face]
###        function lnode_x(lnode)
###            node = nodes[lnode]
###            x = node_x[node]
###        end
###    end
###    prototype = zero(eltype(node_x))
###    accessor(face_lnode_x,prototype)
###end
###
###
#### Untabulated
###function shape_function_accessor(f,space::AbstractSpace)
###    if is_reference_domain(domain(space))
###        shape_function_accessor_reference(f,space)
###    else
###        shape_function_accessor_physical(f,space)
###    end
###end
###
#### Tabulated
###function shape_function_accessor(f,space::AbstractSpace,measure::AbstractQuadrature,disable_tabulation=Val(false))
###    if is_reference_domain(domain(space))
###        # TODO: do we need to remove the tabulation or inline this?
###        shape_function_accessor_reference(f,space,measure)
###    else
###        shape_function_accessor_physical(f,space,measure,disable_tabulation)
###    end
###end
###
#### Untabulated
###function shape_function_accessor_reference(f,space::AbstractSpace)
###    rid_to_dof_s = map(shape_functions,reference_spaces(space))
###    face_to_rid = face_reference_id(space)
###    function face_dof_g(face,face_around=nothing)
###        rid = face_to_rid[face]
###        dof_s = rid_to_dof_s[rid]
###        function dof_g(dof)
###            s = dof_s[dof]
###            function g(x)
###                f(s,x)
###            end
###        end
###    end
###    prototype = first(first(rid_to_dof_s))
###    accessor(face_dof_g,prototype)
###end
###
#### Tabulated
###function shape_function_accessor_reference(f,space::AbstractSpace,measure::AbstractQuadrature)
###    domain_measure = domain(measure)
###    domain_space = domain(space)
###    d = num_dims(domain_measure)
###    D = num_dims(domain_space)
###    faces_around = GT.faces_around(domain_measure)
###    if d == D
###        shape_function_accessor_reference_interior(f,space,measure)
###    elseif d+1==D && faces_around !== nothing
###        shape_function_accessor_reference_boundary(f,space,measure)
###    else
###        shape_function_accessor_reference_skeleton(f,space,measure)
###    end
###end
###
###function shape_function_accessor_reference_interior(f,space::AbstractSpace,measure::AbstractQuadrature)
###    dom = GT.domain(measure)
###    mesh = GT.mesh(dom)
###    d = num_dims(dom)
###    rid_to_point_to_x = map(coordinates,reference_quadratures(measure))
###    # NB the TODOs below can be solved by introducing an extra nesting level
###    # TODO assumes same reference elems for integration and for interpolation
###    rid_to_tab = map(rid_to_point_to_x,reference_spaces(space)) do point_to_x, refface
###        # tabulator(refface)(f,point_to_x)
###        collect(permutedims(tabulator(refface)(f,point_to_x), (2, 1)))
###    end
###    prototype = first(first(rid_to_tab))
###    face_to_rid = face_reference_id(mesh,d)
###    sface_to_face = faces(dom)
###    function face_point_dof_s(sface,face_around=nothing)
###        face = sface_to_face[sface]
###        rid = face_to_rid[face]
###        tab = rid_to_tab[rid]
###        function point_dof_s(point,J=nothing)
###            function dof_s(dof)
###                tab[dof, point]
###                # tab[point, dof]
###            end
###        end
###    end
###    accessor(face_point_dof_s,prototype)
###end
###
###function shape_function_accessor_reference_skeleton(f,space::AbstractSpace,measure::AbstractQuadrature)
###    d = num_dims(domain(measure))
###    D = num_dims(domain(space))
###    mesh = GT.mesh(space)
###    topo = topology(mesh)
###    drid_to_refdface = reference_spaces(mesh,Val(d))
###    Drid_to_refDface = reference_spaces(mesh,Val(D))
###    Drid_to_reffe = reference_spaces(space)
###    drid_to_point_to_x = map(coordinates,reference_quadratures(measure))
###    dface_to_Dfaces, dface_to_ldfaces = GT.face_incidence_ext(topo,d,D)
###    Dface_to_ldface_to_perm = GT.face_permutation_ids(topo,D,d)
###    face_to_dface = faces(domain(measure))
###    Dface_to_Drid = face_reference_id(mesh,D)
###    dface_to_drid = face_reference_id(mesh,d)
###    # NB the TODOs below can be solved by introducing two extra nesting levels
###    # TODO this assumes the same reffes for mesh and quadrature
###
###    point_to_x = drid_to_point_to_x[1]
###    refdface = drid_to_refdface[1]
###    reffe,refDface = Drid_to_reffe[1], Drid_to_refDface[1]
###    lpv = reference_map(refdface,refDface)[1][1]
###    p2q = lpv.(point_to_x)
###    prototype_tab = tabulator(reffe)(f,p2q)
###    prototype = first(prototype_tab)
###    prototype_inner = Vector{Vector{typeof(prototype_tab)}}
###    prototype_outer = Vector{Vector{prototype_inner}}
###
###    drid_Drid_ldface_perm_to_tab::prototype_outer = map(drid_to_point_to_x,drid_to_refdface,1:length(drid_to_point_to_x)) do point_to_x,refdface,_ # fix the type to vector by adding a loop range
###        # TODO this assumes the same reffes for mesh and interpolation
###        map(Drid_to_reffe,Drid_to_refDface,1:length(Drid_to_reffe)) do reffe,refDface,_
###            ldface_perm_varphi = reference_map(refdface,refDface)
###            map(ldface_perm_varphi) do perm_varphi
###                map(perm_varphi) do varphi
###                    point_to_q = varphi.(point_to_x)
###                    # elem::typeof(prototype_tab) = tabulator(reffe)(f,point_to_q)
###                    elem::typeof(prototype_tab) = collect(permutedims(tabulator(reffe)(f,point_to_q), (2, 1)))
###                end
###            end
###        end
###    end
###    
###    function face_point_dof_s(face,face_around)
###        dface = face_to_dface[face]
###        Dfaces = dface_to_Dfaces[dface]
###        ldfaces = dface_to_ldfaces[dface]
###        Dface = Dfaces[face_around]
###        ldface = ldfaces[face_around]
###        drid = dface_to_drid[dface]
###        Drid = Dface_to_Drid[Dface]
###        perm = Dface_to_ldface_to_perm[Dface][ldface]
###        tab = drid_Drid_ldface_perm_to_tab[drid][Drid][ldface][perm]
###        function point_dof_s(point)
###            function dof_s(dof)
###                # tab[point,dof]
###                tab[dof, point]
###            end
###        end
###    end
###    accessor(face_point_dof_s,prototype)
###end
###
###function shape_function_accessor_reference_boundary(f,space::AbstractSpace,measure::AbstractQuadrature)
###    faces_around = GT.faces_around(domain(measure))
###    face_point_dof_s = shape_function_accessor_reference_skeleton(f,space,measure)
###    prototype = GT.prototype(face_point_dof_s)
###    function face_point_dof_b(face,dummy=nothing)
###        face_around=faces_around[face]
###        face_point_dof_s(face,face_around)
###    end
###    accessor(face_point_dof_b,prototype)
###end
###
###function reference_map(refdface::AbstractFaceSpace,refDface::AbstractFaceSpace)
###    d = num_dims(refdface)
###    dof_to_f = shape_functions(refdface)
###    boundary = refDface |> GT.domain |> GT.mesh
###    lface_to_nodes = GT.face_nodes(boundary,d)
###    node_to_coords = GT.node_coordinates(boundary)
###    lface_to_lrefid = GT.face_reference_id(boundary,d)
###    lrefid_to_lrefface = GT.reference_spaces(boundary,d)
###    lrefid_to_perm_to_ids = map(GT.node_permutations,lrefid_to_lrefface)
###    
###    func_template = (dof_to_coeff, dof_to_f, ndofs) -> (x -> sum(dof->dof_to_coeff[dof]*dof_to_f[dof](x),1:ndofs))
###    d2c = node_to_coords[lface_to_nodes[1][ lrefid_to_perm_to_ids[lface_to_lrefid[1]][1] ]]
###    proto = func_template(d2c, dof_to_f, length(d2c))
###
###    result::Vector{Vector{typeof(proto)}} = map(1:GT.num_faces(boundary,d)) do lface
###        lrefid = lface_to_lrefid[lface]
###        nodes = lface_to_nodes[lface]
###        perm_to_ids = lrefid_to_perm_to_ids[lrefid]
###        map(perm_to_ids) do ids
###            #dof_to_coeff = node_to_coords[nodes[ids]] # Wrong
###            dof_to_coeff = similar(node_to_coords[nodes])
###            dof_to_coeff[ids] = node_to_coords[nodes]
###            ndofs = length(dof_to_coeff)
###            result_inner::typeof(proto) = func_template(dof_to_coeff, dof_to_f, ndofs)
###        end
###    end
###    return result
###end
###
###function inv_map(f,x0)
###    function pseudo_inverse_if_not_square(J)
###        m,n = size(J)
###        if m != n
###            pinv(J)
###        else
###            inv(J)
###        end
###    end
###    function invf(fx)
###        x = x0
###        tol = 1.0e-12
###        J = nothing
###        niters = 100
###        for _ in 1:niters
###            J = ForwardDiff.jacobian(f,x)
###            Jinv = pseudo_inverse_if_not_square(J)
###            dx = Jinv*(fx-f(x))
###            x += dx
###            if norm(dx) < tol
###                return x
###            end
###        end
###        error("Max iterations reached")
###        x
###    end
###end
###
###function nodes_accessor(mesh::AbstractMesh,vD,domain::AbstractDomain)
###    D = val_parameter(vD)
###    d = num_dims(domain)
###    faces_around = GT.faces_around(domain)
###    if d == D
###        nodes_accessor_interior(mesh,vD,domain)
###    elseif d+1==D && faces_around !== nothing
###        nodes_accessor_boundary(mesh,vD,domain)
###    else
###        nodes_accessor_skeleton(mesh,vD,domain)
###    end
###end
###
###function nodes_accessor_interior(mesh::AbstractMesh,vD,domain::AbstractDomain)
###    D = val_parameter(vD)
###    Dface_nodes = face_nodes(mesh,D)
###    face_to_Dface = faces(domain)
###    function face_to_nodes(face,face_around=nothing)
###        Dface = face_to_Dface[face]
###        Dface_nodes[Dface]
###    end
###    prototype = Int32[1]
###    accessor(face_to_nodes,prototype)
###end
###
###function nodes_accessor_skeleton(mesh::AbstractMesh,vD,domain::AbstractDomain)
###    D = val_parameter(vD)
###    d = num_dims(domain)
###    Dface_nodes = face_nodes(mesh,D)
###    face_to_dface = faces(domain)
###    topo = topology(mesh)
###    dface_to_Dfaces = face_incidence(topo,d,D)
###    function face_to_nodes(face,face_around)
###        dface = face_to_dface[face]
###        Dfaces = dface_to_Dfaces[dface]
###        Dface = Dfaces[face_around]
###        Dface_nodes[Dface]
###    end
###    prototype = Int32[1]
###    accessor(face_to_nodes,prototype)
###end
###
###function nodes_accessor_boundary(mesh::AbstractMesh,vD,domain::AbstractDomain)
###    face_to_nodes = nodes_accessor_skeleton(mesh,vD,domain)
###    faces_around = GT.faces_around(domain)
###    prototype = GT.prototype(face_to_nodes)
###    function face_to_nodes_2(face,dummy=nothing)
###        face_around=faces_around[face]
###        face_to_nodes(face,face_around)
###    end
###    accessor(face_to_nodes_2,prototype)
###end
###
#### Untabulated
###function physical_map_accessor(f,mesh::AbstractMesh,vD)
###    D = val_parameter(vD)
###    space = mesh_space(mesh,vD)
###    face_lnode_s = shape_function_accessor_reference(value,space)
###    node_to_x = node_coordinates(mesh)
###    Dface_nodes = face_nodes(mesh,D)
###    z = zero(eltype(node_to_x))
###    prototype = q->f(y->GT.prototype(face_lnode_s)(y)*z,q)
###    function face_phi(Dface,face_around=nothing)
###        lnode_to_s = face_lnode_s(Dface,face_around)
###        lnode_to_node = Dface_nodes[Dface]
###        function phi(q)
###            nlnodes = length(lnode_to_node)
###            sum(1:nlnodes) do lnode
###                node = lnode_to_node[lnode]
###                x = node_to_x[node]
###                s = lnode_to_s(lnode)
###                s(q)*x
###            end
###        end
###        q->f(phi,q)
###    end
###    accessor(face_phi,prototype)
###end
###
#### Tabulated
###
###function physical_map_accessor(f,measure::AbstractQuadrature,vD)
###    error("Case not implemented and not needed unless you want higher order derivatives of the physical map.")
###end
###
###function physical_map_accessor(f::typeof(value),measure::AbstractQuadrature,vD)
###    mesh = GT.mesh(domain(measure))
###    space = mesh_space(mesh,vD)
###    face_point_lnode_s = shape_function_accessor_reference(f,space,measure)
###    face_to_nodes = nodes_accessor(mesh,vD,domain(measure))
###    node_to_x = node_coordinates(mesh)
###    z = zero(eltype(node_to_x))
###    prototype = z*GT.prototype(face_point_lnode_s)
###    function face_point_phi(face,face_around=nothing)
###        point_lnode_s = face_point_lnode_s(face,face_around)
###        lnode_to_node = face_to_nodes(face,face_around)
###        function point_s(point)
###            lnode_s = point_lnode_s(point)
###            nlnodes = length(lnode_to_node)
###            sum(1:nlnodes) do lnode
###                node = lnode_to_node[lnode]
###                x = node_to_x[node]
###                s = lnode_s(lnode)
###                x*s
###            end
###        end
###    end
###    accessor(face_point_phi,prototype)
###end
###
###function physical_map_accessor(f::typeof(ForwardDiff.jacobian),measure::AbstractQuadrature,vD)
###    mesh = GT.mesh(domain(measure))
###    space = mesh_space(mesh,vD)
###    face_point_lnode_s = shape_function_accessor_reference(ForwardDiff.gradient,space,measure)
###    face_to_nodes = nodes_accessor(mesh,vD,domain(measure))
###    node_to_x = node_coordinates(mesh)
###    z = zero(eltype(node_to_x))
###    prototype = outer(z,GT.prototype(face_point_lnode_s))
###
###    D = val_parameter(vD)
###    P = typeof(prototype)
###
###    max_n_faces_around = 2
###    jacobian_cache = zeros(P, max_n_faces_around)
###    function face_point_phi(face,face_around=nothing)
###        point_lnode_s = face_point_lnode_s(face,face_around)
###        lnode_to_node = face_to_nodes(face,face_around)
###        index = face_around === nothing ? 1 : face_around
###        nlnodes = length(lnode_to_node)
###        is_simplex = (nlnodes == D + 1) # TODO: find a better way to check whether it is a simplex
###        if is_simplex
###            lnode_s_1 = point_lnode_s(1)
###
###            jacobian_cache[index] = sum(1:nlnodes) do lnode
###                node = lnode_to_node[lnode]
###                x = node_to_x[node]
###                s = lnode_s_1(lnode)
###                outer(x,s)
###            end
###        end
###        function point_s(point) 
###            if is_simplex
###                return jacobian_cache[index]
###            end
###            lnode_s = point_lnode_s(point)
###            # nlnodes = length(lnode_to_node)
###            jacobian_cache[index] = sum(1:nlnodes) do lnode
###                node = lnode_to_node[lnode]
###                x = node_to_x[node]
###                s = lnode_s(lnode)
###                outer(x,s)
###            end
###            return jacobian_cache[index]
###        end
###    end
###    accessor(face_point_phi,prototype)
###end
###
###function coordinate_accessor(measure::AbstractQuadrature)
###    vD=Val(num_dims(domain(measure)))
###    physical_map_accessor(value,measure,vD)
###end
###
###function jacobian_accessor(measure::AbstractQuadrature,vD=Val(num_dims(domain(measure))))
###    physical_map_accessor(ForwardDiff.jacobian,measure,vD)
###end
###
###function num_points_accessor(measure::AbstractQuadrature)
###    rid_to_point_to_w = map(weights,reference_quadratures(measure))
###    rid_to_n = map(length,rid_to_point_to_w)
###    dface_to_rid = face_reference_id(measure)
###    face_to_dface = faces(domain(measure))
###    prototype = 1
###    function face_npoints(face,face_around=nothing)
###        dface = face_to_dface[face]
###        rid = dface_to_rid[dface]
###        rid_to_n[rid]
###    end
###    accessor(face_npoints,prototype)
###end
###
###function weight_accessor(measure::AbstractQuadrature)
###    if is_reference_domain(domain(measure))
###        weight_accessor_reference(measure)
###    else
###        weight_accessor_physical(measure)
###    end
###end
###
###function weight_accessor_reference(measure::AbstractQuadrature)
###    face_point_J = jacobian_accessor(measure)
###    rid_to_point_to_w = map(weights,reference_quadratures(measure))
###    dface_to_rid = face_reference_id(measure)
###    face_to_dface = faces(domain(measure))
###    prototype = first(first(rid_to_point_to_w))
###    function face_point_v(face,face_around=nothing)
###        dface = face_to_dface[face]
###        rid = dface_to_rid[dface]
###        point_to_w = rid_to_point_to_w[rid]
###        function point_v(point, J=nothing)
###            point_to_w[point]
###        end
###    end
###    accessor(face_point_v,prototype)
###end
###
###function weight_accessor_physical(measure::AbstractQuadrature)
###    face_point_v = weight_accessor_reference(measure)
###    face_point_J = jacobian_accessor(measure)
###    prototype = change_of_measure(GT.prototype(face_point_J))*GT.prototype(face_point_v)
###    function face_point_w(face,face_around=nothing)
###        point_v = face_point_v(face)
###        point_J = face_point_J(face)
###        function point_w(point,J = point_J(point))
###            v = point_v(point)
###            change_of_measure(J)*v
###        end
###        return point_w
###    end
###    accessor(face_point_w,prototype)
###end
###
#### push reference shape functions to physical ones
###
###function shape_function_accessor_modifier(f::typeof(value),space::AbstractSpace)
###    function modifier(v,J)
###        v
###    end
###    function face_dof_modifier(face,face_around=nothing)
###        function dof_modifier(dof)
###            return modifier
###        end
###    end
###    accessor(face_dof_modifier,modifier)
###end
###
###function shape_function_accessor_modifier(f::typeof(ForwardDiff.gradient),space::AbstractSpace)
###    function modifier(v,J)
###        transpose(J)\v 
###    end
###    function face_dof_modifier(face,face_around=nothing)
###        function dof_modifier(dof)
###            return modifier
###        end
###    end
###    accessor(face_dof_modifier,modifier)
###end
###
###function shape_function_accessor_modifier(f::typeof(ForwardDiff.jacobian),space::AbstractSpace)
###    function modifier(v,J)
###        v/J
###    end
###    function face_dof_modifier(face,face_around=nothing)
###        function dof_modifier(dof)
###            return modifier
###        end
###    end
###    accessor(face_dof_modifier,modifier)
###end
###
###
###function shape_function_accessor_modifier(f,space::AbstractSpace,domain::AbstractDomain)
###    D = num_dims(space)
###    d = num_dims(domain)
###    faces_around = GT.faces_around(domain)
###    if d == D
###        shape_function_accessor_modifier_interior(f,space,domain)
###    elseif d+1==D && faces_around !== nothing
###        shape_function_accessor_modifier_boundary(f,space,domain)
###    else
###        shape_function_accessor_modifier_skeleton(f,space,domain)
###    end
###end
###
###function shape_function_accessor_modifier_interior(f,space::AbstractSpace,domain::AbstractDomain)
###    Dface_to_modif = shape_function_accessor_modifier(f,space)
###    face_to_Dface = faces(domain)
###    function face_to_modif(face,face_around=nothing)
###        Dface = face_to_Dface[face]
###        Dface_to_modif(Dface)
###    end
###    accessor(face_to_modif,GT.prototype(Dface_to_modif))
###end
###
###function shape_function_accessor_modifier_skeleton(f,space::AbstractSpace,domain::AbstractDomain)
###    D = num_dims(GT.domain(space))
###    d = num_dims(domain)
###    face_to_dface = faces(domain)
###    topo = topology(GT.mesh(space))
###    dface_to_Dfaces = face_incidence(topo,d,D)
###    Dface_to_modif = shape_function_accessor_modifier(f,space)
###    function face_to_modif(face,face_around)
###        dface = face_to_dface[face]
###        Dfaces = dface_to_Dfaces[dface]
###        Dface = Dfaces[face_around]
###        Dface_to_modif(Dface)
###    end
###    accessor(face_to_modif,GT.prototype(Dface_to_modif))
###end
###
###function shape_function_accessor_modifier_boundary(f,space::AbstractSpace,domain::AbstractDomain)
###    face_modif = shape_function_accessor_modifier_skeleton(f,space,domain)
###    faces_around = GT.faces_around(domain)
###    function face_to_modif(face,dummy=nothing)
###        face_around=faces_around[face]
###        face_modif(face,face_around)
###    end
###    accessor(face_to_modif,GT.prototype(face_modif))
###end
###
#### Untabulated
###function shape_function_accessor_physical(f,space::AbstractSpace)
###    domain = GT.domain(space)
###    D = num_dims(domain)
###    mesh = GT.mesh(domain)
###    face_dof_s_ref = shape_function_accessor_reference(value,space)
###    face_dof_modif = shape_function_accessor_modifier(value,space)
###    face_phi = physical_map_accessor(value,mesh,Val(D))
###    face_Dphi = physical_map_accessor(ForwardDiff.jacobian,mesh,Val(D))
###    x0 = zero(SVector{D,Float64})
###    prototype = q->f(x->GT.prototype(face_dof_modif)(GT.prototype(face_dof_s_ref)(x),GT.prototype(face_Dphi)(x)),q)
###    function face_dof_s_phys(face,face_around=nothing)
###        phi = face_phi(face)
###        Dphi = face_Dphi(face)
###        invphi = inv_map(phi,x0)
###        dof_s_ref = face_dof_s_ref(face)
###        dof_modif = face_dof_modif(face)
###        function dof_s_phys(dof)
###            modif = dof_modif(dof)
###            s_ref = dof_s_ref(dof)
###            function s_phys(x)
###                v = s_ref(x)
###                J = Dphi(x)
###                modif(v,J)
###            end
###            x->f(s_phys,x)
###        end
###    end
###    accessor(face_dof_s_phys,prototype)
###end
###
#### This gets tabulated only for particular instances of f.
###function shape_function_accessor_physical(f,space::AbstractSpace,measure::AbstractQuadrature,disable_tabulation)
###    if num_dims(domain(space)) != num_dims(domain(measure))
###        error("case not implemented")
###    end
###    Dface_dof_x_s = shape_function_accessor_physical(f,space)
###    face_point_x = coordinate_accessor(measure)
###    face_dface = faces(domain(measure))
###    function face_point_dof_s(face,face_around=nothing)
###        point_x = face_point_x(face)
###        Dface = face_to_dface[face]
###        dof_x_s = Dface_dof_x_s(Dface)
###        function point_dof_s(point,J=nothing)
###            x = point_x(point)
###            function dof_s(dof)
###                dof_x_s(dof)(x)
###            end
###        end
###    end
###    prototype = GT.prototype(Dface_dof_x_s)(GT.prototype(face_point_x))
###    accessor(face_point_dof_s,prototype)
###end
###
#### Tabulated
###for T in (:value,:(ForwardDiff.gradient),:(ForwardDiff.jacobian))
###    @eval begin
###        function shape_function_accessor_physical(f::typeof($T),space::AbstractSpace,measure::AbstractQuadrature,disable_tabulation::Val{false})
###            face_point_dof_v = shape_function_accessor_reference(f,space,measure)
###            face_ndofs = num_dofs_accessor(space,GT.domain(measure))
###            dface_to_modif = shape_function_accessor_modifier(f,space,domain(measure))
###            D = num_dims(domain(space))
###            face_point_Dphi = jacobian_accessor(measure,Val(D))
###            prototype = GT.prototype(dface_to_modif)(GT.prototype(face_point_dof_v),GT.prototype(face_point_Dphi))
###            P = typeof(prototype)
###            max_n_faces_around = 2
###
###            # face_around_dof_s = fill(zeros(P,max_num_reference_dofs(space)),max_n_faces_around) # incorrect, all elements have the same objectid
###            face_around_dof_s::Vector{Vector{P}} = [zeros(P,max_num_reference_dofs(space)) for _ in 1:max_n_faces_around]
###
###            function face_point_dof_s(face,face_around=nothing)
###                point_dof_v = face_point_dof_v(face,face_around)
###                dof_modif = dface_to_modif(face,face_around)
###                ndofs = face_ndofs(face,face_around)
###                point_Dphi = face_point_Dphi(face,face_around)
###
###                dof_s = if face_around === nothing # define dof_s outside the branches, so that the compiler can inference the type easier
###                    face_around_dof_s[1] 
###                else
###                    face_around_dof_s[face_around]
###                end
###                
###                function point_J_dof_s(point,J)
###                    dof_v = point_dof_v(point)
###                    for dof in 1:ndofs
###                       v = dof_v(dof)
###                       modif = dof_modif(dof)
###                       dof_s[dof] = modif(v,J)
###                    end
###                    function dof_f(dof)
###                        dof_s[dof]
###                    end
###                end
###                function point_dof_s(point,J = nothing)
###                    J2 = J === nothing ? point_Dphi(point) : J
###                    point_J_dof_s(point,J2)
###                end
###                return point_dof_s
###            end
###            accessor(face_point_dof_s,prototype)
###        end # function
###
###
###        function shape_function_accessor_physical(f::typeof($T),space::AbstractSpace,measure::AbstractQuadrature,disable_tabulation::Val{true})
###            face_point_dof_v = shape_function_accessor_reference(f,space,measure)
###            face_ndofs = num_dofs_accessor(space,GT.domain(measure))
###            dface_to_modif = shape_function_accessor_modifier(f,space,domain(measure))
###            D = num_dims(domain(space))
###            face_point_Dphi = jacobian_accessor(measure,Val(D))
###            prototype = GT.prototype(dface_to_modif)(GT.prototype(face_point_dof_v),GT.prototype(face_point_Dphi))
###            
###            function face_point_dof_s(face,face_around=nothing)
###                point_dof_v = face_point_dof_v(face,face_around)
###                dof_modif = dface_to_modif(face,face_around)
###                point_Dphi = face_point_Dphi(face,face_around)
###                
###                function point_J_dof_s(point,J)
###                    dof_v = point_dof_v(point)
###                    function dof_f(dof)
###                        v = dof_v(dof)
###                        modif = dof_modif(dof)
###                        modif(v,J)
###                    end
###                end
###                function point_dof_s(point,J = nothing)
###                    J2 = J === nothing ? point_Dphi(point) : J
###                    point_J_dof_s(point,J2)
###                end
###                return point_dof_s
###            end
###            accessor(face_point_dof_s,prototype)
###        end # function
###    end # @eval
###end # for
###
###function dofs_accessor(space::AbstractSpace,domain::AbstractDomain)
###    D = num_dims(space)
###    d = num_dims(domain)
###    faces_around = GT.faces_around(domain)
###    if d == D
###        dofs_accessor_interior(space,domain)
###    elseif d+1==D && faces_around !== nothing
###        dofs_accessor_boundary(space,domain)
###    else
###        dofs_accessor_skeleton(space,domain)
###    end
###end
###
###function dofs_accessor_interior(space::AbstractSpace,dom::AbstractDomain)
###    face_to_Dface = faces(dom)
###    Dface_to_dofs_ = face_dofs(space)
###    prototype = Int32[1]
###    function face_to_dofs(face,face_around=nothing)
###        Dface = face_to_Dface[face]
###        dofs = Dface_to_dofs_[Dface]
###        dofs
###    end
###    prototype = Int32[1]
###    accessor(face_to_dofs,prototype)
###end
###
###function dofs_accessor_skeleton(space::AbstractSpace,domain::AbstractDomain)
###    D = num_dims(GT.domain(space))
###    d = num_dims(domain)
###    Dface_dofs = face_dofs(space)
###    face_to_dface = faces(domain)
###    topo = topology(GT.mesh(space))
###    dface_to_Dfaces = face_incidence(topo,d,D)
###    function face_to_dofs(face,face_around)
###        dface = face_to_dface[face]
###        Dfaces = dface_to_Dfaces[dface]
###        Dface = Dfaces[face_around]
###        Dface_dofs[Dface]
###    end
###    prototype = Int32[1]
###    accessor(face_to_dofs,prototype)
###end
###
###function dofs_accessor_boundary(space::AbstractSpace,domain::AbstractDomain)
###    face_to_dofs = dofs_accessor_skeleton(space,domain)
###    faces_around = GT.faces_around(domain)
###    prototype = GT.prototype(face_to_dofs)
###    function face_to_dofs_2(face,dummy=nothing)
###        face_around=faces_around[face]
###        face_to_dofs(face,face_around)
###    end
###    accessor(face_to_dofs_2,prototype)
###end
###
###function num_dofs_accessor(space::AbstractSpace,domain::AbstractDomain)
###    D = num_dims(space)
###    d = num_dims(domain)
###    faces_around = GT.faces_around(domain)
###    if d == D
###        num_dofs_accessor_interior(space,domain)
###    elseif d+1==D && faces_around !== nothing
###        num_dofs_accessor_boundary(space,domain)
###    else
###        num_dofs_accessor_skeleton(space,domain)
###    end
###end
###
###function num_dofs_accessor_interior(space::AbstractSpace,domain::AbstractDomain)
###    rid_to_n = map(num_dofs,reference_spaces(space))
###    Dface_to_rid = face_reference_id(space)
###    face_to_Dface = faces(domain)
###    function face_to_n(face,face_around=nothing)
###        Dface = face_to_Dface[face]
###        rid = Dface_to_rid[Dface]
###        rid_to_n[rid]
###    end
###    accessor(face_to_n,1)
###end
###
###function num_dofs_accessor_skeleton(space::AbstractSpace,domain::AbstractDomain)
###    mesh = GT.mesh(space)
###    rid_to_n = map(num_dofs,reference_spaces(space))
###    Dface_to_rid = face_reference_id(space)
###    face_to_dface = faces(domain)
###    topo = topology(mesh)
###    d = num_dims(domain)
###    D = num_dims(GT.domain(space))
###    dface_to_Dfaces = face_incidence(topo,d,D)
###    function face_to_n(face,face_around)
###        dface = face_to_dface[face]
###        Dfaces = dface_to_Dfaces[dface]
###        Dface = Dfaces[face_around]
###        rid = Dface_to_rid[Dface]
###        rid_to_n[rid]
###    end
###    accessor(face_to_n,1)
###end
###
###function num_dofs_accessor_boundary(space::AbstractSpace,domain::AbstractDomain)
###    a = num_dofs_accessor_skeleton(space,domain)
###    faces_around = GT.faces_around(domain)
###    function face_to_n(face,dummy=nothing)
###        face_around=faces_around[face]
###        a(face,face_around)
###    end
###    accessor(face_to_n,GT.prototype(a))
###end
###
###struct DiscreteFieldAccessor{A,B} <: AbstractAccessor
###    update::A
###    accessor::B
###end
###
###prototype(f::DiscreteFieldAccessor) = prototype(f.accessor)
###
###function (f::DiscreteFieldAccessor)(face,face_around=nothing)
###    f.accessor(face,face_around)
###end
###
###function update(f::DiscreteFieldAccessor;discrete_field)
###    uh = discrete_field
###    accessor = f.update(uh)
###    DiscreteFieldAccessor(f.update,accessor)
###end
###
###function discrete_field_accessor(f,uh::DiscreteField,measure::AbstractQuadrature)
###    face_to_dofs = dofs_accessor(GT.space(uh),GT.domain(measure))
###    face_to_point_to_ldof_to_s = shape_function_accessor(f,GT.space(uh),measure)
###    function field_to_accessor(uh)
###        space = GT.space(uh)
###        free_values = GT.free_values(uh)
###        dirichlet_values = GT.dirichlet_values(uh)
###        prototype = zero(eltype(free_values))*GT.prototype(face_to_point_to_ldof_to_s)
###        z = zero(prototype)
###        function face_point_u(face,face_around=nothing)
###            ldof_to_dof = face_to_dofs(face,face_around)
###            point_to_ldof_to_s = face_to_point_to_ldof_to_s(face,face_around)
###            function point_u(point,J=nothing)
###                ldof_to_s = point_to_ldof_to_s(point,J)
###                nldofs = length(ldof_to_dof)
###                sum(1:nldofs;init=z) do ldof
###                    dof = ldof_to_dof[ldof]
###                    s = ldof_to_s(ldof)
###                    if dof > 0
###                        v = free_values[dof]
###                    else
###                        v = dirichlet_values[-dof]
###                    end
###                    v*s
###                end
###            end
###        end
###        GT.accessor(face_point_u,prototype)
###    end
###    accessor = field_to_accessor(uh)
###    DiscreteFieldAccessor(field_to_accessor,accessor)
###end
###
##### TODO not needed anymore
####struct DirichletAccessor{A,B} <: AbstractAccessor
####    update::A
####    accessor::B
####end
####
####prototype(f::DirichletAccessor) = prototype(f.accessor)
####
####function (f::DirichletAccessor)(face)
####    f.accessor(face)
####end
####
####function update(f::DirichletAccessor;discrete_field)
####    uh = discrete_field
####    accessor = f.update(uh)
####    DirichletAccessor(f.update,accessor)
####end
####
####function dirichlet_accessor(uh::DiscreteField,domain::AbstractDomain)
####    face_to_dofs = dofs_accessor(GT.space(uh),domain)
####    function field_to_accessor(uh)
####        space = GT.space(uh)
####        dirichlet_values = GT.dirichlet_values(uh)
####        prototype = nothing
####        function face_dirichlet!(face,face_around=nothing)
####            ldof_to_dof = face_to_dofs(face,face_around)
####            nldofs = length(ldof_to_dof)
####            function dirichlet!(A,b)
####                m = size(A,1)
####                z = zero(eltype(b))
####                for i in 1:m
####                    bi = z
####                    for j in 1:nldofs
####                        dof = ldof_to_dof[j]
####                        if dof < 0
####                            uj = dirichlet_values[-dof]
####                            bi += A[i,j]*uj
####                        end
####                    end
####                    b[i] -= bi
####                end
####                nothing
####            end
####        end
####        GT.accessor(face_dirichlet!,prototype)
####    end
####    accessor = field_to_accessor(uh)
####    DirichletAccessor(field_to_accessor,accessor)
####end
###
###function unit_normal_accessor(measure::AbstractQuadrature)
###    domain = GT.domain(measure)
###    mesh = GT.mesh(domain)
###    D = num_dims(mesh)
###    d = num_dims(domain)
###    if D==d && is_physical_domain(domain)
###        unit_normal_accessor_physical_interior(measure)
###    elseif D==d+1 && is_physical_domain(domain)
###        unit_normal_accessor_physical(measure)
###    elseif D==d+1 && is_reference_domain(domain)
###        unit_normal_accessor_reference(measure)
###    else
###        error()
###    end
###end
###
###function unit_normal_accessor_physical_interior(measure::AbstractQuadrature)
###    error("not implemented")
###end
###
###function unit_normal_accessor_physical(measure::AbstractQuadrature)
###    domain = GT.domain(measure)
###    d = num_dims(domain)
###    D = d+1
###    face_n_ref = unit_normal_accessor_reference(measure)
###    face_point_J = jacobian_accessor(measure,Val(D))
###    function face_point_n(face,face_around=nothing)
###        point_J = face_point_J(face,face_around)
###        n_ref = face_n_ref(face,face_around)
###        function n_phys(point, J=nothing)
###            J2 = J === nothing ? point_J(point) : J
###            map_unit_normal(J2,n_ref)
###        end
###        return n_phys
###    end
###    accessor(face_point_n,GT.prototype(face_n_ref))
###end
###
###
###function unit_normal_accessor_reference(measure::AbstractQuadrature)
###    D = num_dims(mesh(measure))
###    d = num_dims(domain(measure))
###    @assert D == d+1
###    if GT.faces_around(domain(measure)) === nothing
###        unit_normal_accessor_reference_skeleton(measure)
###    else
###        unit_normal_accessor_reference_boundary(measure)
###    end
###end
###
###function unit_normal_accessor_reference_skeleton(measure::AbstractQuadrature)
###    domain = GT.domain(measure)
###    d = num_dims(domain)
###    D = d+1
###    mesh = GT.mesh(domain)
###    topo = topology(mesh)
###    dface_to_Dfaces, dface_to_ldfaces = GT.face_incidence_ext(topo,d,D)
###    Dface_to_Drid = face_reference_id(mesh,D)
###    face_to_dface = faces(domain)
###    Drid_to_ldface_to_n = map(GT.reference_spaces(mesh,Val(D))) do refface
###        boundary = refface |> GT.domain |> GT.mesh
###        boundary |> GT.normals # TODO also rename?
###    end
###    prototype = first(first(Drid_to_ldface_to_n))
###    function face_n(face,face_around)
###        dface = face_to_dface[face]
###        Dfaces = dface_to_Dfaces[dface]
###        ldfaces = dface_to_ldfaces[dface]
###        Dface = Dfaces[face_around]
###        ldface = ldfaces[face_around]
###        Drid = Dface_to_Drid[Dface]
###        n = Drid_to_ldface_to_n[Drid][ldface]
###    end
###    accessor(face_n,prototype)
###end
###
###function unit_normal_accessor_reference_boundary(measure::AbstractQuadrature)
###    face_n = unit_normal_accessor_reference_skeleton(measure)
###    faces_around = GT.faces_around(domain(measure))
###    function face_n_2(face,dummy=nothing)
###        face_around=faces_around[face]
###        face_n(face,face_around)
###    end
###    accessor(face_n_2,GT.prototype(face_n))
###end
###
###
###function face_diameter(domain::AbstractDomain)
###    d = num_dims(domain)
###    dinv = 1/d
###    measure = GT.measure(domain,1)
###    face_point_w = weight_accessor(measure)
###    face_npoints = num_points_accessor(measure)
###    z = zero(prototype(face_point_w))
###    nfaces = num_faces(domain)
###    diams = fill(z,nfaces)
###    for face in 1:nfaces
###        point_w = face_point_w(face)
###        npoints = face_npoints(face)
###        s = z
###        for point in 1:npoints
###            w = point_w(point)
###            s += w
###        end
###        diams[face] =  s^dinv
###    end
###    diams
###end
###
###
####function form_argument_accessor(f,space::AbstractSpace)
####    shape_function_accessor(f,space)
####end
####
#### TODO: type instability when using accessor for customized functions
###function form_argument_accessor(f,space::AbstractSpace,measure::AbstractQuadrature,field=1)
###    face_point_dof_s = shape_function_accessor(f,space,measure)
###    prototype = GT.prototype(face_point_dof_s)
###    the_field = field
###    z = zero(prototype)
###    function face_point_dof_a(face,face_around=nothing)
###        the_face_around = face_around
###        point_dof_s = face_point_dof_s(face,face_around)
###        function point_dof_a(point, J=nothing)
###            dof_s = point_dof_s(point, J)
###            function dof_a(dof,field=1,face_around=nothing)
###                mask = face_around == the_face_around && field == the_field
###                if mask
###                    s = dof_s(dof)
###                    s
###                else
###                    z
###                end
###            end
###        end
###    end
###    accessor(face_point_dof_a,prototype)
###end
###
###function num_faces_around_accesor(domain_space,domain)
###    d = num_dims(domain)
###    D = num_dims(domain_space)
###    faces_around = GT.faces_around(domain)
###    if d == D
###        num_faces_around_accesor_interior(domain_space,space)
###    elseif d+1==D && faces_around !== nothing
###        num_faces_around_accesor_interior(domain_space,space)
###    else
###        num_faces_around_accesor_skeleton(domain_space,space)
###    end
###end
###
###function num_faces_around_accesor_interior(space_domain,domain)
###    function n_faces_around(face,face_around=nothing)
###        1
###    end
###    accessor(n_faces_around,1)
###end
###
###function num_faces_around_accesor_skeleton(space_domain,domain)
###    function n_faces_around(face,face_around=nothing)
###        2
###    end
###    accessor(n_faces_around,1)
###end
###
###
###function max_num_points(measure::AbstractQuadrature)
###    lengths = map(x -> length(weights(x)), reference_quadratures(measure))
###    return reduce(max, lengths)
###end
###
#### remove it as we use quadrature directly
#### function num_points_accessor(measure::Measure)
####     num_points_accessor(quadrature(measure))
#### end
###
#### function coordinate_accessor(measure::Measure)
####     coordinate_accessor(quadrature(measure))
#### end
###
#### function jacobian_accessor(measure::Measure,args...)
####     jacobian_accessor(quadrature(measure),args...)
#### end
###
#### function weight_accessor(measure::Measure)
####     weight_accessor(quadrature(measure))
#### end
###
#### function discrete_field_accessor(f,uh::DiscreteField,measure::Measure)
####     discrete_field_accessor(f,uh,quadrature(measure))
#### end
###
#### function shape_function_accessor(f,space::AbstractSpace,measure::Measure)
####     shape_function_accessor(f,space,quadrature(measure))
#### end
###
#### function form_argument_accessor(f,space::AbstractSpace,measure::Measure)
####     form_argument_accessor(f,space,quadrature(measure))
#### end
###
#### function physical_map_accessor(f,measure::Measure,vD)
####     physical_map_accessor(f,quadrature(measure),vD)
#### end
###
#### function unit_normal_accessor(measure::Measure)
####     unit_normal_accessor(quadrature(measure))
#### end
###
###
#### currently we do 1 step lowering
#### shape_function_accessor_modifier_term: always tabulated
#### function shape_function_accessor_modifier_term(f, v, J)
####     if f == value
####         v
####     elseif f == ForwardDiff.jacobian 
####         :($v/$J)
####     elseif f == ForwardDiff.gradient
####         :(transpose($J)\$v)
####     else
####         error("shape function accessor modifier not supported for this function f: $f")
####     end
#### end
###
###
#### function shape_function_accessor_modifier_interior_term(f, v, J)
####     shape_function_accessor_modifier_term(f, v, J)
#### end
###
#### function shape_function_accessor_modifier_skeleton_term(f, space, dom, face, the_face_around, dof, v, J)
####     # is it correct?
####     shape_function_accessor_modifier_term(f, v, J)
#### end
###
#### function shape_function_accessor_modifier_boundary_term(f, space, dom, face, the_face_around, dof, v, J)
####     face_around = :($(GT.face_around)(dom))
####     shape_function_accessor_modifier_skeleton_term(f, space, dom, face, face_around, dof, v, J)
#### end
###
###
#### function shape_function_accessor_modifier_term(f, space, dom, face, the_face_around, dof, v, J, integral_type)
####     # TODO: inline 1 step further
####     if integral_type == :interior 
####         shape_function_accessor_modifier_interior_term(f, v, J)
####     elseif integral_type == :boundary 
####         shape_function_accessor_modifier_boundary_term(f, space, dom, face, the_face_around, dof, v, J)
####     else
####         shape_function_accessor_modifier_skeleton_term(f, space, dom, face, the_face_around, dof, v, J)
####     end
#### end
###
###
###
#### function shape_function_accessor_physical_term(f, space, measure, face, the_face_around, point, J, dof, integral_type)
####     face_point_dof_v = :(shape_function_accessor_reference($f,$space,$measure)) # TODO: further expand it
###
####     # dface_to_modif = :(shape_function_accessor_modifier($f,$space, $(GT.domain)($measure)))
####     D = :(num_dims($(GT.domain)($space)))
####     face_point_Dphi = :(jacobian_accessor($measure,Val($D)))
###
####     point_dof_v = :($face_point_dof_v($face,$the_face_around))
####     # dof_modif = :($dface_to_modif($face, $the_face_around))
####     point_Dphi = :($face_point_Dphi($face,$the_face_around))
###    
####     # J2 = :(ifelse($J == nothing, point_Dphi($point), $J))
####     # TODO: can we assume that the jacobian is always passed from outside?
####     J2 = if J === nothing
####         :($point_Dphi($point))
####     else
####         J
####     end
####     dof_v = :($point_dof_v($point))
####     v = :($dof_v($dof))
####     # modif = :($dof_modif($dof))
###
####     shape_function = shape_function_accessor_modifier_term(f, space, :($(GT.domain)($measure)), face, the_face_around, dof, v, J2, integral_type)
###
####     shape_function
####     # z, z
#### end # function
###
###
#### function shape_function_accessor_term(f, space, measure, face, the_face_around, point, J, dof, is_reference, integral_type)
####     tabulated_functions = Set([value, ForwardDiff.jacobian, ForwardDiff.gradient])
####     if is_reference
####         # z = :( zero($(GT.prototype)(shape_function_accessor_reference($f,$space,$measure))) )
####         shape_function = :(shape_function_accessor_reference($f,$space,$measure)($face, $the_face_around)($point, $J)($dof))
####     elseif f in tabulated_functions
####         shape_function = shape_function_accessor_physical_term(f, space, measure, face, the_face_around, point, J, dof, integral_type)
####         # z = shape_function_accessor_physical_term(f, space, measure, 1, 1, 1, nothing, 1, integral_type) # TODO: find a better way to get prototype. However, this will probably be removed in the final code.
####     else # Untabulated, keep using accessors
####         # z = :( zero($(GT.prototype)(shape_function_accessor_physical($f,$space,$measure))) )
####         shape_function = :(shape_function_accessor_physical($f,$space,$measure)($face, $the_face_around)($point, $J)($dof))
####     end
####     return shape_function
#### end
###
###
###function form_argument_accessor_term(f,space,measure,the_field, face,the_face_around, point, J, dof,field,face_around, is_reference, integral_type)
###
###    mask = :(($face_around == $the_face_around && $field == $the_field))
###    disable_tabulation = :(Val(true))
###    # shape_function = shape_function_accessor_term(f, space, measure, face, the_face_around, point, J, dof, is_reference, integral_type) # TODO: add an arg to remove tabulation
###    z = :( zero(GT.prototype(shape_function_accessor($f,$space,$measure,$disable_tabulation))) ) # TODO find a better way to do the prototype, maybe inlining 1 step further
###    shape_function = :(shape_function_accessor($f,$space,$measure,$disable_tabulation)($face, $the_face_around)($point, $J)($dof))
###    :(ifelse($mask, $shape_function, $z))
###
###end
###

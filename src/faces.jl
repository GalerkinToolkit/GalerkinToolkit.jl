
abstract type AbstractFaceNew <: AbstractType end
abstract type AbstractMeshFaceNew <: AbstractFaceNew end
abstract type AbstractQuadratureFaceNew <: AbstractFaceNew end
abstract type AbstractSpaceFaceNew <: AbstractFaceNew end
abstract type AbstractDiscreteFieldFaceNew <: AbstractSpaceFaceNew end
abstract type AbstractSamplingFaceNew <: AbstractFaceNew end

abstract type AbstractPointNew <: AbstractType end
#abstract type AbstractQuadraturePointNew <: AbstractPointNew end
#abstract type AbstractSamplingPointNew <: AbstractPointNew end

struct Each{A,B,C}
    update_id::A
    prototype::B
    ids::C
end
Base.length(iter::Each) = length(iter.ids)
Base.isdone(iter::Each,id) = id > length(iter)
function Base.getindex(iter::Each,i::Integer)
    id = iter.ids[i]
    iter.update_id(iter.prototype,id)
end
function Base.iterate(iter::Each,id=1)
    if Base.isdone(iter,id)
        nothing
    else
        accessor = iter[id]
        (accessor,id+1)
    end
end
Base.keys(iter::Each) = LinearIndices((length(iter),))
prototype(x::Each) = x.prototype
function tabulate(f,x::Each)
    Each(x.update_id,tabulate(f,x.prototype),x.ids)
end

struct AnyIndexNew <: AbstractType end
const ANY_INDEX_NEW = AnyIndexNew()

struct MeshFaceNew{A,B,C} <: AbstractMeshFaceNew
    mesh::A
    num_dims::B
    id::C
end
mesh(f::MeshFaceNew) = f.mesh
id(f::MeshFaceNew) = f.id
num_dims(f::MeshFaceNew) = val_parameter(f.num_dims)
function update_id(face::MeshFaceNew,id)
    MeshFaceNew(face.mesh,face.num_dims,id)
end
function update_id(face::MeshFaceNew,face2::MeshFaceNew)
    @assert num_dims(face) == num_dims(face_2)
    id = GT.id(face2)
    MeshFaceNew(face.mesh,face.num_dims,id)
end

function each_face_new(domain::AbstractDomain)
    mesh = GT.mesh(domain)
    d = num_dims(domain)
    num_dims = Val(d)
    id = ANY_INDEX_NEW
    face = MeshFaceNew(mesh,num_dims,id)
    ids = faces(domain)
    Each(update_id,face,ids)
end

function each_face_new(mesh::AbstractMesh,vald)
    d = val_parameter(vald)
    num_dims = Val(d)
    ids = GT.faces(mesh,d)
    face = MeshFaceNew(mesh,num_dims,id)
    Each(update_id,face,ids)
end

function num_faces_around(dface::AbstractMeshFaceNew,valD)
    length(each_face_around_new(dface,valD))
end

function each_face_around_new(dface::AbstractMeshFaceNew,valD)
    D = val_parameter(valD)
    d = GT.num_dims(dface)
    mesh = GT.mesh(dface)
    topo = topology(mesh)
    ids = face_incidence(topo,d,D)[GT.id(dface)]
    num_dims = Val(D)
    id = ANY_INDEX_NEW
    Dface = MeshFaceNew(mesh,num_dims,id)
    Each(update_id,Dface,ids)
end

function each_local_face(Df::MeshFaceNew,vald)
    d = val_parameter(vald)
    D = GT.num_dims(Df)
    l_mesh = GT.mesh(GT.reference_domain(Df))
    l_id = ANY_INDEX_NEW
    l_dface = MeshFaceNew(l_mesh,Val(d),l_id)
    nldfaces = num_faces(l_mesh,d)
    Each(update_id,l_dface,1:nldfaces)
end

function local_face(Df::AbstractMeshFaceNew,df::AbstractMeshFaceNew)
    D = GT.num_dims(Df)
    d = GT.num_dims(df)
    L_faces = each_local_face(Df,Val(d))
    mesh = GT.mesh(Df)
    topo = topology(mesh)
    Df_id = GT.id(Df)
    df_id = GT.id(df)
    Bs = face_incidence(topo,D,d)[Df_id]
    a_id = findfirst(B->B==df_id,Bs)
    L_faces[a_id]
end

function num_local_faces(Df::MeshFaceNew,vald)
    length(each_local_face(Df,vald))
end

function global_face(A::AbstractFaceNew,a::AbstractFaceNew)
    D = GT.num_dims(A)
    d = GT.num_dims(a)
    mesh = GT.mesh(A)
    topo = topology(mesh)
    id = face_incidence(topo,D,d)[GT.id(A)][GT.id(a)]
    MeshFaceNew(mesh,Val(d),id)
end

struct QuadratureFaceNew{A,B} <: AbstractQuadratureFaceNew
    mesh_quadrature::A
    id::B
end
mesh_quadrature(a::QuadratureFaceNew) = a.mesh_quadrature
mesh(a::QuadratureFaceNew) = mesh(domain(mesh_quadrature(a)))
function mesh_face(a::QuadratureFaceNew)
    mesh = GT.mesh(a.mesh_quadrature)
    d = num_dims(GT.domain(a.mesh_quadrature))
    MeshFaceNew(mesh,Val(d),a.id)
end
num_dims(a::QuadratureFaceNew) = num_dims(GT.domain(a.mesh_quadrature))
id(a::QuadratureFaceNew) = a.id
update_id(a::QuadratureFaceNew,id) = QuadratureFaceNew(mesh_quadrature(a),id)
function update_id(a::QuadratureFaceNew,face2::MeshFaceNew)
    q = mesh_quadrature(a)
    @assert num_dims(q) == num_dims(face_2)
    QuadratureFaceNew(q,GT.id(face2))
end


struct PointNew{A,B} <: AbstractPointNew
    parent_face::A
    id::B
end
parent_face(a::PointNew) = a.parent_face
id(a::PointNew) = a.id
update_id(a::PointNew,id) = PointNew(parent_face(a),id)
update_id(a::PointNew,point::AbstractPointNew) = PointNew(parent_face(a),id(point))

function each_point_new(face::AbstractFaceNew)
    id = ANY_INDEX_NEW
    point = PointNew(face,id)
    npoints = num_points(face)
    Each(update_id,point,1:npoints)
end

struct MeshSamplingFaceNew{A,B,C,D,E,F} <: AbstractSamplingFaceNew
    mesh_face::A
    quadrature_face::B
    local_face::C
    tabulated_values::D
    tabulated_gradients::E
    reference_unit_normals::F
end
mesh_space(a::MeshSamplingFaceNew) = mesh_space(mesh(a),Val(num_dims(a)))
mesh_quadrature(a::MeshSamplingFaceNew) = mesh_quadrature(quadrature_face(a))
mesh_face(a::MeshSamplingFaceNew) = a.mesh_face
local_face(a::MeshSamplingFaceNew) = a.local_face
mesh_sampling_face(a::MeshSamplingFaceNew) = a
quadrature_face(a::MeshSamplingFaceNew) = a.quadrature_face
mesh(a::MeshSamplingFaceNew) = mesh(quadrature_face(a))
num_dims(a::MeshSamplingFaceNew) = num_dims(mesh_face(a))
id(a::MeshSamplingFaceNew) = a.id
nodes(a::MeshSamplingFaceNew) = nodes(mesh_face(a))
node_coordinates(a::MeshSamplingFaceNew) = node_coordinates(mesh_face(a))

function replace_tabulators(f::typeof(value),a::MeshSamplingFaceNew,tabulated_values)
    MeshSamplingFaceNew(
                     a.mesh_face,
                     a.quadrature_face,
                     a.local_face,
                     tabulated_values,
                     a.tabulated_gradients,
                     a.reference_unit_normals)
end

function replace_tabulators(f::typeof(ForwardDiff.gradient),a::MeshSamplingFaceNew,tabulated_gradients)
    MeshSamplingFaceNew(
                     a.mesh_face,
                     a.quadrature_face,
                     a.local_face,
                     a.tabulated_values,
                     tabulated_gradients,
                     a.reference_unit_normals)
end

tabulators(f::typeof(value),a::MeshSamplingFaceNew) = a.tabulated_values
tabulators(f::typeof(ForwardDiff.gradient),a::MeshSamplingFaceNew) = a.tabulated_gradients

function replace_reference_unit_normals(a::MeshSamplingFaceNew,reference_unit_normals)
    MeshSamplingFaceNew(
                     a.mesh_face,
                     a.quadrature_face,
                     a.local_face,
                     a.tabulated_values,
                     a.tabulated_gradients,
                     reference_unit_normals)
end

reference_unit_normals(a::MeshSamplingFaceNew) = a.reference_unit_normals

function each_face_new(
    mesh::AbstractMesh,
    valD,
    mesh_quadrature::AbstractQuadrature,
    val_loop_dim=Val(num_dims(domain(mesh_quadrature))))

    loop_dim = val_parameter(val_loop_dim)
    domain_d = GT.domain(mesh_quadrature)
    d = num_dims(domain_d)
    D = val_parameter(valD)
    id = ANY_INDEX_NEW
    if d == loop_dim
        quadrature_face = QuadratureFaceNew(mesh_quadrature,id)
        mesh_face = MeshFaceNew(mesh,Val(D),id)
        local_face = nothing
        ids = GT.faces(domain_d)
        tabulated_values = nothing
        tabulated_gradients = nothing
        reference_unit_normals = nothing
        face0 = MeshSamplingFaceNew(
                                    mesh_face,
                                    quadrature_face,
                                    local_face,
                                    tabulated_values,
                                    tabulated_gradients,
                                    reference_unit_normals)
        face = setup_mesh_sampling_face(face0)
        return Each(update_quadrature_face_id,face,ids)
    end
    if D == loop_dim
        error("not implemented")
    end
    error()
end

function setup_mesh_sampling_face(a0)
    a1 = GT.tabulate(GT.value,a0)
    a2 = GT.tabulate(ForwardDiff.gradient,a1)
    D = num_dims(mesh_face(a0))
    d = num_dims(quadrature_face(a0))
    if d + 1 == D
        a3 = setup_unit_normal(a2)
    else
        a3 = a2
    end
    a3
end

function each_face_around_new(a::MeshSamplingFaceNew)
    n_faces_around = num_faces_around(a)
    ids = 1:n_faces_around
    Each(update_face_around_id,a,ids)
end

function num_faces_around(a::MeshSamplingFaceNew)
    dface = mesh_face(quadrature_face(a))
    D = num_dims(mesh_face(a))
    num_faces_around(dface,Val(D))
end

# NB update_id alone does not makes sense since there are several ids involved.

function update_quadrature_face_id(a::MeshSamplingFaceNew,id)
    quadrature_face_2 = update_id(quadrature_face(a),id)
    d = num_dims(quadrature_face(a))
    D = num_dims(mesh_face(a))
    if d == D
        mesh_face_2 = update_id(a.mesh_face,id)
    else
        mesh_face_2 = a.mesh_face
    end
    face0 = MeshSamplingFaceNew(
                     mesh_face_2,
                     quadrature_face_2,
                     a.local_face,
                     a.tabulated_values,
                     a.tabulated_gradients,
                     a.reference_unit_normals)
    domain = GT.domain(GT.domain(mesh_quadrature(a)))
    faces_around = GT.faces_around(domain)
    if faces_around === nothing
        face1 = face0
    else
        face_around = faces_around[inverse_faces(domain)[id]]
        face1 = update_face_around_id(face0,face_around)
    end
    face1
end

function update_face_around_id(a::MeshSamplingFaceNew,face_around_0)
    mesh = GT.mesh(a)
    topo = GT.topology(mesh)
    d = num_dims(quadrature_face(a))
    D = num_dims(mesh_face(a))
    dface = id(quadrature_face(a))
    domain = GT.domain(GT.domain(mesh_quadrature(a)))
    faces_around_permutation = GT.faces_around_permutation(domain)
    if faces_around_permutation !== nothing
        face = inverse_faces(domain)[dface]
        face_around = faces_around_permutation[face][face_around_0]
    else
        face_around = face_around_0
    end
    Dface = face_incidence(topo,d,D)[dface][face_around]
    mesh_face_2 = update_id(mesh_face(a),Dface)
    local_face_2 = GT.local_face(mesh_face_2,mesh_face(quadrature_face(a)))
    MeshSamplingFaceNew(
                     mesh_face_2,
                     a.quadrature_face,
                     local_face_2,
                     a.tabulated_values,
                     a.tabulated_gradients,
                     a.reference_unit_normals)
end

function update_mesh_face_id(a::MeshSamplingFaceNew,id)
end

function update_local_face_id(a::MeshSamplingFaceNew,local_id)
end

function each_face_new(mesh_quadrature::AbstractQuadrature)
    domain = GT.domain(mesh_quadrature)
    d = num_dims(domain)
    mesh = GT.mesh(domain)
    each_face_new(mesh,Val(d),mesh_quadrature)
end

struct SpaceSamplingFaceNew{A,B,C,D,E} <: AbstractSamplingFaceNew
    mesh_space::A
    mesh_sampling_face::B
    tabulated_values::C
    tabulated_gradients::D
    tabulated_jacobians::E
end
mesh_space(a::SpaceSamplingFaceNew) = a.mesh_space
mesh_sampling_face(a::SpaceSamplingFaceNew) = a.mesh_sampling_face
mesh_face(a::SpaceSamplingFaceNew) = mesh_face(mesh_sampling_face(a))
quadrature_face(a::SpaceSamplingFaceNew) = quadrature_face(mesh_sampling_face(a))
mesh(a::SpaceSamplingFaceNew) = mesh(quadrature_face(a))

function each_face_new(mesh_space::AbstractSpace,mesh_quadrature::AbstractQuadrature)
    mesh = GT.mesh(mesh_space)
    D =  GT.num_dims(mesh_space)
    mesh_faces = each_face_new(mesh,Val(D),mesh_quadrature)
    mesh_sampling_face = prototype(mesh_sampling_face)
    face = SpaceSamplingFaceNew(
        mesh_space,
        mesh_sampling_face,
        tabulated_values,
        tabulated_gradients,
        tabulated_jacobians)
    Each(update_mesh_sampling_face,face,mesh_faces)
end

function each_face_around_new(a::SpaceSamplingFaceNew)
    mesh_sampling_faces_around = each_face_around_new(mesh_sampling_face(a))
    Each(update_mesh_sampling_face,a,mesh_sampling_faces_around)
end

function update_mesh_sampling_face(a::SpaceSamplingFaceNew,mesh_sampling_face)
    msf = mesh_sampling_face(a)
    msf2 = update_id(msf,id)
    SpaceSamplingFaceNew(
        a.mesh_space,
        msf2,
        a.tabulated_values,
        a.tabulated_gradients,
        a.tabulated_jacobians,
       )
end

function replace_tabulators(f::typeof(value),a::SpaceSamplingFaceNew,tabulated_values)
    SpaceSamplingFaceNew(
        a.mesh_space,
        a.mesh_sampling_face,
        tabulated_values,
        a.tabulated_gradients,
        a.tabulated_jacobians,
       )
end

function replace_tabulators(f::typeof(ForwardDiff.gradient),a::SpaceSamplingFaceNew,tabulated_gradients)
    SpaceSamplingFaceNew(
        a.mesh_space,
        a.mesh_sampling_face,
        a.tabulated_values,
        tabulated_gradients,
        a.tabulated_jacobians,
       )
end

function replace_tabulators(f::typeof(ForwardDiff.jacobian),a::SpaceSamplingFaceNew,tabulated_jacobians)
    SpaceSamplingFaceNew(
        a.mesh_space,
        a.mesh_sampling_face,
        a.tabulated_values,
        a.tabulated_gradients,
        tabulated_jacobians,
       )
end

struct DiscreteFieldSamplingFaceNew{A,B} <: AbstractSamplingFaceNew
    mesh_discrete_field::A
    space_sampling_face::B
end
mesh_space(a::DiscreteFieldSamplingFaceNew) = a.mesh_space
mesh_sampling_face(a::DiscreteFieldSamplingFaceNew) = a.mesh_sampling_face
mesh_face(a::DiscreteFieldSamplingFaceNew) = mesh_face(mesh_sampling_face(a))
quadrature_face(a::DiscreteFieldSamplingFaceNew) = quadrature_face(mesh_sampling_face(a))
mesh(a::DiscreteFieldSamplingFaceNew) = mesh(quadrature_face(a))




function reference_id(f::AbstractMeshFaceNew)
    mesh = GT.mesh(f)
    d = num_dims(f)
    f_id = GT.id(f)
    face_reference_id(mesh,d)[f_id]
end

function reference_id(f::AbstractSpaceFaceNew)
    space = GT.mesh_space(f)
    f_id = GT.id(f)
    face_reference_id(space)[f_id]
end

function reference_id(f::AbstractQuadratureFaceNew)
    q = GT.mesh_quadrature(f)
    f_id = GT.id(f)
    face_reference_id(q)[f_id]
end

function reference_space(f::AbstractMeshFaceNew)
    mesh = GT.mesh(f)
    rid = reference_id(f)
    d = num_dims(f)
    reference_spaces(mesh,Val(d))[rid]
end

function reference_space(f::AbstractSpaceFaceNew)
    space = GT.mesh_space(f)
    rid = reference_id(f)
    reference_spaces(space)[rid]
end

function reference_quadrature(f::AbstractFaceNew)
    reference_quadrature(quadrature_face(f))
end

function reference_quadrature(f::AbstractQuadratureFaceNew)
    q = GT.mesh_quadrature(f)
    rid = reference_id(f)
    reference_quadratures(q)[rid]
end

function reference_domain(f::AbstractMeshFaceNew)
    domain(reference_space(f))
end

function reference_domain(f::AbstractSpaceFaceNew)
    domain(reference_space(f))
end

function reference_domain(f::AbstractQuadratureFaceNew)
    domain(reference_quadrature(f))
end

function reference_shape_functions(f::AbstractFaceNew)
    shape_functions(reference_space(f))
end

function shape_functions(f::AbstractFaceNew)
    error("not implemented")
end

function nodes(f::AbstractMeshFaceNew)
    d = num_dims(f)
    mesh = GT.mesh(f)
    f_id = GT.id(f)
    face_nodes(mesh,d)[f_id]
end

function nodes(f::AbstractSpaceFaceNew)
    mesh = GT.space(f)
    f_id = GT.id(f)
    face_nodes(mesh)[f_id]
end

function num_nodes(a::AbstractFaceNew)
    length(nodes(a))
end

function node_coordinates(a::AbstractMeshFaceNew)
    mesh = GT.mesh(a) 
    node_x = node_coordinates(mesh)
    view(node_x,nodes(a))
end

function node_coordinates(a::AbstractSpaceFaceNew)
    mesh = GT.space(a) 
    node_x = node_coordinates(mesh)
    view(node_x,nodes(a))
end

function permutation_id(A::AbstractMeshFaceNew,a::AbstractMeshFaceNew)
    T = GT.topology(GT.mesh(A))
    D = GT.num_dims(A)
    d = GT.num_dims(a)
    face_permutation_ids(T,D,d)[id(A)][id(a)]
end

function node_permutation(A::AbstractMeshFaceNew,a::AbstractMeshFaceNew)
    k = permutation_id(A,a)
    Ps = GT.node_permutations(reference_space(a))
    Ps[k]
end

function barycenter(a::AbstractFaceNew)
    x = node_coordinates(a)
    sum(x)/length(x)
end

function diameter(a::AbstractFaceNew)
    nodes = GT.nodes(a)
    mesh = GT.mesh(a)
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

function tabulate(f,face::AbstractSamplingFaceNew)
    D = num_dims(mesh_face(face))
    d = num_dims(quadrature_face(face))
    aligned = Val(d==D)
    tabulate(f,face,aligned)
end

function tabulate(f,face::AbstractSamplingFaceNew,aligned::Val{true})
    space = GT.mesh_space(face)
    quadrature = GT.mesh_quadrature(face)
    rid_point_x = map(coordinates,reference_quadratures(quadrature))
    # NB the TODOs below can be solved by introducing an extra nesting level
    # TODO assumes same reference elems for integration and for interpolation
    rid_to_tab = map(rid_point_x,reference_spaces(space)) do point_to_x, refface
        collect(permutedims(tabulator(refface)(f,point_to_x))) # TODO fix this globally
    end
    replace_tabulators(f,face,rid_to_tab)
end

function tabulate(f,face::AbstractSamplingFaceNew,aligned::Val{false})
    space = GT.mesh_space(face)
    quadrature = GT.mesh_quadrature(face)
    D = num_dims(space)
    d = num_dims(quadrature)
    mesh = GT.mesh(space)
    topo = topology(mesh)
    drid_refdface = reference_spaces(mesh,Val(d))
    Drid_refDface = reference_spaces(mesh,Val(D))
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
    replace_tabulators(f,face,drid_Drid_ldface_perm_tab)
end

function reference_map(refdface,refDface)
    D = num_dims(refDface)
    d = num_dims(refdface)
    dof_to_f = shape_functions(refdface)
    boundary = refDface |> GT.domain |> GT.mesh
    lface_to_nodes = GT.face_nodes(boundary,d)
    node_to_coords = GT.node_coordinates(boundary)
    lface_to_lrefid = GT.face_reference_id(boundary,d)
    lrefid_to_lrefface = GT.reference_spaces(boundary,d)
    lrefid_to_perm_to_ids = map(GT.node_permutations,lrefid_to_lrefface)
    
    func_template = (dof_to_coeff, dof_to_f, ndofs) -> (x -> sum(dof->dof_to_coeff[dof]*dof_to_f[dof](x),1:ndofs))
    d2c = node_to_coords[lface_to_nodes[1][ lrefid_to_perm_to_ids[lface_to_lrefid[1]][1] ]]
    proto = func_template(d2c, dof_to_f, length(d2c))

    result::Vector{Vector{typeof(proto)}} = map(1:GT.num_faces(boundary,d)) do lface
        lrefid = lface_to_lrefid[lface]
        nodes = lface_to_nodes[lface]
        perm_to_ids = lrefid_to_perm_to_ids[lrefid]
        map(perm_to_ids) do ids
            #dof_to_coeff = node_to_coords[nodes[ids]] # Wrong
            dof_to_coeff = similar(node_to_coords[nodes])
            dof_to_coeff[ids] = node_to_coords[nodes]
            ndofs = length(dof_to_coeff)
            result_inner::typeof(proto) = func_template(dof_to_coeff, dof_to_f, ndofs)
        end
    end
    return result
end

function tabulator(f,face::AbstractSamplingFaceNew)
    D = num_dims(mesh_face(face))
    d = num_dims(quadrature_face(face))
    aligned = Val(d==D)
    tabulator(f,face,aligned)
end

function tabulator(f,face::AbstractSamplingFaceNew,aligned::Val{true})
    tabulators = GT.tabulators(f,face)
    drid_tab = tabulators
    dface = face
    drid = reference_id(quadrature_face(dface))
    tab = drid_tab[drid]
end

function tabulator(f,face::AbstractSamplingFaceNew,aligned::Val{false})
    tabulators = GT.tabulators(f,face)
    Dface = mesh_face(face)
    dface = mesh_face(quadrature_face(face))
    ldface = local_face(face)
    perm = permutation_id(Dface,ldface)
    drid = reference_id(dface)
    Drid = reference_id(Dface)
    drid_Drid_ldface_perm_tab = tabulators
    drid_Drid_ldface_perm_tab[drid][Drid][id(ldface)][perm]
end

function shape_functions(f,point::AbstractPointNew)
    ref_funs = reference_shape_functions(f,point)
    phys_funs = shape_functions_allocation(f,point)
    map_shape_functions!(f,phys_funs,ref_funs)
    phys_funs
end

function reference_shape_functions(f,point::AbstractPointNew)
    tabulator = GT.tabulator(f,parent_face(point))
    p_id = id(point)
    ref_funs = view(tabulator,:,p_id)
    ref_funs
end

function num_points(a::AbstractFaceNew)
    num_points(reference_quadrature(a))
end

function coordinate(a::AbstractPointNew)
    coordinate(GT.value,a)
end

function ForwardDiff.jacobian(a::AbstractPointNew)
    coordinate(ForwardDiff.jacobian,a)
end

function mesh_sampling_point(point::AbstractPointNew)
    face = parent_face(point)
    PointNew(mesh_sampling_face(face),id(point))
end

function coordinate(::typeof(value),point::AbstractPointNew)
    a = mesh_sampling_point(point)
    s = reference_shape_functions(value,a)
    x = node_coordinates(parent_face(a))
    n = length(s)
    sum(i->x[i]*s[i],1:n)
end

function coordinate(::typeof(ForwardDiff.jacobian),point::AbstractPointNew)
    a = mesh_sampling_point(point)
    s = reference_shape_functions(ForwardDiff.gradient,a)
    x = node_coordinates(parent_face(a))
    n = length(s)
    sum(i->outer(x[i],s[i]),1:n)
end

function reference_weight(a::AbstractPointNew)
    quadrature = GT.reference_quadrature(parent_face(a))
    weights(quadrature)[id(a)]
end

function weight(a::AbstractPointNew)
    w = reference_weight(a)
    J = coordinate(ForwardDiff.jacobian,a)
    change_of_measure(J)*w
end

function field(f,point::AbstractPointNew)
    values = GT.values(face(point))
    funs = GT.shape_functions(f,point)
    nvalues = length(values)
    sum(1:nvalues) do i
        values[i]*funs[i]
    end
end

function setup_unit_normal(a::AbstractFaceNew)
    d = num_dims(quadrature_face(a))
    D = num_dims(mesh_face(a))
    mesh = GT.mesh(a)
    @assert d + 1 == D
    Drid_ldface_nref = map(GT.reference_spaces(mesh,D)) do refface
        GT.normals(GT.mesh(GT.domain(refface)))# TODO rename normals?
    end
    replace_reference_unit_normals(a,Drid_ldface_nref)
end

function reference_unit_normal(a::AbstractFaceNew)
    Drid_ldface_nref = reference_unit_normals(a)
    Dface = mesh_face(a)
    ldface = local_face(a)
    Drid = reference_id(Dface)
    Drid_ldface_nref[Drid][id(ldface)]
end

function unit_normal(a::AbstractPointNew)
    nref = reference_unit_normal(parent_face(a))
    J = coordinate(ForwardDiff.jacobian,a)
    map_unit_normal(J,nref)
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



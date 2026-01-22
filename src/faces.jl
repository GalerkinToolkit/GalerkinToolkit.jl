
abstract type AbstractFaceNew <: AbstractType end
abstract type AbstractMeshFaceNew <: AbstractFaceNew end
abstract type AbstractLocalMeshFaceNew <: AbstractMeshFaceNew end
abstract type AbstractQuadratureFaceNew <: AbstractFaceNew end
abstract type AbstractSpaceFaceNew <: AbstractFaceNew end
abstract type AbstractDiscreteFieldFaceNew <: AbstractSpaceFaceNew end
abstract type AbstractSampingFaceNew <: AbstractFaceNew end

abstract type AbstractPointNew <: AbstractType end
abstract type AbstractQuadraturePointNew <: AbstractPointNew end
abstract type AbstractSamplingPointNew <: AbstractPointNew end

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
prototype(x::Each) = x.prototype

struct AnyIndex <: AbstractType end
const ANY_INDEX = AnyIndex()

struct MeshFaceNew <: AbstractMeshFaceNew
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
    id = ANY_INDEX
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

function each_face_around_new(dface::AbstractMeshFaceNew,valD)
    D = val_parameter(valD)
    d = GT.num_dims(dface)
    mesh = GT.mesh(dface)
    topo = topology(mesh)
    ids = face_incidence(topo,d,D)[id(dface)]
    num_dims = Val(D)
    id = ANY_INDEX
    Dface = MeshFaceNew(mesh,num_dims,id)
    Each(update_id,Dface,ids)
end

struct LocalMeshFaceNew{A,B} <: AbstractLocalMeshFaceNew
    local_face::A
    parent_face::B
end
mesh(f::MeshLocalFaceNew) = mesh(local_face(f))
id(f::MeshLocalFaceNew) = id(local_face(f))
num_dims(f::MeshLocalFaceNew) = num_dims(local_face(f))
function update_id(face::MeshLocalFaceNew,id)
    local_face = update_id(GT.local_face(face),id)
    LocalMeshFaceNew(local_face,parent_face(face))
end
local_face(a::LocalMeshFaceNew) = a.local_face
parent_face(a::LocalMeshFaceNew) = a.parent_face
function global_face(a::LocalMeshFaceNew)
    d = num_dims(a)
    ldid = id(a)
    parent_face = GT.parent_face(a)
    D = num_dims(parent_face)
    Did = id(parent_face)
    mesh = GT.mesh(a)
    topo = topology(mesh)
    did = face_incidence(topo,D,d)[Did][ldid]
    MeshFaceNew(mesh,Val(d),did)
end

function each_local_face(Df::MeshFaceNew,vald)
    d = val_parameter(vald)
    D = GT.num_dims(Df)
    l_mesh = GT.reference_domain(Df)
    l_id = ANY_INDEX
    dface = MeshFaceNew(l_mesh,Val(d),l_id)
    l_dface = LocalMeshFaceNew(dface,Df)
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

struct QuadratureFaceNew{A,B} <: AbstractQuadratureFaceNew
    mesh_quadrature::A
    id::B
end
mesh_quadrature(a::QuadratureFaceNew) = a.mesh_quadrature
mesh(a::QuadratureFaceNew) = mesh(domain(mesh_quadrature(a)))
id(a::QuadratureFaceNew) = a.id
update_id(a::QuadratureFaceNew,id) = QuadratureFaceNew(mesh_quadrature(a),id)
function update_id(a::QuadratureFaceNew,face2::MeshFaceNew)
    q = mesh_quadrature(a)
    @assert num_dims(q) == num_dims(face_2)
    QuadratureFaceNew(q,GT.id(face2))
end

struct QuadraturePointNew{A,B} <: AbstractQuadraturePointNew
    parent_face::A
    id::B
end
parent_face(a::QuadraturePointNew) = a.parent_face
id(a::QuadraturePointNew) = a.id
update_id(a::QuadraturePointNew,id) = QuadraturePointNew(parent_face(a),id)
update_id(a::QuadraturePointNew,point::AbstractPointNew) = QuadraturePointNew(parent_face(a),id(point))

function each_point_new(face::QuadratureFaceNew)
    id = ANY_INDEX
    point = QuadraturePointNew(face,id)
    npoints = num_points(face)
    Each(update_id,point,1:npoints)
end

struct MeshSamplingFaceNew{A,B,C,D} <: AbstractSamplingFaceNew
    mesh_face::A
    quadrature_face::B
    tabulated_values::C
    tabulated_gradients::D
end
mesh_face(a::MeshSamplingFaceNew) = a.mesh_face
quadrature_face(a::MeshSamplingFaceNew) = a.quadrature_face
mesh(a::MeshSamplingFaceNew) = mesh(quadrature_face(a))
num_dims(a::MeshSamplingFaceNew) = val_parameter(a.num_dims)
id(a::MeshSamplingFaceNew) = a.id
function update_id(a::MeshSamplingFaceNew,id)
    update_id(a,id,mesh_face(a))
end
function update_id(a::MeshSamplingFaceNew,id,mesh_face::AbstractMeshFaceNew)
    mesh_face_2 = update_id(mesh_face,id)
    quadrature_face_2 = update_id(quadrature_face(a),id)
    MeshSamplingFaceNew(
                     mesh_face_2,
                     quadrature_face_2,
                     a.tabulated_values,
                     a.tabulated_gradients)
end
function update_id(a::MeshSamplingFaceNew,id::AbstractLocalMeshFaceNew,mesh_face::AbstractLocalMeshFaceNew)
    mesh_face_2 = id
    @assert num_dims(mesh_face_2) == num_dims(mesh_face)
    global_face_2 = GT.global_face(mesh_face_2)
    quadrature_face_2 = update_id(quadrature_face(a),global_face_2)
    MeshSamplingFaceNew(
                     mesh_face_2,
                     quadrature_face_2,
                     a.tabulated_values,
                     a.tabulated_gradients)
end

function replace_tabulators(f::typeof(value),a::MeshSamplingFaceNew,tabulated_values)
    MeshSamplingFaceNew(
                     a.mesh_face,
                     a.quadrature_face,
                     tabulated_values,
                     a.tabulated_gradients)
end

function replace_tabulators(f::typeof(ForwardDiff.gradient),a::MeshSamplingFaceNew,tabulated_gradients)
    MeshSamplingFaceNew(
                     a.mesh_face,
                     a.quadrature_face,
                     a.tabulated_values,
                     tabulated_gradients)
end

function each_face_new(mesh::AbstractMesh,valD,mesh_quadrature::AbstractQuadrature)
    d = num_dims(mesh_quadrature)
    D = val_parameter(valD)
    q_faces = each_face_new(mesh_quadrature)
    quadrature_face = prototype(q_faces)
    if d == D
        mesh_face = MeshFaceNew(mesh,Val(D),id(quadrature_face))
        ids = ids(q_faces)
    else
        local_face = MeshFaceNew(mesh,Val(d),ANY_INDEX)
        parent_face = MeshFaceNew(mesh,Val(D),ANY_INDEX)
        mesh_face = LocalMeshFaceNew(local_face,parent_face)
        ids = nothing
    end
    tabulated_values = nothing
    tabulated_grediants = nothing
    face = MeshSamplingFaceNew(mesh_face,quadrature_face,tabulated_values,tabulated_grediants)
    Each(update_id,face,ids)
end

struct MeshSamplingPointNew{A,B} <: AbstractSamplingPointNew
    parent_face::A
    quadrature_point::B
end
parent_face(a::MeshSamplingPointNew) = a.parent_face
quadrature_point(a::MeshSamplingPointNew) = a.quadrature_point
id(a::MeshSamplingPointNew) = id(quadrature_point(a))
function update_id(a::MeshSamplingPointNew,id)
    q_point = update_id(quadrature_point(a),id)
    MeshSamplingPointNew(a.parent_face,q_point)
end

function each_point_new(face::MeshSamplingFaceNew)
    q_points = each_point_new(quadrature_face(face))
    q_point = prototype(q_point)
    ids = GT.ids(q_points)
    point = MeshSamplingPointNew(face,q_point)
    Each(update_id,point,ids)
end

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

function reference_quadrature(f::AbstractQuadratureFaceNew)
    q = GT.mesh_quadrature(f)
    rid = reference_id(f)
    reference_spaces(q)[rid]
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

function barycenter(a::AbstractFaceNew)
    x = node_coordinates(a)
    sum(x)/length(x)
end

function tabulate(f,face::AbstractSamplingFaceNew)
    tabulate(f,face,mesh_face(face))
end

function tabulate(f,face::AbstractSamplingFaceNew,mesh_face::AbstractMeshFaceNew)
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

function tabulate(f,face::AbstractSamplingFaceNew,mesh_face::AbstractLocalMeshFaceNew)
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
    replace_tabulators(f,a,drid_Drid_ldface_perm_tab)
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
    tabulator(f,face,mesh_face(f))
end

function tabulator(f,face::AbstractSamplingFaceNew,mesh_face::AbstractMeshFaceNew)
    drid_tab = tabulators
    dface = face
    drid = reference_id(dface)
    tab = drid_tab[drid]
end

function tabulator(f,face::AbstractSamplingFaceNew,mesh_face::AbstractLocalMeshFaceNew)
    tabulators = GT.tabulators(f,face)
    ldface = face
    Dface = parent_face(ldface)
    dface = global_face(ldface)
    perm = permutation_id(ldface)
    drid = reference_id(dface)
    Drid = reference_id(Dface)
    drid_Drid_ldface_perm_tab = tabulators
    drid_Drid_ldface_perm_tab[drid][Drid][id(ldface)][perm]
end

function shape_functions(f,point::AbstractSamplingPointNew)
    ref_funs = reference_shape_functions(f,point)
    phys_funs = shape_functions_allocation(f,point)
    map_shape_functions!(f,phys_funs,ref_funs)
    phys_funs
end

function reference_shape_functions(f,point::AbstractSamplingPointNew)
    tabulator = GT.tabulator(f,parent_face(point))
    p_id = id(point)
    ref_funs = view(tabulator,:,p)
    ref_funs
end

function num_points(a::AbstractFaceNew)
    num_points(quadrature(a))
end

function coordinate(a::AbstractSamplingPointNew)
    coordinate(GT.value,a)
end

function ForwardDiff.jacobian(a::AbstractSamplingPointNew)
    coordinate(ForwardDiff.jacobian,a)
end

function coordinate(::typeof(value),face::AbstractSamplingPointNew)
    a = mesh_sampling_point(face)
    s = shape_functions(value,a)
    x = node_coordinates(a)
    n = length(s)
    sum(i->x[i]*s[i],1:n)
end

function coordinate(::typeof(ForwardDiff.jacobian),face::AbstractPointNew)
    a = mesh_sampling_point(face)
    s = shape_functions(ForwardDiff.gradient,a)
    x = node_coordinates(a)
    n = num_nodes(a)
    sum(i->outer(x[i],s[i]),1:n)
end

function reference_weight(a::AbstractPointNew)
    quadrature = GT.quadrature(a)
    weights(quadrature)[id(a)]
end

function weight(a::AbstractPointNew)
    w = weight(reference(a))
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

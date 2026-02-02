
abstract type AbstractFaceNew <: AbstractType end
abstract type AbstractMeshFaceNew <: AbstractFaceNew end
abstract type AbstractQuadratureFaceNew <: AbstractFaceNew end
abstract type AbstractSpaceFaceNew <: AbstractFaceNew end
abstract type AbstractDiscreteFieldFaceNew <: AbstractSpaceFaceNew end
abstract type AbstractSamplingFaceNew <: AbstractFaceNew end
abstract type AbstractPointNew <: AbstractType end


struct Each{A,B,C}
    update::A
    prototype::B
    iterator::C #
end
Base.length(iter::Each) = length(iter.iterator)
Base.isdone(iter::Each,id) = id > length(iter)
function Base.getindex(iter::Each,i::Integer)
    item = iter.iterator[i]
    iter.update(iter.prototype,item)
end
function Base.iterate(iter::Each,i=1)
    if Base.isdone(iter,i)
        nothing
    else
        accessor = iter[i]
        (accessor,i+1)
    end
end
Base.keys(iter::Each) = LinearIndices((length(iter),))
prototype(x::Each) = x.prototype
iterator(x::Each) = x.iterator
function tabulate(f,x::Each)
    Each(x.update,tabulate(f,x.prototype),x.iterator)
end

struct AnyIndexNew <: AbstractType end
const ANY_INDEX_NEW = AnyIndexNew()
#function Base.getindex(a::Vector,i::AnyIndexNew)
#    first(a)
#end
#function Base.getindex(a::JaggedArray,i::AnyIndexNew)
#    first(a)
#end
#function Base.getindex(a::Matrix,i::AnyIndexNew,j::AnyIndexNew)
#    first(a)
#end

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
#function update_id(face::MeshFaceNew,face2::MeshFaceNew)
#    @assert num_dims(face) == num_dims(face_2)
#    id = GT.id(face2)
#    MeshFaceNew(face.mesh,face.num_dims,id)
#end

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

struct SpaceFaceNew{A,B} <: AbstractSpaceFaceNew
    mesh_space::A
    id::B
end
mesh_space(a::SpaceFaceNew) = a.mesh_space
id(a::SpaceFaceNew) = a.id

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

allocate_shape_funcions(f,a::MeshSamplingFaceNew) = a

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

struct SpaceSamplingFaceNew{A,B,C,D,E,F,G,H} <: AbstractSamplingFaceNew
    mesh_space::A
    mesh_sampling_face::B
    tabulated_values::C
    tabulated_gradients::D
    tabulated_jacobians::E
    shape_function_values::F
    shape_function_gradients::G
    shape_function_jacobians::H
end
mesh_space(a::SpaceSamplingFaceNew) = a.mesh_space
mesh_sampling_face(a::SpaceSamplingFaceNew) = a.mesh_sampling_face
space_sampling_face(a::SpaceSamplingFaceNew) = a
mesh_quadrature(a::SpaceSamplingFaceNew) = mesh_quadrature(a.mesh_sampling_face)
mesh_face(a::SpaceSamplingFaceNew) = mesh_face(mesh_sampling_face(a))
local_face(a::SpaceSamplingFaceNew) = local_face(mesh_sampling_face(a))
quadrature_face(a::SpaceSamplingFaceNew) = quadrature_face(mesh_sampling_face(a))
mesh(a::SpaceSamplingFaceNew) = mesh(quadrature_face(a))
function space_face(a::SpaceSamplingFaceNew)
    space = GT.mesh_space(a)
    f = mesh_face(a)
    SpaceFaceNew(space,id(f))
end
num_faces_around(a::SpaceSamplingFaceNew) = num_faces_around(mesh_sampling_face(a))

function setup_space_sampling_face(a0::SpaceSamplingFaceNew,tabulate)
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
    a3
end

function each_face_new(mesh_space::AbstractSpace,mesh_quadrature::AbstractQuadrature;tabulate=())
    mesh = GT.mesh(mesh_space)
    D =  GT.num_dims(mesh_space)
    mesh_sampling_faces = each_face_new(mesh,Val(D),mesh_quadrature)
    mesh_sampling_face = prototype(mesh_sampling_faces)
    ids = GT.iterator(mesh_sampling_faces)
    tabulated_values = nothing
    tabulated_gradients = nothing
    tabulated_jacobians = nothing
    shape_function_values = nothing
    shape_function_gradients = nothing
    shape_function_jacobians = nothing
    face0 = SpaceSamplingFaceNew(
                                 mesh_space,
                                 mesh_sampling_face,
                                 tabulated_values,
                                 tabulated_gradients,
                                 tabulated_jacobians,
                                 shape_function_values,
                                 shape_function_gradients,
                                 shape_function_jacobians,
                                )
    face = setup_space_sampling_face(face0,tabulate)
    Each(update_quadrature_face_id,face,ids)
end

function update_quadrature_face_id(a::SpaceSamplingFaceNew,id)
    msf = GT.mesh_sampling_face(a)
    mesh_sampling_face = update_quadrature_face_id(msf,id)
    a1 = SpaceSamplingFaceNew(
        a.mesh_space,
        mesh_sampling_face,
        a.tabulated_values,
        a.tabulated_gradients,
        a.tabulated_jacobians,
        a.shape_function_values,
        a.shape_function_gradients,
        a.shape_function_jacobians,
       )
    domain = GT.domain(GT.domain(mesh_quadrature(msf)))
    faces_around = GT.faces_around(domain)
    if faces_around === nothing
        a2 = a1
    else
        face_around = faces_around[inverse_faces(domain)[id]]
        a2 = update_face_around_id_allocation(a1,face_around)
    end
    a2
end

function each_face_around_new(a::SpaceSamplingFaceNew)
    n_faces_around = num_faces_around(a)
    ids = 1:n_faces_around
    Each(update_face_around_id,a,ids)
end

function update_face_around_id(a::SpaceSamplingFaceNew,face_around_id)
    a1 = update_face_around_id_parent(a,face_around_id)
    a2 = update_face_around_id_allocation(a1,face_around_id)
end

function update_face_around_id_parent(a::SpaceSamplingFaceNew,face_around_id)
    msf = GT.mesh_sampling_face(a)
    mesh_sampling_face = update_face_around_id(msf,face_around_id)
    SpaceSamplingFaceNew(
                         a.mesh_space,
                         mesh_sampling_face,
                         a.tabulated_values,
                         a.tabulated_gradients,
                         a.tabulated_jacobians,
                         a.shape_function_values,
                         a.shape_function_gradients,
                         a.shape_function_jacobians,
                        )
end

function update_face_around_id_allocation(a::SpaceSamplingFaceNew,face_around_id)
    v = a.shape_function_values
    g = a.shape_function_gradients
    j = a.shape_function_jacobians
    SpaceSamplingFaceNew(
        a.mesh_space,
        a.mesh_sampling_face,
        a.tabulated_values,
        a.tabulated_gradients,
        a.tabulated_jacobians,
        v === nothing ? v : v[face_around_id],
        g === nothing ? g : g[face_around_id],
        j === nothing ? j : j[face_around_id],
       )
end

function replace_tabulators(f::typeof(value),a::SpaceSamplingFaceNew,tabulated_values)
    SpaceSamplingFaceNew(
        a.mesh_space,
        a.mesh_sampling_face,
        tabulated_values,
        a.tabulated_gradients,
        a.tabulated_jacobians,
        a.shape_function_values,
        a.shape_function_gradients,
        a.shape_function_jacobians,
       )
end

function replace_tabulators(f::typeof(ForwardDiff.gradient),a::SpaceSamplingFaceNew,tabulated_gradients)
    SpaceSamplingFaceNew(
        a.mesh_space,
        a.mesh_sampling_face,
        a.tabulated_values,
        tabulated_gradients,
        a.tabulated_jacobians,
        a.shape_function_values,
        a.shape_function_gradients,
        a.shape_function_jacobians,
       )
end

function replace_tabulators(f::typeof(ForwardDiff.jacobian),a::SpaceSamplingFaceNew,tabulated_jacobians)
    SpaceSamplingFaceNew(
        a.mesh_space,
        a.mesh_sampling_face,
        a.tabulated_values,
        a.tabulated_gradients,
        tabulated_jacobians,
        a.shape_function_values,
        a.shape_function_gradients,
        a.shape_function_jacobians,
       )
end

function replace_shape_function_allocation(f::typeof(value),a::SpaceSamplingFaceNew,shape_function_values)
    SpaceSamplingFaceNew(
        a.mesh_space,
        a.mesh_sampling_face,
        a.tabulated_values,
        a.tabulated_gradients,
        a.tabulated_jacobians,
        shape_function_values,
        a.shape_function_gradients,
        a.shape_function_jacobians,
       )
end

function replace_shape_function_allocation(f::typeof(ForwardDiff.gradient),a::SpaceSamplingFaceNew,shape_function_gradients)
    SpaceSamplingFaceNew(
        a.mesh_space,
        a.mesh_sampling_face,
        a.tabulated_values,
        a.tabulated_gradients,
        a.tabulated_jacobians,
        a.shape_function_values,
        shape_function_gradients,
        a.shape_function_jacobians,
       )
end

function replace_shape_function_allocation(f::typeof(ForwardDiff.jacobian),a::SpaceSamplingFaceNew,shape_function_jacobians)
    SpaceSamplingFaceNew(
        a.mesh_space,
        a.mesh_sampling_face,
        a.tabulated_values,
        a.tabulated_gradients,
        a.tabulated_jacobians,
        a.shape_function_values,
        a.shape_function_gradients,
        shape_function_jacobians,
       )
end

function tabulators(::typeof(value),a::SpaceSamplingFaceNew)
    a.tabulated_values
end

function tabulators(::typeof(ForwardDiff.gradient),a::SpaceSamplingFaceNew)
    a.tabulated_gradients
end

function tabulators(::typeof(ForwardDiff.jacobian),a::SpaceSamplingFaceNew)
    a.tabulated_jacobians
end

function shape_functions_allocation(::typeof(value),a::SpaceSamplingFaceNew)
    a.shape_function_values
end

function shape_functions_allocation(::typeof(ForwardDiff.gradient),a::SpaceSamplingFaceNew)
    a.shape_function_gradients
end

function shape_functions_allocation(::typeof(ForwardDiff.jacobian),a::SpaceSamplingFaceNew)
    a.shape_function_jacobians
end

struct DiscreteFieldSamplingFaceNew{A,B,C} <: AbstractSamplingFaceNew
    mesh_discrete_field::A
    space_sampling_face::B
    values_allocation::C
end
mesh_space(a::DiscreteFieldSamplingFaceNew) = a.mesh_space
space_sampling_face(a::DiscreteFieldSamplingFaceNew) = a.space_sampling_face
mesh_discrete_field(a::DiscreteFieldSamplingFaceNew) = a.mesh_discrete_field
mesh_sampling_face(a::DiscreteFieldSamplingFaceNew) = mesh_sampling_face(space_sampling_face(a))
mesh_face(a::DiscreteFieldSamplingFaceNew) = mesh_face(mesh_sampling_face(a))
local_face(a::DiscreteFieldSamplingFaceNew) = local_face(mesh_sampling_face(a))
space_face(a::DiscreteFieldSamplingFaceNew) = space_face(space_sampling_face(a))
quadrature_face(a::DiscreteFieldSamplingFaceNew) = quadrature_face(mesh_sampling_face(a))
mesh(a::DiscreteFieldSamplingFaceNew) = mesh(quadrature_face(a))
mesh_space(a::DiscreteFieldSamplingFaceNew) = mesh_space(space_sampling_face(a))
num_faces_around(a::DiscreteFieldSamplingFaceNew) = num_faces_around(mesh_sampling_face(a))

function replace_values_allocation(a::DiscreteFieldSamplingFaceNew,values)
    DiscreteFieldSamplingFaceNew(
                                 a.mesh_discrete_field,
                                 a.space_sampling_face,
                                 values
                                )
end
values_allocation(a::DiscreteFieldSamplingFaceNew) = a.values_allocation
tabulators(f,a::DiscreteFieldSamplingFaceNew) = tabulators(f,space_sampling_face(a))
shape_functions_allocation(f,a::DiscreteFieldSamplingFaceNew) = shape_functions_allocation(f,space_sampling_face(a))

function each_face_new(uh::DiscreteField,mesh_quadrature::AbstractQuadrature;kwargs...)
    V = GT.space(uh)
    V_faces = each_face_new(V,mesh_quadrature;kwargs...)
    space_sampling_face = GT.prototype(V_faces)
    ids = GT.iterator(V_faces)
    values_allocation = nothing
    face0 = DiscreteFieldSamplingFaceNew(uh,space_sampling_face,values_allocation)
    face = allocate_values(face0)
    Each(update_quadrature_face_id,face,ids)
end

function update_quadrature_face_id(a::DiscreteFieldSamplingFaceNew,id)
    ssf = GT.space_sampling_face(a)
    space_sampling_face = update_quadrature_face_id(ssf,id)
    a1 = DiscreteFieldSamplingFaceNew(
                                      a.mesh_discrete_field,
                                      space_sampling_face,
                                      a.values_allocation
                                     )
    domain = GT.domain(GT.domain(mesh_quadrature(ssf)))
    faces_around = GT.faces_around(domain)
    if faces_around === nothing
        a2 = a1
    else
        face_around = faces_around[inverse_faces(domain)[id]]
        a2 = update_face_around_id_allocation(a1,face_around)
    end
    a2
end

function each_face_around_new(a::DiscreteFieldSamplingFaceNew)
    n_faces_around = num_faces_around(a)
    ids = 1:n_faces_around
    Each(update_face_around_id,a,ids)
end

function update_face_around_id(a::DiscreteFieldSamplingFaceNew,face_around_id)
    a1 = update_face_around_id_parent(a,face_around_id)
    a2 = update_face_around_id_allocation(a1,face_around_id)
end

function update_face_around_id_parent(a::DiscreteFieldSamplingFaceNew,face_around_id)
    ssf = GT.space_sampling_face(a)
    space_sampling_face = update_face_around_id(ssf,face_around_id)
    a1 = DiscreteFieldSamplingFaceNew(
                                      a.mesh_discrete_field,
                                      space_sampling_face,
                                      a.values_allocation
                                     )
end

function update_face_around_id_allocation(a::DiscreteFieldSamplingFaceNew,face_around_id)
    v = a.values_allocation
    a1 = DiscreteFieldSamplingFaceNew(
                                      a.mesh_discrete_field,
                                      a.space_sampling_face,
                                      v === nothing ? v : v[face_around_id],
                                     )
end

function reference_id(f::AbstractMeshFaceNew)
    mesh = GT.mesh(f)
    d = num_dims(f)
    f_id = GT.id(f)
    fr = face_reference_id(mesh,d)
    if f_id isa AnyIndexNew
        one(eltype(fr))
    else
        fr[f_id]
    end
end

function reference_id(f::AbstractSpaceFaceNew)
    space = GT.mesh_space(f)
    f_id = GT.id(f)
    fr = face_reference_id(space)
    if f_id isa AnyIndexNew
        one(eltype(fr))
    else
        fr[f_id]
    end
end

function reference_id(f::AbstractQuadratureFaceNew)
    q = GT.mesh_quadrature(f)
    f_id = GT.id(f)
    if f_id isa AnyIndexNew
        one(eltype(face_reference_id(q)))
    else
        face_reference_id(q)[f_id]
    end
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
    fn = face_nodes(mesh,d)
    if f_id isa AnyIndexNew
        T = eltype(eltype(fn))
        T[]
    else
        fn[f_id]
    end
end

function nodes(f::AbstractSpaceFaceNew)
    space = GT.space(f)
    f_id = GT.id(f)
    fn = face_nodes(space)
    if f_id isa AnyIndexNew
        T = eltype(eltype(fn))
        T[]
    else
        fn[f_id]
    end
end

function num_nodes(a::AbstractFaceNew)
    length(nodes(a))
end

function node_coordinates(a::AbstractMeshFaceNew)
    mesh = GT.mesh(a) 
    node_x = node_coordinates(mesh)
    view(node_x,nodes(a))
end

function dofs(a::AbstractFaceNew)
    space = GT.mesh_space(a)
    face_dofs = GT.face_dofs(space)
    i = id(a)
    face_dofs[i]
end

function num_dofs(a::AbstractFaceNew)
    length(dofs(a))
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

function allocate_shape_funcions(f,face::AbstractSamplingFaceNew)
    D = num_dims(mesh_face(face))
    d = num_dims(quadrature_face(face))
    aligned = Val(d==D)
    allocate_shape_funcions(f,face,aligned)
end

function allocate_values(face::AbstractSamplingFaceNew)
    D = num_dims(mesh_face(face))
    d = num_dims(quadrature_face(face))
    aligned = Val(d==D)
    allocate_values(face,aligned)
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
    face1 = replace_tabulators(f,face,rid_to_tab)
    allocate_shape_funcions(f,face1)
end

function allocate_shape_funcions(f,face::AbstractSamplingFaceNew,aligned::Val{true})
    space = GT.mesh_space(face)
    point = prototype(each_point_new(face))
    dof = ANY_INDEX_NEW
    sphys = shape_function(f,point,dof)
    ndofsr = max_num_reference_dofs(space)
    dof_sphys = zeros(typeof(sphys),ndofsr)
    replace_shape_function_allocation(f,face,dof_sphys)
end

function allocate_values(face::AbstractSamplingFaceNew,aligned::Val{true})
    space = GT.mesh_space(face)
    uh = GT.mesh_discrete_field(face)
    ndofsr = max_num_reference_dofs(space)
    T = eltype(GT.free_values(uh))
    alloc = zeros(T,ndofsr)
    replace_values_allocation(face,alloc)
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
    face1 = replace_tabulators(f,face,drid_Drid_ldface_perm_tab)
    allocate_shape_funcions(f,face1)
end

function allocate_shape_funcions(f,face::AbstractSamplingFaceNew,aligned::Val{false})
    space = GT.mesh_space(face)
    point = prototype(each_point_new(face))
    dof = ANY_INDEX_NEW
    sphys = shape_function(f,point,dof)
    ndofsr = max_num_reference_dofs(space)
    max_num_faces_around = 2 # TODO
    nfa = max_num_faces_around
    dof_sphys = [zeros(typeof(sphys),ndofsr) for _ in 1:nfa]
    replace_shape_function_allocation(f,face,dof_sphys)
end

function allocate_values(face::AbstractSamplingFaceNew,aligned::Val{false})
    space = GT.mesh_space(face)
    uh = GT.mesh_discrete_field(face)
    ndofsr = max_num_reference_dofs(space)
    T = eltype(GT.free_values(uh))
    max_num_faces_around = 2 # TODO
    nfa = max_num_faces_around
    alloc = [zeros(T,ndofsr) for _ in 1:nfa]
    replace_values_allocation(face,alloc)
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
    if GT.id(dface) isa AnyIndexNew
        perm = 1
    else
        perm = permutation_id(Dface,ldface)
    end
    if ldface === nothing
        ldid = 1
    else
        ldid = GT.id(ldface)
    end
    drid = reference_id(dface)
    Drid = reference_id(Dface)
    drid_Drid_ldface_perm_tab = tabulators
    drid_Drid_ldface_perm_tab[drid][Drid][ldid][perm]
end

function dofs(face::AbstractSamplingFaceNew)
    sf = space_face(face)
    dofs(sf)
end

function shape_functions_allocation(f,point::AbstractPointNew)
    sf = space_sampling_face(parent_face(point))
    shape_functions_allocation(f,sf)
end

function shape_functions(f,point::AbstractPointNew)
    map_shape_functions!(f,point)
    phys_funs = shape_functions_allocation(f,point)
    ndofs = num_dofs(parent_face(point))
    view(phys_funs,1:ndofs)
end

function shape_function(f,point::AbstractPointNew,dof)
    space = mesh_space(parent_face(point))
    sref = reference_shape_function(f,point,dof)
    point2 = mesh_sampling_point(point)
    sphys = map_shape_function(f,space,dof,point2,sref)
end

function map_shape_functions!(f,point::AbstractPointNew)
    space = mesh_space(parent_face(point))
    point2 = mesh_sampling_point(point)
    ref_funs = reference_shape_functions(f,point)
    phys_funs = shape_functions_allocation(f,point)
    ndofs = length(ref_funs)
    dof = 0
    while dof < ndofs
        dof += 1
        sref = ref_funs[dof]
        sphys = map_shape_function(f,space,dof,point2,sref)
        phys_funs[dof] = sphys
    end
    nothing
end

function reference_shape_functions(f,point::AbstractPointNew)
    tabulator = GT.tabulator(f,parent_face(point))
    p_id = id(point)
    if p_id isa AnyIndexNew
        ref_funs = view(tabulator,:,1)
    else
        ref_funs = view(tabulator,:,p_id)
    end
    ref_funs
end

function reference_shape_function(f,point::AbstractPointNew,dof)
    tabulator = GT.tabulator(f,parent_face(point))
    p_id = id(point)
    if dof isa AnyIndexNew
        zero(eltype(tabulator))
    else
        tabulator[dof,p_id]
    end
end

function values(a::AbstractFaceNew)
    fill_values!(a)
    v = values_allocation(a)
    ndofs = num_dofs(a)
    view(v,1:ndofs)
end

function fill_values!(a::AbstractFaceNew)
    v = values_allocation(a)
    uh = mesh_discrete_field(a)
    fv = free_values(uh)
    dv = dirichlet_values(uh)
    dofs = GT.dofs(a)
    for i in 1:length(dofs)
        g = dofs[i]
        if g < 0
            iv = dv[-g]
        else
            iv = fv[g]
        end
        v[i] = iv
    end
    nothing
end

function field(f,a::AbstractPointNew)
    s = shape_functions(f,a)
    pf = parent_face(a)
    x = values(pf)
    n = num_dofs(pf)
    zi = zero(eltype(x))*zero(eltype(s))
    z = zero(zi+zi)
    sum(i->x[i]*s[i],1:n;init=z)
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

function coordinate(::typeof(ForwardDiff.jacobian),point::AbstractPointNew)
    face = parent_face(point)
    mesh = GT.mesh(face)
    a = mesh_sampling_point(point)
    i_s = reference_shape_functions(ForwardDiff.gradient,a)
    node_x = node_coordinates(mesh)
    i_node = GT.nodes(mesh_face(face))
    n = length(i_node)
    x0 = zero(eltype(node_x))
    s0 = zero(eltype(i_s))
    o0 = outer(x0,s0)
    J0 = zero(o0+o0)
    J = sum_jacobian(node_x,i_node,i_s,n,J0)
    #sum(i->outer(x[i],s[i]),1:n)
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

function map_shape_function(::typeof(GT.value),space,dof,mesh_face,sref)
    sref
end

@inline function map_shape_function(::typeof(ForwardDiff.gradient),space,dof,mesh_face,sref)
    J = jacobian(mesh_face)
    sphys = transpose(J)\sref
end

function map_shape_function(::typeof(ForwardDiff.jacobian),space,dof,mesh_face,sref)
    J = jacobian(mesh_face)
    sphys = sref/J
end


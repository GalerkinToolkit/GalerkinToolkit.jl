abstract type AbstractFaceNew <: AbstractType end
abstract type AbstractPointNew <: AbstractType end
abstract type AbstractTabulation <: AbstractType end

function each_face_new(
    mesh::AbstractMesh, valD, quadrature::AbstractQuadrature)
    state = tabulate(mesh,valD,quadrature)
    each_face_new(state)
end

function each_face_new(quadrature::AbstractQuadrature)
    domain = GT.domain(quadrature)
    mesh = GT.mesh(domain)
    num_dims = GT.num_dims(domain)
    each_face_new(mesh,Val(num_dims),quadrature)
end

struct TabulatedMesh{A,B,C,D,E,F,G,H,I,J} <: AbstractTabulation
    mesh::A
    num_dims::Val{B}
    quadrature::C
    reference_shape_functions_value::D
    reference_shape_functions_gradient::E
    reference_unit_normals::F
    face_nodes::G
    node_coordinates::H
    reference_weights::I
    reference_coordinates::J
end
tabulated_mesh(a::TabulatedMesh) = a
reference_shape_functions(::typeof(value),a::TabulatedMesh) = a.reference_shape_functions_value
reference_shape_functions(::typeof(gradient),a::TabulatedMesh) = a.reference_shape_functions_gradient
face_nodes(a::TabulatedMesh) = a.face_nodes
node_coordinates(a::TabulatedMesh) = a.node_coordinates
reference_weights(a::TabulatedMesh) = a.reference_weights
quadrature(a::TabulatedMesh) = a.quadrature

function at_face_id(state::TabulatedMesh,face_id)
    TabulatedFace(state,face_id,nothing)
end

function at_point_id(state::TabulatedMesh,face,i)
    point = PointNew(face,i,nothing)
    jacobian = GT.jacobian(point)
    workspace = MeshPointWorkspace(jacobian)
    PointNew(face,i,workspace)
end

struct MeshPointWorkspace{A} <: AbstractType
    jacobian::A
end
jacobian(a::MeshPointWorkspace) = a.jacobian

function tabulate(mesh::AbstractMesh, valD, quadrature::AbstractQuadrature)
    D = val_parameter(valD)
    num_dims = Val(D)
    space = mesh_space(mesh,num_dims)
    reference_shape_functions_value = reference_shape_functions(value,space,quadrature)
    reference_shape_functions_gradient = reference_shape_functions(gradient,space,quadrature)
    #TODO a single reference element for the moment
    reference_unit_normals = first(map(normals,reference_spaces(space)))
    face_nodes = GT.face_nodes(mesh,num_dims)
    node_coordinates = GT.node_coordinates(mesh)
    reference_weights = first(map(weights,reference_quadratures(quadrature)))
    reference_coordinates = first(map(coordinates,reference_quadratures(quadrature)))
    TabulatedMesh(
                  mesh,
                  num_dims,
                  quadrature,
                  reference_shape_functions_value,
                  reference_shape_functions_gradient,
                  reference_unit_normals,
                  face_nodes,
                  node_coordinates,
                  reference_weights,
                  reference_coordinates,)
end

struct TabulatedSpace{A,B,C,D,E,F,G,H,I} <: AbstractTabulation
    space::A
    tabulated_mesh::B
    reference_shape_functions_value::C
    reference_shape_functions_gradient::D
    reference_shape_functions_jacobian::E
    shape_functions_value::F
    shape_functions_gradient::G
    shape_functions_jacobian::H
    face_dofs::I
end
quadrature(a::TabulatedSpace) = quadrature(tabulated_mesh(a))
tabulated_mesh(a::TabulatedSpace) = a.tabulated_mesh
tabulated_space(a::TabulatedSpace) = a
tabulated_field(a::TabulatedSpace) = a
reference_shape_functions(::typeof(value),a::TabulatedSpace) = a.reference_shape_functions_value
reference_shape_functions(::typeof(gradient),a::TabulatedSpace) = a.reference_shape_functions_gradient
reference_shape_functions(::typeof(jacobian),a::TabulatedSpace) = a.reference_shape_functions_jacobian
shape_functions(::typeof(value),a::TabulatedSpace) = a.shape_functions_value
shape_functions(::typeof(gradient),a::TabulatedSpace) = a.shape_functions_gradient
shape_functions(::typeof(jacobian),a::TabulatedSpace) = a.shape_functions_jacobian
face_dofs(a::TabulatedSpace) = a.face_dofs

function tabulate(space::AbstractSpace,quadrature::AbstractQuadrature;tabulate=())
    mesh = GT.mesh(space)
    D = GT.num_dims(mesh)
    tabulated_mesh = GT.tabulate(mesh,Val(D),quadrature)
    if GT.value in tabulate
        reference_shape_functions_value, shape_functions_value = allocate_shape_funcions(value,space,tabulated_mesh)
    else
        reference_shape_functions_value = nothing
        shape_functions_value = nothing
    end
    if ForwardDiff.gradient in tabulate
        reference_shape_functions_gradient, shape_functions_gradient = allocate_shape_funcions(gradient,space,tabulated_mesh)
    else
        reference_shape_functions_gradient = nothing
        shape_functions_gradient = nothing
    end
    if ForwardDiff.jacobian in tabulate
        reference_shape_functions_jacobian, shape_functions_jacobian = allocate_shape_funcions(jacobian,space,tabulated_mesh)
    else
        reference_shape_functions_jacobian = nothing
        shape_functions_jacobian = nothing
    end
    face_dofs = GT.face_dofs(space)
    a = TabulatedSpace(
                   space,
                   tabulated_mesh,
                   reference_shape_functions_value,
                   reference_shape_functions_gradient,
                   reference_shape_functions_jacobian,
                   shape_functions_value,
                   shape_functions_gradient,
                   shape_functions_jacobian,
                   face_dofs,
                  )
end

function each_face_new(space::AbstractSpace, quadrature::AbstractQuadrature;tabulate=())
    state = GT.tabulate(space,quadrature;tabulate)
    each_face_new(state)
end

function at_face_id(state::TabulatedSpace,i)
    mesh_face = at_face_id(tabulated_mesh(state),i)
    at_mesh_face(state,mesh_face)
end

function at_mesh_face(state::TabulatedSpace,mesh_face)
    TabulatedFace(state,id(mesh_face),mesh_face)
end

function at_point_id(state::TabulatedSpace,face,i)
    mesh_point = at_point_id(workspace(face),i)
    at_mesh_point(state,face,mesh_point)
end

function at_mesh_point(state::TabulatedSpace,face,mesh_point)
    point = PointNew(face,id(mesh_point),mesh_point)
    if state.shape_functions_value !== nothing
        map_shape_functions!(value,point)
    end
    if state.shape_functions_gradient !== nothing
        map_shape_functions!(gradient,point)
    end
    if state.shape_functions_jacobian !== nothing
        map_shape_functions!(jacobian,point)
    end
    point
end

struct TabulatedField{A,B,C,D} <: AbstractTabulation
    space::A
    tabulated_space::B
    free_values::C
    dirichlet_values::D
end
quadrature(a::TabulatedField) = quadrature(tabulated_mesh(a))
tabulated_mesh(a::TabulatedField) = tabulated_mesh(tabulated_space(a))
tabulated_space(a::TabulatedField) = a.tabulated_space
tabulated_field(a::TabulatedField) = a
reference_shape_functions(f,a::TabulatedField) = reference_shape_functions(f,tabulated_space(a))
shape_functions(f,a::TabulatedField) = shape_functions(f,tabulated_space(a))
face_dofs(a::TabulatedField) = face_dofs(tabulated_space(a))
free_values(a::TabulatedField) = a.free_values
dirichlet_values(a::TabulatedField) = a.dirichlet_values

function tabulate(uh::DiscreteField,quadrature::AbstractQuadrature;tabulate=())
    space = GT.space(uh)
    tabulated_space = GT.tabulate(space,quadrature;tabulate)
    free_values = GT.free_values(uh)
    dirichlet_values = GT.dirichlet_values(uh)
    a = TabulatedField(
                   space,
                   tabulated_space,
                   free_values,
                   dirichlet_values,
                  )
end

function each_face_new(uh::DiscreteField, quadrature::AbstractQuadrature;tabulate=())
    state = GT.tabulate(uh,quadrature;tabulate)
    each_face_new(state)
end

function at_face_id(state::TabulatedField,id)
    mesh_face = at_face_id(tabulated_mesh(state),id)
    at_mesh_face(state,mesh_face)
end

function at_mesh_face(state::TabulatedField,mesh_face)
    space_face = at_mesh_face(tabulated_space(state),mesh_face)
    TabulatedFace(state,id(mesh_face),space_face)
end

function at_point_id(state::TabulatedField,face,i)
    mesh_point = at_point_id(workspace(workspace(face)),i)
    at_mesh_point(state,face,mesh_point)
end

function at_mesh_point(state::TabulatedField,face,mesh_point)
    space_point = at_mesh_point(tabulated_space(state),workspace(face),mesh_point)
    point = PointNew(face,id(mesh_point),space_point)
end

struct Each{A,B,C} <: AbstractType
    update::A
    state::B
    iterator::C
end
Base.length(iter::Each) = length(iter.iterator)
Base.isdone(iter::Each,id) = id > length(iter)
function Base.getindex(iter::Each,i::Integer)
    item = iter.iterator[i]
    iter.update(iter.state,item)
end
function Base.getindex(iter::Each,mesh_point::AbstractPointNew)
    at_mesh_point(iter.state,mesh_point)
end
function Base.getindex(iter::Each,mesh_face::AbstractFaceNew)
    at_mesh_face(iter.state,mesh_face)
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
state(x::Each) = x.state
iterator(x::Each) = x.iterator

struct TabulatedFace{A,B,C} <: AbstractFaceNew
    parent::A
    id::B
    workspace::C
end
parent(a::TabulatedFace) = a.parent
id(a::TabulatedFace) = a.id
tabulated_mesh(a::TabulatedFace) = tabulated_mesh(parent(a))
tabulated_space(a::TabulatedFace) = tabulated_space(parent(a))
tabulated_field(a::TabulatedFace) = tabulated_field(parent(a))
workspace(a::TabulatedFace) = a.workspace
at_face_id(a,id) = TabulatedFace(a,id,nothing)


function each_face_new(state::AbstractTabulation)
    quadrature = GT.quadrature(state)
    ids = faces(domain(quadrature))
    Each(at_face_id,state,ids)
end

struct PointNew{A,B,C} <: AbstractPointNew
    parent::A
    id::B
    workspace::C
end
parent(a::PointNew) = a.parent
id(a::PointNew) = a.id
workspace(a::PointNew) = a.workspace
at_point_id(a,id) = at_point_id(parent(a),a,id)
at_point_id(s,a,id) = PointNew(a,id)

function each_point_new(face::TabulatedFace)
    Each(at_point_id,face,1:num_points(face))
end

function at_mesh_point(state::AbstractTabulation,mesh_point::AbstractPointNew)
    mesh_face = parent(mesh_point)
    face = at_mesh_face(state,mesh_face)
    at_mesh_point(state,face,mesh_point)
end

function Base.getindex(face::TabulatedFace,mesh_point::AbstractPointNew)
    state = parent(face)
    at_mesh_point(state,face,mesh_point)
end

function reference_shape_functions(f,space::AbstractSpace,quadrature::AbstractQuadrature)
    rid_point_x = map(coordinates,reference_quadratures(quadrature))
    # NB the TODOs below can be solved by introducing an extra nesting level
    # TODO assumes same reference elems for integration and for interpolation
    rid_to_tab = map(rid_point_x,reference_spaces(space)) do point_to_x, refface
        collect(permutedims(tabulator(refface)(f,point_to_x))) # TODO fix this globally
    end
    # For the moment we assume a single element type
    first(rid_to_tab)
end

function allocate_shape_funcions(f,space,tabulated_mesh)
    quadrature = GT.quadrature(tabulated_mesh)
    face = each_face_new(tabulated_mesh)[1]
    point = each_point_new(face)[1]
    i_ref_fun = reference_shape_functions(f,space,quadrature)
    ref_fun = i_ref_fun[1]
    fun = map_shape_function(f,space,1,point,ref_fun)
    ndofs = max_num_reference_dofs(space)
    i_ref_fun, zeros(typeof(fun),ndofs)
end

function num_points(face::AbstractFaceNew)
    alloc = tabulated_mesh(face)
    p_w = reference_weights(alloc)
    length(p_w)
end

function reference_weight(point::AbstractPointNew)
    alloc = tabulated_mesh(parent(point))
    p_w = reference_weights(alloc)
    p = id(point)
    p_w[p]
end

function weight(a::AbstractPointNew)
    w = reference_weight(a)
    J = jacobian(a)
    change_of_measure(J)*w
end

function coordinate(point::AbstractPointNew)
    alloc = tabulated_mesh(parent(point))
    i_p_ref_fun = reference_shape_functions(value,alloc)
    node_x = node_coordinates(alloc)
    f_i_node = face_nodes(alloc)
    f = id(parent(point))
    p = id(point)
    x0 = zero(eltype(node_x))
    s0 = zero(eltype(i_p_ref_fun))
    o0 = outer(x0,s0)
    J0 = zero(o0+o0)
    sum_coordinate_no_views(node_x,f_i_node,i_p_ref_fun,p,f,J0)
end

@noinline function sum_coordinate_no_views(node_x,f_i_node,i_p_ref_fun,p,f,J0)
    J = zero(J0)
    i = 0
    n = size(i_p_ref_fun,1)
    node_ptr = f_i_node.ptrs[f]
    while i < n
        i += 1
        node = f_i_node.data[node_ptr]
        node_ptr += 1
        x = node_x[node]
        ref_fun = i_p_ref_fun[i,p]
        J += x*ref_fun
    end
    J
    ##TODO sum leads to much faster than hand-written loop, but why?
    #the answer was the noinline above. But why?
    #J = sum(i->outer(node_x[i_node[i]],i_s[i]),1:n;init=zero(J0))
end

function jacobian(point::AbstractPointNew)
    jacobian(point,workspace(point))
end

function jacobian(point::AbstractPointNew,workspace::Nothing)
    alloc = tabulated_mesh(parent(point))
    i_p_ref_fun = reference_shape_functions(gradient,alloc)
    node_x = node_coordinates(alloc)
    f_i_node = face_nodes(alloc)
    f = id(parent(point))
    p = id(point)
    x0 = zero(eltype(node_x))
    s0 = zero(eltype(i_p_ref_fun))
    o0 = outer(x0,s0)
    J0 = zero(o0+o0)
    sum_jacobian_no_views(node_x,f_i_node,i_p_ref_fun,p,f,J0)
end

function jacobian(point::AbstractPointNew,workspace)
    jacobian(workspace)
end

@noinline function sum_jacobian_no_views(node_x,f_i_node,i_p_ref_fun,p,f,J0)
    J = zero(J0)
    i = 0
    n = size(i_p_ref_fun,1)
    node_ptr = f_i_node.ptrs[f]
    while i < n
        i += 1
        node = f_i_node.data[node_ptr]
        node_ptr += 1
        x = node_x[node]
        ref_fun = i_p_ref_fun[i,p]
        J += outer(x,ref_fun)
    end
    J
    ##TODO sum leads to much faster than hand-written loop, but why?
    #the answer was the noinline above. But why?
    #J = sum(i->outer(node_x[i_node[i]],i_s[i]),1:n;init=zero(J0))
end

function shape_functions(f,point::AbstractPointNew)
    alloc = GT.tabulated_space(parent(point))
    i_fun = shape_functions(f,alloc)
end

function map_shape_functions!(f,point::AbstractPointNew)
    alloc = GT.tabulated_space(parent(point))
    i_fun = shape_functions(f,alloc)
    i_p_ref_fun = reference_shape_functions(f,alloc)
    p = id(point)
    for i in 1:size(i_p_ref_fun,1)
        ref_fun = i_p_ref_fun[i,p]
        fun = map_shape_function(f,space,i,point,ref_fun)
        i_fun[i] = fun
    end
    nothing
end

function dofs(face::AbstractFaceNew)
    alloc = GT.tabulated_space(face)
    f_i_dof = face_dofs(alloc)
    f_i_dof[id(face)]
end

function num_dofs(face::AbstractFaceNew)
    alloc = GT.tabulated_space(face)
    f_i_dof = face_dofs(alloc)
    f = id(face)
    f_i_dof.ptrs[f+1] - f_i_dof.ptrs[f]
end

function field(f,point::AbstractPointNew)
    alloc = GT.tabulated_field(parent(point))
    i_fun = shape_functions(f,alloc)
    @assert i_fun !== nothing
    free_dof_val = free_values(alloc)
    diri_dof_val = dirichlet_values(alloc)
    f_i_dof = face_dofs(alloc)
    f = id(parent(point))
    x0 = zero(eltype(free_dof_val))
    s0 = zero(eltype(i_fun))
    o0 = x0*s0
    J0 = zero(o0+o0)
    sum_field_no_views(free_dof_val,diri_dof_val,f_i_dof,i_fun,f,J0)
end

@noinline function sum_field_no_views(free_dof_val,diri_dof_val,f_i_dof,i_fun,f,J0)
    J = zero(J0)
    i = 0
    n = length(i_fun)
    dof_ptr = f_i_dof.ptrs[f]
    while i < n
        i += 1
        dof = f_i_dof.data[dof_ptr]
        dof_ptr += 1
        if dof < 0
            val = diri_dof_val[-dof]
        else
            val = free_dof_val[dof]
        end
        fun = i_fun[i]
        J += val*fun
    end
    J
    ##TODO sum leads to much faster than hand-written loop, but why?
    #the answer was the noinline above. But why?
    #J = sum(i->outer(node_x[i_node[i]],i_s[i]),1:n;init=zero(J0))
end



function device_layout(fs::EachFace;kwargs...)
    accessor = device_layout(fs.accessor;kwargs...)
    EachFace(accessor)
end

function device_layout(mesh_face::MeshFace;tile_size=nothing)
    @assert tile_size === nothing "not implemented yet"
    mesh = GT.mesh(mesh_face)
    D = num_dims(mesh_face)
    face_nodes = GT.face_nodes_coalesced(mesh,D)
    node_coordinates = GT.node_coordinates(mesh)
    face_reference_id = GT.face_reference_id(mesh,D)
    tabulated_values = mesh_face.space_face.workspace.values
    tabulated_gradients = mesh_face.space_face.workspace.gradients
    quad = mesh_face.space_face.quadrature
    reference_weights = map(weights,GT.reference_quadratures(quad))
    Dface = 0
    point = 0
    MeshFaceWithDeviceLayout(
                     Dface,
                     point,
                     node_coordinates,
                     face_nodes,
                     face_reference_id,
                     reference_weights,
                     tabulated_values,
                     tabulated_gradients)
end

function face_nodes_coalesced(mesh,d)
    face_nodes = GT.face_nodes(mesh,d)
    nlnodes = max_num_reference_nodes(mesh,d)
    nfaces = num_faces(mesh,d)
    T = Int32
    face_nodes_coalesced = zeros(T,nfaces,nlnodes)
    for face in 1:nfaces
        nodes = face_nodes[face]
        for lnode in 1:length(nodes)
            node = nodes[lnode]
            face_nodes_coalesced[face,lnode] = node
        end
    end
    face_nodes_coalesced
end

struct MeshFaceWithDeviceLayout{A,B,C,D,E,F} <: AbstractMeshFace
    Dface::Int
    point::Int
    node_coordinates::A
    face_nodes::B
    face_reference_id::C
    reference_weights::D
    tabulated_values::E
    tabulated_gradients::F
end

function at_face(x::MeshFaceWithDeviceLayout,face)
    MeshFaceWithDeviceLayout(
                     face,
                     x.point,
                     x.node_coordinates,
                     x.face_nodes,
                     x.face_reference_id,
                     x.reference_weights,
                     x.tabulated_values,
                     x.tabulated_gradients)
end

function num_faces(x::MeshFaceWithDeviceLayout)
    length(x.face_reference_id)
end

function at_point(x::MeshFaceWithDeviceLayout,point)
    MeshFaceWithDeviceLayout(
                     x.Dface,
                     point,
                     x.node_coordinates,
                     x.face_nodes,
                     x.face_reference_id,
                     x.reference_weights,
                     x.tabulated_values,
                     x.tabulated_gradients)
end

function num_points(x::MeshFaceWithDeviceLayout)
    face = x.Dface
    rid = x.face_reference_id[face]
    values = x.tabulated_values[rid]
    size(values,2)
end

function coordinate(f::MeshFaceWithDeviceLayout)
    face = f.Dface
    rid = f.face_reference_id[face]
    point = f.point
    values = f.tabulated_values[rid]
    n = size(values,1)
    sum(1:n) do i
        s = values[i,point]
        node = f.face_nodes[face,i]
        x = f.node_coordinates[node]
        s*x
    end
end

function ForwardDiff.jacobian(f::MeshFaceWithDeviceLayout)
    face = f.Dface
    rid = f.face_reference_id[face]
    point = f.point
    values = f.tabulated_gradients[rid]
    n = size(values,1)
    sum(1:n) do i
        s = values[i,point]
        node = f.face_nodes[face,i]
        x = f.node_coordinates[node]
        outer(x,s)
    end
end

function weight(f::MeshFaceWithDeviceLayout)
    J = jacobian(f)
    face = f.Dface
    rid = f.face_reference_id[face]
    point = f.point
    w = f.reference_weights[rid][point]
    dJ =change_of_measure(J) 
    w*dJ
end

function Adapt.adapt_structure(to,x::MeshFaceWithDeviceLayout)
    MeshFaceWithDeviceLayout(
                     Adapt.adapt(to,x.Dface),
                     Adapt.adapt(to,x.point),
                     Adapt.adapt(to,x.node_coordinates),
                     Adapt.adapt(to,x.face_nodes),
                     Adapt.adapt(to,x.face_reference_id),
                     Adapt.adapt(to,x.reference_weights),
                     Adapt.adapt(to,x.tabulated_values),
                     Adapt.adapt(to,x.tabulated_gradients)
                    )
end


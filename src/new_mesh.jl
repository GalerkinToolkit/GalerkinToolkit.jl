
"""
    abstract type AbstractMesh

# Basic queries

- [`node_coordinates`](@ref)
- [`face_nodes`](@ref)
- [`face_reference_id`](@ref)
- [`reference_spaces`](@ref)
- [`periodic_nodes`](@ref)
- [`physical_faces`](@ref)
- [`outward_normals`](@ref)

# Basic constructors

- [`mesh`](@ref)
- [`mesh_from_gmsh`](@ref)
- [`cartesian_mesh`](@ref)

"""
abstract type AbstractMesh end

num_dims(m::AbstractMesh) = length(reference_spaces(m))-1
num_ambient_dims(m::AbstractMesh) = length(eltype(node_coordinates(m)))
options(m::AbstractMesh) = options(first(last(reference_spaces(m))))
num_faces(m::AbstractMesh,d) = length(face_reference_id(m,d))
num_nodes(fe::AbstractMesh) = length(node_coordinates(fe))

function label_faces_in_dim!(m::AbstractMesh,d;physical_name="__$d-FACES__")
    groups = physical_faces(m,d)
    if haskey(groups,physical_name)
        return m
    end
    Ti = int_type(options(m))
    faces = collect(Ti,1:num_faces(m,d))
    groups[physical_name] = faces
    m
end

function domain(mesh::AbstractMesh,d;is_reference_domain=Val(false),physical_name="__$d-FACES__")
    label_faces_in_dim!(mesh,d;physical_name)
    mesh_domain(;
        mesh,
        num_dims=Val(val_parameter(d)),
        physical_names=[physical_name],
        is_reference_domain)
end

function reference_domains(a::AbstractMesh,d)
    map(domain,reference_spaces(a,d))
end

function reference_domains(a::AbstractMesh)
    D = num_dims(a)
    ntuple(d->reference_domains(a,d-1),Val(D+1))
end

function remove_interior(mesh::AbstractMesh)
    GT.mesh(;
         node_coordinates = node_coordinates(mesh),
         face_nodes = face_nodes(mesh)[1:end-1],
         face_reference_id = face_reference_id(mesh)[1:end-1],
         reference_spaces = reference_spaces(mesh)[1:end-1],
         physical_faces = physical_faces(mesh)[1:end-1],
         outward_normals = outward_normals(mesh),
         is_cell_complex = Val(is_cell_complex(mesh)),
         periodic_nodes = periodic_nodes(mesh)
        )
end

struct Mesh{A} <: AbstractMesh
    contents::A
end

function mesh(;
        node_coordinates,
        face_nodes,
        face_reference_id,
        reference_spaces,
        periodic_nodes = default_periodic_nodes(reference_spaces),
        physical_faces = default_physical_faces(reference_spaces),
        outward_normals = nothing,
        is_cell_complex = Val(false),
        workspace = nothing,
    )
    contents = (;
                node_coordinates,
                face_nodes,
                face_reference_id,
                reference_spaces,
                periodic_nodes,
                physical_faces,
                outward_normals,
                is_cell_complex,
                workspace,
               )
    mesh = Mesh(contents)
end

function replace_workspace(mesh::Mesh,workspace)
    contents = (;
                node_coordinates=node_coordinates(mesh),
                face_nodes=face_nodes(mesh),
                face_reference_id=face_reference_id(mesh),
                reference_spaces=reference_spaces(mesh),
                periodic_nodes=periodic_nodes(mesh),
                physical_faces=physical_faces(mesh),
                outward_normals=outward_normals(mesh),
                is_cell_complex=Val(is_cell_complex(mesh)),
                workspace,
               )
    Mesh(contents)
end

node_coordinates(m::Mesh) = m.contents.node_coordinates
face_nodes(m::Mesh) = m.contents.face_nodes
face_nodes(m::Mesh,d) = m.contents.face_nodes[d+1]
face_reference_id(m::Mesh) = m.contents.face_reference_id
face_reference_id(m::Mesh,d) = m.contents.face_reference_id[d+1]
reference_spaces(m::Mesh) = m.contents.reference_spaces
reference_spaces(m::Mesh,d) = m.contents.reference_spaces[d+1]
physical_faces(m::Mesh) = m.contents.physical_faces
physical_faces(m::Mesh,d) = m.contents.physical_faces[d+1]
periodic_nodes(m::Mesh) = m.contents.periodic_nodes
is_cell_complex(m::Mesh) = val_parameter(m.contents.is_cell_complex)
workspace(m::Mesh) = m.contents.workspace

function default_physical_faces(reference_spaces)
    [ Dict{String,Vector{int_type(options(first(last(reference_spaces))))}}() for _ in 1:length(reference_spaces) ]
end

function default_periodic_nodes(reference_spaces)
    Ti = int_type(options(first(last(reference_spaces))))
    Ti[] => Ti[]
end

function outward_normals(m::Mesh)
    m.contents.outward_normals
end

"""
"""
function physical_names(mesh,d)
    groups = physical_faces(mesh,d)
    Set(keys(groups))
end

function physical_names(mesh;merge_dims=Val(false))
    D = num_dims(mesh)
    d_to_names = [ physical_names(mesh,d) for d in 0:D]
    if val_parameter(merge_dims) == false
        return d_to_names
    end
    reduce(union,d_to_names)
end

abstract type AbstractChain <: AbstractType end

num_dims(m::AbstractChain) = num_dims(domain(first(reference_spaces(m))))
num_ambient_dims(m::AbstractChain) = length(eltype(node_coordinates(m)))
options(m::AbstractChain) = options(first(reference_spaces(m)))
num_faces(m::AbstractChain) = length(face_reference_id(m))

"""
"""
function mesh(chain::AbstractChain)
    D = num_dims(chain)
    cell_nodes = face_nodes(chain)
    cell_reference_id = face_reference_id(chain)
    reference_cells = reference_spaces(chain)
    node_coords = node_coordinates(chain)
    face_to_nodes = Vector{typeof(cell_nodes)}(undef,D+1)
    face_to_refid = Vector{typeof(cell_reference_id)}(undef,D+1)
    for d in 0:D-1
        face_to_nodes[d+1] = Vector{Int}[]
        face_to_refid[d+1] = Int[]
    end
    face_to_nodes[end] = cell_nodes
    face_to_refid[end] = cell_reference_id
    ref_cell = first(reference_cells)
    ref_faces = reference_spaces(mesh_from_space(ref_cell))
    refid_to_refface = push(ref_faces[1:end-1],reference_cells)
    cell_groups = physical_faces(chain)
    groups = [ typeof(cell_groups)() for d in 0:D]
    groups[end] = cell_groups
    pnodes = periodic_nodes(chain)
    onormals = outward_normals(chain)
    GT.mesh(;
      node_coordinates = node_coords,
      face_nodes = face_to_nodes,
      face_reference_id = face_to_refid,
      reference_spaces = refid_to_refface,
      periodic_nodes = pnodes,
      physical_faces = groups,
      outward_normals = onormals)
end

struct Chain{A} <: AbstractChain
    contents::A
end

function chain(;
        node_coordinates,
        face_nodes,
        face_reference_id,
        reference_spaces,
        periodic_nodes = default_periodic_nodes((reference_spaces,)),
        physical_faces = default_physical_faces((reference_spaces,))[end],
        outward_normals = nothing,
    )
    contents = (;
                node_coordinates,
                face_nodes,
                face_reference_id,
                reference_spaces,
                periodic_nodes,
                physical_faces,
                outward_normals,
               )
    Chain(contents)
end

node_coordinates(m::Chain) = m.contents.node_coordinates
face_nodes(m::Chain) = m.contents.face_nodes
face_reference_id(m::Chain) = m.contents.face_reference_id
reference_spaces(m::Chain) = m.contents.reference_spaces
physical_faces(m::Chain) = m.contents.physical_faces
periodic_nodes(m::Chain) = m.contents.periodic_nodes
outward_normals(m::Chain) = m.contents.outward_normals

function chain(mesh::AbstractMesh,D=Val(num_dims(mesh)))
    d = val_parameter(D)
    chain(;
          node_coordinates=node_coordinates(mesh),
          face_nodes=face_nodes(mesh,d),
          face_reference_id=face_reference_id(mesh,d),
          reference_spaces=reference_spaces(mesh,d),
          periodic_nodes=periodic_nodes(mesh),
          physical_faces=physical_faces(mesh,d),
          outward_normals=outward_normals(mesh),
         )
end

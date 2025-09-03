

domain(a::AbstractDomain) = a

is_n_cube(geo::AbstractDomain) = false
is_simplex(geo::AbstractDomain) = false
is_axis_aligned(geo::AbstractDomain) = false
is_unitary(geo::AbstractDomain) = false
num_ambient_dims(a::AbstractDomain) = num_ambient_dims(mesh(a))
is_physical_domain(a::AbstractDomain) = ! is_reference_domain(a)
options(a::AbstractDomain) = options(mesh(a))
num_faces(a::AbstractDomain) = length(faces(a))

function is_unit_n_cube(geo::AbstractDomain)
    is_n_cube(geo) && is_unitary(geo)
end

function is_unit_simplex(geo::AbstractDomain)
    is_simplex(geo) && is_unitary(geo)
end

function is_boundary(dom::AbstractDomain)
    face_around(dom) !== nothing && (num_dims(dom) + 1) == num_dims(mesh(dom))
end

function max_num_faces_around(interpolation_domain::AbstractDomain,integration_domain::AbstractDomain)
    D = num_dims(interpolation_domain)
    d = num_dims(integration_domain)
    face_around = GT.face_around(integration_domain)
    if d == D
        1
    elseif d+1==D && face_around !== nothing
        1
    else
        2
    end
end

num_faces(geo::AbstractFaceDomain) = 1
faces(geo::AbstractFaceDomain) = [1]
inverse_faces(geo::AbstractFaceDomain) = [1]
face_around(geo::AbstractFaceDomain) = nothing
geometries(geo::AbstractFaceDomain,d) = 1:num_faces(mesh(geo),d)
num_geometries(geo::AbstractFaceDomain,d) = length(geometries(geo,d))

function geometries(geo::AbstractFaceDomain)
    m = mesh(geo)
    D = num_dims(m)
    [ 1:num_faces(m,d) for d in 0:D]
end

function vertex_permutations(geo::AbstractFaceDomain)
    ## Admissible if the following map is admissible
    # phi_i(x) = sum_i x_perm[i] * fun_i(x)
    # This map sends vertex i to vertex perm[i]
    Ti = int_type(options(geo))
    D = num_dims(geo)
    if D == 0
        return [[Ti(1)]]
    end
    geo_mesh = mesh(geo)
    vertex_to_geo_nodes = face_nodes(geo_mesh,0)
    vertex_to_geo_node = map(first,vertex_to_geo_nodes)
    nvertices = length(vertex_to_geo_node)
    # TODO compute this more lazily
    # so that we never compute it for 3d 
    # since it is not needed
    if D > 2
        return [collect(Ti,1:nvertices)]
    end
    permutations = Combinatorics.permutations(1:nvertices)
    if is_simplex(geo)
        return collect(Vector{Ti},permutations)
    end
    admissible_permutations = Vector{Ti}[]
    ref_face = lagrange_space(geo,1)
    fun_mesh = complexify(ref_face)
    geo_node_coords = node_coordinates(geo_mesh)
    fun_node_coords = node_coordinates(fun_mesh)
    vertex_coords = geo_node_coords[vertex_to_geo_node]
    degree = 1
    quad = quadrature(geo,degree)
    q = coordinates(quad)
    Tx = eltype(vertex_coords)
    TJ = typeof(zero(Tx)*zero(Tx)')
    A = tabulator(ref_face)(ForwardDiff.gradient,q)
    refvol = compute_volume(vertex_coords,TJ,A)
    perm_vertex_coords = similar(vertex_coords)
    for permutation in permutations
        for (j,cj) in enumerate(permutation)
          perm_vertex_coords[j] = vertex_coords[cj]
        end
        vol2 = compute_volume(perm_vertex_coords,TJ,A)
        if (refvol + vol2) ≈ (2*refvol)
            push!(admissible_permutations,permutation)
        end
    end
    admissible_permutations
end

function compute_volume(vertex_coords,TJ,A)
    vol = zero(eltype(TJ))
    for iq in 1:size(A,1)
        J = zero(TJ)
        for fun_node in 1:size(A,2)
            vertex = fun_node # TODO we are assuming that the vertices and nodes match
            g = A[iq,fun_node]
            x = vertex_coords[vertex]
            J += g*x'
        end
        vol += abs(det(J))
    end
    vol
end

function topology(geom::AbstractFaceDomain)
    D = num_dims(geom)
    if D != 0
        boundary = topology(GT.remove_interior(mesh(geom)))
    else
        boundary = nothing
    end
    vertex_permutations = GT.vertex_permutations(geom)
    face_topology(;boundary,vertex_permutations)
end

complexify(geom::AbstractFaceDomain) = mesh(geom)

"""
    unit_n_cube(d)
    unit_n_cube(Val(d))

Return an object representing a unit `d`-cube.
"""
function unit_n_cube(d;options=GT.options())
    D = val_parameter(d)
    num_dims = Val{D}()
    UnitNCube(num_dims,options)
end

"""
    struct UnitNCube{...} <: AbstractFaceDomain

"""
struct UnitNCube{D,A} <: AbstractFaceDomain
    num_dims::Val{D}
    options::A
end

options(geo::UnitNCube) = geo.options
num_dims(geo::UnitNCube{D}) where D = D
is_n_cube(geo::UnitNCube) = true
is_simplex(geo::UnitNCube) = false
is_simplex(geo::UnitNCube{0}) = true
is_simplex(geo::UnitNCube{1}) = true
is_axis_aligned(geo::UnitNCube) = true
is_unitary(geo::UnitNCube) = true
is_reference_domain(geo::UnitNCube) = true
is_physical_domain(geo::UnitNCube) = true
reference_domain(geo::UnitNCube) = geo

function bounding_box(geo::UnitNCube)
    D = num_dims(geo)
    Tv = real_type(options(geo))
    p0 = ntuple(i->zero(real_type(options(geo))),Val(D)) |> SVector{D,Tv}
    p1 = ntuple(i->one(real_type(options(geo))),Val(D)) |> SVector{D,Tv}
    (p0,p1)
end

function mesh(geom::UnitNCube{0})
    Tv = real_type(options(geom))
    Ti = int_type(options(geom))
    Tr = reference_int_type(options(geom))
    node_coordinates = [SVector{0,Tv}()]
    face_nodes = [Ti[1]]
    face_reference_id = [Tr[1]]
    space = lagrange_space(geom,1)
    reference_spaces = ((space,),)
    mesh(;
         node_coordinates,
         face_nodes,
         face_reference_id,
         reference_spaces,
         is_cell_complex = Val(true),
        )
end

function mesh(geom::UnitNCube{1})
    Tv = real_type(options(geom))
    Ti = int_type(options(geom))
    Tr = reference_int_type(options(geom))
    geom0 = if is_unit_n_cube(geom)
        unit_n_cube(Val(0);options=options(geom))
    else
        unit_simplex(Val(0);options=options(geom))
    end
    geom1 = geom
    space0 = lagrange_space(geom0,1)
    space1 = lagrange_space(geom1,1)
    node_coordinates = SVector{1,Tv}[(0,),(1,)]
    face_nodes = Vector{Vector{Ti}}[[[1],[2]],[[1,2]]]
    face_reference_id = Vector{Tr}[[1,1],[1]]
    reference_spaces = ((space0,),(space1,))
    normals = SVector{1,Tv}[(-1,),(1,)]
    mesh(;
         node_coordinates,
         face_nodes,
         face_reference_id,
         reference_spaces,
         normals,
         is_cell_complex = Val(true),
        )
end

function mesh(geom::UnitNCube{2})
    Tv = real_type(options(geom))
    Ti = int_type(options(geom))
    Tr = reference_int_type(options(geom))
    geom0 = unit_n_cube(Val(0);options=options(geom))
    geom1 = unit_n_cube(Val(1);options=options(geom))
    geom2 = geom
    space0 = lagrange_space(geom0,1)
    space1 = lagrange_space(geom1,1)
    space2 = lagrange_space(geom2,1)
    node_coordinates = SVector{2,Tv}[(0,0),(1,0),(0,1),(1,1)]
    face_nodes = Vector{Vector{Ti}}[[[1],[2],[3],[4]],[[1,2],[3,4],[1,3],[2,4]],[[1,2,3,4]]]
    face_reference_id = Vector{Tr}[[1,1,1,1],[1,1,1,1],[1]]
    normals = SVector{2,Tv}[(0,-1),(0,1),(-1,0),(1,0)]
    reference_spaces = ((space0,),(space1,),(space2,))
    mesh(;
         node_coordinates,
         face_nodes,
         face_reference_id,
         reference_spaces,
         is_cell_complex = Val(true),
         normals
        )
end

function mesh(geom::UnitNCube{3})
    Tv = real_type(options(geom))
    Ti = int_type(options(geom))
    Tr = reference_int_type(options(geom))
    geom0 = unit_n_cube(Val(0);options=options(geom))
    geom1 = unit_n_cube(Val(1);options=options(geom))
    geom2 = unit_n_cube(Val(2);options=options(geom))
    geom3 = geom
    space0 = lagrange_space(geom0,1)
    space1 = lagrange_space(geom1,1)
    space2 = lagrange_space(geom2,1)
    space3 = lagrange_space(geom3,1)
    node_coordinates = SVector{3,Tv}[(0,0,0),(1,0,0),(0,1,0),(1,1,0),(0,0,1),(1,0,1),(0,1,1),(1,1,1)]
    face_nodes = [
                  Vector{Ti}[[1],[2],[3],[4],[5],[6],[7],[8]],
                  Vector{Ti}[[1,2],[3,4],[1,3],[2,4],[5,6],[7,8],[5,7],[6,8],[1,5],[3,7],[2,6],[4,8]],
                  Vector{Ti}[[1,2,3,4],[5,6,7,8],[1,2,5,6],[3,4,7,8],[1,3,5,7],[2,4,6,8]],
                  Vector{Ti}[[1,2,3,4,5,6,7,8]],
                 ]
    face_reference_id = [ones(Tr,8),ones(Tr,12),ones(Tr,6),ones(Tr,1)]
    normals = SVector{3,Tv}[(0,0,-1),(0,0,1),(0,-1,0),(0,1,0),(-1,0,0),(1,0,0)]
    reference_spaces = ((space0,),(space1,),(space2,),(space3,))
    mesh(;
         node_coordinates,
         face_nodes,
         face_reference_id,
         reference_spaces,
         is_cell_complex = Val(true),
         normals
        )
end

function simplexify(geo::UnitNCube)
    @assert is_unit_n_cube(geo)
    D = num_dims(geo)
    if D in (0,1)
        return GT.mesh(geo)
    end
    simplex = unit_simplex(Val(D))
    ref_cell = lagrange_space(simplex,1)
    node_coords = node_coordinates(GT.mesh(geo))
    cell_nodes = simplex_nodes(geo)
    ncells = length(cell_nodes)
    cell_reference_id = fill(Int8(1),ncells)
    reference_cells = (ref_cell,)
    chain = GT.chain(
        node_coordinates = node_coords,
        face_nodes = cell_nodes,
        face_reference_id = cell_reference_id,
        reference_spaces = reference_cells,
       )
    mesh = GT.mesh(chain)
    mesh_complex = complexify(mesh)
    groups = [Dict{String,Vector{Int32}}() for d in 0:D]
    for d in 0:(D-1)
        sface_to_nodes = face_nodes(mesh_complex,d)
        cface_to_nodes = face_nodes(GT.mesh(geo),d)
        nsfaces = length(sface_to_nodes)
        ncfaces = length(cface_to_nodes)
        sface_touched = fill(false,nsfaces)
        for cface in 1:ncfaces
            fill!(sface_touched,false)
            nodes_c = cface_to_nodes[cface]
            for sface in 1:nsfaces
                nodes_s = sface_to_nodes[sface]
                if all(map(n->n in nodes_c,nodes_s))
                    sface_touched[sface] = true
                end
            end
            sfaces_in_group = findall(sface_touched)
            group_name = "$d-face-$cface"
            groups[d+1][group_name] = sfaces_in_group
        end
    end
    groups[end]["interior"] = 1:(num_faces(mesh_complex,D))
    sface_is_boundary = fill(false,num_faces(mesh_complex,D-1))
    for (_,sfaces) in groups[end-1]
        sface_is_boundary[sfaces] .= true
    end
    groups[end-1]["boundary"] = findall(sface_is_boundary)
    group_faces(mesh_complex) .= groups
    mesh_complex
end

function simplex_nodes(geo::UnitNCube)
    D = num_dims(geo)
    if D == 0
        [[1]]
    elseif D == 1
        [[1,2]]
    elseif D == 2
        [[1,2,3],[4,3,2]] # The second is oriented for proper shadowing with Makie
    elseif D ==3
        [[1,2,3,7], [1,2,5,7], [2,3,4,7],
        [2,4,7,8], [2,5,6,7], [2,6,7,8]]
    else
        error("case not implemented")
    end
end

"""
    unit_simplex(d)
    unit_simplex(Val(d))

Return an object representing a unit simplex of dimension `d`.
"""
function unit_simplex(d;options=GT.options())
    D = val_parameter(d)
    num_dims = Val{D}()
    UnitSimplex(num_dims,options)
end

"""
    struct UnitSimplex{...} <: AbstractFaceDomain
"""
struct UnitSimplex{D,A} <: AbstractFaceDomain
    num_dims::Val{D}
    options::A
end

options(geo::UnitSimplex) = geo.options
num_dims(geo::UnitSimplex{D}) where D = D
is_simplex(geo::UnitSimplex) = true
is_n_cube(geo::UnitSimplex) = false
is_n_cube(geo::UnitSimplex{0}) = true
is_n_cube(geo::UnitSimplex{1}) = true
is_axis_aligned(geo::UnitSimplex) = true
is_unitary(geo::UnitSimplex) = true
bounding_box(geo::UnitSimplex) = bounding_box(unit_n_cube(Val(num_dims(geo)),options=options(geo)))
is_reference_domain(geo::UnitSimplex) = true
is_physical_domain(geo::UnitSimplex) = true
reference_domain(geo::UnitSimplex) = geo

function mesh(geom::UnitSimplex{0})
    cube = unit_n_cube(Val(0);options=options(geom))
    mesh(cube)
end

function mesh(geom::UnitSimplex{1})
    cube = unit_n_cube(Val(1);options=options(geom))
    mesh(cube)
end

function mesh(geom::UnitSimplex{2})
    Tv = real_type(options(geom))
    Ti = int_type(options(geom))
    Tr = reference_int_type(options(geom))
    geom0 = unit_simplex(Val(0);options=options(geom))
    geom1 = unit_simplex(Val(1);options=options(geom))
    geom2 = geom
    space0 = lagrange_space(geom0,1)
    space1 = lagrange_space(geom1,1)
    space2 = lagrange_space(geom2,1)
    node_coordinates = SVector{2,Tv}[(0,0),(1,0),(0,1)]
    face_nodes = Vector{Vector{Ti}}[[[1],[2],[3]],[[1,2],[1,3],[2,3]],[[1,2,3]]]
    face_reference_id = Vector{Tr}[[1,1,1],[1,1,1],[1]]
    reference_spaces = ((space0,),(space1,),(space2,))
    n1 = sqrt(2)/2
    normals = SVector{2,Tv}[(0,-1),(-1,0),(n1,n1)]
    mesh(;
         node_coordinates,
         face_nodes,
         face_reference_id,
         reference_spaces,
         normals,
         is_cell_complex = Val(true),
        )
end

function mesh(geom::UnitSimplex{3})
    Tv = real_type(options(geom))
    Ti = int_type(options(geom))
    Tr = reference_int_type(options(geom))
    geom0 = unit_simplex(Val(0);options=options(geom))
    geom1 = unit_simplex(Val(1);options=options(geom))
    geom2 = unit_simplex(Val(2);options=options(geom))
    geom3 = geom
    space0 = lagrange_space(geom0,1)
    space1 = lagrange_space(geom1,1)
    space2 = lagrange_space(geom2,1)
    space3 = lagrange_space(geom3,1)
    node_coordinates = SVector{3,Tv}[(0,0,0),(1,0,0),(0,1,0),(0,0,1)]
    face_nodes = [
                  Vector{Ti}[[1],[2],[3],[4]],
                  Vector{Ti}[[1,2],[1,3],[2,3],[1,4],[2,4],[3,4]],
                  Vector{Ti}[[1,2,3],[1,2,4],[1,3,4],[2,3,4]],
                  Vector{Ti}[[1,2,3,4]]
                 ]
    face_reference_id = [ones(Tr,4),ones(Tr,6),ones(Tr,4),ones(Tr,1)]
    n1 = sqrt(3)/3
    normals = SVector{3,Tv}[(0,0,-1),(0,-1,0),(-1,0,0),(n1,n1,n1)]
    reference_spaces = ((space0,),(space1,),(space2,),(space3,))
    mesh(;
         node_coordinates,
         face_nodes,
         face_reference_id,
         reference_spaces,
         is_cell_complex = Val(true),
         normals
        )
end

function simplexify(geo::UnitSimplex)
    mesh(geo)
end

"""
    abstract type AbstractMeshDomain{A} <: AbstractDomain{A} end
"""
abstract type AbstractMeshDomain{A} <: AbstractDomain{A} end

function face_reference_id(a::AbstractMeshDomain)
    d = num_dims(a)
    face_reference_id(mesh(a),d)
end

function reference_domains(a::AbstractMeshDomain)
    d = num_dims(a)
    reference_domains(mesh(a),d)
end

function Base.:(==)(a::AbstractMeshDomain,b::AbstractMeshDomain)
    flag = true
    flag = flag && (GT.mesh_id(a) == GT.mesh_id(b))
    flag = flag && (GT.group_names(a) == GT.group_names(b))
    flag = flag && (GT.num_dims(a) == GT.num_dims(b))
    flag = flag && (GT.is_reference_domain(a) == GT.is_reference_domain(b))
    flag
end

struct PhysicalDomain{A,B} <: AbstractMeshDomain{A}
    mesh::A
    contents::B
end
is_reference_domain(a::PhysicalDomain) = false

struct ReferenceDomain{A,B} <: AbstractMeshDomain{A}
    mesh::A
    contents::B
end
is_reference_domain(a::ReferenceDomain) = true

const MeshDomain = Union{PhysicalDomain{A,B},ReferenceDomain{A,B}} where {A,B} 

function mesh_domain(mesh;
    mesh_id = objectid(mesh),
    num_dims = Val(GT.num_dims(mesh)),
    group_names=GT.group_names(mesh,num_dims),
    is_reference_domain = Val(false),
    face_around=nothing,
    workspace=nothing,
    setup = Val(true),
    )

    contents = (;
                mesh_id,
                group_names,
                num_dims = Val(val_parameter(num_dims)),
                face_around,
                workspace,
               )
    domain = ReferenceDomain(mesh,contents)
    if val_parameter(setup)
        domain2 = setup_domain(domain)
    else
        domain2 = domain
    end
    if val_parameter(is_reference_domain)
        domain3 = domain2
    else
        domain3 = physical_domain(domain2)
    end
    domain3
end

function reference_domain(domain::PhysicalDomain)
    (;mesh, contents) = domain
    ReferenceDomain(mesh,contents)
end

function reference_domain(domain::ReferenceDomain)
    domain
end

function physical_domain(domain::ReferenceDomain)
    (;mesh, contents) = domain
    PhysicalDomain(mesh,contents)
end

function physical_domain(domain::PhysicalDomain)
    domain
end

mesh(a::MeshDomain) = a.mesh
mesh_id(a::MeshDomain) = a.contents.mesh_id
group_names(a::MeshDomain) = a.contents.group_names
face_around(a::MeshDomain) = a.contents.face_around
num_dims(a::MeshDomain) = GT.val_parameter(a.contents.num_dims)
workspace(a::MeshDomain) = a.contents.workspace

function setup_domain(domain::MeshDomain)
    if GT.workspace(domain) !== nothing
        return domain
    end
    faces = GT.faces(domain)
    inverse_faces = GT.inverse_faces(domain)
    workspace = (;faces,inverse_faces)
    replace_workspace(domain,workspace)
end

function PartitionedArrays.partition(pdomain::MeshDomain)
    if GT.workspace(pdomain) !== nothing
        return GT.workspace(pdomain).domain_partition
    end
    pmesh = GT.mesh(pdomain)
    p_mesh = partition(pmesh)
    domain_partition = map(p_mesh) do mesh
        mesh_domain(mesh;
                    num_dims = Val(GT.num_dims(pdomain)),
                    group_names=GT.group_names(pdomain),
                    is_reference_domain = Val(is_reference_domain(pdomain)),
                    face_around=face_around(pdomain),
                    setup = Val(false))
    end
end

function setup_domain(pdomain::MeshDomain{<:AbstractPMesh})
    if GT.workspace(pdomain) !== nothing
        return pdomain
    end
    p_domain = partition(pdomain)
    domain_partition = map(setup_domain,p_domain)
    workspace = (;domain_partition)
    replace_workspace(pdomain,workspace)
end

function replace_workspace(domain::MeshDomain,workspace)
    mesh = GT.mesh(domain)
    num_dims = GT.num_dims(domain)
    mesh_id = GT.mesh_id(domain)
    group_names = GT.group_names(domain)
    is_reference_domain = GT.is_reference_domain(domain)
    face_around = GT.face_around(domain)
    GT.mesh_domain(mesh;num_dims,mesh_id,group_names,is_reference_domain,face_around,workspace)
end

function faces(domain::MeshDomain)
    if workspace(domain) !== nothing
        return workspace(domain).faces
    end
    Ti = int_type(options(domain))
    mesh = domain |> GT.mesh
    D = GT.num_dims(domain)
    Dface_to_tag = zeros(Ti,GT.num_faces(mesh,D))
    tag_to_name = GT.group_names(domain)
    fill!(Dface_to_tag,zero(eltype(Dface_to_tag)))
    group_faces = GT.group_faces(mesh,D)
    for (tag,name) in enumerate(tag_to_name)
        for (name2,faces) in group_faces
            if name != name2
                continue
            end
            Dface_to_tag[faces] .= tag
        end
    end
    physical_Dfaces = findall(i->i!=0,Dface_to_tag)
    Ti.(physical_Dfaces)
end

function inverse_faces(domain::MeshDomain)
    if workspace(domain) !== nothing
        return workspace(domain).inverse_faces
    end
    Ti = int_type(options(domain))
    d = num_dims(domain)
    ndfaces = num_faces(mesh(domain),d)
    dface_to_face = zeros(Ti,ndfaces)
    face_to_dface = faces(domain)
    dface_to_face[face_to_dface] = 1:length(face_to_dface)
    dface_to_face
end


# TODO there are functions that depending on the input
# return objects on different interfaces
# topology
# boundary
# reference_faces
#

abstract type AbstractType end
function Base.show(io::IO,data::GT.AbstractType)
    print(io,"GalerkinToolkit.$(nameof(typeof(data)))(…)")
end

function push(a::AbstractVector,x)
    b = copy(a)
    push!(b,x)
    b
end

function push(a::Tuple,x)
    (a...,x)
end


val_parameter(a) = a
val_parameter(::Val{a}) where a = a

"""
    num_dims(a)

Return the number of space dimensions of `a`. Defaults to `a.num_dims`.
"""
num_dims(a) = val_parameter(a.num_dims)
"""
"""
node_coordinates(a) = a.node_coordinates
"""
"""
reference_faces(a) = a.reference_faces
"""
"""
face_nodes(a) = a.face_nodes
"""
"""
face_incidence(a) = a.face_incidence
"""
"""
face_reference_id(a) = a.face_reference_id
"""
"""
physical_faces(a) = a.physical_faces
has_physical_faces(a) = hasproperty(a,:physical_faces) && a.physical_faces !== nothing
"""
"""
periodic_nodes(a) = a.periodic_nodes
# The uncommented code would also ensure that an empty periodic node list 
# corresponds to mesh with non-periodic BCs
has_periodic_nodes(a) = hasproperty(a,:periodic_nodes) && 
    a.periodic_nodes !== nothing # && length(a.periodic_nodes.first) != 0
"""
"""
geometry(a) = a.geometry
"""
"""
boundary(a) = a.boundary
"""
"""
is_n_cube(a) = hasproperty(a,:is_n_cube) ? val_parameter(a.is_n_cube) : false
"""
"""
is_simplex(a) = hasproperty(a,:is_simplex) ? val_parameter(a.is_simplex) : false

"""
"""
is_axis_aligned(a) = a.is_axis_aligned
"""
"""
bounding_box(a) = a.bounding_box
"""
"""
vertex_permutations(a) = a.vertex_permutations
face_own_dofs(a) = a.face_own_dofs
face_own_dof_permutations(a) = a.face_own_dof_permutations
node_to_dofs(a) = a.node_to_dofs
dof_to_node(a) = a.dof_to_node
dof_to_index(a) = a.dof_to_index
num_dofs(a) = a.num_dofs
"""
"""
coordinates(a) = a.coordinates
"""
"""
weights(a) = a.weights
order_per_dir(a) = a.order_per_dir
"""
"""
monomial_exponents(a) = a.monomial_exponents
"""
"""
primal_basis(a) = a.primal_basis
"""
"""
dual_basis(a) = a.dual_basis
"""
"""
lib_to_user_nodes(a) = a.lib_to_user_nodes
"""
"""
interior_nodes(a) = a.interior_nodes
"""
"""
face_permutation_ids(a) = a.face_permutation_ids
face_permutation_ids(a,m,n) = face_permutation_ids(a)[m+1,n+1]
local_nodes(a) = a.local_nodes
local_node_colors(a) = a.local_node_colors
"""
"""
real_type(a) = a.real_type
"""
"""
int_type(a) = a.int_type

# TODO rename for unit_normals
"""
"""
outwards_normals(a) = a.outwards_normals

function repeat_per_dir(geo,a)
    D = num_dims(geo)
    ntuple(i->a,Val(D))
end
repeat_per_dir(geo,a::NTuple) = a

reference_faces(a,d) = reference_faces(a)[val_parameter(d)+1]
face_nodes(a,d) = face_nodes(a)[val_parameter(d)+1]
face_incidence(a,d1,d2) = face_incidence(a)[val_parameter(d1)+1,val_parameter(d2)+1]
function face_local_faces(topo,d,D)
    dface_to_Dfaces = JaggedArray(GT.face_incidence(topo,d,D))
    Dface_to_dfaces = GT.face_incidence(topo,D,d)
    dface_to_lfaces_data = similar(dface_to_Dfaces.data) 
    fill!(dface_to_lfaces_data,0)
    dface_to_lfaces_ptrs = dface_to_Dfaces.ptrs
    dface_to_lfaces = JaggedArray(dface_to_lfaces_data,dface_to_lfaces_ptrs)
    for (dface1,Dfaces) in enumerate(dface_to_Dfaces)
        for (lface1,Dface) in enumerate(Dfaces)
            dfaces = Dface_to_dfaces[Dface]
            for (lface2,dface2) in enumerate(dfaces)
                if dface1 == dface2
                    dface_to_lfaces[dface1][lface1] = lface2
                end
            end
        end
    end
    dface_to_lfaces
end
face_reference_id(a,d) = face_reference_id(a)[val_parameter(d)+1]
num_faces(a) = map(length,face_reference_id(a))
num_faces(a,d) = length(face_reference_id(a,d))
physical_faces(a,d) = physical_faces(a)[val_parameter(d)+1]
"""
"""
num_nodes(a) = length(node_coordinates(a))
num_ambient_dims(a) = length(eltype(node_coordinates(a)))
function face_offset(a)
    D = num_dims(a)
    offsets = zeros(Int,D+1)
    for d in 1:D
        offsets[d+1] = offsets[d] + num_faces(a,d-1)
    end
    offsets
end
function face_offset(a,d)
    face_offset(a)[d+1]
end
function face_range(a,d)
    o = face_offset(a,d)
    o .+ (1:num_faces(a,d))
end
function face_range(a)
    D = num_dims(a)
    map(d->face_range(a,d),0:D)
end
function face_dim(a,d)
    n = num_faces(a,d)
    fill(d,n)
end
function face_dim(a)
    D = num_dims(a)
    reduce(vcat,map(d->face_dim(a,d),0:D))
end

"""
    abstract type AbstractFaceGeometry

# Basic queries

- [`num_dims`](@ref)
- [`is_axis_aligned`](@ref)
- [`is_simplex`](@ref)
- [`is_n_cube`](@ref)
- [`real_type`](@ref)
- [`int_type`](@ref)
- [`is_unit_n_cube`](@ref)
- [`is_unit_simplex`](@ref)
- [`is_unitary`](@ref)
- [`bounding_box`](@ref)
- [`boundary`](@ref)
- [`vertex_permutations`](@ref)

# Basic constructors

- [`unit_simplex`](@ref)
- [`unit_n_cube`](@ref)

# Supertype hierarchy

    AbstractFaceGeometry <: GT.AbstractType

"""
abstract type AbstractFaceGeometry <: GT.AbstractType end

struct ExtrusionPolytope{D,Tv,Ti} <: AbstractFaceGeometry
    extrusion::NTuple{D,Bool}
    bounding_box::NTuple{2,SVector{D,Tv}}
    int_type::Type{Ti}
end

"""
"""
function unit_simplex(num_dims;real_type=Float64,int_type=Int)
    D = val_parameter(num_dims)
    extrusion = ntuple(i->false,Val(D))
    Tv = real_type
    p0 = ntuple(i->zero(real_type),Val(D)) |> SVector{D,Tv}
    p1 = ntuple(i->one(real_type),Val(D)) |> SVector{D,Tv}
    bounding_box = (p0,p1)
    ExtrusionPolytope(extrusion,bounding_box,int_type)
end

"""
"""
function unit_n_cube(num_dims;real_type=Float64,int_type=Int)
    D = val_parameter(num_dims)
    extrusion = ntuple(i->true,Val(D))
    Tv = real_type
    p0 = ntuple(i->zero(real_type),Val(D)) |> SVector{D,Tv}
    p1 = ntuple(i->one(real_type),Val(D)) |> SVector{D,Tv}
    bounding_box = (p0,p1)
    ExtrusionPolytope(extrusion,bounding_box,int_type)
end


num_dims(p::ExtrusionPolytope{D}) where D = D
is_axis_aligned(p::ExtrusionPolytope) = true
is_simplex(geom::ExtrusionPolytope) = all(i->i==false,geom.extrusion)
is_n_cube(geom::ExtrusionPolytope) = all(geom.extrusion)
real_type(p::ExtrusionPolytope{D,Tv}) where {D,Tv} = Tv
int_type(p::ExtrusionPolytope{D,Tv,Ti}) where {D,Tv,Ti} = Ti

"""
"""
function is_unit_n_cube(geo)
    is_n_cube(geo) && is_unitary(geo)
end

"""
"""
function is_unit_simplex(geo)
    is_simplex(geo) && is_unitary(geo)
end

"""
"""
function is_unitary(geom)
    ! is_axis_aligned(geom) && return false
    my_bounding_box = bounding_box(geom)
    all(i->i==0,first(my_bounding_box)) && all(i->i==1,last(my_bounding_box))
end

abstract type AbstractFiniteElement <: AbstractType end

"""
    abstract type AbstractMeshFace

# Basic queries

- [`geometry`](@ref)
- [`num_dims`](@ref)
- [`num_nodes`](@ref)
- [`node_coordinates`](@ref)
- [`monomial_exponents`](@ref)
- [`primal_basis`](@ref)
- [`dual_basis`](@ref)
- [`lib_to_user_nodes`](@ref)
- [`shape_functions`](@ref)
- [`tabulator`](@ref)
- [`boundary`](@ref)
- [`interior_nodes`](@ref)
- [`interior_node_permutations`](@ref)

# Basic constructors

- [`lagrange_mesh_face`](@ref)

"""
abstract type AbstractMeshFace <: AbstractFiniteElement end

num_dims(f::AbstractMeshFace) = num_dims(geometry(f))

num_dofs(f::AbstractMeshFace) = num_nodes(f)

abstract type AbstractLagrangeMeshFace <: AbstractMeshFace end

AutoHashEquals.@auto_hash_equals struct GenericLagrangeMeshFace{A,B,C} <: AbstractLagrangeMeshFace
    geometry::A
    order_per_dir::B
    space::Symbol
    lib_to_user_nodes::C
end

lagrange_mesh_face(args...) = GenericLagrangeMeshFace(args...)

"""
"""
function lagrange_mesh_face(geometry,order;
        space = default_space(geometry),
        lib_to_user_nodes = int_type(geometry)[])
    D = num_dims(geometry)
    order_per_dir = repeat_per_dir(geometry,order)
    lagrange_mesh_face(
               geometry,
               order_per_dir,
               space,
               lib_to_user_nodes)
end

order(fe::AbstractLagrangeMeshFace) = maximum(order_per_dir(fe),init=0)

function lib_to_user_nodes(fe::AbstractLagrangeMeshFace)
    if length(fe.lib_to_user_nodes) == 0
        nnodes = num_nodes(fe)
        Ti = int_type(geometry(fe))
        collect(Ti.(1:nnodes))
    else
        fe.lib_to_user_nodes
    end
end

function monomial_exponents(fe::AbstractLagrangeMeshFace)
    monomial_exponents_from_space(fe.space,fe.order_per_dir,fe.geometry |> int_type)
end

function node_coordinates(fe::AbstractLagrangeMeshFace)
    if order(fe) == 0 && boundary(geometry(fe)) !== nothing
        x = node_coordinates(boundary(geometry(fe)))
        return [ sum(x)/length(x) ]
    end
    @assert fe |> geometry |> is_unitary
    mexps = monomial_exponents(fe)
    lib_node_to_coords = node_coordinates_from_monomials_exponents(mexps,fe.order_per_dir,fe.geometry |> real_type)
    if length(fe.lib_to_user_nodes) == 0
        return lib_node_to_coords
    end
    user_node_to_coords = similar(lib_node_to_coords)
    user_node_to_coords[fe.lib_to_user_nodes] = lib_node_to_coords
    user_node_to_coords
end

function primal_basis(fe::AbstractLagrangeMeshFace)
    map(e->(x-> prod(x.^e)),fe|>monomial_exponents)
end

function dual_basis(fe::AbstractLagrangeMeshFace)
    node_coordinates_reffe = fe|>node_coordinates
    map(x->(f->f(x)),node_coordinates_reffe)
end

function default_space(geom)
    if is_simplex(geom)
        :P
    elseif is_n_cube(geom)
        :Q
    else
        error("Not implemented")
    end
end

function node_coordinates_from_monomials_exponents(monomial_exponents,order_per_dir,::Type{Tv}) where Tv
    if length(monomial_exponents) == 0
        return SVector{length(order_per_dir),Tv}[]
    end
    node_coordinates = map(monomial_exponents) do exponent
        map(exponent,order_per_dir) do e,order
            if order != 0
               Tv(e/order)
            else
               Tv(e)
            end
        end |> SVector{length(order_per_dir),Tv}
    end
end

function monomial_exponents_from_space(space,args...)
    if space == :Q
        monomial_exponents_from_filter((e,o)->true,args...)
    elseif space == :P
        monomial_exponents_from_filter((e,o)->sum(e)<=maximum(o,init=0),args...)
    else
        error("Case not implemented (yet)")
    end
end

function monomial_exponents_from_filter(f,order_per_dir,::Type{Ti}) where Ti
    terms_per_dir = Tuple(map(d->d+1,order_per_dir))
    D = length(terms_per_dir)
    cis = CartesianIndices(terms_per_dir)
    m = count(ci->f(SVector{length(terms_per_dir),Ti}(Tuple(ci) .- 1),order_per_dir),cis)
    result = zeros(SVector{D,Ti},m)
    li = 0
    for ci in cis
        t = SVector{D,Ti}(Tuple(ci) .- 1)
        if f(t,order_per_dir)
            li += 1
            result[li] = t
        end
    end
    result
end

inner(a,b) = sum(map(*,a,b))
value(f,x) = f(x)

"""
"""
function shape_functions(fe)
    primal = primal_basis(fe)
    dual = dual_basis(fe)
    primal_t = permutedims(primal)
    A = value.(dual,primal_t)
    B = A\I
    n = length(primal)
    map(1:n) do i
        Bi = view(B,:,i)
        x->begin
            primal_t_x = map(f->f(x),primal_t)
            (primal_t_x*Bi)[1,1]#TODO
        end
    end
end

"""
"""
function tabulator(fe)
    primal = primal_basis(fe)
    dual = dual_basis(fe)
    primal_t = permutedims(primal)
    A = value.(dual,primal_t)
    B = A\I
    (f,x) -> begin
        C = broadcast(f,primal_t,x)
        C*B
    end
end

"""
    abstract type AbstractMesh

# Basic queries

- [`node_coordinates`](@ref)
- [`face_nodes`](@ref)
- [`face_reference_id`](@ref)
- [`reference_faces`](@ref)
- [`periodic_nodes`](@ref)
- [`physical_faces`](@ref)
- [`physical_nodes`](@ref)
- [`outwards_normals`](@ref)

# Basic constructors

- [`mesh_from_arrays`](@ref)
- [`mesh_from_gmsh`](@ref)
- [`cartesian_mesh`](@ref)
- [`mesh_from_chain`](@ref)

"""
abstract type AbstractMesh <: GT.AbstractType end

# TODO rename to Mesh
struct GenericMesh{A,B,C,D,E,F,G} <: AbstractMesh
    node_coordinates::A
    face_nodes::B
    face_reference_id::C
    reference_faces::D
    periodic_nodes::E
    physical_faces::F
    outwards_normals::G
end

function mesh_from_arrays(args...)
    GenericMesh(args...)
end

function replace_node_coordinates(mesh::AbstractMesh,node_coordinates)
    GenericMesh(
                node_coordinates,
                mesh.face_nodes,
                mesh.face_reference_id,
                mesh.reference_faces,
                mesh.periodic_nodes,
                mesh.physical_faces,
                mesh.outwards_normals)
end

"""
"""
function mesh_from_arrays(
    node_coordinates,
    face_nodes,
    face_reference_id,
    reference_faces;
    periodic_nodes = eltype(eltype(face_reference_id))[] => eltype(eltype(face_reference_id))[],
    physical_faces = map(i->Dict{String,Vector{eltype(eltype(face_reference_id))}}(),face_reference_id),
    outwards_normals = nothing
    )
    GT.mesh_from_arrays(
            node_coordinates,
            face_nodes,
            face_reference_id,
            reference_faces,
            periodic_nodes,
            physical_faces,
            outwards_normals)
end

num_dims(mesh::AbstractMesh) = length(reference_faces(mesh))-1

const INVALID_ID = 0

function default_gmsh_options()
    [
     "General.Terminal"=>1,
     "Mesh.SaveAll"=>1,
     "Mesh.MedImportGroupsOfNodes"=>1
    ]
end

function with_gmsh(f;options=default_gmsh_options())
    gmsh.initialize()
    for (k,v) in options
        gmsh.option.setNumber(k,v)
    end
    try
        return f()
    finally
        gmsh.finalize()
    end
end

"""
"""
function mesh_from_gmsh(file;complexify=true,renumber=true,kwargs...)
    @assert ispath(file) "File not found: $(file)"
    with_gmsh(;kwargs...) do
        gmsh.open(file)
        renumber && gmsh.model.mesh.renumberNodes()
        renumber && gmsh.model.mesh.renumberElements()
        mesh_from_gmsh_module(;complexify)
    end
end

function mesh_from_gmsh_module(;complexify=true)
    entities = gmsh.model.getEntities()
    nodeTags, coord, parametricCoord = gmsh.model.mesh.getNodes()

    # find num_dims
    ddim = -1
    for e in entities
        ddim = max(ddim,e[1])
    end
    if ddim == -1
        error("No entities in the msh file.")
    end
    D = ddim

    # find embedded_dimension
    dtouched = [false,false,false]
    for node in nodeTags
        if !(coord[(node-1)*3+1] + 1 ≈ 1)
            dtouched[1] = true
        end
        if !(coord[(node-1)*3+2] + 1 ≈ 1)
            dtouched[2] = true
        end
        if !(coord[(node-1)*3+3] + 1 ≈ 1)
            dtouched[3] = true
        end
    end
    if dtouched[3]
        adim = 3
    elseif dtouched[2]
        adim = 2
    elseif dtouched[1]
        adim = 1
    else
        adim = 0
    end

    # Setup node coords
    nmin = minimum(nodeTags)
    nmax = maximum(nodeTags)
    nnodes = length(nodeTags)
    if !(nmax == nnodes && nmin == 1)
        error("Only consecutive node tags allowed.")
    end
    my_node_to_coords = zeros(SVector{adim,Float64},nnodes)
    m = zero(MVector{adim,Float64})
    for node in nodeTags
        for j in 1:adim
            k = (node-1)*3 + j
            xj = coord[k]
            m[j] = xj
        end
        my_node_to_coords[node] = m
    end

    # Setup face nodes
    offsets = zeros(Int32,D+1)
    my_face_nodes = Vector{JaggedArray{Int32,Int32}}(undef,D+1)
    for d in 0:D
        elemTypes, elemTags, nodeTags = gmsh.model.mesh.getElements(d)
        ndfaces = 0
        for t in 1:length(elemTypes)
            ndfaces += length(elemTags[t])
        end
        if ndfaces != 0
            nmin::Int = minimum( minimum, elemTags )
            nmax::Int = maximum( maximum, elemTags )
            if !( (nmax-nmin+1) == ndfaces)
                error("Only consecutive elem tags allowed.")
            end
            offsets[d+1] = nmin-1
        end
        ptrs = zeros(Int32,ndfaces+1)
        dface = 0
        for t in 1:length(elemTypes)
            elementName, dim, order, numNodes, nodeCoord =
            gmsh.model.mesh.getElementProperties(elemTypes[t])
            for e in 1:length(elemTags[t])
                dface += 1
                ptrs[dface+1] = numNodes
            end
        end
        length_to_ptrs!(ptrs)
        ndata = ptrs[end]-1
        data = zeros(Int32,ndata)
        dface = 1
        for t in 1:length(elemTypes)
            p = ptrs[dface]-Int32(1)
            for (i,node) in enumerate(nodeTags[t])
                data[p+i] = node
            end
            dface += length(elemTags[t])
        end
        my_face_nodes[d+1] = JaggedArray(data,ptrs)
    end

    # Setup face_reference_id
    my_face_reference_id = Vector{Vector{Int32}}(undef,D+1)
    for d in 0:D
        elemTypes, elemTags, nodeTags = gmsh.model.mesh.getElements(d)
        ndfaces = length(my_face_nodes[d+1])
        dface_to_refid = zeros(Int8,ndfaces)
        refid = 0
        dface = 0
        for t in 1:length(elemTypes)
            refid += 1
            for e in 1:length(elemTags[t])
                dface += 1
                dface_to_refid[dface] = refid
            end
        end
        my_face_reference_id[d+1] = dface_to_refid
    end

    # Setup reference faces
    my_reference_faces = ()
    for d in D:-1:0
        elemTypes, elemTags, nodeTags = gmsh.model.mesh.getElements(d)
        refdfaces = ()
        for t in 1:length(elemTypes)
            refface = reference_face_from_gmsh_eltype(elemTypes[t])
            refdfaces = (refdfaces...,refface)
        end
        if refdfaces == ()
            refdfaces = reference_faces(boundary(first(first(my_reference_faces))),d)
        end
        my_reference_faces = (refdfaces,my_reference_faces...)
    end

    ## Setup periodic nodes
    node_to_master_node = fill(Int32(INVALID_ID),nnodes)
    for (dim,tag) in entities
        tagMaster, nodeTags, nodeTagsMaster, = gmsh.model.mesh.getPeriodicNodes(dim,tag)
        for i in 1:length(nodeTags)
            node = nodeTags[i]
            master_node = nodeTagsMaster[i]
            node_to_master_node[node] = master_node
        end
    end
    pnode_to_node = Int32.(findall(i->i!=INVALID_ID,node_to_master_node))
    pnode_to_master = node_to_master_node[pnode_to_node]
    periodic_nodes = pnode_to_node => pnode_to_master

    # Setup physical groups
    my_groups = [ Dict{String,Vector{Int32}}() for d in 0:D]
    for d in 0:D
        offset = Int32(offsets[d+1])
        dimTags = gmsh.model.getPhysicalGroups(d)
        for (dim,tag) in dimTags
            @boundscheck @assert dim == d
            g_entities = gmsh.model.getEntitiesForPhysicalGroup(dim,tag)
            ndfaces_in_physical_group = 0
            for entity in g_entities
                elemTypes, elemTags, nodeTags = gmsh.model.mesh.getElements(dim,entity)
                for t in 1:length(elemTypes)
                    ndfaces_in_physical_group += length(elemTags[t])
                end
            end
            dfaces_in_physical_group = zeros(Int32,ndfaces_in_physical_group)
            ndfaces_in_physical_group = 0
            for entity in g_entities
                elemTypes, elemTags, nodeTags = gmsh.model.mesh.getElements(dim,entity)
                for t in 1:length(elemTypes)
                    for etag in elemTags[t]
                        ndfaces_in_physical_group += 1
                        dfaces_in_physical_group[ndfaces_in_physical_group] = Int32(etag)-offset
                    end
                end
            end
            groupname = gmsh.model.getPhysicalName(dim,tag)
            my_groups[d+1][groupname] = dfaces_in_physical_group
        end
    end
    mesh = GT.mesh_from_arrays(
            my_node_to_coords,
            my_face_nodes,
            my_face_reference_id,
            my_reference_faces;
            physical_faces = my_groups,
            periodic_nodes,)

    if complexify
        mesh = GT.complexify(mesh)
    end
    mesh
end

function reference_face_from_gmsh_eltype(eltype)
    if eltype == 1
        order = 1
        geom = unit_n_cube(Val(1))
        lib_to_gmsh = [1,2]
    elseif eltype == 2
        order = 1
        geom = unit_simplex(Val(2))
        lib_to_gmsh = [1,2,3]
    elseif eltype == 3
        order = 1
        geom = unit_n_cube(Val(2))
        lib_to_gmsh = [1,2,4,3]
    elseif eltype == 4
        order = 1
        geom = unit_simplex(Val(3))
        lib_to_gmsh = [1,2,3,4]
    elseif eltype == 5
        order = 1
        lib_to_gmsh = [1,2,4,3,5,6,8,7]
    elseif eltype == 15
        order = 1
        geom = unit_n_cube(Val(0))
        lib_to_gmsh = [1]
    elseif eltype == 8
        order = 2
        geom = unit_n_cube(Val(1))
        lib_to_gmsh = [1,3,2]
    elseif eltype == 9
        order = 2
        geom = unit_simplex(Val(2))
        lib_to_gmsh = [1,4,2,6,5,3]
    else
        en, = gmsh.model.mesh.getElementProperties(eltype)
        error("Unsupported element type. elemType: $eltype ($en)")
    end
    lagrange_mesh_face(geom,order;lib_to_user_nodes=lib_to_gmsh)
end

# WriteVTK prototype functions to be defined inside 'ext/GalerkinToolkitWriteVTKExt.jl' extension module.

function vtk_points end
function vtk_points! end
function vtk_cells end
function vtk_cells! end
function vtk_args end
function vtk_args! end
function vtk_physical_faces end
function vtk_physical_faces! end
function vtk_physical_nodes end
function vtk_physical_nodes! end
function vtk_mesh_cell end
function vtk_mesh_cell! end

opposite_faces(geom,d) = opposite_faces(geom)[d+1]

# TODO use the same convention than in Gridap
# allow the user to customize the boundary object ids
function boundary(geom::ExtrusionPolytope{0})
    nothing
end

function opposite_faces(geom::ExtrusionPolytope{0})
    [[1]]
end

function boundary(geom::ExtrusionPolytope{1})
    Tv = real_type(geom)
    Ti = int_type(geom)
    if is_unit_simplex(geom)
        vertex = unit_simplex(0;real_type=Tv,int_type=Ti)
    elseif is_unit_n_cube(geom)
        vertex = unit_n_cube(0;real_type=Tv,int_type=Ti)
    else
        error("Not implemented")
    end
    order = 1
    fe = lagrange_mesh_face(vertex,order)
    node_coordinates = SVector{1,Tv}[(0,),(1,)]
    face_nodes = [Vector{Ti}[[1],[2]]]
    face_reference_id = [Ti[1,1]]
    reference_faces = ([fe],)
    outwards_normals = SVector{1,Tv}[(-1,),(1,)]
    GT.mesh_from_arrays(
            node_coordinates,
            face_nodes,
            face_reference_id,
            reference_faces;
            outwards_normals
           )
end

function opposite_faces(geom::ExtrusionPolytope{1})
    [[2,1],[1]]
end

function boundary(geom::ExtrusionPolytope{2})
    Tv = real_type(geom)
    Ti = int_type(geom)
    if is_unit_simplex(geom)
        n1 = sqrt(2)/2
        g0 = unit_simplex(0;real_type=Tv,int_type=Ti)
        g1 = unit_simplex(1;real_type=Tv,int_type=Ti)
        node_coordinates = SVector{2,Tv}[(0,0),(1,0),(0,1)]
        face_nodes = [Vector{Ti}[[1],[2],[3]],Vector{Ti}[[1,2],[1,3],[2,3]]]
        face_reference_id = [Ti[1,1,1],Ti[1,1,1]]
        outwards_normals = SVector{2,Tv}[(0,-1),(-1,0),(n1,n1)]
    elseif is_unit_n_cube(geom)
        g0 = unit_n_cube(0;real_type=Tv,int_type=Ti)
        g1 = unit_n_cube(1;real_type=Tv,int_type=Ti)
        node_coordinates = SVector{2,Tv}[(0,0),(1,0),(0,1),(1,1)]
        face_nodes = [Vector{Ti}[[1],[2],[3],[4]],Vector{Ti}[[1,2],[3,4],[1,3],[2,4]]]
        face_reference_id = [Ti[1,1,1,1],Ti[1,1,1,1]]
        outwards_normals = SVector{2,Tv}[(0,-1),(0,1),(-1,0),(1,0)]
    else
        error("Not implemented")
    end
    order = 1
    fe0 = lagrange_mesh_face(g0,order)
    fe1 = lagrange_mesh_face(g1,order)
    reference_faces = ([fe0],[fe1])
    GT.mesh_from_arrays(
            node_coordinates,
            face_nodes,
            face_reference_id,
            reference_faces;
            outwards_normals
           )
end

function opposite_faces(geom::ExtrusionPolytope{2})
    @assert is_n_cube(geom)
    [[4,3,2,1],[2,1,4,3],[1]]
end

function boundary(geom::ExtrusionPolytope{3})
    Tv = real_type(geom)
    Ti = int_type(geom)
    if is_unit_simplex(geom)
        g0 = unit_simplex(0;real_type=Tv,int_type=Ti)
        g1 = unit_simplex(1;real_type=Tv,int_type=Ti)
        g2 = unit_simplex(2;real_type=Tv,int_type=Ti)
        node_coordinates = SVector{3,Tv}[(0,0,0),(1,0,0),(0,1,0),(0,0,1)]
        face_nodes = [
            Vector{Ti}[[1],[2],[3],[4]],
            Vector{Ti}[[1,2],[1,3],[2,3],[1,4],[2,4],[3,4]],
            Vector{Ti}[[1,2,3],[1,2,4],[1,3,4],[2,3,4]],
           ]
        face_reference_id = [ones(Ti,4),ones(Ti,6),ones(Ti,4)]
        n1 = sqrt(3)/3
        outwards_normals = SVector{3,Tv}[(0,0,-1),(0,-1,0),(-1,0,0),(n1,n1,n1)]
    elseif is_unit_n_cube(geom)
        g0 = unit_n_cube(0;real_type=Tv,int_type=Ti)
        g1 = unit_n_cube(1;real_type=Tv,int_type=Ti)
        g2 = unit_n_cube(2;real_type=Tv,int_type=Ti)
        node_coordinates = SVector{3,Tv}[(0,0,0),(1,0,0),(0,1,0),(1,1,0),(0,0,1),(1,0,1),(0,1,1),(1,1,1)]
        face_nodes = [
                      Vector{Ti}[[1],[2],[3],[4],[5],[6],[7],[8]],
                      Vector{Ti}[[1,2],[3,4],[1,3],[2,4],[5,6],[7,8],[5,7],[6,8],[1,5],[3,7],[2,6],[4,8]],
                      Vector{Ti}[[1,2,3,4],[5,6,7,8],[1,2,5,6],[3,4,7,8],[1,3,5,7],[2,4,6,8]],
           ]
        face_reference_id = [ones(Ti,8),ones(Ti,12),ones(Ti,6)]
        outwards_normals = SVector{3,Tv}[(0,0,-1),(0,0,1),(0,-1,0),(0,1,0),(-1,0,0),(1,0,0)]
    else
        error("Not implemented")
    end
    order = 1
    fe0 = lagrange_mesh_face(g0,order)
    fe1 = lagrange_mesh_face(g1,order)
    fe2 = lagrange_mesh_face(g2,order)
    reference_faces = ([fe0,],[fe1],[fe2])
    GT.mesh_from_arrays(
            node_coordinates,
            face_nodes,
            face_reference_id,
            reference_faces;
            outwards_normals
           )
end

function opposite_faces(geom::ExtrusionPolytope{3})
    @assert is_n_cube(geom)
    [[8,7,6,5,4,3,2,1],[6,5,8,7,2,1,4,3,12,11,10,9],[2,1,4,3,6,5],[1]]
end

function boundary(refface::AbstractMeshFace)
    boundary_from_mesh_face(refface)
end

function boundary_from_mesh_face(refface)
    geom = geometry(refface)
    node_coordinates_inter = node_coordinates(refface)
    order_inter = order(refface)
    D = num_dims(geom)
    if D == 0
        return nothing
    end
    mesh_geom = boundary(geom)
    node_coordinates_geom = node_coordinates(mesh_geom)
    face_nodes_inter = Vector{Vector{Vector{Int}}}(undef,D)
    node_coordinates_aux = map(xi->map(xii->round(Int,order_inter*xii),xi),node_coordinates_inter)
    ref_faces_geom = reference_faces(mesh_geom)
    ref_faces = map(ref_faces_geom) do ref_faces_geom_d
        map(r->lagrange_mesh_face(geometry(r),order_inter),ref_faces_geom_d)
    end
    face_ref_id_geom = face_reference_id(mesh_geom)
    for d in 0:(D-1)
        s_ref = map(ref_faces_geom[d+1],ref_faces[d+1]) do r_geom,r
            m = num_nodes(r)
            n = num_nodes(r_geom)
            x = node_coordinates(r)
            tabulator(r_geom)(value,x)
        end
        face_nodes_geom = face_nodes(mesh_geom,d)
        nfaces = length(face_ref_id_geom[d+1])
        face_nodes_inter_d = Vector{Vector{Int}}(undef,nfaces)
        for face in 1:nfaces
            ref_id_geom = face_ref_id_geom[d+1][face]
            s = s_ref[ref_id_geom]
            nodes_geom = face_nodes_geom[face]
            nnodes, nnodes_geom = size(s)
            x_mapped = map(1:nnodes) do i
                x = zero(eltype(node_coordinates_inter))
                for k in 1:nnodes_geom
                    x += node_coordinates_geom[nodes_geom[k]]*s[i,k]
                end
                map(xi->round(Int,order_inter*xi),x)
            end
            my_nodes = indexin(x_mapped,node_coordinates_aux)
            face_nodes_inter_d[face] = my_nodes
        end
        face_nodes_inter[d+1] = face_nodes_inter_d
    end
    GT.mesh_from_arrays(
        node_coordinates_inter,
        face_nodes_inter,
        face_ref_id_geom,
        ref_faces)
end

function face_nodes(fe::AbstractMeshFace,d)
    face_nodes_from_mesh_face(fe,d)
end

function face_nodes_from_mesh_face(fe,d)
    D = num_dims(fe)
    if d == D
        [collect(Int32,1:GT.num_nodes(fe))]
    else
        boundary = GT.boundary(fe)
        GT.face_nodes(boundary,d)
    end
end

function face_interior_nodes(fe::AbstractMeshFace,d)
    face_interior_nodes_from_mesh_face(fe,d)
end

function face_interior_nodes_from_mesh_face(fe,d)
    D = num_dims(fe)
    if  d == D
        [GT.interior_nodes(fe)]
    else
        boundary = GT.boundary(fe)
        dface_to_lnode_to_node = GT.face_nodes(boundary,d)
        dface_to_ftype = GT.face_reference_id(boundary,d)
        ftype_to_refdface = GT.reference_faces(boundary,d)
        ftype_to_lnodes = map(GT.interior_nodes,ftype_to_refdface)
        map(dface_to_ftype,dface_to_lnode_to_node) do ftype,lnode_to_node
            lnodes = ftype_to_lnodes[ftype]
            lnode_to_node[lnodes]
        end
    end
end

function num_interior_nodes(fe::AbstractMeshFace)
    length(interior_nodes(fe))
end

function face_interior_node_permutations(fe::AbstractMeshFace,d)
    face_interior_node_permutations_from_mesh_face(fe,d)
end

function face_interior_node_permutations_from_mesh_face(fe,d)
    D = num_dims(fe)
    if  d == D
        [[ collect(1:num_interior_nodes(fe)) ]]
    else
        boundary = GT.boundary(fe)
        dface_to_ftype = GT.face_reference_id(boundary,d)
        ftype_to_refdface = GT.reference_faces(boundary,d)
        ftype_to_perms = map(GT.interior_node_permutations,ftype_to_refdface)
        map(dface_to_ftype) do ftype
            perms = ftype_to_perms[ftype]
        end
    end
end

function interior_nodes(fe::AbstractMeshFace)
    interior_nodes_from_mesh_face(fe)
end

function interior_nodes_from_mesh_face(fe)
    nnodes = num_nodes(fe)
    D = num_dims(fe|>geometry)
    if D == 0
        return collect(1:nnodes)
    else
        node_is_touched = fill(true,nnodes)
        mesh = boundary(fe)
        for d in 0:(D-1)
            face_to_nodes = face_nodes(fe,d)
            for nodes in face_to_nodes
                node_is_touched[nodes] .= false
            end
        end
        return findall(node_is_touched)
    end
end

function vertex_permutations(geom::AbstractFaceGeometry)
    vertex_permutations_from_face_geometry(geom)
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

## Admissible if the following map is admissible
# phi_i(x) = sum_i x_perm[i] * fun_i(x)
# This map sends vertex i to vertex perm[i]
function vertex_permutations_from_face_geometry(geo)
    D = num_dims(geo)
    if D == 0
        return [[1]]
    end
    geo_mesh = boundary(geo)
    vertex_to_geo_nodes = face_nodes(geo_mesh,0)
    vertex_to_geo_node = map(first,vertex_to_geo_nodes)
    nvertices = length(vertex_to_geo_node)
    # TODO compute this more lazily
    # so that we never compute it for 3d 
    # since it is not needed
    if D > 2
        return [collect(1:nvertices)]
    end
    permutations = Combinatorics.permutations(1:nvertices)
    if is_simplex(geo)
        return collect(permutations)
    end
    admissible_permutations = Vector{Int}[]
    order = 1
    ref_face = lagrange_mesh_face(geo,order)
    fun_mesh = boundary(ref_face)
    geo_node_coords = node_coordinates(geo_mesh)
    fun_node_coords = node_coordinates(fun_mesh)
    vertex_coords = geo_node_coords[vertex_to_geo_node]
    degree = 1
    quad = default_quadrature(geo,degree)
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

"""
"""
function interior_node_permutations(fe::AbstractMeshFace)
    interior_ho_nodes = interior_nodes(fe)
    node_permutations_from_mesh_face(fe,interior_ho_nodes)
end

"""
"""
function node_permutations(fe::AbstractMeshFace)
    interior_ho_nodes = 1:num_nodes(fe)
    node_permutations_from_mesh_face(fe,interior_ho_nodes)
end

function node_permutations_from_mesh_face(refface,interior_ho_nodes)
    ho_nodes_coordinates = node_coordinates(refface)
    geo = geometry(refface)
    vertex_perms = vertex_permutations(geo)
    if length(interior_ho_nodes) == 0
        return map(i->Int[],vertex_perms)
    end
    if length(vertex_perms) == 1
        return map(i->collect(1:length(interior_ho_nodes)),vertex_perms)
    end
    if order(refface) == 0 # TODO ugly. It assumes the hack above for node coordinates of faces of order 0
        return map(i->collect(1:length(interior_ho_nodes)),vertex_perms)
    end
    geo_mesh = boundary(geo)
    vertex_to_geo_nodes = face_nodes(geo_mesh,0)
    vertex_to_geo_node = map(first,vertex_to_geo_nodes)
    ref_face = lagrange_mesh_face(geo,1)
    fun_mesh = boundary(ref_face)
    geo_node_coords = node_coordinates(geo_mesh)
    fun_node_coords = node_coordinates(fun_mesh)
    vertex_coords = geo_node_coords[vertex_to_geo_node]
    q = ho_nodes_coordinates[interior_ho_nodes]
    Tx = eltype(vertex_coords)
    A = zeros(Float64,length(q),length(fun_node_coords))
    A = tabulator(ref_face)(value,q)
    perm_vertex_coords = similar(vertex_coords)
    node_perms = similar(vertex_perms)
    for (iperm,permutation) in enumerate(vertex_perms)
        for (j,cj) in enumerate(permutation)
          perm_vertex_coords[j] = vertex_coords[cj]
        end
        node_to_pnode = fill(INVALID_ID,length(interior_ho_nodes))
        for iq in 1:size(A,1)
            y = zero(Tx)
            for fun_node in 1:size(A,2)
                vertex = fun_node # TODO we are assuming that the vertices and nodes match
                g = A[iq,fun_node]
                x = perm_vertex_coords[vertex]
                y += g*x
            end
            pnode = findfirst(i->(norm(i-y)+1)≈1,q)
            if !isnothing(pnode)
               node_to_pnode[iq] = pnode
            end
        end
        node_perms[iperm] = node_to_pnode
    end
    node_perms
end

"""
    abstract type AbstractMeshTopology

# Basic queries

- [`face_incidence`](@ref)
- [`face_reference_id`](@ref)
- [`face_permutation_ids`](@ref)
- [`reference_faces`](@ref)

# Basic constructors

- [`topology`](@ref)

"""
abstract type AbstractMeshTopology <: GT.AbstractType end

struct GenericMeshTopology{A,B,C,D} <: AbstractMeshTopology
    face_incidence::A
    face_reference_id::B
    face_permutation_ids::C
    reference_faces::D
end

function mesh_topology(args...)
    GenericMeshTopology(args...)
end

"""
"""
function topology(mesh::AbstractMesh)
    topology_from_mesh(mesh)
end

function topology_from_mesh(mesh)
    # Assumes that the input is a cell complex
    T = JaggedArray{Int32,Int32}
    D = num_dims(mesh)
    my_face_incidence = Matrix{T}(undef,D+1,D+1)
    my_face_reference_id  = [ face_reference_id(mesh,d) for d in 0:D ]
    my_reference_faces = Tuple([ map(topology,reference_faces(mesh,d)) for d in 0:D ])
    my_face_permutation_ids = Matrix{T}(undef,D+1,D+1)
    topo = mesh_topology(
        my_face_incidence,
        my_face_reference_id,
        my_face_permutation_ids,
        my_reference_faces,
       )
    for d in 0:D
        fill_face_interior_mesh_topology!(topo,mesh,d)
    end
    for d in 1:D
        fill_face_vertices_mesh_topology!(topo,mesh,d)
        fill_face_coboundary_mesh_topology!(topo,mesh,d,0)
    end
    for d in 1:(D-1)
        for n in (D-d):-1:1
            m = n+d
            fill_face_boundary_mesh_topology!(topo,mesh,m,n)
            fill_face_coboundary_mesh_topology!(topo,mesh,m,n)
        end
    end
    for d in 0:D
        for n in 0:d
            fill_face_permutation_ids!(topo,d,n)
        end
    end
    topo
end

function fill_face_interior_mesh_topology!(mesh_topology,mesh,d)
    n = num_faces(mesh,d)
    ptrs = collect(Int32,1:(n+1))
    data = collect(Int32,1:n)
    mesh_topology.face_incidence[d+1,d+1] = JaggedArray(data,ptrs)
end

function generate_face_coboundary(nface_to_mfaces,nmfaces)
    ptrs = zeros(Int32,nmfaces+1)
    nnfaces = length(nface_to_mfaces)
    for nface in 1:nnfaces
        mfaces = nface_to_mfaces[nface]
        for mface in mfaces
            ptrs[mface+1] += Int32(1)
        end
    end
    length_to_ptrs!(ptrs)
    ndata = ptrs[end]-1
    data = zeros(Int32,ndata)
    for nface in 1:nnfaces
        mfaces = nface_to_mfaces[nface]
        for mface in mfaces
            p = ptrs[mface]
            data[p] = nface
            ptrs[mface] += Int32(1)
        end
    end
    rewind_ptrs!(ptrs)
    mface_to_nfaces = JaggedArray(data,ptrs)
    mface_to_nfaces
end

function fill_face_coboundary_mesh_topology!(topology,mesh,n,m)
    nmfaces = num_faces(mesh,m)
    nface_to_mfaces = face_incidence(topology,n,m)
    face_incidence(topology)[m+1,n+1] = generate_face_coboundary(nface_to_mfaces,nmfaces)
end

function fill_face_vertices_mesh_topology!(topo,mesh,d)
    function barrier(nnodes,vertex_to_nodes,dface_to_nodes,dface_to_refid,refid_to_lvertex_to_lnodes)
        vertex_to_node = JaggedArray(vertex_to_nodes).data
        #if vertex_to_node == 1:nnodes
        #    return dface_to_nodes
        #end
        node_to_vertex = zeros(Int32,nnodes)
        nvertices = length(vertex_to_nodes)
        node_to_vertex[vertex_to_node] = 1:nvertices
        ndfaces = length(dface_to_nodes)
        dface_to_vertices_ptrs = zeros(Int32,ndfaces+1)
        for dface in 1:ndfaces
            refid = dface_to_refid[dface]
            nlvertices = length(refid_to_lvertex_to_lnodes[refid])
            dface_to_vertices_ptrs[dface+1] = nlvertices
        end
        length_to_ptrs!(dface_to_vertices_ptrs)
        ndata = dface_to_vertices_ptrs[end]-1
        dface_to_vertices_data = zeros(Int32,ndata)
        for dface in 1:ndfaces
            refid = dface_to_refid[dface]
            lvertex_to_lnodes = refid_to_lvertex_to_lnodes[refid]
            nlvertices = length(lvertex_to_lnodes)
            lnode_to_node = dface_to_nodes[dface]
            offset = dface_to_vertices_ptrs[dface]-1
            for lvertex in 1:nlvertices
                lnode = first(lvertex_to_lnodes[lvertex])
                node = lnode_to_node[lnode]
                vertex = node_to_vertex[node]
                dface_to_vertices_data[offset+lvertex] = vertex
            end
        end
        dface_to_vertices = JaggedArray(dface_to_vertices_data,dface_to_vertices_ptrs)
    end
    nnodes = num_nodes(mesh)
    vertex_to_nodes = face_nodes(mesh,0)
    dface_to_nodes = face_nodes(mesh,d)
    dface_to_refid = face_reference_id(mesh,d)
    refid_refface = reference_faces(mesh,d)
    refid_to_lvertex_to_lnodes = map(refface->face_nodes(boundary(refface),0),refid_refface)
    face_incidence(topo)[d+1,0+1] = barrier(nnodes,vertex_to_nodes,dface_to_nodes,dface_to_refid,refid_to_lvertex_to_lnodes)
end

function fill_face_boundary_mesh_topology!(topo,mesh,D,d)
    function barrier(
            Dface_to_vertices,
            vertex_to_Dfaces,
            dface_to_vertices,
            vertex_to_dfaces,
            Dface_to_refid,
            Drefid_to_ldface_to_lvertices)

        # Count
        ndfaces = length(dface_to_vertices)
        nDfaces = length(Dface_to_vertices)
        # Allocate output
        ptrs = zeros(Int32,nDfaces+1)
        for Dface in 1:nDfaces
            Drefid = Dface_to_refid[Dface]
            ldface_to_lvertices = Drefid_to_ldface_to_lvertices[Drefid]
            ptrs[Dface+1] = length(ldface_to_lvertices)
        end
        length_to_ptrs!(ptrs)
        ndata = ptrs[end]-1
        data = fill(Int32(INVALID_ID),ndata)
        Dface_to_dfaces = JaggedArray(data,ptrs)
        for Dface in 1:nDfaces
            Drefid = Dface_to_refid[Dface]
            ldface_to_lvertices = Drefid_to_ldface_to_lvertices[Drefid]
            lvertex_to_vertex = Dface_to_vertices[Dface]
            ldface_to_dface = Dface_to_dfaces[Dface]
            for (ldface,lvertices) in enumerate(ldface_to_lvertices)
                # Find the global d-face for this local d-face
                dface2 = Int32(INVALID_ID)
                vertices = view(lvertex_to_vertex,lvertices)
                for (i,lvertex) in enumerate(lvertices)
                    vertex = lvertex_to_vertex[lvertex]
                    dfaces = vertex_to_dfaces[vertex]
                    for dface1 in dfaces
                        vertices1 = dface_to_vertices[dface1]
                        if same_valid_ids(vertices,vertices1)
                            dface2 = dface1
                            break
                        end
                    end
                    if dface2 != Int32(INVALID_ID)
                        break
                    end
                end
                @boundscheck begin
                    msg = """


                    Error in: topology_from_mesh

                    The given mesh is provably not a cell complex.
                    """
                    @assert dface2 != Int32(INVALID_ID) msg
                end
                ldface_to_dface[ldface] = dface2
            end # (ldface,lvertices)
        end # Dface
        Dface_to_dfaces
    end
    Dface_to_vertices = face_incidence(topo,D,0)
    vertex_to_Dfaces = face_incidence(topo,0,D)
    dface_to_vertices = face_incidence(topo,d,0)
    vertex_to_dfaces = face_incidence(topo,0,d)
    Dface_to_refid = face_reference_id(topo,D)
    refid_refface = reference_faces(topo,D)
    Drefid_to_ldface_to_lvertices = map(refface->face_incidence(boundary(refface),d,0),refid_refface)
    Dface_to_dfaces = barrier(
            Dface_to_vertices,
            vertex_to_Dfaces,
            dface_to_vertices,
            vertex_to_dfaces,
            Dface_to_refid,
            Drefid_to_ldface_to_lvertices)
    topo.face_incidence[D+1,d+1] = Dface_to_dfaces
end

function fill_face_permutation_ids!(top,D,d)
    function barrier!(
            cell_to_lface_to_pindex,
            cell_to_lface_to_face,
            cell_to_cvertex_to_vertex,
            cell_to_ctype,
            ctype_to_lface_to_cvertices,
            face_to_fvertex_to_vertex,
            face_to_ftype,
            ftype_to_pindex_to_cfvertex_to_fvertex)

        ncells = length(cell_to_lface_to_face)
        for cell in 1:ncells
            ctype = cell_to_ctype[cell]
            lface_to_cvertices = ctype_to_lface_to_cvertices[ctype]
            a = cell_to_lface_to_face.ptrs[cell]-1
            c = cell_to_cvertex_to_vertex.ptrs[cell]-1
            for (lface,cfvertex_to_cvertex) in enumerate(lface_to_cvertices)
                face = cell_to_lface_to_face.data[a+lface]
                ftype = face_to_ftype[face]
                b = face_to_fvertex_to_vertex.ptrs[face]-1
                pindex_to_cfvertex_to_fvertex = ftype_to_pindex_to_cfvertex_to_fvertex[ftype]
                pindexfound = false
                for (pindex, cfvertex_to_fvertex) in enumerate(pindex_to_cfvertex_to_fvertex)
                    found = true
                    for (cfvertex,fvertex) in enumerate(cfvertex_to_fvertex)
                        vertex1 = face_to_fvertex_to_vertex.data[b+fvertex]
                        cvertex = cfvertex_to_cvertex[cfvertex]
                        vertex2 = cell_to_cvertex_to_vertex.data[c+cvertex]
                        if vertex1 != vertex2
                            found = false
                            break
                        end
                    end
                    if found
                        cell_to_lface_to_pindex.data[a+lface] = pindex
                        pindexfound = true
                        break
                    end
                end
                @assert pindexfound "Valid pindex not found"
            end
        end
    end
    @assert D >= d
    cell_to_lface_to_face = JaggedArray(face_incidence(top,D,d))
    data = similar(cell_to_lface_to_face.data,Int8)
    ptrs = cell_to_lface_to_face.ptrs
    cell_to_lface_to_pindex = JaggedArray(data,ptrs)
    if d == D || d == 0
        fill!(cell_to_lface_to_pindex.data,Int8(1))
        top.face_permutation_ids[D+1,d+1] = cell_to_lface_to_pindex
        return top
    end
    face_to_fvertex_to_vertex = JaggedArray(face_incidence(top,d,0))
    face_to_ftype = face_reference_id(top,d)
    ref_dfaces = reference_faces(top,d)
    ftype_to_pindex_to_cfvertex_to_fvertex = map(vertex_permutations,ref_dfaces)
    cell_to_cvertex_to_vertex = JaggedArray(face_incidence(top,D,0))
    cell_to_ctype = face_reference_id(top,D)
    ref_Dfaces = reference_faces(top,D)
    ctype_to_lface_to_cvertices = map(p->face_incidence(boundary(p),d,0),ref_Dfaces)
    barrier!(
             cell_to_lface_to_pindex,
             cell_to_lface_to_face,
             cell_to_cvertex_to_vertex,
             cell_to_ctype,
             ctype_to_lface_to_cvertices,
             face_to_fvertex_to_vertex,
             face_to_ftype,
             ftype_to_pindex_to_cfvertex_to_fvertex)
    top.face_permutation_ids[D+1,d+1] = cell_to_lface_to_pindex
    top
end

function intersection!(a,b,na,nb)
  function findeq!(i,a,b,nb)
    for j in 1:nb
      if a[i] == b[j]
        return
      end
    end
    a[i] = INVALID_ID
    return
  end
  for i in 1:na
    if a[i] == INVALID_ID
      continue
    end
    findeq!(i,a,b,nb)
  end
end

function find_eq(v,b)
    for vs in b
        if v == vs
            return true
        end
    end
    return false
end

function is_subset(a,b)
    for i in 1:length(a)
        v = a[i]
        if v == INVALID_ID
            continue
        end
        if !find_eq(v,b)
            return false
        end
    end
    return true
end

function same_valid_ids(a,b)
    if !is_subset(a,b)
        return false
    end
    if !is_subset(b,a)
        return false
    end
    return true
end

## TODO AbstractFaceTopology <: AbstractMeshTopology
"""
    abstract type AbstractFaceTopology

# Basic queries

- [`boundary`](@ref)
- [`vertex_permutations`](@ref)

# Basic constructors

- [`topology`](@ref)

"""
abstract type AbstractFaceTopology <: GT.AbstractType end

struct GenericFaceTopology{A,B} <: AbstractFaceTopology
    boundary::A
    vertex_permutations::B
end

function face_topology(args...)
    GenericFaceTopology(args...)
end

num_dims(a::AbstractFaceTopology) = num_dims(boundary(a))+1

function topology(fe::AbstractMeshFace)
    topology_from_mesh_face(fe)
end

function topology_from_mesh_face(refface)
    geom = geometry(refface)
    D = num_dims(geom)
    if D != 0
        myboundary = geom |> boundary |> topology
    else
        myboundary = nothing
    end
    myperms = vertex_permutations(geom)
    face_topology(myboundary,myperms)
end

"""
"""
function complexify(mesh::AbstractMesh;kwargs...)
    complexify_mesh(mesh;kwargs...)
end

function complexify_mesh(mesh;glue=Val(false))
    Ti = Int32
    T = JaggedArray{Ti,Ti}
    D = num_dims(mesh)
    oldface_to_newvertices = Vector{T}(undef,D+1)
    newvertex_to_oldfaces = Vector{T}(undef,D+1)
    newface_incidence = Matrix{T}(undef,D+1,D+1)
    nnewfaces = zeros(Int,D+1)
    newface_refid = Vector{Vector{Ti}}(undef,D+1)
    newreffaces = Vector{Any}(undef,D+1)
    newface_nodes = Vector{T}(undef,D+1)
    old_to_new = Vector{Vector{Ti}}(undef,D+1)
    node_to_newvertex, n_new_vertices = find_node_to_vertex(mesh) # Optimizable for linear meshes
    for d in 0:D
        oldface_to_newvertices[d+1] = fill_face_vertices(mesh,d,node_to_newvertex) # Optimizable for linear meshes
        newvertex_to_oldfaces[d+1] = generate_face_coboundary(oldface_to_newvertices[d+1],n_new_vertices) # Optimizable for linear meshes
    end
    newface_incidence[D+1,0+1] = oldface_to_newvertices[D+1]
    newface_incidence[0+1,D+1] = newvertex_to_oldfaces[D+1]
    nnewfaces[D+1] = length(oldface_to_newvertices[D+1])
    newface_refid[D+1] = face_reference_id(mesh,D)
    newreffaces[D+1] = reference_faces(mesh,D)
    newface_nodes[D+1] = face_nodes(mesh,D)
    old_to_new[D+1] = collect(Ti,1:length(newface_nodes[D+1]))
    # TODO optimize for d==0
    for d in (D-1):-1:0
        n = d+1
        new_nface_to_new_vertices = newface_incidence[n+1,0+1]
        new_vertex_to_new_nfaces = newface_incidence[0+1,n+1]
        old_dface_to_new_vertices = oldface_to_newvertices[d+1]
        new_vertex_to_old_dfaces = newvertex_to_oldfaces[d+1]
        new_nface_to_nrefid = newface_refid[n+1]
        old_dface_to_drefid = face_reference_id(mesh,d)
        drefid_to_ref_dface = reference_faces(mesh,d)
        old_dface_to_nodes = face_nodes(mesh,d)
        new_nface_to_nodes = newface_nodes[n+1]
        nrefid_to_ldface_to_lvertices = map(a->face_incidence(topology(boundary(geometry(a))),d,0),newreffaces[n+1])
        nrefid_to_ldface_to_lnodes = map(a->face_nodes(boundary(a),d),newreffaces[n+1])
        nrefid_to_ldface_to_drefrefid = map(a->face_reference_id(boundary(a),d),newreffaces[n+1])
        nrefid_to_drefrefid_to_ref_dface = map(a->reference_faces(boundary(a),d),newreffaces[n+1])
        new_nface_to_new_dfaces, n_new_dfaces, old_dface_to_new_dface = generate_face_boundary(
            new_nface_to_new_vertices,
            new_vertex_to_new_nfaces,
            old_dface_to_new_vertices,
            new_vertex_to_old_dfaces,
            new_nface_to_nrefid,
            nrefid_to_ldface_to_lvertices)
        new_dface_to_new_vertices = generate_face_vertices(
            n_new_dfaces,
            old_dface_to_new_dface,
            old_dface_to_new_vertices,
            new_nface_to_new_vertices,
            new_nface_to_new_dfaces,
            new_nface_to_nrefid,
            nrefid_to_ldface_to_lvertices)
        new_vertex_to_new_dfaces = generate_face_coboundary(new_dface_to_new_vertices,n_new_vertices)
        new_dface_to_nodes = generate_face_vertices(
            n_new_dfaces,
            old_dface_to_new_dface,
            old_dface_to_nodes,
            new_nface_to_nodes,
            new_nface_to_new_dfaces,
            new_nface_to_nrefid,
            nrefid_to_ldface_to_lnodes)
        new_dface_to_new_drefid, new_refid_to_ref_dface = generate_reference_faces(
            n_new_dfaces,
            old_dface_to_new_dface,
            old_dface_to_drefid,
            drefid_to_ref_dface,
            new_nface_to_new_dfaces,
            new_nface_to_nrefid,
            nrefid_to_ldface_to_drefrefid,
            nrefid_to_drefrefid_to_ref_dface)
        newface_incidence[n+1,d+1] = new_nface_to_new_dfaces
        newface_incidence[d+1,0+1] = new_dface_to_new_vertices
        newface_incidence[0+1,d+1] = new_vertex_to_new_dfaces
        newface_refid[d+1] = new_dface_to_new_drefid
        newreffaces[d+1] = new_refid_to_ref_dface
        nnewfaces[d+1] = n_new_dfaces
        newface_nodes[d+1] = new_dface_to_nodes
        old_to_new[d+1] = old_dface_to_new_dface
    end
    node_to_coords = node_coordinates(mesh)
    old_physical_faces = physical_faces(mesh)
    new_physical_faces = [ Dict{String,Vector{Int32}}() for d in 0:D] # TODO hardcoded
    for d in 0:D
        old_groups = old_physical_faces[d+1]
        for (group_name,old_group_faces) in old_groups
            new_group_faces = similar(old_group_faces)
            new_group_faces .= old_to_new[d+1][old_group_faces]
            new_physical_faces[d+1][group_name] = new_group_faces
        end
    end
    new_mesh = GT.mesh_from_arrays(
            node_to_coords,
            newface_nodes,
            newface_refid,
            Tuple(newreffaces);
            physical_faces = new_physical_faces,
            periodic_nodes = periodic_nodes(mesh),
            outwards_normals = outwards_normals(mesh)
           )
    if val_parameter(glue)
        new_mesh, old_to_new
    else
        new_mesh
    end
end

function generate_face_vertices(
    n_new_dfaces,
    old_dface_to_new_dface,
    old_dface_to_new_vertices,
    new_nface_to_new_vertices,
    new_nface_to_new_dfaces,
    new_nface_to_nrefid,
    nrefid_to_ldface_to_lvertices
    )

    Ti = eltype(eltype(old_dface_to_new_vertices))
    new_dface_to_touched = fill(false,n_new_dfaces)
    new_dface_to_new_vertices_ptrs = zeros(Ti,n_new_dfaces+1)
    n_old_dfaces = length(old_dface_to_new_dface)
    for old_dface in 1:n_old_dfaces
        new_dface = old_dface_to_new_dface[old_dface]
        new_vertices = old_dface_to_new_vertices[old_dface]
        new_dface_to_new_vertices_ptrs[new_dface+1] = length(new_vertices)
        new_dface_to_touched[new_dface] = true
    end
    n_new_nfaces = length(new_nface_to_new_vertices)
    for new_nface in 1:n_new_nfaces
        nrefid = new_nface_to_nrefid[new_nface]
        ldface_to_lvertices = nrefid_to_ldface_to_lvertices[nrefid]
        ldface_to_new_dface = new_nface_to_new_dfaces[new_nface]
        n_ldfaces = length(ldface_to_new_dface)
        for ldface in 1:n_ldfaces
            new_dface = ldface_to_new_dface[ldface]
            if new_dface_to_touched[new_dface]
                continue
            end
            lvertices = ldface_to_lvertices[ldface]
            new_dface_to_new_vertices_ptrs[new_dface+1] = length(lvertices)
            new_dface_to_touched[new_dface] = true
        end

    end
    length_to_ptrs!(new_dface_to_new_vertices_ptrs)
    ndata = new_dface_to_new_vertices_ptrs[end]-1
    new_dface_to_new_vertices_data = zeros(Ti,ndata)
    new_dface_to_new_vertices = JaggedArray(new_dface_to_new_vertices_data,new_dface_to_new_vertices_ptrs)
    fill!(new_dface_to_touched,false)
    for old_dface in 1:n_old_dfaces
        new_dface = old_dface_to_new_dface[old_dface]
        new_vertices_in = old_dface_to_new_vertices[old_dface]
        new_vertices_out = new_dface_to_new_vertices[new_dface]
        for i in 1:length(new_vertices_in)
            new_vertices_out[i] = new_vertices_in[i]
        end
        new_dface_to_touched[new_dface] = true
    end
    for new_nface in 1:n_new_nfaces
        nrefid = new_nface_to_nrefid[new_nface]
        ldface_to_lvertices = nrefid_to_ldface_to_lvertices[nrefid]
        ldface_to_new_dface = new_nface_to_new_dfaces[new_nface]
        n_ldfaces = length(ldface_to_new_dface)
        new_vertices_in = new_nface_to_new_vertices[new_nface]
        for ldface in 1:n_ldfaces
            new_dface = ldface_to_new_dface[ldface]
            if new_dface_to_touched[new_dface]
                continue
            end
            new_vertices_out = new_dface_to_new_vertices[new_dface]
            lvertices = ldface_to_lvertices[ldface]
            for i in 1:length(lvertices)
                new_vertices_out[i] = new_vertices_in[lvertices[i]]
            end
            new_dface_to_touched[new_dface] = true
        end

    end
    new_dface_to_new_vertices
end

function generate_reference_faces(
        n_new_dfaces,
        old_dface_to_new_dface,
        old_dface_to_drefid,
        drefid_to_ref_dface,
        new_nface_to_new_dfaces,
        new_nface_to_nrefid,
        nrefid_to_ldface_to_drefrefid,
        nrefid_to_drefrefid_to_ref_dface)

    i_to_ref_dface = collect(Any,drefid_to_ref_dface)
    drefid_to_i = collect(1:length(drefid_to_ref_dface))
    i = length(i_to_ref_dface)
    Ti = Int32
    nrefid_to_drefrefid_to_i = map(a->zeros(Ti,length(a)),nrefid_to_drefrefid_to_ref_dface)
    for (nrefid,drefrefid_to_ref_dface) in enumerate(nrefid_to_drefrefid_to_ref_dface)
        for (drefrefid, ref_dface) in enumerate(drefrefid_to_ref_dface)
            push!(i_to_ref_dface,ref_dface)
            i += 1
            nrefid_to_drefrefid_to_i[nrefid][drefrefid] = i
        end
    end
    u_to_ref_dface = unique(i_to_ref_dface)
    i_to_u = indexin(i_to_ref_dface,u_to_ref_dface)
    new_dface_to_u = zeros(Ti,n_new_dfaces)
    new_dface_to_touched = fill(false,n_new_dfaces)
    for (old_dface,new_dface) in enumerate(old_dface_to_new_dface)
        drefid = old_dface_to_drefid[old_dface]
        i = drefid_to_i[drefid]
        u = i_to_u[i]
        new_dface_to_u[new_dface] = u
        new_dface_to_touched[new_dface] = true
    end
    for (new_nface,nrefid) in enumerate(new_nface_to_nrefid)
        ldface_to_drefrefid = nrefid_to_ldface_to_drefrefid[nrefid]
        ldface_to_new_dface = new_nface_to_new_dfaces[new_nface]
        drefrefid_to_i = nrefid_to_drefrefid_to_i[nrefid]
        for (ldface,new_dface) in enumerate(ldface_to_new_dface)
            new_dface = ldface_to_new_dface[ldface]
            if new_dface_to_touched[new_dface]
                continue
            end
            drefrefid = ldface_to_drefrefid[ldface]
            i = drefrefid_to_i[drefrefid]
            u = i_to_u[i]
            new_dface_to_u[new_dface] = u
            new_dface_to_touched[new_dface] = true
        end
    end
    new_dface_to_u, Tuple(u_to_ref_dface)
end

function generate_face_boundary(
    Dface_to_vertices,
    vertex_to_Dfaces,
    dface_to_vertices,
    vertex_to_dfaces,
    Dface_to_refid,
    Drefid_to_ldface_to_lvertices)

    # Count
    ndfaces = length(dface_to_vertices)
    nDfaces = length(Dface_to_vertices)
    nvertices = length(vertex_to_Dfaces)
    maxldfaces = 0
    for ldface_to_lvertices in Drefid_to_ldface_to_lvertices
        maxldfaces = max(maxldfaces,length(ldface_to_lvertices))
    end
    maxDfaces = 0
    for vertex in 1:length(vertex_to_Dfaces)
        Dfaces = vertex_to_Dfaces[vertex]
        maxDfaces = max(maxDfaces,length(Dfaces))
    end
    # Allocate output
    ptrs = zeros(Int32,nDfaces+1)
    for Dface in 1:nDfaces
        Drefid = Dface_to_refid[Dface]
        ldface_to_lvertices = Drefid_to_ldface_to_lvertices[Drefid]
        ptrs[Dface+1] = length(ldface_to_lvertices)
    end
    length_to_ptrs!(ptrs)
    ndata = ptrs[end]-1
    data = fill(Int32(INVALID_ID),ndata)
    Dface_to_dfaces = GenericJaggedArray(data,ptrs)
    # Main loop
    Dfaces1 = fill(Int32(INVALID_ID),maxDfaces)
    Dfaces2 = fill(Int32(INVALID_ID),maxDfaces)
    ldfaces1 = fill(Int32(INVALID_ID),maxDfaces)
    nDfaces1 = 0
    nDfaces2 = 0
    newdface = Int32(ndfaces)
    old_to_new = collect(Int32,1:ndfaces)
    for Dface in 1:nDfaces
        Drefid = Dface_to_refid[Dface]
        ldface_to_lvertices = Drefid_to_ldface_to_lvertices[Drefid]
        lvertex_to_vertex = Dface_to_vertices[Dface]
        ldface_to_dface = Dface_to_dfaces[Dface]
        for (ldface,lvertices) in enumerate(ldface_to_lvertices)
            # Do nothing if this local face has already been processed by
            # a neighbor
            if ldface_to_dface[ldface] != Int32(INVALID_ID)
                continue
            end
            # Find if there is already a global d-face for this local d-face
            # if yes, then use the global id of this d-face
            # if not, create a new one
            dface2 = Int32(INVALID_ID)
            fill!(Dfaces1,Int32(INVALID_ID))
            fill!(Dfaces2,Int32(INVALID_ID))
            vertices = view(lvertex_to_vertex,lvertices)
            for (i,lvertex) in enumerate(lvertices)
                vertex = lvertex_to_vertex[lvertex]
                dfaces = vertex_to_dfaces[vertex]
                for dface1 in dfaces
                    vertices1 = dface_to_vertices[dface1]
                    if same_valid_ids(vertices,vertices1)
                        dface2 = dface1
                        break
                    end
                end
                if dface2 != Int32(INVALID_ID)
                    break
                end
            end
            if dface2 == Int32(INVALID_ID)
                newdface += Int32(1)
                dface2 = newdface
            end
            # Find all D-faces around this local d-face
            for (i,lvertex) in enumerate(lvertices)
                vertex = lvertex_to_vertex[lvertex]
                Dfaces = vertex_to_Dfaces[vertex]
                if i == 1
                    copyto!(Dfaces1,Dfaces)
                    nDfaces1 = length(Dfaces)
                else
                    copyto!(Dfaces2,Dfaces)
                    nDfaces2 = length(Dfaces)
                    intersection!(Dfaces1,Dfaces2,nDfaces1,nDfaces2)
                end
            end
            # Find their correspondent local d-face and set the d-face
            for Dface1 in Dfaces1
                if Dface1 != INVALID_ID
                    Drefid1 = Dface_to_refid[Dface1]
                    lvertex1_to_vertex1 = Dface_to_vertices[Dface1]
                    ldface1_to_lvertices1 = Drefid_to_ldface_to_lvertices[Drefid1]
                    ldface2 = Int32(INVALID_ID)
                    for (ldface1,lvertices1) in enumerate(ldface1_to_lvertices1)
                        vertices1 = view(lvertex1_to_vertex1,lvertices1)
                        if same_valid_ids(vertices,vertices1)
                            ldface2 = ldface1
                            break
                        end
                    end
                    @boundscheck @assert ldface2 != INVALID_ID # TODO: issue gmsh quad eles
                    Dface_to_dfaces[Dface1][ldface2] = dface2
                end
            end
        end # (ldface,lvertices)
    end # Dface
    Dface_to_dfaces, newdface, old_to_new
end

function fill_face_vertices(mesh,d,node_to_vertex)
    function barrier(node_to_vertex,face_to_nodes,face_to_refid,refid_to_lvertex_to_lnodes)
        Ti = eltype(node_to_vertex)
        nfaces = length(face_to_nodes)
        face_to_vertices_ptrs = zeros(Ti,nfaces+1)
        for face in 1:nfaces
            nodes = face_to_nodes[face]
            refid = face_to_refid[face]
            lvertex_to_lnodes = refid_to_lvertex_to_lnodes[refid]
            nlvertices = length(lvertex_to_lnodes)
            face_to_vertices_ptrs[face+1] = nlvertices
        end
        length_to_ptrs!(face_to_vertices_ptrs)
        ndata = face_to_vertices_ptrs[end]-1
        face_to_vertices_data = zeros(Ti,ndata)
        face_to_vertices = JaggedArray(face_to_vertices_data,face_to_vertices_ptrs)
        for face in 1:nfaces
            vertices = face_to_vertices[face]
            nodes = face_to_nodes[face]
            refid = face_to_refid[face]
            lvertex_to_lnodes = refid_to_lvertex_to_lnodes[refid]
            for (lvertex,lnodes) in enumerate(lvertex_to_lnodes)
                lnode = first(lnodes)
                vertex = node_to_vertex[nodes[lnode]]
                @boundscheck @assert vertex != INVALID_ID
                vertices[lvertex] = vertex
            end
        end
        face_to_vertices
    end
    face_to_nodes = face_nodes(mesh,d)
    face_to_refid = face_reference_id(mesh,d)
    refid_to_lvertex_to_lnodes = map(reference_faces(mesh,d)) do a
        if num_dims(geometry(a)) != 0
            face_nodes(boundary(a),0)
        else
            [interior_nodes(a)]
        end
    end
    barrier(node_to_vertex,face_to_nodes,face_to_refid,refid_to_lvertex_to_lnodes)
end

function find_node_to_vertex(mesh)
    function barrier!(node_to_vertex,face_to_nodes,face_to_refid,refid_to_lvertex_to_lnodes)
        valid_id = one(eltype(node_to_vertex))
        for face in eachindex(face_to_nodes)
            nodes = face_to_nodes[face]
            refid = face_to_refid[face]
            lvertex_to_lnodes = refid_to_lvertex_to_lnodes[refid]
            for lnodes in lvertex_to_lnodes
                lnode = first(lnodes)
                node = nodes[lnode]
                node_to_vertex[node] = valid_id
            end
        end
    end
    Ti = Int32
    nnodes = num_nodes(mesh)
    node_to_vertex = zeros(Ti,nnodes)
    fill!(node_to_vertex,Ti(INVALID_ID))
    D = num_dims(mesh)
    for d in 0:D
        face_to_nodes = face_nodes(mesh,d)
        if length(face_to_nodes) == 0
            continue
        end
        face_to_refid = face_reference_id(mesh,d)
        refid_to_lvertex_to_lnodes = map(reference_faces(mesh,d)) do a
            if num_dims(geometry(a)) != 0
                face_nodes(boundary(a),0)
            else
                [interior_nodes(a)]
            end
        end
        barrier!(node_to_vertex,face_to_nodes,face_to_refid,refid_to_lvertex_to_lnodes)
    end
    vertex = Ti(0)
    for node in eachindex(node_to_vertex)
        if node_to_vertex[node] != Ti(INVALID_ID)
            vertex += Ti(1)
            node_to_vertex[node] = vertex
        end
    end
    node_to_vertex, vertex
end

function physical_nodes(mesh,d)
    nnodes = num_nodes(mesh)
    node_to_touched = fill(false,nnodes)
    node_groups = Dict{String,Vector{Int32}}()
    face_to_nodes = face_nodes(mesh,d)
    for (name,faces) in physical_faces(mesh,d)
        fill!(node_to_touched,false)
        for face in faces
            nodes = face_to_nodes[face]
            node_to_touched[nodes] .= true
        end
        node_groups[name] = findall(node_to_touched)
    end
    node_groups
end

"""
"""
function physical_nodes(mesh;
    merge_dims=Val(false),
    disjoint=Val(false),
    name_priority=nothing)

    if val_parameter(disjoint) == true && val_parameter(merge_dims) == false
        error("disjoint=true requires merge_dims=true")
    end
    D = num_dims(mesh)
    d_to_groups = [ physical_nodes(mesh,d) for d in 0:D ]
    if val_parameter(merge_dims) == false
        return d_to_groups
    end
    names = physical_names(mesh;merge_dims)
    nnodes = num_nodes(mesh)
    node_groups = Dict{String,Vector{Int32}}()

    if val_parameter(disjoint) == false
        node_to_touched = fill(false,nnodes)
        for name in names
            fill!(node_to_touched,false)
            for groups in d_to_groups
                for (name2,nodes) in groups
                    if name != name2
                        continue
                    end
                    node_to_touched[nodes] .= true
                end
            end
            node_groups[name] = findall(node_to_touched)
        end
    else
        if name_priority === nothing
            tag_to_name = reverse(sort(collect(names)))
        else
            tag_to_name = name_priority
        end
        node_to_tag = zeros(Int32,nnodes)
        classify_mesh_nodes!(node_to_tag,mesh,tag_to_name)
        for (tag,name) in enumerate(tag_to_name)
            node_groups[name] = findall(t->t==tag,node_to_tag)
        end
    end
    node_groups
end

"""
"""
function classify_mesh_nodes!(node_to_tag,mesh,tag_to_name,dmax=num_dims(mesh))
    fill!(node_to_tag,zero(eltype(node_to_tag)))
    for d in dmax:-1:0
        face_to_nodes = face_nodes(mesh,d)
        face_groups = physical_faces(mesh,d)
        for (tag,name) in enumerate(tag_to_name)
            for (name2,faces) in face_groups
                if name != name2
                    continue
                end
                for face in faces
                    nodes = face_to_nodes[face]
                    node_to_tag[nodes] .= tag
                end
            end
        end
    end
    node_to_tag
end

function classify_mesh_faces!(dface_to_tag,mesh,d,tag_to_name)
    fill!(dface_to_tag,zero(eltype(dface_to_tag)))
    face_groups = physical_faces(mesh,d)
    for (tag,name) in enumerate(tag_to_name)
        for (name2,faces) in face_groups
            if name != name2
                continue
            end
            dface_to_tag[faces] .= tag
        end
    end
    dface_to_tag
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

function label_interior_faces!(mesh::AbstractMesh;physical_name="interior")
    D = num_dims(mesh)
    d = D-1
    topo = topology(mesh)
    face_to_cells = face_incidence(topo,d,D)
    groups = physical_faces(mesh,d)
    faces = findall(cells->length(cells)==2,face_to_cells)
    groups[physical_name] = faces
    mesh
end

function label_boundary_faces!(mesh::AbstractMesh;physical_name="boundary")
    D = num_dims(mesh)
    d = D-1
    topo = topology(mesh)
    face_to_cells = face_incidence(topo,d,D)
    groups = physical_faces(mesh,d)
    faces = findall(cells->length(cells)==1,face_to_cells)
    groups[physical_name] = faces
    mesh
end

"""
abstract type AbstractChain

# Basic queries

- [`node_coordinates`](@ref)
- [`face_nodes`](@ref)
- [`face_reference_id`](@ref)
- [`reference_faces`](@ref)
- [`periodic_nodes`](@ref)
- [`physical_faces`](@ref)
- [`physical_nodes`](@ref)
- [`outwards_normals`](@ref)

# Basic constructors

- [`chain_from_arrays`](@ref)

"""
abstract type AbstractChain <: GT.AbstractType end

struct GenericChain{A,B,C,D,E,F,G} <: AbstractChain
    node_coordinates::A
    face_nodes::B
    face_reference_id::C
    reference_faces::D
    periodic_nodes::E
    physical_faces::F
    outwards_normals::G
end

function chain_from_arrays(args...)
    GenericChain(args...)
end

"""
"""
function chain_from_arrays(
    node_coordinates,
    face_nodes,
    face_reference_id,
    reference_faces;
    periodic_nodes = eltype(eltype(face_reference_id))[] => eltype(eltype(face_reference_id))[],
    physical_faces = Dict{String,Vector{eltype(eltype(face_reference_id))}}(),
    outwards_normals = nothing
    )
    chain_from_arrays(
            node_coordinates,
            face_nodes,
            face_reference_id,
            reference_faces,
            periodic_nodes,
            physical_faces,
            outwards_normals)
end

num_dims(mesh::AbstractChain) = num_dims(first(reference_faces(mesh)))

function mesh(chain::AbstractChain)
    mesh_from_chain(chain)
end

"""
"""
function mesh_from_chain(chain)
    D = num_dims(chain)
    cell_nodes = face_nodes(chain)
    cell_reference_id = face_reference_id(chain)
    reference_cells = reference_faces(chain)
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
    ref_faces = reference_faces(boundary(ref_cell))
    refid_to_refface = push(ref_faces,reference_cells)
    cell_groups = physical_faces(chain)
    groups = [ typeof(cell_groups)() for d in 0:D]
    groups[end] = cell_groups
    pnodes = periodic_nodes(chain)
    onormals = outwards_normals(chain)
    GT.mesh_from_arrays(
      node_coords,
      face_to_nodes,
      face_to_refid,
      refid_to_refface;
      periodic_nodes = pnodes,
      physical_faces = groups,
      outwards_normals = onormals)
end

# TODO simplexify mesh
"""
"""
function simplexify(geo::AbstractFaceGeometry)
    simplexify_face_geometry(geo)
end

function simplexify_face_geometry(geo)
    if is_unit_simplex(geo)
        simplexify_unit_simplex(geo)
    elseif is_unit_n_cube(geo)
        simplexify_unit_n_cube(geo)
    else
        error("case not implemented")
    end
end

function simplexify_unit_simplex(geo)
    @assert is_unit_simplex(geo)
    refface = lagrange_mesh_face(geo,1)
    mesh_from_reference_face(refface)
end

function simplexify_unit_n_cube(geo)
    @assert is_unit_n_cube(geo)
    D = num_dims(geo)
    if D in (0,1)
        return simplexify_unit_simplex(geo)
    end
    simplex = unit_simplex(Val(D))
    order = 1
    ref_cell = lagrange_mesh_face(simplex,order)
    node_coords = node_coordinates(boundary(geo))
    cell_nodes = simplex_node_ids_n_cube(geo)
    ncells = length(cell_nodes)
    cell_reference_id = fill(Int8(1),ncells)
    reference_cells = [ref_cell]
    chain = chain_from_arrays(
        node_coords,
        cell_nodes,
        cell_reference_id,
        reference_cells,
       )
    mesh = mesh_from_chain(chain)
    mesh_complex = complexify(mesh)
    groups = [Dict{String,Vector{Int32}}() for d in 0:D]
    for d in 0:(D-1)
        sface_to_nodes = face_nodes(mesh_complex,d)
        cface_to_nodes = face_nodes(boundary(geo),d)
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
    physical_faces(mesh_complex) .= groups
    mesh_complex
end

function simplex_node_ids_n_cube(geo)
    D = num_dims(geo)
    # TODO check orientation of nodes
    # this assumes lexicographic ordering
    # for 3d nodes ids are carefully selected
    # such that opposite faces match.
    # This allows one to complexify meshes
    # of ncubes with oriented faces
    if D == 0
        [[1]]
    elseif D == 1
        [[1,2]]
    elseif D == 2
        [[1,2,3],[2,3,4]]
    elseif D ==3
        [[1,2,3,7], [1,2,5,7], [2,3,4,7],
        [2,4,7,8], [2,5,6,7], [2,6,7,8]]
    else
        error("case not implemented")
    end
end

function simplexify_reference_face(ref_face)
    mesh_geom  = simplexify(geometry(ref_face))
    D = num_dims(mesh_geom)
    node_coordinates_geom = node_coordinates(mesh_geom)
    ref_faces_geom = reference_faces(mesh_geom,D)
    face_nodes_geom = face_nodes(mesh_geom,D)
    face_ref_id_geom = face_reference_id(mesh_geom,D)
    nfaces = length(face_ref_id_geom)
    # We need the same order in all directions
    # for this to make sense
    my_order = order(ref_face)
    ref_faces_inter = map(r_geom->lagrange_mesh_face(geometry(r_geom),my_order),ref_faces_geom)
    s_ref = map(ref_faces_geom,ref_faces_inter) do r_geom,r
        m = num_nodes(r)
        n = num_nodes(r_geom)
        x = node_coordinates(r)
        tabulator(r_geom)(value,x)
    end
    node_coordinates_inter = node_coordinates(ref_face)
    node_coordinates_aux = map(xi->map(xii->round(Int,my_order*xii),xi),node_coordinates_inter)
    face_nodes_inter = Vector{Vector{Int}}(undef,nfaces)
    for face in 1:nfaces
        ref_id_geom = face_ref_id_geom[face]
        s = s_ref[ref_id_geom]
        nodes_geom = face_nodes_geom[face]
        nnodes, nnodes_geom = size(s)
        x_mapped = map(1:nnodes) do i
            x = zero(eltype(node_coordinates_inter))
            for k in 1:nnodes_geom
                x += node_coordinates_geom[nodes_geom[k]]*s[i,k]
            end
            map(xi->round(Int,my_order*xi),x)
        end
        my_nodes = indexin(x_mapped,node_coordinates_aux)
        face_nodes_inter[face] = my_nodes
    end
    ref_inter = 
    chain = chain_from_arrays(;
        num_dims=Val(D),
        node_coordinates=node_coordinates_inter,
        face_nodes = face_nodes_inter,
        face_reference_id = face_ref_id_geom,
        reference_faces = ref_faces_inter,
    )
    mesh = mesh_from_chain(chain)
    mesh_complex = complexify(mesh)
    pg = physical_faces(mesh_complex)
    pg .= physical_faces(mesh_geom)
    mesh_complex
end

"""
"""
function cartesian_mesh(
    domain,cells_per_dir;
    boundary=true,
    complexify=true,
    simplexify=false,
    parts_per_dir=nothing,
    parts = parts_per_dir === nothing ? nothing : LinearIndices((prod(parts_per_dir),)),
    partition_strategy = GT.partition_strategy()
    )
    mesh = if boundary
        if simplexify
            structured_simplex_mesh_with_boundary(domain,cells_per_dir)
        else
            cartesian_mesh_with_boundary(domain,cells_per_dir)
        end
    else
        if simplexify
            chain = structured_simplex_chain(domain,cells_per_dir)
        else
            chain = cartesian_chain(domain,cells_per_dir)
        end
        mesh_from_chain(chain)
    end
    if complexify
        mesh = GT.complexify(mesh)
    end
    if parts_per_dir === nothing
        return mesh
    end
    # TODO this can (and should!) be heavily optimized
    @assert partition_strategy !== nothing
    graph_nodes = partition_strategy.graph_nodes
    graph_edges = partition_strategy.graph_edges
    ghost_layers = partition_strategy.ghost_layers
    @assert graph_nodes === :cells
    @assert graph_edges === :nodes
    np = prod(parts_per_dir)
    graph = mesh_graph(mesh;partition_strategy)
    parts_seq = LinearIndices((np,))
    graph_partition = zeros(Int,prod(cells_per_dir))
    for ids in uniform_partition(LinearIndices((np,)),parts_per_dir,cells_per_dir)
        graph_partition[local_to_global(ids)] = local_to_owner(ids)
    end
    partition_mesh(mesh,np;partition_strategy,parts,graph_partition)
end

function bounding_box_from_domain(domain)
    l = length(domain)
    D = div(l,2)
    pmin = SVector(ntuple(d->domain[2*(d-1)+1],Val(D)))
    pmax = SVector(ntuple(d->domain[2*d],Val(D)))
    (pmin,pmax)
end

function domain_from_bounding_box(box)
    l = sum(length,box)
    ntuple(Val(l)) do i
        vector = mod(i-1,2)+1
        component = div(i-1,2)+1
        box[vector][component]
    end
end

function cartesian_mesh_with_boundary(domain,cells_per_dir)
    if any(i->i!=1,cells_per_dir) && any(i->i<2,cells_per_dir)
        error("At least 2 cells in any direction (or 1 cell in all directions)")
    end
    function barrier(
      D,
      cell_to_nodes,
      nnodes,
      d_to_ldface_to_lnodes)

      node_to_n = zeros(Int32,nnodes)
      for nodes in cell_to_nodes
        for node in nodes
          node_to_n[node] += Int32(1)
        end
      end
      J = typeof(JaggedArray(Vector{Int32}[]))
      face_to_nodes = Vector{J}(undef,D)
      groups = [ Dict{String,Vector{Int32}}() for d in 0:(D-1) ]
      ngroups = 0
      for d in 0:(D-1)
        nmax = 2^d
        ldface_to_lnodes = d_to_ldface_to_lnodes[d+1]
        ndfaces = 0
        for nodes in cell_to_nodes
          for (ldface,lnodes) in enumerate(ldface_to_lnodes)
            isboundary = true
            for lnode in lnodes
              node = nodes[lnode]
              if node_to_n[node] > nmax
                isboundary = false
                break
              end
            end
            if isboundary
              ndfaces += 1
            end
          end
        end
        ptrs = zeros(Int32,ndfaces+1)
        for dface in 1:ndfaces
          ptrs[dface+1] += Int32(nmax)
        end
        length_to_ptrs!(ptrs)
        ndata = ptrs[end]-1
        data = zeros(Int32,ndata)
        dface_to_physical_group = zeros(Int32,ndfaces)
        ndfaces = 0
        for nodes in cell_to_nodes
          for (ldface,lnodes) in enumerate(ldface_to_lnodes)
            isboundary = true
            for lnode in lnodes
              node = nodes[lnode]
              if node_to_n[node] > nmax
                isboundary = false
                break
              end
            end
            if isboundary
              ndfaces += 1
              group = ngroups + ldface
              dface_to_physical_group[ndfaces] = group
              p = ptrs[ndfaces]-Int32(1)
              for (i,lnode) in enumerate(lnodes)
                node = nodes[lnode]
                data[p+i] = node
              end
            end
          end
        end
        nldfaces = length(ldface_to_lnodes)
        face_to_nodes[d+1] = JaggedArray(data,ptrs)
        for ldface in 1:nldfaces
          group = ngroups + ldface
          group_name = "$(d)-face-$ldface"
          faces_in_physical_group = findall(g->g==group,dface_to_physical_group)
          groups[d+1][group_name] = faces_in_physical_group
        end
        ngroups += nldfaces
        if d == (D-1)
            groups[d+1]["boundary"] = 1:length(dface_to_physical_group)
        end
      end # d
      ngroups += 1
      groups, face_to_nodes
    end # barrier
    chain = cartesian_chain(domain,cells_per_dir)
    interior_mesh = mesh_from_chain(chain)
    D = num_dims(interior_mesh)
    cell_to_nodes = face_nodes(interior_mesh,D)
    reference_cells = reference_faces(interior_mesh,D)
    node_coords = node_coordinates(interior_mesh)
    ref_cell = first(reference_cells)
    refid_to_refface = reference_faces(boundary(ref_cell))
    nnodes = num_nodes(interior_mesh)
    d_to_ldface_to_lnodes = [face_nodes(boundary(ref_cell),d) for d in 0:(D-1)]
    groups, face_to_nodes = barrier(
      D,
      cell_to_nodes,
      nnodes,
      d_to_ldface_to_lnodes)
    face_to_refid = [ ones(Int8,length(face_to_nodes[d+1]))  for d in 0:(D-1)]
    mesh_face_nodes = push(face_to_nodes,face_nodes(interior_mesh,D))
    mesh_face_reference_id = push(face_to_refid,face_reference_id(interior_mesh,D))
    mesh_reference_faces = push(refid_to_refface,reference_faces(interior_mesh,D))
    mesh_groups = push(groups,physical_faces(interior_mesh,D))
    GT.mesh_from_arrays(
     node_coords,
     mesh_face_nodes,
     mesh_face_reference_id,
     mesh_reference_faces;
     physical_faces=mesh_groups,
    )
end

function cartesian_chain(domain,cells_per_dir)
    box = bounding_box_from_domain(domain)
    D = length(cells_per_dir)
    nodes_per_dir = cells_per_dir .+ 1
    pmin = first(box)
    pmax = last(box)
    extent_per_dir = pmax .- pmin
    h_per_dir = SVector(extent_per_dir ./ cells_per_dir)
    nnodes = prod(nodes_per_dir)
    ncells = prod(cells_per_dir)
    nlnodes = 2^D
    cell_nodes_ptrs = fill(Int32(nlnodes),ncells+1)
    cell_nodes_ptrs[1] = 0
    length_to_ptrs!(cell_nodes_ptrs)
    cell_nodes_data = zeros(Int32,ncells*nlnodes)
    cell_nodes = JaggedArray(cell_nodes_data,cell_nodes_ptrs)
    cell_cis = CartesianIndices(cells_per_dir)
    cell_lis = LinearIndices(cells_per_dir)
    node_cis = CartesianIndices(nodes_per_dir)
    node_lis = LinearIndices(nodes_per_dir)
    lnode_cis = CartesianIndices(ntuple(i->0:1,Val(D)))
    for (cell_li,cell_ci) in enumerate(cell_cis)
        nodes = cell_nodes[cell_li]
        for (lnode_li,lnode_ci) in enumerate(lnode_cis)
            node_ci = CartesianIndex(Tuple(cell_ci) .+ Tuple(lnode_ci))
            node_li = node_lis[node_ci]
            nodes[lnode_li] = node_li
        end
    end
    node_coords = zeros(SVector{D,Float64},nnodes)
    for (node_li,node_ci) in enumerate(node_cis)
        anchor = SVector(Tuple(node_ci) .- 1)
        x = pmin .+ h_per_dir .* anchor
        node_coords[node_li] = x
    end
    cell_reference_id = fill(Int32(1),ncells)
    cell_geometry = unit_n_cube(Val(D))
    order = 1
    ref_cell = lagrange_mesh_face(cell_geometry,order)
    reference_cells = [ref_cell]
    interior_cells = collect(Int32,1:length(cell_nodes))
    groups = Dict(["interior"=>interior_cells,"$D-face-1"=>interior_cells])
    chain = chain_from_arrays(
        node_coords,
        cell_nodes,
        cell_reference_id,
        reference_cells;
        physical_faces=groups,
       )
    chain
end

function structured_simplex_chain(domain,cells_per_dir)
    box = bounding_box_from_domain(domain)
    D = length(cells_per_dir)
    nodes_per_dir = cells_per_dir .+ 1
    pmin = first(box)
    pmax = last(box)
    extent_per_dir = pmax .- pmin
    h_per_dir = SVector(extent_per_dir ./ cells_per_dir)
    nnodes = prod(nodes_per_dir)
    nlnodes = 2^D
    cell_geometry = unit_n_cube(Val(D))
    ref_simplex_mesh = simplexify(cell_geometry)
    rscell_to_lnodes = face_nodes(ref_simplex_mesh,D)
    nrscells = length(rscell_to_lnodes)
    nslnodes = length(first(rscell_to_lnodes))
    ncells = prod(cells_per_dir)*nrscells
    cell_nodes_ptrs = fill(Int32(nslnodes),ncells+1)
    cell_nodes_ptrs[1] = 0
    length_to_ptrs!(cell_nodes_ptrs)
    ndata = cell_nodes_ptrs[end]-1
    cell_nodes_data = zeros(Int32,ndata)
    cell_nodes = JaggedArray(cell_nodes_data,cell_nodes_ptrs)
    cell_cis = CartesianIndices(cells_per_dir)
    cell_lis = LinearIndices(cells_per_dir)
    node_cis = CartesianIndices(nodes_per_dir)
    node_lis = LinearIndices(nodes_per_dir)
    lnode_cis = CartesianIndices(ntuple(i->0:1,Val(D)))
    clnodes = zeros(Int,nlnodes)
    scell = 0
    for (cell_li,cell_ci) in enumerate(cell_cis)
        for (lnode_li,lnode_ci) in enumerate(lnode_cis)
            node_ci = CartesianIndex(Tuple(cell_ci) .+ Tuple(lnode_ci))
            node_li = node_lis[node_ci]
            clnodes[lnode_li] = node_li
        end
        for srcell in 1:nrscells
            scell += 1
            nodes = cell_nodes[scell]
            lnodes = rscell_to_lnodes[srcell]
            for (i,lnode) in enumerate(lnodes)
                nodes[i] = clnodes[lnode]
            end
        end
    end
    node_coords = zeros(SVector{D,Float64},nnodes)
    for (node_li,node_ci) in enumerate(node_cis)
        anchor = SVector(Tuple(node_ci) .- 1)
        x = pmin .+ h_per_dir .* anchor
        node_coords[node_li] = x
    end
    cell_reference_id = fill(Int32(1),ncells)
    order = 1
    reference_cells = reference_faces(ref_simplex_mesh,D)
    interior_cells = collect(Int32,1:length(cell_nodes))
    groups = Dict(["interior"=>interior_cells,"$D-face-1"=>interior_cells])
    chain = chain_from_arrays(
        node_coords,
        cell_nodes,
        cell_reference_id,
        reference_cells;
        physical_faces=groups,
       )
    chain
end

function structured_simplex_mesh_with_boundary(domain,cells_per_dir)
    if any(i->i!=1,cells_per_dir) && any(i->i<2,cells_per_dir)
        error("At least 2 cells in any direction (or 1 cell in all directions)")
    end
    function barrier(
      D,
      cell_to_nodes,
      nnodes,
      d_to_ldface_to_lnodes,
      d_to_ldface_to_sldface_to_lnodes)

        node_to_n = zeros(Int32,nnodes)
        for nodes in cell_to_nodes
            for node in nodes
                node_to_n[node] += Int32(1)
            end
        end
        J = typeof(JaggedArray(Vector{Int32}[]))
        face_to_nodes = Vector{J}(undef,D)
        groups = [ Dict{String,Vector{Int32}}() for d in 0:(D-1) ]
        ngroups = 0
        for d in 0:(D-1)
            nmax = 2^d
            ldface_to_lnodes = d_to_ldface_to_lnodes[d+1]
            ldface_to_sldface_to_lnodes = d_to_ldface_to_sldface_to_lnodes[d+1]
            ndfaces = 0
            for nodes in cell_to_nodes
                for (ldface,lnodes) in enumerate(ldface_to_lnodes)
                    isboundary = true
                    for lnode in lnodes
                        node = nodes[lnode]
                        if node_to_n[node] > nmax
                            isboundary = false
                            break
                        end
                    end
                    if isboundary
                        ndfaces += length(ldface_to_sldface_to_lnodes[ldface])
                    end
                end
            end
            nslnodes = length(ldface_to_sldface_to_lnodes[begin][begin])
            ptrs = zeros(Int32,ndfaces+1)
            for dface in 1:ndfaces
                ptrs[dface+1] += Int32(nslnodes)
            end
            length_to_ptrs!(ptrs)
            ndata = ptrs[end]-1
            data = zeros(Int32,ndata)
            dface_to_physical_group = zeros(Int32,ndfaces)
            ndfaces = 0
            for nodes in cell_to_nodes
                for (ldface,lnodes) in enumerate(ldface_to_lnodes)
                    isboundary = true
                    for lnode in lnodes
                        node = nodes[lnode]
                        if node_to_n[node] > nmax
                            isboundary = false
                            break
                        end
                    end
                    if isboundary
                        group = ngroups + ldface
                        sldface_to_lnodes = ldface_to_sldface_to_lnodes[ldface]
                        nsldfaces = length(sldface_to_lnodes)
                        for sldface in 1:nsldfaces
                            ndfaces += 1
                            dface_to_physical_group[ndfaces] = group
                            p = ptrs[ndfaces]-Int32(1)
                            mylnodes = sldface_to_lnodes[sldface]
                            for (i,lnode) in enumerate(mylnodes)
                                node = nodes[lnode]
                                data[p+i] = node
                            end
                        end
                    end
                end
            end
            nldfaces = length(ldface_to_lnodes)
            face_to_nodes[d+1] = JaggedArray(data,ptrs)
            for ldface in 1:nldfaces
                group = ngroups + ldface
                group_name = "$(d)-face-$ldface"
                faces_in_physical_group = findall(g->g==group,dface_to_physical_group)
                groups[d+1][group_name] = faces_in_physical_group
            end
            ngroups += nldfaces
            if d == (D-1)
                groups[d+1]["boundary"] = 1:length(dface_to_physical_group)
            end
        end # d
        ngroups += 1
        groups, face_to_nodes
    end # barrier
    chain = cartesian_chain(domain,cells_per_dir)
    interior_mesh = mesh_from_chain(chain)
    D = num_dims(interior_mesh)
    cell_to_nodes = face_nodes(interior_mesh,D)
    reference_cells = reference_faces(interior_mesh,D)
    ref_cell = first(reference_cells)
    refid_to_refface = reference_faces(boundary(ref_cell))
    nnodes = num_nodes(interior_mesh)

    cell_geometry = unit_n_cube(Val(D))
    ref_simplex_mesh = simplexify(cell_geometry)
    d_to_ldface_to_sldface_to_lnodes = [
      [ face_nodes(ref_simplex_mesh,d)[physical_faces(ref_simplex_mesh,d)["$d-face-$ldface"]]
      for ldface in 1:num_faces(boundary(ref_cell),d) ] for d in 0:(D-1)]
    d_to_ldface_to_lnodes = [face_nodes(boundary(ref_cell),d) for d in 0:(D-1)]
    groups, face_to_nodes = barrier(
      D,
      cell_to_nodes,
      nnodes,
      d_to_ldface_to_lnodes,
      d_to_ldface_to_sldface_to_lnodes
     )
    simplex_chain = structured_simplex_chain(domain,cells_per_dir)
    node_coords = node_coordinates(simplex_chain)
    face_to_refid = [ ones(Int8,length(face_to_nodes[d+1]))  for d in 0:(D-1)]
    mesh_face_nodes = push(face_to_nodes,face_nodes(simplex_chain))
    mesh_face_reference_id = push(face_to_refid,face_reference_id(simplex_chain))
    mesh_reference_faces = reference_faces(ref_simplex_mesh)
    mesh_groups = push(groups,physical_faces(simplex_chain))
    GT.mesh_from_arrays(
     node_coords,
     mesh_face_nodes,
     mesh_face_reference_id,
     mesh_reference_faces;
     physical_faces=mesh_groups,
    )
end

function visualization_mesh(mesh::AbstractMesh,args...;kwargs...)
    visualization_mesh_from_mesh(mesh,args...;kwargs...)
end

function visualization_mesh_from_mesh(mesh,dim,ids=num_faces(mesh,dim);order=nothing,refinement=nothing)
    function barrier(
            refid_to_tabulation,
            refid_to_scell_to_snodes,
            refid_to_scell_to_srefid,
            refid_to_srefid_to_oid,
            refid_to_srefid_to_vrefface,
            refid_to_snode_to_coords,
            node_to_coords,
            cell_to_nodes,
            cell_to_refid,
            ::Val{Dn}) where Dn

        ncells = length(cell_to_refid)
        nvnodes = 0
        nvcells = 0
        for cell in 1:ncells
            refid = cell_to_refid[cell]
            nvnodes += length(refid_to_snode_to_coords[refid])
            nvcells += length(refid_to_scell_to_srefid[refid])

        end
        nrefids = length(refid_to_srefid_to_oid)
        i_to_oid = reduce(vcat,refid_to_srefid_to_oid)
        i_to_vrefface = reduce(vcat,refid_to_srefid_to_vrefface)
        refid_to_srefid_to_i = Vector{Vector{Int}}(undef,nrefids)
        i = 0
        for refid in 1:nrefids
            srefid_to_oid = refid_to_srefid_to_oid[refid]
            nsrefids = length(srefid_to_oid)
            srefid_to_i = zeros(Int,nsrefids)
            for srefid in 1:nsrefids
                i += 1
                srefid_to_i[srefid] = i
            end
            refid_to_srefid_to_i[refid] = srefid_to_i
        end
        vrefid_to_oid = unique(i_to_oid)
        i_to_vrefid = indexin(i_to_oid,vrefid_to_oid)
        vrefid_to_i = indexin(vrefid_to_oid,i_to_oid)
        vrefid_to_vrefface = i_to_vrefface[vrefid_to_i]
        Tx = SVector{Dn,Float64}
        vnode_to_coords = zeros(Tx,nvnodes)
        vcell_to_vnodes_ptrs = zeros(Int32,nvcells+1)
        vcell_to_vrefid = zeros(Int32,nvcells)
        vcell_to_cell = zeros(Int32,nvcells)
        cell_to_vnodes = fill(0:1,ncells)
        vcell = 0
        for cell in 1:ncells
            refid = cell_to_refid[cell]
            scell_to_snodes = refid_to_scell_to_snodes[refid]
            nscells = length(scell_to_snodes)
            scell_to_srefid = refid_to_scell_to_srefid[refid]
            srefid_to_i = refid_to_srefid_to_i[refid]
            for scell in 1:nscells
                srefid = scell_to_srefid[scell]
                i = srefid_to_i[srefid]
                vrefid = i_to_vrefid[i]
                snodes = scell_to_snodes[scell]
                vcell += 1
                vcell_to_vnodes_ptrs[vcell+1] = length(snodes)
                vcell_to_vrefid[vcell] = vrefid
                vcell_to_cell[vcell] = cell
            end
        end
        length_to_ptrs!(vcell_to_vnodes_ptrs)
        ndata = vcell_to_vnodes_ptrs[end]-1
        vcell_to_vnodes_data = zeros(Int32,ndata)
        vcell = 0
        vnode = 0
        vnode_prev = 1
        for cell in 1:ncells
            refid = cell_to_refid[cell]
            scell_to_snodes = refid_to_scell_to_snodes[refid]
            nscells = length(scell_to_snodes)
            for scell in 1:nscells
                snodes = scell_to_snodes[scell]
                vcell += 1
                p = vcell_to_vnodes_ptrs[vcell]
                for (i,snode) in enumerate(snodes)
                    vcell_to_vnodes_data[p-1+i] = snode + vnode
                end
            end
            tabulation = refid_to_tabulation[refid]
            nsnodes = size(tabulation,1)
            nodes = cell_to_nodes[cell]
            for snode in 1:nsnodes
                y = zero(Tx)
                for (i,node) in enumerate(nodes)
                    coeff = tabulation[snode,i]
                    x = node_to_coords[node]
                    y += coeff*x
                end
                vnode += 1
                vnode_to_coords[vnode] = y
            end
            cell_to_vnodes[cell] = vnode_prev:vnode
            vnode_prev = vnode + 1
        end
        vcell_to_vnodes = JaggedArray(vcell_to_vnodes_data,vcell_to_vnodes_ptrs)
        vchain = chain_from_arrays(
                        vnode_to_coords,
                        vcell_to_vnodes,
                        vcell_to_vrefid,
                        vrefid_to_vrefface)
        vmesh = mesh_from_chain(vchain)
        vglue = (;parent_face=vcell_to_cell,
                 reference_coordinates=refid_to_snode_to_coords,
                 face_fine_nodes = cell_to_vnodes,
                 num_dims=Val(dim))
        vmesh, vglue
    end # barrier
    refid_to_refface = reference_faces(mesh,dim)
    refid_to_refmesh = map(refid_to_refface) do ref_face
        if order === nothing && refinement === nothing
            # Use the given cells as visualization cells
            mesh_from_reference_face(ref_face)
        elseif order !== nothing && refinement === nothing
            # Use cells of given order as visualization cells
            geo = geometry(ref_face)
            ref_face_ho = lagrange_mesh_face(geo,order)
            mesh_from_reference_face(ref_face_ho)
        elseif order === nothing && refinement !== nothing
            # Use linear sub-cells with $refinement per direction per direction
            geom = geometry(ref_face)
            refine_reference_geometry(geom,refinement)
        else
            error("order and refinement kw-arguments can not be given at the same time")
        end
    end
    refid_to_tabulation = map(refid_to_refface,refid_to_refmesh) do refface,refmesh
        x = node_coordinates(refmesh)
        tabulator(refface)(value,x)
    end
    refid_to_scell_to_snodes = map(refmesh->face_nodes(refmesh,dim),refid_to_refmesh)
    refid_to_scell_to_srefid = map(refmesh->face_reference_id(refmesh,dim),refid_to_refmesh)
    refid_to_srefid_to_oid = map(refmesh->map(objectid,reference_faces(refmesh,dim)),refid_to_refmesh)
    refid_to_srefid_to_vrefface = map(refmesh->reference_faces(refmesh,dim),refid_to_refmesh)
    refid_to_snode_to_coords = map(node_coordinates,refid_to_refmesh)
    node_to_coords = node_coordinates(mesh)
    cell_to_nodes = view(face_nodes(mesh,dim),ids)
    cell_to_refid = view(face_reference_id(mesh,dim),ids)
    Dn = num_ambient_dims(mesh)
    barrier(
            refid_to_tabulation,
            refid_to_scell_to_snodes,
            refid_to_scell_to_srefid,
            refid_to_srefid_to_oid,
            refid_to_srefid_to_vrefface,
            refid_to_snode_to_coords,
            node_to_coords,
            cell_to_nodes,
            cell_to_refid,
            Val(Dn))
end

function refine_reference_geometry(geo,resolution)
    function refine_n_cube_aligned(geo,n)
        box = bounding_box(geo)
        domain = domain_from_bounding_box(box)
        D = num_dims(geo)
        cells = ntuple(i->n,Val(D))
        cartesian_mesh(domain,cells)
    end
    function refine_unit_triangle(geo,n)
        # Copyed + adapted from Gridap
        tri_num(n) = n*(n+1)÷2
        v(n,i,j) = tri_num(n) - tri_num(n-i+1) + j
        D = 2
        quad_to_tris = ((1,2,3),(2,4,3))
        quad = CartesianIndices( (0:1,0:1) )
        Tp = SVector{2,Float64}
        n_verts = tri_num(n+1)
        n_cells = tri_num(n)+tri_num(n-1)
        n_verts_x_cell = 3
        X = zeros(Tp,n_verts)
        T = [ zeros(Int,n_verts_x_cell) for i in 1:n_cells ]
        for i in 1:n+1
          for j in 1:n+1-i+1
            vert = v(n+1,i,j)
            X[vert] = SVector((i-1)/n,(j-1)/n)
          end
        end
        for i in 1:n
          for j in 1:n-(i-1)
            verts = ntuple( lv-> v(n+1, (i,j).+quad[lv].I ...), Val{2^D}() )
            cell = v(n,i,j)
            T[cell] .= map(i->verts[i],quad_to_tris[1])
            if (i-1)+(j-1) < n-1
              cell = tri_num(n) + v(n-1,i,j)
              T[cell] .= map(i->verts[i],quad_to_tris[2])
            end
          end
        end
        refface = lagrange_mesh_face(geo,1)
        chain = chain_from_arrays(
                       X,
                       T,
                       fill(1,length(T)),
                       [refface]
                      )
        mesh_from_chain(chain)
    end
    function refine_unit_tet(geo,n)
        # Copyed + adapted from Gridap
        tri_num(n) = n*(n+1)÷2
        tet_num(n) = n*(n+1)*(n+2)÷6
        v(n,i,j) = tri_num(n) - tri_num(n-i+1) + j
        v(n,i,j,k) = tet_num(n) - tet_num(n-i+1) + v(n-i+1,j,k)
        D = 3
        cube_to_tets = ((1,2,3,5),(2,4,3,6),(3,5,7,6),(2,3,5,6),(3,4,7,6),(4,6,7,8))
        cube = CartesianIndices( (0:1,0:1,0:1) )
        n_core_tets = length(cube_to_tets)-2
        Tp = SVector{3,Float64}
        n_verts = tet_num(n+1)
        n_cells = tet_num(n)+n_core_tets*tet_num(n-1)+tet_num(n-2)
        n_verts_x_cell = 4
        X = zeros(Tp,n_verts)
        T = [ zeros(Int,n_verts_x_cell) for i in 1:n_cells ]
        for i in 1:n+1
          for j in 1:n+1-(i-1)
            for k in 1:n+1-(i-1)-(j-1)
              vert = v(n+1,i,j,k)
              X[vert] = SVector((i-1)/n,(j-1)/n,(k-1)/n)
            end
          end
        end
        for i in 1:n
          for j in 1:n-(i-1)
            for k in 1:n-(i-1)-(j-1)
              verts = ntuple( lv-> v(n+1, (i,j,k).+cube[lv].I ...), Val{2^D}() )
              cell = v(n,i,j,k)
              T[cell] .= map(i->verts[i],cube_to_tets[1])
              if (i-1)+(j-1)+(k-1) < n-1
                cell = tet_num(n) + (v(n-1,i,j,k)-1)*n_core_tets
                for t in 1:n_core_tets
                  T[cell+t] .= map(i->verts[i],cube_to_tets[t+1])
                end
              end
              if (i-1)+(j-1)+(k-1) < n-2
                cell = tet_num(n) + n_core_tets*tet_num(n-1) + v(n-2,i,j,k)
                T[cell] .= map(i->verts[i],cube_to_tets[end])
              end
            end
          end
        end
        refface = lagrange_mesh_face(geo,1)
        chain = chain_from_arrays(
                       X,
                       T,
                       fill(1,length(T)),
                       [refface],
                      )
        mesh_from_chain(chain)
    end
    if is_n_cube(geo) && is_axis_aligned(geo)
        refine_n_cube_aligned(geo,resolution)
    elseif is_unit_simplex(geo) && num_dims(geo) == 2
        refine_unit_triangle(geo,resolution)
    elseif is_unit_simplex(geo) && num_dims(geo) == 3
        refine_unit_tet(geo,resolution)
    else
        error("Case not implemented (yet)")
    end
end

function mesh(refface::AbstractMeshFace)
    mesh_from_reference_face(refface)
end

function mesh_from_reference_face(ref_face)
    boundary_mesh = boundary(ref_face)
    D = num_dims(geometry(ref_face))
    nnodes = num_nodes(ref_face)
    face_to_nodes = push(face_nodes(boundary_mesh),[collect(1:nnodes)])
    face_to_refid = push(face_reference_id(boundary_mesh),[1])
    refid_refface = push(reference_faces(boundary_mesh),[ref_face])
    node_to_coords = node_coordinates(ref_face)
    groups = [ Dict{String,Vector{Int32}}() for d in 0:D]
    for d in 0:D
        for face in 1:length(face_to_refid[d+1])
            groups[d+1]["$d-face-$face"] = [face]
        end
    end
    groups[end-1]["boundary"] = 1:length(face_to_refid[D-1+1])
    groups[end]["interior"] = [1]
    mesh = GT.mesh_from_arrays(
                node_to_coords,
                face_to_nodes,
                face_to_refid,
                refid_refface)
end

partition_from_mask(a) = partition_from_mask(identity,a)

function partition_from_mask(f,node_to_mask)
    T = Vector{Int32}
    free_nodes = convert(T,findall(f,node_to_mask))
    dirichlet_nodes = convert(T,findall(i->!f(i),node_to_mask))
    nfree = length(free_nodes)
    ndiri = length(dirichlet_nodes)
    permutation = T(undef,nfree+ndiri)
    permutation[free_nodes] = 1:nfree
    permutation[dirichlet_nodes] = (1:ndiri) .+ nfree
    TwoWayPartition(free_nodes,dirichlet_nodes,permutation)
end

struct TwoWayPartition{A} <: AbstractVector{A}
    first::A
    last::A
    permutation::A
end

permutation(a::TwoWayPartition) = a.permutation
Base.size(a::TwoWayPartition) = (2,)
Base.IndexStyle(::Type{<:TwoWayPartition}) = IndexLinear()
function Base.getindex(a::TwoWayPartition,i::Int)
    @boundscheck @assert i in (1,2)
    if i == 1
        a.first
    else
        a.last
    end
end

function restrict(mesh::AbstractMesh,args...)
    restrict_mesh(mesh,args...)
end

function restrict_mesh(mesh,lnode_to_node,lface_to_face_mesh)
    nnodes = num_nodes(mesh)
    node_to_lnode = zeros(Int32,nnodes)
    node_to_lnode[lnode_to_node] = 1:length(lnode_to_node)
    lnode_to_coords = node_coordinates(mesh)[lnode_to_node]
    lface_to_lnodes_mesh = map(lface_to_face_mesh,face_nodes(mesh)) do lface_to_face,face_to_nodes
        lface_to_nodes = view(face_to_nodes,lface_to_face)
        lface_to_lnodes = JaggedArray(lface_to_nodes)
        f = node->node_to_lnode[node]
        lface_to_lnodes.data .= f.(lface_to_lnodes.data)
        lface_to_lnodes
    end
    lface_to_refid_mesh = map((a,b)->b[a],lface_to_face_mesh,face_reference_id(mesh))
    D = num_dims(mesh)
    lgroups_mesh = map(lface_to_face_mesh,num_faces(mesh),physical_faces(mesh)) do lface_to_face, nfaces, groups
        lgroups = Dict{String,Vector{Int32}}()
        face_to_lface = zeros(Int32,nfaces)
        face_to_lface[lface_to_face] = 1:length(lface_to_face)
        for (k,faces) in groups
            lgroups[k] = filter(i->i!=0,face_to_lface[faces])
        end
        lgroups
    end
    pnode_to_node,pnode_to_master = periodic_nodes(mesh)
    plnode_to_lnode = filter(i->i!=0,node_to_lnode[pnode_to_node])
    plnode_to_lmaster = filter(i->i!=0,node_to_lnode[pnode_to_master])
    if outwards_normals(mesh) !== nothing
        lnormals = outwards_normals(mesh)[lface_to_face_mesh[end]]
    else
        lnormals = nothing
    end

    lmesh = GT.mesh_from_arrays(
        lnode_to_coords,
        lface_to_lnodes_mesh,
        lface_to_refid_mesh,
        reference_faces(mesh);
        physical_faces = lgroups_mesh,
        periodic_nodes = (plnode_to_lnode=>plnode_to_lmaster),
        outwards_normals = lnormals
        )

    lmesh
end

struct PartitionStrategy{A,B} <: GT.AbstractType
    graph_nodes::Symbol
    graph_edges::Symbol
    graph_nodes_dim::A
    graph_edges_dim::B
    ghost_layers::Int
end

function partition_strategy(;
    graph_nodes=:cells,
    graph_edges=:nodes,
    ghost_layers=1,
    graph_nodes_dim=nothing,
    graph_edges_dim=nothing)

    @assert graph_nodes in (:cells,:nodes,:faces)
    @assert graph_edges in (:cells,:nodes,:faces)

    if graph_nodes ∉ (:cells,:nodes) && graph_nodes_dim === nothing
        error("graph_nodes_dim needs to be defined")
    end

    if graph_edges ∉ (:cells,:nodes) && graph_edges_dim === nothing
        error("graph_edges_dim needs to be defined")
    end

    PartitionStrategy(
                      graph_nodes,
                      graph_edges,
                      graph_nodes_dim,
                      graph_edges_dim,
                      ghost_layers)

end

function mesh_graph(mesh::AbstractMesh;
    partition_strategy=GT.partition_strategy())
    graph_nodes = partition_strategy.graph_nodes
    graph_edges = partition_strategy.graph_edges
    graph_nodes_dim = partition_strategy.graph_nodes_dim
    graph_edges_dim = partition_strategy.graph_edges_dim
    function barrier(nnodes,d_to_cell_to_nodes)
        ndata = 0
        for cell_to_nodes in d_to_cell_to_nodes
            ncells = length(cell_to_nodes)
            for cell in 1:ncells
                nodes = cell_to_nodes[cell]
                nlnodes = length(nodes)
                ndata += nlnodes*nlnodes
            end
        end
        I = zeros(Int32,ndata)
        J = zeros(Int32,ndata)
        p = 0
        for cell_to_nodes in d_to_cell_to_nodes
            ncells = length(cell_to_nodes)
            for cell in 1:ncells
                nodes = cell_to_nodes[cell]
                nlnodes = length(nodes)
                for j in 1:nlnodes
                    for i in 1:nlnodes
                        p += 1
                        I[p] = nodes[i]
                        J[p] = nodes[j]
                    end
                end
            end
        end
        V = ones(Int8,ndata)
        g = sparse(I,J,V,nnodes,nnodes)
        fill!(g.nzval,Int8(1))
        g
    end

    if graph_nodes === :nodes
        nnodes = num_nodes(mesh)
        if graph_edges === :cells
            face_nodes_mesh = face_nodes(mesh,num_dims(mesh))
            return barrier(nnodes,[face_nodes_mesh])
        elseif graph_edges === :faces && graph_edges_dim === :all
            face_nodes_mesh = face_nodes(mesh)
            return barrier(nnodes,face_nodes_mesh)
        elseif graph_edges === :faces
            face_nodes_mesh = face_nodes(mesh,graph_edges_dim)
            return barrier(nnodes,[face_nodes_mesh])
        else
            error("case not implemented")
        end
    elseif graph_nodes === :cells
        D = num_dims(mesh)
        ndfaces = num_faces(mesh,D)
        if graph_edges === :nodes
            dface_to_nodes = face_nodes(mesh,D)
            nnodes = num_nodes(mesh)
            node_to_dfaces = generate_face_coboundary(dface_to_nodes,nnodes)
            return barrier(ndfaces,[node_to_dfaces])
        elseif graph_edges === :faces && graph_edges_dim !== :all
            topo = topology(mesh)
            node_to_dfaces = face_incidence(topo,graph_edges_dim,D)
            return barrier(ndfaces,[node_to_dfaces])
        else
            error("case not implemented")
        end
    else
        error("case not implemented")
    end
end

function simplexify(mesh::AbstractMesh;glue=Val(false))

    # TODO add attributes to mesh to figure out if we really
    # need to simplexify and/or complexify
    #mesh, = complexify(mesh)

    D = num_dims(mesh)
    refid_to_refcell = reference_faces(mesh,D)
    #@assert length(refid_to_refcell) == 1
    refcell = first(refid_to_refcell)
    #@assert order(refcell) == 1
    ltmesh = simplexify(geometry(refcell))

    ltcell_to_lnodes = face_nodes(ltmesh,D)
    cell_to_nodes = face_nodes(mesh,D)

    ncells = length(cell_to_nodes)
    nltcells = length(ltcell_to_lnodes)
    ntcells = ncells * nltcells
    tcell_to_nodes_ptrs = zeros(Int32,ntcells+1)
    tcell_to_cell = zeros(Int32,ntcells)

    tcell = 0
    for cell in 1:ncells
        for lnodes in ltcell_to_lnodes
            tcell +=1
            tcell_to_nodes_ptrs[tcell+1] = length(lnodes)
            tcell_to_cell[tcell] = cell
        end
    end
    length_to_ptrs!(tcell_to_nodes_ptrs)
    ndata = tcell_to_nodes_ptrs[end]-1
    T = eltype(eltype(cell_to_nodes))
    tcell_to_nodes_data = zeros(T,ndata)
    k = 1
    for cell in 1:ncells
        nodes = cell_to_nodes[cell]
        for lnodes in ltcell_to_lnodes
            for lnode in lnodes
                node = nodes[lnode]
                tcell_to_nodes_data[k] = node
                k += 1
            end
        end
    end
    tcell_to_nodes = JaggedArray(tcell_to_nodes_data,tcell_to_nodes_ptrs)
    reftcell = first(reference_faces(ltmesh,D))
    tcell_to_refid = fill(Int8(1),ntcells)
    refid_to_reftcell = [reftcell]
    node_to_coords = node_coordinates(mesh)

    tchain = chain_from_arrays(
        node_to_coords,
        tcell_to_nodes,
        tcell_to_refid,
        refid_to_reftcell,)
    tmesh = complexify(mesh_from_chain(tchain))

    topo = topology(mesh)
    ttopo = topology(tmesh)
    ltopo = topology(ltmesh)
    d_to_tface_to_face = map(0:(D-1)) do d
        cell_to_faces = face_incidence(topo,D,d)
        tcell_to_tfaces = JaggedArray(face_incidence(ttopo,D,d))
        lface_to_lnodes = face_nodes(boundary(refcell),d)
        ltface_to_ltnodes = face_nodes(boundary(reftcell),d)
        if length(tcell_to_tfaces.data) != 0
            ntfaces = maximum(tcell_to_tfaces.data)
        else
            ntfaces = 0
        end
        tface_to_face = simplexify_generate_tface_to_face(
                                                          cell_to_faces,
                                                          tcell_to_tfaces,
                                                          ltcell_to_lnodes,
                                                          ltface_to_ltnodes,
                                                          lface_to_lnodes,
                                                          ntfaces)

    end
    push!(d_to_tface_to_face,tcell_to_cell)
    for d in 0:D
        groups = physical_faces(mesh,d)
        tgroups = physical_faces(tmesh,d)
        nfaces = num_faces(mesh,d)
        ntfaces = num_faces(tmesh,d)
        face_to_mask = fill(false,nfaces)
        tface_to_mask = fill(false,ntfaces)
        tface_to_face = d_to_tface_to_face[d+1]
        for (k,v) in groups
            fill!(face_to_mask,false)
            fill!(tface_to_mask,false)
            face_to_mask[v] .= true
            for tface in 1:ntfaces
                face = tface_to_face[tface]
                if face != 0
                    tface_to_mask[tface] = face_to_mask[face]
                end
            end
            tgroups[k] = findall(tface_to_mask)
        end
    end
    if val_parameter(glue)
        tmesh, d_to_tface_to_face
    else
        tmesh
    end
end

function simplexify_generate_tface_to_face(
        cell_to_faces,
        tcell_to_tfaces,
        ltcell_to_lnodes,
        ltface_to_ltnodes,
        lface_to_lnodes,
        ntfaces)

    nltcells = length(ltcell_to_lnodes)
    nltfaces = length(ltface_to_ltnodes)
    nlfaces = length(lface_to_lnodes)
    ltcell_ltface_to_lface = [ zeros(Int32,nltfaces) for ltcell in 1:nltcells  ]
    for ltcell in 1:nltcells
        ltnode_to_lnode = ltcell_to_lnodes[ltcell]
        for ltface in 1:nltfaces
            ltnodes = ltface_to_ltnodes[ltface]
            for lface in 1:nlfaces
                lnodes = lface_to_lnodes[lface]
                allin = true
                for ltnode in ltnodes
                    lnode = ltnode_to_lnode[ltnode]
                    if !(lnode in lnodes)
                        allin = false
                        break
                    end
                end
                if allin
                    ltcell_ltface_to_lface[ltcell][ltface] = lface
                    break
                end
            end
        end
    end
    ltcell_to_lfaces = ltcell_ltface_to_lface
    tface_to_face = zeros(Int32,ntfaces)
    ncells = length(cell_to_faces)
    nltcells = length(ltcell_to_lfaces)
    tcell = 1
    for cell in 1:ncells
        faces = cell_to_faces[cell]
        for ltcell in 1:nltcells
            tfaces = tcell_to_tfaces[tcell]
            nltfaces = length(tfaces)
            ltface_to_lface = ltcell_to_lfaces[ltcell]
            for ltface in 1:nltfaces
                tface = tfaces[ltface]
                lface = ltface_to_lface[ltface]
                if lface != 0
                    face = faces[lface]
                    tface_to_face[tface] = face
                end
            end
            tcell += 1
        end
    end
    tface_to_face
end


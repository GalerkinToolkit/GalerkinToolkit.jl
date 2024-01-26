module GalerkinToolkit

using WriteVTK
using FastGaussQuadrature
using StaticArrays
using LinearAlgebra
using ForwardDiff
using Gmsh
using PartitionedArrays
using Combinatorics
using SparseArrays

abstract type GalerkinToolkitDataType end
function Base.show(io::IO,data::GalerkinToolkitDataType)
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
num_dims(a) = val_parameter(a.num_dims)
node_coordinates(a) = a.node_coordinates
reference_faces(a) = a.reference_faces
face_nodes(a) = a.face_nodes
face_incidence(a) = a.face_incidence
face_reference_id(a) = a.face_reference_id
physical_groups(a) = a.physical_groups
has_physical_groups(a) = hasproperty(a,:physical_groups) && a.physical_groups !== nothing
periodic_nodes(a) = a.periodic_nodes
has_periodic_nodes(a) = hasproperty(a,:periodic_nodes) && a.periodic_nodes !== nothing
geometry(a) = a.geometry
boundary(a) = a.boundary
is_n_cube(a) = hasproperty(a,:is_n_cube) ? val_parameter(a.is_n_cube) : false
is_simplex(a) = hasproperty(a,:is_simplex) ? val_parameter(a.is_simplex) : false
is_axis_aligned(a) = a.is_axis_aligned
bounding_box(a) = a.bounding_box
vertex_permutations(a) = a.vertex_permutations
face_own_dofs(a) = a.face_own_dofs
face_own_dof_permutations(a) = a.face_own_dof_permutations
node_to_dofs(a) = a.node_to_dofs
dof_to_node(a) = a.dof_to_node
dof_to_index(a) = a.dof_to_index
num_dofs(a) = a.num_dofs
coordinates(a) = a.coordinates
weights(a) = a.weights
order_per_dir(a) = a.order_per_dir
monomial_exponents(a) = a.monomial_exponents
lib_to_user_nodes(a) = a.lib_to_user_nodes
interior_nodes(a) = a.interior_nodes
face_permutation_ids(a) = a.face_permutation_ids
face_permutation_ids(a,m,n) = face_permutation_ids(a)[m+1,n+1]
local_nodes(a) = a.local_nodes
local_node_colors(a) = a.local_node_colors
real_type(a) = a.real_type
int_type(a) = a.int_type
outwards_normals(a) = a.outwards_normals

reference_faces(a,d) = reference_faces(a)[val_parameter(d)+1]
face_nodes(a,d) = face_nodes(a)[val_parameter(d)+1]
face_incidence(a,d1,d2) = face_incidence(a)[val_parameter(d1)+1,val_parameter(d2)+1]
face_reference_id(a,d) = face_reference_id(a)[val_parameter(d)+1]
num_faces(a) = map(length,face_reference_id(a))
num_faces(a,d) = length(face_reference_id(a,d))
physical_groups(a,d) = physical_groups(a)[val_parameter(d)+1]
num_nodes(a) = length(node_coordinates(a))
num_ambient_dims(a) = length(eltype(node_coordinates(a)))
function face_offsets(a)
    D = num_dims(a)
    offsets = zeros(Int,D+1)
    for d in 1:D
        offsets[d+1] = offsets[d] + num_faces(a,d-1)
    end
    offsets
end
function face_dim(a,d)
    n = num_faces(a,d)
    fill(d,n)
end
function face_dim(a)
    D = num_dims(a)
    reduce(vcat,map(d->face_dim(a,d),0:D))
end

abstract type AbstractFaceGeometry <: GalerkinToolkitDataType end

struct ExtrusionPolytope{D,Tv,Ti} <: AbstractFaceGeometry
    extrusion::NTuple{D,Bool}
    bounding_box::NTuple{2,SVector{D,Tv}}
    int_type::Type{Ti}
end

function unit_simplex(num_dims;real_type=Float64,int_type=Int)
    D = val_parameter(num_dims)
    extrusion = ntuple(i->false,Val(D))
    Tv = real_type
    p0 = ntuple(i->zero(real_type),Val(D)) |> SVector{D,Tv}
    p1 = ntuple(i->one(real_type),Val(D)) |> SVector{D,Tv}
    bounding_box = (p0,p1)
    ExtrusionPolytope(extrusion,bounding_box,int_type)
end

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

function is_unit_n_cube(geo)
    is_n_cube(geo) && is_unitary(geo)
end

function is_unit_simplex(geo)
    is_simplex(geo) && is_unitary(geo)
end

function is_unitary(geom)
    ! is_axis_aligned(geom) && return false
    my_bounding_box = bounding_box(geom)
    all(i->i==0,first(my_bounding_box)) && all(i->i==1,last(my_bounding_box))
end

struct GenericCuadrature{A,B} <: GalerkinToolkitDataType
    coordinates::A
    weights::B
end
struct Cuadrature{D,T} <: GalerkinToolkitDataType
    coordinates::Vector{SVector{D,T}}
    weights::Vector{T}
end
function quadrature(coordinates,weights)
    GenericCuadrature(coordinates,weights)
end
function quadrature(coordinates::Vector{SVector{D,T}},weights::Vector{T}) where {D,T}
    Cuadrature(coordinates,weights)
end

function duffy_quadrature(geo,degree)
    @assert is_unit_simplex(geo)
    D = num_dims(geo)
    if D == 0
        x = zeros(SVector{0,real_type},1)
        w = ones(real_type,1)
        return quadrature(x,w)
    end
    function map_to(a,b,(points,weights))
      points_ab = similar(points)
      weights_ab = similar(weights)
      points_ab .= 0.5*(b-a)*points .+ 0.5*(a+b)
      weights_ab .= 0.5*(b-a)*weights
      (points_ab, weights_ab)
    end
    function duffy_map(q)
        D = length(q)
        a = 1.0
        m = ntuple(Val(D)) do i
            if i == 1
                q[i]
            else
                a *= (1-q[i-1])
                a*q[i]
            end
        end
        typeof(q)(m)
    end
    n = ceil(Int, (degree + 1.0) / 2.0 )
    beta = 0
    dim_to_quad_1d = map(1:(D-1)) do d
        alpha = (D-1)-(d-1)
        map_to(0,1,gaussjacobi(n,alpha,beta))
    end
    quad_1d = map_to(0,1,gausslegendre(n))
    push!(dim_to_quad_1d,quad_1d)
    coords_per_dir = map(first,dim_to_quad_1d)
    weights_per_dir =  map(last,dim_to_quad_1d)
    a = 0.5
    for d in (D-1):-1:1
        ws_1d = weights_per_dir[d]
        ws_1d[:] *= a
        a *= 0.5
    end
    m = prod(map(length,weights_per_dir))
    Tv = real_type(geo)
    w = zeros(Tv,m)
    x = zeros(SVector{D,Tv},m)
    tensor_product!(identity,x,coords_per_dir)
    tensor_product!(prod,w,weights_per_dir)
    x .= duffy_map.(x)
    quadrature(x,w)
end

function tensor_product_quadrature(geo,degree_per_dir)
    @assert is_n_cube(geo)
    @assert is_axis_aligned(geo)
    my_bounding_box = bounding_box(geo)
    D = num_dims(geo)
    limits_per_dir = ntuple(i->(my_bounding_box[1][i],my_bounding_box[2][i]),Val(D))
    n_per_dir = map(d->ceil(Int,(d+1)/2),degree_per_dir)
    function quadrature_1d(n,limits)
        x,w = FastGaussQuadrature.gausslegendre(n)
        a,b = limits
        x .= (0.5*(b-a)) .*x .+ (0.5*(b+a))
        w .*= 0.5*(b-a)
        quadrature(x,w)
    end
    quad_per_dir = map(quadrature_1d,n_per_dir,limits_per_dir)
    coords_per_dir = map(coordinates,quad_per_dir)
    weights_per_dir = map(weights,quad_per_dir)
    m = prod(map(length,weights_per_dir))
    Tv = real_type(geo)
    w = zeros(Tv,m)
    x = zeros(SVector{D,Tv},m)
    tensor_product!(identity,x,coords_per_dir)
    tensor_product!(prod,w,weights_per_dir)
    quadrature(x,w)
end

function tensor_product!(f,result,values_per_dir)
    shape = Tuple(map(length,values_per_dir))
    cis = CartesianIndices(shape)
    lis = LinearIndices(cis)
    for ci in cis
        li = lis[ci]
        result[li] = f(map((q,i)->q[i],values_per_dir,Tuple(ci)))
    end
    result
end

function repeat_per_dir(geo,a)
    D = num_dims(geo)
    ntuple(i->a,Val(D))
end
repeat_per_dir(geo,a::NTuple) = a

function default_quadrature(geo,degree)
    if is_n_cube(geo) && is_axis_aligned(geo)
        D = num_dims(geo)
        degree_per_dir = repeat_per_dir(geo,degree)
        tensor_product_quadrature(geo,degree_per_dir)
    elseif is_unit_simplex(geo)
        duffy_quadrature(geo,degree)
    else
        error("Not implemented")
    end
end

abstract type AbstractMeshFace <: GalerkinToolkitDataType end

num_dims(f::AbstractMeshFace) = num_dims(geometry(f))

abstract type AbstractLagrangeFE <: AbstractMeshFace end

struct GenericLagrangeFE{A,B,C,D} <: AbstractLagrangeFE
    geometry::A
    order_per_dir::B
    space::Symbol
    lib_to_user_nodes::C
    major::Symbol
    shape::D
end

struct ScalarShape end
const SCALAR_SHAPE = ScalarShape()

lagrangian_fe(args...) = GenericLagrangeFE(args...)

function lagrangian_fe(geometry,order;
        space = default_space(geometry),
        lib_to_user_nodes = int_type(geometry)[],
        major = :component,
        shape = SCALAR_SHAPE)
    D = num_dims(geometry)
    order_per_dir = repeat_per_dir(geometry,order)
    lagrangian_fe(
               geometry,
               order_per_dir,
               space,
               lib_to_user_nodes,
               major,
               shape)
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

order(fe::AbstractLagrangeFE) = maximum(order_per_dir(fe),init=0)

function lib_to_user_nodes(fe::AbstractLagrangeFE)
    if length(fe.lib_to_user_nodes) == 0
        nnodes = num_nodes(fe)
        Ti = int_type(geometry(fe))
        collect(Ti.(1:nnodes))
    else
        fe.lib_to_user_nodes
    end
end

function monomial_exponents(fe::AbstractLagrangeFE)
    monomial_exponents_from_space(fe.space,fe.order_per_dir,fe.geometry |> int_type)
end

function node_coordinates(fe::AbstractLagrangeFE)
    @assert fe |> geometry |> is_unitary
    mexps = monomial_exponents(fe)
    node_coordinates_from_monomials_exponents(mexps,fe.order_per_dir,fe.geometry |> real_type)
end

function node_coordinates_from_monomials_exponents(monomial_exponents,order_per_dir,real_type)
    D = length(order_per_dir)
    Tv = real_type
    if length(monomial_exponents) == 0
        return  
    end
    node_coordinates = map(monomial_exponents) do exponent
        map(exponent,order_per_dir) do e,order
            if order != 0
                real_type(e/order)
            else
                real_type(e)
            end
        end |> SVector{D,Tv}
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

function monomial_exponents_from_filter(f,order_per_dir,int_type)
    Ti = int_type
    terms_per_dir = Tuple(map(d->d+1,order_per_dir))
    D = length(terms_per_dir)
    cis = CartesianIndices(terms_per_dir)
    m = count(ci->f(SVector{D,Ti}(Tuple(ci) .- 1),order_per_dir),cis)
    result = zeros(SVector{D,int_type},m)
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

function tensor_basis(fe::AbstractLagrangeFE)
    Tv = fe |> geometry |> real_type
    if fe.shape == SCALAR_SHAPE
        return Tv(1)
    else
        cis = CartesianIndices(val_parameter(fe.shape))
        l = prod(val_parameter(fe.shape))
        init_tensor = SArray{Tuple{val_parameter(fe.shape)...},Tv}
        cis_flat = cis[:]
        return map(cis_flat) do ci
            init_tensor(ntuple(j->cis[j]==ci ? 1 : 0 ,Val(l)))
        end
    end
end

function primal_basis(fe::AbstractLagrangeFE)
    scalar_basis = map(e->(x-> prod(x.^e)),fe|>monomial_exponents)
    if fe.shape == SCALAR_SHAPE
        return scalar_basis
    else
        primal_nested = map(scalar_basis) do monomial
            map(tensor_basis(fe))  do e
                x -> monomial(x)*e
            end
        end
        return reduce(vcat,primal_nested)
    end
end

function dual_basis(fe::AbstractLagrangeFE)
    node_coordinates_reffe = fe|>node_coordinates
    scalar_basis = map(x->(f->f(x)),node_coordinates_reffe)
    if fe.shape == SCALAR_SHAPE
        return scalar_basis
    else
        if fe.major === :component
            dual_nested = map(node_coordinates_reffe) do x
                map(tensor_basis(fe)) do e
                    f->inner(e,f(x))
                end
            end
        elseif fe.major === :node
            dual_nested = map(tensor_basis(fe)) do e
                map(node_coordinates_reffe) do x
                    f->inner(e,f(x))
                end
            end
        else
            error("Not Implemented")
        end
        return reduce(vcat,dual_nested)
    end
end

inner(a,b) = sum(map(*,a,b))
value(f,x) = f(x)

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

abstract type AbstractFEMesh <: GalerkinToolkitDataType end

struct GenericFEMesh{A,B,C,D,E,F,G} <: AbstractFEMesh
    node_coordinates::A
    face_nodes::B
    face_reference_id::C
    reference_faces::D
    periodic_nodes::E
    physical_groups::F
    outwards_normals::G
end

function fe_mesh(args...)
    GenericFEMesh(args...)
end

function fe_mesh(
    node_coordinates,
    face_nodes,
    face_reference_id,
    reference_faces;
    periodic_nodes = eltype(eltype(face_reference_id))[],
    physical_groups = map(i->Dict{String,Vector{eltype(eltype(face_reference_id))}}(),face_reference_id),
    outwards_normals = nothing
    )
    fe_mesh(
            node_coordinates,
            face_nodes,
            face_reference_id,
            reference_faces,
            periodic_nodes,
            physical_groups,
            outwards_normals)
end

num_dims(mesh::AbstractFEMesh) = length(reference_faces(mesh))-1

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
    mesh = fe_mesh(
            my_node_to_coords,
            my_face_nodes,
            my_face_reference_id,
            my_reference_faces;
            physical_groups = my_groups,
            periodic_nodes,)

    if complexify
        mesh, _ = GalerkinToolkit.complexify(mesh)
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
    lagrangian_fe(geom,order;lib_to_user_nodes=lib_to_gmsh)
end

function vtk_points(mesh)
    function barrirer(coords)
        nnodes = length(coords)
        points = zeros(3,nnodes)
        for node in 1:nnodes
            coord = coords[node]
            for i in 1:length(coord)
                points[i,node] = coord[i]
            end
        end
        points
    end
    coords = node_coordinates(mesh)
    barrirer(coords)
end

function vtk_cells(mesh,d)
    function barrirer(face_to_refid,face_to_nodes,refid_mesh_cell)
        cells = map(face_to_refid,face_to_nodes) do refid, nodes
            mesh_cell = refid_mesh_cell[refid]
            if mesh_cell === nothing
                msg = """
                Not enough information to visualize this mesh via vtk:
                vtk_mesh_cell returns nothing for the reference face in position $refid in dimension $d.
                """
                error(msg)
            end
            mesh_cell(nodes)
        end
        cells
    end
    face_to_nodes = face_nodes(mesh,d)
    face_to_refid = face_reference_id(mesh,d)
    refid_refface = reference_faces(mesh,d)
    refid_mesh_cell = map(vtk_mesh_cell,refid_refface)
    barrirer(face_to_refid,face_to_nodes,refid_mesh_cell)
end

"""
    args = vtk_args(mesh,d)

Return the arguments `args` to be passed in final position
to functions like `WriteVTK.vtk_grid`.
"""
function vtk_args(mesh,d)
    points = vtk_points(mesh)
    cells = vtk_cells(mesh,d)
    points, cells
end

function vtk_args(mesh)
    points = vtk_points(mesh)
    D = num_dims(mesh)
    allcells = [vtk_cells(mesh,d) for d in 0:D if num_faces(mesh,d) != 0]
    cells = reduce(vcat,allcells)
    points, cells
end

function vtk_physical_groups!(vtk,mesh,d;physical_groups=physical_groups(mesh,d))
    ndfaces = num_faces(mesh,d)
    for group in physical_groups
        name,faces = group
        face_mask = zeros(Int,ndfaces)
        face_mask[faces] .= 1
        vtk[name,WriteVTK.VTKCellData()] = face_mask
    end
    vtk
end

function vtk_physical_groups!(vtk,mesh;physical_groups=physical_groups(mesh))
    nfaces = sum(num_faces(mesh))
    offsets = face_offsets(mesh)
    D = num_dims(mesh)
    data = Dict{String,Vector{Int}}()
    for d in 0:D
        for group in physical_groups[d+1]
            name, = group
            if !haskey(data,name)
                face_mask = zeros(Int,nfaces)
                data[name] = face_mask
            end
        end
    end
    for d in 0:D
        for group in physical_groups[d+1]
            offset = offsets[d+1]
            name,faces = group
            face_mask = data[name]
            face_mask[faces.+offset] .= 1
        end
    end
    for (name,face_mask) in data
        vtk[name,WriteVTK.VTKCellData()] = face_mask
    end
    vtk
end

function vtk_mesh_cell(ref_face)
    geom = geometry(ref_face)
    d = num_dims(geom)
    nnodes = num_nodes(ref_face)
    lib_to_user = lib_to_user_nodes(ref_face)
    if d == 0 && nnodes == 1
        cell_type = WriteVTK.VTKCellTypes.VTK_VERTEX
        vtk_to_lib = [1]
    elseif d == 1 && (is_simplex(geom) || is_n_cube(geom)) && nnodes == 2
        cell_type = WriteVTK.VTKCellTypes.VTK_LINE
        vtk_to_lib = [1,2]
    elseif d == 1 && (is_simplex(geom) || is_n_cube(geom)) && nnodes == 3
        cell_type = WriteVTK.VTKCellTypes.VTK_QUADRATIC_EDGE
        vtk_to_lib = [1,3,2]
    elseif d == 2 && is_n_cube(geom) && nnodes == 4
        cell_type = WriteVTK.VTKCellTypes.VTK_QUAD
        vtk_to_lib = [1,2,4,3]
    elseif d == 2 && is_simplex(geom) && nnodes == 3
        cell_type = WriteVTK.VTKCellTypes.VTK_TRIANGLE
        vtk_to_lib = [1,2,3]
    elseif d == 2 && is_simplex(geom) && nnodes == 6
        cell_type = WriteVTK.VTKCellTypes.VTK_QUADRATIC_TRIANGLE
        vtk_to_lib = [1,3,6,2,5,4]
    elseif d == 3 && is_n_cube(geom) && nnodes == 8
        cell_type = WriteVTK.VTKCellTypes.VTK_HEXAHEDRON
        vtk_to_lib = [1,2,4,3,5,6,8,7]
    elseif d == 3 && is_simplex(geom) && nnodes == 4
        cell_type = WriteVTK.VTKCellTypes.VTK_TETRA
        vtk_to_lib = [1,2,3,4]
    elseif d == 3 && is_simplex(geom) && nnodes == 10
        cell_type = WriteVTK.VTKCellTypes.VTK_QUADRATIC_TETRA
        vtk_to_lib = [1,3,6,10,2,5,4,7,8,9]
    else
        return nothing
    end
    nodes -> WriteVTK.MeshCell(cell_type,(nodes[lib_to_user])[vtk_to_lib])
end

# TODO use the same convention than in Gridap
# allow the user to customize the boundary object ids
function boundary(geom::ExtrusionPolytope{0})
    nothing
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
    fe = lagrangian_fe(vertex,order)
    node_coordinates = SVector{1,Tv}[(0,),(1,)]
    face_nodes = [Vector{Ti}[[1],[2]]]
    face_reference_id = [Ti[1,1]]
    reference_faces = ([fe],)
    outwards_normals = SVector{1,Tv}[(-1,),(1,)]
    fe_mesh(
            node_coordinates,
            face_nodes,
            face_reference_id,
            reference_faces;
            outwards_normals
           )
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
    fe0 = lagrangian_fe(g0,order)
    fe1 = lagrangian_fe(g1,order)
    reference_faces = ([fe0],[fe1])
    fe_mesh(
            node_coordinates,
            face_nodes,
            face_reference_id,
            reference_faces;
            outwards_normals
           )
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
    fe0 = lagrangian_fe(g0,order)
    fe1 = lagrangian_fe(g1,order)
    fe2 = lagrangian_fe(g2,order)
    reference_faces = ([fe0,],[fe1],[fe2])
    fe_mesh(
            node_coordinates,
            face_nodes,
            face_reference_id,
            reference_faces;
            outwards_normals
           )
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
        map(r->lagrangian_fe(geometry(r),order_inter),ref_faces_geom_d)
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
    fe_mesh(
        node_coordinates_inter,
        face_nodes_inter,
        face_ref_id_geom,
        ref_faces)
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
        for d in 0:D
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
    ref_face = lagrangian_fe(geo,order)
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
    function compute_volume(vertex_coords)
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
    refvol = compute_volume(vertex_coords)
    perm_vertex_coords = similar(vertex_coords)
    for permutation in permutations
        for (j,cj) in enumerate(permutation)
          perm_vertex_coords[j] = vertex_coords[cj]
        end
        vol2 = compute_volume(perm_vertex_coords)
        if (refvol + vol2) ≈ (2*refvol)
            push!(admissible_permutations,permutation)
        end
    end
    admissible_permutations
end

function interior_node_permutations(fe::AbstractMeshFace)
    interior_node_permutations_from_mesh_face(fe)
end

function interior_node_permutations_from_mesh_face(refface)
    interior_ho_nodes = interior_nodes(refface)
    ho_nodes_coordinates = node_coordinates(refface)
    geo = geometry(refface)
    vertex_perms = vertex_permutations(geo)
    if length(interior_ho_nodes) == 0
        return map(i->Int[],vertex_perms)
    end
    if length(vertex_perms) == 1
        return map(i->collect(1:length(interior_ho_nodes)),vertex_perms)
    end
    geo_mesh = boundary(geo)
    vertex_to_geo_nodes = face_nodes(geo_mesh,0)
    vertex_to_geo_node = map(first,vertex_to_geo_nodes)
    ref_face = lagrangian_reference_face(geo)
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
            if pnode != nothing
               node_to_pnode[iq] = pnode
            end
        end
        node_perms[iperm] = node_to_pnode
    end
    node_perms
end

abstract type AbstractMeshTopology <: GalerkinToolkitDataType end

struct GenericMeshTopology{A,B,C,D} <: AbstractMeshTopology
    face_incidence::A
    face_reference_id::B
    face_permutation_ids::C
    reference_faces::D
end

function mesh_topology(args...)
    GenericMeshTopology(args...)
end

function topology(mesh::AbstractFEMesh)
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

function same_valid_ids(a,b)
  function is_subset(a,b)
    for i in 1:length(a)
      v = a[i]
      if v == INVALID_ID
        continue
      end
      c = find_eq(v,b)
      if c == false; return false; end
    end
    return true
  end
  function find_eq(v,b)
    for vs in b
      if v == vs
        return true
      end
    end
    return false
  end
  c = is_subset(a,b)
  if c == false; return false; end
  c = is_subset(b,a)
  if c == false; return false; end
  return true
end

## TODO AbstractFaceTopology <: AbstractMeshTopology
abstract type AbstractFaceTopology <: GalerkinToolkitDataType end

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

function complexify(mesh::AbstractFEMesh)
    complexify_mesh(mesh)
end

function complexify_mesh(mesh)
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
    old_physical_groups = physical_groups(mesh)
    new_physical_groups = [ Dict{String,Vector{Int32}}() for d in 0:D] # TODO hardcoded
    for d in 0:D
        old_groups = old_physical_groups[d+1]
        for (group_name,old_group_faces) in old_groups
            new_group_faces = similar(old_group_faces)
            new_group_faces .= old_to_new[d+1][old_group_faces]
            new_physical_groups[d+1][group_name] = new_group_faces
        end
    end
    new_mesh = fe_mesh(
            node_to_coords,
            newface_nodes,
            newface_refid,
            Tuple(newreffaces);
            physical_groups = new_physical_groups,
            periodic_nodes = periodic_nodes(mesh),
            outwards_normals = outwards_normals(mesh)
           )
    new_mesh, old_to_new
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
                    @boundscheck @assert ldface2 != INVALID_ID
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

abstract type AbstractFEChain <: GalerkinToolkitDataType end

struct GenericFEChain{A,B,C,D,E,F,G} <: AbstractFEChain
    node_coordinates::A
    face_nodes::B
    face_reference_id::C
    reference_faces::D
    periodic_nodes::E
    physical_groups::F
    outwards_normals::G
end

function fe_chain(args...)
    GenericFEChain(args...)
end

function fe_chain(
    node_coordinates,
    face_nodes,
    face_reference_id,
    reference_faces;
    periodic_nodes = eltype(eltype(face_reference_id))[],
    physical_groups = map(i->Dict{String,Vector{eltype(eltype(face_reference_id))}}(),face_reference_id),
    outwards_normals = nothing
    )
    fe_chain(
            node_coordinates,
            face_nodes,
            face_reference_id,
            reference_faces,
            periodic_nodes,
            physical_groups,
            outwards_normals)
end

num_dims(mesh::AbstractFEChain) = num_dims(first(reference_faces(mesh)))

function fe_mesh(chain::AbstractFEChain)
    mesh_from_chain(chain)
end

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
    cell_groups = physical_groups(chain)
    groups = [ typeof(cell_groups)() for d in 0:D]
    groups[end] = cell_groups
    pnodes = periodic_nodes(chain)
    onormals = outwards_normals(chain)
    fe_mesh(
      node_coords,
      face_to_nodes,
      face_to_refid,
      refid_to_refface,
      pnodes,
      groups,
      onormals)
end

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
    refface = lagrangian_reference_face(geo)
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
    ref_cell = lagrangian_fe(simplex,order)
    node_coords = node_coordinates(boundary(geo))
    cell_nodes = simplex_node_ids_n_cube(geo)
    ncells = length(cell_nodes)
    cell_reference_id = fill(Int8(1),ncells)
    reference_cells = [ref_cell]
    chain = fe_chain(
        node_coords,
        cell_nodes,
        cell_reference_id,
        reference_cells,
       )
    mesh = mesh_from_chain(chain)
    mesh_complex, = complexify_mesh(mesh)
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
    physical_groups(mesh_complex) .= groups
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
    ref_faces_inter = map(r_geom->lagrangian_reference_face(geometry(r_geom),order=my_order),ref_faces_geom)
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
    chain = Chain(;
        num_dims=Val(D),
        node_coordinates=node_coordinates_inter,
        face_nodes = face_nodes_inter,
        face_reference_id = face_ref_id_geom,
        reference_faces = ref_faces_inter,
    )
    mesh = mesh_from_chain(chain)
    mesh_complex, = complexify_mesh(mesh)
    pg = physical_groups(mesh_complex)
    pg .= physical_groups(mesh_geom)
    mesh_complex
end

function cartesian_mesh(domain,cells_per_dir;boundary=true,complexify=true,simplexify=false)
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
        mesh, = complexify_mesh(mesh)
    end
    mesh
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
    mesh_groups = push(groups,physical_groups(interior_mesh,D))
    fe_mesh(
     node_coords,
     mesh_face_nodes,
     mesh_face_reference_id,
     mesh_reference_faces;
     physical_groups=mesh_groups,
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
    ref_cell = lagrangian_fe(cell_geometry,order)
    reference_cells = [ref_cell]
    interior_cells = collect(Int32,1:length(cell_nodes))
    groups = Dict(["interior"=>interior_cells,"$D-face-1"=>interior_cells])
    chain = fe_chain(
        node_coords,
        cell_nodes,
        cell_reference_id,
        reference_cells;
        physical_groups=groups,
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
    chain = fe_chain(
        node_coords,
        cell_nodes,
        cell_reference_id,
        reference_cells;
        physical_groups=groups,
       )
    chain
end

function structured_simplex_mesh_with_boundary(domain,cells_per_dir)
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
      [ face_nodes(ref_simplex_mesh,d)[physical_groups(ref_simplex_mesh,d)["$d-face-$ldface"]]
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
    mesh_groups = push(groups,physical_groups(simplex_chain))
    fe_chain(
     node_coords,
     mesh_face_nodes,
     mesh_face_reference_id,
     mesh_reference_faces;
     physical_groups=mesh_groups,
    )
end


end # module

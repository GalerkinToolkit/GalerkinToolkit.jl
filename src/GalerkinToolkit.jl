module GalerkinToolkit

using WriteVTK
using FastGaussQuadrature
using StaticArrays
using LinearAlgebra
using ForwardDiff
using Gmsh
using PartitionedArrays

# TEMPORARY helper in the design process
# Using NamedTuples is starting to be too verbouse
# in stack traces.
struct AnonymousObject
    __data__::NamedTuple
    AnonymousObject(;kwargs...) = new((;kwargs...))
end

function Base.propertynames(x::AnonymousObject, private::Bool=false)
    propertynames(x.__data__,private)
end

function Base.getproperty(x::AnonymousObject,sym::Symbol)
    if sym === :__data__
        Base.getfield(x,sym)
    else
        Base.getproperty(x.__data__,sym)
    end
end

function set(a::AnonymousObject;kwargs...)
    AnonymousObject(;a.__data__...,kwargs...)
end

val_parameter(a) = a
val_parameter(::Val{a}) where a = a

num_dims(a) = val_parameter(a.num_dims)
node_coordinates(a) = a.node_coordinates
reference_faces(a) = a.reference_faces
face_nodes(a) = a.face_nodes
face_incidence(a) = a.face_incidence
face_reference_id(a) = a.face_reference_id
vtk_mesh_cell(a) = a.vtk_mesh_cell
physical_groups(a) = a.physical_groups
has_physical_groups(a) = hasproperty(a,:physical_groups)
geometry(a) = a.geometry
topology(a) = a.topology
boundary(a) = a.boundary
is_n_cube(a) = hasproperty(a,:is_n_cube) ? val_parameter(a.is_n_cube) : false
is_simplex(a) = hasproperty(a,:is_simplex) ? val_parameter(a.is_simplex) : false
is_axis_aligned(a) = a.is_axis_aligned
bounding_box(a) = a.bounding_box
coordinates(a) = a.coordinates
weights(a) = a.weights
interpolation(a) = a.interpolation
shape_functions(a) = a.shape_functions
gradient!(a) = a.gradient!
value!(a) = a.value!
order(a) = a.order
monomial_exponents(a) = a.monomial_exponents
lib_to_user_nodes(a) = a.lib_to_user_nodes
interior_nodes(a) = a.interior_nodes

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

function is_unit_n_cube(geom)
    !(is_axis_aligned(geom) && is_n_cube(geom)) && return false
    my_bounding_box = bounding_box(geom)
    all(i->i==0,first(my_bounding_box)) && all(i->i==1,last(my_bounding_box))
end

function is_unit_simplex(geom)
    !(is_axis_aligned(geom) && is_simplex(geom)) && return false
    my_bounding_box = bounding_box(geom)
    all(i->i==0,first(my_bounding_box)) && all(i->i==1,last(my_bounding_box))
end

function push(a::AbstractVector,x)
    b = copy(a)
    push!(b,x)
    b
end

function push(a::Tuple,x)
    (a...,x)
end

function mesh_from_reference_face(ref_face)
    boundary_mesh = boundary(interpolation(ref_face))
    D = num_dims(geometry(ref_face))
    nnodes = num_nodes(interpolation(ref_face))
    face_to_nodes = push(face_nodes(boundary_mesh),[collect(1:nnodes)])
    face_to_refid = push(face_reference_id(boundary_mesh),[1])
    refid_refface = push(reference_faces(boundary_mesh),[ref_face])
    node_to_coords = node_coordinates(interpolation(ref_face))
    AnonymousObject(;
      num_dims=Val(D),
      node_coordinates=node_to_coords,
      face_nodes=face_to_nodes,
      face_reference_id=face_to_refid,
      reference_faces=refid_refface)
end

evaluate(f,x) = f(x)

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

function classify_mesh_nodes!(node_to_tag,mesh,tag_to_name,dmax=num_dims(mesh))
    fill!(node_to_tag,zero(eltype(node_to_tag)))
    for d in dmax:-1:0
        face_to_nodes = face_nodes(mesh,d)
        face_groups = physical_groups(mesh,d)
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
    TwoPartPartition(free_nodes,dirichlet_nodes,permutation)
end

struct TwoPartPartition{A} <: AbstractVector{A}
    first::A
    last::A
    permutation::A
end

permutation(a::TwoPartPartition) = a.permutation
Base.size(a::TwoPartPartition) = (2,)
Base.IndexStyle(::Type{<:TwoPartPartition}) = IndexLinear()
function Base.getindex(a::TwoPartPartition,i::Int)
    @boundscheck @assert i in (1,2)
    if i == 1
        a.first
    else
        a.last
    end
end

function tensor_product_quadrature(
    degree_per_dir,
    limits_per_dir;
    weight_type=Float64,
    coordinate_type=SVector{length(limits_per_dir),weight_type},
    allocator=zeros)
    D = length(degree_per_dir)
    n_per_dir = map(d->ceil(Int,(d+1)/2),degree_per_dir)
    quad_per_dir = map(n_per_dir,limits_per_dir) do n,limits
        x,w = FastGaussQuadrature.gausslegendre(n)
        a,b = limits
        x .= (0.5*(b-a)) .*x .+ (0.5*(b+a))
        w .*= 0.5*(b-a)
        (;coordinates=x,weights=w)
    end
    coords_per_dir = map(coordinates,quad_per_dir)
    weights_per_dir = map(weights,quad_per_dir)
    m = prod(map(length,weights_per_dir))
    w = allocator(weight_type,m)
    x = allocator(coordinate_type,m)
    tensor_product!(identity,x,coords_per_dir)
    tensor_product!(prod,w,weights_per_dir)
    (;coordinates=x,weights=w)
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

function quadrature(geom,degree;kwargs...)
    if is_n_cube(geom)
        if is_axis_aligned(geom)
            my_bounding_box = bounding_box(geom)
            D = length(first(my_bounding_box))
            limits_per_dir = ntuple(i->(my_bounding_box[1][i],my_bounding_box[2][i]),Val(D))
            degree_per_dir = ntuple(i->degree,Val(D))
            tensor_product_quadrature(degree_per_dir,limits_per_dir;kwargs...)
        else
        error("Not implemented")
        end
    else
        error("Not implemented")
    end
end

function lagrange_monomial_exponents(f,order,D;kwargs...)
    function tensor_monomial_exponents(f,degree_per_dir; monomial_type=SVector{length(degree_per_dir),Int}, allocator=zeros)
        terms_per_dir = Tuple(map(d->d+1,degree_per_dir))
        D = length(terms_per_dir)
        cis = CartesianIndices(terms_per_dir)
        m = count(ci->f(Tuple(ci) .- 1,degree_per_dir),cis)
        li = 0
        result = zeros(monomial_type,m)
        for ci in cis
            t = Tuple(ci) .- 1
            if f(t,degree_per_dir)
                li += 1
                result[li] = t
            end
        end
        result
    end
    d = val_parameter(D)
    order_per_dir = ntuple(i->order,Val(d))
    tensor_monomial_exponents((e,os)->f(e,order),order_per_dir;kwargs...)
end

function lagrange_shape_functions(lagrange_monomial_exponents,node_coordinates)
    monomials = map(e->(x-> prod(x.^e)),lagrange_monomial_exponents)
    monomials_t = permutedims(monomials)
    A = evaluate.(monomials_t,node_coordinates)
    B = A\I
    function value!(r,x;scratch=similar(r))
        C = broadcast!(evaluate,scratch,monomials_t,x)
        mul!(r,C,B)
        r
    end
    function gradient!(r,x;scratch=similar(r,eltype(x)))
        C = broadcast!(ForwardDiff.gradient,scratch,monomials_t,x)
        mul!(r,C,B)
        r
    end
    (;value!,gradient!)
end

function lagrange_interpolation(
    geom,order;
    type=:default,
    lib_to_user_nodes = nothing,
    coordinate_type=SVector{num_dims(geom),Float64},
    kwargs...)
    f = if is_unit_n_cube(geom)
        if type === :default
            (e,o)->true
        elseif type === :Q
            (e,o)->true
        else
            error("Not implemented")
        end
    elseif is_unit_simplex(geom)
        if type === :default
            (e,o)->sum(e)<=o
        elseif type === :P
            (e,o)->sum(e)<=o
        else
            error("Not implemented")
        end
    else
        error("Not implemented")
    end
    d = num_dims(geom)
    monomial_exponents = lagrange_monomial_exponents(f,order,Val(d);kwargs...)
    T = eltype(coordinate_type)
    node_coordinates = map(monomial_exponents) do exponent
        c = map(exponent) do e
            if order != 0
                T(e/order)
            else
                T(e)
            end
        end
        coordinate_type(c)
    end
    nnodes = length(node_coordinates)
    lib_to_user_nodes = lib_to_user_nodes !== nothing ? lib_to_user_nodes : collect(1:nnodes)
    node_coordinates[lib_to_user_nodes] = node_coordinates
    shape_functions = lagrange_shape_functions(monomial_exponents,node_coordinates)
    boundary, interior_nodes = lagrange_reference_face_boundary(geom,node_coordinates,order)
    AnonymousObject(;shape_functions,node_coordinates,order,monomial_exponents,lib_to_user_nodes,boundary,interior_nodes)
end

function lagrange_reference_face_boundary(geom,node_coordinates_inter,order)
    D = num_dims(geom)
    if D == 0
        return nothing, collect(1:length(node_coordinates_inter))
    end
    mesh_geom = boundary(geom)
    node_coordinates_geom = node_coordinates(mesh_geom)
    face_nodes_inter = Vector{Vector{Vector{Int}}}(undef,D)
    node_coordinates_aux = map(xi->map(xii->round(Int,order*xii),xi),node_coordinates_inter)
    ref_faces_geom = reference_faces(mesh_geom)
    ref_faces = map(ref_faces_geom) do ref_faces_geom_d
        map(r->lagrange_reference_face(geometry(r),order),ref_faces_geom_d)
    end
    face_ref_id_geom = face_reference_id(mesh_geom)
    node_is_touched = fill(true,length(node_coordinates_inter))
    for d in 0:(D-1)
        s_ref = map(ref_faces_geom[d+1],ref_faces[d+1]) do r_geom,r
            m = num_nodes(interpolation(r))
            n = num_nodes(interpolation(r_geom))
            f = shape_functions(interpolation(r_geom))
            x = node_coordinates(interpolation(r))
            f.value!(zeros(m,n),x)
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
                map(xi->round(Int,order*xi),x)
            end
            my_nodes = indexin(x_mapped,node_coordinates_aux)
            node_is_touched[my_nodes] .= false
            face_nodes_inter_d[face] = my_nodes
        end
        face_nodes_inter[d+1] = face_nodes_inter_d
    end
    mesh_inter = AnonymousObject(;
        num_dims = Val(D-1),
        node_coordinates=node_coordinates_inter,
        face_nodes=face_nodes_inter,
        face_reference_id=face_ref_id_geom,
        reference_faces=ref_faces
    )
    interior_nodes = findall(node_is_touched)
    mesh_inter, interior_nodes
end

function lagrange_reference_face(geometry,args...;kwargs...)
    interpolation = lagrange_interpolation(geometry,args...;kwargs...)
    vtk_mesh_cell = vtk_mesh_cell_from_geometry(geometry,interpolation)
    AnonymousObject(;geometry,interpolation,vtk_mesh_cell)
end

function vtk_mesh_cell_from_geometry(geom,interpolation)
    d = num_dims(geom)
    nnodes = num_nodes(interpolation)
    lib_to_user = lib_to_user_nodes(interpolation)
    if d == 0 && nnodes == 1
        cell_type = WriteVTK.VTKCellTypes.VTK_VERTEX
        vtk_to_lex = [1]
    elseif d == 1 && (is_simplex(geom) || is_n_cube(geom)) && nnodes == 2
        cell_type = WriteVTK.VTKCellTypes.VTK_LINE
        vtk_to_lex = [1,2]
    elseif d == 1 && (is_simplex(geom) || is_n_cube(geom)) && nnodes == 3
        cell_type = WriteVTK.VTKCellTypes.VTK_QUADRATIC_EDGE
        vtk_to_lex = [1,3,2]
    elseif d == 2 && is_n_cube(geom) && nnodes == 4
        cell_type = WriteVTK.VTKCellTypes.VTK_QUAD
        vtk_to_lex = [1,2,4,3]
    elseif d == 2 && is_simplex(geom) && nnodes == 3
        cell_type = WriteVTK.VTKCellTypes.VTK_TRIANGLE
        vtk_to_lex = [1,2,3]
    elseif d == 2 && is_simplex(geom) && nnodes == 6
        cell_type = WriteVTK.VTKCellTypes.VTK_QUADRATIC_TRIANGLE
        vtk_to_lex = [1,3,6,2,4,5]
    elseif d == 3 && is_n_cube(geom) && nnodes == 8
        cell_type = WriteVTK.VTKCellTypes.VTK_HEXAHEDRON
        vtk_to_lex = [1,2,4,3,5,6,8,7]
    elseif d == 3 && is_simplex(geom) && nnodes == 4
        cell_type = WriteVTK.VTKCellTypes.VTK_TETRA
        vtk_to_lex = [1,2,3,4]
    else
        return nothing
    end
    nodes -> WriteVTK.MeshCell(cell_type,nodes[lib_to_user][vtk_to_lex])
end

function unit_n_cube(D;kwargs...)
    d = val_parameter(D)
    num_dims = Val(d)
    is_n_cube=true
    is_axis_aligned=true
    bounding_box=SVector{d,Float64}[ntuple(i->0,Val(d)),ntuple(i->1,Val(d))]
    is_simplex = d in (0,1)
    boundary = unit_n_cube_boundary(D;kwargs...)
    geometry = AnonymousObject(;num_dims,is_n_cube,is_simplex,is_axis_aligned,bounding_box,boundary)
    geometry
end

function unit_simplex(D;kwargs...)
    d = val_parameter(D)
    num_dims = Val(d)
    is_simplex=true
    is_axis_aligned=true
    is_n_cube = d in (0,1)
    bounding_box=SVector{d,Float64}[ntuple(i->0,Val(d)),ntuple(i->1,Val(d))]
    boundary = unit_simplex_boundary(D;kwargs...)
    geometry = AnonymousObject(;num_dims,is_n_cube,is_simplex,is_axis_aligned,bounding_box,boundary)
    geometry
end

function unit_n_cube_boundary(
    D;
    reference_face=nothing)

    # TODO order the boundary in a more standard way
    # allow the user define an alternative ordering

    d = val_parameter(D)
    if reference_face === nothing && d != 0
        reference_face = lagrange_reference_face(unit_n_cube(d-1),1)
    end
    my_boundary = if d == 0
        nothing
    elseif d == 1
        node_coordinates = SVector{1,Float64}[(0,),(1,)]
        face_nodes = [[[1],[2]]]
        face_reference_id = [[1,1]]
        vertex = reference_face
        my_reference_faces = ([vertex],)
        AnonymousObject(;num_dims=Val(0),node_coordinates,face_nodes,face_reference_id,reference_faces=my_reference_faces)
    elseif d == 2
        node_coordinates = SVector{2,Float64}[(0,0),(1,0),(0,1),(1,1)]
        face_nodes = [[[1],[2],[3],[4]],[[1,2],[3,4],[1,3],[2,4]]]
        face_reference_id = [[1,1,1,1],[1,1,1,1]]
        segment = reference_face
        vertex = first(reference_faces(boundary(geometry(segment)),0))
        my_reference_faces = ([vertex],[segment])
        AnonymousObject(;num_dims=Val(1),node_coordinates,face_nodes,face_reference_id,reference_faces=my_reference_faces)
    elseif d == 3
        node_coordinates = SVector{3,Float64}[(0,0,0),(1,0,0),(0,1,0),(1,1,0),(0,0,1),(1,0,1),(0,1,1),(1,1,1)]
        face_nodes = [
            [[1],[2],[3],[4],[5],[6],[7],[8]],
            [[1,2],[3,4],[1,3],[2,4],[5,6],[7,8],[5,7],[6,8],[1,5],[3,7],[2,6],[4,8]],
            [[1,2,3,4],[5,6,7,8],[1,2,5,6],[3,4,7,8],[1,3,5,7],[2,4,6,8]],
           ]
        face_reference_id = [ones(Int,8),ones(Int,12),ones(Int,6)]
        quad = reference_face
        segment = first(reference_faces(boundary(geometry(quad)),1))
        vertex = first(reference_faces(boundary(geometry(segment)),0))
        my_reference_faces = ([vertex],[segment],[quad])
        AnonymousObject(;num_dims=Val(2),node_coordinates,face_nodes,face_reference_id,reference_faces=my_reference_faces)
    else
        @error "Case not implemented"
    end
    if my_boundary !== nothing
        vtk_grid("debug_$d",vtk_args(my_boundary)...) |> vtk_save
        topology = topology_from_mesh(my_boundary)
        set(my_boundary;topology)
    else
        (;topology=nothing)
    end
end

function unit_simplex_boundary(
    D;
    reference_face=nothing)

    d = val_parameter(D)
    if reference_face === nothing && d != 0
        reference_face = lagrange_reference_face(unit_simplex(Val(d-1)),1)
    end
    my_boundary = if d == 0
        nothing
    elseif d == 1
        node_coordinates = SVector{1,Float64}[(0,),(1,)]
        face_nodes = [[[1],[2]]]
        face_reference_id = [[1,1]]
        vertex = reference_face
        my_reference_faces = ([vertex],)
        AnonymousObject(;num_dims=Val(0),node_coordinates,face_nodes,face_reference_id,reference_faces=my_reference_faces)
    elseif d == 2
        node_coordinates = SVector{2,Float64}[(0,0),(1,0),(0,1)]
        face_nodes = [[[1],[2],[3]],[[1,2],[1,3],[2,3]]]
        face_reference_id = [[1,1,1],[1,1,1]]
        segment = reference_face
        vertex = first(reference_faces(boundary(geometry(segment)),0))
        my_reference_faces = ([vertex],[segment])
        AnonymousObject(;num_dims=Val(1),node_coordinates,face_nodes,face_reference_id,reference_faces=my_reference_faces)
    elseif d == 3
        node_coordinates = SVector{3,Float64}[(0,0,0),(1,0,0),(0,1,0),(0,0,1)]
        face_nodes = [
            [[1],[2],[3],[4]],
            [[1,2],[1,3],[2,3],[1,4],[2,4],[3,4]],
            [[1,2,3],[2,3,4],[3,2,4],[3,1,4]],
           ]
        face_reference_id = [ones(Int,4),ones(Int,6),ones(Int,4)]
        tri = reference_face
        segment = first(reference_faces(boundary(geometry(tri)),1))
        vertex = first(reference_faces(boundary(geometry(segment)),0))
        my_reference_faces = ([vertex],[segment],[tri])
        AnonymousObject(;num_dims=Val(2),node_coordinates,face_nodes,face_reference_id,reference_faces=my_reference_faces)
    else
        @error "case not implemented"
    end
    if my_boundary !== nothing
        topology = topology_from_mesh(my_boundary)
        set(my_boundary;topology)
    else
        (;topology=nothing)
    end
end

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
            refdfaces = reference_faces(first(first(my_reference_faces)),d)
        end
        my_reference_faces = (refdfaces,my_reference_faces...)
    end

    ## Setup periodic nodes
    #node_to_main_node = fill(Int32(INVALID_ID),nnodes)
    #for (dim,tag) in entities
    #    tagMaster, nodeTags, nodeTagsMaster, = gmsh.model.mesh.getPeriodicNodes(dim,tag)
    #    for i in 1:length(nodeTags)
    #        node = nodeTags[i]
    #        main_node = nodeTagsMaster[i]
    #        node_to_main_node[node] = main_node
    #    end
    #end
    #my_free_and_periodic_nodes = partition_from_mask(i->i==INVALID_ID,node_to_main_node)
    #my_periodic_nodes = last(my_free_and_periodic_nodes)
    #my_periodic_to_master = node_to_main_node[my_periodic_nodes]
    #my_periodic_to_coeff = ones(length(my_periodic_to_master))

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

    mesh = AnonymousObject(;
            num_dims=Val(D),
            node_coordinates = my_node_to_coords,
            face_nodes = my_face_nodes,
            face_reference_id = my_face_reference_id,
            reference_faces = my_reference_faces,
            physical_groups = my_groups)

    if complexify
        mesh, _ = complexify_mesh(mesh)
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
    lagrange_reference_face(geom,order;lib_to_user_nodes=lib_to_gmsh)
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
        nrefid_to_ldface_to_lnodes = map(a->face_nodes(boundary(interpolation(a)),d),newreffaces[n+1])
        nrefid_to_ldface_to_drefrefid = map(a->face_reference_id(boundary(interpolation(a)),d),newreffaces[n+1])
        nrefid_to_drefrefid_to_ref_dface = map(a->reference_faces(boundary(interpolation(a)),d),newreffaces[n+1])
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
    new_mesh = AnonymousObject(;
        num_dims = Val(D),
        node_coordinates=node_to_coords,
        face_nodes=newface_nodes,
        face_reference_id=newface_refid,
        reference_faces=Tuple(newreffaces)
       )
    if has_physical_groups(mesh)
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
        new_mesh = set(new_mesh,physical_groups=new_physical_groups)
    end
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
            face_nodes(boundary(interpolation(a)),0)
        else
            [interior_nodes(interpolation(a))]
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
                face_nodes(boundary(interpolation(a)),0)
            else
                [interior_nodes(interpolation(a))]
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

function topology_from_mesh(mesh)
    # Assumes that the input is a cell complex
    T = JaggedArray{Int32,Int32}
    D = num_dims(mesh)
    my_face_incidence = Matrix{T}(undef,D+1,D+1)
    my_face_reference_id  = [ face_reference_id(mesh,d) for d in 0:D ]
    my_reference_faces = Tuple([ map(f->f|>geometry|>boundary|>topology,reference_faces(mesh,d)) for d in 0:D ])
    topo = AnonymousObject(;
        face_incidence=my_face_incidence,
        face_reference_id=my_face_reference_id,
        reference_faces=my_reference_faces
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
    refid_to_lvertex_to_lnodes = map(refface->face_nodes(boundary(interpolation(refface)),0),refid_refface)
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
    Drefid_to_ldface_to_lvertices = map(refface->face_incidence(refface,d,0),refid_refface)
    Dface_to_dfaces = barrier(
            Dface_to_vertices,
            vertex_to_Dfaces,
            dface_to_vertices,
            vertex_to_dfaces,
            Dface_to_refid,
            Drefid_to_ldface_to_lvertices)
    topo.face_incidence[D+1,d+1] = Dface_to_dfaces
end

function bounding_box_from_domain(domain)
    l = length(domain)
    D = div(l,2)
    pmin = SVector(ntuple(d->domain[2*(d-1)+1],Val(D)))
    pmax = SVector(ntuple(d->domain[2*d],Val(D)))
    (pmin,pmax)
end

function cartesian_mesh(domain,cells_per_dir)
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
    ref_cell = lagrange_reference_face(cell_geometry,order)
    reference_cells = [ref_cell]
    chain = (;
        num_dims=Val(D),
        node_coordinates=node_coords,
        face_nodes=cell_nodes,
        face_reference_id=cell_reference_id,
        reference_faces=reference_cells)
    mesh = mesh_from_chain(chain)
    mesh
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
    refid_to_refface = (ntuple(d->[],Val(D))...,reference_cells)
    AnonymousObject(;
      num_dims=Val(D),
      node_coordinates=node_coords,
      face_nodes=face_to_nodes,
      face_reference_id=face_to_refid,
      reference_faces=refid_to_refface,
      )
end




end # module

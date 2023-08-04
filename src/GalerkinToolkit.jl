module GalerkinToolkit

using WriteVTK
using FastGaussQuadrature
using StaticArrays
using LinearAlgebra
using ForwardDiff
using Gmsh
using PartitionedArrays
using Combinatorics

# Inspited by PropertyUtils.jl
# Use a dict directly instead? (provably not, we would not be type stable even if we wanted,really?)
struct Object
    item::Dict{Symbol,Any}
end
Object(;kwargs...) = Object(Dict(kwargs))
Base.propertynames(i::Object) = collect(keys(getfield(i, :item)))
Base.getproperty(i::Object, x::Symbol) = getindex(getfield(i, :item), x)
Base.setproperty!(i::Object, name::Symbol, x) = setindex!(getfield(i, :item), x, name)

# TODO setproperty instead?
function setproperties(a::Object;kwargs...)
    dict = getfield(a, :item)
    item = merge(dict,kwargs)
    Object(item)
end

function Base.show(io::IO,a::Object)
    print(io,"GalerkinToolkit.Object(")
    dict = getfield(a, :item)
    for (i,(k,v)) in enumerate(dict)
        if i != 1
            print(io,", ")
        end
        print(io,k," = ")
        if isa(v,Object)
            print(io,"GalerkinToolkit.Object(...)")
        elseif isa(v,AbstractArray) && length(v) > 10
            print(io,"[...]")
        else
            print(io,v)
        end
    end
    print(io,")")
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
periodic_nodes(a) = a.periodic_nodes
has_periodic_nodes(a) = hasproperty(a,:periodic_nodes)
geometry(a) = a.geometry
topology(a) = a.topology
boundary(a) = a.boundary
is_n_cube(a) = hasproperty(a,:is_n_cube) ? val_parameter(a.is_n_cube) : false
is_simplex(a) = hasproperty(a,:is_simplex) ? val_parameter(a.is_simplex) : false
is_axis_aligned(a) = a.is_axis_aligned
bounding_box(a) = a.bounding_box
vertex_permutations(a) = a.vertex_permutations
coordinates(a) = a.coordinates
weights(a) = a.weights
shape_functions(a) = a.shape_functions
tabulation_matrix!(a) = a.tabulation_matrix!
tabulation_matrix(a) = a.tabulation_matrix
order(a) = a.order
monomial_exponents(a) = a.monomial_exponents
lib_to_user_nodes(a) = a.lib_to_user_nodes
interior_nodes(a) = a.interior_nodes
face_permutation_ids(a) = a.face_permutation_ids

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

function mesh_from_reference_face(ref_face;physical_groups=Val(true))
    boundary_mesh = boundary(ref_face)
    D = num_dims(geometry(ref_face))
    nnodes = num_nodes(ref_face)
    face_to_nodes = push(face_nodes(boundary_mesh),[collect(1:nnodes)])
    face_to_refid = push(face_reference_id(boundary_mesh),[1])
    refid_refface = push(reference_faces(boundary_mesh),[ref_face])
    node_to_coords = node_coordinates(ref_face)
    mesh = Object(;
      num_dims=Val(D),
      node_coordinates=node_to_coords,
      face_nodes=face_to_nodes,
      face_reference_id=face_to_refid,
      reference_faces=refid_refface)
    if val_parameter(physical_groups)
        groups = [ Dict{String,Vector{Int32}}() for d in 0:D]
        for d in 0:D
            for face in 1:num_faces(mesh,d)
                groups[d+1]["$d-face-$face"] = [face]
            end
        end
        groups[end-1]["boundary"] = 1:num_faces(mesh,D-1)
        groups[end]["interior"] = [1]
        mesh = setproperties(mesh,physical_groups=groups)
    end
    mesh
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

function monomial_exponents_from_filter(
    f,order_per_dir; monomial_type=SVector{length(order_per_dir),Int}, allocator=zeros)
    terms_per_dir = Tuple(map(d->d+1,order_per_dir))
    D = length(terms_per_dir)
    cis = CartesianIndices(terms_per_dir)
    m = count(ci->f(Tuple(ci) .- 1,order_per_dir),cis)
    li = 0
    result = zeros(monomial_type,m)
    for ci in cis
        t = Tuple(ci) .- 1
        if f(t,order_per_dir)
            li += 1
            result[li] = t
        end
    end
    result
end

function monomial_exponents_from_space(space,order_per_dir;kwargs...)
    filter = if space == :Q
        (e,o)->true
    elseif space == :P
        (e,o)->sum(e)<=maximum(o)
    else
        error("Case not implemented (yet)")
    end
    monomial_exponents_from_filter(filter,order_per_dir;kwargs...)
end

function shape_functions_and_dofs_from_bases(primal,dual)
    primal_t = permutedims(primal)
    function dofs_from_dual_basis(dual)
        function tabulation_matrix(primal)
            A = value.(dual,primal_t)
        end
        (;tabulation_matrix)
    end
    dofs = dofs_from_dual_basis(dual)
    A = dofs.tabulation_matrix(primal)
    B = A\I
    function tabulation_matrix!(f,r,x;scratch=nothing)
        C = if scratch !== nothing
            broadcast!(f,scratch,primal_t,x)
        else
            broadcast(f,primal_t,x)
        end
        mul!(r,C,B)
        r
    end
    function tabulation_matrix(f,x)
        C = broadcast(f,primal_t,x)
        C*B
    end
    # We want this to be type stable for
    # performance reasons
    shape_functions = (;tabulation_matrix!,tabulation_matrix)
    shape_functions, dofs
end


value(f,x) = f(x)

function default_polynomial_space(geo)
    if is_n_cube(geo) 
        :Q
    elseif is_simplex(geo)
        :P
    else
        error("Case not implemented")
    end
end

function reference_face_from_geometry(
    geometry;
    order=1,
    order_per_dir = ntuple(i->order,Val(num_dims(geometry))),
    space=:default,
    interior_node_permutations=Val(true),
    lib_to_user_nodes = nothing,
    coordinate_type=SVector{num_dims(geometry),Float64},
    kwargs...)

    D = num_dims(geometry)
    if space === :default
        space = default_polynomial_space(geometry)
    end
    monomial_exponents =
        monomial_exponents_from_space(space,order_per_dir,kwargs...)
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
    lib_to_user_nodes = if lib_to_user_nodes !== nothing
        lib_to_user_nodes 
    else
        collect(1:nnodes)
    end
    node_coordinates[lib_to_user_nodes] = node_coordinates
    monomials = map(e->(x-> prod(x.^e)),monomial_exponents)
    dofs = map(x->(f->f(x)),node_coordinates)
    shape_functions, dofs = shape_functions_and_dofs_from_bases(monomials,dofs)
    refface = Object(;
        geometry,
        shape_functions,
        dofs,
        node_coordinates,
        order= D==0 ? 1 : maximum(order_per_dir),
        order_per_dir,
        monomial_exponents,
        lib_to_user_nodes)
    boundary, interior_nodes = reference_face_boundary_from_reference_face(refface)
    vtk_mesh_cell = vtk_mesh_cell_from_reference_face(refface)
    refface1 = setproperties(refface;boundary,interior_nodes,vtk_mesh_cell)
    if val_parameter(interior_node_permutations)
        interior_node_permutations = interior_node_permutations_from_reference_face(refface1)
    else
        interior_node_permutations = nothing
    end
    setproperties(refface1;interior_node_permutations)
end

function reference_face_boundary_from_reference_face(refface)
    geom = geometry(refface)
    node_coordinates_inter = node_coordinates(refface)
    order_inter = order(refface)
    D = num_dims(geom)
    if D == 0
        return nothing, collect(1:length(node_coordinates_inter))
    end
    mesh_geom = boundary(geom)
    node_coordinates_geom = node_coordinates(mesh_geom)
    face_nodes_inter = Vector{Vector{Vector{Int}}}(undef,D)
    node_coordinates_aux = map(xi->map(xii->round(Int,order_inter*xii),xi),node_coordinates_inter)
    ref_faces_geom = reference_faces(mesh_geom)
    ref_faces = map(ref_faces_geom) do ref_faces_geom_d
        map(r->reference_face_from_geometry(geometry(r);order=order_inter),ref_faces_geom_d)
    end
    face_ref_id_geom = face_reference_id(mesh_geom)
    node_is_touched = fill(true,length(node_coordinates_inter))
    for d in 0:(D-1)
        s_ref = map(ref_faces_geom[d+1],ref_faces[d+1]) do r_geom,r
            m = num_nodes(r)
            n = num_nodes(r_geom)
            f = shape_functions(r_geom)
            x = node_coordinates(r)
            tabulation_matrix(f)(value,x)
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
            node_is_touched[my_nodes] .= false
            face_nodes_inter_d[face] = my_nodes
        end
        face_nodes_inter[d+1] = face_nodes_inter_d
    end
    mesh_inter = Object(;
        num_dims = Val(D-1),
        node_coordinates=node_coordinates_inter,
        face_nodes=face_nodes_inter,
        face_reference_id=face_ref_id_geom,
        reference_faces=ref_faces
    )
    interior_nodes = findall(node_is_touched)
    mesh_inter, interior_nodes
end

function vtk_mesh_cell_from_reference_face(ref_face)
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

function unit_n_cube(D;kwargs...)
    d = val_parameter(D)
    num_dims = Val(d)
    is_n_cube=true
    is_axis_aligned=true
    bounding_box=SVector{d,Float64}[ntuple(i->0,Val(d)),ntuple(i->1,Val(d))]
    is_simplex = d in (0,1)
    boundary = unit_n_cube_boundary(D;kwargs...)
    geometry = Object(;num_dims,is_n_cube,is_simplex,is_axis_aligned,bounding_box,boundary)
    vertex_permutations = vertex_permutations_from_geometry(geometry)
    setproperties(geometry;vertex_permutations)
end

function unit_simplex(D;kwargs...)
    d = val_parameter(D)
    num_dims = Val(d)
    is_simplex=true
    is_axis_aligned=true
    is_n_cube = d in (0,1)
    bounding_box=SVector{d,Float64}[ntuple(i->0,Val(d)),ntuple(i->1,Val(d))]
    boundary = unit_simplex_boundary(D;kwargs...)
    geometry = Object(;num_dims,is_n_cube,is_simplex,is_axis_aligned,bounding_box,boundary)
    vertex_permutations = vertex_permutations_from_geometry(geometry)
    setproperties(geometry;vertex_permutations)
end

function unit_n_cube_boundary(
    D;
    reference_face=nothing)

    # TODO use the same convention than in Gridap
    # allow the user define an alternative ordering

    d = val_parameter(D)
    if d == 0
        return nothing
    end
    if reference_face === nothing
        reference_face = reference_face_from_geometry(unit_n_cube(d-1))
    end
    my_boundary = if d == 1
        node_coordinates = SVector{1,Float64}[(0,),(1,)]
        face_nodes = [[[1],[2]]]
        face_reference_id = [[1,1]]
        vertex = reference_face
        my_reference_faces = ([vertex],)
        Object(;num_dims=Val(0),node_coordinates,face_nodes,face_reference_id,reference_faces=my_reference_faces)
    elseif d == 2
        node_coordinates = SVector{2,Float64}[(0,0),(1,0),(0,1),(1,1)]
        face_nodes = [[[1],[2],[3],[4]],[[1,2],[3,4],[1,3],[2,4]]]
        face_reference_id = [[1,1,1,1],[1,1,1,1]]
        segment = reference_face
        vertex = first(reference_faces(boundary(geometry(segment)),0))
        my_reference_faces = ([vertex],[segment])
        Object(;num_dims=Val(1),node_coordinates,face_nodes,face_reference_id,reference_faces=my_reference_faces)
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
        Object(;num_dims=Val(2),node_coordinates,face_nodes,face_reference_id,reference_faces=my_reference_faces)
    else
        @error "Case not implemented"
    end
    topology = topology_from_mesh(my_boundary)
    setproperties(my_boundary;topology)
end

function unit_simplex_boundary(
    # TODO use the same convention than in Gridap
    # allow the user to customize the boundary object ids
    D;
    reference_face=nothing)

    d = val_parameter(D)
    if d == 0
        return nothing
    end
    if reference_face === nothing
        reference_face = reference_face_from_geometry(unit_simplex(Val(d-1)))
    end
    my_boundary = if d == 1
        node_coordinates = SVector{1,Float64}[(0,),(1,)]
        face_nodes = [[[1],[2]]]
        face_reference_id = [[1,1]]
        vertex = reference_face
        my_reference_faces = ([vertex],)
        Object(;num_dims=Val(0),node_coordinates,face_nodes,face_reference_id,reference_faces=my_reference_faces)
    elseif d == 2
        node_coordinates = SVector{2,Float64}[(0,0),(1,0),(0,1)]
        face_nodes = [[[1],[2],[3]],[[1,2],[1,3],[2,3]]]
        face_reference_id = [[1,1,1],[1,1,1]]
        segment = reference_face
        vertex = first(reference_faces(boundary(geometry(segment)),0))
        my_reference_faces = ([vertex],[segment])
        Object(;num_dims=Val(1),node_coordinates,face_nodes,face_reference_id,reference_faces=my_reference_faces)
    elseif d == 3
        node_coordinates = SVector{3,Float64}[(0,0,0),(1,0,0),(0,1,0),(0,0,1)]
        face_nodes = [
            [[1],[2],[3],[4]],
            [[1,2],[1,3],[2,3],[1,4],[2,4],[3,4]],
            [[1,2,3],[1,2,4],[1,3,4],[2,3,4]],
           ]
        face_reference_id = [ones(Int,4),ones(Int,6),ones(Int,4)]
        tri = reference_face
        segment = first(reference_faces(boundary(geometry(tri)),1))
        vertex = first(reference_faces(boundary(geometry(segment)),0))
        my_reference_faces = ([vertex],[segment],[tri])
        Object(;num_dims=Val(2),node_coordinates,face_nodes,face_reference_id,reference_faces=my_reference_faces)
    else
        error("case not implemented")
    end
    topology = topology_from_mesh(my_boundary)
    setproperties(my_boundary;topology)
end

function vertex_permutations_from_geometry(geo)
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
    ref_face = reference_face_from_geometry(geo,interior_node_permutations=Val(false))
    fun_mesh = boundary(ref_face)
    geo_node_coords = node_coordinates(geo_mesh)
    fun_node_coords = node_coordinates(fun_mesh)
    vertex_coords = geo_node_coords[vertex_to_geo_node]
    shape_funs = shape_functions(ref_face)
    degree = 1
    quad = quadrature(geo,degree)
    q = coordinates(quad)
    Tx = eltype(vertex_coords)
    TJ = typeof(zero(Tx)*zero(Tx)')
    A = tabulation_matrix(shape_funs)(ForwardDiff.gradient,q)
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

function interior_node_permutations_from_reference_face(refface)
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
    ref_face = reference_face_from_geometry(geo)
    fun_mesh = boundary(ref_face)
    geo_node_coords = node_coordinates(geo_mesh)
    fun_node_coords = node_coordinates(fun_mesh)
    vertex_coords = geo_node_coords[vertex_to_geo_node]
    shape_funs = shape_functions(ref_face)
    q = ho_nodes_coordinates[interior_ho_nodes]
    Tx = eltype(vertex_coords)
    A = zeros(Float64,length(q),length(fun_node_coords))
    A = tabulation_matrix(shape_funs)(value,q)
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

function simplexify_reference_geometry(geo)
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
    refface = reference_face_from_geometry(geo)
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
    ref_cell = reference_face_from_geometry(simplex)
    node_coords = node_coordinates(boundary(geo))
    cell_nodes = simplex_node_ids_n_cube(geo)
    ncells = length(cell_nodes)
    cell_reference_id = fill(Int8(1),ncells)
    reference_cells = [ref_cell]
    chain = (;
        num_dims=Val(D),
        node_coordinates=node_coords,
        face_nodes=cell_nodes,
        face_reference_id=cell_reference_id,
        reference_faces=reference_cells,
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
    mesh_complex = setproperties(mesh_complex,physical_groups=groups)
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
    mesh_geom  = simplexify_reference_geometry(geometry(ref_face))
    D = num_dims(mesh_geom)
    node_coordinates_geom = node_coordinates(mesh_geom)
    ref_faces_geom = reference_faces(mesh_geom,D)
    face_nodes_geom = face_nodes(mesh_geom,D)
    face_ref_id_geom = face_reference_id(mesh_geom,D)
    nfaces = length(face_ref_id_geom)
    # We need the same order in all directions
    # for this to make sense
    my_order = order(ref_face)
    ref_faces_inter = map(r_geom->reference_face_from_geometry(geometry(r_geom),order=my_order),ref_faces_geom)
    s_ref = map(ref_faces_geom,ref_faces_inter) do r_geom,r
        m = num_nodes(r)
        n = num_nodes(r_geom)
        f = shape_functions(r_geom)
        x = node_coordinates(r)
        tabulation_matrix(f)(value,x)
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
    chain = Object(;
        num_dims=Val(D),
        node_coordinates=node_coordinates_inter,
        face_nodes = face_nodes_inter,
        face_reference_id = face_ref_id_geom,
        reference_faces = ref_faces_inter,
    )
    mesh = mesh_from_chain(chain)
    mesh_complex, = complexify_mesh(mesh)
    if has_physical_groups(mesh_geom)
        mesh_complex = setproperties(mesh_complex,physical_groups=physical_groups(mesh_geom))
    end
    mesh_complex
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

function mesh_from_gmsh(file;complexify=true,topology=true,renumber=true,kwargs...)
    @assert ispath(file) "File not found: $(file)"
    with_gmsh(;kwargs...) do
        gmsh.open(file)
        renumber && gmsh.model.mesh.renumberNodes()
        renumber && gmsh.model.mesh.renumberElements()
        mesh_from_gmsh_module(;complexify,topology)
    end
end

function mesh_from_gmsh_module(;complexify=true,topology=true)
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
    mesh = Object(;
            num_dims=Val(D),
            node_coordinates = my_node_to_coords,
            face_nodes = my_face_nodes,
            face_reference_id = my_face_reference_id,
            reference_faces = my_reference_faces,
            physical_groups = my_groups,
            periodic_nodes,
           )

    if complexify
        mesh, _ = complexify_mesh(mesh)
        if topology
            topo = topology_from_mesh(mesh)
            mesh = setproperties(mesh,topology=topo)
        end
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
    reference_face_from_geometry(geom;order,lib_to_user_nodes=lib_to_gmsh)
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
    new_mesh = Object(;
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
        new_mesh = setproperties(new_mesh,physical_groups=new_physical_groups)
    end
    if has_periodic_nodes(mesh)
        new_mesh = setproperties(new_mesh;periodic_nodes=periodic_nodes(mesh))
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

function reference_topology_from_reference_face(refface)
    geom = geometry(refface)
    D = num_dims(geom)
    if D != 0
        myboundary = geom |> boundary |> topology
    else
        myboundary = nothing
    end
    myperms = vertex_permutations(geom)
    Object(;num_dims=Val(D),boundary=myboundary,vertex_permutations=myperms)
end

function topology_from_mesh(mesh)
    # Assumes that the input is a cell complex
    T = JaggedArray{Int32,Int32}
    D = num_dims(mesh)
    my_face_incidence = Matrix{T}(undef,D+1,D+1)
    my_face_reference_id  = [ face_reference_id(mesh,d) for d in 0:D ]
    my_reference_faces = Tuple([ map(reference_topology_from_reference_face,reference_faces(mesh,d)) for d in 0:D ])
    my_face_permutation_ids = Matrix{T}(undef,D+1,D+1)
    topo = Object(;
        face_incidence=my_face_incidence,
        face_reference_id=my_face_reference_id,
        reference_faces=my_reference_faces,
        face_permutation_ids=my_face_permutation_ids,
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
        return cell_to_lface_to_pindex
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

function cartesian_mesh(domain,cells_per_dir;boundary=true,complexify=true,simplexify=false,topology=true)
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
        if topology
            topo = topology_from_mesh(mesh)
            mesh = setproperties(mesh,topology=topo)
        end
    end
    mesh
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
    Object(;
     num_dims=Val(D),
     node_coordinates=node_coords,
     face_nodes=mesh_face_nodes,
     face_reference_id=mesh_face_reference_id,
     reference_faces=mesh_reference_faces,
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
    ref_cell = reference_face_from_geometry(cell_geometry)
    reference_cells = [ref_cell]
    interior_cells = collect(Int32,1:length(cell_nodes))
    groups = Dict(["interior"=>interior_cells,"$D-face-1"=>interior_cells])
    chain = Object(;
        num_dims=Val(D),
        node_coordinates=node_coords,
        face_nodes=cell_nodes,
        face_reference_id=cell_reference_id,
        reference_faces=reference_cells,
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
    ref_simplex_mesh = simplexify_reference_geometry(cell_geometry)
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
    chain = Object(;
        num_dims=Val(D),
        node_coordinates=node_coords,
        face_nodes=cell_nodes,
        face_reference_id=cell_reference_id,
        reference_faces=reference_cells,
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
    ref_simplex_mesh = simplexify_reference_geometry(cell_geometry)
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
    Object(;
     num_dims=Val(D),
     node_coordinates=node_coords,
     face_nodes=mesh_face_nodes,
     face_reference_id=mesh_face_reference_id,
     reference_faces=mesh_reference_faces,
     physical_groups=mesh_groups,
    )
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
    mesh = Object(;
      num_dims=Val(D),
      node_coordinates=node_coords,
      face_nodes=face_to_nodes,
      face_reference_id=face_to_refid,
      reference_faces=refid_to_refface,
      )
    if has_physical_groups(chain)
        cell_groups = physical_groups(chain)
        groups = [ typeof(cell_groups)() for d in 0:D]
        groups[end] = cell_groups
        mesh = setproperties(mesh,physical_groups=groups)
    end
    mesh
end

function visualization_mesh_from_mesh(mesh,dim=num_dims(mesh);order=nothing,resolution=nothing)
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
        vchain = Object(;
                        num_dims=Val(dim),
                        node_coordinates=vnode_to_coords,
                        face_nodes=vcell_to_vnodes,
                        face_reference_id=vcell_to_vrefid,
                        reference_faces=vrefid_to_vrefface)
        vmesh = mesh_from_chain(vchain)
        vglue = (;parent_face=vcell_to_cell,
                 reference_coordinates=refid_to_snode_to_coords,
                 face_fine_nodes = cell_to_vnodes,
                 num_dims=Val(dim))
        vmesh, vglue
    end # barrier
    refid_to_refface = reference_faces(mesh,dim)
    refid_to_refmesh = map(refid_to_refface) do ref_face
        if order === nothing && resolution === nothing
            # Use the given cells as visualization cells
            mesh_from_reference_face(ref_face)
        elseif order !== nothing && resolution === nothing
            # Use cells of given order as visualization cells
            geo = geometry(ref_face)
            ref_face_ho = reference_face_from_geometry(geo;order)
            mesh_from_reference_face(ref_face_ho)
        elseif order === nothing && resolution !== nothing
            # Use linear sub-cells with $resolution per direction per direction
            geom = geometry(ref_face)
            refine_reference_geometry(geom,resolution)
        else
            error("order and resolution kw-arguments can not be given at the same time")
        end
    end
    refid_to_tabulation = map(refid_to_refface,refid_to_refmesh) do refface,refmesh
        shape_funs = shape_functions(refface)
        x = node_coordinates(refmesh)
        tabulation_matrix(shape_funs)(value,x)
    end
    refid_to_scell_to_snodes = map(refmesh->face_nodes(refmesh,dim),refid_to_refmesh)
    refid_to_scell_to_srefid = map(refmesh->face_reference_id(refmesh,dim),refid_to_refmesh)
    refid_to_srefid_to_oid = map(refmesh->map(objectid,reference_faces(refmesh,dim)),refid_to_refmesh)
    refid_to_srefid_to_vrefface = map(refmesh->reference_faces(refmesh,dim),refid_to_refmesh)
    refid_to_snode_to_coords = map(node_coordinates,refid_to_refmesh)
    node_to_coords = node_coordinates(mesh)
    cell_to_nodes = face_nodes(mesh,dim)
    cell_to_refid = face_reference_id(mesh,dim)
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
        refface = reference_face_from_geometry(geo)
        chain = Object(;
                       num_dims=Val(2),
                       node_coordinates = X,
                       face_nodes = T,
                       face_reference_id = fill(1,length(T)),
                       reference_faces = [refface]
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
        refface = reference_face_from_geometry(geo)
        chain = Object(;
                       num_dims=Val(3),
                       node_coordinates = X,
                       face_nodes = T,
                       face_reference_id = fill(1,length(T)),
                       reference_faces = [refface],
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


# This should replace reference_face_from_geometry
# TODO
# two different versions
# shape = :scalar, (), (1,) etc
# :scalar should be different than ()
function lagrangian_reference_element(
    geom;
    space=:default,
    order=1,
    shape=Val(()),
    init_tensor=SArray{Tuple{val_parameter(shape)...}},
    major=:component,
    inner=(a,b)->sum(map(*,a,b)),
    order_per_dir=ntuple(i->order,Val(num_dims(geom))),
    lib_to_user_nodes = nothing,
    init_coordinate=SVector,
    interior_node_permutations=Val(true),
    )

    if space === :default
        space = default_polynomial_space(geom)
    end
    monomial_exponents = monomial_exponents_from_space(space,order_per_dir)
    monomials = map(e->(x-> prod(x.^e)),monomial_exponents)
    myshape = val_parameter(shape)
    cis = CartesianIndices(myshape)
    cis_flat = cis[:]
    l = prod(myshape)
    if myshape == ()
        tensor_basis = [1]
    else
        tensor_basis = map(cis_flat) do ci
            init_tensor(ntuple(j->cis[j]==ci ? 1 : 0 ,Val(l)))
        end
    end
    if myshape == ()
        primal = monomials
    else
        primal_nested = map(monomials) do monomial
            map(tensor_basis)  do e
                x -> monomial(x)*e
            end
        end
        primal = reduce(vcat,primal_nested)
    end
    T = Float64
    node_coordinates = map(monomial_exponents) do exponent
        c = map(exponent) do e
            if order != 0
                T(e/order)
            else
                T(e)
            end
        end
        init_coordinate(c)
    end
    nnodes = length(node_coordinates)
    lib_to_user_nodes = if lib_to_user_nodes !== nothing
        lib_to_user_nodes 
    else
        collect(1:nnodes)
    end
    node_coordinates[lib_to_user_nodes] = node_coordinates
    if myshape == ()
       dual = map(x->(f->f(x)),node_coordinates)
       dof_to_node = 1:nnodes
       dof_to_index = fill(CartesianIndex(()),nnodes)
    else
        if major === :component
            dual_nested = map(node_coordinates) do x
                map(tensor_basis) do e
                    f->inner(e,f(x))
                end
            end
            dof_to_node_nested = map(1:nnodes) do inode
                fill(inode,length(cis))
            end
            dof_to_index_nested = map(1:nnodes) do inode
                cis[:]
            end
        elseif major === :node
            dual_nested = map(tensor_basis) do e
                map(node_coordinates) do x
                    f->inner(e,f(x))
                end
            end
            dof_to_node_nested = map(cis_flat) do ci
                collect(1:nnodes)
            end
            dof_to_index_nested = map(cis_flat) do ci
                fill(ci,nnodes)
            end
        end
        dual = reduce(vcat,dual_nested)
        dof_to_node = reduce(vcat,dof_to_node_nested)
        dof_to_index = reduce(vcat,dof_to_index_nested)
    end
    node_to_dofs = zeros(eltype(tensor_basis),nnodes)
    ndofs = length(dof_to_node)
    lis = LinearIndices(cis)
    for dof in 1:ndofs
        node = dof_to_node[dof]
        ci = dof_to_index[dof]
        node_to_dofs[node] += dof*tensor_basis[lis[ci]]
    end
    shape_functions, dofs = shape_functions_and_dofs_from_bases(primal,dual)
    D = num_dims(geom)
    refface = Object(;
        geometry=geom,
        num_dofs = ndofs,
        shape_functions,
        dofs,
        node_coordinates,
        order= D==0 ? 1 : maximum(order_per_dir),
        order_per_dir,
        monomial_exponents,
        lib_to_user_nodes,
        dof_to_node,
        dof_to_index,
        node_to_dofs,
       )
    boundary, interior_nodes_refface = reference_face_boundary_from_reference_face(refface)
    if length(interior_nodes_refface) !=0 && myshape !=()
        own_dofs_nested = map(interior_nodes_refface) do node
            collect(node_to_dofs[node][:])
        end
        own_dofs = reduce(vcat,own_dofs_nested)
    else
        own_dofs = interior_nodes_refface
    end
    face_own_dofs = Vector{Vector{Vector{Int}}}(undef,D+1)
    face_own_dofs[end] = [own_dofs]
    for d in 0:(D-1)
        face_to_refid = face_reference_id(boundary,d)
        face_to_nodes = face_nodes(boundary,d)
        refid_to_refface = reference_faces(boundary,d)
        refid_to_interior_nodes = map(interior_nodes,refid_to_refface)
        ndfaces = num_faces(boundary,d)
        face_own_dofs[d+1] = Vector{Vector{Int}}(undef,ndfaces)
        for face in 1:ndfaces
            refid = face_to_refid[face]
            my_interior_nodes = refid_to_interior_nodes[refid]
            mynodes = face_to_nodes[face][my_interior_nodes]
            if length(mynodes) !=0 && myshape != ()
                own_dofs_nested = map(mynodes) do node
                    collect(node_to_dofs[node][:])
                end
                face_own_dofs[d+1][face] = reduce(vcat,own_dofs_nested)
            else
                face_own_dofs[d+1][face] = mynodes
            end
        end
    end
    vtk_mesh_cell = vtk_mesh_cell_from_reference_face(refface)
    refface1 = setproperties(refface;boundary,interior_nodes=interior_nodes_refface,vtk_mesh_cell)
    if val_parameter(interior_node_permutations)
        interior_node_permutations = interior_node_permutations_from_reference_face(refface1)
    else
        interior_node_permutations = nothing
    end
    #TODO interior dofs on the boundary and permutations
    setproperties(refface1;interior_node_permutations,face_own_dofs)
end


end # module

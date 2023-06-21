module GalerkinToolkit

using WriteVTK
using FastGaussQuadrature
using StaticArrays
using LinearAlgebra
using ForwardDiff

export vtk_args

val_parameter(a) = a
val_parameter(::Val{a}) where a = a

num_dims(a) = val_parameter(a.num_dims)
node_coordinates(a) = a.node_coordinates
reference_faces(a) = a.reference_faces
face_nodes(a) = a.face_nodes
face_reference_id(a) = a.face_reference_id
vtk_mesh_cell(a) = a.vtk_mesh_cell
physical_groups(a) = a.physical_groups
geometry(a) = a.geometry
is_n_cube(a) = a.is_n_cube
is_axis_aligned(a) = a.is_axis_aligned
bounding_box(a) = a.bounding_box
coordinates(a) = a.coordinates
weights(a) = a.weights
interpolation(a) = a.interpolation
shape_functions(a) = a.shape_functions
gradient!(a) = a.gradient!
value!(a) = a.value!

reference_faces(a,d) = reference_faces(a)[val_parameter(d)+1]
face_nodes(a,d) = face_nodes(a)[val_parameter(d)+1]
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
    allcells = [vtk_cells(mesh,d) for d in 0:D]
    cells = reduce(vcat,allcells)
    points, cells
end

function vtk_physical_groups!(vtk,geo,d;physical_groups=physical_groups(geo,d))
    ndfaces = num_faces(geo,d)
    for group in physical_groups
        name,faces = group
        face_mask = zeros(Int,ndfaces)
        face_mask[faces] .= 1
        vtk[name,WriteVTK.VTKCellData()] = face_mask
    end
    vtk
end

function vtk_physical_groups!(vtk,geo;physical_groups=physical_groups(geo))
    nfaces = sum(num_faces(geo))
    offsets = face_offsets(geo)
    D = num_dims(geo)
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

function quadrature(geom,degree_per_dir;kwargs...)
    if is_n_cube(geom)
        if is_axis_aligned(geom)
            my_bounding_box = bounding_box(geom)
            D = length(first(my_bounding_box))
            limits_per_dir = ntuple(i->(my_bounding_box[1][i],my_bounding_box[2][i]),Val(D))
            tensor_product_quadrature(degree_per_dir,limits_per_dir;kwargs...)
        else
        error("Not implemented")
        end
    else
        error("Not implemented")
    end
end

function monomial_exponents(f,degree_per_dir; monomial_type=SVector{length(degree_per_dir),Int}, allocator=zeros)
  terms_per_dir = Tuple(map(d->d+1,degree_per_dir))
  D = length(terms_per_dir)
  cis = CartesianIndices(terms_per_dir)
  lis = LinearIndices(cis)
  m = count(ci->f(Tuple(ci) .- 1),cis)
  result = zeros(monomial_type,m)
  for ci in cis
      t = Tuple(ci) .- 1
      if f(t)
          li = lis[ci]
          result[li] = t
      end
  end
  result
end

function lagrange_shape_functions(monomial_exponents,node_coordinates)
    monomials = map(e->(x-> prod(x.^e)),monomial_exponents)
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
    geom,degree_per_dir;
    type=:default,
    coordinate_type=SVector{length(degree_per_dir),Float64},
    kwargs...)
    f = if is_unit_n_cube(geom)
        if type === :default
            e->true
        elseif type === :Q
            e->true
        else
            error("Not implemented")
        end
    elseif is_unit_simplex(geom)
        if type === :default
            e->true
        elseif type === :Q
            e->true
        else
            error("Not implemented")
        end
    else
        error("Not implemented")
    end
    exponents = monomial_exponents(f,degree_per_dir;kwargs...)
    node_coordinates = map(coordinate_type,exponents)
    shape_functions = lagrange_shape_functions(exponents,node_coordinates)
    (;shape_functions,node_coordinates)
end

end # module

module MeshInterfaceTest

using FastGaussQuadrature
using StaticArrays
using LinearAlgebra
using ForwardDiff

get_coordinates(a) = a.coordinates
get_weights(a) = a.weights

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
    coords_per_dir = map(get_coordinates,quad_per_dir)
    weights_per_dir = map(get_weights,quad_per_dir)
    m = prod(map(length,weights_per_dir))
    weights = allocator(weight_type,m)
    coordinates = allocator(coordinate_type,m)
    tensor_product!(identity,coordinates,coords_per_dir)
    tensor_product!(prod,weights,weights_per_dir)
    (;coordinates,weights)
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


get_val_parameter(a) = a
get_val_parameter(::Val{a}) where a = a

has_is_n_cube(a) = true
get_is_n_cube(a) = hasproperty(a,:is_n_cube) ? get_val_parameter(a.is_n_cube) : false

has_is_simplex(a) = true
get_is_simplex(a) = hasproperty(a,:is_simplex) ? get_val_parameter(a.is_simplex) : false

has_is_axis_aligned(a) = true
get_is_axis_aligned(a) = hasproperty(a,:is_axis_aligned) ? get_val_parameter(a.is_axis_aligned) : false

has_bounding_box(a) = hasproperty(a,:bounding_box)
get_bounding_box(a) = a.bounding_box

has_geometry(a) = hasproperty(a,:geometry)
get_geometry(a) = a.geometry

function has_num_dims(a)
    hasproperty(a,:num_dims) && return true
    has_geometry(a) && return has_num_dims(get_geometry(a))
    false
end

function get_num_dims(a)
    hasproperty(a,:num_dims) && return get_val_parameter(a.num_dims)
    (has_geometry(a) && has_num_dims(get_geometry(a))) && return get_num_dims(get_geometry(a))
    error("num_dims not found in data")
end

function has_num_ambient_dims(a)
    hasproperty(a,:num_ambient_dims) && return true
    (has_geometry(a) && has_num_ambient_dims(get_geometry(a))) && return true
    has_node_coordinates(a) && return true
    false
end

function get_num_ambient_dims(a)
    hasproperty(a,:num_ambient_dims) && return get_val_parameter(a.num_ambient_dims)
    (has_geometry(a) && has_num_ambient_dims(get_geometry(a))) && return get_num_ambient_dims(get_geometry(a))
    has_node_coordinates(a) && return length(eltype(get_node_coordinates(a)))
    error("num_ambient_dims not found in data")
end

has_node_coordinates(a) = hasproperty(a,:node_coordinates)
get_node_coordinates(a) = a.node_coordinates

has_face_nodes(a) = hasproperty(a,:face_nodes)
get_face_nodes(a) = a.face_nodes
get_face_nodes(a,d) = get_face_nodes(a)[get_val_parameter(d)+1]

has_face_reference_id(a) = hasproperty(a,:face_reference_id)
get_face_reference_id(a) = a.face_reference_id
get_face_reference_id(a,d) = get_face_reference_id(a)[get_val_parameter(d)+1]

has_reference_faces(a) = hasproperty(a,:reference_faces)
get_reference_faces(a) = a.reference_faces
get_reference_faces(a,d) = get_reference_faces(a)[get_val_parameter(d)+1]

function quadrature(geom,degree_per_dir;kwargs...)
    has_geometry(geom) && return quadrature(get_geometry(geom),degree_per_dir;kwargs...)
    if get_is_n_cube(geom)
        if get_is_axis_aligned(geom)
            bounding_box = get_bounding_box(geom)
            D = length(first(bounding_box))
            limits_per_dir = ntuple(i->(bounding_box[1][i],bounding_box[2][i]),Val(D))
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

evaluate(f,x) = f(x)

function shape_functions_from_coordinates(monomial_exponents,node_coordinates)
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

# look first for a field?
function get_is_unit_n_cube(geom)
    !(get_is_axis_aligned(geom) && get_is_n_cube(geom)) && return false
    bounding_box = get_bounding_box(geom)
    all(i->i==0,first(bounding_box)) && all(i->i==1,last(bounding_box))
end

function get_is_unit_simplex(geom)
    !(get_is_axis_aligned(geom) && get_is_simplex(geom)) && return false
    bounding_box = get_bounding_box(geom)
    all(i->i==0,first(bounding_box)) && all(i->i==1,last(bounding_box))
end

function interpolation_from_geometry(
    geom,degree_per_dir;
    type=:default,
    coordinate_type=SVector{length(degree_per_dir),Float64},
    kwargs...)
    f = if get_is_unit_n_cube(geom)
        if type === :default
            e->true
        elseif type === :Q
            e->true
        else
            error("Not implemented")
        end
    elseif get_is_unit_simplex(geom)
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
    shape_functions = shape_functions_from_coordinates(exponents,node_coordinates)
    (;shape_functions,node_coordinates)
end



ms = monomial_exponents(i->true,(2,2))

quad = tensor_product_quadrature((2,2),((0,1),(0,1)))

geo = (;is_n_cube=Val(true), is_axis_aligned=true, bounding_box=[(0,0,0),(2,2,2)])

quad = quadrature(geo,(2,2,3))

display(get_coordinates(quad))
display(get_weights(quad))

geometry = (;num_dims=Val(0))
boundary = nothing
vertex = (;geometry,boundary)

geometry = (;num_dims = 1, is_n_cube=true, is_simplex=true, is_axis_aligned=true, bounding_box=SVector{1,Float64}[(0,),(1,)])

node_coordinates = SVector{1,Float64}[(0,),(1,)]
face_nodes = [[[1],[2]]]
face_reference_id = [[1,1]]
reference_faces = [[vertex]]
boundary = (;node_coordinates,face_nodes,face_reference_id,reference_faces)


exponents = monomial_exponents(i->true,(2,))
node_coordinates = 1.0*exponents
shape_functions = shape_functions_from_coordinates(exponents,node_coordinates)

n = length(node_coordinates)
sq = shape_functions.value!(zeros(n,n),node_coordinates)

display(sq)

sq = shape_functions.gradient!(zeros(eltype(node_coordinates),n,n),node_coordinates)
display(sq)

interpolation = interpolation_from_geometry(geometry,(2,2))

display(interpolation)

segment = (;geometry,boundary,interpolation)

quad = quadrature(segment,(1,))

display(quad)


end # module

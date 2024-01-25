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
has_physical_groups(a) = hasproperty(a,:physical_groups) && a.physical_groups !== nothing
periodic_nodes(a) = a.periodic_nodes
has_periodic_nodes(a) = hasproperty(a,:periodic_nodes) && a.periodic_nodes !== nothing
geometry(a) = a.geometry
topology(a) = a.topology
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
shape_functions(a) = a.shape_functions
tabulation_matrix!(a) = a.tabulation_matrix!
tabulation_matrix(a) = a.tabulation_matrix
order(a) = a.order
monomial_exponents(a) = a.monomial_exponents
lib_to_user_nodes(a) = a.lib_to_user_nodes
interior_nodes(a) = a.interior_nodes
face_permutation_ids(a) = a.face_permutation_ids
face_permutation_ids(a,m,n) = face_permutation_ids(a)[m+1,n+1]
local_nodes(a) = a.local_nodes
local_node_colors(a) = a.local_node_colors
real_type(a) = a.real_type
int_type(a) = a.int_type

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

struct ExtrusionPolytope{D,Tv,Ti} <: GalerkinToolkitDataType
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

abstract type AbstractLagrangeFE <: GalerkinToolkitDataType end

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

#struct LagrangeFE{D,Tv,Ti} <: AbstractLagrangeFE
#    geometry::ExtrusionPolytope{D,Tv,Ti}
#    order_per_dir::NTuple{D,Ti}
#    space::Symbol
#    lib_to_user_nodes::Vector{Ti}
#    major::Symbol
#    shape::ScalarShape
#end
#
#struct TensorValuedLagrangeFE{D,S,Tv,Ti} <: AbstractLagrangeFE
#    geometry::ExtrusionPolytope{D,Tv,Ti}
#    order_per_dir::NTuple{D,Ti}
#    space::Symbol
#    lib_to_user_nodes::Vector{Ti}
#    major::Symbol
#    shape::NTuple{S,Ti}
#end

lagrangian_fe(args...) = GenericLagrangeFE(args...)

#function lagrangian_fe(
#        geometry::ExtrusionPolytope{D,Tv,Ti},
#        order_per_dir::NTuple{D,Ti},
#        space::Symbol,
#        lib_to_user_nodes::Vector{Ti},
#        major::Symbol,
#        shape::ScalarShape) where {D,Tv,Ti}
#
#    LagrangeFE(
#               geometry,
#               order_per_dir,
#               space,
#               lib_to_user_nodes,
#               major,
#               shape)
#end
#
#function lagrangian_fe(
#        geometry::ExtrusionPolytope{D,Tv,Ti},
#        order_per_dir::NTuple{D,Ti},
#        space::Symbol,
#        lib_to_user_nodes::Vector{Ti},
#        major::Symbol,
#        shape::NTuple{S,Ti}) where {D,S,Tv,Ti}
#
#    TensorValuedLagrangeFE(
#               geometry,
#               order_per_dir,
#               space,
#               lib_to_user_nodes,
#               major,
#               shape)
#end

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
    monomial_exponents_from_space(fe.space,fe.order_per_dir,fe.geometry |> real_type)
end

function node_coordinates(fe::AbstractLagrangeFE)
    @assert fe |> geometry |> is_unitary
    mexps = monomial_exponents(fe)
    node_coordinates_from_monomials_exponents(mexps,fe.order_per_dir,fe.geometry |> real_type)
end

function node_coordinates_from_monomials_exponents(monomial_exponents,order_per_dir,real_type)
    node_coordinates = map(monomial_exponents) do exponent
        map(exponent,order_per_dir) do e,order
            if order != 0
                real_type(e/order)
            else
                real_type(e)
            end
        end
    end
end

function monomial_exponents_from_space(space,args...)
    if space == :Q
        monomial_exponents_from_filter((e,o)->true,args...)
    elseif space == :P
        monomial_exponents_from_filter((e,o)->sum(e)<=maximum(o),args...)
    else
        error("Case not implemented (yet)")
    end
end

function monomial_exponents_from_filter(f,order_per_dir,int_type)
    terms_per_dir = Tuple(map(d->d+1,order_per_dir))
    D = length(terms_per_dir)
    cis = CartesianIndices(terms_per_dir)
    m = count(ci->f(Tuple(ci) .- 1,order_per_dir),cis)
    li = 0
    result = zeros(SVector{D,int_type},m)
    for ci in cis
        t = Tuple(ci) .- 1
        if f(t,order_per_dir)
            li += 1
            result[li] = t
        end
    end
    result
end

function tensor_basis(fe::AbstractLagrangeFE)
    Ti = fe |> geometry |> int_type
    if fe.shape == SCALAR_SHAPE
        return Ti(1)
    else
        cis = CartesianIndices(val_parameter(fe.shape))
        l = prod(val_parameter(fe.shape))
        init_tensor = SArray{Tuple{val_parameter(fe.shape)...},Ti}
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


end # module

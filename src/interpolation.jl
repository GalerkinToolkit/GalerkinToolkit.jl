
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

"""
"""
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

abstract type AbstractSpace <: gk.AbstractType end

@enum FreeOrDirichlet FREE=1 DIRICHLET=2

function dofs(a,free_or_diri::FreeOrDirichlet)
    if free_or_diri == FREE
        free_dofs(a)
    else
        dirichlet_dofs(a)
    end
end

function values(a,free_or_diri::FreeOrDirichlet)
    if free_or_diri == FREE
        free_values(a)
    else
        dirichlet_values(a)
    end
end

struct FieldIndex{B}
    field::Int
    dof::B
end

function monolithic_field_major_strategy(field_to_dofs)
    offset = 0
    field_to_offset = map(field_to_dofs) do dofs
        myoffset = offset
        offset += length(dofs)
        return myoffset
    end
    ndofs = offset
    dof_range = Base.OneTo(ndofs)
    allocate_values(::Type{T},r::Base.OneTo) where T = Vector{T}(undef,length(r))
    restrict_values(v,field) = view(v,field_to_dofs[field])
    (;dofs,allocate_values,restrict_values)
end

function cartesian_product(spaces::AbstractSpace...;
        free_values_strategy = gk.monolithic_field_major_strategy(map(gk.free_dofs,spaces)),
        dirichlet_values_strategy = gk.monolithic_field_major_strategy(map(gk.dirichlet_dofs,spaces)),
    )
    CartesianProductSpace(
        spaces,
        free_values_strategy,
        dirichlet_values_strategy,)
end

struct CartesianProductSpace{A,B,C} <: gk.AbstractSpace
    spaces::A
    free_values_strategy::B
    dirichlet_values_strategy::C
end

function ×(a::AbstractSpace,b::AbstractSpace)
    cartesian_product(a,b)
end

function ×(a::CartesianProductSpace,b::AbstractSpace)
    cartesian_product(a.spaces...,b)
end

function ×(a::AbstractSpace,b::CartesianProductSpace)
    cartesian_product(a,b.spaces...)
end

function components(a::CartesianProductSpace)
    a.spaces
end

function component(u::CartesianProductSpace,feild::Integer)
    a.spaces[field]
end

function num_fields(a::CartesianProductSpace)
    length(a.spaces)
end

function domain(a::CartesianProductSpace)
    domains = map(gk.domain,a.spaces)
    domain = first(domains)
    if all(dom == domain,domains)
        domain
    else
        # TODO generalize the concept of AbstractQuantity
        # by introducing PiecewiseQuantity and PiecewiseDomain
        # include also domain id in index
        error("Case not yet implemented")
    end
end

function dof_map(a::CartesianProductSpace,dim)
    terms = map(s->gk.dof_map(s,dim),a.spaces)
    index -> begin
        i = index.dof_per_dim[dim]
        field = i.field
        dof = i.dof
        dof_per_dim = set_at_dim(index.dof_per_dim,Val(dim),dof)
        index2 = replace_dof_per_dim(index,dof_per_dim)
        terms[field](index2)
    end
end

function basis(a::CartesianProductSpace,dim)
    bases = map(s->gk.basis(s,dim),a.spaces)
    prototypes = map(gk.prototype,bases)
    field_to_field = ntuple(identity,gk.num_fields(a))
    domain = gk.domain(a)
    prototype = x->begin
        # TODO if we want to support recursive Cartesian spaces we need a 
        # "MathTuple" instead of a Tuple so that we can call zero(g(x))
        map(field_to_field) do field
            g = prototypes[field]
            g(x)
        end
    end
    gk.quantity(prototype,domain) do index
        i = index.dof_per_dim[dim]
        @assert isa(i,FieldIndex)
        field = i.field
        dof = i.dof
        dof_per_dim = set_at_dim(index.dof_per_dim,Val(dim),dof)
        index2 = replace_dof_per_dim(index,dof_per_dim)
        f = terms[field](index2)
        # TODO if we want to support recursive Cartesian spaces we need a 
        # "MathTuple" instead of a Tuple so that we can call zero(g(x))
        x -> begin
            map(field_to_field) do field_bis
                if field == field_bis
                    f(x)
                else
                    g = prototypes[field_bis]
                    zero(g(x))
                end
            end
        end
    end
end

function free_dofs(a::CartesianProductSpace)
    gk.dofs(a.free_values_strategy)
end

function dirichlet_dofs(a::CartesianProductSpace)
    gk.dofs(a.dirichlet_values_strategy)
end

function discrete_field(space,free_values,dirichlet_values)
    DiscreteField(space,free_values,dirichlet_values)
end

struct DiscreteField{A,B,C} <: gk.AbstractQuantity
    space::A
    free_values::B
    dirichlet_values::C
end

function undef_field(::Type{T},space::AbstractSpace) where T
    free_values = allocate_values(T,free_values_strategy(space))
    dirichet_values = allocate_values(T,dirichlet_values_strategy(space))
    discrete_field(space,free_values,dirichlet_values)
end

function zero_field(::Type{T},space::AbstractSpace) where T
    fill_field(zero(T),space)
end

function fill_field!(u::DiscreteField,z)
    fill!(free_values(u),z)
    fill!(dirichlet_values(u),z)
    u
end

function fill_field(z,space::AbstractSpace)
    u = undef_field(T,space)
    fill_field!(u,z)
end

function num_fields(u::DiscreteField)
    num_fields(gk.space(u))
end

function components(u)
    ntuple(1:gk.num_fields(u)) do field
        component(u,field)
    end
end

function component(u::DiscreteField,feild::Integer)
    space = gk.space(u)
    space_field = component(space,field)
    free_values_field = restrict_values(free_values_strategy(u))(free_values(u),field)
    dirichlet_values_field = restrict_values(dirichlet_values_strategy(u))(dirichlet_values(u),field)
    discrete_field(space_field,free_values_field,dirichlet_values_field)
end

function free_values(x::DiscreteField)
    x.free_values
end

function dirichlet_values(x::DiscreteField)
    x.dirichlet_values
end

function domain(u::DiscreteField)
    gk.domain(gk.space(u))
end

function prototype(u::DiscreteField)
    space = gk.space(u)
    basis = gk.basis(space)
    f = gk.prototype(basis)
    T = eltype(gk.free_values(u))
    z = zero(T)
    # TODO use also dirichlet values type?
    x -> f(x)*z + f(x)*z
end

function term(u::DiscreteField)
    field_to_u = components(u)
    field_to_free_vals = map(free_values,field_to_u)
    field_to_diri_vals = map(dirichlet_values,field_to_u)
    dof_map = gk.dof_map(space)
    dim = 2
    basis = gk.basis(space,dim)
    basis_term = gk.term(basis)
    nfields = num_fields(u)
    index -> begin
        x -> begin
            sum(1:nfields) do field
                sum(1:ndofs) do dof
                    i = FieldIndex(field,dof)
                    dof_per_dim = (nothing,i)
                    index2 = replace_dof_per_dim(index,dof_per_dim)
                    dof = dof_map(index2)
                    f = basis_term(index2)
                    coeff =  dof > 0 ? field_to_free_vals[field][dof] : field_to_diri_vals[field][-dof]
                    f(x)*coeff
                end
            end
        end
    end
end


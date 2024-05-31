
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

Base.iterate(m::AbstractSpace) = iterate(components(m))
Base.iterate(m::AbstractSpace,state) = iterate(components(m),state)
Base.getindex(m::AbstractSpace,field::Integer) = component(m,field)
Base.length(m::AbstractSpace) = num_fields(m)

function free_dofs(a::AbstractSpace)
    gk.dofs(free_values_strategy(a))
end

function dirichlet_dofs(a::AbstractSpace)
    gk.dofs(dirichlet_values_strategy(a))
end

# Single field by default

function num_fields(a::AbstractSpace)
    1
end

function components(a::AbstractSpace)
    (a,)
end

function component(a::AbstractSpace,field)
    @assert field == 1
    a
end

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

struct VectorStrategy{A,B,C}
    dofs::A
    allocate_values::B
    restrict_values::C
end

dofs(a::VectorStrategy) = a.dofs
allocate_values(a::VectorStrategy,args...) = a.allocate_values(args...)
restrict_values(a::VectorStrategy,args...) = a.restrict_values(args...)

function monolithic_field_major_strategy(field_to_dofs)
    offset = 0
    field_to_offset = map(field_to_dofs) do dofs
        myoffset = offset
        offset += length(dofs)
        return myoffset
    end
    ndofs = offset
    dofs = Base.OneTo(ndofs)
    allocate_values(::Type{T}) where T = Vector{T}(undef,ndofs)
    restrict_values(v,field) = view(v,field_to_dofs[field])
    VectorStrategy(dofs,allocate_values,restrict_values)
end

function cartesian_product(spaces::AbstractSpace...;
        free_values_strategy = gk.monolithic_field_major_strategy,
        dirichlet_values_strategy = gk.monolithic_field_major_strategy,
    )
    CartesianProductSpace(
        spaces,
        free_values_strategy(map(gk.free_dofs,spaces)),
        dirichlet_values_strategy(map(gk.dirichlet_dofs,spaces)),)
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

function free_values_strategy(a::CartesianProductSpace)
    a.free_values_strategy
end

function dirichlet_values_strategy(a::CartesianProductSpace)
    a.dirichlet_values_strategy
end

function components(a::CartesianProductSpace)
    a.spaces
end

function component(u::CartesianProductSpace,field::Integer)
    u.spaces[field]
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
        # I am not sure about this. Not needed in practice
        error("Case not yet implemented")
    end
end

function domain(a::CartesianProductSpace,field)
    gk.domain(gk.component(a,field))
end

function num_face_dofs(a::CartesianProductSpace,dim)
    terms = map(s->gk.num_face_dofs(s,dim),a.spaces)
    index -> begin
        field = index.field_per_dim[dim]
        @assert field !== nothing
        index2 = replace_field_per_dim(index,nothing)
        terms[field](index2)
    end
end

function dof_map(a::CartesianProductSpace,dim)
    terms = map(s->gk.dof_map(s,dim),a.spaces)
    index -> begin
        field = index.field_per_dim[dim]
        @assert field !== nothing
        index2 = replace_field_per_dim(index,nothing)
        terms[field](index2)
    end
end

function shape_functions(a::CartesianProductSpace,dim)
    fields = ntuple(identity,gk.num_fields(a))
    map(fields) do field
        gk.shape_functions(a,dim,field)
    end
end

function shape_functions(a::CartesianProductSpace,dim,field)
    f = component(a,field)
    shape_functions(f,dim)
end

function dual_basis(a::CartesianProductSpace)
    error("Not implemented yet. Not needed in practice.")
end

function discrete_field(space,free_values,dirichlet_values)
    DiscreteField(space,free_values,dirichlet_values)
end

struct DiscreteField{A,B,C} <: gk.AbstractQuantity
    space::A
    free_values::B
    dirichlet_values::C
end

Base.iterate(m::DiscreteField) = iterate(components(m))
Base.iterate(m::DiscreteField,state) = iterate(components(m),state)
Base.getindex(m::DiscreteField,field::Integer) = component(m,field)
Base.length(m::DiscreteField) = num_fields(m)

function undef_field(::Type{T},space::AbstractSpace) where T
    free_values = allocate_values(free_values_strategy(space),T)
    dirichlet_values = allocate_values(dirichlet_values_strategy(space),T)
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
    u = undef_field(typeof(z),space)
    fill_field!(u,z)
end

function num_fields(u::DiscreteField)
    num_fields(gk.space(u))
end

function components(u)
    ntuple(gk.num_fields(u)) do field
        component(u,field)
    end
end

function component(u::DiscreteField,field::Integer)
    space = gk.space(u)
    space_field = component(space,field)
    free_values_field = restrict_values(free_values_strategy(space),free_values(u),field)
    dirichlet_values_field = restrict_values(dirichlet_values_strategy(space),dirichlet_values(u),field)
    discrete_field(space_field,free_values_field,dirichlet_values_field)
end

function space(x::DiscreteField)
    x.space
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
    dim = 2
    basis = gk.shape_functions(space,dim)
    f = gk.prototype(basis)
    T = eltype(gk.free_values(u))
    z = zero(T)
    # TODO use also dirichlet values type?
    x -> f(x)*z + f(x)*z
end

function term(u::DiscreteField)
    space = gk.space(u)
    free_vals = gk.free_values(u)
    diri_vals = gk.dirichlet_values(u)
    dim = 2
    dof_map = gk.dof_map(space,dim)
    num_face_dofs = gk.num_face_dofs(space,dim)
    basis = gk.shape_functions(space,dim)
    basis_term = gk.term(basis)
    index -> begin
        x -> begin
            index2 = gk.index(
                              face=index.face,
                              local_face=index.local_face,
                              face_around=index.face_around)
            ndofs = num_face_dofs(index2)
            sum(1:ndofs) do dof
                dof_per_dim = (nothing,dof)
                index3 = replace_dof_per_dim(index2,dof_per_dim)
                dof = dof_map(index3)
                f = basis_term(index3)
                coeff =  dof > 0 ? free_vals[dof] : diri_vals[-dof]
                f(x)*coeff
            end
        end
    end
end

# TODO
# With an extra param we can allow interpolation
# in another domain
# what about coBoundary?
# Or should we detect the domain of f?
#gk.interpolate!(uref,v;domain=gk.domain(uref))
#gk.interpolate!(uref,v;domain_map=)
#gk.interpolate!(uref,v;domain_glue=)
#gk.interpolate!(u,v;domain=Ω)
# We could automatic discover the domain of f
# but we can ask for an addition kwarg domain
# to confirm that we want the domain of f
# and we can also provide another one witht he glue?

function interpolate!(f,u::DiscreteField)
    interpolate!(f,u,nothing)
end

function interpolate_free!(f,u::DiscreteField)
    interpolate!(f,u,FREE)
end

function interpolate_dirichlet!(f,u::DiscreteField)
    # TODO for dirichlet we perhaps want to allow integrating on boundaries?
    interpolate!(f,u,DIRICHLET)
end

function interpolate!(f,u::DiscreteField,free_or_diri::Union{Nothing,FreeOrDirichlet})
    free_vals = gk.free_values(u)
    diri_vals = gk.dirichlet_values(u)
    space = gk.space(u)
    dim = 1
    sigma = gk.dual_basis(space,dim)
    dof_map = gk.dof_map(space,dim)
    vals = sigma(f)
    vals_term = gk.term(vals)
    domain = gk.domain(space)
    num_face_dofs = gk.num_face_dofs(space,dim)
    for face in 1:gk.num_faces(domain)
        index = gk.index(;face)
        for dof in 1:num_face_dofs(index)
            index2 = replace_dof_per_dim(index,(dof,))
            v = vals_term(index2)
            dof = dof_map(index2)
            if dof > 0
                if free_or_diri != DIRICHLET
                    free_vals[dof] = v
                end
            else
                if free_or_diri != FREE
                    diri_vals[-dof] = v
                end
            end
        end
    end
end

# This space is not suitable for assembling problems in general, as there might be gaps in the dofs

function iso_parametric_space(dom::AbstractDomain;
        dirichlet_boundary=nothing,
        free_values_strategy = gk.monolithic_field_major_strategy,
        dirichlet_values_strategy = gk.monolithic_field_major_strategy,
    )
    style = domain_style(dom)
    @assert isa(style,GlobalDomain{true}) "iso_parametric_space only accepts domains in the reference space"
    IsoParametricSpace(dom,dirichlet_boundary,free_values_strategy,dirichlet_values_strategy)
end

struct IsoParametricSpace{A,B,C,D} <: AbstractSpace
    domain::A
    dirichlet_boundary::B
    free_values_strategy::C
    dirichlet_values_strategy::D
end

function free_and_dirichlet_nodes(V::IsoParametricSpace)
    free_and_dirichlet_nodes(V,V.dirichlet_boundary)
end

function free_and_dirichlet_nodes(V::IsoParametricSpace,dirichlet_boundary::Nothing)
    n = V |> gk.domain |> gk.mesh |> gk.num_nodes
    gk.partition_from_mask(fill(true,n))
end

function free_and_dirichlet_nodes(V::IsoParametricSpace,dirichlet_boundary::AbstractDomain)
    mesh = V |> gk.domain |> gk.mesh
    n = mesh |> gk.num_nodes
    node_to_tag = fill(Int32(0),n)
    tag_to_name = dirichlet_boundary |> gk.physical_names
    gk.classify_mesh_nodes!(node_to_tag,mesh,tag_to_name)
    gk.partition_from_mask(i->i==0,node_to_tag)
end

function face_dofs(V::IsoParametricSpace)
    face_dofs(V,V.dirichlet_boundary)
end

function face_dofs(V::IsoParametricSpace,dirichlet_boundary::Nothing)
    domain = gk.domain(V)
    d = gk.face_dim(domain)
    faces = gk.faces(domain)
    face_nodes(gk.mesh(domain),d)
    face_to_nodes = view(face_nodes(gk.mesh(domain),d),faces)
    face_to_nodes
end

function face_dofs(V::IsoParametricSpace,dirichlet_boundary::AbstractDomain)
    face_to_nodes = JaggedArray(gk.face_dofs(V,nothing))
    free_and_dirichlet_nodes = gk.free_and_dirichlet_nodes(V)
    node_to_free_node = gk.permutation(free_and_dirichlet_nodes)
    n_free = length(first(free_and_dirichlet_nodes))
    function node_to_dof(node)
        free_node = node_to_free_node[node]
        if free_node <= n_free
            return Int(free_node)
        end
        Int(-(free_node-n_free))
    end
    face_to_dofs_data = node_to_dof.(face_to_nodes.data)
    face_to_dofs = JaggedArray(face_to_dofs_data,face_to_nodes.ptrs)
    face_to_dofs
end

function reference_fes(V::IsoParametricSpace)
    domain = gk.domain(V)
    d = gk.face_dim(domain)
    # TODO LagrangianMesh face needs to be a AbstractElement
    reference_faces(gk.mesh(domain),d)
end

function face_reference_id(V::IsoParametricSpace)
    domain = gk.domain(V)
    faces = gk.faces(domain)
    d = gk.face_dim(domain)
    face_refid = face_reference_id(gk.mesh(domain),d)
    view(face_refid,faces)
end

function free_values_strategy(a::IsoParametricSpace)
    n = length(first(free_and_dirichlet_nodes(a)))
    a.free_values_strategy( (Base.OneTo(n),) )
end

function dirichlet_values_strategy(a::IsoParametricSpace)
    n = length(last(free_and_dirichlet_nodes(a)))
    a.dirichlet_values_strategy( (Base.OneTo(n),) )
end

function domain(a::IsoParametricSpace)
    a.domain
end

function domain(a::IsoParametricSpace,field)
    @assert field == 1
    a.domain
end

function num_face_dofs(a::IsoParametricSpace,dim)
    face_to_dofs = face_dofs(a)
    index -> begin
        #TODO
        #index.field_per_dim !== nothing && @assert index.field_per_dim[dim] == 1
        face = index.face
        length(face_to_dofs[face])
    end
end

function dof_map(a::IsoParametricSpace,dim)
    face_to_dofs = face_dofs(a)
    index -> begin
        #TODO
        #index.field_per_dim !== nothing && @assert index.field_per_dim[dim] == 1
        face = index.face
        dof = index.dof_per_dim[dim]
        face_to_dofs[face][dof]
    end
end

function shape_functions(a::IsoParametricSpace,dim)
    field=1
    gk.shape_functions(a,dim,field)
end

function shape_function_mask(f,face_around_per_dim,face_around,dim)
    x -> face_around_per_dim[dim] == face_around ? f(x) : zero(f(x))
end

function shape_function_mask(f,face_around_per_dim::Nothing,face_around,dim)
    f
end

function shape_functions(a::IsoParametricSpace,dim,field)
    face_to_refid = face_reference_id(a)
    refid_to_reffes = reference_fes(a)
    refid_to_funs = map(gk.shape_functions,refid_to_reffes)
    domain = gk.domain(a)
    prototype = first(first(refid_to_funs))
    gk.quantity(prototype,domain) do index
        face = index.face
        refid = face_to_refid[face]
        dof = index.dof_per_dim[dim]
        f = refid_to_funs[refid][dof]
        g = shape_function_mask(f,index.face_around_per_dim,index.face_around,dim)
        shape_function_mask(g,index.field_per_dim,field,dim)
    end
end

function dual_basis(a::IsoParametricSpace,dim)
    face_to_refid = face_reference_id(a)
    refid_to_reffes = reference_fes(a)
    refid_to_funs = map(gk.dual_basis,refid_to_reffes)
    domain = gk.domain(a)
    prototype = first(first(refid_to_funs))
    gk.quantity(prototype,domain) do index
        face = index.face
        refid = face_to_refid[face]
        dof = index.dof_per_dim[dim]
        refid_to_funs[refid][dof]
    end
end


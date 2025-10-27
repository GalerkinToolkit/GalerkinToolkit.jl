
struct Field{A,B} <: AbstractField
    quantity::A
    domain::B
end

domain(a::AbstractField) = a.domain

function term(f::AbstractField,i)
    term(f.quantity,i)
end

function term(f::AbstractField)
    term(f.quantity)
end

function analytical_field(f,dom::AbstractDomain;piecewise=Val(false))
    D = num_dims(dom)
    if val_parameter(piecewise)
        mesh = GT.mesh(dom)
        d = num_dims(dom)
        nfaces = num_faces(mesh,d)
        location = zeros(Int32,nfaces)
        group_faces = GT.group_faces(mesh,d)
        for (i,name) in enumerate(group_names(dom))
            faces = group_faces[name]
            location[faces] .= i
        end
        f2(name) = x->f(x,name)
        q = quantity() do opts
            domain = GT.domain(opts)
            index = GT.index(opts)
            domain_face = leaf_term(domain_face_index(index))
            location_term = leaf_term(location)
            domain_face_to_face = call_term(map(leaf_term,(GT.faces,domain))...)
            face = RefTerm(domain_face_to_face,domain_face)
            loc = RefTerm(location_term,face)
            names_term = call_term(map(leaf_term,(GT.group_names,dom))...)
            name = RefTerm(names_term,loc)
            call_term(leaf_term(f2),name)
        end
    else
        q = quantity() do opts
            leaf_term(f)
        end
        location = nothing
    end
    AnalyticalField(f,q,dom,location)
end

struct AnalyticalField{A,B,C,D} <: AbstractField
    definition::A
    quantity::B
    domain::C
    location::D
end

@inline is_piecewise(a::AnalyticalField) = a.location !== nothing

function face_constant_field(data,dom::AbstractDomain)
    mesh = GT.mesh(dom)
    d = num_dims(dom)
    nfaces = num_faces(mesh,d)
    T = eltype(data)
    mesh_data = zeros(T,nfaces)
    mesh_data[GT.faces(dom)] = data
    q = face_quantity(mesh_data,mesh,d)
    q2 = call(a->(x->a),q)
    Field(q2,dom)
end

function piecewise_field(fields::AbstractField...)
    PiecewiseField(fields)
end

struct PiecewiseField{A} <: AbstractType
    fields::A
end

function domain(u::PiecewiseField)
    domains = map(GT.domain,u.fields)
    PiecewiseDomain(domains)
end

function piecewise_domain(domains::AbstractDomain...)
    PiecewiseDomain(domains)
end

struct PiecewiseDomain{A} <: AbstractType
    domains::A
end

function discrete_field(space::AbstractSpace,free_values)
    @assert num_dirichlet_dofs(space) == 0 "This space has Dirichlet DOFs. You should provide dirichlet_values."
    dirichlet_values = similar(free_values,dirichlet_dofs(space))
    discrete_field(space,free_values,dirichlet_values)
end

function discrete_field(space::AbstractSpace,free_values,dirichlet_values)
    mesh = space |> GT.mesh
    workspace = setup_discrete_field_workspace(space,free_values,dirichlet_values)
    DiscreteField(mesh,space,free_values,dirichlet_values,workspace)
end

function setup_discrete_field_workspace(space,fv,dv)
    nothing
end

struct DiscreteField{A,B,C,D,W} <: GT.AbstractField
    mesh::A
    space::B
    free_values::C
    dirichlet_values::D
    workspace::W
end

# term is defined in compiler.jl

#Moved to CartesianProductSpace
#Base.iterate(m::DiscreteField) = iterate(fields(m))
#Base.iterate(m::DiscreteField,state) = iterate(fields(m),state)
#Base.getindex(m::DiscreteField,field::Integer) = field(m,field)
#Base.length(m::DiscreteField) = num_fields(m)

# @enum FreeOrDirichlet FREE=1 DIRICHLET=2 FREE_AND_DIRICHLET=3
abstract type FreeOrDirichlet <: AbstractType end
struct EnumFree <: FreeOrDirichlet end # fix type instability by making the enum different types
struct EnumDirichlet <: FreeOrDirichlet end
struct EnumFreeAndDirichlet <: FreeOrDirichlet end
const FREE = EnumFree()
const DIRICHLET = EnumDirichlet()
const FREE_AND_DIRICHLET = EnumFreeAndDirichlet()

function values(a::DiscreteField,free_or_diri::FreeOrDirichlet)
    @assert free_or_diri != FREE_AND_DIRICHLET
    if free_or_diri == FREE
        free_values(a)
    else
        dirichlet_values(a)
    end
end

function allocate_values(::Type{T},dofs) where T
    n = length(dofs)
    Vector{T}(undef,n)
end

function allocate_values(::Type{T},dofs::PRange) where T
    pzeros(T,dofs)
end

function allocate_values(::Type{T},dofs::BRange) where T
    map(dofs.blocks) do dofs_i
        allocate_values(T,dofs_i)
    end |> BVector
end

function rand_values(::Type{T},dofs) where T
    n = length(dofs)
    rand(T,n)
end

function rand_values(::Type{T},dofs::BRange) where T
    map(dofs.blocks) do dofs_i
        rand_values(T,dofs_i)
    end |> BVector
end

function constant_values(f,dofs)
    n = length(dofs)
    Fill(f,n)
end

function constant_values(f,dofs::BRange)
    map(dofs.blocks) do dofs_i
        constant_values(f,dofs_i)
    end |> BVector
end

function undef_field(::Type{T},space::AbstractSpace) where T
    free_values = allocate_values(T,GT.free_dofs(space))
    dirichlet_values = allocate_values(T,GT.dirichlet_dofs(space))
    discrete_field(space,free_values,dirichlet_values)
end

function rand_field(::Type{T},space::AbstractSpace) where T
    free_values = rand_values(T,GT.free_dofs(space))
    dirichlet_values = zeros(T,GT.dirichlet_dofs(space))
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

function undef_dirichlet_field(::Type{T},space::AbstractSpace) where T
    free_values = constant_values(zero(T),GT.free_dofs(space))
    dirichlet_values = allocate_values(T,GT.dirichlet_dofs(space))
    discrete_field(space,free_values,dirichlet_values)
end

function zero_dirichlet_field(::Type{T},space::AbstractSpace) where T
    uh = undef_dirichlet_field(T,space)
    fill!(dirichlet_values(uh),zero(T))
    uh
end

function dirichlet_field(space::AbstractSpace,dirichlet_values::AbstractVector)
    T = eltype(dirichlet_values)
    free_values = constant_values(zero(T),GT.free_dofs(space))
    discrete_field(space,free_values,dirichlet_values)
end

function num_fields(u::DiscreteField)
    num_fields(GT.space(u))
end

function fields(u)
    ntuple(GT.num_fields(u)) do field
        GT.field(u,field)
    end
end

function field(u::DiscreteField,field::Integer)
    space = GT.space(u)
    space_field = GT.field(space,field)
    free_values_field = view(free_values(u),Block(field))
    dirichlet_values_field = view(dirichlet_values(u),Block(field))
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
    GT.domain(GT.space(u))
end

#function prototype(u::DiscreteField)
#    GT.prototype(u.quantity)
#    #space = GT.space(u)
#    #dim = 2
#    #basis = GT.shape_functions(space,dim)
#    #f = GT.prototype(basis)
#    #T = eltype(GT.free_values(u))
#    #z = zero(T)
#    ## TODO use also dirichlet values type?
#    #x -> f(x)*z + f(x)*z
#end

#function term(u::DiscreteField,i...)
#    GT.term(u.quantity,i...)
#    #space = GT.space(u)
#    #free_vals = GT.free_values(u)
#    #diri_vals = GT.dirichlet_values(u)
#    #space = GT.space(u)
#    #free_vals = GT.free_values(u)
#    #diri_vals = GT.dirichlet_values(u)
#    #dim = 2
#    #dof_map = GT.dof_map(space,dim)
#    #num_face_dofs = GT.num_face_dofs(space,dim)
#    #basis = GT.shape_functions(space,dim)
#    #basis_term = GT.term(basis)
#    #index -> begin
#    #    x -> begin
#    #        index2 = GT.index(
#    #                          face=index.face,
#    #                          local_face=index.local_face,
#    #                          face_around=index.face_around)
#    #        ndofs = num_face_dofs(index2)
#    #        sum(1:ndofs) do dof
#    #            dof_per_dim = (nothing,dof)
#    #            index3 = replace_dof_per_dim(index2,dof_per_dim)
#    #            dof = dof_map(index3)
#    #            f = basis_term(index3)
#    #            coeff =  dof > 0 ? free_vals[dof] : diri_vals[-dof]
#    #            f(x)*coeff
#    #        end
#    #    end
#    #end
#end

# TODO
# With an extra param we can allow interpolation
# in another domain
# what about coBoundary?
# Or should we detect the domain of f?
#GT.interpolate!(uref,v;domain=GT.domain(uref))
#GT.interpolate!(uref,v;domain_map=)
#GT.interpolate!(uref,v;domain_glue=)
#GT.interpolate!(u,v;domain=Ω)
# We could automatic discover the domain of f
# but we can ask for an addition kwarg domain
# to confirm that we want the domain of f
# and we can also provide another one witht he glue?
# Solution: Two options: u is defined either on the domain of the space
# or on the Dirichlet boundary


#function interpolate!(f,u::DiscreteField;free_or_dirichlet=FREE_AND_DIRICHLET)
#    interpolate!(f,u;free_or_dirichlet)
#end


function interpolate(f,space::AbstractSpace;free_or_dirichlet=FREE_AND_DIRICHLET)
    sigma = GT.dual_basis_quantity(space)
    vals = sigma(f)
    index = GT.index(Val(1))
    domain = GT.domain(space)
    opts = QuantityOptions(domain,index)
    t = term(vals,opts)
    T = typeof(prototype(t))
    u = zero_field(T,space)
    interpolate!(f,u;free_or_dirichlet)
    u
end

#function interpolate!(f,u::DiscreteField,field)
#    interpolate!(f,u,nothing,field)
#end

function interpolate_free!(f,u::DiscreteField)
    interpolate!(f,u;free_or_dirichlet = FREE)
end

function interpolate_free(f,space::AbstractSpace)
    interpolate(f,space;free_or_dirichlet = FREE)
end

#function interpolate_free!(f,u::DiscreteField,field)
#    interpolate!(f,u,FREE,field)
#end

function interpolate_dirichlet!(f,u::DiscreteField)
    # TODO for dirichlet we perhaps want to allow integrating on boundaries?
    interpolate!(f,u;free_or_dirichlet = DIRICHLET)
end

function interpolate_dirichlet(f,space::AbstractSpace)
    interpolate(f,space;free_or_dirichlet = DIRICHLET)
end

#function interpolate_dirichlet!(f,u::DiscreteField,field)
#    # TODO for dirichlet we perhaps want to allow integrating on boundaries?
#    interpolate!(f,u,DIRICHLET,field)
#end

function interpolate!(f,u::DiscreteField;free_or_dirichlet=FREE_AND_DIRICHLET)
    interpolate_impl!(f,u,space(u),free_or_dirichlet)
end

#function interpolate!(f,u::DiscreteField,free_or_diri::Union{Nothing,FreeOrDirichlet},field)
#    ui = GT.field(u,field)
#    interpolate_impl!(f,ui,free_or_diri)
#    u
#end

#function interpolate_impl!(f,u,space,free_or_diri;location=1)
#    free_vals = GT.free_values(u)
#    diri_vals = GT.dirichlet_values(u)
#    dof = gensym("fe-dof")
#    sigma = GT.dual_basis_quantity(space,dof)
#    face_to_dofs = GT.face_dofs(space)
#    vals = sigma(f)
#    domain = GT.domain(space)
#    dirichlet_dof_location = GT.dirichlet_dof_location(space)
#    index = generate_index(domain)
#    sface_to_face = get_symbol!(index,faces(domain),"sface_to_face")
#    t = term(vals,index)
#    expr_qty = t |> expression |> simplify
#    D = num_dims(domain)
#    face = face_index(index,D)
#    s_qty = GT.topological_sort(expr_qty,(face,dof))
#    expr = quote
#        (args,storage) -> begin
#            (;face_to_dofs,free_vals,diri_vals,nfaces,dirichlet_dof_location,location,free_or_diri) = args
#            $(unpack_index_storage(index,:storage))
#            $(s_qty[1])
#            for sface in 1:nfaces
#                $face = $sface_to_face[sface]
#                $(s_qty[2])
#                dofs = face_to_dofs[$face]
#                ndofs = length(dofs)
#                for $dof in 1:ndofs
#                    v = $(s_qty[3])
#                    gdof = dofs[$dof]
#                    if gdof > 0
#                        if free_or_diri != DIRICHLET
#                            free_vals[gdof] = v
#                        end
#                    else
#                        diri_dof = -gdof
#                        if free_or_diri != FREE && dirichlet_dof_location[diri_dof] == location
#                            diri_vals[diri_dof] = v
#                        end
#                    end
#                end
#            end
#        end
#    end
#    loop! = eval(expr)
#    storage = GT.index_storage(index)
#    nfaces = GT.num_faces(domain)
#    args = (;face_to_dofs,free_vals,diri_vals,nfaces,dirichlet_dof_location,location,free_or_diri)
#    invokelatest(loop!,args,storage)
#    u
#end

function interpolate_impl!(f::PiecewiseField,u,space,free_or_diri)
    for (location,field) in f.fields |> enumerate
        interpolate_impl!(field,u,space,free_or_diri;location)
    end
    u
end

function semi_discrete_field(T,V::AbstractSpace)
    semi_discrete_field(T,V) do t,uh
        uh
    end
end

function semi_discrete_field(f,T,V::AbstractSpace)
    uh = zero_field(T,V)
    semi_discrete_field(f,uh)
end

function semi_discrete_field(uh::DiscreteField)
    semi_discrete_field(uh) do t,uh
        uh
    end
end

function semi_discrete_field(f,uh::DiscreteField)
    SemiDiscreteField(f,uh)
end

struct SemiDiscreteField{A,B}
    update::A
    discrete_field::B
end

free_values(uh::SemiDiscreteField) = free_values(uh.discrete_field)
dirichlet_values(uh::SemiDiscreteField) = dirichlet_values(uh.discrete_field)

function (u::SemiDiscreteField)(t)
    u.update(t,u.discrete_field)
    u.discrete_field
end

function space(u::SemiDiscreteField)
    space(u.discrete_field)
end

function discrete_field(u::SemiDiscreteField)
    u.discrete_field
end

function face_diameter_field(Ω::AbstractDomain)
    dΩ = GT.quadrature(Ω,0)
    dims = map(diameter,each_face(dΩ))
    face_constant_field(dims,Ω)
end

function interpolate!(f,u::SemiDiscreteField;free_or_diri=FREE_AND_DIRICHLET)
    interpolate!(f,u.discrete_field;free_or_dirichlet)
end

function interpolate_free!(f,u::SemiDiscreteField)
    interpolate_free!(f,u.discrete_field)
end

function interpolate_dirichlet!(f,u::SemiDiscreteField)
    interpolate_dirichlet!(f,u.discrete_field)
end



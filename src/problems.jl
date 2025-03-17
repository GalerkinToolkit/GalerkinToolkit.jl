
# Measure

function measure(domain::AbstractDomain,degree)
    measure(quadrature,domain,degree)
end

function measure(f,domain::AbstractFaceDomain,degree)
    f(domain,degree)
end

function measure(f,domain::AbstractMeshDomain,degree)
    mesh = GT.mesh(domain)
    d = GT.num_dims(domain)
    dface_to_drefid = GT.face_reference_id(mesh,d)
    drefid_refdface = GT.reference_spaces(mesh,Val(d))
    refid_to_quad = map(drefid_refdface) do refdface
        geo = GT.domain(refdface)
        f(geo,degree)
    end
    mesh_quadrature(;
        domain,
        face_reference_id=dface_to_drefid,
        reference_quadratures=refid_to_quad)
end

# Integrals

function integrate(f,quadrature::AbstractQuadrature)
    coefficient = 1
    contribution = DomainContribution(f,quadrature,coefficient)
    integral(contribution)
end
const ∫ = integrate

struct DomainContribution{A,B,C}
    integrand::A
    quadrature::B
    coefficient::C
end

domain(a::DomainContribution) = domain(a.quadrature)
quadrature(a::DomainContribution) = a.quadrature
coefficient(a::DomainContribution) = a.coefficient

function prototype(a::DomainContribution)
    prototype(optimize(term(a,index(Val(0)))))
end

function replace_coefficient(a::DomainContribution,coefficient)
    DomainContribution(a.integrand,a.quadrature,coefficient)
end

function quantity(contribution::DomainContribution)
    (;integrand,quadrature) = contribution
    domain = GT.domain(contribution)
    x = coordinate_quantity(quadrature)
    dV = weight_quantity(quadrature)
    coefficient = GT.coefficient(contribution)
    alpha = uniform_quantity(coefficient)
    alpha*integrand(x)*dV
end

function term(contribution::DomainContribution,index)
    domain = GT.domain(contribution)
    opts = QuantityOptions(domain,index)
    term(quantity(contribution),opts)
end

function integral(contributions...)
    Integral(contributions)
end

struct Integral{A}
    contributions::A
end

contributions(i::Integral) = contributions

function contribution(i::Integral,domain::AbstractDomain)
    for contribution in contributions(i)
        if domain == GT.comain(contribution)
            return contribution
        end
    end
    error("domain is not found")
end

function Base.:+(int1::Integral,int2::Integral)
    # TODO merge contributions on the same domain?
    contribs = (GT.contributions(int1)...,GT.contributions(int2)...)
    Integral(contribs)
end

function Base.:*(v::Number,int::Integral)
    contribs = map(GT.contributions(int)) do contribution
        coefficient = v*GT.coefficient(contribution)
        replace_coefficient(contribution,coefficient)
    end
    Integral(contribs)
end

function Base.:-(int1::Integral,int2::Integral)
    int1 + (-1)*int2
end

function Base.:-(int1::Real,int2::Integral)
    @assert int1 == 0
    (-1)*int2
end

function Base.:*(int::Integral,v::Number)
    v*int
end

function Base.:/(int::Integral,v::Number)
    (1/v)*int
end

# Sample on the quadrature points. This is needed to visualize fields

function sample_accessor(f,quadrature::AbstractQuadrature;
        parameters = (),
        reuse = isempty(parameters) ? Val(false) : Val(true),
    )
    g = generate_sample_accessor(f,quadrature;parameters)
    accessor = g(parameters...)
    if val_parameter(reuse)
        accessor, g
    else
        accessor
    end
end

function update_sample_accessor(g;parameters)
    g(parameters...)
end

# 0-forms

function assemble_scalar(integral::Integral;
        parameters = (),
        reuse = isempty(parameters) ? Val(false) : Val(true),
    )
    contributions = GT.contributions(integral)
    pairs = map(contributions) do contribution
        init = zero(prototype(contribution))
        params_loop = generate_assemble_scalar(contribution;parameters)
        loop = params_loop(parameters...)
        loop(init), (params_loop, init)
    end
    b = sum(first,pairs)
    loops = map(last,pairs)
    if val_parameter(reuse)
        b, loops
    else
        b
    end
end

function update_scalar(loops;parameters)
    sum(loops) do (params_loop, init)
        loop = params_loop(parameters...)
        loop(init)
    end
end

function Base.sum(int::Integral)
    assemble_scalar(int)
end


function assemble_face_contribution!(face_to_val,contribution::DomainContribution;
        parameters = (),
        reuse = isempty(parameters) ? Val(false) : Val(true),
    )
    params_loop = generate_assemble_face_contribution(contribution;parameters)
    loop! = params_loop(parameters...)
    loop!(face_to_val)
    if val_parameter(reuse)
        face_to_val, params_loop
    else
        face_to_val
    end
end

function update_face_contribution!(face_to_val, params_loop;parameters)
    loop! = params_loop(parameters...)
    loop!(face_to_val)
    face_to_val
end

function face_contribution(int::Integral,domain::AbstractDomain)
    contribution = GT.contribution(int,domain)
    nfaces = GT.domain(contribution)
    T = typeof(prototype(contribution))
    face_to_val = zeros(T,nfaces)
    assemble_face_contribution!(face_to_val,contribution)
end

function face_contribution!(a,int::Integral,domain::AbstractDomain)
    contribution = GT.contribution(int,domain)
    assemble_face_contribution!(face_to_val,contribution)
end

# 1-forms

function assemble_vector(f,::Type{T},space;
    parameters = (),
    reuse = isempty(parameters) ? Val(false) : Val(true),
    vector_strategy = monolithic_vector_assembly_strategy(),
    free_or_dirichlet = FREE,
    ) where T

    arg = 1
    dv = form_argument_quantity(space,arg)
    integral = f(dv)
    contributions = GT.contributions(integral)
    domains = map(GT.domain,contributions)
    alloc = allocate_vector(T,space,domains...;vector_strategy,free_or_dirichlet)
    loops = map(contributions) do contribution
        params_loop = generate_assemble_vector(contribution,space;parameters)
        loop! = params_loop(parameters...)
        loop!(alloc)
        params_loop
    end
    b, vector_cache = compress(alloc;reuse=Val(true))
    cache = (;loops,alloc,vector_cache)
    if val_parameter(reuse)
        b, cache
    else
        b
    end
end

function update_vector!(b,cache;parameters)
    cache = (;loops,alloc,vector_cache)
    reset!(alloc)
    map(loops) do params_loop
        loop! = params_loop(parameters...)
        loop!(alloc)
    end
    compress!(alloc,b,vector_cache)
    b
end

# 1-forms (dual)

function assemble_values!(f,uh::DiscreteField,space::AbstractSpace;
    parameters = (),
    reuse = isempty(parameters) ? Val(false) : Val(true),
    free_or_dirichlet = FREE_AND_DIRICHLET, 
    location=1,
    )

    ds = dual_operator_quantity(space)
    qty = ds(f)
    params_loop = generate_assemble_values(space;parameters,free_or_dirichlet,location) do index
        opts = QuantityOptions(domain,index)
        term(qty,opts)
    end
    loop! = params_loop(parameters...)
    loop!(uh)
    if val_parameter(reuse)
        uh, params_loop
    else
        uh
    end
end

function update_values!(uh,params_loop;parameters)
    loop! = params_loop(parameters...)
    loop!(uh)
    uh
end

# 2-forms

function assemble_matrix(f,::Type{T},trial_space,test_space;
    parameters = (),
    reuse = isempty(parameters) ? Val(false) : Val(true),
    matrix_strategy = monolithic_matrix_assembly_strategy(),
    free_or_dirichlet = (FREE,FREE),
    ) where T
    arg_test = 1
    arg_trial = 2
    du = form_argument_quantity(space,arg_trial)
    dv = form_argument_quantity(space,arg_test)
    integral = f(du,dv)
    contributions = GT.contributions(integral)
    domains = map(GT.domain,contributions)
    alloc = allocate_matrix(T,space,domains...;matrix_strategy,free_or_dirichlet)
    loops = map(contributions) do contribution
        params_loop = generate_assemble_matrix(contribution,space;parameters)
        loop! = params_loop(parameters...)
        loop!(alloc)
        params_loop
    end
    b, matrix_cache = compress(alloc;reuse=Val(true))
    cache = (;loops,alloc,matrix_cache)
    if val_parameter(reuse)
        b, cache
    else
        b
    end
end

function update_matrix!(A,cache;parameters)
    cache = (;loops,alloc,vector_cache)
    reset!(alloc)
    map(loops) do params_loop
        loop! = params_loop(parameters...)
        loop!(alloc)
    end
    compress!(alloc,b,vector_cache)
    b
end

function assemble_matrix_with_free_and_dirichlet_columns(f,::Type{T},trial_space,test_space;
    parameters = (),
    reuse = isempty(parameters) ? Val(false) : Val(true),
    matrix_strategy = monolithic_matrix_assembly_strategy(),
    free_or_dirichlet = FREE, # Rows
    ) where T
    # We are doing two calls, this can be optimized
    A,Acache = assemble_matrix(f,T,trial_space,test_space;
        parameters,reuse=Val(true),matrix_strategy,free_or_dirichlet=(free_or_dirichlet,FREE))
    Ad,Adcache = assemble_matrix(f,T,trial_space,test_space;
        parameters,reuse=Val(true),matrix_strategy,free_or_dirichlet=(free_or_dirichlet,DIRICHLET))
    cache = (Acache,Adcache)
    if val_parameter(reuse)
        A,Ad,cache
    else
        A,Ad
    end
end

function update_matrix_with_free_and_dirichlet_columns!(A,Ad,cache;parameters)
    (Acache,Adcache) = cache
    update_matrix!(A,Acache;parameters)
    update_matrix!(Ad,Adcache;parameters)
    A,Ad
end

# 1 and 2 forms.

function assemble_matrix_and_vector(a,l,::Type{T},U,V;
    parameters = (),
    reuse = isempty(parameters) ? Val(false) : Val(true),
    matrix_strategy = monolithic_matrix_assembly_strategy(),
    vector_strategy = monolithic_vector_assembly_strategy(),
    ) where T
    A,matrix_cache = assemble_matrix(a,T,U,V;parameters,reuse=Val(true),matrix_strategy)
    b,vector_cache = assemble_vector(l,T,V;parameters,reuse=Val(true),vector_strategy)
    cache = (;matrix_cache,vector_cache)
    if val_parameter(reuse)
        A,b,cache
    else
        A,b
    end
end

function update_matrix_and_vector!(A,b,cache;parameters)
    (;matrix_cache,vector_cache) = cache
    update_matrix!(A,matrix_cache;parameters)
    update_vector!(b,vector_cache;parameters)
    A,b
end

function assemble_matrix_and_vector_with_free_and_dirichlet_columns(a,l,::Type{T},U,V;
    parameters = (),
    reuse = isempty(parameters) ? Val(false) : Val(true),
    matrix_strategy = monolithic_matrix_assembly_strategy(),
    vector_strategy = monolithic_vector_assembly_strategy(),
    free_or_dirichlet = FREE, #Rows
    ) where T
    A,Ad,matrix_cache = assemble_matrix_with_free_and_dirichlet_columns(
        a,T,U,V;parameters,reuse=Val(true),matrix_strategy,free_or_dirichlet)
    b,vector_cache = assemble_vector(l,T,V;parameters,reuse=Val(true),vector_strategy,free_or_dirichlet)
    cache = (;matrix_cache,vector_cache)
    if val_parameter(reuse)
        A,Ad,b,cache
    else
        A,Ad,b
    end
end

function update_matrix_and_vector_with_free_and_dirichlet_columns!(A,Ad,b,cache;parameters)
    (;matrix_cache,vector_cache) = cache
    update_matrix_with_free_and_dirichlet_columns!(A,Ad,matrix_cache;parameters)
    update_vector!(b,vector_cache;parameters)
    A,b
end

# Linear problems

function linear_problem(uhd::DiscreteField,a,l,V=GT.space(uhd);
        matrix_strategy = monolithic_matrix_assembly_strategy(),
        vector_strategy = monolithic_vector_assembly_strategy(),
    )
    U = GT.space(uhd)
    xd = dirichlet_values(uhd)
    T = eltype(xd)
    A,Ad,b = assemble_matrix_and_vector_with_free_and_dirichlet_columns(a,l,T,U,V;matrix_strategy,vector_strategy)
    mul!(b,Ad,xd,-1,1)
    x = similar(b,axes(A,2))
    PS.linear_problem(x,A,b)
end

function linear_problem(::Type{T},U::AbstractSpace,a,l,V=U;
        matrix_strategy = monolithic_matrix_assembly_strategy(),
        vector_strategy = monolithic_vector_assembly_strategy(),
    ) where T
    A,b = assemble_matrix_and_vector(a,l,T,U,V;matrix_strategy,vector_strategy)
    x = similar(b,axes(A,2))
    PS.linear_problem(x,A,b)
end

# Nonlinear problems

function nonlinear_problem(uh0::DiscreteField,r,j,V=GT.space(uh0),
        matrix_strategy = monolithic_matrix_assembly_strategy(),
        vector_strategy = monolithic_vector_assembly_strategy(),
    )
    U = GT.space(uh0)
    #b0,residual_cache = assemble_residual(r,uh,V;reuse=Val(true),vector_strategy)
    #A0,jacobian_cache = assemble_jacobian(j,uh,V,reuse=Val(true),matrix_strategy)
    parameters = (uh,)
    b0,residual_cache = assemble_vector(r(uh),V;parameters,vector_strategy)
    A0,jacobian_cache = assemble_matrix(j(uh),U,V;parameters,matrix_strategy)
    workspace = (;uh0,residual_cache,jacobian_cache)
    PS.nonlinear_problem(nonlinear_problem_update,x0,b0,A0,workspace)
end

function nonlinear_problem_update(p)
    (;uh0,residual_cache,jacobian_cache) = PS.workspace(p)
    x = PS.solution(p)
    uh = solution_field(uh0,x)
    A = PS.jacobian(p)
    b = PS.residual(p)
    parameters = (uh,)
    if b !== nothing && A !== nothing
        update_vector!(b,residual_cache;parameters)
        update_matrix!(A,jacobian_cache;parameters)
    elseif b !== nothing && A === nothing
        update_vector!(b,residual_cache;parameters)
    elseif b === nothing && A !== nothing
        update_matrix!(A,jacobian_cache;parameters)
    end
    p = PS.update(p,solution=x,residual=b,jacobian=A)
end

# Solution fields

function solution_field(U::AbstractSpace,x::AbstractVector)
    T = eltype(x)
    uhd = zero_dirichlet_field(T,U)
    solution_field(uhd,x)
end

function solution_field(uhd::DiscreteField,x::AbstractVector)
    diri_vals = dirichlet_values(uhd)
    U = GT.space(uhd)
    free_vals = free_values_from_solution(x,free_dofs(U))
    discrete_field(U,free_vals,diri_vals)
end

function solution_field!(uh::DiscreteField,x::AbstractVector)
    U = GT.space(uh)
    free_vals_src = free_values_from_solution(x,free_dofs(U))
    free_vals_dest = free_values(uh)
    copyto!(free_vals_dest,free_vals_src)
    uh
end

function solution_field(U::AbstractSpace,p::PS.AbstractProblem)
    solution_field(U,PS.solution(p))
end

function solution_field(U::DiscreteField,p::PS.AbstractProblem)
    solution_field(U,PS.solution(p))
end

function solution_field(U::AbstractSpace,p::PS.AbstractSolver)
    solution_field(U,PS.solution(p))
end

function solution_field(U::DiscreteField,p::PS.AbstractSolver)
    solution_field(U,PS.solution(p))
end

function solution_field!(uh::DiscreteField,p::PS.AbstractProblem)
    solution_field!(uh,PS.solution(p))
end

function solution_field!(uh::DiscreteField,p::PS.AbstractSolver)
    solution_field!(uh,PS.solution(p))
end

function free_values_from_solution(x,dofs)
    x
end

function free_values_from_solution(x,dofs::BRange)
    nfields = blocklength(dofs)
    map(1:nfields) do field
        pend = blocklasts(dofs)[field]
        pini = 1 + pend - length(blocks(dofs)[field])
        view(x,pini:pend)
    end |> BVector
end

function free_values_from_solution(x::BVector,dofs::BRange)
    x
end


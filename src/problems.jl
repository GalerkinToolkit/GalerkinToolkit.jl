
function integrate(f,quadrature::AbstractQuadrature)
    contribution = DomainContribution(f,quadrature)
    integral(contribution)
end
const ∫ = integrate

struct DomainContribution{A,B}
    integrand::A
    quadrature::B
end

domain(a::DomainContribution) = domain(a.quadrature)

function term(contribution::DomainContribution,index)
    (;integrand,quadrature) = contribution
    domain = GT.domain(contribution)
    x = coordinate_quantity(quadrature)
    dV =weight_quantity(quadrature)
    quantity = call(*,integrand(x),dV)
    opts = QuantityOptions(domain,index)
    term(quantity,opts)
end

function integral(contributions...)
    Integral(contributions)
end

struct Integral{A}
    contributions::A
end

contributions(i::Integral) = contributions

function nonlinear_problem(uh0::DiscreteField,r,j,V=GT.space(uh0),
        matrix_strategy = monolithic_matrix_assembly_strategy(),
        vector_strategy = monolithic_vector_assembly_strategy(),
    )
    U = GT.space(uh0)
    b0,residual_cache = assemble_residual(r,uh,V;reuse=Val(true),vector_strategy)
    A0,jacobian_cache = assemble_jacobian(j,uh,V,reuse=Val(true),matrix_strategy)
    workspace = (;uh0,residual_cache,jacobian_cache)
    PS.nonlinear_problem(nonlinear_problem_update,x0,b0,A0,workspace)
end

function nonlinear_problem_update(p)
    (;uh0,residual_cache,jacobian_cache) = PS.workspace(p)
    x = PS.solution(p)
    uh = solution_field(uh0,x)
    A = PS.jacobian(p)
    b = PS.residual(p)
    if b !== nothing && A !== nothing
        assemble_residual!(b,uh,residual_cache)
        assemble_jacobian!(A,uh,jacobian_cache)
    elseif b !== nothing && A === nothing
        assemble_residual!(b,uh,residual_cache)
    elseif b === nothing && A !== nothing
        assemble_jacobian!(A,uh,jacobian_cache)
    end
    p = PS.update(p,solution=x,residual=b,jacobian=A)
end

function assemble_residual(r,uh::DiscreteField,V=GT.space(uh),
        vector_strategy = monolithic_vector_assembly_strategy(),
        reuse = Val(false),
    )
    l = r(uh)
    contributions = GT.contributions(l)
    domains = map(GT.domain,contributions)
    x = free_values(uh)
    T = eltype(x)
    alloc = allocate_vector(T,V,domains..;vector_strategy)
    loops = map(contributions) do contribution
        uh_loop = generate_vector_assembly_loop(contribution,V,uh)
        loop! = uh_loop(uh)
        loop!(alloc)
        uh_loop
    end
    b, vector_cache = compress(alloc;reuse=Val(true))
    cache = (;loops,alloc,vector_cache)
    if val_parameter(reuse)
        b, cache
    else
        b
    end
end

function assemble_residual!(b,uh,cache)
    cache = (;loops,alloc,vector_cache)
    reset!(alloc)
    map(loops) do uh_loop
        loop! = uh_loop(uh)
        loop!(alloc)
    end
    compress!(alloc,b,vector_cache)
    b
end


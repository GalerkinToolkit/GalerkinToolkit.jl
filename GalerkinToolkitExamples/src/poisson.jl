"""
This implements a Poisson solver with several methods
"""
module Poisson

import GalerkinToolkit as GT
using GalerkinToolkit: ∫
import PartitionedSolvers as ps
import ForwardDiff
using StaticArrays
using LinearAlgebra
using SparseArrays
using PartitionedArrays
using TimerOutputs

function main(params_in)
    params_default = default_params()
    params = add_default_params(params_in,params_default)
    if params[:implementation] === :automatic
        main_automatic(params)
    elseif params[:implementation] === :hand_written
        main_hand_written(params)
    else
        error(" implementation should be either :automatic or :hand_written")
    end
end

gradient(u) = x->ForwardDiff.gradient(u,x)
jacobian(u) = x->ForwardDiff.jacobian(u,x)
laplacian(u) = x-> tr(ForwardDiff.jacobian(y->ForwardDiff.gradient(u,y),x))

Δ(u) = GT.call(laplacian,u)
∇(u) = GT.call(gradient,u)
Δ(u,x) = Δ(u)(x)
∇(u,x) = ∇(u)(x)

mean(u,x) = 0.5*(u(x)[+]+u(x)[-])
jump(u,n,x) = u(x)[+]*n[+](x) + u(x)[-]*n[-](x)

function main_automatic(params)
    timer = params[:timer]
    results = Dict{Symbol,Any}()

    mesh = params[:mesh]
    D = GT.num_dims(mesh)
    Ω = GT.interior(mesh,physical_names=params[:domain_tags])
    Γd = GT.boundary(mesh;physical_names=params[:dirichlet_tags])
    Γn = GT.boundary(mesh;physical_names=params[:neumann_tags])
    integration_degree = params[:integration_degree]
    dΩ = GT.measure(Ω,integration_degree)
    dΓn = GT.measure(Γn,integration_degree)

    u = GT.analytical_field(params[:u],Ω)
    f(x) = -Δ(u,x)
    n_Γn = GT.unit_normal(Γn,Ω)
    g(x) = n_Γn(x)⋅∇(u,x)

    interpolation_degree = params[:interpolation_degree]
    γ = integration_degree*(integration_degree+1)
    γ = γ/10.0

    @assert params[:discretization_method] in (:continuous_galerkin,:interior_penalty)

    if params[:discretization_method] !== :continuous_galerkin
        conformity = :L2
        GT.label_interior_faces!(mesh;physical_name="__INTERIOR_FACES__")
        Λ = GT.skeleton(mesh;physical_names=["__INTERIOR_FACES__"])
        dΛ = GT.measure(Λ,integration_degree)
        n_Λ = GT.unit_normal(Λ,Ω)
        h_Λ = GT.face_diameter_field(Λ)
    else
        conformity = :default
    end

    if params[:dirichlet_method] === :strong
        V = GT.lagrange_space(Ω,interpolation_degree;conformity,dirichlet_boundary=Γd)
        uhd = GT.dirichlet_field(Float64,V)
        GT.interpolate_dirichlet!(u,uhd)
    else
        n_Γd = GT.unit_normal(Γd,Ω)
        h_Γd = GT.face_diameter_field(Γd)
        dΓd = GT.measure(Γd,integration_degree)
        V = GT.lagrange_space(Ω,interpolation_degree;conformity)
    end

    if params[:dirichlet_method] === :multipliers
        Q = GT.lagrange_space(Γd,interpolation_degree-1;conformity=:L2)
        VxQ = V × Q
    end

    function a(u,v)
        r = ∫( x->∇(u,x)⋅∇(v,x), dΩ)
        if params[:dirichlet_method] === :nitsche
            r += ∫( x->
                   (γ/h_Γd(x))*v(x)*u(x)-v(x)*n_Γd(x)⋅∇(u,x)-n_Γd(x)⋅∇(v,x)*u(x), dΓd)
        end
        if params[:discretization_method] === :interior_penalty
            r += ∫( x->
                   (γ/h_Λ(x))*jump(v,n_Λ,x)⋅jump(u,n_Λ,x)-jump(v,n_Λ,x)⋅mean(∇(u),x)-mean(∇(v),x)⋅jump(u,n_Λ,x), dΛ)
        end
        r
    end
    function l(v)
        r =∫( x->f(x)*v(x), dΩ) + ∫( x->g(x)*v(x), dΓn)
        if params[:dirichlet_method] === :nitsche
            r += ∫( x->
                (γ/h_Γd(x))*v(x)*u(x)-n_Γd(x)⋅∇(v,x)*u(x), dΓd)
        end
        r
    end

    if params[:dirichlet_method] === :multipliers
        function A((u,p),(v,q))
            r = a(u,v)
            r += ∫(x->(u(x)+p(x))*(v(x)+q(x))-u(x)*v(x)-p(x)*q(x), dΓd)
            r
        end
        function L((v,q))
            r = l(v)
            r += ∫(x->u(x)*q(x), dΓd)
            r
        end
    end

    @timeit timer "assembly" begin
        if params[:dirichlet_method] === :strong
            # TODO give a hint when dirichlet BCS are homogeneous
            x,K,b = GT.linear_problem(uhd,a,l)
        elseif params[:dirichlet_method] === :multipliers
            x,K,b = GT.linear_problem(Float64,VxQ,A,L)
        else
            x,K,b = GT.linear_problem(Float64,V,a,l)
        end
    end

    @timeit timer "solver" begin
        P = ps.setup(params[:solver],x,K,b)
        # TODO the solver should tell if the initial guess will be actually used
        fill!(x,0) # set initial guess
        ps.solve!(x,P,b)
        if params[:dirichlet_method] === :strong
            uh = GT.solution_field(uhd,x)
        elseif params[:dirichlet_method] === :multipliers
            uh_ph = GT.solution_field(VxQ,x)
            uh, ph = uh_ph
        else
            uh = GT.solution_field(V,x)
        end
    end

    @timeit timer "error_norms" begin
        eh(x) = u(x) - uh(x)
        ∇eh(x) = ∇(u,x) - ∇(uh,x)
        el2 = ∫( x->abs2(eh(x)), dΩ) |> sum |> sqrt
        eh1 = ∫( x->∇eh(x)⋅∇eh(x), dΩ) |> sum |> sqrt
        results[:error_h1_norm] = eh1
        results[:error_l2_norm] = el2
    end

    @timeit timer "vtk" if params[:export_vtu]
        GT.vtk_plot(params[:example_path]*"_Ω",Ω;refinement=4) do plt
            GT.plot!(plt,u;label="u")
            GT.plot!(plt,f;label="f")
            GT.plot!(plt,uh;label="uh")
        end
        GT.vtk_plot(params[:example_path]*"_Γn",Γn;refinement=4) do plt
            GT.plot!(plt,n_Γn;label="n")
            GT.plot!(plt,g;label="g")
        end
    end

    if params[:verbose]
        display(timer)
    end

    results
end

function main_hand_written(params)
    timer = TimerOutput()
    results = Dict{Symbol,Any}()
    error("Not implemented yet")
    results
end

function add_default_params(params_in,params_default)
    UmD = setdiff(keys(params_in),keys(params_default))
    for k in UmD
        @warn "Parameter with key :$k is unused"
    end
    UiD = intersect(keys(params_default),keys(params_in))
    DmU = setdiff(keys(params_default),keys(params_in))
    a = [ k=>params_in[k] for k in UiD]
    b = [ k=>params_default[k] for k in DmU]
    Dict(vcat(a,b))
end

function default_params()
    params = Dict{Symbol,Any}()
    params[:implementation] = :automatic
    params[:discretization_method] = :continuous_galerkin
    params[:verbose] = true
    params[:timer] = TimerOutput()
    params[:mesh] = GT.cartesian_mesh((0,1,0,1),(10,10))
    params[:u] = (x) -> sum(x)
    params[:domain_tags] = ["interior"]
    params[:dirichlet_tags] = ["boundary"]
    params[:neumann_tags] = String[]
    params[:dirichlet_method] = :strong
    params[:integration_degree] = 1
    params[:interpolation_degree] = 2*params[:integration_degree]
    params[:solver] = ps.lu_solver()
    params[:example_path] = joinpath(mkpath(joinpath(@__DIR__,"..","output")),"example_001")
    params[:export_vtu] = true
    params
end

end #module

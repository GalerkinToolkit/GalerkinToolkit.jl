"""
This implements a Poisson solver with several methods
"""
module Poisson

import GalerkinToolkit as gk
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

Δ(u) = gk.call(laplacian,u)
∇(u) = gk.call(gradient,u)
Δ(u,x) = Δ(u)(x)
∇(u,x) = ∇(u)(x)

mean(u,x) = 0.5*(u(x)[1]+u(x)[2])
jump(u,n,x) = u(x)[1]*n[1](x) + u(x)[2]*n[2](x)

function main_automatic(params)
    timer = params[:timer]
    results = Dict{Symbol,Any}()

    mesh = params[:mesh]
    D = gk.num_dims(mesh)
    Ω = gk.interior(mesh,physical_names=params[:domain_tags])
    Γd = gk.boundary(mesh;physical_names=params[:dirichlet_tags])
    Γn = gk.boundary(mesh;physical_names=params[:neumann_tags])
    integration_degree = params[:integration_degree]
    dΩ = gk.measure(Ω,integration_degree)
    dΓn = gk.measure(Γn,integration_degree)

    u = gk.analytical_field(params[:u],Ω)
    f(x) = -Δ(u,x)
    n_Γn = gk.unit_normal(Γn,Ω)
    g(x) = n_Γn(x)⋅∇(u,x)

    interpolation_degree = params[:interpolation_degree]

    @assert params[:discretization_method] in (:continuous_galerkin,:interior_penalty)

    if params[:discretization_method] !== :continuous_galerkin
        conformity = :L2
        gk.label_interior_faces!(mesh;physical_name="__INTERIOR_FACES__")
        Λ = gk.skeleton(mesh;physical_names=["__INTERIOR_FACES__"])
        dΛ = gk.measure(Λ,integration_degree)
        n_Λ = gk.unit_normal(Λ,Ω)
        h_Λ = gk.face_diameter_field(Λ)
    else
        conformity = :default
    end

    if params[:dirichlet_method] === :strong
        V = gk.lagrange_space(Ω,interpolation_degree;conformity,dirichlet_boundary=Γd)
        uh = gk.zero_field(Float64,V)
        gk.interpolate_dirichlet!(u,uh)
    else
        n_Γd = gk.unit_normal(Γd,Ω)
        h_Γd = gk.face_diameter_field(Γd)
        dΓd = gk.measure(Γd,integration_degree)
        γ = integration_degree*(integration_degree+1)
        γ = γ/10.0
        V = gk.lagrange_space(Ω,interpolation_degree;conformity)
        uh = gk.zero_field(Float64,V)
    end

    if params[:dirichlet_method] === :multipliers
        Q = gk.lagrange_space(Γd,interpolation_degree-1;conformity=:L2)
        VxQ = V × Q
        uh_qh = gk.zero_field(Float64,VxQ)
        uh, qh = uh_qh
    end

    function a(u,v)
        r = ∫( x->∇(u,x)⋅∇(v,x), dΩ)
        if params[:dirichlet_method] === :nitsche
            r += ∫( x->
                   (γ/h_Γd(x))*v(x)*u(x)-v(x)*n_Γd(x)⋅∇(u,x)-n_Γd(x)⋅∇(v,x)*u(x), dΓd)
        end
        if params[:discretization_method] === :interior_penalty
            r += ∫( x->
                   (γ/h_Γd(x))*jump(v,n_Λ,x)⋅jump(u,n_Λ,x)-jump(v,n_Λ,x)⋅mean(∇(u),x)-mean(∇(v),x)⋅jump(u,n_Λ,x), dΛ)
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
        Uh = uh_qh
    else
        A = a
        L = l
        Uh = uh
    end

    @timeit timer "assembly" begin
        # TODO give a hint when dirichlet BCS are homogeneous or not present
        x,A,b = gk.linear_problem(Uh,A,L)
    end

    @timeit timer "solver" begin
        P = ps.setup(params[:solver],x,A,b)
        ps.solve!(x,P,b)
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
        gk.vtk_plot(params[:example_path]*"_Ω",Ω;refinement=4) do plt
            gk.plot!(plt,u;label="u")
            gk.plot!(plt,f;label="f")
            gk.plot!(plt,uh;label="uh")
        end
        gk.vtk_plot(params[:example_path]*"_Γn",Γn;refinement=4) do plt
            gk.plot!(plt,n_Γn;label="n")
            gk.plot!(plt,g;label="g")
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
    params[:mesh] = gk.cartesian_mesh((0,1,0,1),(10,10))
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

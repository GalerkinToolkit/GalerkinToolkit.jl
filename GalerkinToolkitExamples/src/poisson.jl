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

Δ(u,x) = gk.call(laplacian,u)(x)
∇(u,x) = gk.call(gradient,u)(x)

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
    if params[:dirichlet_method] === :strong
        V = gk.lagrange_space(Ω,interpolation_degree;dirichlet_boundary=Γd)
        uh = gk.zero_field(Float64,V)
        gk.interpolate_dirichlet!(u,uh)
    else
        n_Γd = gk.unit_normal(Γd,Ω)
        h_Γd = gk.face_diameter_field(Γd)
        dΓd = gk.measure(Γd,integration_degree)
        γ = integration_degree*(integration_degree+1)
        γ = γ/10.0
        V = gk.lagrange_space(Ω,interpolation_degree)
        uh = gk.zero_field(Float64,V)
    end

    function a(u,v)
        r = ∫( x->∇(u,x)⋅∇(v,x), dΩ)
        if params[:dirichlet_method] === :nitsche
            r += ∫( x->
                   (γ/h_Γd(x))*v(x)*u(x)-v(x)*n_Γd(x)⋅∇(u,x)-n_Γd(x)⋅∇(v,x)*u(x), dΓd)
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

    @timeit timer "assembly" begin
        # TODO give a hint when dirichlet BCS are homogeneous or not present
        x,A,b = gk.linear_problem(uh,a,l)
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
    params[:verbose] = true
    params[:timer] = TimerOutput()
    params[:mesh] = gk.cartesian_mesh((0,1,0,1),(10,10))
    params[:u] = (x) -> sum(x)
    params[:domain_tags] = ["interior"]
    params[:dirichlet_tags] = ["boundary"]
    params[:neumann_tags] = String[]
    params[:dirichlet_method] = :strong
    params[:integration_degree] = 2
    params[:interpolation_degree] = 1
    params[:solver] = ps.lu_solver()
    params[:example_path] = joinpath(mkpath(joinpath(@__DIR__,"..","output")),"example_001")
    params[:export_vtu] = true
    params
end

end #module

"""
This implements a Poisson solver with several methods
"""
module Example001

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

function main_automatic(params)
    timer = params[:timer]
    results = Dict{Symbol,Any}()

    mesh = params[:mesh]
    D = gk.num_dims(mesh)
    Ω = gk.domain(mesh,physical_names=params[:domain_tags])
    Ωr = gk.reference_domain(Ω)
    Γd = gk.domain(mesh;face_dim=D-1,physical_names=params[:dirichlet_tags])
    Γn = gk.domain(mesh;face_dim=D-1,physical_names=params[:neumann_tags])
    Γnr = gk.reference_domain(Γn)
    ϕ_Ωr_Ω = gk.domain_map(Ωr,Ω)
    ϕ_Γnr_Γn = gk.domain_map(Γnr,Γn)
    ϕ_Γnr_Ωr = gk.domain_map(Γnr,Ωr;face_around=1)
    integration_degree = params[:integration_degree]
    dΩr = gk.measure(Ωr,integration_degree)
    dΓnr = gk.measure(Γnr,integration_degree)

    u = gk.analytical_field(params[:u],Ω)
    f(x) = -gk.call(Δ,u,x) # TODO
    nn = gk.unit_normal(Γnr,Ω;face_around=1)
    g(q) = nn(q)⋅ForwardDiff.gradient(u,ϕ_Γnr_Γn(q))

    function ∇(u,q)
        J = ForwardDiff.jacobian(ϕ_Ωr_Ω,q)
        grad = ForwardDiff.gradient(u,q)
        J\grad
    end
    function dV(q)
        J = ForwardDiff.jacobian(ϕ_Ωr_Ω,q)
        abs(det(J))
    end
    function dSn(q)
        J = ForwardDiff.jacobian(ϕ_Γnr_Γn,q)
        Jt = transpose(J)
        sqrt(det(Jt*J))
    end
    Δ(u,x) = tr(ForwardDiff.jacobian(y->ForwardDiff.gradient(u,y),x))

    if params[:dirichlet_method] == :nitsche
        Γdr = gk.reference_domain(Γd)
        ϕ_Γdr_Γd = gk.domain_map(Γdr,Γd)
        ϕ_Γdr_Ωr = gk.domain_map(Γdr,Ωr;face_around=1)
        nd = gk.unit_normal(Γdr,Ω;face_around=1)
        hd = gk.face_diameter_field(Γd)
        dΓdr = gk.measure(Γdr,integration_degree)
        γ = integration_degree*(integration_degree+1)
        γ = γ/10.0
        function dSd(q)
            J = ForwardDiff.jacobian(ϕ_Γdr_Γd,q)
            Jt = transpose(J)
            sqrt(det(Jt*J))
        end
    end

    interpolation_degree = params[:interpolation_degree]
    @assert params[:dirichlet_method] in (:strong,:nitsche)
    if params[:dirichlet_method] == :nitsche
        V = gk.lagrange_space(Ωr,interpolation_degree)
    else
        V = gk.lagrange_space(Ωr,interpolation_degree;dirichlet_boundary=Γd)
    end
    uh = gk.zero_field(Float64,V)
    gk.interpolate_dirichlet!(u∘ϕ_Ωr_Ω,uh)

    function a(u,v)
        r = ∫( q->∇(u,q)⋅∇(v,q)*dV(q), dΩr)
        if params[:dirichlet_method] == :nitsche
            r += ∫( q->
                   (γ/hd(q))*v(ϕ_Γdr_Ωr(q))*u(ϕ_Γdr_Ωr(q))*dSd(q)-
                   v(ϕ_Γdr_Ωr(q))*(nd(q)⋅∇(u,ϕ_Γdr_Ωr(q)))*dSd(q)-
                   (nd(q)⋅∇(v,ϕ_Γdr_Ωr(q)))*u(ϕ_Γdr_Ωr(q))*dSd(q), dΓdr)
        end
        r
    end
    function l(v)
        r =∫( q->f(ϕ_Ωr_Ω(q))*v(q)*dV(q) , dΩr) +
        ∫( q->g(q)*v(ϕ_Γnr_Ωr(q))*dSn(q) , dΓnr)
        if params[:dirichlet_method] == :nitsche
            r += ∫( q->
                (γ/hd(q))*v(ϕ_Γdr_Ωr(q))*u(ϕ_Γdr_Γd(q))*dSd(q)-
                (nd(q)⋅∇(v,ϕ_Γdr_Ωr(q)))*u(ϕ_Γdr_Γd(q))*dSd(q), dΓdr)
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
        eh(q) = u(ϕ_Ωr_Ω(q)) - uh(q)
        ∇eh(q) = ForwardDiff.gradient(u,ϕ_Ωr_Ω(q)) - ∇(uh,q)
        el2 = ∫( q->abs2(eh(q))*dV(q), dΩr) |> sum |> sqrt
        eh1 = ∫( q->∇eh(q)⋅∇eh(q)*dV(q), dΩr) |> sum |> sqrt
        results[:error_h1_norm] = eh1
        results[:error_l2_norm] = el2
    end

    @timeit timer "vtk" if params[:export_vtu]
        gk.vtk_plot(params[:example_path]*"_Ωr",Ωr;refinement=4) do plt
            gk.plot!(plt,q->u(ϕ_Ωr_Ω(q));label="u")
            gk.plot!(plt,q->f(ϕ_Ωr_Ω(q));label="f")
            gk.plot!(plt,uh;label="uh")
        end
        gk.vtk_plot(params[:example_path]*"_Γnr",Γnr;refinement=4) do plt
            gk.plot!(plt,nn;label="n")
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

"""
This implements a vanilla Poisson solver
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

    # Geometry
    @timeit timer "geometry" begin
        mesh = params[:mesh]
        D = gk.num_dims(mesh)
        Ω = gk.domain(mesh,physical_names=params[:domain_tags])
        Ωref = gk.reference_domain(Ω)
        Γdiri = gk.domain(mesh;face_dim=D-1,physical_names=params[:dirichlet_tags])
        Γ = gk.domain(mesh;face_dim=D-1,physical_names=params[:neumann_tags])
        Γref = gk.reference_domain(Γ)
    end

    # Geometrical map
    ϕ_Ωref_Ω = gk.domain_map(Ωref,Ω)
    ϕ_Γref_Γ = gk.domain_map(Γref,Γ)
    ϕ_Γref_Ωref = gk.domain_map(Γref,Ωref)
    function ∇(u,q)
        J = ForwardDiff.jacobian(ϕ_Ωref_Ω,q)
        g = ForwardDiff.gradient(u,q)
        J\g
    end
    function dV(q)
        J = ForwardDiff.jacobian(ϕ_Ωref_Ω,q)
        abs(det(J))
    end
    function dS(q)
        J = ForwardDiff.jacobian(ϕ_Γref_Γ,q)
        Jt = transpose(J)
        sqrt(det(Jt*J))
    end
    Δ(u,x) = tr(ForwardDiff.jacobian(y->ForwardDiff.gradient(u,y),x))

    # Manufactured functions
    u = gk.analytical_field(params[:u],Ω)
    f(x) = -gk.call(Δ,u,x) # TODO

    # Interpolation
    @timeit timer "interpolation" begin
        interpolation_degree = params[:interpolation_degree]
        V = gk.lagrange_space(Ωref,interpolation_degree;dirichlet_boundary=Γdiri)
        uh = gk.zero_field(Float64,V)
        gk.interpolate_dirichlet!(u∘ϕ_Ωref_Ω,uh)
    end

    # Integration
    integration_degree = params[:integration_degree]
    dΩref = gk.measure(Ωref,integration_degree)
    dΓref = gk.measure(Γref,integration_degree)

    # Weak form
    a(u,v) = ∫( q->∇(u,q)⋅∇(v,q)*dV(q), dΩref)
    l(v) =
        ∫( q->f(ϕ_Ωref_Ω(q))*v(q)*dV(q) , dΩref) #+
        #∫( q->g(ϕ_Γref_Γ(q))*v(ϕ_Γref_Ωref(q))*dS(q) , dΓref)

    # FE assembly
    @timeit timer "assembly" begin
        x,A,b = gk.linear_problem(uh,a,l)
    end

    # Linear solve
    @timeit timer "solver" begin
        P = ps.setup(params[:solver],x,A,b)
        ps.solve!(x,P,b)
    end

    # Error norms
    @timeit timer "error_norms" begin
        eh(q) = u(ϕ_Ωref_Ω(q)) - uh(q)
        ∇eh(q) = ForwardDiff.gradient(u,ϕ_Ωref_Ω(q)) - ∇(uh,q)
        el2 = ∫( q->abs2(eh(q))*dV(q), dΩref) |> sum |> sqrt
        eh1 = ∫( q->∇eh(q)⋅∇eh(q)*dV(q), dΩref) |> sum |> sqrt
        results[:error_h1_norm] = eh1
        results[:error_l2_norm] = el2
    end

    # Vtk output
    @timeit timer "vtk" if ! params[:export_vtu]
        gk.vtk_plot(params[:example_path],Ωref;refinement=4) do plt
            gk.plot!(plt,q->u(ϕ_Ωref_Ω(q));label="u")
            gk.plot!(plt,uh;label="u")
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
    params[:integration_degree] = 2
    params[:interpolation_degree] = 1
    params[:solver] = ps.lu_solver()
    params[:example_path] = joinpath(mkpath(joinpath(@__DIR__,"..","output")),"example_001")
    params[:export_vtu] = true
    params
end

end #module

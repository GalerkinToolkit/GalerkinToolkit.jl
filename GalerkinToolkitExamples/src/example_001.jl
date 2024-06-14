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
using IterativeSolvers: cg!

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
    end

    # Geometrical map
    ϕ = gk.domain_map(Ωref,Ω)
    function ∇(u,q)
        J = ForwardDiff.jacobian(ϕ,q)
        g = ForwardDiff.gradient(u,q)
        J\g
    end
    function dV(q)
        J = ForwardDiff.jacobian(ϕ,q)
        abs(det(J))
    end

    # Manufactured functions
    u = gk.analytical_field(params[:u],Ω)
    f = gk.analytical_field(params[:f],Ω)

    # Interpolation
    @timeit timer "interpolation" begin
        interpolation_degree = params[:interpolation_degree]
        V = gk.lagrange_space(Ωref,interpolation_degree;dirichlet_boundary=Γdiri)
        uh = gk.zero_field(Float64,V)
        gk.interpolate_dirichlet!(u∘ϕ,uh)
    end

    # Integration
    integration_degree = params[:integration_degree]
    dΩref = gk.measure(Ωref,integration_degree)


    # Weak form
    a(u,v) = ∫( q->∇(u,q)⋅∇(v,q)*dV(q), dΩref)
    l(v) = ∫( q->f(q)*dV(q), dΩref)

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
        eh(q) = u(ϕ(q)) - uh(q)
        ∇eh(q) = ForwardDiff.gradient(u,ϕ(q)) - ∇(uh,q)
        el2 = ∫( q->abs2(eh(q))*dV(q), dΩref) |> sum |> sqrt
        eh1 = ∫( q->∇eh(q)⋅∇eh(q)*dV(q), dΩref) |> sum |> sqrt
        results[:error_h1_norm] = eh1
        results[:error_l2_norm] = el2
    end

    # Vtk output
    @timeit timer "vtk" if ! params[:export_vtu]
        gk.vtk_plot(params[:example_path],Ωref;refinement=4) do plt
            gk.plot!(plt,q->u(ϕ(q));label="u")
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
    params[:f] = (x) -> 0.0
    params[:g] = (x) -> 0.0
    params[:domain_tags] = ["interior"]
    params[:dirichlet_tags] = ["boundary"]
    params[:neumann_tags] = String[]
    params[:integration_degree] = 2
    params[:interpolation_degree] = 1
    params[:solver] = ps.lu_solver()
    params[:example_path] = joinpath(mkpath("GalerkinToolkitExamples_output"),"example_001")
    params[:export_vtu] = true
    params
end

end #module

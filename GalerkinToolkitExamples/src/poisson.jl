"""
This implements a Poisson solver with several methods
"""
module Poisson

import GalerkinToolkit as GT
using GalerkinToolkit: ∫
import PartitionedSolvers as PS
import ForwardDiff
using StaticArrays
using LinearAlgebra
using SparseArrays
using PartitionedArrays
using TimerOutputs
using WriteVTK
using AbstractTrees


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

laplacian(u,x) = tr(ForwardDiff.jacobian(y->ForwardDiff.gradient(u,y),x))
# TODO hide this
laplacian(u::GT.AbstractQuantity,x::GT.AbstractQuantity) = GT.call(laplacian,u,x)

const ∇ = ForwardDiff.gradient
const Δ = laplacian

mean(f,u,x) = 0.5*(f(u[1],x)+f(u[2],x))
jump(u,n,x) = u[2](x)*n[2](x) + u[1](x)*n[1](x)

function main_automatic(params)
    timer = params[:timer]
    results = Dict{Symbol,Any}()

    mesh = params[:mesh]
    D = GT.num_dims(mesh)
    Ω = GT.interior(mesh,group_names=params[:domain_tags])
    Γd = GT.boundary(mesh;group_names=params[:dirichlet_tags])
    Γn = GT.boundary(mesh;group_names=params[:neumann_tags])
    integration_degree = params[:integration_degree]
    dΩ = GT.measure(Ω,integration_degree)
    dΓn = GT.measure(Γn,integration_degree)

    u = GT.analytical_field(params[:u],Ω)
    f(x) = -Δ(u,x)
    n_Γn = GT.unit_normal(mesh,D-1)
    g(x) = n_Γn(x)⋅∇(u,x)

    interpolation_degree = params[:interpolation_degree]
    γ = GT.uniform_quantity(integration_degree*(integration_degree+1)/10.0)

    @assert params[:discretization_method] in (:continuous_galerkin,:interior_penalty)

    if params[:discretization_method] !== :continuous_galerkin
        conformity = :L2
        GT.group_interior_faces!(mesh;group_name="__INTERIOR_FACES__")
        Λ = GT.skeleton(mesh;group_names=["__INTERIOR_FACES__"])
        dΛ = GT.measure(Λ,integration_degree)
        n_Λ = GT.unit_normal(mesh,D-1)
        h_Λ = GT.face_diameter_field(Λ)
    else
        conformity = :default
    end

    if params[:dirichlet_method] === :strong
        V = GT.lagrange_space(Ω,interpolation_degree;conformity,dirichlet_boundary=Γd)
        uhd = GT.zero_dirichlet_field(Float64,V)
        GT.interpolate_dirichlet!(u,uhd)
    else
        n_Γd = GT.unit_normal(mesh,D-1)
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
                   (γ/h_Λ(x))*jump(v,n_Λ,x)⋅jump(u,n_Λ,x)-jump(v,n_Λ,x)⋅mean(∇,u,x)-mean(∇,v,x)⋅jump(u,n_Λ,x), dΛ)
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
            p = GT.PartitionedSolvers_linear_problem(uhd,a,l)
        elseif params[:dirichlet_method] === :multipliers
            p = GT.PartitionedSolvers_linear_problem(Float64,VxQ,A,L)
        else
            p = GT.PartitionedSolvers_linear_problem(Float64,V,a,l)
        end
    end

    @timeit timer "solver" begin
        s = params[:solver](p)
        # TODO the solver should tell if the initial guess will be actually used
        x = PS.solution(p)
        fill!(x,0) # set initial guess
        s = PS.update(s,solution=x)
        s = PS.solve(s)
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
        vtk_grid(params[:example_path]*"_Ω",Ω;plot_params=(;refinement=4)) do plt
            GT.plot!(plt,u;label="u")
            GT.plot!(plt,f;label="f")
            GT.plot!(plt,uh;label="uh")
        end
        vtk_grid(params[:example_path]*"_Γn",Γn;plot_params=(;refinement=4)) do plt
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
    timer = params[:timer]
    results = Dict{Symbol,Any}()

    mesh = params[:mesh]
    D = GT.num_dims(mesh)
    Ω = GT.interior(mesh,group_names=params[:domain_tags])
    Γd = GT.boundary(mesh;group_names=params[:dirichlet_tags])
    Γn = GT.boundary(mesh;group_names=params[:neumann_tags])
    integration_degree = params[:integration_degree]
    dΩ = GT.measure(Ω,integration_degree)
    dΓn = GT.measure(Γn,integration_degree)

    u = GT.analytical_field(params[:u],Ω)
    f(x) = -Δ(u,x)
    n_Γn = GT.unit_normal(mesh,D-1)
    g(x) = n_Γn(x)⋅∇(u,x)

    interpolation_degree = params[:interpolation_degree]
    V = GT.lagrange_space(Ω,interpolation_degree;dirichlet_boundary=Γd)
    uhd = GT.zero_dirichlet_field(Float64,V)
    GT.interpolate_dirichlet!(u,uhd)
    T = Float64

    @timeit timer "assembly" begin
        A_alloc = GT.allocate_matrix(T,V,V,Ω)
        free_or_dirichlet=(GT.FREE,GT.DIRICHLET)
        Ad_alloc = GT.allocate_matrix(T,V,V,Ω;free_or_dirichlet)
        assemble_matrix!(A_alloc,Ad_alloc,V,dΩ)
        A = GT.compress(A_alloc)
        Ad = GT.compress(Ad_alloc)
        xd = GT.dirichlet_values(uhd)
        b = -Ad*xd
        x = similar(b)
        fill!(x,0) # set initial guess
        p = PS.linear_problem(x,A,b)
    end

    @timeit timer "solver" begin
        s = params[:solver](p)
        s = PS.solve(s)
        uh = GT.solution_field(uhd,x)
    end

    @timeit timer "error_norms" begin
        el2,eh1 = integrate_error(u,uh,dΩ)
        results[:error_h1_norm] = eh1
        results[:error_l2_norm] = el2
    end

    @timeit timer "vtk" if params[:export_vtu]
        vtk_grid(params[:example_path]*"_Ω",Ω;plot_params=(;refinement=4)) do plt
            GT.plot!(plt,u;label="u")
            GT.plot!(plt,f;label="f")
            GT.plot!(plt,uh;label="uh")
        end
    end

    if params[:verbose]
        display(timer)
    end
    results
end


function assemble_matrix!(A_alloc,Ad_alloc,V,dΩ)

    #Accessors to the quantities on the
    #integration points
    Ω = GT.domain(dΩ)
    face_point_J = GT.jacobian_accessor(dΩ)
    face_point_dV = GT.weight_accessor(dΩ)
    face_npoints = GT.num_points_accessor(dΩ)
    face_dofs = GT.dofs_accessor(V,Ω)
    ∇ = ForwardDiff.gradient
    face_point_dof_∇s = GT.shape_function_accessor(∇,V,dΩ)

    #Temporaries
    n = GT.max_num_reference_dofs(V)
    T = Float64
    Auu = zeros(T,n,n)

    #Numerical integration loop
    for face in 1:GT.num_faces(Ω)

        #Get quantities at current face
        npoints = face_npoints(face)
        point_J = face_point_J(face)
        point_dV = face_point_dV(face)
        point_dof_∇s = face_point_dof_∇s(face)
        dofs = face_dofs(face)

        #Reset face matrix
        fill!(Auu,zero(T))

        #Loop over integration points
        for point in 1:npoints

            #Get quantities at current integration point
            J = point_J(point)
            dV = point_dV(point,J)
            dof_∇s = point_dof_∇s(point,J)

            #Fill in face matrix
            for (i,dofi) in enumerate(dofs)
                ∇v = dof_∇s(i)
                for (j,dofj) in enumerate(dofs)
                    ∇u = dof_∇s(j)
                    Auu[i,j] += ∇v⋅∇u*dV
                end
            end
        end

        #Add face contribution to the
        #global allocations
        GT.contribute!(A_alloc,Auu,dofs,dofs)
        GT.contribute!(Ad_alloc,Auu,dofs,dofs)
    end
end

function norm2(a)
    a⋅a
end

function integrate_error(g,uh,dΩ)

    #Accessors to the quantities on the
    #integration points
    face_point_x = GT.coordinate_accessor(dΩ)
    face_point_J = GT.jacobian_accessor(dΩ)
    face_point_dV = GT.weight_accessor(dΩ)
    face_npoints = GT.num_points_accessor(dΩ)
    face_point_uhx = GT.discrete_field_accessor(GT.value,uh,dΩ)
    face_point_∇uhx = GT.discrete_field_accessor(ForwardDiff.gradient,uh,dΩ)

    #Numerical integration loop
    s = 0.0
    s1 = 0.0
    Ω = GT.domain(dΩ)
    for face in 1:GT.num_faces(Ω)

        #Get quantities at current face
        npoints = face_npoints(face)
        point_x = face_point_x(face)
        point_J = face_point_J(face)
        point_dV = face_point_dV(face)
        point_uhx = face_point_uhx(face)
        point_∇uhx = face_point_∇uhx(face)

        for point in 1:npoints

            #Get quantities at current integration point
            x = point_x(point)
            J = point_J(point)
            dV = point_dV(point,J)
            uhx = point_uhx(point,J)
            ∇uhx = point_∇uhx(point,J)

            #Add contribution
            s += abs2(uhx-g.definition(x))*dV
            s1 += norm2(∇uhx-ForwardDiff.gradient(g.definition,x))
        end
    end

    #Compute the final norms
    el2 = sqrt(s)
    eh1 = sqrt(s1)

    (el2,eh1)
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
    params[:u] = sum
    params[:domain_tags] = ["interior"]
    params[:dirichlet_tags] = ["boundary"]
    params[:neumann_tags] = String[]
    params[:dirichlet_method] = :strong
    params[:interpolation_degree] = 1
    params[:integration_degree] = 2*params[:interpolation_degree]
    params[:solver] = PS.LinearAlgebra_lu
    params[:example_path] = joinpath(mkpath(joinpath(@__DIR__,"..","output")),"example_001")
    params[:export_vtu] = true
    params
end

end #module

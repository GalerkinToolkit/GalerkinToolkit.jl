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
laplacian(u,x) = tr(ForwardDiff.jacobian(y->ForwardDiff.gradient(u,y),x))
laplacian(u) = x->laplacian(u,x)

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
            p = GT.linear_problem(uhd,a,l)
        elseif params[:dirichlet_method] === :multipliers
            p = GT.linear_problem(Float64,VxQ,A,L)
        else
            p = GT.linear_problem(Float64,V,a,l)
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
    state0 = (;results,params)
    @timeit timer "setup" state1 = setup(state0)
    @timeit timer "assemble_system" state2 = assemble_system(state1)
    @timeit timer "solve_system" state3 = solve_system(state2)
    @timeit timer "integrate_error_norms" integrate_error_norms(state3)
    @timeit timer "export_results" export_results(state3)
    if params[:verbose]
        display(timer)
    end
    results
end

function setup(state0)
    state1 = setup_domains(state0)
    state2 = setup_integration(state1)
    state3 = setup_interpolation(state2)
    state4 = setup_dirichlet(state3)
    return state4
end

function assemble_system(state0)
    state1 = assemble_allocate(state0)
    state2 = assemble_fill(state1)
    state3 = assemble_compress(state2)
    return state3
end

function export_results(state)
    (;params,domains,interpolation) = state
    (;x_diri,x_free,face_to_rid_V,face_to_dofs_V) = interpolation.workspace
    (;Ω) = domains
    (;V) = interpolation.spaces

    T = Float64
    z = zero(T)

    rid_to_reffe = GT.reference_fes(V)
    plt = GT.plot(Ω;refinement=4)
    vis_mesh, vis_glue = GT.visualization_mesh(plt,glue=Val(true))
    n_vis_nodes = GT.num_nodes(vis_mesh)
    vis_node_to_uh = zeros(T,n_vis_nodes)
    rid_to_ref_vis_xs = vis_glue.reference_coordinates
    rid_to_vis_cache = map(rid_to_ref_vis_xs,rid_to_reffe) do vis_xs,reffe
        ref_vals = GT.tabulator(reffe)(GT.value,vis_xs)
        npoints,ndofs = size(ref_vals)
        ue = zeros(T,ndofs)
        (;ref_vals,ndofs,ue,npoints)
    end

    face_to_cell = GT.faces(Ω)
    nfaces = length(face_to_cell)

    vis_node = 0
    for face in 1:nfaces
        rid = face_to_rid_V[face]
        cache = rid_to_vis_cache[rid]
        dofs = face_to_dofs_V[face]
        ndofs = cache.ndofs
        # Fill values
        for i in 1:ndofs
            dof = dofs[i]
            if dof < 0
                cache.ue[i] = x_diri[-dof]
            else
                cache.ue[i] = x_free[dof]
            end
        end
        # Compute FE solution for each visualization point
        for point in 1:cache.npoints
            uhx = z
            for i in 1:ndofs
                ui = cache.ue[i]
                uhx += cache.ref_vals[point,i]*ui
            end
            vis_node += 1
            vis_node_to_uh[vis_node] = uhx
        end
    end
    # Set visualization data
    GT.node_data(plt)["uh"] = vis_node_to_uh
    # Write it to vtk
    vtk_grid(params[:example_path]*"_Ω",plt) |> WriteVTK.close
end

function integrate_error_norms(state)
    (;results,params,domains,interpolation,integration) = state
    (;x_diri,x_free,rid_to_cache_V,face_to_rid_V,face_to_dofs_V) = interpolation.workspace
    (;rid_to_cache_dΩ,nodes_to_coords) = integration.workspace
    (;Ω) = domains
    mesh = params[:mesh]
    u = params[:u]

    T = Float64
    el2 = zero(T)
    eh1 = zero(T)
    z = zero(T)

    face_to_cell = GT.faces(Ω)
    nfaces = length(face_to_cell)

    D = GT.num_dims(mesh)
    cell_to_nodes = GT.face_nodes(mesh,D)

    Tgrad = SVector{D,T}
    Tx = Tgrad
    TJ = typeof(zero(Tx)*transpose(zero(Tx)))

    for face in 1:nfaces
        rid = face_to_rid_V[face]
        cache_V = rid_to_cache_V[rid]
        cache_dΩ = rid_to_cache_dΩ[rid]
        dofs = face_to_dofs_V[face]
        cell = face_to_cell[face]
        nodes = cell_to_nodes[cell]
        ndofs = cache_V.ndofs
        # Fill values
        for i in 1:ndofs
            dof = dofs[i]
            if dof < 0
                cache_V.u[i] = x_diri[-dof]
            else
                cache_V.u[i] = x_free[dof]
            end
        end
        # Integration loop
        for point in 1:cache_dΩ.npoints
            # Compute integration coordinate and Jacobian transpose
            x = zero(Tx)
            Jt = zero(TJ)
            for (inode,node) in enumerate(nodes)
                coords = nodes_to_coords[node]
                x += cache_dΩ.ref_vals[point,inode]*coords
                Jt += coords*transpose(cache_dΩ.ref_grads[point,inode])
            end
            # Integration weight
            dV = abs(det(Jt))*cache_dΩ.weights[point]
            # Compute physical gradients
            invJt = inv(Jt)
            for i in 1:ndofs
                cache_V.phys_grads[i] = invJt*cache_V.ref_grads[point,i]
            end
            # Compute FE solution and its gradients
            uhx = z
            ∇uhx = zero(Tx)
            for i in 1:ndofs
                ui = cache_V.u[i]
                uhx += cache_V.ref_vals[point,i]*ui
                ∇uhx += cache_V.phys_grads[i]*ui
            end
            # Compute and integrate the error
            ex = u(x)-uhx
            ∇ex = ForwardDiff.gradient(u,x)-∇uhx
            el2 += ex*ex*dV
            eh1 += ∇ex⋅∇ex*dV
        end
    end

    el2 = sqrt(el2)
    eh1 = sqrt(eh1)

    results[:error_h1_norm] = eh1
    results[:error_l2_norm] = el2
end

function solve_system(state)
    (;params,linear_system,interpolation) = state
    (;matrix,rhs) = linear_system
    (;uh) = interpolation
    (;x_free) = interpolation.workspace
    p = PS.linear_problem(x_free,matrix,rhs)
    s = params[:solver](p)
    s = PS.solve(s)
    x_free = PS.solution(s)
    uh = GT.solution_field(uh,x_free)
    (;interpolation,linear_system,state...)
end

function assemble_compress(state)
    (;assemble_strategy) = state
    (;matrix_strategy,matrix_strategy_init,matrix_strategy_alloc) = assemble_strategy.matrix
    (;vector_strategy,vector_strategy_init,vector_strategy_alloc) = assemble_strategy.vector
    rhs, rhs_cache = vector_strategy.compress(vector_strategy_alloc,vector_strategy_init)
    matrix, matrix_cache = matrix_strategy.compress(matrix_strategy_alloc,matrix_strategy_init)
    linear_system = (;matrix,rhs,matrix_cache,rhs_cache)
    (;linear_system,state...)
end

function assemble_fill(state)
    (;params,domains,interpolation,integration,assemble_strategy) = state
    (;matrix_strategy,matrix_strategy_init,matrix_strategy_counter,matrix_strategy_alloc) = assemble_strategy.matrix
    (;vector_strategy,vector_strategy_init,vector_strategy_counter,vector_strategy_alloc) = assemble_strategy.vector
    (;Ω) = domains
    (;x_diri,rid_to_cache_V,face_to_rid_V,face_to_dofs_V) = interpolation.workspace
    (;rid_to_cache_dΩ,nodes_to_coords) = integration.workspace
    mesh = params[:mesh]
    u = params[:u]
    f(x) = -laplacian(u,x)
    T = Float64

    matrix_strategy_counter = matrix_strategy.counter(matrix_strategy_init)
    vector_strategy_counter = vector_strategy.counter(vector_strategy_init)
    z = zero(T)
    field_i = field_j = 1

    face_to_cell = GT.faces(Ω)
    nfaces = length(face_to_cell)

    D = GT.num_dims(mesh)
    cell_to_nodes = GT.face_nodes(mesh,D)

    Tgrad = SVector{D,T}
    Tx = Tgrad
    TJ = typeof(zero(Tx)*transpose(zero(Tx)))

    for face in 1:nfaces
        rid = face_to_rid_V[face]
        cache_V = rid_to_cache_V[rid]
        cache_dΩ = rid_to_cache_dΩ[rid]
        
        # Face degrees of freedom
        dofs = face_to_dofs_V[face]

        # Set face index
        cell = face_to_cell[face]
        
        # Definition of the face in terms of node ids
        nodes = cell_to_nodes[cell]

        fill!(cache_V.A,z)
        fill!(cache_V.b,z)
        fill!(cache_V.u,z)

        ndofs = cache_V.ndofs
        
        # Fill in Dirichlet values
        # Use face dofs to assign the initial values of uh at the boundary
        for i in 1:ndofs
            dof = dofs[i]
            if dof < 0
                cache_V.u[i] = x_diri[-dof]
            end
        end

        # Integration loop
        for point in 1:cache_dΩ.npoints  #-> loop over quadrature points
            # Compute integration coordinate and Jacobian transpose
            x = zero(Tx)
            Jt = zero(TJ)
            for (inode,node) in enumerate(nodes)
                coords = nodes_to_coords[node]
                x += cache_dΩ.ref_vals[point,inode]*coords
                Jt += coords*transpose(cache_dΩ.ref_grads[point,inode])
            end
            # Integration weight
            dV = abs(det(Jt))*cache_dΩ.weights[point]
            # Compute physical gradients
            invJt = inv(Jt)
            for i in 1:ndofs
                cache_V.phys_grads[i] = invJt*cache_V.ref_grads[point,i]
            end
            # Integrate matrix and vector
            fx = f(x)
            for i in 1:ndofs
                cache_V.b[i] += cache_V.ref_vals[i]*fx*dV
                for j in 1:ndofs
                    Aij = cache_V.phys_grads[i]⋅cache_V.phys_grads[j]*dV
                    cache_V.A[i,j] += Aij
                    cache_V.b[i] -= Aij*cache_V.u[j]
                end
            end
        end
        # Assemble the matrix and vector
        for (i,dof_i) in enumerate(dofs)
            vector_strategy_counter =
            vector_strategy.set!(vector_strategy_alloc,
                                 vector_strategy_counter,
                                 vector_strategy_init,
                                 cache_V.b[i],
                                 dof_i,field_i)
            for (j,dof_j) in enumerate(dofs)
                matrix_strategy_counter =
                matrix_strategy.set!(matrix_strategy_alloc,
                                     matrix_strategy_counter,
                                     matrix_strategy_init,
                                     cache_V.A[i,j],
                                     dof_i,dof_j,field_i,field_j)
            end
        end
    end
    (;assemble_strategy,state...)
end

function assemble_allocate(state)
    (;domains,interpolation) = state
    (;V) = interpolation.spaces
    (;face_to_dofs_V) = interpolation.workspace
    (;Ω) = domains
    free_dofs = GT.free_dofs(V)
    T = Float64

    matrix_strategy =  GT.monolithic_matrix_assembly_strategy()
    vector_strategy = GT.monolithic_vector_assembly_strategy()

    matrix_strategy_init = matrix_strategy.init(free_dofs,free_dofs,T)
    vector_strategy_init = vector_strategy.init(free_dofs,T)
    field_i = field_j = 1

    matrix_strategy_counter = matrix_strategy.counter(matrix_strategy_init)
    vector_strategy_counter = vector_strategy.counter(vector_strategy_init)

    face_to_cell = GT.faces(Ω)
    nfaces = length(face_to_cell)

    for face in 1:nfaces
        dofs = face_to_dofs_V[face]
        for dof_i in dofs
            vector_strategy_counter =
            vector_strategy.count(vector_strategy_counter,
                                  vector_strategy_init,
                                  dof_i,field_i)
            for dof_j in dofs
                matrix_strategy_counter =
                matrix_strategy.count(
                                      matrix_strategy_counter,
                                      matrix_strategy_init,
                                      dof_i,dof_j,field_i,field_j)
            end
        end
    end

    matrix_strategy_alloc = matrix_strategy.allocate(matrix_strategy_counter,matrix_strategy_init)
    vector_strategy_alloc = vector_strategy.allocate(vector_strategy_counter,vector_strategy_init)

    matrix = (;matrix_strategy,matrix_strategy_init,matrix_strategy_counter,matrix_strategy_alloc)
    vector = (;vector_strategy,vector_strategy_init,vector_strategy_counter,vector_strategy_alloc)

    assemble_strategy = (;matrix,vector)
    (;assemble_strategy,state...)
end

function setup_domains(state)
    (;params,) = state
    mesh = params[:mesh]
    Ω = GT.interior(mesh,physical_names=params[:domain_tags])
    Γd = GT.boundary(mesh;physical_names=params[:dirichlet_tags])
    domains = (;Ω,Γd)
    (;domains,state...)
end

function setup_integration(state)
    (;params,domains) = state
    (;Ω) = domains
    integration_degree = params[:integration_degree]
    mesh = params[:mesh]
    dΩ = GT.measure(Ω,integration_degree)
    D = GT.num_dims(mesh)
    rid_to_refface = GT.reference_faces(mesh,D)
    rid_to_refquad = GT.reference_quadratures(dΩ)
    rid_to_cache_dΩ = map(rid_to_refface,rid_to_refquad) do refface,refquad
            x = GT.coordinates(refquad)
            npoints = length(x)
            weights = GT.weights(refquad)
            ref_vals = GT.tabulator(refface)(GT.value,x)
            ref_grads = GT.tabulator(refface)(ForwardDiff.gradient,x)
            (;x,weights,ref_vals,ref_grads,npoints)
    end
    face_to_rid_dΩ = nothing #face_reference_id(dΩ)
    Dface_to_nodes = GT.face_nodes(mesh,D)
    nodes_to_coords = GT.node_coordinates(mesh)
    workspace = (;rid_to_cache_dΩ,face_to_rid_dΩ,Dface_to_nodes,nodes_to_coords)
    measures = (;dΩ,)
    integration = (;measures,workspace)
    (;integration,state...)
end

function setup_interpolation(state)
    (;params,domains,integration) = state
    (;dΩ) = integration.measures
    (;Ω,Γd) = domains
    interpolation_degree = params[:interpolation_degree]
    mesh = params[:mesh]
    V = GT.lagrange_space(Ω,interpolation_degree;dirichlet_boundary=Γd)
    rid_to_reffe = GT.reference_fes(V)
    rid_to_refquad = GT.reference_quadratures(dΩ)
    T = Float64
    D = GT.num_dims(mesh)
    Tgrad = SVector{D,T}
    rid_to_cache_V = map(rid_to_reffe,rid_to_refquad) do reffe,refquad
            x = GT.coordinates(refquad)
            ndofs = GT.num_dofs(reffe)
            A = zeros(T,ndofs,ndofs)
            b = zeros(T,ndofs)
            u = zeros(T,ndofs)
            ref_vals = GT.tabulator(reffe)(GT.value,x)
            ref_grads = GT.tabulator(reffe)(ForwardDiff.gradient,x)
            phys_grads = zeros(Tgrad,ndofs)
            (;A,b,u,ref_vals,ref_grads,ndofs,phys_grads)
    end
    face_to_rid_V = GT.face_reference_id(V)
    uh = GT.zero_field(Float64,V)
    x_free = GT.free_values(uh)
    x_diri = GT.dirichlet_values(uh)
    face_to_dofs_V = GT.face_dofs(V)
    interpolation = (;uh,spaces=(;V),workspace=(;rid_to_cache_V,face_to_rid_V,face_to_dofs_V,x_free,x_diri))
    (;interpolation,state...)
end

function setup_dirichlet(state)
    (;params,interpolation) = state
    (;V) = interpolation.spaces
    (;x_diri) = interpolation.workspace
    u = params[:u]
    vnode_to_coord = GT.node_coordinates(V)
    diri_dof_to_vnode = GT.dirichlet_dof_node(V)
    x_diri .= u.(vnode_to_coord[diri_dof_to_vnode])
    (;interpolation,state...)
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
    params[:integration_degree] = 1
    params[:interpolation_degree] = 2*params[:integration_degree]
    params[:solver] = PS.LinearAlgebra_lu
    params[:example_path] = joinpath(mkpath(joinpath(@__DIR__,"..","output")),"example_001")
    params[:export_vtu] = true
    params
end

end #module

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

function main_hand_written_new(params)
    results = Dict{Symbol,Any}()
    state0 = (;results,params)
    @timeit timer "setup" state1 = setup(state0)
    @timeit timer "assemble_system" state2 = assemble_system(state1)
    @timeit timer "solve_system" state3 = solve_system(state2)
    @timeit timer "integrate_error_norms" integrate_error_norms(state3)
    @timeit timer "export_results" export_results(state3)
    display(timer)
    results
end

function setup(state0)
    state1 = setup_domains(state0)
    state2 = setup_integration(state1)
    state3 = setup_interpolation(state3)
    state4 = setup_dirichlet(state3)
end

function assemble_system(state0)
    state1 = assemble_allocate(state)
    state2 = assemble_fill(state1)
    state3 = assemble_compress(state2)
end

function solve_system(state0)
end

function integrate_error_norms(state0)
end

function export_results(state0)
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
    integration_degree = params[:integration_degree]
    dΩ = GT.measure(Ω,integration_degree)
    rid_to_cache_dΩ = map(rid_to_refface,rid_to_refquad) do refface,refquad
            x = GT.coordinates(refquad)
            weights = GT.weights(refquad)
            ref_vals = GT.tabulator(refface)(GT.value,x)
            ref_grads = GT.tabulator(refface)(ForwardDiff.gradient,x)
    end
    face_to_rid_dΩ = face_reference_id(dΩ)
    Dface_to_nodes = face_nodes(mesh,D)
    nodes_to_coords = node_coordinates(mesh)
    workspace = (;rid_to_cache_dΩ,face_to_rid_dΩ,Dface_to_nodes,node_to_coords)
    measures = (;dΩ,)
    integration = (;measures,workspace)
    (;integration,state...)
end

function setup_interpolation(state)
    (;domains,integration) = state
    (;Ω,Γd) = domains
    interpolation_degree = params[:interpolation_degree]
    V = GT.lagrange_space(Ω,interpolation_degree;dirichlet_boundary=Γd)
    rid_to_cache_dΩ = integration.workspace.rid_to_cache_dΩ
    rid_to_cache_V = map(rid_to_reffe,rid_to_cache_dΩ) do reffe,cache
            ndofs = GT.num_dofs(reffe)
            A = zeros(T,ndofs,ndofs)
            b = zeros(T,ndofs)
            u = zeros(T,ndofs)
            ref_vals = GT.tabulator(reffe)(GT.value,x)
            ref_grads = GT.tabulator(reffe)(ForwardDiff.gradient,x)
    end
    face_to_rid_V = face_reference_id(V)
    uh = GT.zero_field(Float64,V)
    x_free = GT.free_values(uh)
    x_diri = GT.dirichlet_values(uh)
    face_to_dofs_V = face_dofs(V)
    interpolation = (;uh,spaces=(;V),workspace=(;rid_to_cache_V,face_to_rid_V,face_to_dofs_V,x_free,x_diri))
    (;interpolation,state...)
end

function setup_dirichlet(state)
    vnode_to_coord = GT.node_coordinates(V)
    diri_dof_to_vnode = GT.dirichlet_dof_node(V)
    x_diri .= u.(vnode_to_coord[diri_dof_to_vnode])
end


function main_hand_written(params)

    @assert params[:discretization_method] == :continuous_galerkin
    @assert params[:dirichlet_method] === :strong

    timer = params[:timer]
    results = Dict{Symbol,Any}()

    mesh = params[:mesh]
    D = GT.num_dims(mesh)

    # Set computational domains
    Ω = GT.interior(mesh,physical_names=params[:domain_tags])
    Γd = GT.boundary(mesh;physical_names=params[:dirichlet_tags])
    
    interpolation_degree = params[:interpolation_degree]
    integration_degree = params[:integration_degree]

    # Set up numerical quadratures
    dΩ = GT.measure(Ω,integration_degree)

    u = params[:u]
    f(x) = -laplacian(u,x)

    # Set interpolation space -> Lagrange polynomial basis function, the element type and geometry
    V = GT.lagrange_space(Ω,interpolation_degree;dirichlet_boundary=Γd)

    # Initializing DiscreteField with zero interior and boundary values
    uh = GT.zero_field(Float64,V)

    # Extract uh values associated with out-of-the-boundary nodes
    x_free = GT.free_values(uh)
    # Extract uh values associated with nodes at boundary
    x_diri = GT.dirichlet_values(uh)

    # Calculate faces degrees of freedom (dofs) -> negative dofs are associated with nodes at boundary
    face_to_dofs = GT.face_dofs(V)

    # Get node coordinates in the computational FE space
    vnode_to_coord = GT.node_coordinates(V)

    # Calculate degree of freedom for each boundary node in the FE space
    diri_dof_to_vnode = GT.dirichlet_dof_node(V)

    # Calculate the value of the guess u function at each boundary node
    x_diri .= u.(vnode_to_coord[diri_dof_to_vnode])

    # Describe interpolation between nodes
    rid_to_reffe = GT.reference_fes(V)
    face_to_rid = GT.face_reference_id(V)

    # rid_to_refquad -> `Cuadrature`
    rid_to_refquad = GT.reference_quadratures(dΩ)

    # Gather node indexes of each face
    cell_to_nodes = GT.face_nodes(mesh,D)
    
    # Get node coordinates in the mesh
    node_to_coords = GT.node_coordinates(mesh)
    
    # rid_to_refface -> GenericLagrangeMeshFace
    rid_to_refface = GT.reference_faces(mesh,D)
    cell_to_rid = GT.face_reference_id(mesh,D)
    face_to_cell = GT.faces(Ω)
    nfaces = length(face_to_cell)

    # Returns a series of callables containing recipes on how to build the matrix and vector
    # TODO the strategies need some refactoring
    matrix_strategy =  GT.monolithic_matrix_assembly_strategy()
    vector_strategy = GT.monolithic_vector_assembly_strategy()

    # Consider only interpolation points which are not in the boundary
    free_dofs = GT.free_dofs(V)
    T = Float64

    # Initialize stiffness matrix and load vector
    matrix_strategy_init = matrix_strategy.init(free_dofs,free_dofs,T)
    vector_strategy_init = vector_strategy.init(free_dofs,T)
    field_i = field_j = 1

    # Count
    matrix_strategy_counter = matrix_strategy.counter(matrix_strategy_init)
    vector_strategy_counter = vector_strategy.counter(vector_strategy_init)

    # Count the number of terms in the matrix -> important for defining storage sparse matrix format later
    @timeit timer "count_terms" begin
        for face in 1:nfaces
            dofs = face_to_dofs[face]
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
    end

    # Allocate
    matrix_strategy_alloc = matrix_strategy.allocate(matrix_strategy_counter,matrix_strategy_init)
    vector_strategy_alloc = vector_strategy.allocate(vector_strategy_counter,vector_strategy_init)

    # Pre compute
    # NB we are assuming that the number of reference elements in the mesh is the same as in the
    # finite element space
    # In general we would need 3 kinds of reference objects
    # 1) reference geometries, 2) reference mesh faces, and 3) reference fes
    # All these concepts are in the code, but we are not using them now (maybe we should?)
    T = Float64
    Tgrad = SVector{D,T}
    Tx = Tgrad
    TJ = typeof(zero(Tx)*transpose(zero(Tx)))

    @timeit timer "cache_quantities" begin
        rid_to_cache = map(rid_to_reffe,rid_to_refface,rid_to_refquad) do reffe,refface,refquad
            
            # Set matrix and vector dimensions as the number of out-of-boundary nodes in the reference FE space
            ndofs = GT.num_dofs(reffe)

            A = zeros(T,ndofs,ndofs)
            b = zeros(T,ndofs)
            ue = zeros(T,ndofs)

            # Get the quadrature points
            x = GT.coordinates(refquad)

            # Calculate the gradient of the basis functions at the quadrature points
            ref_grads = GT.tabulator(reffe)(ForwardDiff.gradient,x)

            # Calculate the value of the basis functions at the quadrature points
            ref_vals = GT.tabulator(reffe)(GT.value,x)

            # Evaluate the values and gradients of basis functions for the mesh
            ref_vals_geom = GT.tabulator(refface)(GT.value,x)
            ref_grads_geom = GT.tabulator(refface)(ForwardDiff.gradient,x)

            phys_grads = zeros(Tgrad,ndofs)
            npoints = length(x)
            weights = GT.weights(refquad)
            (;A,b,ue,ref_vals_geom,ref_grads_geom,phys_grads,ref_vals,ref_grads,npoints,ndofs,weights)
        end
    end

    # Fill
    matrix_strategy_counter = matrix_strategy.counter(matrix_strategy_init)
    vector_strategy_counter = vector_strategy.counter(vector_strategy_init)
    z = zero(T)

    # Calculate matrix and vector for each face first and then summed over the terms of all faces
    @timeit timer "fill_and_assembly" begin
        for face in 1:nfaces
            rid = face_to_rid[face]
            cache = rid_to_cache[rid]
            
            # Face degrees of freedom
            dofs = face_to_dofs[face]

            # Set face index
            cell = face_to_cell[face]
            
            # Definition of the face in terms of node ids
            nodes = cell_to_nodes[cell]

            fill!(cache.A,z)
            fill!(cache.b,z)
            fill!(cache.ue,z)

            ndofs = cache.ndofs
            
            # Fill in Dirichlet values
            # Use face dofs to assign the initial values of uh at the boundary
            for i in 1:ndofs
                dof = dofs[i]
                if dof < 0
                    cache.ue[i] = x_diri[-dof]
                end
            end

            # Integration loop
            for point in 1:cache.npoints  #-> loop over quadrature points
                # Compute integration coordinate and Jacobian transpose
                x = zero(Tx)
                Jt = zero(TJ)
                for (inode,node) in enumerate(nodes)
                    coords = node_to_coords[node]
                    x += cache.ref_vals_geom[point,inode]*coords
                    Jt += coords*transpose(cache.ref_grads_geom[point,inode])
                end
                # Integration weight
                dV = abs(det(Jt))*cache.weights[point]
                # Compute physical gradients
                invJt = inv(Jt)
                for i in 1:ndofs
                    cache.phys_grads[i] = invJt*cache.ref_grads[point,i]
                end
                # Integrate matrix and vector
                fx = f(x)
                for i in 1:ndofs
                    cache.b[i] += cache.ref_vals[i]*fx*dV
                    for j in 1:ndofs
                        Aij = cache.phys_grads[i]⋅cache.phys_grads[j]*dV
                        cache.A[i,j] += Aij
                        cache.b[i] -= Aij*cache.ue[j]
                    end
                end
            end
            # Assemble the matrix and vector
            for (i,dof_i) in enumerate(dofs)
                vector_strategy_counter =
                vector_strategy.set!(vector_strategy_alloc,
                                     vector_strategy_counter,
                                     vector_strategy_init,
                                     cache.b[i],
                                     dof_i,field_i)
                for (j,dof_j) in enumerate(dofs)
                    matrix_strategy_counter =
                    matrix_strategy.set!(matrix_strategy_alloc,
                                         matrix_strategy_counter,
                                         matrix_strategy_init,
                                         cache.A[i,j],
                                         dof_i,dof_j,field_i,field_j)
                end
            end
        end

    end
    # Build the linear system and solve it
    # TODO do not return rhs_cache, and matrix_cache by default
    rhs, rhs_cache = vector_strategy.compress(vector_strategy_alloc,vector_strategy_init)
    matrix, matrix_cache = matrix_strategy.compress(matrix_strategy_alloc,matrix_strategy_init)
    p = PS.linear_problem(x_free,matrix,rhs)


    @timeit timer "solver" begin
        s = params[:solver](p)
        s = PS.solve(s)
        x_free = PS.solution(s)
        uh = GT.solution_field(uh,x_free)
    end

    # Integrate error norms
    el2 = zero(T)
    eh1 = zero(T)

    @timeit timer "error_norms" begin
        for face in 1:nfaces
            rid = face_to_rid[face]
            cache = rid_to_cache[rid]
            dofs = face_to_dofs[face]
            cell = face_to_cell[face]
            nodes = cell_to_nodes[cell]
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
            # Integration loop
            for point in 1:cache.npoints
                # Compute integration coordinate and Jacobian transpose
                x = zero(Tx)
                Jt = zero(TJ)
                for (inode,node) in enumerate(nodes)
                    coords = node_to_coords[node]
                    x += cache.ref_vals_geom[point,inode]*coords
                    Jt += coords*transpose(cache.ref_grads_geom[point,inode])
                end
                # Integration weight
                dV = abs(det(Jt))*cache.weights[point]
                # Compute physical gradients
                invJt = inv(Jt)
                for i in 1:ndofs
                    cache.phys_grads[i] = invJt*cache.ref_grads[point,i]
                end
                # Compute FE solution and its gradients
                uhx = z
                ∇uhx = zero(Tx)
                for i in 1:ndofs
                    ui = cache.ue[i]
                    uhx += cache.ref_vals[point,i]*ui
                    ∇uhx += cache.phys_grads[i]*ui
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

    end
    results[:error_h1_norm] = eh1
    results[:error_l2_norm] = el2

    @timeit timer "vtk" if params[:export_vtu]
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
        vis_node = 0
        for face in 1:nfaces
            rid = face_to_rid[face]
            cache = rid_to_vis_cache[rid]
            dofs = face_to_dofs[face]
            cell = face_to_cell[face]
            nodes = cell_to_nodes[cell]
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

    if params[:verbose]
        display(timer)
    end

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

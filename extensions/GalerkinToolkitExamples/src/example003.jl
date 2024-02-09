module Example003

import GalerkinToolkit as gk
import ForwardDiff
using StaticArrays
using LinearAlgebra
using SparseArrays
using WriteVTK
using PartitionedArrays: JaggedArray, sparse_matrix, sparse_matrix!
using TimerOutputs
using NLsolve
using Random
using GalerkinToolkitExamples: Example001

using Preconditioners
using IterativeSolvers: cg!

# This one implements a vanilla sequential iso-parametric p-Laplacian solver by only
# using the mesh interface.

function main(params_in)

    # Process params
    params_default = default_params()
    params = add_default_params(params_in,params_default)
    results = Dict{Symbol,Any}()
    timer = params[:timer]

    # Setup main data structures
    @timeit timer "setup" state = setup(params)
    add_basic_info(results,params,state)

    @timeit timer "solve_problem" x = solve_problem(params,state)

    # Post process
    @timeit timer "setup_uh" uh = setup_uh(x,state)
    @timeit timer "integrate_error_norms" integrate_error_norms(results,uh,state)
    @timeit timer "export_results" export_results(uh,params,state)

    print_timer(timer,allocations=false)

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
    outdir = mkpath(joinpath(@__DIR__,"..","output"))
    params = Dict{Symbol,Any}()
    params[:mesh] = gk.cartesian_mesh((0,1,0,1),(10,10))
    params[:u] = (x) -> sum(x)
    params[:f] = (x) -> 0.0
    params[:g] = (x) -> 0.0
    params[:dirichlet_tags] = ["boundary"]
    params[:neumann_tags] = String[]
    params[:integration_degree] = 2
    params[:p] = 2
    params[:solver] = nlsolve_solver()
    params[:example_path] = joinpath(outdir,"example003")
    params[:export_vtu] = true
    params[:timer] = TimerOutput()
    params
end

function nlsolve_solver(;linear_solver=Example001.lu_solver(),timer=TimerOutput(),options...)
    function setup(x0,nlp)
        @timeit timer "linearize" r0,J0,cache = nlp.linearize(x0)
        dx = similar(r0,axes(J0,2)) # TODO is there any way of reusing this in the nonlinear solve?
        @timeit timer "linear_solver_setup" ls_setup = linear_solver.setup(dx,J0,r0)
        function linsolve(x,A,b)
            # TODO we dont need to re-setup for the first
            # linear solve.
            # This can be avoided with a Ref{Bool} shared
            # between this function and j!
            @timeit timer "linear_solver_setup!"  linear_solver.setup!(ls_setup,A)
            @timeit timer "linear_solver_solve!" linear_solver.solve!(x,ls_setup,b)
            x
        end
        f!(r,x) = nlp.residual!(r,x,cache)
        j!(J,x) = nlp.jacobian!(J,x,cache)
        df = OnceDifferentiable(f!,j!,x0,r0,J0)
        (;df,linsolve,J0,cache,ls_setup,linear_solver)
    end
    function solve!(x,setup)
        (;df,linsolve) = setup
        result = nlsolve(df,x;linsolve,options...)
        x .= result.zero
        nothing
    end
    function setup!(setup,x0,nlp)
        error("todo")
    end
    function finalize!(setup)
        setup.linear_solver.finalize!(setup)
    end
    (;setup,solve!,setup!,finalize!)
end

function setup_dirichlet_bcs(params)
    mesh = params[:mesh]
    u = params[:u]
    dirichlet_tags = params[:dirichlet_tags]
    node_to_tag = zeros(gk.num_nodes(mesh))
    tag_to_name = dirichlet_tags
    gk.classify_mesh_nodes!(node_to_tag,mesh,tag_to_name)
    free_and_dirichlet_nodes = gk.partition_from_mask(i->i==0,node_to_tag)
    node_to_x = gk.node_coordinates(mesh)
    dirichlet_nodes = last(free_and_dirichlet_nodes)
    x_dirichlet = view(node_to_x,dirichlet_nodes)
    u_dirichlet = u.(x_dirichlet)
    dirichlet_bcs = (;free_and_dirichlet_nodes,u_dirichlet,node_to_tag)
end

function setup_integration(params,objects)
    mesh = params[:mesh]
    degree = params[:integration_degree]
    D = gk.num_dims(mesh)
    if objects === :cells
        d = D
    elseif objects === :faces
        d = D-1
    else
        error("")
    end
    # Integration
    ref_cells = gk.reference_faces(mesh,d)
    face_to_rid = gk.face_reference_id(mesh,d)
    integration_rules = map(ref_cells) do ref_cell
        gk.default_quadrature(gk.geometry(ref_cell),degree)
    end
    rid_to_weights = map(gk.weights,integration_rules)
    rid_to_coords = map(gk.coordinates,integration_rules)
    integration = (;rid_to_weights,rid_to_coords,face_to_rid,d)
end

function setup_isomap(params,integration)
    rid_to_coords = integration.rid_to_coords
    face_to_rid = integration.face_to_rid
    d = integration.d
    mesh = params[:mesh]
    ref_cells = gk.reference_faces(mesh,d)
    shape_funs = map(rid_to_coords,ref_cells) do q,ref_cell
        shape_vals = gk.tabulator(ref_cell)(gk.value,q)
        shape_grads = gk.tabulator(ref_cell)(ForwardDiff.gradient,q)
        shape_vals, shape_grads
    end
    rid_to_shape_vals = map(first,shape_funs)
    rid_to_shape_grads = map(last,shape_funs)
    face_to_nodes = JaggedArray(gk.face_nodes(mesh,d))
    node_to_coords = gk.node_coordinates(mesh)
    isomap = (;face_to_nodes,node_to_coords,face_to_rid,rid_to_shape_vals,rid_to_shape_grads,d)
end

function setup_neumann_bcs(params)
    mesh = params[:mesh]
    neumann_tags = params[:neumann_tags]
    neum_face_to_face = Int[]
    if length(neumann_tags) != 0
        D = gk.num_dims(mesh)
        tag_to_groups = gk.physical_faces(mesh,D-1)
        neum_face_to_face = reduce(union,map(tag->tag_to_groups[tag],neumann_tags))
    end
    (;neum_face_to_face)
end

function setup_user_funs(params)
    u = params[:u]
    f = params[:f]
    g = params[:g]
    user_funs = (;u,f,g)
end

function setup_dofs(params,dirichlet_bcs,cell_isomap,face_isomap)
    free_and_dirichlet_nodes = dirichlet_bcs.free_and_dirichlet_nodes
    node_to_free_node = gk.permutation(free_and_dirichlet_nodes)
    n_free = length(first(free_and_dirichlet_nodes))
    n_dofs = n_free
    cell_to_nodes = cell_isomap.face_to_nodes
    face_to_nodes = face_isomap.face_to_nodes
    function node_to_dof(node)
        free_node = node_to_free_node[node]
        if free_node <= n_free
            return Int(free_node)
        end
        Int(-(free_node-n_free))
    end
    cell_to_dofs_data = node_to_dof.(cell_to_nodes.data)
    face_to_dofs_data = node_to_dof.(face_to_nodes.data)
    cell_to_dofs = JaggedArray(cell_to_dofs_data,cell_to_nodes.ptrs)
    face_to_dofs = JaggedArray(face_to_dofs_data,face_to_nodes.ptrs)
    dofs = (;cell_to_dofs,face_to_dofs,n_dofs)
    dofs
end

function setup(params)
    timer = params[:timer]
    @timeit timer "dirichlet_bcs" dirichlet_bcs = setup_dirichlet_bcs(params)
    @timeit timer "neumann_bcs" neumann_bcs = setup_neumann_bcs(params)
    @timeit timer "cell_integration" cell_integration = setup_integration(params,:cells)
    @timeit timer "face_integration" face_integration = setup_integration(params,:faces)
    @timeit timer "cell_isomap" cell_isomap = setup_isomap(params,cell_integration)
    @timeit timer "face_isomap" face_isomap = setup_isomap(params,face_integration)
    @timeit timer "dofs" dofs = setup_dofs(params,dirichlet_bcs,cell_isomap,face_isomap)
    @timeit timer "user_funs" user_funs = setup_user_funs(params)
    solver = params[:solver]
    p = params[:p]
    state = (;p,timer,solver,dirichlet_bcs,neumann_bcs,cell_integration,cell_isomap,face_integration,face_isomap,user_funs,dofs)
end

function add_basic_info(results,params,state)
    node_to_x = state.cell_isomap.node_to_coords
    free_and_dirichlet_nodes = state.dirichlet_bcs.free_and_dirichlet_nodes
    cell_to_rid = state.cell_isomap.face_to_rid
    nfree = length(first(free_and_dirichlet_nodes))
    nnodes = length(node_to_x)
    ncells = length(cell_to_rid)
    results[:nnodes] = nnodes
    results[:nfree] = nfree
    results[:ncells] = ncells
    results
end

function nonlinear_problem(state)
    timer = state.timer
    function initial()
        n = state.dofs.n_dofs
        Random.seed!(1)
        rand(Float64,n)
    end
    function linearize(u_dofs)
        @timeit timer "assemble_sybmolic" r,J,V,K = assemble_sybmolic(state)
        cache = (V,K)
        @timeit timer "residual_cells!" residual_cells!(r,u_dofs,state)
        @timeit timer "residual_faces!" residual_faces!(r,u_dofs,state)
        @timeit timer "jacobian_cells!" jacobian_cells!(V,u_dofs,state)
        @timeit timer "sparse_matrix!" sparse_matrix!(J,V,K)
        r,J,cache
    end
    function residual!(r,u_dofs,cache)
        fill!(r,0)
        @timeit timer "residual_cells!" residual_cells!(r,u_dofs,state)
        @timeit timer "residual_faces!" residual_faces!(r,u_dofs,state)
        r
    end
    function jacobian!(J,u_dofs,cache)
        (V,K) = cache 
        @timeit timer "jacobian_cells!" jacobian_cells!(V,u_dofs,state)
        @timeit timer "sparse_matrix!" sparse_matrix!(J,V,K)
        J
    end
    (;initial,linearize,residual!,jacobian!)
end

function assemble_sybmolic(state)
    cell_to_dofs = state.dofs.cell_to_dofs
    n_dofs = state.dofs.n_dofs

    n_coo = 0
    ncells = length(cell_to_dofs)
    for cell in 1:ncells
        dofs = cell_to_dofs[cell]
        ndofs = length(dofs)
        n_coo += ndofs*ndofs
    end

    I_coo = Vector{Int32}(undef,n_coo)
    J_coo = Vector{Int32}(undef,n_coo)
    V_coo = Vector{Float64}(undef,n_coo)
    r = zeros(Float64,state.dofs.n_dofs)

    n_coo = 0
    for cell in 1:ncells
        dofs = cell_to_dofs[cell]
        ndofs = length(dofs)
        for i in 1:ndofs
            dofs_i = dofs[i]
            for j in 1:ndofs
                n_coo += 1
                I_coo[n_coo] = dofs_i
                J_coo[n_coo] = dofs[j]
            end
        end
    end

    J,K = sparse_matrix(I_coo,J_coo,V_coo,n_dofs,n_dofs;reuse=true)
    r,J,V_coo,K
end

@inline function p_laplace_residual(∇u,∇dv,dv,f,p)
    ∇dv⋅((norm(∇u)^(p-2))*∇u) - f*dv
end

@inline function p_laplace_jacobian(∇u,∇du,∇dv,p)
    pm2 = p-2
    pm4 = p-4
    ∇dv⋅((norm(∇u)^pm2)*∇du) + ∇dv⋅(pm2*(norm(∇u)^pm4)*(∇u⋅∇du)*∇u)
end

function residual_cells!(r,u_dofs,state)

    cell_to_dofs = state.dofs.cell_to_dofs
    cell_to_nodes = state.cell_isomap.face_to_nodes
    cell_to_rid = state.cell_isomap.face_to_rid
    free_and_dirichlet_nodes = state.dirichlet_bcs.free_and_dirichlet_nodes
    node_to_x = state.cell_isomap.node_to_coords
    ∇s = state.cell_isomap.rid_to_shape_grads
    s = state.cell_isomap.rid_to_shape_vals
    u_dirichlet = state.dirichlet_bcs.u_dirichlet
    f = state.user_funs.f
    w = state.cell_integration.rid_to_weights
    d = state.cell_integration.d
    u_dirichlet = state.dirichlet_bcs.u_dirichlet
    ncells = length(cell_to_nodes)
    p = state.p

    # Allocate auxiliary buffers
    ∇x = map(i->similar(i,size(i,2)),∇s)
    Aes = map(i->zeros(size(i,2),size(i,2)),∇s)
    fes = map(i->zeros(size(i,2)),∇s)
    ues = map(i->zeros(size(i,2)),∇s)
    ∇st = map(m->collect(permutedims(m)),∇s)
    st = map(m->collect(permutedims(m)),s)

    Tx = eltype(node_to_x)
    TJ = typeof(zero(Tx)*zero(Tx)')
    T = eltype(Tx)

    i_coo = 0
    for cell in 1:ncells
        rid = cell_to_rid[cell]
        nodes = cell_to_nodes[cell]
        dofs = cell_to_dofs[cell]
        Ae = Aes[rid]
        fe = fes[rid]
        ue = ues[rid]
        nl = length(nodes)
        ∇ste = ∇st[rid]
        ∇xe = ∇x[rid]
        ste = st[rid]
        we = w[rid]
        nq = length(we)
        fill!(Ae,zero(eltype(Ae)))
        fill!(fe,zero(eltype(Ae)))
        # Integrate
        for iq in 1:nq
            Jt = zero(TJ) 
            xint = zero(Tx)
            for k in 1:nl
                x = node_to_x[nodes[k]]
                ∇sqx = ∇ste[k,iq]
                sqx = ste[k,iq]
                Jt += ∇sqx*x'
                xint += sqx*x
            end
            detJt = det(Jt)
            invJt = inv(Jt)
            dV = abs(detJt)*we[iq]
            for k in 1:nl
                dof_k = dofs[k]
                if dof_k < 1
                    uk = u_dirichlet[-dof_k]
                    ue[k] = uk
                    continue
                end
                ue[k] = u_dofs[dof_k]
            end
            ∇u = zero(Tx)
            for k in 1:nl
                uek = ue[k] 
                ∇xek = invJt*∇ste[k,iq]
                ∇u += ∇xek*uek
                ∇xe[k] = ∇xek
            end
            fx = f(xint)
            for k in 1:nl
                dv = ste[k,iq]
                ∇dv = ∇xe[k]
                fe[k] += ( p_laplace_residual(∇u,∇dv,dv,fx,p) )*dV
            end
        end
        for i in 1:nl
            if !(dofs[i]>0)
                continue
            end
            r[dofs[i]] += fe[i]
        end
    end
end

function jacobian_cells!(V_coo,u_dofs,state)

    cell_to_dofs = state.dofs.cell_to_dofs
    cell_to_nodes = state.cell_isomap.face_to_nodes
    cell_to_rid = state.cell_isomap.face_to_rid
    free_and_dirichlet_nodes = state.dirichlet_bcs.free_and_dirichlet_nodes
    node_to_x = state.cell_isomap.node_to_coords
    ∇s = state.cell_isomap.rid_to_shape_grads
    s = state.cell_isomap.rid_to_shape_vals
    u_dirichlet = state.dirichlet_bcs.u_dirichlet
    f = state.user_funs.f
    w = state.cell_integration.rid_to_weights
    d = state.cell_integration.d
    u_dirichlet = state.dirichlet_bcs.u_dirichlet
    ncells = length(cell_to_nodes)
    p = state.p

    # Allocate auxiliary buffers
    ∇x = map(i->similar(i,size(i,2)),∇s)
    Aes = map(i->zeros(size(i,2),size(i,2)),∇s)
    fes = map(i->zeros(size(i,2)),∇s)
    ues = map(i->zeros(size(i,2)),∇s)
    ∇st = map(m->collect(permutedims(m)),∇s)
    st = map(m->collect(permutedims(m)),s)

    Tx = eltype(node_to_x)
    TJ = typeof(zero(Tx)*zero(Tx)')
    T = eltype(Tx)

    i_coo = 0
    for cell in 1:ncells
        rid = cell_to_rid[cell]
        nodes = cell_to_nodes[cell]
        dofs = cell_to_dofs[cell]
        Ae = Aes[rid]
        fe = fes[rid]
        ue = ues[rid]
        nl = length(nodes)
        ∇ste = ∇st[rid]
        ∇xe = ∇x[rid]
        ste = st[rid]
        we = w[rid]
        nq = length(we)
        fill!(Ae,zero(eltype(Ae)))
        fill!(fe,zero(eltype(Ae)))
        # Integrate
        for iq in 1:nq
            Jt = zero(TJ) 
            for k in 1:nl
                x = node_to_x[nodes[k]]
                ∇sqx = ∇ste[k,iq]
                sqx = ste[k,iq]
                Jt += ∇sqx*x'
            end
            detJt = det(Jt)
            invJt = inv(Jt)
            dV = abs(detJt)*we[iq]
            for k in 1:nl
                dof_k = dofs[k]
                if dof_k < 1
                    uk = u_dirichlet[-dof_k]
                    ue[k] = uk
                    continue
                end
                ue[k] = u_dofs[dof_k]
            end
            ∇u = zero(Tx)
            for k in 1:nl
                uek = ue[k]
                ∇xek = invJt*∇ste[k,iq]
                ∇u += ∇xek*uek
                ∇xe[k] = ∇xek
            end
            for j in 1:nl
                ∇du = ∇xe[j]
                for i in 1:nl
                    ∇dv = ∇xe[i]
                    Ae[i,j] += ( p_laplace_jacobian(∇u,∇du,∇dv,p) )*dV
                end
            end
        end
        # Set the result in the output array
        for i in 1:nl
            for j in 1:nl
                i_coo += 1
                V_coo[i_coo] = Ae[i,j]
            end
        end
    end
end

function residual_faces!(r,u_dofs,state)

    neum_face_to_face = state.neumann_bcs.neum_face_to_face
    if length(neum_face_to_face) == 0
        return nothing
    end

    face_to_rid = state.face_isomap.face_to_rid
    face_to_dofs = state.dofs.face_to_dofs
    face_to_nodes = state.face_isomap.face_to_nodes
    ∇s_f = state.face_isomap.rid_to_shape_grads
    s_f = state.face_isomap.rid_to_shape_vals
    node_to_x = state.face_isomap.node_to_coords
    g = state.user_funs.g
    w_f = state.face_integration.rid_to_weights
    d = state.cell_integration.d

    fes = map(i->zeros(size(i,2)),s_f)
    Tx = SVector{d,Float64}
    Ts = typeof(zero(SVector{d-1,Float64})*zero(Tx)')

    for face in neum_face_to_face
        rid = face_to_rid[face]
        nodes = face_to_nodes[face]
        dofs = face_to_dofs[face]
        ∇se = ∇s_f[rid]
        se = s_f[rid]
        we = w_f[rid]
        nq = length(we)
        fe = fes[rid]
        fill!(fe,zero(eltype(fe)))
        nl = length(nodes)
        for iq in 1:nq
            Jt = zero(Ts) 
            xint = zero(Tx)
            for k in 1:nl
                x = node_to_x[nodes[k]]
                ∇sqx = ∇se[iq,k]
                sqx = se[iq,k]
                Jt += ∇sqx*x'
                xint += sqx*x
            end
            J = transpose(Jt)
            dS = sqrt(det(Jt*J))*we[iq]
            gx = g(xint)
            for k in 1:nl
                fe[k] +=  gx*se[iq,k]*dS
            end
        end
        for i in 1:nl
            if !(dofs[i]>0)
                continue
            end
            r[dofs[i]] -= fe[i]
        end
    end
end

function solve_problem(params,state)
    timer = state.timer
    problem = nonlinear_problem(state)
    solver = params[:solver]
    x = problem.initial()
    @timeit timer "solver.setup" setup = solver.setup(x,problem)
    @timeit timer "solver.solve!" solver.solve!(x,setup)
    @timeit timer "solver.finalize!" solver.finalize!(setup)
    x
end

function setup_uh(x,state)
    free_and_dirichlet_nodes = state.dirichlet_bcs.free_and_dirichlet_nodes
    node_to_x = state.cell_isomap.node_to_coords
    u_free = x
    u_dirichlet = state.dirichlet_bcs.u_dirichlet
    nnodes = length(node_to_x)
    uh = zeros(nnodes)
    uh[first(free_and_dirichlet_nodes)] = u_free
    uh[last(free_and_dirichlet_nodes)] = u_dirichlet
    uh
end

function integrate_error_norms_loop(uh,state)

    cell_to_nodes = state.cell_isomap.face_to_nodes
    cell_to_rid = state.cell_isomap.face_to_rid
    free_and_dirichlet_nodes = state.dirichlet_bcs.free_and_dirichlet_nodes
    node_to_x = state.cell_isomap.node_to_coords
    ∇s = state.cell_isomap.rid_to_shape_grads
    s = state.cell_isomap.rid_to_shape_vals
    u_dirichlet = state.dirichlet_bcs.u_dirichlet
    u = state.user_funs.u
    w = state.cell_integration.rid_to_weights
    d = state.cell_integration.d

    eh1 = 0.0
    el2 = 0.0
    ncells = length(cell_to_rid)
    ues = map(i->zeros(size(i,2)),∇s)
    ∇x = map(i->similar(i,size(i,2)),∇s)
    Tx = eltype(node_to_x)
    TJ = typeof(zero(Tx)*zero(Tx)')
    for cell in 1:ncells
        rid = cell_to_rid[cell]
        nodes = cell_to_nodes[cell]
        ue = ues[rid]
        nl = length(nodes)
        ∇se = ∇s[rid]
        ∇xe = ∇x[rid]
        se = s[rid]
        we = w[rid]
        nq = length(we)
        for k in 1:nl
            nk = nodes[k]
            ue[k] = uh[nk]
        end
        for iq in 1:nq
            Jt = zero(TJ) 
            xint = zero(Tx)
            for k in 1:nl
                x = node_to_x[nodes[k]]
                ∇sqx = ∇se[iq,k]
                sqx = se[iq,k]
                Jt += ∇sqx*x'
                xint += sqx*x
            end
            detJt = det(Jt)
            invJt = inv(Jt)
            dV = abs(detJt)*we[iq]
            for k in 1:nl
                ∇xe[k] = invJt*∇se[iq,k]
            end
            ux = u(xint)
            ∇ux = ForwardDiff.gradient(u,xint)
            ∇uhx = zero(∇ux)
            uhx = zero(ux)
            for k in 1:nl
                uek = ue[k]
                ∇uhx += uek*∇xe[k]
                uhx += uek*se[iq,k]
            end
            ∇ex = ∇ux - ∇uhx
            ex =  ux - uhx
            eh1 += (∇ex⋅∇ex + ex*ex)*dV
            el2 += ex*ex*dV
        end
    end
    eh1, el2
end

function integrate_error_norms(results,uh,state)
    eh1², el2² = integrate_error_norms_loop(uh,state)
    eh1 = sqrt(eh1²)
    el2 = sqrt(el2²)
    results[:eh1] = eh1
    results[:el2] = el2
    results
end

function export_results(uh,params,state)
    if ! params[:export_vtu]
        return nothing
    end
    mesh = params[:mesh]
    example_path = params[:example_path]
    node_to_tag = state.dirichlet_bcs.node_to_tag
    vtk_grid(example_path,gk.vtk_args(mesh)...) do vtk
        gk.vtk_physical_faces!(vtk,mesh)
        vtk["tag"] = node_to_tag
        vtk["uh"] = uh
        vtk["dim"] = gk.face_dim(mesh)
    end
    nothing
end

end # module

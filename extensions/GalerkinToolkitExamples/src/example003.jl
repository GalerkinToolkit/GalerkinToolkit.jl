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

using CUDA
using Serialization
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

    @timeit timer "solve_problem" x = solve_problem(params,state,results)

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
    params[:solver] = nlsolve_solver(;method=:newton)
    params[:example_path] = joinpath(outdir,"example003")
    params[:export_vtu] = true
    params[:autodiff] = :hand
    params[:timer] = TimerOutput()
    params[:jacobian] = :original
    params[:parallelization_level] = :cell
    params[:mem_layout] = :cell_major
    params
end

function nlsolve_solver(;linear_solver=Example001.lu_solver(),timer=TimerOutput(),options...)
    function setup(x0,nlp,params)
        @timeit timer "linearize" r0,J0,cache = nlp.linearize(x0, params)
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
        j!(J,x) = nlp.jacobian!(J,x,cache) # This line calls the jacobian function (specified in _tests).
        df = OnceDifferentiable(f!,j!,x0,r0,J0)
        (;df,linsolve,J0,cache,ls_setup,linear_solver)
    end
    function solve!(x,setup)
        (;df,linsolve) = setup
        result = nlsolve(df,x;linsolve,options...)
        x .= result.zero
        result.f_calls
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
# two versions of the setup, one with double precision. One for residual and one for jacobian. 
# Two parts to it, change the values of cell_to_nodes etc in the setup and then pass a parameter that sets up xe etc. as Float32
function setup_xe(cell_isomap, cell_integration, params)
    cell_to_nodes = cell_isomap.face_to_nodes
    node_to_coords = cell_isomap.node_to_coords
    mem_layout = params[:mem_layout]
    ncells = length(cell_to_nodes)
    nl = length(cell_to_nodes[1])
    d = length(eltype(node_to_coords))
    # dof as the inner loop, since dofs per cell are looped over in assembly.
    xe = Vector{SVector{d, Float64}}(undef, ncells * nl) 
    if mem_layout == :cell_major
        for cell in 1:ncells
            for dof in 1:nl
                index = (cell - 1) * nl + dof
                node = cell_to_nodes[cell][dof]
                xe[index] = node_to_coords[node]
            end
        end
    elseif mem_layout == :dof_major # coalesced for gpu memory
        for dof in 1:nl
            for cell in 1:ncells
                index = (dof - 1) * ncells + cell
                node = cell_to_nodes[cell][dof]
                xe[index] = node_to_coords[node]
            end
        end
    end
    setup_xe = (;xe)
end

function setup_ue!(ue, cell_to_dofs, u_dirichlet, u_dofs, ncells, nl, mem_layout)
    if mem_layout == :cell_major
        index_function = (cell, dof, nl, ncells) -> (cell - 1) * nl + dof
    elseif mem_layout == :dof_major # coalesced for gpu access.
        index_function = (cell, dof, nl, ncells) -> (dof - 1) * ncells + cell
    end

    for cell in 1:ncells
        for dof in 1:nl
            index = index_function(cell, dof, nl, ncells)
            # Then fill in the ue based on dof_k type.
            dof_k = cell_to_dofs[cell][dof]
            # Here only fill in non-dirichlet values.
            if dof_k < 1
                uk = u_dirichlet[-dof_k]
                ue[index] = uk
                continue
            end
            ue[index] = u_dofs[dof_k]
        end
    end
end

function ue_allocation(cell_isomap)
    cell_to_nodes = cell_isomap.face_to_nodes
    ncells = length(cell_to_nodes) # number of elements
    nl = length(cell_to_nodes[1]) # number of nodes per element
    
    # Allocate memory for the element u
    ue = Vector{Float64}(undef, ncells * nl)
    ue_alloc = (;ue)
end

function cpu_gpu_transfer(V_coo, ∇ste, w, xe, Jt, ∇u, ue, ncells, nq, nl, p)
    V_coo_d = cu(V_coo) 
    ∇ste_d = cu(cu.(∇ste))
    w_d = cu(w)
    xe_d = cu(xe) 
    TJ_d = cu(Jt) 
    Tx_d = cu(cu.(∇u)) 
    ue_d = cu(ue) 
    ncells_d = cu(ncells)
    nq_d = cu(nq)
    nl_d = cu(nl)
    p_d = cu(p)
    V_coo_d,∇ste_d,w_d,xe_d,ue_d,TJ_d,Tx_d,ncells_d,nl_d,nq_d,p_d
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
    @timeit timer "ue_allocation" ue_alloc = ue_allocation(cell_isomap)
    @timeit timer "node_setup" xe_setup = setup_xe(cell_isomap, cell_integration, params)
    solver = params[:solver]
    p = params[:p]
    if params[:autodiff] === :hand
        flux = flux_hand
        dflux = dflux_hand
    elseif params[:autodiff] === :flux
        flux = flux_hand
        dflux = dflux_from_flux
    elseif params[:autodiff] === :energy
        flux = flux_from_energy
        dflux = dflux_from_energy
    else
        error("not implemented: autodiff == $autodiff")
    end

    if params[:jacobian] === :original
        jacobian_cells! = jacobian_cells_original!
    elseif params[:jacobian] === :cpu_extension
        jacobian_cells! = jacobian_cells_cpu_extension!
    elseif params[:jacobian] === :gpu_extension
        jacobian_cells! = jacobian_cells_gpu!
    else
        error("unspecified jacobian function")
    end

    prl_level = params[:parallelization_level]
    mem_layout = params[:mem_layout]
    state = (;p,timer,flux,dflux,solver,dirichlet_bcs,neumann_bcs,cell_integration,cell_isomap,face_integration,face_isomap,user_funs,dofs,xe_setup,ue_alloc,jacobian_cells!,prl_level,mem_layout)
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
    function linearize(u_dofs, params)
        @timeit timer "assemble_sybmolic" r,J,V,K = assemble_sybmolic(state, params)
        cache = (V,K)
        @timeit timer "residual_cells!" residual_cells!(r,u_dofs,state)
        @timeit timer "residual_faces!" residual_faces!(r,u_dofs,state)
        @timeit timer "jacobian_cells!" state.jacobian_cells!(V,u_dofs,state)
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
        @timeit timer "jacobian_cells!" state.jacobian_cells!(V,u_dofs,state)
        @timeit timer "sparse_matrix!" sparse_matrix!(J,V,K)
        J
    end
    (;initial,linearize,residual!,jacobian!)
end

function assemble_sybmolic(state,params)
    cell_to_dofs = state.dofs.cell_to_dofs
    n_dofs = state.dofs.n_dofs
    mem_layout = params[:mem_layout]

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
    if mem_layout == :cell_major
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
    # You need to assume all cells are the same type here to make the outer loops work.
    elseif mem_layout == :dof_major # coalesced for gpu memory
        ndofs = length(cell_to_dofs[1])
        for j in 1:ndofs
            for i in 1:ndofs
                for cell in 1:ncells
                    dofs = cell_to_dofs[cell]
                    n_coo += 1
                    I_coo[n_coo] = dofs[i]
                    J_coo[n_coo] = dofs[j]
                end
            end
        end
    end
    J,K = sparse_matrix(I_coo,J_coo,V_coo,n_dofs,n_dofs;reuse=true)
    println(J, " ", K)
    r,J,V_coo,K
end

@inline flux_hand(∇u,p) = (norm(∇u)^(p-2))*∇u
@inline function dflux_hand(∇u,∇du,p)
    pm2 = p-2
    pm4 = p-4
    ((norm(∇u)^pm2)*∇du) + (pm2*(norm(∇u)^pm4)*(∇u⋅∇du)*∇u)
end
@inline function dflux_from_flux(∇u,∇du,p)
    f(x) = flux_hand(x,p)
    x = ∇u
    dx = ∇du
    dfdx = ForwardDiff.jacobian(f,x)
    dfdx*dx
end

@inline energy(∇u,p) = (1/p)*(norm(∇u)^p)
@inline function flux_from_energy(∇u,p)
    f(x) = energy(x,p)
    x = ∇u
    dfdx = ForwardDiff.gradient(f,x)
    dfdx
end
@inline function dflux_from_energy(∇u,∇du,p)
    f(x) = energy(x,p)
    x = ∇u
    dx = ∇du
    dfdx = ForwardDiff.hessian(f,x)
    dfdx*dx
end

function residual_cells!(r,u_dofs,state)
    flux = state.flux
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
                fe[k] += ( ∇dv⋅flux(∇u,p) - fx*dv )*dV
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

# Check which arrays go into which parts of the memory.
# profile on the gpu.
@inline function kernel_generic!(V_coo, cell, nq, nl, Jt, xe, ∇ste, we, ∇u, ue, dflux, p)
    # Integrate
    for iq in 1:nq 
        Jt_local = Jt # These are here so that the threads do not accumulate to same variable in global memory.
        ∇u_local = ∇u

        for k in 1:nl
            k_index = (cell - 1) * nl + k
            x = xe[k_index]
            ∇sqx = ∇ste[k,iq] 
            Jt_local += ∇sqx*x'
        end
        detJt = det(Jt_local)
        invJt = inv(Jt_local) # For the tiling this is tile x 1 _. make T x T (duplicate computation).
        dV = abs(detJt)*we[iq] 

        for k in 1:nl
            k_index = (cell - 1) * nl + k
            uek = ue[k_index]
            ∇xek = invJt*∇ste[k,iq] 
            ∇u_local += ∇xek*uek
        end
        # Next you need ∇du, ∇dv and dV, ∇u
        for j in 1:nl
            ∇du = invJt*∇ste[j,iq]
            for i in 1:nl
                ∇dv = invJt*∇ste[i,iq]
                i_coo = (cell - 1) * nl^2 + ((j - 1) * nl + i) 
                V_coo[i_coo] += (∇dv⋅dflux(∇u_local,∇du,p)) * dV
            end
        end
    end
end
# Maybe the logic should be to get xe + ue, then loop through the iq in parallel. Then the last loop is nl*nq.
# So in this case you have ncells * nl threads that are launched. In this version the ∇ste will never be a 
# problem because threads in a block is approx. 1024 while ∇ste << 1024. 

# So there are two parts to this -> loading values which is according to nl and computing
# which is done according to nq.
@inline function kernel_tiled!(V_coo, tid, nq, nl, Jt, xe, ∇ste, we, ∇u, ue, dflux, p, ncells)
    # variables to load into local memory: Jt, ∇u, ∇ste, xe, we, ue.
    # Think about this as how many values does each thread store (max 12 for 1024 threads per block).
    # Specify the shared memory per threadblock.
    # Assume (for now) nq == nl and that d = 2.
    # So the gain comes from not having to recompute invJt etc. for all iq points.
    # Note here that each shmem array should be 32 or a multiple of it (related to warp execution).
    max_threads = 1024
    Txe = CuStaticSharedArray(Float32, (2,max_threads)) # Two per thread (for the coordinates).
    Tue = CuStaticSharedArray(Float32, max_threads) # One per thread
    T∇ste = CuStaticSharedArray(Float32, (max_threads,2)) # This one is flattened using 2 indices (nq and nl)
    T∇u = CuStaticSharedArray(Float32, (2,max_threads))
    T∇dux = CuStaticSharedArray(Float32, (2,max_threads))
    TdV = CuStaticSharedArray(Float32, max_threads)
    
    # identify your position in nl block, nq and cell. tid is the grid value.
    cell = Int((ceil(tid/nl^2)-1) % ncells)+1 
    j = Int(((ceil(tid/nl)-1) % nl)+1) # Column 
    i = ((tid-1)%nl)+1 # Row
    iq = ((tid-1)%nq)+1 
    node_index = (cell - 1) * nl + i # Here you get node to use in ue/xe
    xe_col = tid % 2 + 1

    # here until row 638 is done for each iq in nq (NOTE THIS!)
    # Load the variables into shared memory from global memory
    TV_coo = 0 # auxiliary for each thread
    # If you want the vector inside xe next to each other then use div(tid + 1, 2)
    Txe[:,tid] = xe[node_index][:] # use global id to load into local id - so which i,j to operate from.
    Tue[tid] = ue[node_index]
    T∇ste[tid,:] = ∇ste[j,iq][:] # This works but very complicated.
    sync_threads() # sync to make sure all data (per threadblock) is loaded into memory

    # ---------------- SWITCH FROM NL INDEXING TO NQ INDEXING ---------------
    # Calculate your part of Jt -> determinants etc. + u_local accumulation
    for k in 1:nl 
        # So you want to round down to the first node in your column
        index = Int32(floor(tid / nl^2) * nl^2) + k # For the first index in the cell
        x = Txe[:, index] # Then you can get the k-th node in cell 
        ∇sqx = T∇ste[(k-1) * nl + iq,:] 
        Jt += ∇sqx*x'
    end
    detJt = det(Jt_local)
    invJt = inv(Jt_local) 
    # Both Jt and ∇u are calculated per iq point and per cell. 

    for k in 1:nl
        index = Int32(floor(tid / nl^2) * nl^2) + k # For the first index in the cell
        uek = Tue[index]
        ∇xek = invJt*T∇ste[(k-1) * nl + iq,:] 
        ∇u += ∇xek*uek # per iq and per cell.
    end

    T∇u[:,tid] = ∇u # per iq and per cell. Maybe a cyclic one here as well. 
    TdV[tid] = abs(detJt)*we[iq] # This one is per iq. So pick first value in cell and loop nq times.
    # If you precompute this then you need to store nq x nl values.
    ∇dux = invJt * T∇ste[(k-1) * nl + iq,:] # To keep the order stable. 1 1 1 2 2 2 for nl and 1 2 3 1 2 3 for iq.
    T∇dux[:,tid] = ∇dux # Precompute these into shared instead of storing invJt matrices.
    sync_threads() # sync so that all data is computed before adding to global memory.

    # Loop over the iq points and accumulate into an intermediate variable before writing to global mem.
    for iq in 1:nq # Each i,j index will loop over all the nq points
        index = Int32(floor(tid / nl^2) * nl^2) + iq # For first index in cell = first iq point
        dV = TdV[index]
        ∇u_local = T∇u[:,index]
        i_index = (cell - 1) * nl + nl * i + iq
        j_index = (cell - 1) * nl + nl * j + iq
        ∇du = T∇dux[:,i_index] # Still need to think about how these are stored.
        ∇dv = T∇dux[:,j_index] # So you jump nl steps into the cell and + the iq point.
        TV_coo += (∇dv⋅dflux(∇u_local,∇du,p)) * dV
    end 
    # Write the total onto global memory
    i_coo = (cell - 1) * nl^2 + ((j - 1) * nl + i) 
    V_coo[i_coo] = TV_coo
end

@inline function kernel_coalesced!(V_coo, cell, nq, nl, Jt, xe, ∇ste, we, ∇u, ue, dflux, p, ncells)
    for iq in 1:nq 
        Jt_local = Jt
        ∇u_local = ∇u

        for k in 1:nl
            k_index = (k - 1) * ncells + cell # switch the cell and k here to get right index.
            x = xe[k_index]
            ∇sqx = ∇ste[k,iq]
            Jt_local += ∇sqx*x'
        end
        detJt = det(Jt_local)
        invJt = inv(Jt_local)
        dV = abs(detJt)*we[iq]

        for k in 1:nl
            k_index = (k - 1) * ncells + cell # switch the cell and k here to get right index.
            uek = ue[k_index]
            ∇xek = invJt*∇ste[k,iq]
            ∇u_local += ∇xek*uek
        end

        for j in 1:nl
            ∇du = invJt*∇ste[j,iq]
            for i in 1:nl
                ∇dv = invJt*∇ste[i,iq]
                i_coo = cell + (j-1) * ncells * nl + (i-1) * ncells # There is also a change here to the calculation
                V_coo[i_coo] += (∇dv⋅dflux(∇u_local,∇du,p)) * dV
            end
        end
    end
end

@inline function kernel_quad!(V_coo, thread, nq, nl, Jt, xe, ∇ste, we, ∇u, ue, dflux, p, ncells)
    # First calc all cells for iq = 1, then all cells for iq = 2 etc...
    cell = ((thread-1)%ncells)+1
    iq = Int(((ceil(thread/ncells)-1) % ncells)+1)

    Jt_local = Jt
    ∇u_local = ∇u

    for k in 1:nl
        k_index = (cell - 1) * nl + k
        x = xe[k_index]
        ∇sqx = ∇ste[k,iq] # On a gpu, all threads are going for the same values here. 
        Jt_local += ∇sqx*x'
    end
    detJt = det(Jt_local)
    invJt = inv(Jt_local)
    dV = abs(detJt)*we[iq] # Same with w here.

    for k in 1:nl
        k_index = (cell - 1) * nl + k
        uek = ue[k_index] # Also make the other option here. 
        ∇xek = invJt*∇ste[k,iq] # here again all threads are reading the same array from memory at the same time. 
        ∇u_local += ∇xek*uek
    end

    for j in 1:nl
        ∇du = invJt*∇ste[j,iq]
        for i in 1:nl
            ∇dv = invJt*∇ste[i,iq]
            i_coo = (cell - 1) * nl^2 + ((j - 1) * nl + i)
            V_coo[i_coo] += (∇dv⋅dflux(∇u_local,∇du,p)) * dV # atomic add to avoid race conditions. 
            # CUDA.@atomic
        end
    end

end
# Here there is a race condition if you launch nq * ncells of threads. 
# add an atomic operation for the V_coo.
@inline function kernel_elem_j!(V_coo, thread, nq, nl, Jt, xe, ∇ste, we, ∇u, ue, dflux, p)
    cell = Int(ceil(thread/nl))
    j = ((thread-1)%nl)+1

    for iq in 1:nq
        Jt_local = Jt
        ∇u_local = ∇u
        for k in 1:nl
            k_index = (cell - 1) * nl + k
            x = xe[k_index]
            ∇sqx = ∇ste[k,iq]
            Jt_local += ∇sqx*x'
        end
        detJt = det(Jt_local)
        invJt = inv(Jt_local)
        dV = abs(detJt)*we[iq]

        for k in 1:nl
            k_index = (cell - 1) * nl + k
            uek = ue[k_index]
            ∇xek = invJt*∇ste[k,iq]
            ∇u_local += ∇xek*uek
        end

        ∇du = invJt*∇ste[j,iq]
        for i in 1:nl
            ∇dv = invJt*∇ste[i,iq]
            i_coo = (cell - 1) * nl^2 + ((j - 1) * nl + i)
            V_coo[i_coo] += (∇dv⋅dflux(∇u_local,∇du,p)) * dV # Almost all time spent here according to @ProfileView
        end
    end
end

@inline function kernel_elem_ij!(V_coo, thread, nq, nl, Jt, xe, ∇ste, we, ∇u, ue, dflux, p)
    # Here it would be possible to define one integer, accumulate and then one global memory write.
    # Set the indices of cell, i and j.
    cell = Int(ceil(thread/(nl^2)))
    i = ((thread-1)%nl)+1 
    j = Int(((ceil(thread/nl)-1) % nl)+1) 
    i_coo = (cell - 1) * nl^2 + ((j - 1) * nl + i) # Define this one only once. 
    # parallelize this loop and add atomic V_coo.
    for iq in 1:nq # quadrature points. 
        Jt_local = Jt
        ∇u_local = ∇u

        for k in 1:nl
            k_index = (cell - 1) * nl + k
            x = xe[k_index]
            ∇sqx = ∇ste[k,iq]
            Jt_local += ∇sqx*x'
        end
        detJt = det(Jt_local)
        invJt = inv(Jt_local)
        dV = abs(detJt)*we[iq]

        for k in 1:nl
            k_index = (cell - 1) * nl + k
            uek = ue[k_index]
            ∇xek = invJt*∇ste[k,iq]
            ∇u_local += ∇xek*uek
        end

        ∇du = invJt*∇ste[j,iq]
        ∇dv = invJt*∇ste[i,iq]
        V_coo[i_coo] += (∇dv⋅dflux(∇u_local,∇du,p)) * dV # Almost all time spent here according to @ProfileView
    end
end

@inline function kernel_full_prl!(V_coo, thread, nq, nl, Jt, xe, ∇ste, we, ∇u, ue, dflux, p, ncells)
    # Set the indices of cell, iq, i and j.
    cell = Int((ceil(thread/nl^2)-1) % ncells)+1
    i = ((thread-1)%nl)+1 
    j = Int(((ceil(thread/nl)-1) % nl)+1) 
    iq = Int(ceil(thread/(nl^2*ncells)))
    i_coo = (cell - 1) * nl^2 + ((j - 1) * nl + i) # Define this one only once. 

    Jt_local = Jt
    ∇u_local = ∇u

    for k in 1:nl
        k_index = (cell - 1) * nl + k
        x = xe[k_index]
        ∇sqx = ∇ste[k,iq]
        Jt_local += ∇sqx*x'
    end
    detJt = det(Jt_local)
    invJt = inv(Jt_local)
    dV = abs(detJt)*we[iq]

    for k in 1:nl
        k_index = (cell - 1) * nl + k
        uek = ue[k_index]
        ∇xek = invJt*∇ste[k,iq]
        ∇u_local += ∇xek*uek
    end

    ∇du = invJt*∇ste[j,iq]
    ∇dv = invJt*∇ste[i,iq]
    V_coo[i_coo] += (∇dv⋅dflux(∇u_local,∇du,p)) * dV # Almost all time spent here according to @ProfileView
    # CUDA.@atomic
end

function jacobian_cells_cpu_extension!(V_coo,u_dofs,state)
    dflux = state.dflux
    cell_to_dofs = state.dofs.cell_to_dofs
    cell_to_nodes = state.cell_isomap.face_to_nodes
    ∇s = state.cell_isomap.rid_to_shape_grads
    s = state.cell_isomap.rid_to_shape_vals
    u_dirichlet = state.dirichlet_bcs.u_dirichlet
    w = state.cell_integration.rid_to_weights
    d = state.cell_integration.d
    ncells = length(cell_to_nodes)
    p = state.p
    xe = state.xe_setup.xe
    ue = state.ue_alloc.ue
    prl_level = state.prl_level
    mem_layout = state.mem_layout

    ∇ste = map(m->collect(permutedims(m)),∇s)[1]
    we = w[1]
    nl = length(cell_to_nodes[1]) 
    nq = length(we) # number of quadrature points.
    Tx = eltype(xe)
    TJ = typeof(zero(Tx)*zero(Tx)')

    Jt = zero(TJ)
    ∇u = zero(Tx)

    fill!(V_coo,zero(Float64)) # Fill V_coo with Float64 zeros to then later accumulate to it
    setup_ue!(ue, cell_to_dofs, u_dirichlet, u_dofs, ncells, nl, mem_layout) # update ue

    # Almost all time is spent inside these kernels (as expected).
    if prl_level == :cell
        for cell in 1:ncells # Call the generic kernel function
            kernel_generic!(V_coo, cell, nq, nl, Jt, xe, ∇ste, we, ∇u, ue, dflux, p)
        end
    elseif prl_level == :coalesce
        for cell in 1:ncells
            kernel_coalesced!(V_coo, cell, nq, nl, Jt, xe, ∇ste, we, ∇u, ue, dflux, p, ncells)
        end
    elseif prl_level == :elem_j
        threads = ncells * nl
        for thread in 1:threads
            kernel_elem_j!(V_coo, thread, nq, nl, Jt, xe, ∇ste, we, ∇u, ue, dflux, p)
        end
    elseif prl_level == :elem_ij
        threads = ncells * nl^2
        for thread in 1:threads
            kernel_elem_ij!(V_coo, thread, nq, nl, Jt, xe, ∇ste, we, ∇u, ue, dflux, p)
        end
    elseif prl_level == :quad
        threads = ncells * nq
        for thread in 1:threads
            kernel_quad!(V_coo, thread, nq, nl, Jt, xe, ∇ste, we, ∇u, ue, dflux, p, ncells)
        end
    elseif prl_level == :full
        threads = ncells * nq * nl^2
        for thread in 1:threads
            kernel_full_prl!(V_coo, thread, nq, nl, Jt, xe, ∇ste, we, ∇u, ue, dflux, p, ncells)
        end
    end
end

function jacobian_cells_original!(V_coo,u_dofs,state)

    dflux = state.dflux
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
                    Ae[i,j] += ( ∇dv⋅dflux(∇u,∇du,p) )*dV
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


function jacobian_cells_gpu!(V_coo,u_dofs,state)
    dflux = state.dflux
    cell_to_dofs = state.dofs.cell_to_dofs
    cell_to_nodes = state.cell_isomap.face_to_nodes
    ∇s = state.cell_isomap.rid_to_shape_grads
    s = state.cell_isomap.rid_to_shape_vals
    u_dirichlet = state.dirichlet_bcs.u_dirichlet
    w = state.cell_integration.rid_to_weights
    d = state.cell_integration.d
    ncells = length(cell_to_nodes)
    p = state.p
    xe = state.xe.node_to_x
    ue = state.ue_alloc.ue
    prl_level = state.prl_level
    
    ∇ste = map(m->collect(permutedims(m)),∇s)[1]
    we = w[1]
    nl = length(cell_to_nodes[1]) 
    nq = length(we)
    Tx = eltype(xe)
    TJ = typeof(zero(Tx)*zero(Tx)')

    Jt = zero(TJ)
    ∇u = zero(Tx)
    # This also needs to be done on the gpu.
    fill!(V_coo,zero(Float64)) # Move this such that each thread does it individually, check how that goes for the kernels i and ij
    setup_ue!(ue, cell_to_dofs, u_dirichlet, u_dofs, ncells, nl, mem_layout) # update ue

    # The _d denotes that it's an array on the device.
    V_coo_d,∇ste_d,w_d,xe_d,ue_d,Jt_d,∇u_d,ncells_d,nl_d,nq_d,p_d = cpu_gpu_transfer(V_coo, ∇ste, 
            we, xe, Jt, ∇u, ue, ncells, nq, nl, p) # This should be timed

    # Maybe set all of this in a separate function with config and kernel launch.
    # Here the configuration is set for the kernel (threads, blocks etc.)
    ckernel=@cuda launch=false assemble_cell_gpu(V_coo_d,∇ste_d,w_d,xe_d,ue_d,Jt_d,∇u_d,ncells_d,p_d,dflux,nl_d,nq_d)
    config = launch_configuration(ckernel.fun)
    threads = min(ncells, config.threads)
    blocks =  cld(ncells, threads)
    
    CUDA.@sync @cuda threads=threads blocks=blocks assemble_cell_gpu(V_coo_d,∇ste_d,w_d,xe_d,ue_d,Jt_d,∇u_d,
                                    ncells_d,p_d,dflux,nl_d,nq_d) # This should be timed.
    # Then you need to copy back V_coo back to cpu memory
    copyto!(V_coo, V_coo_d) # This copies over to already existing memory.
end

function assemble_cell_gpu(V_coo_d,∇ste_d,w_d,xe_d,ue_d,Jt_d,∇u_d,ncells_d,p_d,dflux,nl_d,nq_d)
    idx = threadIdx().x + (blockIdx().x - Int32(1)) * blockDim().x
    if idx <= ncells_d
        kernel_generic!(V_coo_d, idx, nq_d, nl_d, Jt_d, xe_d, ∇ste_d, w_d, ∇u_d, ue_d, dflux, p_d)
    end
    nothing
end

function assemble_cell_gpu_coalesced(V_coo_d,∇ste_d,w_d,xe_d,ue_d,Jt_d,∇u_d,ncells_d,p_d,dflux,nl_d,nq_d)
    idx = threadIdx().x + (blockIdx().x - Int32(1)) * blockDim().x
    if idx <= ncells_d
        kernel_coalesced!(V_coo_d, idx, nq_d, nl_d, Jt_d, xe_d, ∇ste_d, w_d, ∇u_d, ue_d, dflux, p_d, ncells_d)
    end
    nothing
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

function solve_problem(params,state,results)
    timer = state.timer
    problem = nonlinear_problem(state)
    solver = params[:solver]
    x = problem.initial()
    @timeit timer "solver.setup" setup = solver.setup(x,problem,params)
    @timeit timer "solver.solve!" iterations = solver.solve!(x,setup)
    @timeit timer "solver.finalize!" solver.finalize!(setup)
    results[:iterations] = iterations
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

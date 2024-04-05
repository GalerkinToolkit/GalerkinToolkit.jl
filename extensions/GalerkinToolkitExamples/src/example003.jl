module Example003

import GalerkinToolkit as gk
import ForwardDiff
using StaticArrays
using LinearAlgebra
using SparseArrays
using WriteVTK
using PartitionedArrays: JaggedArray, sparse_matrix, sparse_matrix!
using NLsolve
using Random
using GalerkinToolkitExamples: Example001

using Preconditioners
using IterativeSolvers: cg!

using CUDA
# This one implements a vanilla sequential iso-parametric p-Laplacian solver by only
# using the mesh interface.

global jacobian_time = 0.0
global gpu_transfer_time = 0.0

function main(params_in)
    # Process params
    global jacobian_time = 0.0 # set again to zero (for between simulations)
    global gpu_transfer_time = 0.0

    params_default = default_params()
    params = add_default_params(params_in,params_default)
    results = Dict{Symbol,Any}()

    # Setup main data structures
    state, gpu_setup_time = setup(params)
    add_basic_info(results,params,state)

    x = solve_problem(params,state,results)

    # Post process
    uh = setup_uh(x,state)
    integrate_error_norms(results,uh,state)
    export_results(uh,params,state)

    results, jacobian_time, gpu_transfer_time, gpu_setup_time, x
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
    params[:jacobian_implementation] = :original
    params[:float_type] = Dict(:Float => Float64, :Int => Int32)
    params
end

function nlsolve_solver(;linear_solver=Example001.cg_amg_solver(),options...)
    function setup(x0,nlp,params)
        r0,J0,cache = nlp.linearize(x0, params)
        dx = similar(r0,axes(J0,2)) # TODO is there any way of reusing this in the nonlinear solve?
        ls_setup = linear_solver.setup(dx,J0,r0)
        function linsolve(x,A,b)
            # TODO we dont need to re-setup for the first
            # linear solve.
            # This can be avoided with a Ref{Bool} shared
            # between this function and j!
            linear_solver.setup!(ls_setup,A)
            linear_solver.solve!(x,ls_setup,b)
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
        result.f_calls, x
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
    u_dirichlet = u_dirichlet = convert(Vector{params[:float_type][:Float]}, u.(x_dirichlet))
    dirichlet_bcs = (;free_and_dirichlet_nodes,u_dirichlet,node_to_tag)
end

function setup_integration(params,objects,type_setup)
    mesh = params[:mesh]
    degree = params[:integration_degree]
    D = gk.num_dims(mesh)
    float_type = params[:float_type][:Float]
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
    if type_setup == "cell"
        rid_to_weights = convert(Vector{Vector{float_type}}, map(gk.weights,integration_rules))
    elseif type_setup == "face"
        rid_to_weights = convert(Tuple{Vector{float_type}, Vector{float_type}}, map(gk.weights,integration_rules))
    end
    rid_to_coords = map(gk.coordinates,integration_rules)
    integration = (;rid_to_weights,rid_to_coords,face_to_rid,d)
end

function setup_isomap(params,integration, type_setup)
    rid_to_coords = integration.rid_to_coords
    face_to_rid = integration.face_to_rid
    d = Int32(integration.d)
    float_type = params[:float_type][:Float]
    mesh = params[:mesh]
    ref_cells = gk.reference_faces(mesh,d)
    shape_funs = map(rid_to_coords,ref_cells) do q,ref_cell
        shape_vals = gk.tabulator(ref_cell)(gk.value,q)
        shape_grads = gk.tabulator(ref_cell)(ForwardDiff.gradient,q)
        shape_vals, shape_grads
    end
    if type_setup == "cell"
        rid_to_shape_vals = convert(Vector{Matrix{float_type}}, map(first,shape_funs))
    elseif type_setup == "face"
        rid_to_shape_vals = convert(Tuple{Matrix{float_type}, Matrix{float_type}}, map(first,shape_funs))
    end
    rid_to_shape_grads = map(last,shape_funs)
    rid_to_shape_grads = [[float_type.(float_type.(vec(m))) for m in rid_to_shape_grads[1]]]
    face_to_nodes = convert(JaggedArray{Int32, Int32}, JaggedArray(gk.face_nodes(mesh,d)))
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
    cell_to_dofs = convert(JaggedArray{Int32,Int32}, JaggedArray(cell_to_dofs_data,cell_to_nodes.ptrs))
    face_to_dofs = JaggedArray(face_to_dofs_data,face_to_nodes.ptrs)
    dofs = (;cell_to_dofs,face_to_dofs,n_dofs)
    dofs
end

function setup_xe(cell_isomap, cell_integration, params)
    cell_to_nodes = cell_isomap.face_to_nodes
    node_to_coords = cell_isomap.node_to_coords
    jacobian = params[:jacobian_implementation]
    ncells = length(cell_to_nodes)
    nl = length(cell_to_nodes[1])
    d = length(eltype(node_to_coords))
    # dof as the inner loop, since dofs per cell are looped over in assembly.
    xe = Vector{SVector{d, params[:float_type][:Float]}}(undef, ncells * nl) 
    if jacobian == :cpu_v1 || jacobian == :gpu_v1 || jacobian == :original
        for cell in 1:ncells
            for dof in 1:nl
                index = (cell - 1) * nl + dof
                node = cell_to_nodes[cell][dof]
                xe[index] = node_to_coords[node]
            end
        end
    else # coalesced for gpu memory
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

function setup_ue!(ue, cell_to_dofs, u_dirichlet, u_dofs, ncells, nl, jacobian)
    # Here you also want jacobian implementation to determine.
    if jacobian == :cpu_v1 || jacobian == :gpu_v1 || jacobian == :original
        index_function = (cell, dof, nl, ncells) -> (cell - 1) * nl + dof
    else # coalesced for gpu access.
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

function ue_allocation(cell_isomap, params)
    cell_to_nodes = cell_isomap.face_to_nodes
    ncells = length(cell_to_nodes) # number of elements
    nl = length(cell_to_nodes[1]) # number of nodes per element

    # Allocate memory for the element u
    ue = Vector{params[:float_type][:Float]}(undef, ncells * nl)
    ue_alloc = (;ue)
end

function setup_Jt(cell_integration, xe_setup, cell_isomap, float_type, jacobian)
    # Precompute the Jt values. ncells * nq of static matrices  
    w = cell_integration.rid_to_weights
    xe = xe_setup.xe
    cell_to_nodes = cell_isomap.face_to_nodes
    ∇s = cell_isomap.rid_to_shape_grads

    ncells = Int32(length(cell_to_nodes))
    we = float_type[:Float].(w[1])
    nq = Int32(length(we))
    nl = Int32(length(cell_to_nodes[1]))
    Tx = eltype(xe)
    TJ = typeof(zero(Tx)*zero(Tx)')
    ∇ste = map(m->collect(permutedims(m)),∇s)[1]

    zero_Jt = zero(TJ)
    Jt_local = zero_Jt
    Jt = Matrix{typeof(zero_Jt)}(undef, (ncells,nq))

    # Here you need to consider the memory layout. 
    if jacobian == :cpu_v1 || jacobian == :gpu_v1 || jacobian == :original
        for cell in 1:ncells
            for iq in 1:nq 
                Jt_local = zero_Jt
                for k in 1:nl
                    k_index = (cell - 1) * nl + k
                    x = xe[k_index]
                    ∇sqx = ∇ste[k,iq] 
                    Jt_local += ∇sqx*x'
                end
                Jt[cell, iq] = Jt_local
            end
        end
    else
        for cell in 1:ncells
            for iq in 1:nq 
                Jt_local = zero_Jt
                for k in 1:nl
                    k_index = (k - 1) * ncells + cell
                    x = xe[k_index]
                    ∇sqx = ∇ste[k,iq] 
                    Jt_local += ∇sqx*x'
                end
                Jt[cell, iq] = Jt_local
            end
        end
    end
    Jt_precompute = (;Jt,zero_Jt)
end

function cpu_gpu_transfer(params,cell_isomap,dofs,cell_integration,p,xe_setup,ue_alloc,Jt_precompute)
    cell_to_nodes = cell_isomap.face_to_nodes
    cell_to_dofs = dofs.cell_to_dofs
    dofs = cell_to_dofs[1]
    ∇s = cell_isomap.rid_to_shape_grads
    we = params[:float_type][:Float].(cell_integration.rid_to_weights[1])
    d = cell_integration.d
    ncells = Int32(length(cell_to_nodes))
    xe = xe_setup.xe
    ue = ue_alloc.ue
    Jt = Jt_precompute.Jt
    float_type = params[:float_type][:Float]
    ∇ste = map(m->collect(permutedims(m)),∇s)[1]
    nl = Int32(length(cell_to_nodes[1]))
    nq = Int32(length(we))
    Tx = eltype(xe)
    TJ = typeof(zero(Tx)*zero(Tx)')
    ndofs = length(dofs)

    zero_Jt = zero(TJ)
    ∇u = zero(Tx)
    n_coo = ncells * ndofs^2
    # cu() copies default to 32 float_type
    setup_gpu_time  = CUDA.@elapsed CUDA.@sync begin
        V_coo_d = CuVector{float_type}(undef, n_coo)
        ∇ste_d = cu(cu.(∇ste))
        w_d = CuArray{float_type}(we)
        xe_d = CuArray(xe) 
        Jt_d = CuArray(Jt) 
        ∇u_d = cu(cu.(∇u)) 
        ue_d = CuArray{float_type}(ue) 
        ncells_d = cu(ncells)
        nq_d = cu(nq)
        nl_d = cu(nl)
        p_d = cu(Int32(p))
    end
    [V_coo_d,∇ste_d,w_d,xe_d,ue_d,Jt_d,∇u_d,ncells_d,nl_d,nq_d,p_d], setup_gpu_time
end

function setup(params)
    dirichlet_bcs = setup_dirichlet_bcs(params)
    neumann_bcs = setup_neumann_bcs(params)
    cell_integration = setup_integration(params,:cells,"cell")
    face_integration = setup_integration(params,:faces,"face")
    cell_isomap = setup_isomap(params,cell_integration,"cell")
    face_isomap = setup_isomap(params,face_integration,"face")
    dofs = setup_dofs(params,dirichlet_bcs,cell_isomap,face_isomap)
    user_funs = setup_user_funs(params)
    ue_alloc = ue_allocation(cell_isomap, params)
    xe_setup = setup_xe(cell_isomap, cell_integration, params)
    Jt_precompute = setup_Jt(cell_integration,xe_setup,cell_isomap,params[:float_type],params[:jacobian_implementation])

    solver = params[:solver]
    p = Int32(params[:p])
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

    if occursin("original", string(params[:jacobian_implementation]))
        jacobian_cells! = jacobian_cells_original!
    elseif occursin("cpu", string(params[:jacobian_implementation]))
        jacobian_cells! = jacobian_cells_cpu_extension!
    elseif occursin("gpu", string(params[:jacobian_implementation]))
        jacobian_cells! = jacobian_cells_gpu!
        # Allocate memory on gpu
        gpu_pointers, gpu_setup_time = cpu_gpu_transfer(params,cell_isomap,dofs,cell_integration,p,xe_setup,ue_alloc,Jt_precompute)
    else
        error("unspecified jacobian function")
    end

    jacobian_impl = params[:jacobian_implementation]
    float_type = params[:float_type]

    if occursin("gpu", string(params[:jacobian_implementation]))
        return state = (;p,flux,dflux,solver,dirichlet_bcs,neumann_bcs,cell_integration,cell_isomap,
                    face_integration,face_isomap,user_funs,dofs,xe_setup,ue_alloc,jacobian_cells!,
                    gpu_pointers,float_type,jacobian_impl,Jt_precompute), gpu_setup_time
    else
        return state = (;p,flux,dflux,solver,dirichlet_bcs,neumann_bcs,cell_integration,cell_isomap,
                    face_integration,face_isomap,user_funs,dofs,xe_setup,ue_alloc,jacobian_cells!,float_type,
                    jacobian_impl,Jt_precompute), 0
    end
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
    function initial()
        n = state.dofs.n_dofs
        Random.seed!(1)
        rand(Float64,n)
    end
    function linearize(u_dofs, params)
        r,J,V,K = assemble_sybmolic(state, params)
        cache = (V,K)
        residual_cells!(r,u_dofs,state)
        residual_faces!(r,u_dofs,state)
        state.jacobian_cells!(V,u_dofs,state) 
        sparse_matrix!(J,V,K)
        r,J,cache
    end
    function residual!(r,u_dofs,cache)
        fill!(r,0)
        residual_cells!(r,u_dofs,state)
        residual_faces!(r,u_dofs,state)
        r
    end
    function jacobian!(J,u_dofs,cache)
        (V,K) = cache 
        jac_time_iter = @elapsed gpu_data_transfer = state.jacobian_cells!(V,u_dofs,state)
        global jacobian_time += jac_time_iter
        global gpu_transfer_time += gpu_data_transfer
        sparse_matrix!(J,V,K)
        J
    end
    (;initial,linearize,residual!,jacobian!)
end

function assemble_sybmolic(state,params)
    cell_to_dofs = state.dofs.cell_to_dofs
    n_dofs = state.dofs.n_dofs
    jacobian = params[:jacobian_implementation]

    n_coo = 0
    ncells = length(cell_to_dofs)
    for cell in 1:ncells
        dofs = cell_to_dofs[cell]
        ndofs = length(dofs)
        n_coo += ndofs*ndofs
    end

    I_coo = Vector{Int32}(undef,n_coo)
    J_coo = Vector{Int32}(undef,n_coo)
    V_coo = Vector{params[:float_type][:Float]}(undef,n_coo)
    r = zeros(Float64,state.dofs.n_dofs)

    n_coo = 0
    if jacobian == :cpu_v1 || jacobian == :gpu_v1 || jacobian == :original
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
    else # coalesced for gpu memory
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

@inline function kernel_generic!(V_coo, cell, nq, nl, Jt, xe, ∇ste, we, ∇u, ue, dflux, p)
    # Integrate
    for iq in 1:nq 
        Jt_local = Jt[cell,iq]
        ∇u_local = ∇u

        detJt = det(Jt_local)
        invJt = inv(Jt_local)
        dV = abs(detJt)*we[iq] 

        for k in 1:nl # nq * ncells threads
            k_index = (cell - 1) * nl + k
            uek = ue[k_index]
            ∇xek = invJt*∇ste[k,iq] 
            ∇u_local += ∇xek*uek
        end

        for j in 1:nl # ncells * nl^2 threads.
            ∇du = invJt*∇ste[j,iq]
            for i in 1:nl
                ∇dv = invJt*∇ste[i,iq]
                i_coo = (cell - 1) * nl^2 + ((j - 1) * nl + i) 
                V_coo[i_coo] += (∇dv⋅dflux(∇u_local,∇du,p)) * dV
            end
        end
    end
end

@inline function kernel_coalesced!(V_coo, cell, nq, nl, Jt, xe, ∇ste, we, ∇u, ue, dflux, p, ncells)
    for iq in 1:nq 
        ∇u_local = ∇u
        Jt_local = Jt[cell,iq]

        detJt = det(Jt_local)
        invJt = inv(Jt_local)
        dV = abs(detJt)*we[iq]

        for k in 1:nl
            k_index = (k - 1) * ncells + cell
            uek = ue[k_index]
            ∇xek = invJt*∇ste[k,iq]
            ∇u_local += ∇xek*uek
        end

        for j in 1:nl
            ∇du = invJt*∇ste[j,iq]
            for i in 1:nl
                ∇dv = invJt*∇ste[i,iq]
                i_coo = cell + (j-1) * ncells * nl + (i-1) * ncells
                V_coo[i_coo] += (∇dv⋅dflux(∇u_local,∇du,p)) * dV
            end
        end
    end
end

@inline function kernel_generic_v3!(V_coo, cell, i, j, nq, nl, Jt, xe, ∇ste, we, ∇u, ue, dflux, p, ncells)
    V_coo_local = 0
    for iq in 1:nq
        ∇u_local = ∇u
        Jte = Jt[cell,iq]

        detJt = det(Jte)
        invJt = inv(Jte)
        dV = abs(detJt)*we[iq] 

        for k in 1:nl
            k_index = (k - 1) * ncells + cell
            uek = ue[k_index]
            ∇xek = invJt*∇ste[k,iq] 
            ∇u_local += ∇xek*uek
        end

        ∇du = invJt*∇ste[j,iq]
        ∇dv = invJt*∇ste[i,iq]
        V_coo_local += (∇dv⋅dflux(∇u_local,∇du,p)) * dV
    end
    i_coo = cell + (j-1) * ncells * nl + (i-1) * ncells
    V_coo[i_coo] = V_coo_local
end

@inline function kernel_∇ue_tiled!(V_coo, idx, cell, i, j, nq, nl, Jt, ∇ste, we, ∇u, ue, dflux, p)
    # In this kernel you have ncells * nl^2 threads
    # Use here the double modulus counter to get a loop over Jt and we for each cell and iq. Can you assume that nl^2 > nq?
    # First define the data structures in local and shared memory
    #TJt = CuStaticSharedArray(typeof(Jt[1,1]), 1024)
    Twe = CuStaticSharedArray(Float32, 1024)
    Tue = CuStaticSharedArray(Float32, 1024)
    # identify your position in the blocks
    tid = threadIdx().x
    node_index = (cell - 1) * nl + i
    # Then collectively load it into shared memory
    #TJt[tid] = Jt[cell, i] # load the i-th (or iq-th) value.
    Twe[tid] = we[(tid-1)%nq+1] # change this to the double modulus counter. Then you don't need to assume nq == nl
    Tue[tid] = ue[node_index]
    sync_threads()

    V_coo_local = 0
    for iq in 1:nq
        ∇u_local = ∇u
        Jte = Jt[cell,iq]

        detJt = det(Jte)
        invJt = inv(Jte)
        dV = abs(detJt) * Twe[iq] 

        for k in 1:nl
            k_index = Int32(floor(idx / nl^2) * nl^2) + k 
            uek = Tue[k_index]
            ∇xek = invJt*∇ste[k,iq] 
            ∇u_local += ∇xek*uek
        end

        ∇du = invJt*∇ste[j,iq]
        ∇dv = invJt*∇ste[i,iq]
        V_coo_local += (∇dv⋅dflux(∇u_local,∇du,p)) * dV
    end
    i_coo = (cell - 1) * nl^2 + ((j - 1) * nl + i) 
    V_coo[i_coo] = V_coo_local
end

# Maybe the logic should be to get xe + ue, then loop through the iq in parallel. Then the last loop is nl*nq.
# In this case you have ncells * nl threads that are launched. In this version the ∇ste will never be a 
# problem because threads in a block is approx. 1024 while ∇ste << 1024. 
# So there are two parts to this -> loading values which is according to nl and computing
# which is done according to nq.
@inline function kernel_tiled!(V_coo, tid, nq, nl, Jt, xe, ∇ste, we, ∇u, ue, dflux, p, ncells)
    # variables to load into shared memory: Jt, ∇u, ∇ste, xe, we, ue.
    # Think about this as how many values does each thread store (max 12 for 1024 threads per block).
    # Specify the shared memory per threadblock.
    # Assume (for now) nq == nl and that d = 2.
    # So the gain comes from not having to recompute invJt etc. for all iq points.
    # Note here that each shmem array should be 32 or a multiple of it (related to warp execution).
    max_threads = 1024 # This one is then defined outside the kernel as a constant. 
    Txe = CuStaticSharedArray(Float32, (2,max_threads)) # Two per thread (for the coordinates).
    Tue = CuStaticSharedArray(Float32, max_threads) # One per thread
    T∇ste = CuStaticSharedArray(Float32, (max_threads,2)) # This one is flattened using 2 indices (nq and nl)
    T∇u = CuStaticSharedArray(Float32, (2,max_threads))
    T∇dux = CuStaticSharedArray(Float32, (2,max_threads))
    TdV = CuStaticSharedArray(Float32, max_threads)
    # Maybe use static arrays here. 
    # identify your position in nl block, nq and cell. tid is the grid value.
    cell = Int((ceil(tid/nl^2)-1) % ncells)+1 
    j = Int(((ceil(tid/nl)-1) % nl)+1) # Column 
    i = ((tid-1)%nl)+1 # Row
    iq = ((tid-1)%nq)+1 # Identical to the i index.
    node_index = (cell - 1) * nl + i # Here you get node to use in ue/xe
    xe_col = tid % 2 + 1 # Can be used for th xe to alternatve columns of which value to store.     

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

function assemble_v3(V_coo, nq, nl, Jt, xe, ∇ste, we, ∇u, ue, dflux, p, ncells)
    idx = threadIdx().x + (blockIdx().x - Int32(1)) * blockDim().x
    cell = Int((ceil(idx/nl^2)-1) % ncells)+1 
    i = ((idx-1)%nl)+1
    j = Int(((ceil(idx/nl)-1) % nl)+1) 
    if idx <= ncells * nl^2
        kernel_generic_v3!(V_coo, cell, i, j, nq, nl, Jt, xe, ∇ste, we, ∇u, ue, dflux, p, ncells)
    end
    nothing
end

function assemble_kernel_∇ue_tiled(V_coo, nq, nl, Jt, xe, ∇ste, we, ∇u, ue, dflux, p, ncells)
    idx = threadIdx().x + (blockIdx().x - Int32(1)) * blockDim().x
    cell = Int((ceil(idx/nl^2)-1) % ncells)+1 
    i = ((idx-1)%nl)+1
    j = Int(((ceil(idx/nl)-1) % nl)+1)
    if idx <= ncells * nl^2 
        kernel_∇ue_tiled!(V_coo, idx, cell, i, j, nq, nl, Jt, ∇ste, we, ∇u, ue, dflux, p)
    end
    nothing
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

function jacobian_cells_cpu_extension!(V_coo,u_dofs,state)
    dflux = state.dflux
    cell_to_dofs = state.dofs.cell_to_dofs
    cell_to_nodes = state.cell_isomap.face_to_nodes
    ∇s = state.cell_isomap.rid_to_shape_grads
    s = state.cell_isomap.rid_to_shape_vals
    u_dirichlet = state.dirichlet_bcs.u_dirichlet
    w = state.cell_integration.rid_to_weights
    ncells = length(cell_to_nodes)
    p = state.p
    xe = state.xe_setup.xe
    ue = state.ue_alloc.ue
    jacobian = state.jacobian_impl
    float_type = state.float_type
    Jt = state.Jt_precompute.Jt
    zero_Jt = state.Jt_precompute.zero_Jt

    ∇ste = map(m->collect(permutedims(m)),∇s)[1]
    we = float_type[:Float].(w[1])
    nl = Int32(length(cell_to_nodes[1])) 
    nq = Int32(length(we)) 
    Tx = eltype(xe)
    ∇u = zero(Tx)

    fill!(V_coo, 0.0) # reset the jacobian
    setup_ue!(ue, cell_to_dofs, u_dirichlet, u_dofs, ncells, nl, jacobian) # update ue

    if jacobian == :cpu_v1
        for cell in 1:ncells
            kernel_generic!(V_coo, cell, nq, nl, Jt, xe, ∇ste, we, ∇u, ue, dflux, p)
        end
    elseif jacobian == :cpu_v2
        for cell in 1:ncells
            kernel_coalesced!(V_coo, cell, nq, nl, Jt, xe, ∇ste, we, ∇u, ue, dflux, p, ncells)
        end
    elseif jacobian == :cpu_v3
        threads = ncells * nl^2
        for thread in 1:threads
            cell = Int((ceil(thread/nl^2)-1) % ncells)+1
            i = ((thread-1)%nl)+1
            j = Int(((ceil(thread/nl)-1) % nl)+1)
            kernel_generic_v3!(V_coo, cell, i, j, nq, nl, Jt, xe, ∇ste, we, ∇u, ue, dflux, p, ncells)
        end
    end
    gpu_transfer = 0
    gpu_transfer
end

function jacobian_cells_gpu!(V_coo,u_dofs,state)
    dflux = state.dflux
    cell_to_dofs = state.dofs.cell_to_dofs
    cell_to_nodes = state.cell_isomap.face_to_nodes
    ncells = length(cell_to_nodes)
    nl = length(cell_to_nodes[1]) 
    u_dirichlet = state.dirichlet_bcs.u_dirichlet
    xe = state.xe_setup.xe
    ue = state.ue_alloc.ue
    jacobian = state.jacobian_impl

    V_coo_d,∇ste_d,w_d,xe_d,ue_d,Jt_d,∇u_d,ncells_d,nl_d,nq_d,p_d = state.gpu_pointers

    setup_ue!(ue, cell_to_dofs, u_dirichlet, u_dofs, ncells, nl, jacobian) # update ue
    gpu_transfer1 = CUDA.@elapsed CUDA.@sync begin
        copyto!(ue_d, ue)
        fill!(V_coo_d,0.0)
    end
    # max_threads = attribute(CuDevice(0),CUDA.DEVICE_ATTRIBUTE_MAX_THREADS_PER_BLOCK)

    if jacobian == :gpu_v1 # baseline
        # Here the configuration is set for the kernel (threads, blocks etc.)
        ckernel=@cuda launch=false assemble_cell_gpu(V_coo_d,∇ste_d,w_d,xe_d,ue_d,Jt_d,∇u_d,ncells_d,p_d,dflux,nl_d,nq_d)
        config = launch_configuration(ckernel.fun)
        threads = min(ncells, config.threads)
        blocks =  cld(ncells, threads)
        CUDA.@sync @cuda threads=threads blocks=blocks assemble_cell_gpu(V_coo_d,∇ste_d,w_d,xe_d,ue_d,Jt_d,∇u_d,
                                        ncells_d,p_d,dflux,nl_d,nq_d)
    elseif jacobian == :gpu_v2 # coalesced
        ckernel=@cuda launch=false assemble_cell_gpu_coalesced(V_coo_d,∇ste_d,w_d,xe_d,ue_d,Jt_d,∇u_d,ncells_d,p_d,dflux,nl_d,nq_d)
        config = launch_configuration(ckernel.fun)
        threads = min(ncells, config.threads)
        blocks =  cld(ncells, threads)
        
        CUDA.@sync @cuda threads=threads blocks=blocks assemble_cell_gpu_coalesced(V_coo_d,∇ste_d,w_d,xe_d,ue_d,Jt_d,∇u_d,
                                        ncells_d,p_d,dflux,nl_d,nq_d)
    
    elseif jacobian == :gpu_v3 # version 3 ncells * nl^2
        ckernel=@cuda launch=false assemble_v3(V_coo_d,nq_d,nl_d,Jt_d,xe_d,∇ste_d,w_d,∇u_d,ue_d,dflux,p_d,ncells_d)
        config = launch_configuration(ckernel.fun)
        threads = min(ncells * nl^2, config.threads)
        blocks =  cld(ncells * nl^2, threads)
        CUDA.@sync @cuda threads=threads blocks=blocks assemble_v3(V_coo_d,nq_d,nl_d,Jt_d,xe_d,∇ste_d,w_d,∇u_d,ue_d,dflux,p_d,ncells_d)
        
    elseif jacobian == :gpu_v4
        # ncells * nq^2 threads
        threads = fld(512, nl^2) * nl^2
        blocks = Int(cld((ncells*nl^2)/nl, threads))
        CUDA.@sync @cuda threads=threads blocks=blocks assemble_kernel_∇ue_tiled(V_coo_d, nq_d, nl_d, Jt_d, xe_d, ∇ste_d, w_d, ∇u_d, ue_d, dflux, p_d, ncells_d)
    end

    # Then you need to copy back V_coo back to cpu memory
    gpu_transfer2 = CUDA.@elapsed CUDA.@sync begin
        copyto!(V_coo, V_coo_d) # This copies over to already existing memory.
    end

    gpu_transfer1 + gpu_transfer2
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
    problem = nonlinear_problem(state)
    solver = params[:solver]
    x = problem.initial()
    setup = solver.setup(x,problem,params)
    iterations,x = solver.solve!(x,setup)
    solver.finalize!(setup)
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

module Poisson

import GalerkinToolkit as glk
import ForwardDiff
using StaticArrays
using LinearAlgebra
using SparseArrays
using WriteVTK

function main(params_in)

    # Process params
    params_default = default_params_poisson()
    params = add_default_params(params_in,params_default)

    # Setup main data structures
    state = setup(params)

    # Assemble system and solve it
    A,b = assemble_system(state)
    add_neumann_bcs!(b,state)
    uh = solve_system(A,b,state)

    # Post process
    results = Dict{Symbol,Any}()
    add_basic_info!(results,state)
    integrate_error_norms!(results,uh,state)
    export_results(uh,params,state)

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

function default_params_poisson()
    outdir = mkpath(joinpath(@__DIR__,"..","output"))
    params = Dict{Symbol,Any}()
    params[:mesh] = glk.cartesian_mesh((0,1,0,1),(10,10))
    params[:u] = (x) -> sum(x)
    params[:f] = (x) -> 0.0
    params[:g] = (x) -> 0.0
    params[:dirichlet_tags] = ["boundary"]
    params[:neumann_tags] = String[]
    params[:integration_degree] = 2
    params[:solve] = \
    params[:example_path] = joinpath(outdir,"poisson")
    params
end

function setup_dirichlet_bcs(params)
    mesh = params[:mesh]
    u = params[:u]
    dirichlet_tags = params[:dirichlet_tags]
    node_to_tag = zeros(glk.num_nodes(mesh))
    tag_to_name = dirichlet_tags
    glk.classify_mesh_nodes!(node_to_tag,mesh,tag_to_name)
    free_and_dirichlet_nodes = glk.partition_from_mask(i->i==0,node_to_tag)
    node_to_x = glk.node_coordinates(mesh)
    dirichlet_nodes = last(free_and_dirichlet_nodes)
    x_dirichlet = view(node_to_x,dirichlet_nodes)
    u_dirichlet = u.(x_dirichlet)
    dirichlet_bcs = (;free_and_dirichlet_nodes,u_dirichlet,node_to_tag)
end

function setup_integration(params,objects)
    mesh = params[:mesh]
    degree = params[:integration_degree]
    D = glk.num_dims(mesh)
    if objects === :cells
        d = D
    elseif objects === :faces
        d = D-1
    else
        error("")
    end
    # Integration
    ref_cells = glk.reference_faces(mesh,d)
    face_to_rid = glk.face_reference_id(mesh,d)
    integration_rules = map(ref_cells) do ref_cell
        glk.quadrature(glk.geometry(ref_cell),degree)
    end
    rid_to_weights = map(glk.weights,integration_rules)
    rid_to_coords = map(glk.coordinates,integration_rules)
    integration = (;rid_to_weights,rid_to_coords,face_to_rid,d)
end

function setup_isomap(params,integration)
    rid_to_coords = integration.rid_to_coords
    face_to_rid = integration.face_to_rid
    d = integration.d
    mesh = params[:mesh]
    ref_cells = glk.reference_faces(mesh,d)
    shape_funs = map(rid_to_coords,ref_cells) do q,ref_cell
        shape_functions = glk.shape_functions(ref_cell)
        shape_vals = glk.tabulation_matrix(shape_functions)(glk.value,q)
        shape_grads = glk.tabulation_matrix(shape_functions)(ForwardDiff.gradient,q)
        shape_vals, shape_grads
    end
    rid_to_shape_vals = map(first,shape_funs)
    rid_to_shape_grads = map(last,shape_funs)
    face_to_nodes = glk.face_nodes(mesh,d)
    node_to_coords = glk.node_coordinates(mesh)
    isomap = (;face_to_nodes,node_to_coords,face_to_rid,rid_to_shape_vals,rid_to_shape_grads,d)
end

function setup_neumann_bcs(params)
    mesh = params[:mesh]
    neumann_tags = params[:neumann_tags]
    neum_face_to_face = Int[]
    if length(neumann_tags) != 0
        D = glk.num_dims(mesh)
        tag_to_groups = glk.physical_groups(mesh,D-1)
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

function setup(params)
    dirichlet_bcs = setup_dirichlet_bcs(params)
    neumann_bcs = setup_neumann_bcs(params)
    cell_integration = setup_integration(params,:cells)
    face_integration = setup_integration(params,:faces)
    cell_isomap = setup_isomap(params,cell_integration)
    face_isomap = setup_isomap(params,face_integration)
    user_funs = setup_user_funs(params)
    solve = params[:solve]
    state = (;solve,dirichlet_bcs,neumann_bcs,cell_integration,cell_isomap,face_integration,face_isomap,user_funs)
end

function assemble_system(state)

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

    # Count coo entries
    n_coo = 0
    ncells = length(cell_to_nodes)
    node_to_free_node = glk.permutation(free_and_dirichlet_nodes)
    nfree = length(first(free_and_dirichlet_nodes))
    for cell in 1:ncells
        nodes = cell_to_nodes[cell]
        for nj in nodes
            if node_to_free_node[nj] > nfree
                continue
            end
            for ni in nodes
                if node_to_free_node[ni] > nfree
                    continue
                end
                n_coo += 1
            end
        end
    end

    # Allocate coo values
    I_coo = zeros(Int32,n_coo)
    J_coo = zeros(Int32,n_coo)
    V_coo = zeros(Float64,n_coo)
    b = zeros(nfree)

    # Allocate auxiliary buffers
    ∇x = map(i->similar(i,size(i,2)),∇s)
    Aes = map(i->zeros(size(i,2),size(i,2)),∇s)
    fes = map(i->zeros(size(i,2)),∇s)
    ues = map(i->zeros(size(i,2)),∇s)

    Tx = SVector{d,Float64}
    TJ = typeof(zero(Tx)*zero(Tx)')

    # Fill coo values
    i_coo = 0
    for cell in 1:ncells
        rid = cell_to_rid[cell]
        nodes = cell_to_nodes[cell]
        Ae = Aes[rid]
        ue = ues[rid]
        fe = fes[rid]
        nl = length(nodes)
        ∇se = ∇s[rid]
        ∇xe = ∇x[rid]
        se = s[rid]
        we = w[rid]
        nq = length(we)
        fill!(Ae,zero(eltype(Ae)))
        fill!(fe,zero(eltype(Ae)))
        fill!(ue,zero(eltype(Ae)))
        for k in 1:nl
            nk = nodes[k]
            gk = node_to_free_node[nk]
            if  gk <= nfree
                continue
            end
            ue[k] = u_dirichlet[gk-nfree]
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
            for j in 1:nl
                for i in 1:nl
                    Ae[i,j] += ∇xe[i]⋅∇xe[j]*dV
                end
            end
            fx = f(xint)
            for k in 1:nl
                fe[k] +=  fx*se[iq,k]*dV
            end
        end
        for i in 1:nl
            for j in 1:nl
                fe[i] -= Ae[i,j]*ue[j]
            end
        end
        for i in 1:nl
            ni = nodes[i]
            gi = node_to_free_node[ni]
            if gi > nfree
                continue
            end
            b[gi] += fe[i]
        end
        for j in 1:nl
            nj = nodes[j]
            gj = node_to_free_node[nj]
            if  gj > nfree
                continue
            end
            for i in 1:nl
                ni = nodes[i]
                gi = node_to_free_node[ni]
                if gi > nfree
                    continue
                end
                i_coo += 1
                I_coo[i_coo] = gi
                J_coo[i_coo] = gj
                V_coo[i_coo] = Ae[i,j]
            end
        end
    end

    A = sparse(I_coo,J_coo,V_coo,nfree,nfree)

    A,b
end

function add_neumann_bcs!(b,state)

    neum_face_to_face = state.neumann_bcs.neum_face_to_face
    if length(neum_face_to_face) == 0
        return nothing
    end
    face_to_rid = state.face_isomap.face_to_rid
    face_to_nodes = state.face_isomap.face_to_nodes
    ∇s_f = state.face_isomap.rid_to_shape_grads
    s_f = state.face_isomap.rid_to_shape_vals
    node_to_x = state.face_isomap.node_to_coords
    free_and_dirichlet_nodes = state.dirichlet_bcs.free_and_dirichlet_nodes
    g = state.user_funs.g
    w_f = state.face_integration.rid_to_weights
    d = state.cell_integration.d

    fes = map(i->zeros(size(i,2)),s_f)
    Tx = SVector{d,Float64}
    Ts = typeof(zero(SVector{d-1,Float64})*zero(Tx)')
    node_to_free_node = glk.permutation(free_and_dirichlet_nodes)
    nfree = length(first(free_and_dirichlet_nodes))

    for face in neum_face_to_face
        rid = face_to_rid[face]
        nodes = face_to_nodes[face]
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
            ni = nodes[i]
            gi = node_to_free_node[ni]
            if gi > nfree
                continue
            end
            b[gi] += fe[i]
        end
    end

    nothing
end

function solve_system(A,b,state)
    solve = state.solve
    free_and_dirichlet_nodes = state.dirichlet_bcs.free_and_dirichlet_nodes
    u_dirichlet = state.dirichlet_bcs.u_dirichlet
    node_to_x = state.cell_isomap.node_to_coords

    nnodes = length(node_to_x)
    u_free = solve(A,b)
    uh = zeros(nnodes)
    uh[first(free_and_dirichlet_nodes)] = u_free
    uh[last(free_and_dirichlet_nodes)] = u_dirichlet
    uh
end

function add_basic_info!(results,state)
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

function integrate_error_norms!(results,uh,state)

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
    Tx = SVector{d,Float64}
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

    eh1 = sqrt(eh1)
    el2 = sqrt(el2)
    results[:eh1] = eh1
    results[:el2] = el2
    results
end

function export_results(uh,params,state)
    mesh = params[:mesh]
    example_path = params[:example_path]
    node_to_tag = state.dirichlet_bcs.node_to_tag
    vtk_grid(example_path,glk.vtk_args(mesh)...) do vtk
        glk.vtk_physical_groups!(vtk,mesh)
        vtk["tag"] = node_to_tag
        vtk["uh"] = uh
        vtk["dim"] = glk.face_dim(mesh)
    end
end



function diag!(diagA,A::SparseMatrixCSC)
    # TODO improve with binary search
    colptr = A.colptr
    rowval = A.rowval
    nzval = A.nzval
    ncols = length(colptr)-1
    for j in 1:ncols
        pini = colptr[j]
        pend = colptr[j+1]-1
        for p in pini:pend
            i = rowval[p]
            if j != i
                continue
            end
            diagA[j] = nzval[p]
        end
    end
    diagA
end

function spmv!(b,A,x)
    T = eltype(b)
    alpha = zero(T)
    beta = one(T)
    spmv!(b,A,x,alpha,beta)
    b
end

function spmv!(b,A::SparseMatrixCSC,x,alpha,beta)
    colptr = A.colptr
    rowval = A.rowval
    nzval = A.nzval
    nrows,ncols = size(A)
    if alpha == zero(eltype(b))
        b .= alpha
    else
        b .= alpha .* b
    end
    for j in 1:ncols
        xj = x[j]
        pini = colptr[j]
        pend = colptr[j+1]-1
        for p in pini:pend
            i = rowval[p]
            aij = nzval[p]
            b[i] += beta*aij*xj
        end
    end
    b
end

function spmm(A::SparseMatrixCSC,B::SparseMatrixCSC)
    # TODO implementation
    A*B
end

const SparseMatrixCSCT = Transpose{R,SparseMatrixCSC{R,S}} where {R,S}

function sprap(R::SparseMatrixCSCT,A::SparseMatrixCSC,P::SparseMatrixCSC)
    # TODO implementation
    R*A*P
end

function jacobi!(x,A,b;kwargs...)
    setup = jacobi_setup(x,A,b;kwargs...)
    jacobi!(x,A,b,setup)
end

function jacobi_setup(x,A,b;maxiters,omega=1)
    options = (;maxiters,omega)
    jacobi_setup(x,A,b,options)
end

function jacobi_setup(x,A::SparseMatrixCSC,b,options)
  xnew = similar(x)
  diagA = similar(x)
  diag!(diagA,A)
  (;xnew,diagA,options)
end

function jacobi_setup!(setup,A::SparseMatrixCSC)
    diag!(setup.diagA,A)
    setup
end

function jacobi!(xin,A::SparseMatrixCSC,b,setup)
    maxiters = setup.options.maxiters
    omega = setup.options.omega
    x = xin
    xnew = setup.xnew
    diagA = setup.diagA
    colptr = A.colptr
    rowval = A.rowval
    nzval = A.nzval
    nrows,ncols = size(A)
    for iter in 1:maxiters
        xnew .= b .+ diagA .* x ./ omega
        for j in 1:ncols
            xj = x[j]
            pini = colptr[j]
            pend = colptr[j+1]-1
            for p in pini:pend
                i = rowval[p]
                aij = nzval[p]
                xnew[i] -= aij*xj
            end
        end
        xnew .= omega .* xnew ./ diagA
        x,xnew = xnew,x
    end
    xin .= x
    xin
end

function aggregate(A;epsilon)
    options = (;epsilon)
    aggregate(A,options)
end

function aggregate(A::SparseMatrixCSC,options)

    epsi = options.epsilon
    typeof_aggregate = Int32
    typeof_strength = eltype(A.nzval)

    nnodes = size(A,1)
    pending = typeof_aggregate(0)
    isolated = typeof_aggregate(-1)
    
    node_to_aggregate = fill(pending,nnodes)
    node_to_old_aggregate = similar(node_to_aggregate)

    diagA = zeros(eltype(A),nnodes)
    diag!(diagA,A)

    node_to_neigs = jagged_array(A.rowval,A.colptr)
    node_to_vals = jagged_array(A.nzval,A.colptr)
    strongly_connected = (node,ineig) -> begin
        neig = node_to_neigs[node][ineig]
        aii = diagA[node]
        ajj = diagA[neig]
        aij = node_to_vals[node][ineig]
        abs(aij) > epsi*sqrt(aii*ajj)
    end
    coupling_strength = (node,ineig) -> begin
        abs(node_to_vals[node][ineig])
    end

    # Initialization
    for node in 1:nnodes
        neigs = node_to_neigs[node]
        isolated_node = count(i->i!=node,neigs) == 0
        if isolated_node
            node_to_aggregate[node] = isolated
        end
    end

    # Step 1
    aggregate = typeof_aggregate(0)
    for node in 1:nnodes
        if node_to_aggregate[node] != pending
            continue
        end
        neigs = node_to_neigs[node]
        nneigs = length(neigs)
        all_pending = true
        for ineig in 1:nneigs
            neig = neigs[ineig]
            if neig == node || !strongly_connected(node,ineig)
                continue
            end
            all_pending &= (node_to_aggregate[neig] == pending)
        end
        if !all_pending
            continue
        end
        aggregate += typeof_aggregate(1)
        node_to_aggregate[node] = aggregate
        for ineig in 1:nneigs
            neig = neigs[ineig]
            if neig == node || !strongly_connected(node,ineig)
                continue
            end
            node_to_aggregate[neig] = aggregate
        end
    end

    # Step 2
    copy!(node_to_old_aggregate,node_to_aggregate)
    for node in 1:nnodes
        if node_to_aggregate[node] != pending
            continue
        end
        strength = zero(typeof_strength)
        neigs = node_to_neigs[node]
        nneigs = length(neigs)
        for ineig in 1:nneigs
            neig = neigs[ineig]
            if neig == node || !strongly_connected(node,ineig)
                continue
            end
            neig_aggregate = node_to_old_aggregate[neig]
            if neig_aggregate != pending && neig_aggregate != isolated
                neig_strength = coupling_strength(node,ineig)
                if neig_strength > strength
                    strength = neig_strength
                    node_to_aggregate[node] = neig_aggregate
                end
            end
        end
    end

    # Step 3
    for node in 1:nnodes
        if node_to_aggregate[node] != pending
            continue
        end
        aggregate += typeof_aggregate(1)
        node_to_aggregate[node] = aggregate
        neigs = node_to_neigs[node]
        nneigs = length(neigs)
        for ineig in 1:nneigs
            neig = neigs[ineig]
            if neig == node || !strongly_connected(node,ineig)
                continue
            end
            neig_aggregate = node_to_old_aggregate[neig]
            if neig_aggregate == pending || neig_aggregate == isolated
                node_to_aggregate[neig] = aggregate
            end
        end
    end
    naggregates = aggregate

    ## Compression
    aggregate_to_nodes_ptrs = zeros(Int,naggregates+1)
    for node in 1:nnodes
        agg = node_to_aggregate[node]
        if agg == pending
            continue
        end
        aggregate_to_nodes_ptrs[agg+1] += 1
    end
    length_to_ptrs!(aggregate_to_nodes_ptrs)
    ndata = aggregate_to_nodes_ptrs[end]-1
    aggregate_to_nodes_data = zeros(Int,ndata)
    for node in 1:nnodes
        agg = node_to_aggregate[node]
        if agg == pending
            continue
        end
        p = aggregate_to_nodes_ptrs[agg]
        aggregate_to_nodes_data[p] = node
        aggregate_to_nodes_ptrs[agg] += 1
    end
    rewind_ptrs!(aggregate_to_nodes_ptrs)
    aggregate_to_nodes = JaggedArray(aggregate_to_nodes_data,aggregate_to_nodes_ptrs)

    aggregate_to_nodes
end

function prolongator(A,agg_to_nodes;omega=1)
    options = (;omega)
    prolongator(A,agg_to_nodes,options)
end

function prolongator(A::SparseMatrixCSC,agg_to_nodes,options)
    # TODO Optimize
    omega = options.omega
    nrows,ncols = size(A)
    diagA = zeros(eltype(A.nzval),nrows)
    diag!(diagA,A) # TODO We are taking the diagonal to many times
    nagg = length(agg_to_nodes)
    P0 = SparseMatrixCSC(ncols,nagg,agg_to_nodes.ptrs,agg_to_nodes.data,ones(length(agg_to_nodes.data)))
    Dinv = sparse(1:nrows,1:nrows,1 ./ diagA,nrows,nrows)
    Id = sparse(1:nrows,1:nrows,ones(nrows),nrows,nrows)
    P = (I-omega*Dinv*A)*P0
end

function smooth_aggregation(A;epsilon,omega)
    aggrs = aggregate(A;epsilon)
    P = prolongator(A,aggrs;omega)
    R = transpose(P)
    Ac = sprap(R,A,P)
    Ac,R,P
end

function amg!(x,A,b;kwargs...)
    setup = amg_setup(x,A,b;kwargs...)
    amg!(x,A,b,setup)
end

function amg_setup(x,A,b;fine_params,coarse_solver)
    nlevels = length(fine_params)
    fine_levels_setup =  map(fine_params) do fine_level
        (;pre_smoother,pos_smoother,coarsen,cycle) = fine_level
        pre_setup = pre_smoother.setup(x,A,b;pre_smoother.params...)
        pos_setup = pos_smoother.setup(x,A,b;pos_smoother.params...)
        pre_smoother! = (x,b) -> pre_smoother.solve(x,A,b,pre_setup)
        pos_smoother! = (x,b) -> pos_smoother.solve(x,A,b,pos_setup)
        Ac,R,P = coarsen.setup(A;coarsen.params...)
        r = similar(b)
        e = similar(x)
        rc = similar(r,axes(Ac,1))
        ec = similar(e,axes(Ac,2))
        level_setup = (;R,A,Ac,P,r,e,rc,ec,pre_smoother!,pos_smoother!,cycle)
        A = Ac
        b = rc
        x = ec
        level_setup
    end
    coarse_setup = coarse_solver.setup(x,A,b,coarse_solver.params...)
    coarse_solver! = (x,b)->coarse_solver.solve(x,A,b,coarse_setup)
    (;nlevels,fine_levels_setup,coarse_solver!)
end

function amg!(x,A,b,setup)
    level = 1
    amg_cycle!(x,b,setup,level)
    x
end

function amg_cycle!(x,b,setup,level)
    if level == setup.nlevels+1
        return setup.coarse_solver!(x,b)
    end
    level_setup = setup.fine_levels_setup[level]
    (;R,A,P,r,e,rc,ec,pre_smoother!,pos_smoother!,cycle) = level_setup
    pre_smoother!(x,b)
    copy!(r,b)
    spmv!(r,A,x,-1,1) # TODO mul! better? Provably, yes
    mul!(rc,R,r) # TODO spmv!
    fill!(ec,zero(eltype(ec)))
    cycle(ec,rc,setup,level+1)
    spmv!(e,P,ec)
    x .-= e
    pos_smoother!(x,b)
    x
end

function v_cycle!(args...)
    amg_cycle!(args...)
end

function w_cycle!(args...)
    amg_cycle!(args...)
    amg_cycle!(args...)
end


end # module



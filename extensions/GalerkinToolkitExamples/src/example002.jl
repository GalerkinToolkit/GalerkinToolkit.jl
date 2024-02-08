module Example002

import GalerkinToolkit as gk
import ForwardDiff
using StaticArrays
using LinearAlgebra
using SparseArrays
using WriteVTK
using PartitionedArrays
using TimerOutputs
using PetscCall

function main(params_in)

    # Dict to collect results
    results = Dict{Symbol,Any}()

    # Process params
    params_default = default_params()
    params = add_default_params(params_in,params_default)
    parts = linear_indices(partition(params[:mesh]))

    timer = TimerOutput()

    # Setup main data structures
    @timeit timer "setup"  state = setup(params,timer)
    add_basic_info(results,params,state)

    # Assemble system and solve it
    @timeit timer "assemble_system" A,b = assemble_system(state)
    @timeit timer "solve_system" x = solve_system(A,b,params,state)

    # Post process
    @timeit timer "setup_uh"  uh = setup_uh(x,state)
    @timeit timer "integrate_error_norms" integrate_error_norms(results,uh,state)
    @timeit timer "export_results"  export_results(uh,params,state)

    map_main(state.local_states) do local_state
        display(local_state.timer)
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
    np = 4
    parts = DebugArray(LinearIndices((np,)))
    domain = (0,1,0,1)
    cells_per_dir = (10,10)
    parts_per_dir = (2,2)
    ghost_layers = 0
    mesh = gk.cartesian_mesh(domain,cells_per_dir,parts_per_dir;parts,ghost_layers)
    outdir = mkpath(joinpath(@__DIR__,"..","output"))
    params = Dict{Symbol,Any}()
    params[:u] = (x) -> sum(x)
    params[:f] = (x) -> 0.0
    params[:g] = (x) -> 0.0
    params[:dirichlet_tags] = ["boundary"]
    params[:neumann_tags] = String[]
    params[:integration_degree] = 2
    params[:solver] = lu_solver()
    params[:export_vtu] = true
    params[:mesh] = mesh
    params[:example_path] = joinpath(outdir,"example002")
    params
end

# Not a parallel solver
# just for debugging purposes
function lu_solver()
    setup(x,A,b) = lu(A)
    setup! = lu!
    solve! = ldiv!
    finalize!(S) = nothing
    (;setup,setup!,solve!,finalize!)
end

function ksp_solver()
    setup = PetscCall.ksp_setup
    setup! = PetscCall.ksp_setup!
    solve! = PetscCall.ksp_solve!
    finalize! = PetscCall.ksp_finalize!
    (;setup,setup!,solve!,finalize!)
end

function setup(params,timer)
    mesh = params[:mesh]
    local_states_0 = map(partition(mesh)) do mesh
        local_params = copy(params)
        local_params[:mesh] = mesh
        local_setup(local_params,timer)
    end
    @timeit timer "setup_dofs" dofs,local_dofs = setup_dofs(params,local_states_0)
    local_states = map((a,dofs)->(;dofs,a...),local_states_0,local_dofs)
    state = (;local_states,dofs,timer)
    state
end

function local_setup(params,timer)
    @timeit timer "dirichlet_bcs" dirichlet_bcs = setup_dirichlet_bcs(params)
    @timeit timer "neumann_bcs" neumann_bcs = setup_neumann_bcs(params)
    @timeit timer "cell_integration" cell_integration = setup_integration(params,:cells)
    @timeit timer "face_integration" face_integration = setup_integration(params,:faces)
    @timeit timer "cell_isomap" cell_isomap = setup_isomap(params,cell_integration)
    @timeit timer "face_isomap" face_isomap = setup_isomap(params,face_integration)
    @timeit timer "user_funs" user_funs = setup_user_funs(params)
    solver = params[:solver]
    state = (;timer,solver,dirichlet_bcs,neumann_bcs,cell_integration,cell_isomap,face_integration,face_isomap,user_funs)
end

function setup_dofs(params,local_states)
    mesh = params[:mesh]
    node_partition = gk.node_partition(mesh)
    global_node_to_mask = pfill(false,node_partition)
    function fillmask!(node_to_mask,state)
        free_nodes = first(state.dirichlet_bcs.free_and_dirichlet_nodes)
        node_to_mask[free_nodes] .= true
    end
    map(fillmask!,partition(global_node_to_mask),local_states)
    global_dof_to_node, global_node_to_dof = find_local_indices(global_node_to_mask)
    dof_partition = partition(axes(global_dof_to_node,1))
    function setup_local_dofs(state,node_to_dof)
        free_and_dirichlet_nodes = state.dirichlet_bcs.free_and_dirichlet_nodes
        node_to_free_node = gk.permutation(free_and_dirichlet_nodes)
        n_free = length(first(free_and_dirichlet_nodes))
        n_dofs = n_free
        cell_to_nodes = state.cell_isomap.face_to_nodes
        face_to_nodes = state.face_isomap.face_to_nodes
        function map_node_to_dof(node)
            free_node = node_to_free_node[node]
            if free_node > n_free
                return Int(-(free_node-n_free))
            end
            dof = node_to_dof[node]
            Int(dof)
        end
        cell_to_dofs_data = map_node_to_dof.(cell_to_nodes.data)
        face_to_dofs_data = map_node_to_dof.(face_to_nodes.data)
        cell_to_dofs = JaggedArray(cell_to_dofs_data,cell_to_nodes.ptrs)
        face_to_dofs = JaggedArray(face_to_dofs_data,face_to_nodes.ptrs)
        (;cell_to_dofs,face_to_dofs,n_dofs)
    end
    local_dofs = map(setup_local_dofs,local_states,partition(global_node_to_dof))
    (;dof_partition,global_node_to_dof,global_dof_to_node), local_dofs
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

function add_basic_info(results,params,state)
    mesh = params[:mesh]
    nfree = length(state.dofs.global_dof_to_node)
    nnodes = gk.num_nodes(mesh)
    ncells = gk.num_faces(mesh,gk.num_dims(mesh))
    results[:nfree] = nfree
    results[:nnodes] = nnodes
    results[:ncells] = ncells
    results
end

function assemble_system(state)
    timer = state.timer
    dof_partition = partition(axes(state.dofs.global_dof_to_node,1))
    @timeit timer "assemble_system_coo" I,J,V,b_partition = map(assemble_system_coo,state.local_states) |> tuple_of_arrays
    @timeit timer "psparse" A = psparse(I,J,V,dof_partition,dof_partition;subassembled=true,indices=:local) |> fetch
    row_partition = partition(axes(A,1))
    b = PVector(b_partition,dof_partition)
    @timeit timer "assemble_b" b = assemble(b,row_partition) |> fetch
    A,b
end

function assemble_system_coo(state)
    timer = state.timer
    @timeit timer "assemble_sybmolic" args = assemble_sybmolic(state)
    @timeit timer "assemble_numeric_cells!" assemble_numeric_cells!(args...,state)
    @timeit timer "assemble_numeric_faces!" assemble_numeric_faces!(args...,state)
    args
end

function assemble_sybmolic(state)
    cell_to_dofs = state.dofs.cell_to_dofs

    n_coo = 0
    ncells = length(cell_to_dofs)
    for cell in 1:ncells
        dofs = cell_to_dofs[cell]
        ndofs = length(dofs)
        for i in 1:ndofs
            if !(dofs[i]>0)
                continue
            end
            for j in 1:ndofs
                if !(dofs[j]>0)
                    continue
                end
                n_coo += 1
            end
        end
    end

    I_coo = Vector{Int32}(undef,n_coo)
    J_coo = Vector{Int32}(undef,n_coo)
    V_coo = Vector{Float64}(undef,n_coo)
    b = zeros(Float64,state.dofs.n_dofs)

    n_coo = 0
    for cell in 1:ncells
        dofs = cell_to_dofs[cell]
        ndofs = length(dofs)
        for i in 1:ndofs
            if !(dofs[i]>0)
                continue
            end
            dofs_i = dofs[i]
            for j in 1:ndofs
                if !(dofs[j]>0)
                    continue
                end
                n_coo += 1
                I_coo[n_coo] = dofs_i
                J_coo[n_coo] = dofs[j]
            end
        end
    end

    I_coo,J_coo,V_coo,b
end

function assemble_numeric_cells!(I_coo,J_coo,V_coo,b,state)

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

    # Allocate auxiliary buffers
    ∇x = map(i->similar(i,size(i,2)),∇s)
    Aes = map(i->zeros(size(i,2),size(i,2)),∇s)
    fes = map(i->zeros(size(i,2)),∇s)
    ∇st = map(m->collect(permutedims(m)),∇s)
    st = map(m->collect(permutedims(m)),s)

    Tx = eltype(node_to_x)
    TJ = typeof(zero(Tx)*zero(Tx)')

    i_coo = 0
    for cell in 1:ncells
        rid = cell_to_rid[cell]
        nodes = cell_to_nodes[cell]
        dofs = cell_to_dofs[cell]
        Ae = Aes[rid]
        fe = fes[rid]
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
                ∇xe[k] = invJt*∇ste[k,iq]
            end
            for j in 1:nl
                ∇xej = ∇xe[j]
                for i in 1:nl
                    Ae[i,j] += ∇xe[i]⋅∇xej*dV
                end
            end
            fx = f(xint)
            for k in 1:nl
                fe[k] +=  fx*ste[k,iq]*dV
            end
        end
        # Set the result in the output array
        for i in 1:nl
            if !(dofs[i]>0)
                continue
            end
            b[dofs[i]] += fe[i]
            for j in 1:nl
                if !(dofs[j]>0)
                    uj = u_dirichlet[-dofs[j]]
                    b[dofs[i]] -= Ae[i,j]*uj
                    continue
                end
                i_coo += 1
                V_coo[i_coo] = Ae[i,j]
            end
        end
    end

end

function assemble_numeric_faces!(I_coo,J_coo,V_coo,b,state)

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
            b[dofs[i]] += fe[i]
        end
    end
end

function solve_system(A,b,params,state)
    timer = state.timer
    solver = params[:solver]
    x = similar(b,axes(A,2))
    @timeit timer "solver.setup" setup = solver.setup(x,A,b)
    @timeit timer "solver.solve!" solver.solve!(x,setup,b)
    @timeit timer "solver.finalize!" solver.finalize!(setup)
    x
end

function setup_uh(x,state)
    dofs = axes(state.dofs.global_dof_to_node,1)
    local_states = state.local_states
    global_dof_to_x = similar(x,dofs)
    global_dof_to_x .= x
    consistent!(global_dof_to_x) |> wait
    function setup_local_uh(dof_to_x,state,node_to_dof)
        free_and_dirichlet_nodes = state.dirichlet_bcs.free_and_dirichlet_nodes
        node_to_x = state.cell_isomap.node_to_coords
        free_node_to_node = first(free_and_dirichlet_nodes)
        @views u_free = dof_to_x[node_to_dof[free_node_to_node]]
        u_dirichlet = state.dirichlet_bcs.u_dirichlet
        nnodes = length(node_to_x)
        local_uh = zeros(nnodes)
        local_uh[first(free_and_dirichlet_nodes)] = u_free
        local_uh[last(free_and_dirichlet_nodes)] = u_dirichlet
        local_uh
    end
    uh = map(setup_local_uh,partition(global_dof_to_x),local_states,partition(state.dofs.global_node_to_dof))
    uh
end

function integrate_error_norms(results,uh,state)
    local_states = state.local_states
    eh1², el2² = map(integrate_error_norms_loop,uh,local_states) |> tuple_of_arrays
    eh1 = sqrt(sum(eh1²))
    el2 = sqrt(sum(el2²))
    results[:eh1] = eh1
    results[:el2] = el2
    results
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

function export_results(uh,params,state)
    if ! params[:export_vtu]
        return nothing
    end
    example_path = params[:example_path]
    mesh = params[:mesh]
    ranks = linear_indices(gk.node_partition(mesh))
    np = length(ranks)
    map(partition(mesh),ranks,uh,gk.index_partition(mesh)) do mesh,rank,uh,ids
        pvtk_grid(example_path,gk.vtk_args(mesh)...;part=rank,nparts=np) do vtk
            gk.vtk_physical_faces!(vtk,mesh)
            vtk["piece",VTKCellData()] = fill(rank,sum(gk.num_faces(mesh)))
            vtk["owner",VTKPointData()] = local_to_owner(gk.node_indices(ids))
            vtk["uh",VTKPointData()] = uh
        end
    end
end

end # module

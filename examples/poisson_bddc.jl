module PoissonBDDC

import GalerkinToolkit as glk
import ForwardDiff
using StaticArrays
using LinearAlgebra
using SparseArrays
using WriteVTK
using PartitionedArrays
using Metis
using IterativeSolvers

const CORNER = 0
const EDGE = 1
const FACE = 2

struct DDOperator{A,B<:PRange}
    matrices::A
    ids::B
end

Base.eltype(A::DDOperator) = eltype(eltype(A.matrices))
Base.eltype(::Type{<:DDOperator{A}}) where A = eltype(eltype(A))
Base.size(A::DDOperator) = (length(A.ids),length(A.ids))
Base.axes(A::DDOperator) = (A.ids,A.ids)
Base.size(A::DDOperator,i::Integer) = length(A.ids)
Base.axes(A::DDOperator,i::Integer) = A.ids

function Base.:*(A::DDOperator,x::PVector)
    b = similar(x)
    mul!(b,A,x)
    b
end

function LinearAlgebra.mul!(b::PVector,A::DDOperator,x::PVector)
    fill!(b,zero(eltype(b)))
    consistent!(x) |> wait # TODO we can add latency hiding here
    map(mul!,local_values(b),A.matrices,local_values(x))
    assemble!(b) |> wait
    b
end

struct BDDC{T<:DDOperator,S}
    A::T
    workspace::S
end

function LinearAlgebra.ldiv!(x::PVector,M::BDDC,r::PVector)
    consistent!(r) |> wait
    workspace = M.workspace
    # Coarse correction
    rci = map(workspace,local_values(r)) do workspace,ri
        Wi = workspace.Wi
        Phii = workspace.Phii
        PhiiT = transpose(Phii)
        PhiiT*(Wi .* ri)# TODO allocations here
    end
    rank_to_rci = gather(rci,destination=MAIN) # TODO allocations here
    rank_to_v1ci = map(rank_to_rci,workspace) do rank_to_rci,workspace
        rank = workspace.rank
        if rank == MAIN
            ncdofs = workspace.ncdofs
            Pc = workspace.Pc
            rank_to_cdofs = workspace.rank_to_cdofs
            rc = zeros(ncdofs)#TODO allocation here
            nparts = length(rank_to_cdofs)
            for part in 1:nparts
                cdofs = rank_to_cdofs[part]
                rci = rank_to_rci[part]
                rc[cdofs] += rci
            end
            v1ic = Pc\rc # TODO allocation here
            JaggedArray([ v1ic[cdofs] for cdofs in rank_to_cdofs])# TODO allocation here
        else
            JaggedArray([Float64[]])# TODO allocation here
        end
    end
    v1ci = scatter(rank_to_v1ci) # TODO allocation here
    v1 = similar(x) # TODO allocation here
    map(workspace,local_values(v1),v1ci) do workspace,v1i,v1ci
        Wi = workspace.Wi
        Phii = workspace.Phii
        v1i .= Wi.*(Phii*v1ci) # TODO allocations here
    end
    assemble!(v1) |> wait # TODO group communications

    # Neumann correction
    v2 = similar(x) # TODO allocation here
    map(workspace,local_values(v2),local_values(r)) do workspace,v2i,r
        Wi = workspace.Wi
        nldofs = length(Wi)
        Pn = workspace.Pn
        ncldofs = workspace.ncldofs
        rhs = vcat((Wi.*r),zeros(ncldofs)) # TODO allocations here
        sol = Pn\rhs # TODO allocation here
        zi = view(sol,1:nldofs)
        v2i .= Wi .* zi
    end
    assemble!(v2) |> wait # TODO group communications

    # Dirichlet correction
    A = M.A
    r1 = r - A*(v1+v2) # TODO allocations here
    consistent!(r1) |> wait
    v3 = similar(x) # TODO allocation here
    map(workspace,local_values(v3),local_values(r1)) do workspace, v3i, r1i
        RIi = workspace.RIi
        Pd = workspace.Pd
        fill!(v3i,zero(eltype(v3i)))
        v3i[RIi] .= Pd\r1i[RIi] # TODO allocations here
    end
    #assemble!(v3) |> wait # TODO really needed, boundary is 0? TODO group communications

    # Final correction
    x .= v1 .+ v2 .+ v3
    x
end

function main(params_in)

    # Dict to collect results
    results = Dict{Symbol,Any}()

    # Process params
    params_default = default_params()
    params = add_default_params(params_in,params_default)

    # Setup main data structures
    state = setup(params)

    # Assemble system and solve it
    A,b = assemble_system(state)
    uh = solve_system(results,A,b,state)

    # Post process
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

function default_params()
    mesh = glk.cartesian_mesh((0,1,0,1),(10,10),complexify=false)
    np = 4
    ranks = DebugArray(LinearIndices((np,)))
    pmesh = glk.partition_mesh(Metis.partition,ranks,mesh,via=:cells)
    outdir = mkpath(joinpath(@__DIR__,"..","output"))
    params = Dict{Symbol,Any}()
    params[:pmesh] = glk.cartesian_mesh((0,1,0,1),(10,10))
    params[:u] = (x) -> sum(x)
    params[:f] = (x) -> 0.0
    params[:g] = (x) -> 0.0
    params[:dirichlet_tags] = ["boundary"]
    params[:neumann_tags] = String[]
    params[:integration_degree] = 2
    params[:example_path] = joinpath(outdir,"poisson_bddc")
    params[:bddc_type] = (CORNER,EDGE,FACE)
    params
end

function setup_dirichlet_bcs(mesh,params)
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

function setup_integration(mesh,params,objects)
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

function setup_isomap(mesh,params,integration)
    rid_to_coords = integration.rid_to_coords
    face_to_rid = integration.face_to_rid
    d = integration.d
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

function setup_neumann_bcs(mesh,params)
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
    pmesh = params[:pmesh]
    bddc_type = params[:bddc_type]
    ranks = linear_indices(pmesh)
    state1 = map(ranks,pmesh) do rank,mesh
        dirichlet_bcs = setup_dirichlet_bcs(mesh,params)
        neumann_bcs = setup_neumann_bcs(mesh,params)
        cell_integration = setup_integration(mesh,params,:cells)
        face_integration = setup_integration(mesh,params,:faces)
        cell_isomap = setup_isomap(mesh,params,cell_integration)
        face_isomap = setup_isomap(mesh,params,face_integration)
        user_funs = setup_user_funs(params)
        lnode_to_colors = glk.local_node_colors(mesh)
        lnodes = glk.local_nodes(mesh)
        (;
         dirichlet_bcs,
         neumann_bcs,
         cell_integration,
         cell_isomap,
         face_integration,
         face_isomap,
         user_funs,
         rank,
         bddc_type,
         lnode_to_colors,
         lnodes)
    end
    state2 = setup_dof_partition(state1)
    state3 = setup_coarse_dofs(state2)
    state3
end

function setup_coarse_dofs(state1)
    state2 = map(state1) do state1
        rank = state1.rank
        lnode_to_colors = state1.lnode_to_colors
        dirichlet_bcs = state1.dirichlet_bcs
        bddc_type = state1.bddc_type
        free_and_dirichlet_lnodes = dirichlet_bcs.free_and_dirichlet_nodes
        nldofs = length(first(free_and_dirichlet_lnodes))
        lnode_to_ldof = glk.permutation(free_and_dirichlet_lnodes)
        nlnodes = length(lnode_to_colors)
        interface_node_to_lnode = findall(1:nlnodes) do lnode
            colors = lnode_to_colors[lnode]
            ldof = lnode_to_ldof[lnode]
            length(colors)!=1 && ldof <= nldofs
        end
        interface_node_to_colors = view(lnode_to_colors,interface_node_to_lnode)
        # assumes sorted colors
        clnode_to_colors = unique(interface_node_to_colors)
        clnode_to_owner = map(maximum,clnode_to_colors)
        nclnodes = length(clnode_to_colors)
        interface_node_to_clnode = indexin(interface_node_to_colors,clnode_to_colors)
        clnode_to_interface_nodes = glk.inverse_index_map(interface_node_to_clnode,nclnodes)
        f = interface_node->lnode_to_ldof[interface_node_to_lnode[interface_node]]
        clnode_to_interface_nodes.data .= f.(clnode_to_interface_nodes.data)
        clnode_to_ldofs = clnode_to_interface_nodes
        clnode_to_type = fill(EDGE,nclnodes)
        clnode_to_type[ map(i->length(i) == 2,clnode_to_colors) ] .= FACE
        clnode_to_type[ map(i->length(i) == 1,clnode_to_ldofs) ] .= CORNER
        mask = map(i->i in bddc_type,clnode_to_type)
        cldof_to_ldofs = JaggedArray(clnode_to_ldofs[mask])
        cldof_to_owner = clnode_to_owner[mask]
        nocdofs = count(owner->owner==rank,cldof_to_owner)
        (;cldof_to_ldofs,cldof_to_owner,nocdofs)
    end
    nown = map(i->i.nocdofs,state2)
    ncdofs = sum(nown)
    cdof_partition = variable_partition(nown,ncdofs)
    node_partition = map(i->i.lnodes,state1)
    v = PVector{Vector{Int}}(undef,node_partition)
    map(state1,state2,local_values(v),cdof_partition) do state1,state2,lnode_to_v,codofs
        rank = state1.rank
        cldof_to_owner = state2.cldof_to_owner
        cldof_to_ldofs = state2.cldof_to_ldofs
        dirichlet_bcs = state1.dirichlet_bcs
        ldof_to_lnode = first(dirichlet_bcs.free_and_dirichlet_nodes)
        codof_to_cdof = own_to_global(codofs)
        codof = 0
        for (cldof,owner) in enumerate(cldof_to_owner)
            if owner != rank
                continue
            end
            codof += 1
            cdof = codof_to_cdof[codof]
            ldofs = cldof_to_ldofs[cldof]
            lnodes = view(ldof_to_lnode,ldofs)
            lnode_to_v[lnodes] .= cdof
        end
    end
    consistent!(v) |> wait
    cldof_to_cdof = map(state1,state2,local_values(v)) do state1, state2, lnode_to_v
        rank = state1.rank
        cldof_to_ldofs = state2.cldof_to_ldofs
        dirichlet_bcs = state1.dirichlet_bcs
        ldof_to_lnode = first(dirichlet_bcs.free_and_dirichlet_nodes)
        ncldofs = length(cldof_to_ldofs)
        my_cldof_to_cdof = zeros(Int,ncldofs)
        for cldof in 1:ncldofs
            ldofs = cldof_to_ldofs[cldof]
            cdof = lnode_to_v[ldof_to_lnode[first(ldofs)]]
            my_cldof_to_cdof[cldof] = cdof
        end
        my_cldof_to_cdof
    end
    rank_to_cdofs = gather(cldof_to_cdof,destination=MAIN)
    state3 = map(state1,state2,rank_to_cdofs) do state1, state2, rank_to_cdofs
        cldof_to_ldofs = state2.cldof_to_ldofs
        coarse_dofs = (;cldof_to_ldofs, rank_to_cdofs,ncdofs)
        (;coarse_dofs,state1...)
    end
    state3
end

function setup_dof_partition(state1)
    nodofs = map(state1) do state1
        rank = state1.rank
        dirichlet_bcs = state1.dirichlet_bcs
        free_and_dirichlet_lnodes = dirichlet_bcs.free_and_dirichlet_nodes
        lnodes = state1.lnodes
        lnode_to_owner = local_to_owner(lnodes)
        ldof_to_lnode = first(free_and_dirichlet_lnodes)
        count(lnode->lnode_to_owner[lnode]==rank,ldof_to_lnode)
    end
    ndofs = sum(nodofs)
    dof_partition = variable_partition(nodofs,ndofs)
    node_partition = map(i->i.lnodes,state1)
    v = PVector{Vector{Int}}(undef,node_partition)
    fill!(v,0)
    map(local_values(v),dof_partition,state1) do lnode_to_dof, odofs, state1
        rank = state1.rank
        lnodes = state1.lnodes
        dirichlet_bcs = state1.dirichlet_bcs
        ldof_to_lnode = first(dirichlet_bcs.free_and_dirichlet_nodes)
        odof_to_dof = own_to_global(odofs)
        lnode_to_owner = local_to_owner(lnodes)
        nldofs = length(ldof_to_lnode)
        odof = 0
        for ldof in 1:nldofs
            lnode = ldof_to_lnode[ldof]
            if lnode_to_owner[lnode] == rank
                odof += 1
                dof = odof_to_dof[odof]
                lnode_to_dof[lnode] = dof
            end
        end
    end
    consistent!(v) |> wait
    state2 = map(local_values(v),state1) do lnode_to_dof,state1
        rank = state1.rank
        lnodes = state1.lnodes
        dirichlet_bcs = state1.dirichlet_bcs
        lnode_to_owner = local_to_owner(lnodes)
        ldof_to_lnode = first(dirichlet_bcs.free_and_dirichlet_nodes)
        ldof_to_dof = lnode_to_dof[ldof_to_lnode]
        ldof_to_owner = lnode_to_owner[ldof_to_lnode]
        ldofs = LocalIndices(ndofs,rank,ldof_to_dof,ldof_to_owner)
        (;ldofs,state1...)
    end
    state2
end

function assemble_system(state)
    Alocal,blocal = map(assemble_system_locally,state) |> tuple_of_arrays
    map(add_neumann_bcs_locally!,blocal,state)
    dof_partition = map(i->i.ldofs,state)
    b = PVector(blocal,dof_partition)
    assemble!(b) |> wait # do not forget this
    A = DDOperator(Alocal,axes(b,1))
    A,b
end

function assemble_system_locally(state)

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
    I_coo = zeros(Int,n_coo)
    J_coo = zeros(Int,n_coo)
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

function add_neumann_bcs_locally!(b,state)

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

function solve_system(results,A,b,state)
    M = setup_BDDC(A,state)
    x = similar(b)
    fill!(x,zero(eltype(b)))
    IterativeSolvers.cg!(copy(x),A,b,Pl=M,verbose=i_am_main(state))
    ranks = linear_indices(state)
    t = PTimer(ranks)
    tic!(t,barrier=true)
    IterativeSolvers.cg!(x,A,b,Pl=M)
    toc!(t,"cg!")
    display(t)
    consistent!(x) |> wait
    node_partition = map(i->i.lnodes,state)
    uh = pzeros(node_partition)
    map(state,local_values(uh),local_values(x)) do state, uh, u_free
        free_and_dirichlet_nodes = state.dirichlet_bcs.free_and_dirichlet_nodes
        u_dirichlet = state.dirichlet_bcs.u_dirichlet
        uh[first(free_and_dirichlet_nodes)] = u_free
        uh[last(free_and_dirichlet_nodes)] = u_dirichlet
    end
    uh
end

function setup_BDDC(A,state)
    workspace1 = map(A.matrices,state) do Ai,state
        # Setup Neumann problems
        cldof_to_ldofs = state.coarse_dofs.cldof_to_ldofs
        nldofs = size(Ai,1)
        ncldofs = length(cldof_to_ldofs)
        cldof_to_coefs = map(cldof_to_ldofs) do ldofs
            fill(1/length(ldofs),length(ldofs))
        end
        mycolptr = cldof_to_ldofs.ptrs
        myrowval = cldof_to_ldofs.data
        mynzval = reduce(vcat,cldof_to_coefs)
        CiT = SparseMatrixCSC(nldofs,ncldofs,mycolptr,myrowval,mynzval)
        Ci = transpose(CiT)
        Zi = spzeros(ncldofs,ncldofs)
        Ti = [Ai CiT; Ci Zi]
        Pn = lu(Ti)# TODO with some work one could use Cholesky
        rhs = [zeros(nldofs,ncldofs);Matrix(I,ncldofs,ncldofs)]
        Phii = (Pn\rhs)[1:nldofs,:]

        # Setup Dirichlet problems
        lnode_to_colors = state.lnode_to_colors
        free_and_dirichlet_lnodes = state.dirichlet_bcs.free_and_dirichlet_nodes
        ldof_to_lnode = first(free_and_dirichlet_lnodes)
        nldofs = length(ldof_to_lnode)
        lnode_to_ldof = glk.permutation(free_and_dirichlet_lnodes)
        nlnodes = length(lnode_to_colors)
        idof_to_lnode = findall(1:nlnodes) do lnode
            colors = lnode_to_colors[lnode]
            ldof = lnode_to_ldof[lnode]
            length(colors)==1 && ldof <= nldofs
        end
        idof_to_ldof = lnode_to_ldof[idof_to_lnode]
        RIi = idof_to_ldof
        AIi = Ai[RIi,RIi]
        Pd = cholesky(AIi)

        # Setup weight
        ldof_to_colors = lnode_to_colors[ldof_to_lnode]
        ldof_to_w = map(colors->1/length(colors),ldof_to_colors)
        Wi = ldof_to_w

        rank_to_cdofs = state.coarse_dofs.rank_to_cdofs
        rank = state.rank
        ncdofs = state.coarse_dofs.ncdofs
        (;Pn,Phii,Pd,RIi,Wi,ncldofs,rank_to_cdofs,rank,ncdofs)
    end

    # Setup coarse solver
    Aci = map(workspace1,A.matrices) do workspace1,Ai
        Phii = workspace1.Phii
        PhiiT = transpose(Phii)
        PhiiT*Ai*Phii
    end
    rank_to_Aci = gather_matrix(Aci,destination=MAIN)
    workspace2 = map(rank_to_Aci,workspace1) do rank_to_Aci,workspace1
        rank_to_cdofs = workspace1.rank_to_cdofs
        rank = workspace1.rank
        ncdofs = workspace1.ncdofs
        if rank == MAIN
            nparts = length(rank_to_cdofs)
            ncoo = 0
            for part in 1:nparts
                cdofs = rank_to_cdofs[part]
                ncoo += length(cdofs)*length(cdofs)
            end
            Icoo = zeros(Int,ncoo)
            Jcoo = zeros(Int,ncoo)
            Vcoo = zeros(ncoo)
            ncoo = 0
            for part in 1:nparts
                cdofs = rank_to_cdofs[part]
                myAci = rank_to_Aci[part]
                for lj in 1:size(myAci,2)
                    for li in 1:size(myAci,1)
                        ncoo += 1
                        Icoo[ncoo] = cdofs[li]
                        Jcoo[ncoo] = cdofs[lj]
                        Vcoo[ncoo] = myAci[li,lj]
                    end
                end
            end
            Ac = sparse(Icoo,Jcoo,Vcoo,ncdofs,ncdofs)
            Ac = 0.5*(Ac+transpose(Ac)) # TODO cholesky requires exactly symmetric matrices (no rounding errors allowed)
            Pc = cholesky(Ac)
        else
            Pc = nothing
        end
        (;Pc,workspace1...)
    end
    BDDC(A,workspace2)
end

# TODO this should be handled by PartitionedArrays
function gather_matrix(Aci;destination=MAIN)
    # TODO, we can pack all the data in the same message
    s = gather(map(i->size(i),Aci);destination)
    data = gather(map(i->i[:],Aci);destination)
    ranks = linear_indices(Aci)
    map(ranks,data,s) do rank,data,s
        if rank == destination
            map(reshape,data,s)
        else
            nothing
        end
    end
end

function export_results(uh,params,state)
    example_path = params[:example_path]
    pmesh = params[:pmesh]
    ranks = linear_indices(pmesh)
    np = length(ranks)
    map(pmesh,ranks,local_values(uh)) do mesh,rank,uh
        pvtk_grid(example_path,glk.vtk_args(mesh)...;part=rank,nparts=np) do vtk
            glk.vtk_physical_groups!(vtk,mesh)
            vtk["piece",VTKCellData()] = fill(rank,sum(glk.num_faces(mesh)))
            vtk["owner",VTKPointData()] = local_to_owner(glk.local_nodes(mesh))
            vtk["uh",VTKPointData()] = uh
            vtk["interface",VTKPointData()] = map(colors->Int(length(colors)!=1),glk.local_node_colors(mesh))
        end
    end
end

function integrate_error_norms!(results,uh,state)
    errs = map(integrate_error_norms_local,local_values(uh),state)
    eh1 = sum(map(i->i.eh1,errs)) |> sqrt
    el2 = sum(map(i->i.el2,errs)) |> sqrt
    results[:eh1] = eh1
    results[:el2] = el2
    results
end

function integrate_error_norms_local(uh,state)

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

    (;eh1,el2)
end


end # module

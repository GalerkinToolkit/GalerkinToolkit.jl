
const DD_CORNER = 0
const DD_EDGE = 1
const DD_FACE = 2

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

bddc_cef() = (DD_CORNER,DD_EDGE,DD_FACE)
bddc_c() = (DD_CORNER,)
bddc_cf() = (DD_CORNER,DD_FACE)
bddc_ce() = (DD_CORNER,DD_EDGE)

function bddc_preconditioner(A::DDOperator;bddc_type=bddc_cef())
    dof_partition = partition(A.ids)
    coarse_dofs = generate_coarse_objects(bddc_type,dof_partition)
    setup_BDDC(A,coarse_dofs)
end

function generate_coarse_objects(
    bddc_type,
    dof_partition,
    ldof_to_colors=collect_parts_around(dof_partition),
    )
    ranks = linear_indices(dof_partition)
    state = map(ranks,ldof_to_colors) do rank,ldof_to_colors
        interface_dof_to_ldof = findall(colors->length(colors)!=1,ldof_to_colors)
        interface_dof_to_colors = view(ldof_to_colors,interface_dof_to_ldof)
        cldof_to_colors = unique(interface_dof_to_colors)
        cldof_to_owner = map(maximum,cldof_to_colors)
        ncldofs = length(cldof_to_colors)
        interface_dof_to_cldof = indexin(interface_dof_to_colors,cldof_to_colors)
        cldof_to_interface_dofs = inverse_index_map(interface_dof_to_cldof,ncldofs)
        f = interface_dof->interface_dof_to_ldof[interface_dof]
        cldof_to_interface_dofs.data .= f.(cldof_to_interface_dofs.data)
        cldof_to_ldofs = cldof_to_interface_dofs
        cldof_to_type = fill(DD_EDGE,ncldofs)
        cldof_to_type[ map(i->length(i) == 2,cldof_to_colors) ] .= DD_FACE
        cldof_to_type[ map(i->length(i) == 1,cldof_to_ldofs) ] .= DD_CORNER
        mask = map(i->i in bddc_type,cldof_to_type)
        cldof_to_colors = JaggedArray(cldof_to_colors[mask])
        cldof_to_ldofs = JaggedArray(cldof_to_ldofs[mask])
        cldof_to_owner = cldof_to_owner[mask]
        ncodofs = count(owner->owner==rank,cldof_to_owner)
        (;ldof_to_colors,cldof_to_colors,cldof_to_ldofs,cldof_to_owner,ncodofs,rank)
    end
    nown = map(i->i.ncodofs,state)
    ncdofs = sum(nown)
    codof_partition = variable_partition(nown,ncdofs)
    v = PVector{Vector{Int}}(undef,dof_partition)
    map(state,local_values(v),codof_partition) do state, ldof_to_cdof, codofs
        rank = state.rank
        cldof_to_owner = state.cldof_to_owner
        cldof_to_ldofs = state.cldof_to_ldofs
        codof_to_cdof = own_to_global(codofs)
        codof = 0
        for (cldof,owner) in enumerate(cldof_to_owner)
            if owner != rank
                continue
            end
            codof += 1
            cdof = codof_to_cdof[codof]
            ldofs = cldof_to_ldofs[cldof]
            ldof_to_cdof[ldofs] .= cdof
        end
    end
    consistent!(v) |> wait
    cldof_to_cdof_all = map(state,local_values(v)) do state, ldof_to_cdof
        rank = state.rank
        cldof_to_ldofs = state.cldof_to_ldofs
        ncldofs = length(cldof_to_ldofs)
        cldof_to_cdof = zeros(Int,ncldofs)
        for cldof in 1:ncldofs
            ldofs = cldof_to_ldofs[cldof]
            cdof = ldof_to_cdof[first(ldofs)]
            cldof_to_cdof[cldof] = cdof
        end
        cldof_to_cdof
    end
    rank_to_cdofs = gather(cldof_to_cdof_all,destination=MAIN)
    map(state,rank_to_cdofs) do state, rank_to_cdofs
        cldof_to_ldofs = state.cldof_to_ldofs
        (;cldof_to_ldofs, rank_to_cdofs,ncdofs, state...)
    end
end

function setup_BDDC(A,state)
    workspace1 = map(A.matrices,state) do Ai,state
        # Setup Neumann problems
        cldof_to_ldofs = state.cldof_to_ldofs
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
        ldof_to_colors = state.ldof_to_colors
        idof_to_ldof = findall(colors->length(colors)==1,ldof_to_colors)
        RIi = idof_to_ldof
        AIi = Ai[RIi,RIi]
        Pd = cholesky(AIi)

        # Setup weight
        ldof_to_w = map(colors->1/length(colors),ldof_to_colors)
        Wi = ldof_to_w

        rank_to_cdofs = state.rank_to_cdofs
        rank = state.rank
        ncdofs = state.ncdofs
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

function collect_parts_around(dof_partition)
    v = PVector{Vector{Int32}}(undef,dof_partition)
    map(i->fill!(i,Int32(1)),local_values(v))
    assemble!(v) |> wait
    consistent!(v) |> wait
    ranks = linear_indices(dof_partition)
    ldof_to_parts_around = map(ranks,local_values(v)) do rank,ldof_to_n_parts_around
        nldofs = length(ldof_to_n_parts_around)
        ptrs = zeros(Int32,nldofs+1)
        for ldof in 1:nldofs
            ptrs[ldof+1] = ldof_to_n_parts_around[ldof]
        end
        length_to_ptrs!(ptrs)
        ndata = ptrs[end]-1
        data = fill(Int32(rank),ndata)
        JaggedArray(data,ptrs)
    end
    map(fill!,local_values(v),ranks)
    cache = v.cache
    vector_partition = partition(v)
    buffer_snd = map(vector_partition,cache) do values,cache
        local_indices_snd = cache.local_indices_snd
        for (p,lid) in enumerate(local_indices_snd.data)
            cache.buffer_snd.data[p] = values[lid]
        end
        cache.buffer_snd
    end
    neighbors_snd, neighbors_rcv, buffer_rcv = map(cache) do cache
        cache.neighbors_snd, cache.neighbors_rcv, cache.buffer_rcv
    end |> tuple_of_arrays
    graph = ExchangeGraph(neighbors_snd,neighbors_rcv)
    exchange!(buffer_rcv,buffer_snd,graph) |> wait
    map(ldof_to_parts_around,vector_partition,cache) do ldof_to_parts_around,values,cache
        local_indices_rcv = cache.local_indices_rcv
        ldof_to_parts_around_ptrs = copy(ldof_to_parts_around.ptrs)
        for (p,lid) in enumerate(local_indices_rcv.data)
            q = ldof_to_parts_around_ptrs[lid]
            ldof_to_parts_around.data[q] = cache.buffer_rcv.data[p]
            ldof_to_parts_around_ptrs[lid] += 1
        end
        map(sort!,ldof_to_parts_around)
    end
    w = PVector(ldof_to_parts_around,dof_partition)
    consistent!(w) |> wait
    ldof_to_parts_around
end

function inverse_index_map(a_to_b,nb)
    b_to_as_ptrs = zeros(Int32,nb+1)
    for b in a_to_b
        b_to_as_ptrs[b+1] += 1
    end
    length_to_ptrs!(b_to_as_ptrs)
    ndata = b_to_as_ptrs[end]-1
    b_to_as_data = zeros(Int32,ndata)
    for (a,b) in enumerate(a_to_b)
        p = b_to_as_ptrs[b]
        b_to_as_data[p] = a
        b_to_as_ptrs[b] += 1
    end
    rewind_ptrs!(b_to_as_ptrs)
    JaggedArray(b_to_as_data,b_to_as_ptrs)
end


function monolithic_vector_assembly_strategy()
    function init(dofs,::Type{T}) where T
        n_total_dofs = length(dofs)
        offsets = [0]
        (;offsets,n_total_dofs,T_val=Val(T))
    end
    function init(block_dofs::BRange,::Type{T}) where T
        n_total_dofs = length(block_dofs)
        offsets = blocklasts(block_dofs) .- map(length,blocks(block_dofs))
        (;offsets,n_total_dofs,T_val=Val(T))
    end
    function scalar_type(setup)
        val_parameter(setup.T_val)
    end
    function counter(setup)
        0
    end
    function count(n,setup,i,field_i,free_or_diri)
        if free_or_diri == FREE && i>0
            n+=1
        end
        if free_or_diri == DIRICHLET && i<0
            n+=1
        end
        n
    end
    function allocate(n,setup)
        Ti = Int32
        T = val_parameter(setup.T_val)
        I = zeros(Ti,n)
        V = zeros(T,n)
        (;I,V)
    end
    function set!(alloc,n,setup,v,i,field,free_or_diri)
        if free_or_diri == FREE && i>0
            n+=1
            alloc.I[n] = i+setup.offsets[field]
            alloc.V[n] = v
        end
        if free_or_diri == DIRICHLET && i<0
            n+=1
            alloc.I[n] = setup.offsets[field]-i
            alloc.V[n] = v
        end
        n
    end
    function compress(alloc,setup)
        I = alloc.I
        V = alloc.V
        n_total_dofs = setup.n_total_dofs
        vec = PartitionedArrays.dense_vector(I,V,n_total_dofs)
        cache = nothing
        (vec, cache)
    end
    function compress!(vec,cache,alloc,setup)
        I = alloc.I
        V = alloc.V
        PartitionedArrays.dense_vector!(vec,I,V)
        vec
    end
    (;init,scalar_type,counter,count,allocate,set!,compress,compress!)
end

function monolithic_matrix_assembly_strategy(;matrix_type=nothing)
    function init(dofs_test,dofs_trial,::Type{T}) where T
        n_total_rows = length(dofs_test)
        n_total_cols = length(dofs_trial)
        offsets_rows = [0]
        offsets_cols = [0]
        (;offsets_rows,offsets_cols,n_total_rows,n_total_cols,T_val=Val(T))
    end
    function init(dofs_test::BRange,dofs_trial,::Type{T}) where T
        n_total_rows = length(dofs_test)
        n_total_cols = length(dofs_trial)
        offsets_rows = blocklasts(dofs_test) .- map(length,blocks(dofs_test))
        offsets_cols = blocklasts(dofs_trial) .- map(length,blocks(dofs_trial))
        (;offsets_rows,offsets_cols,n_total_rows,n_total_cols,T_val=Val(T))
    end
    function scalar_type(setup)
        val_parameter(setup.T_val)
    end
    function counter(setup)
        0
    end
    function count(n,setup,i,j,field_i,field_j,free_or_diri_i,free_or_diri_j)
        if free_or_diri_i == FREE && free_or_diri_j == FREE && i>0 && j>0
            n+=1
        end
        if free_or_diri_i == FREE && free_or_diri_j == DIRICHLET && i>0 && j<0
            n+=1
        end
        if free_or_diri_i == DIRICHLET && free_or_diri_j == FREE && i<0 && j>0
            n+=1
        end
        if free_or_diri_i == DIRICHLET && free_or_diri_j == DIRICHLET && i<0 && j<0
            n+=1
        end
        n
    end
    function allocate(n,setup)
        Ti = Int32
        T = val_parameter(setup.T_val)
        I = zeros(Ti,n)
        J = zeros(Ti,n)
        V = zeros(T,n)
        (;I,J,V)
    end
    function set!(alloc,n,setup,v,i,j,field_i,field_j,free_or_diri_i,free_or_diri_j)
        if free_or_diri_i == FREE && free_or_diri_j == FREE && i>0 && j>0
            n+=1
            alloc.I[n] = i+setup.offsets_rows[field_i]
            alloc.J[n] = j+setup.offsets_cols[field_j]
            alloc.V[n] = v
        end
        if free_or_diri_i == FREE && free_or_diri_j == DIRICHLET && i>0 && j<0
            n+=1
            alloc.I[n] = i+setup.offsets_rows[field_i]
            alloc.J[n] = setup.offsets_cols[field_j]-j
            alloc.V[n] = v
        end
        if free_or_diri_i == DIRICHLET && free_or_diri_j == FREE && i<0 && j>0
            n+=1
            alloc.I[n] = setup.offsets_rows[field_i]-i
            alloc.J[n] = j+setup.offsets_cols[field_j]
            alloc.V[n] = v
        end
        if free_or_diri_i == DIRICHLET && free_or_diri_j == DIRICHLET && i<0 && j<0
            n+=1
            alloc.I[n] = setup.offsets_rows[field_i]-i
            alloc.J[n] = setup.offsets_cols[field_j]-j
            alloc.V[n] = v
        end
        n
    end
    function compress(alloc,setup)
        I = alloc.I
        J = alloc.J
        V = alloc.V
        n_total_rows = setup.n_total_rows
        n_total_cols = setup.n_total_cols
        T = val_parameter(setup.T_val)
        M = matrix_type === nothing ? SparseMatrixCSC{T,Int} : matrix_type
        A,cache = sparse_matrix(M,I,J,V,n_total_rows,n_total_cols;reuse=Val(true))
        (A, cache)
    end
    function compress!(A,cache,alloc,setup)
        V = alloc.V
        sparse_matrix!(A,V,cache)
        A
    end
    (;init,scalar_type,counter,count,allocate,set!,compress,compress!)
end

function allocate_vector(
        ::Type{T},space_test::AbstractSpace,domains::AbstractDomain...;
        vector_strategy = monolithic_vector_assembly_strategy(),
        block_coupling = map(domain->fill(true,num_fields(space_test)),domains),
        free_or_dirichlet=FREE,
    ) where T
    setup = vector_strategy.init(GT.dofs(space_test,free_or_dirichlet),T)
    counter = vector_strategy.counter(setup)
    counter = increment(counter,domains,block_coupling) do counter1,domain,field_test_to_mask
        fields_test = ntuple(identity,Val(num_fields(space_test)))
        field_test_face_dofs = map(V->dofs_accessor(V,domain),fields(space_test))
        increment(counter1,fields_test,field_test_face_dofs) do counter2,field_test,face_dofs
            mask = field_test_to_mask[field_test]
            if mask
                nfaces = num_faces(domain)
                for face in 1:nfaces
                    dofs = face_dofs(face)
                    for dof_test in dofs
                        counter2 = vector_strategy.count(counter2,setup,dof_test,field_test,free_or_dirichlet)
                    end
                end
            end
            counter2
        end
    end
    coo = vector_strategy.allocate(counter,setup)
    counter_ref = Ref(vector_strategy.counter(setup))
    data = (;setup,coo,counter_ref,vector_strategy,free_or_dirichlet)
    VectorAllocation(data)
end

function increment(f,i,a::Tuple{Any},b::Tuple{Any})
    f(i,a[1],b[1])
end

function increment(f,i,a::Tuple{Any,Any},b::Tuple{Any,Any})
    i=f(i,a[1],b[1])
    f(i,a[2],b[2])
end

function increment(f,i,a::Tuple,b::Tuple)
    a1 = a[1]; b1 = b[1]
    at = Base.tail(a); bt = Base.tail(b)
    i=f(i,a[1],b[1])
    increment(f,i,at,bt)
end

struct VectorAllocation{A} <: AbstractType
    data::A
end

function Base.eltype(a::VectorAllocation)
    a.data.vector_strategy.scalar_type(a.data.setup)
end

function reset!(alloc::VectorAllocation)
    (;counter_ref,setup,vector_strategy) = alloc.data
    counter_ref[] = vector_strategy.counter(setup)
    alloc
end

function contribute!(alloc::VectorAllocation,b,dofs_test,field_test=1)
    (;setup,coo,counter_ref,vector_strategy,free_or_dirichlet) = alloc.data
    counter = counter_ref[]
    for (i,dof_test) in enumerate(dofs_test)
        counter = vector_strategy.set!(coo,counter,setup,b[i],dof_test,field_test,free_or_dirichlet)
    end
    counter_ref[] = counter
    alloc
end

function compress(alloc::VectorAllocation;reuse=Val(false))
    (;setup,coo,vector_strategy) = alloc.data
    b,cache = vector_strategy.compress(coo,setup)
    if val_parameter(reuse)
        b,cache
    else
        b
    end
end

function compress!(alloc::VectorAllocation,A,A_cache)
    (;setup,coo,vector_strategy) = alloc.data
    vector_strategy.compress!(A,A_cache,coo,setup)
    A
end

function allocate_matrix(
        ::Type{T},space_test::AbstractSpace,space_trial,domains::AbstractDomain...;
        matrix_strategy = monolithic_matrix_assembly_strategy(),
        block_coupling = map(domain->fill(true,num_fields(space_test),num_fields(space_trial)),domains),
        free_or_dirichlet = (FREE,FREE)
    ) where T
    setup = matrix_strategy.init(GT.dofs(space_test,free_or_dirichlet[1]),GT.dofs(space_trial,free_or_dirichlet[2]),T)
    counter = matrix_strategy.counter(setup)
    counter = increment(counter,domains,block_coupling) do counter1,domain,field_test_trial_to_mask
        field_test_face_dofs = map(V->dofs_accessor(V,domain),fields(space_test))
        field_trial_face_dofs = map(V->dofs_accessor(V,domain),fields(space_trial))
        fields_test = ntuple(identity,Val(num_fields(space_test)))
        fields_trial = ntuple(identity,Val(num_fields(space_trial)))
        increment(counter1,fields_trial,field_trial_face_dofs) do counter2,field_test,face_dofs_test
            increment(counter2,fields_trial,field_trial_face_dofs) do counter3,field_trial,face_dofs_trial
                nfaces = num_faces(domain)
                mask = field_test_trial_to_mask[field_test,field_trial]
                if mask
                    for face in 1:nfaces
                        dofs_test = face_dofs_test(face)
                        dofs_trial = face_dofs_trial(face)
                        for dof_test in dofs_test
                            for dof_trial in dofs_trial
                                counter3 = matrix_strategy.count(counter3,setup,dof_test,dof_trial,field_test,field_trial,free_or_dirichlet...)
                            end
                        end
                    end
                end
                counter3
            end
        end
    end
    coo = matrix_strategy.allocate(counter,setup)
    counter_ref = Ref(matrix_strategy.counter(setup))
    data = (;setup,coo,counter_ref,matrix_strategy,free_or_dirichlet)
    MatrixAllocation(data)
end

struct MatrixAllocation{A} <: AbstractType
    data::A
end

function Base.eltype(a::MatrixAllocation)
    a.data.matrix_strategy.scalar_type(a.data.setup)
end

function reset!(alloc::MatrixAllocation)
    (;counter_ref,setup,matrix_strategy) = alloc.data
    counter_ref[] = matrix_strategy.counter(setup)
    alloc
end

function contribute!(alloc::MatrixAllocation,b,dofs_test,dofs_trial,field_test=1,field_trial=1)
    (;setup,coo,counter_ref,matrix_strategy,free_or_dirichlet) = alloc.data
    counter = counter_ref[]
    for (i,dof_test) in enumerate(dofs_test)
        for (j,dof_trial) in enumerate(dofs_trial)
            counter = matrix_strategy.set!(coo,counter,setup,b[i,j],dof_test,dof_trial,field_test,field_trial,free_or_dirichlet...)
        end
    end
    counter_ref[] = counter
    alloc
end

function compress(alloc::MatrixAllocation;reuse=Val(false))
    (;setup,coo,matrix_strategy) = alloc.data
    b,cache = matrix_strategy.compress(coo,setup)
    if val_parameter(reuse)
        b,cache
    else
        b
    end
end

function compress!(alloc::MatrixAllocation,A,A_cache)
    (;setup,coo,matrix_strategy) = alloc.data
    matrix_strategy.compress!(A,A_cache,coo,setup)
    A
end

#function assemble_vector(f,::Type{T},space;kwargs...) where T
#    axis = 1
#    dv = GT.form_argument_quantity(space,axis)
#    integral = f(dv)
#    assemble_vector(integral,T,space;kwargs...)
#end
#
#function assemble_vector(integral::Number,::Type{T},space;
#    reuse = Val(false),
#    vector_strategy = monolithic_vector_assembly_strategy(),
#    ) where T
#    @assert integral == 0
#    free_dofs = GT.free_dofs(space)
#    setup = vector_strategy.init(free_dofs,T)
#    counter = vector_strategy.counter(setup)
#    alloc = vector_strategy.allocate(counter,setup)
#    vec,cache = vector_strategy.compress(alloc,setup)
#    if val_parameter(reuse) == false
#        vec
#    else
#        vec, (cache,alloc,setup,vector_strategy)
#    end
#end
#
#function assemble_vector(integral::Integral,::Type{T},space;
#    reuse = Val(false),
#    vector_strategy = monolithic_vector_assembly_strategy(),
#    free_or_dirichlet = FREE,
#    ) where T
#
#    setup = vector_strategy.init(GT.free_dofs(space),T)
#    state0 = (;space,vector_strategy,setup,free_or_dirichlet)
#    state1 = assemble_vector_count(integral,state0)
#    state2 = assemble_vector_allocate(state1)
#    state3 = assemble_vector_fill!(integral,state2)
#    b, vector_cache = assemble_vector_compress(state3)
#    cache = (;vector_cache,state2)
#    if val_parameter(reuse) == false
#        b
#    else
#        b, cache
#    end
#end
#
#function assemble_vector!(f,space,b,cache)
#    dim = 1
#    dv = GT.form_argument_quantity(space,dim)
#    integral = f(dv)
#    assemble_vector!(integral,b,cache)
#end
#
#function assemble_vector!(integral::Integral,b,cache)
#    (;vector_cache,state2) = cache
#    state3 = assemble_vector_fill!(integral,state2)
#    assemble_vector_compress!(b,vector_cache,state3)
#    b
#end
#
#function assemble_vector!(integral::Integral,V,b,cache)
#    (;vector_cache,state2) = cache
#    state3 = assemble_vector_fill!(integral,state2)
#    assemble_vector_compress!(b,vector_cache,state3)
#    b
#end
#
#function assemble_vector_count(integral,state)
#    (;space,vector_strategy,setup,free_or_dirichlet) = state
#    contributions = GT.contributions(integral)
#    nfields = GT.num_fields(space)
#    field_to_domain = map(field->domain(GT.field(space,field)),1:nfields)
#    field_to_Dface_to_dofs = map(field->face_dofs(GT.field(space,field)),1:nfields)
#    field_to_Dface_to_sDface = map(inverse_faces,field_to_domain) # TODO faces or inverse_faces??
#    field_to_D = map(num_dims,field_to_domain)
#    counter = vector_strategy.counter(setup)
#    for measure_and_contribution in contributions
#        measure, = measure_and_contribution
#        domain = GT.domain(measure)
#        sface_to_face = faces(domain)
#        topo = topology(mesh(domain))
#        d = num_dims(domain)
#        nsfaces = length(sface_to_face)
#        for sface in 1:nsfaces
#            face = sface_to_face[sface]
#            for field in 1:nfields
#                Dface_to_sDface = field_to_Dface_to_sDface[field]
#                D = field_to_D[field]
#                if D < d
#                    continue
#                end
#                face_to_Dfaces = face_incidence(topo,d,D)
#                Dface_to_dofs = field_to_Dface_to_dofs[field]
#                Dfaces = face_to_Dfaces[face]
#                for Dface in Dfaces
#                    sDface = Dface_to_sDface[Dface]
#                    if sDface == 0
#                        continue
#                    end
#                    dofs = Dface_to_dofs[Dface]
#                    for dof in dofs
#                        counter = vector_strategy.count(counter,setup,dof,field,free_or_dirichlet)
#                    end
#                end
#            end
#        end
#    end
#    (;counter,state...)
#end
#
#function assemble_vector_allocate(state)
#    (;counter,vector_strategy,setup) = state
#    alloc = vector_strategy.allocate(counter,setup)
#    (;alloc,state...)
#end
#
#function assemble_vector_fill!(integral,state)
#    (;space,vector_strategy,alloc,setup,free_or_dirichlet) = state
#    contributions = GT.contributions(integral)
#    nfields = GT.num_fields(space)
#    field_to_domain = map(field->domain(GT.field(space,field)),1:nfields)
#    field_to_Dface_to_dofs = map(field->face_dofs(GT.field(space,field)),1:nfields)
#    field_to_Dface_to_sDface = map(inverse_faces,field_to_domain)
#    field_to_D = map(num_dims,field_to_domain)
#    counter = vector_strategy.counter(setup)
#    T = vector_strategy.scalar_type(setup)
#    b = zeros(T,max_num_reference_dofs(space))
#    for measure_and_contribution in contributions
#        measure,contribution = measure_and_contribution
#        domain = GT.domain(measure)
#        sface_to_face = faces(domain)
#        topo = topology(mesh(domain))
#        d = num_dims(domain)
#        form_arity = 1
#        index = GT.generate_index(domain,form_arity)
#        t = GT.term(contribution,index)
#        term_npoints = GT.num_points(measure)
#        expr_qty = t |> expression |> simplify
#        expr_npoints = term_npoints(index) |> expression |> simplify
#        face = face_index(index,d)
#        point = point_index(index)
#        axis = 1
#        idof = dof_index(index,axis)
#        Dface_around = face_around_index(index,axis)
#        field = field_index(index,axis)
#        s_qty = GT.topological_sort(expr_qty,(face,field,Dface_around,point,idof))
#        s_npoints = GT.topological_sort(expr_npoints,(face,))
#        expr = quote
#            (counter,args,storage) -> begin
#                $(unpack_index_storage(index,:storage))
#                $(s_qty[1])
#                $(s_npoints[1])
#                b = args.b
#                nsfaces = length(args.sface_to_face)
#                for sface in 1:nsfaces
#                    $face = args.sface_to_face[sface]
#                    $(s_qty[2])
#                    npoints = $(s_npoints[2])
#                    for $field in 1:args.nfields
#                        $(s_qty[3])
#                        Dface_to_sDface = args.field_to_Dface_to_sDface[$field]
#                        D = args.field_to_D[$field]
#                        if D < args.d
#                            continue
#                        end
#                        face_to_Dfaces = face_incidence(args.topo,args.d,D)
#                        Dface_to_dofs = args.field_to_Dface_to_dofs[$field]
#                        Dfaces = face_to_Dfaces[$face]
#                        for ($Dface_around,Dface) in enumerate(Dfaces)
#                            $(s_qty[4])
#                            fill!(b,zero(eltype(b)))
#                            sDface = Dface_to_sDface[Dface]
#                            if sDface == 0
#                                continue
#                            end
#                            dofs = Dface_to_dofs[Dface]
#                            for $point in 1:npoints
#                                $(s_qty[5])
#                                for ($idof,dof) in enumerate(dofs)
#                                    b[$idof] += $(s_qty[6])
#                                end
#                            end
#                            for ($idof,dof) in enumerate(dofs)
#                                counter = args.vector_strategy.set!(args.alloc,counter,args.setup,b[$idof],dof,$field,args.free_or_dirichlet)
#                            end
#                        end
#                    end
#                end
#                counter
#            end
#        end
#        loop! = eval(expr)
#        storage = index_storage(index)
#        args = (;d,nfields,alloc,vector_strategy,setup,b,topo,sface_to_face,field_to_D,field_to_Dface_to_sDface,field_to_Dface_to_dofs,free_or_dirichlet)
#        counter = invokelatest(loop!,counter,args,storage)
#    end
#    state
#end
#
#function assemble_vector_compress(state)
#    (;alloc,setup,vector_strategy) = state
#    b, vector_cache = vector_strategy.compress(alloc,setup)
#    b, vector_cache
#end
#
#function assemble_vector_compress!(b,vector_cache,state)
#    (;alloc,setup,vector_strategy) = state
#    vector_strategy.compress!(b,vector_cache,alloc,setup)
#    b
#end
#
#function assemble_matrix(f,::Type{T},trial_space,test_space;kwargs...) where T
#    test_dim = 1
#    trial_dim = 2
#    dv = GT.form_argument_quantity(test_space,test_dim)
#    du = GT.form_argument_quantity(trial_space,trial_dim)
#    integral = f(du,dv)
#    assemble_matrix(integral,T,trial_space,test_space;kwargs...)
#end
#
#function assemble_matrix(integral::Integral,::Type{T},trial_space,test_space;
#        reuse=false,
#        matrix_strategy = monolithic_matrix_assembly_strategy(),
#        free_or_dirichlet = (FREE,FREE),
#    ) where T
#    setup = matrix_strategy.init(dofs(test_space,free_or_dirichlet[1]),dofs(trial_space,free_or_dirichlet[2]),T)
#    state0 = (;test_space,trial_space,matrix_strategy,setup,free_or_dirichlet)
#    state1 = assemble_matrix_count(integral,state0)
#    state2 = assemble_matrix_allocate(state1)
#    state3 = assemble_matrix_fill!(integral,state2)
#    A, matrix_cache = assemble_matrix_compress(state3)
#    cache = (;matrix_cache,state2)
#    if reuse == false
#        A
#    else
#        A, cache
#    end
#end
#
#function assemble_matrix!(f,trial_space,test_space,A,cache)
#    test_dim = 1
#    trial_dim = 2
#    dv = GT.form_argument_quantity(test_space,test_dim)
#    du = GT.form_argument_quantity(trial_space,trial_dim)
#    integral = f(du,dv)
#    assemble_matrix!(integral,A,cache)
#end
#
#function assemble_matrix!(integral::Integral,A,cache)
#    (;matrix_cache,state2) = cache
#    state3 = assemble_matrix_fill!(integral,state2)
#    assemble_matrix_compress!(A,matrix_cache,state3)
#    A
#end
#
#function assemble_matrix!(integral::Integral,U,V,A,cache)
#    (;matrix_cache,state2) = cache
#    state3 = assemble_matrix_fill!(integral,state2)
#    assemble_matrix_compress!(A,matrix_cache,state3)
#    A
#end
#
#function assemble_matrix_count(integral,state)
#    (;test_space,trial_space,matrix_strategy,setup,free_or_dirichlet) = state
#    contributions = GT.contributions(integral)
#    test, trial = 1, 2
#    axis_to_space = (test_space,trial_space)
#    axis_to_nfields = map(num_fields,axis_to_space)
#    axis_to_field_to_domain = map(space->map(field->domain(GT.field(space,field)),1:num_fields(space)),axis_to_space)
#    axis_to_field_to_Dface_to_dofs = map(space->map(field->face_dofs(GT.field(space,field)),1:num_fields(space)),axis_to_space)
#    axis_to_field_to_Dface_to_sDface = map(field_to_domain->map(inverse_faces,field_to_domain),axis_to_field_to_domain) # TODO: faces or inverse_faces??
#    axis_to_field_to_D = map(field_to_domain->map(num_dims,field_to_domain),axis_to_field_to_domain)
#    counter = matrix_strategy.counter(setup)
#    for measure_and_contribution in contributions
#        measure, = measure_and_contribution
#        domain = GT.domain(measure)
#        sface_to_face = faces(domain)
#        topo = topology(mesh(domain))
#        d = num_dims(domain)
#        nsfaces = length(sface_to_face)
#        for sface in 1:nsfaces
#            face = sface_to_face[sface]
#            for field_test in 1:axis_to_nfields[test]
#                test_Dface_to_sDface = axis_to_field_to_Dface_to_sDface[test][field_test]
#                test_D = axis_to_field_to_D[test][field_test]
#                if test_D < d
#                    continue
#                end
#                test_face_to_Dfaces = face_incidence(topo,d,test_D)
#                test_Dface_to_dofs = axis_to_field_to_Dface_to_dofs[test][field_test]
#                test_Dfaces = test_face_to_Dfaces[face]
#                for field_trial in 1:axis_to_nfields[trial]
#                    trial_Dface_to_sDface = axis_to_field_to_Dface_to_sDface[trial][field_trial]
#                    trial_D = axis_to_field_to_D[trial][field_trial]
#                    if trial_D < d
#                        continue
#                    end
#                    trial_face_to_Dfaces = face_incidence(topo,d,trial_D)
#                    trial_Dface_to_dofs = axis_to_field_to_Dface_to_dofs[trial][field_trial]
#                    trial_Dfaces = trial_face_to_Dfaces[face]
#                    for test_Dface in test_Dfaces 
#                        test_sDface = test_Dface_to_sDface[test_Dface]
#                        if test_sDface == 0
#                            continue
#                        end
#                        test_dofs = test_Dface_to_dofs[test_Dface]
#                        for trial_Dface in trial_Dfaces
#                            trial_sDface = trial_Dface_to_sDface[trial_Dface]
#                            if trial_sDface == 0
#                                continue
#                            end
#                            trial_dofs = trial_Dface_to_dofs[trial_Dface]
#                            for test_dof in test_dofs
#                                for trial_dof in trial_dofs
#                                    counter = matrix_strategy.count(counter,setup,test_dof,trial_dof,field_test,field_trial,free_or_dirichlet...)
#                                end
#                            end
#                        end
#                    end
#                end
#            end
#        end
#    end
#    (;counter,state...)
#end
#
#function assemble_matrix_allocate(state)
#    (;matrix_strategy,counter,setup) = state
#    alloc = matrix_strategy.allocate(counter,setup)
#    (;alloc,state...)
#end
#
#function assemble_matrix_fill!(integral,state)
#    (;test_space,trial_space,alloc,matrix_strategy,setup,free_or_dirichlet) = state 
#
#    contributions = GT.contributions(integral)
#    test, trial = 1, 2
#    axis_to_space = (test_space,trial_space)
#    axis_to_nfields = map(num_fields,axis_to_space)
#    axis_to_field_to_domain = map(space->map(field->domain(GT.field(space,field)),1:num_fields(space)),axis_to_space)
#    axis_to_field_to_Dface_to_dofs = map(space->map(field->face_dofs(GT.field(space,field)),1:num_fields(space)),axis_to_space)
#    axis_to_field_to_Dface_to_sDface = map(field_to_domain->map(inverse_faces,field_to_domain),axis_to_field_to_domain) # TODO: faces or inverse_faces??
#    axis_to_field_to_D = map(field_to_domain->map(num_dims,field_to_domain),axis_to_field_to_domain)
#    counter = matrix_strategy.counter(setup)
#    T = matrix_strategy.scalar_type(setup)
#    b = zeros(T,max_num_reference_dofs(trial_space), max_num_reference_dofs(test_space))
#
#    for measure_and_contribution in contributions
#        measure,contribution = measure_and_contribution
#        domain = GT.domain(measure)
#        sface_to_face = faces(domain)
#        topo = topology(mesh(domain))
#        d = num_dims(domain)
#        form_arity = 2
#        index = GT.generate_index(domain,form_arity)
#        t = GT.term(contribution,index)
#        term_npoints = GT.num_points(measure)
#        expr_qty = t |> expression |> simplify
#        expr_npoints = term_npoints(index) |> expression |> simplify
#        face = face_index(index,d)
#        point = point_index(index)
#        axes = (test, trial)
#
#        idof_test, idof_trial = map(axis -> dof_index(index,axis), axes) # TODO: test, trial
#        test_Dface_around, trial_Dface_around = map(axis -> face_around_index(index,axis), axes)
#        field_test, field_trial = map(axis -> field_index(index,axis), axes)
#
#        s_qty = GT.topological_sort(expr_qty,(face, field_test, field_trial, test_Dface_around, trial_Dface_around, point, idof_test, idof_trial)) # 8 deps
#        s_npoints = GT.topological_sort(expr_npoints,(face,))
#        
#        expr = quote
#            function loop!(counter,args,storage)
#                $(unpack_index_storage(index,:storage))
#                $(s_qty[1])
#                $(s_npoints[1])
#                b = args.b
#                nsfaces = length(args.sface_to_face)
#                for sface in 1:nsfaces
#                    $face = args.sface_to_face[sface]
#                    $(s_qty[2])
#                    npoints = $(s_npoints[2])
#                    for $field_test in 1:args.axis_to_nfields[$test]
#                        $(s_qty[3])
#                        test_Dface_to_sDface = args.axis_to_field_to_Dface_to_sDface[$test][$field_test]
#                        test_D = args.axis_to_field_to_D[$test][$field_test]
#                        if test_D < args.d
#                            continue
#                        end
#                        test_face_to_Dfaces = face_incidence(args.topo,args.d,test_D)
#                        test_Dface_to_dofs = args.axis_to_field_to_Dface_to_dofs[$test][$field_test]
#                        test_Dfaces = test_face_to_Dfaces[$face]
#                        for $field_trial in 1:args.axis_to_nfields[$trial]
#                            $(s_qty[4])
#                            trial_Dface_to_sDface = args.axis_to_field_to_Dface_to_sDface[$trial][$field_trial]
#                            trial_D = args.axis_to_field_to_D[$trial][$field_trial]
#                            if trial_D < args.d
#                                continue
#                            end
#                            trial_face_to_Dfaces = face_incidence(args.topo,args.d,trial_D)
#                            trial_Dface_to_dofs = args.axis_to_field_to_Dface_to_dofs[$trial][$field_trial]
#                            trial_Dfaces = trial_face_to_Dfaces[$face]
#                            for ($test_Dface_around,test_Dface) in enumerate(test_Dfaces)
#                                $(s_qty[5])
#                                test_sDface = test_Dface_to_sDface[test_Dface]
#                                if test_sDface == 0
#                                    continue
#                                end
#                                test_dofs = test_Dface_to_dofs[test_Dface]
#                                for ($trial_Dface_around,trial_Dface) in enumerate(trial_Dfaces)
#                                    $(s_qty[6])
#                                    trial_sDface = trial_Dface_to_sDface[trial_Dface]
#                                    if trial_sDface == 0
#                                        continue
#                                    end
#                                    trial_dofs = trial_Dface_to_dofs[trial_Dface]
#                                    fill!(b,zero(eltype(b)))
#                                    for $point in 1:npoints
#                                        $(s_qty[7])
#                                        for ($idof_test, dof_test) in enumerate(test_dofs)
#                                            $(s_qty[8])
#                                            for ($idof_trial, dof_trial) in enumerate(trial_dofs)
#                                                b[$idof_test, $idof_trial] += $(s_qty[9])
#                                            end
#                                        end
#                                    end
#                                    for ($idof_test, dof_test) in enumerate(test_dofs)
#                                        for ($idof_trial, dof_trial) in enumerate(trial_dofs)
#                                            counter = args.matrix_strategy.set!(args.alloc,counter,args.setup,b[$idof_test, $idof_trial],dof_test,dof_trial,$field_test,$field_trial,args.free_or_dirichlet...)
#                                        end
#                                    end
#                                end
#                            end
#                        end
#                    end
# 
#                end
#                counter
#            end
#        end
#        loop! = eval(expr)
#        storage = index_storage(index)
#        args = (;d,axis_to_nfields,alloc,matrix_strategy,setup,b,topo,sface_to_face,axis_to_field_to_D,axis_to_field_to_Dface_to_sDface,axis_to_field_to_Dface_to_dofs,free_or_dirichlet) 
#        counter = invokelatest(loop!,counter,args,storage)
#    end
#    state
#end
#
#function assemble_matrix_compress(state)
#    (;alloc,matrix_strategy,setup) = state
#    A, matrix_cache = matrix_strategy.compress(alloc,setup)
#    A, matrix_cache
#end
#
#function assemble_matrix_compress!(A,matrix_cache,state)
#    (;alloc,matrix_strategy,setup) = state
#    matrix_strategy.compress!(A,matrix_cache,alloc,setup)
#    A
#end
#
#function assemble_matrix_and_vector(a,l,::Type{T},U,V;
#        reuse = Val(false),
#        matrix_strategy = monolithic_matrix_assembly_strategy(),
#        vector_strategy = monolithic_vector_assembly_strategy(),
#    ) where T
#    A,Acache = assemble_matrix(a,T,U,V;reuse=Val(true),matrix_strategy)
#    b,bcache = assemble_vector(l,T,V;reuse=Val(true),vector_strategy)
#    cache = (;Acache,bcache)
#    if val_parameter(reuse)
#        A,b,cache
#    else
#        A,b
#    end
#end
#
#function assemble_matrix_and_vector!(a,l,U,V,A,b,cache)
#    (;matrix_cache,vector_cache) = cache
#    assemble_matrix!(a,U,V,A,matrix_cache)
#    assemble_vector!(l,V,b,vector_cache)
#    A,b
#end
#
#function assemble_matrix_and_vector_with_dirichlet(a,l,U,V,dirichlet_values;kwargs...)
#    T = eltype(dirichlet_values)
#    free_values = constant_values(zero(T),GT.free_dofs(U))
#    ud = discrete_field(U,free_values,dirichlet_values)
#    l2(v) = l(v) - a(ud,v)
#    assemble_matrix_and_vector(a,l2,T,U,V;kwargs...)
#end
#
#function assemble_matrix_and_vector_with_dirichlet!(a,l,U,V,dirichlet_values,A,b,cache)
#    T = eltype(diri_vals)
#    free_values = constant_values(zero(T),GT.free_dofs(U))
#    ud = discrete_field(U,free_values,dirichlet_values)
#    l2(v) = l(v) - a(ud,v)
#    assemble_matrix_and_vector!(a,l2,U,V,A,b,cache)
#end
#
#function linear_problem(uhd::DiscreteField,a,l,V=GT.space(uhd);
#        matrix_strategy = monolithic_matrix_assembly_strategy(),
#        vector_strategy = monolithic_vector_assembly_strategy(),
#    )
#    U = GT.space(uhd)
#    A,b = assemble_matrix_and_vector_with_dirichlet(a,l,U,V,dirichlet_values(uhd);matrix_strategy,vector_strategy)
#    x = similar(b,axes(A,2))
#    PS.linear_problem(x,A,b)
#    #T = eltype(dirichlet_values(uhd))
#    #A = assemble_matrix(a,T,U,V)
#    #b = assemble_vector(l,T,V)
#    #d = assemble_vector(v->a(uhd,v),T,V)
#    #b .= b .- d
#    #x = similar(b,axes(A,2))
#    #PS.linear_problem(x,A,b)
#end
#
#function linear_problem(::Type{T},U::AbstractSpace,a,l,V=U;
#        matrix_strategy = monolithic_matrix_assembly_strategy(),
#        vector_strategy = monolithic_vector_assembly_strategy(),
#    ) where T
#    A,b = assemble_matrix_and_vector(a,l,T,U,V;matrix_strategy,vector_strategy)
#    x = similar(b,axes(A,2))
#    PS.linear_problem(x,A,b)
#end
#
#function solution_field(U::AbstractSpace,x::AbstractVector)
#    T = eltype(x)
#    uhd = zero_dirichlet_field(T,U)
#    solution_field(uhd,x)
#end
#
#function solution_field(uhd::DiscreteField,x::AbstractVector)
#    diri_vals = dirichlet_values(uhd)
#    U = GT.space(uhd)
#    free_vals = free_values_from_solution(x,free_dofs(U))
#    discrete_field(U,free_vals,diri_vals)
#end
#
#function solution_field!(uh::DiscreteField,x::AbstractVector)
#    U = GT.space(uh)
#    free_vals_src = free_values_from_solution(x,free_dofs(U))
#    free_vals_dest = free_values(uh)
#    copyto!(free_vals_dest,free_vals_src)
#    uh
#end
#
#function solution_field(U::AbstractSpace,p::PS.AbstractProblem)
#    solution_field(U,PS.solution(p))
#end
#
#function solution_field(U::DiscreteField,p::PS.AbstractProblem)
#    solution_field(U,PS.solution(p))
#end
#
#function solution_field(U::AbstractSpace,p::PS.AbstractSolver)
#    solution_field(U,PS.solution(p))
#end
#
#function solution_field(U::DiscreteField,p::PS.AbstractSolver)
#    solution_field(U,PS.solution(p))
#end
#
#function solution_field(uts,p::PS.AbstractODEProblem)
#    x = PS.solution(p)
#    t = x[1]
#    uas = x[2:end]
#    uhs = map((ua,ut)->solution_field(ut(t),ua),uas,uts)
#    t, uhs
#end
#
#function solution_field(uts,p::PS.AbstractODESolver)
#    solution_field(uts,PS.problem(p))
#end
#
#function solution_field!(uh::DiscreteField,p::PS.AbstractProblem)
#    solution_field!(uh,PS.solution(p))
#end
#
#function solution_field!(uh::DiscreteField,p::PS.AbstractSolver)
#    solution_field!(uh,PS.solution(p))
#end
#
#function free_values_from_solution(x,dofs)
#    x
#end
#
#function free_values_from_solution(x,dofs::BRange)
#    nfields = blocklength(dofs)
#    map(1:nfields) do field
#        pend = blocklasts(dofs)[field]
#        pini = 1 + pend - length(blocks(dofs)[field])
#        view(x,pini:pend)
#    end |> BVector
#end
#
#function free_values_from_solution(x::BVector,dofs::BRange)
#    x
#end
#
#function nonlinear_problem(uh0::DiscreteField,r,j,V=GT.space(uh0);
#        matrix_strategy = monolithic_matrix_assembly_strategy(),
#        vector_strategy = monolithic_vector_assembly_strategy(),
#    )
#    a0 = j(uh0)
#    l0 = r(uh0)
#    U = GT.space(uh0)
#    x0 = free_values(uh0)
#    T = eltype(x0)
#    A0,b0,cache = assemble_matrix_and_vector(a0,l0,T,U,V;
#       reuse=Val(true),matrix_strategy,vector_strategy)
#    PS.nonlinear_problem(x0,b0,A0,cache) do p
#        x = PS.solution(p)
#        # TODO call to consistent! in parallel code
#        uh = solution_field(uh0,x)
#        a = j(uh)
#        l = r(uh)
#        A = PS.jacobian(p)
#        b = PS.residual(p)
#        ws = PS.workspace(p)
#        if b !== nothing && A !== nothing
#            assemble_matrix_and_vector!(a,l,U,V,A,b,ws)
#        elseif b !== nothing && A === nothing
#            assemble_vector!(l,V,b,ws.vector_cache)
#        elseif b === nothing && A !== nothing
#            assemble_matrix!(a,U,V,A,ws.matrix_cache)
#        end
#        p = PS.update(p,solution=x,residual=b,jacobian=A)
#    end
#end

#function semi_discrete_field(T,V::AbstractSpace)
#    semi_discrete_field(T,V) do t,uh
#        uh
#    end
#end
#
#function semi_discrete_field(f,T,V::AbstractSpace)
#    uh = zero_field(T,V)
#    semi_discrete_field(f,uh)
#end
#
#function semi_discrete_field(uh::DiscreteField)
#    semi_discrete_field(uh) do t,uh
#        uh
#    end
#end
#
#function semi_discrete_field(f,uh::DiscreteField)
#    SemiDiscreteField(f,uh)
#end
#
#struct SemiDiscreteField{A,B}
#    update::A
#    discrete_field::B
#end
#
#free_values(uh::SemiDiscreteField) = free_values(uh.discrete_field)
#dirichlet_values(uh::SemiDiscreteField) = dirichlet_values(uh.discrete_field)
#
#function (u::SemiDiscreteField)(t)
#    u.update(t,u.discrete_field)
#    u.discrete_field
#end
#
#function space(u::SemiDiscreteField)
#    space(u.discrete_field)
#end
#
#function discrete_field(u::SemiDiscreteField)
#    u.discrete_field
#end
#
#function nonlinear_ode(
#    uts,tspan,res,jacs,V=space(uts[1]);
#    matrix_strategy = monolithic_matrix_assembly_strategy(),
#    vector_strategy = monolithic_vector_assembly_strategy(),
#    )
#
#    U = space(uts[1])
#    test=1
#    trial=2
#    v = form_argument_quantity(V,test)
#    du = form_argument_quantity(U,trial)
#    t0 = first(tspan)
#    uhs_0 = map(uh->uh(t0),uts)
#    coeffs_0 = map(uh_0 -> one(eltype(free_values(uh_0))),uhs_0) 
#    r0_int = res(t0,uhs_0)(v)
#    js0_int = map((jac,coeff)->coeff*jac(t0,uhs_0)(du,v),jacs,coeffs_0)
#    j0_int = sum(js0_int)
#    T = eltype(free_values(uhs_0[1]))
#    A0,b0,cache = assemble_matrix_and_vector(j0_int,r0_int,T,U,V;
#       reuse=Val(true),matrix_strategy,vector_strategy)
#    x0 = (t0,map(free_values,uhs_0)...)
#    PS.ode_problem(x0,b0,A0,tspan,coeffs_0,cache) do p
#        coeffs = PS.coefficients(p)
#        x = PS.solution(p)
#        A = PS.jacobian(p)
#        b = PS.residual(p)
#        ws = PS.workspace(p)
#        t = x[1]
#        uhs_t = map((uh,y)->solution_field(uh(t),y),uts,x[2:end])
#        r_int = res(t,uhs_t)(v)
#        js_int = map((jac,coeff)->coeff*jac(t,uhs_t)(du,v),jacs,coeffs)
#        j_int = sum(js_int)
#        if b !== nothing && A !== nothing
#            assemble_matrix_and_vector!(j_int,r_int,U,V,A,b,ws)
#        elseif b !== nothing && A === nothing
#            assemble_vector!(r_int,V,b,ws.vector_cache)
#        elseif b === nothing && A !== nothing
#            assemble_matrix!(j_int,U,V,A,ws.matrix_cache)
#        end
#        p = PS.update(p,solution=x,residual=b,jacobian=A)
#    end
#end
#
### TODO It is being used currently?
##function vector_allocation(vector_strategy,space_test,domains...)
##    b_setup = vector_strategy.init(GT.free_dofs(space_test),T)
##    b_counter = vector_strategy.counter(b_setup)
##    for domain in domains
##        for field in 1:nfields
##            for sface in 1:nsfaces
##                for face_around in 1:nfaces_around
##                    for dof in dofs
##                        b_counter = vector_strategy.count(b_counter,b_setup,dof_test,field_test)
##                    end
##                end
##            end
##        end
##    end
##    b_alloc = vector_strategy.allocate(b_counter,b_setup)
##end
#
##function allocate_vector(
##        ::Type{T},space_test::AbstractSpace,domains::AbstractDomain...;
##        vector_strategy = monolithic_vector_assembly_strategy(),
##    ) where T
##    setup = vector_strategy.init(GT.free_dofs(space_test),T)
##    counter = vector_strategy.counter(setup)
##    @assert num_fields(space_test) == 1
##    field_test = 1
##    for domain in domains
##        face_dofs = dofs_accessor(space_test,domain)
##        nfaces = num_faces(domain)
##        for face in 1:nfaces
##            dofs = face_dofs(face)
##            for dof_test in dofs
##                counter = vector_strategy.count(counter,setup,dof_test,field_test)
##            end
##        end
##    end
##    coo = vector_strategy.allocate(counter,setup)
##    counter_ref = Ref(vector_strategy.counter(setup))
##    data = (;setup,coo,counter_ref,vector_strategy)
##    VectorAllocation(data)
##end
#
##function allocate_vector(
##        ::Type{T},space_test::AbstractSpace,domains::AbstractDomain...;
##        vector_strategy = monolithic_vector_assembly_strategy(),
##        block_coupling = map(domain->fill(true,num_fields(space_test)),domains),
##    ) where T
##    setup = vector_strategy.init(GT.free_dofs(space_test),T)
##    counter_ref = Ref(vector_strategy.counter(setup))
##    fields_test = ntuple(identity,Val(num_fields(space_test)))
##    let domain = domains[1], field_test_to_mask = block_coupling[1]
##    #map(block_coupling,domains) do field_test_to_mask,domain
##        field_test_face_dofs = map(V->dofs_accessor(V,domain),fields(space_test))
##        nfaces = num_faces(domain)
##        let field_test = fields_test[1], face_dofs = field_test_face_dofs[1]
##        #map(fields_test,field_test_face_dofs) do field_test,face_dofs
##            mask = field_test_to_mask[field_test]
##            if mask
##                counter = counter_ref[]
##                for face in 1:nfaces
##                    dofs = face_dofs(face)
##                    for dof_test in dofs
##                        counter = vector_strategy.count(counter,setup,dof_test,field_test)
##                    end
##                end
##                counter_ref[] = counter
##            end
##        end
##    end
##    coo = vector_strategy.allocate(counter_ref[],setup)
##    counter_ref = Ref(vector_strategy.counter(setup))
##    data = (;setup,coo,counter_ref,vector_strategy)
##    VectorAllocation(data)
##end


#function allocate_vector_domains(domains::Tuple{Any},block_coupling,counter,setup,vector_strategy,space_test)
#    domain = domains[1]
#    field_test_to_mask = block_coupling[1]
#    allocate_vector_domain(domain,field_test_to_mask,counter,setup,vector_strategy,space_test)
#end

#function allocate_vector_domain(domain,field_test_to_mask,counter,setup,vector_strategy,space_test)
#    fields_test = ntuple(identity,Val(num_fields(space_test)))
#    field_test_face_dofs = map(V->dofs_accessor(V,domain),fields(space_test))
#    nfaces = num_faces(domain)
#    let field_test = fields_test[1], face_dofs = field_test_face_dofs[1]
#        #map(fields_test,field_test_face_dofs) do field_test,face_dofs
#        mask = field_test_to_mask[field_test]
#        if mask
#            for face in 1:nfaces
#                dofs = face_dofs(face)
#                for dof_test in dofs
#                    counter = vector_strategy.count(counter,setup,dof_test,field_test)
#                end
#            end
#        end
#    end
#    counter
#end


#function vector_assembly_loop!(fs,alloc,faces)
#    for face in faces
#        map(fs) do f
#            b,dofs,field = f(face)
#            contribute!(alloc,b,dofs,field)
#            nothing
#        end
#    end
#    alloc
#end


#function allocate_matrix(
#        ::Type{T},space_test::AbstractSpace,space_trial,domains::AbstractDomain...;
#        matrix_strategy = monolithic_matrix_assembly_strategy(),
#    ) where T
#    setup = matrix_strategy.init(GT.free_dofs(space_test),GT.free_dofs(space_trial),T)
#    counter = matrix_strategy.counter(setup)
#    @assert num_fields(space_test) == 1
#    @assert num_fields(space_trial) == 1
#    field_test = 1
#    field_trial = 1
#    for domain in domains
#        face_dofs_test = dofs_accessor(space_test,domain)
#        face_dofs_trial = dofs_accessor(space_trial,domain)
#        nfaces = num_faces(domain)
#        for face in 1:nfaces
#            dofs_test = face_dofs_test(face)
#            dofs_trial = face_dofs_trial(face)
#            for dof_test in dofs_test
#                for dof_trial in dofs_trial
#                    counter = matrix_strategy.count(counter,setup,dof_test,dof_trial,field_test,field_trial)
#                end
#            end
#        end
#    end
#    coo = matrix_strategy.allocate(counter,setup)
#    counter_ref = Ref(matrix_strategy.counter(setup))
#    data = (;setup,coo,counter_ref,matrix_strategy)
#    MatrixAllocation(data)
#end


#function matrix_assembly_loop!(fs,alloc,faces)
#
#
#function matrix_assembly_loop!(fs,alloc,faces)
#    for face in faces
#        map(fs) do f
#            A,dofs_i,dofs_j,field_i,field_j = f(face)
#            contribute!(alloc,A,dofs_i,dofs_j,field_i,field_j)
#            nothing
#        end
#    end
#    alloc
#end




#function system_assembler(uhd::DiscreteField,space_trial,space_test,domains...;
#        vector_strategy=monolithic_vector_assembly_strategy(),
#        matrix_strategy=monolithic_matrix_assembly_strategy(),
#    )
#    T = eltype(free_values(uh))
#    b_setup = vector_strategy.init(GT.free_dofs(space_test),T)
#    A_setup = matrix_strategy.init(free_dofs(space_test),free_dofs(space_trial),T)
#    b_counter = vector_strategy.counter(b_setup)
#    A_counter = matrix_strategy.counter(A_setup)
#    rface_to_dofs_test = face_dofs(space_test)
#    rface_to_dofs_trial = face_dofs(space_trial)
#    face_to_rface_test = inverse_faces(GT.domain(space_test))
#    face_to_rface_trial = inverse_faces(GT.domain(space_trial))
#    @assert num_fields(space_test) == 1
#    @assert num_fields(space_trial) == 1
#    field_test = 1 # TODO
#    field_trial = 1 # TODO
#    for domain in domains
#        sface_to_face = faces(domain)
#        for face in sface_to_face
#            rface_test = face_to_rface_test[face]
#            rface_trial = face_to_rface_trial[face]
#            dofs_test = rface_to_dofs_test[rface_test]
#            dofs_trial = rface_to_dofs_trial[rface_trial]
#            for dof_test in dofs_test
#                b_counter = vector_strategy.count(b_counter,b_setup,dof_test,field_test)
#                for dof_trial in dofs_trial
#                    A_counter = matrix_strategy.count(A_counter,A_setup,dof_test,dof_trial,field_test,field_trial)
#                end
#            end
#        end
#    end
#    b_alloc = vector_strategy.allocate(b_counter,b_setup)
#    A_alloc = matrix_strategy.allocate(A_counter,A_setup)
#    b_counter_ref = vector_strategy.counter(b_setup) |> Ref
#    A_counter_ref = matrix_strategy.counter(A_setup) |> Ref
#    map(domains) do domain
#        sface_to_face = faces(domain)
#        function set_domain!(sface,A,b)
#            face = sface_to_face[sface]
#            rface_test = face_to_rface_test[face]
#            rface_trial = face_to_rface_trial[face]
#            dofs_test = rface_to_dofs_test[rface_test]
#            dofs_trial = rface_to_dofs_trial[rface_trial]
#
#            for dof_test in dofs_test
#                b_counter = vector_strategy.count(b_counter,b_setup,dof_test,field_test)
#                for dof_trial in dofs_trial
#                    A_counter = matrix_strategy.count(A_counter,A_setup,dof_test,dof_trial,field_test,field_trial)
#                end
#            end
#        end
#    end
#
#
#end

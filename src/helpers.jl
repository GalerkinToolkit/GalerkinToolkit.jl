# Helpers for hand-written assembly

function coordinate_accessor(measure::Measure)
    mesh = measure.mesh
    dom = measure.domain
    @assert !is_reference_domain(dom)
    d = num_dims(dom)
    rid_to_point_to_x = map(coordinates,reference_quadratures(measure))
    rid_to_tab = map(rid_to_point_to_x,reference_faces(mesh,d)) do point_to_x, refface
        tabulator(refface)(value,point_to_x)
    end
    face_to_rid = face_reference_id(mesh,d)
    face_to_nodes = face_nodes(mesh,d)
    node_to_x = node_coordinates(mesh)
    sface_to_face = faces(dom)
    function face_point_x(sface)
        face = sface_to_face[sface]
        rid = face_to_rid[face]
        tab = rid_to_tab[rid]
        nodes = face_to_nodes[face]
        function point_x(point)
            nnodes = length(nodes)
            sum(1:nnodes) do i
                node = nodes[i]
                x = node_to_x[node]
                x*tab[point,i]
            end
        end
    end
end

outer(a,b) = a*transpose(b)

function jacobian_accessor(measure::Measure)
    mesh = measure.mesh
    dom = measure.domain
    @assert !is_reference_domain(dom)
    d = num_dims(dom)
    rid_to_point_to_x = map(coordinates,reference_quadratures(measure))
    rid_to_tab = map(rid_to_point_to_x,reference_faces(mesh,d)) do point_to_x, refface
        tabulator(refface)(ForwardDiff.gradient,point_to_x)
    end
    face_to_rid = face_reference_id(mesh,d)
    face_to_nodes = face_nodes(mesh,d)
    node_to_x = node_coordinates(mesh)
    sface_to_face = faces(dom)
    function face_point_x(sface)
        face = sface_to_face[sface]
        rid = face_to_rid[face]
        tab = rid_to_tab[rid]
        nodes = face_to_nodes[face]
        function point_x(point)
            nnodes = length(nodes)
            sum(1:nnodes) do i
                node = nodes[i]
                x = node_to_x[node]
                outer(x,tab[point,i])
            end
        end
    end
end

function num_points_accessor(measure::Measure)
    mesh = measure.mesh
    dom = measure.domain
    d = num_dims(dom)
    rid_to_point_to_w = map(weights,reference_quadratures(measure))
    face_to_rid = face_reference_id(mesh,d)
    sface_to_face = faces(dom)
    function face_npoints(sface)
        face = sface_to_face[sface]
        rid = face_to_rid[face]
        point_to_w = rid_to_point_to_w[rid]
        length(point_to_w)
    end
end

function weight_accessor(measure::Measure)
    mesh = measure.mesh
    dom = measure.domain
    @assert !is_reference_domain(dom)
    d = num_dims(dom)
    rid_to_point_to_w = map(weights,reference_quadratures(measure))
    face_to_rid = face_reference_id(mesh,d)
    sface_to_face = faces(dom)
    function face_point_x(sface)
        face = sface_to_face[sface]
        rid = face_to_rid[face]
        point_to_w = rid_to_point_to_w[rid]
        function point_x(point,J)
            w = point_to_w[point]
            change_of_measure(J)*w
        end
    end
end

function shape_function_accessor(f::typeof(value),space::AbstractSpace,measure::Measure)
    mesh = measure.mesh
    dom = measure.domain
    @assert !is_reference_domain(dom)
    d = num_dims(dom)
    @assert num_dims(domain(space)) == d
    # TODO assumes same reference elems for integration and for interpolation
    rid_to_point_to_x = map(coordinates,reference_quadratures(measure))
    rid_to_tab = map(rid_to_point_to_x,reference_fes(space)) do point_to_x, refface
        tabulator(refface)(f,point_to_x)
    end
    face_to_rid = face_reference_id(mesh,d)
    sface_to_face = faces(dom)
    function face_point_dof_s(sface)
        face = sface_to_face[sface]
        rid = face_to_rid[face]
        tab = rid_to_tab[rid]
        function point_dof_s(point,J=nothing)
            function dof_s(dof)
                tab[point,dof]
            end
        end
    end
end

function shape_function_accessor(f::typeof(ForwardDiff.gradient),space::AbstractSpace,measure::Measure)
    mesh = measure.mesh
    dom = measure.domain
    @assert !is_reference_domain(dom)
    d = num_dims(dom)
    @assert num_dims(domain(space)) == d
    rid_to_point_to_x = map(coordinates,reference_quadratures(measure))
    rid_to_tab = map(rid_to_point_to_x,reference_fes(space)) do point_to_x, refface
        tabulator(refface)(f,point_to_x)
    end
    rid_to_dof_to_phys = map(rid_to_tab) do tab
        T = eltype(tab)
        ndofs = size(tab,2)
        # This assumes that the jacobian is a square matrix
        zeros(T,ndofs)
    end
    face_to_rid = face_reference_id(mesh,d)
    sface_to_face = faces(dom)
    function face_point_dof_s(sface)
        face = sface_to_face[sface]
        rid = face_to_rid[face]
        tab = rid_to_tab[rid]
        dof_to_phys = rid_to_dof_to_phys[rid]
        ndofs = length(dof_to_phys)
        function point_dof_s(point,J)
            # NB you cannot evaluate this at more than one point at once
            for dof in 1:ndofs
                dof_to_phys[dof] = transpose(J)\tab[point,dof]
            end
            function dof_s(dof)
                dof_to_phys[dof]
            end
        end
    end
end

function dofs_accessor(space::AbstractSpace,dom::AbstractDomain)
    @assert num_fields(space) == 1
    d = num_dims(dom)
    dom2 = domain(space)
    @assert num_dims(dom2) == d
    sface_to_face = faces(dom)
    face_to_rface = inverse_faces(dom)
    rface_to_dofs = face_dofs(space)
    function face_to_dofs(sface)
        face = sface_to_face[sface]
        rface = face_to_rface[face]
        dofs = rface_to_dofs[rface]
        dofs
    end
end

function discrete_field_accessor(f,uh::DiscreteField,measure::Measure)
    dom = domain(measure)
    space = GT.space(uh)
    d = num_dims(dom)
    dom2 = domain(space)
    @assert num_dims(dom2) == d
    sface_to_face = faces(dom)
    face_to_rface = inverse_faces(dom)
    rface_to_dofs = face_dofs(space)
    rface_to_rid = face_reference_id(space)
    free_vals = free_values(uh)
    diri_vals = dirichlet_values(uh)
    face_point_dof_s = shape_function_accessor(f,space,measure)
    function face_point_val(sface)
        face = sface_to_face[sface]
        rface = face_to_rface[face]
        dofs = rface_to_dofs[rface]
        rid = rface_to_rid[rface]
        point_dof_s = face_point_dof_s(sface)
        ndofs = length(dofs)
        function point_val(point,J)
            dof_s = point_dof_s(point,J)
            sum(1:ndofs) do i
                dof = dofs[i]
                s = dof_s(i)
                if dof > 0
                    v = free_vals[dof]
                else
                    v = diri_vals[-dof]
                end
                v*s
            end
        end
    end
end

function dirichlet_accessor(uh::DiscreteField,dom::AbstractDomain)
    space = GT.space(uh)
    d = num_dims(dom)
    dom2 = domain(space)
    @assert num_dims(dom2) == d
    sface_to_face = faces(dom)
    face_to_rface = inverse_faces(dom)
    rface_to_dofs = face_dofs(space)
    rid_to_u = map(reference_fes(space)) do fe
        T = eltype(dirichlet_values(uh))
        zeros(T,num_dofs(fe))
    end
    rface_to_rid = face_reference_id(space)
    diri_vals = dirichlet_values(uh)
    function face_dirichlet!(sface)
        face = sface_to_face[sface]
        rface = face_to_rface[face]
        dofs = rface_to_dofs[rface]
        rid = rface_to_rid[rface]
        u = rid_to_u[rid]
        fill!(u,zero(eltype(u)))
        for (i,dof) in enumerate(dofs)
            if dof < 0
                u[i] = diri_vals[-dof]
            end
        end
        function dirichlet!(A,b)
            m,n = size(A)
            z = zero(eltype(b))
            for i in 1:m
                bi = z
                for j in 1:n
                    bi += A[i,j]*u[j]
                end
                b[i] -= bi
            end
        end
    end
end

function allocate_vector(
        ::Type{T},space_test::AbstractSpace,domains::AbstractDomain...;
        vector_strategy = monolithic_vector_assembly_strategy(),
    ) where T
    setup = vector_strategy.init(GT.free_dofs(space_test),T)
    counter = vector_strategy.counter(setup)
    @assert num_fields(space_test) == 1
    field_test = 1
    for domain in domains
        face_dofs = dofs_accessor(space_test,domain)
        nfaces = num_faces(domain)
        for face in 1:nfaces
            dofs = face_dofs(face)
            for dof_test in dofs
                counter = vector_strategy.count(counter,setup,dof_test,field_test)
            end
        end
    end
    coo = vector_strategy.allocate(counter,setup)
    counter_ref = Ref(vector_strategy.counter(setup))
    data = (;setup,coo,counter_ref,vector_strategy)
    VectorAllocation(data)
end

struct VectorAllocation{A} <: AbstractType
    data::A
end

function contribute!(alloc::VectorAllocation,b,dofs_test,field_test=1)
    (;setup,coo,counter_ref,vector_strategy) = alloc.data
    counter = counter_ref[]
    for (i,dof_test) in enumerate(dofs_test)
        counter = vector_strategy.set!(coo,counter,setup,b[i],dof_test,field_test)
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

function allocate_matrix(
        ::Type{T},space_test::AbstractSpace,space_trial,domains::AbstractDomain...;
        matrix_strategy = monolithic_matrix_assembly_strategy(),
    ) where T
    setup = matrix_strategy.init(GT.free_dofs(space_test),GT.free_dofs(space_trial),T)
    counter = matrix_strategy.counter(setup)
    @assert num_fields(space_test) == 1
    @assert num_fields(space_trial) == 1
    field_test = 1
    field_trial = 1
    for domain in domains
        face_dofs_test = dofs_accessor(space_test,domain)
        face_dofs_trial = dofs_accessor(space_trial,domain)
        nfaces = num_faces(domain)
        for face in 1:nfaces
            dofs_test = face_dofs_test(face)
            dofs_trial = face_dofs_trial(face)
            for dof_test in dofs_test
                for dof_trial in dofs_trial
                    counter = matrix_strategy.count(counter,setup,dof_test,dof_trial,field_test,field_trial)
                end
            end
        end
    end
    coo = matrix_strategy.allocate(counter,setup)
    counter_ref = Ref(matrix_strategy.counter(setup))
    data = (;setup,coo,counter_ref,matrix_strategy)
    MatrixAllocation(data)
end

struct MatrixAllocation{A} <: AbstractType
    data::A
end

function contribute!(alloc::MatrixAllocation,b,dofs_test,dofs_trial,field_test=1,field_trial=1)
    (;setup,coo,counter_ref,matrix_strategy) = alloc.data
    counter = counter_ref[]
    for (i,dof_test) in enumerate(dofs_test)
        for (j,dof_trial) in enumerate(dofs_trial)
            counter = matrix_strategy.set!(coo,counter,setup,b[i,j],dof_test,dof_trial,field_test,field_trial)
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

# TODO remove after merging branch assembly

function inverse_faces(domain::AbstractDomain)
    d = num_dims(domain)
    ndfaces = num_faces(mesh(domain),d)
    dface_to_face = zeros(Int32,ndfaces)
    face_to_dface = faces(domain)
    dface_to_face[face_to_dface] = 1:length(face_to_dface)
    dface_to_face
end

function inverse_faces(domain::AbstractDomain{<:PMesh})
    map(GT.inverse_faces,partition(domain))
end

function analytical_field_tmp(callee,domain)
    AnalyticalField(mesh(domain),callee,domain)
end

struct AnalyticalField{A,B,C} <: AbstractQuantity{A}
    mesh::A
    definition::B
    domain::C
end

prototype(a::AnalyticalField) = a.definition
domain(a::AnalyticalField) = a.domain
term(a::AnalyticalField) = term(constant_quantity(a.definition,a.domain))







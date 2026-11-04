# TODO: check topologies
# TODO: reuse data for multiple kernels, including boundary/skeleon and integration over the same domain but different element types
# TODO: can we fix some tpyes to reduce compilation time?
# TODO: is the GPU layout really GPU friendly? other libraries like cuda-dolfinx and firedrake uses the original layout and load them into private/shared memory
struct MeshIntegrationData{L,A,B,C,D}
    layout::L # CPULayout or GPULayout
    node_coordinates::A # cpu layout: Vector{StaticArraysCore.SVector{d, Float64}} indexed by node_coordinates[node]
                        # gpu layout: Matrix{Float64}  indexed by node_coordinates[node, d]
    face_nodes::B # cpu layout: Vector{StaticArray{}} indexed by  face_nodes[face][node]
                  # gpu layout: Matrix{}  indexed by face_nodes[face, node]
    reference_gradients::C # cpu layout: Matrix{StaticArraysCore.SVector{d, Float64}}   indexed by reference_gradients[node, point] 
                           # we don't have a gpu layout for this, as it is usually on the constant memory
                           # TODO: if simplex: only store 1 quadrature point
    reference_weights::D # Vector{Float64}   indexed by  reference_weights[point]
end

node_coordinates(a::MeshIntegrationData) = a.node_coordinates
face_nodes(a::MeshIntegrationData) = a.face_nodes
reference_gradients(a::MeshIntegrationData) = a.reference_gradients
reference_weights(a::MeshIntegrationData) = a.reference_weights

# lower the mesh to plain arrays with a fixed data layout
# TODO: support the boundary/skeleton, and reuse the memory.
function mesh_integration_data(mesh::AbstractMesh, Q::AbstractQuadrature)
    # TODO: get a subset, like boundary, dimension (now we assume that it is complexified), etc.
    dom = GT.domain(Q)
    d = num_dims(dom)
    D = num_dims(mesh)
    node_coordinates = GT.node_coordinates(mesh)
    layout = CPULayout() # TODO: build the function on the GPU or with GPU layout directly 
    face_nodes = GT.face_nodes(mesh, d)  

    reference_ids = face_reference_id(mesh)
    rid = reference_ids[d+1][1] # assuming that there is only 1 element type in the mesh or the selected subset.
    @assert all(x -> x == rid, reference_ids[d+1]) 
    # Maybe we need to remove this layer of abstraction of MeshQuadrature and FaceQuadrature
    
    # remove jagged array. use the same array size for all
    
    static_face_nodes = SVector{length(first(face_nodes))}.(face_nodes)


    # gradients tabulation from the reference accessors. 
    # TODO: check correctness. Also check the number of quadrature points / DOFs per simplex element
    ref_quadrature = reference_quadratures(Q)[rid]
    reference_weights = weights(ref_quadrature)
    point_to_x = coordinates(ref_quadrature)
    
    # m2 = GT.mesh(domain(Q))
    # space = mesh_space(m2,Val(D))

    # m2 = GT.mesh(domain(Q))
    # space = mesh_space(m2,Val(2))

    # mesh_domain = GT.mesh_domain(m; num_dims=Val(D))
    # space = mesh_space(m,D)
    # refface = reference_spaces(space)[rid]
    refface = reference_spaces(mesh, Val(D))[rid]
    reference_gradients = collect(permutedims(tabulator(refface)(ForwardDiff.gradient,point_to_x), (2, 1))) 

    MeshIntegrationData(layout, node_coordinates, static_face_nodes, reference_gradients, reference_weights)
end


# check whether it is possible to load the layout and type to the GPU, although these are only used for construction and type inference
struct SpaceIntegrationData{F,A}
    # layout::L # CPULayout or GPULayout
    tabulate_function::F # tabulated functions over the reference element. maybe not needed
    # type::Type{T} # element type of the CPU tabulation, including its static shape
    reference_tabulator::A # layout: Matrix{T}   indexed by reference_gradients[dof, point] 
end


reference_tabulator(a::SpaceIntegrationData) = a.reference_tabulator

function space_integration_data(f, V::AbstractSpace, Q::AbstractQuadrature)
    dom = GT.domain(Q)
    mesh = GT.mesh(dom)
    d = num_dims(dom)

    reference_ids = face_reference_id(mesh)
    rid = reference_ids[d+1][1] # assuming that there is only 1 element type in the mesh or the selected subset.
    @assert all(x -> x == rid, reference_ids[d+1]) 
    

    # tabulation from the reference accessors. 
    ref_quadrature = reference_quadratures(Q)[rid]
    # reference_weights = weights(ref_quadrature)
    point_to_x = coordinates(ref_quadrature)
    refface = reference_spaces(V)[rid]
    reference_tabulator = collect(permutedims(tabulator(refface)(f,point_to_x), (2, 1))) 

    SpaceIntegrationData(f, reference_tabulator)
end


# View packed StaticArrays as scalars without copying. Flatten all static
# dimensions into one component axis, including single-component/empty arrays.
function integration_scalar_view(values::AbstractArray{T}) where {T<:StaticArray}
    reshape(reinterpret(eltype(T), vec(values)), length(T), size(values)...)
end

integration_gpu_layout(values::AbstractArray) = permutedims(values, (2, 1))

function integration_gpu_layout(values::AbstractArray{T}) where {T<:StaticArray}
    scalars = integration_scalar_view(values)
    permutedims(scalars, reverse(ntuple(identity, ndims(scalars))))
end

function integration_cpu_layout(values::AbstractArray, ::Type{T}) where {T}
    permutedims(values, (2, 1))
end

function integration_cpu_layout(values::AbstractArray, ::Type{T}) where {T<:StaticArray}
    # Allocate the final packed buffer and permute directly into its scalar view.
    # Reinterpreting an intermediate permutation would leave a wrapper instead
    # of the plain arrays expected by the CPU layout.
    result = Array{T}(undef, reverse(Base.front(size(values))))
    permutedims!(integration_scalar_view(result), values,
        reverse(ntuple(identity, ndims(values))))
    result
end

# function to_layout(data::SpaceIntegrationData, layout)
#     if isa(layout, GPULayout) && isa(data.layout, CPULayout)
#         reference_tabulator = integration_gpu_layout(data.reference_tabulator)
#     elseif isa(layout, CPULayout) && isa(data.layout, GPULayout)
#         reference_tabulator = integration_cpu_layout(data.reference_tabulator, data.type)
#     else
#         return data
#     end
#     SpaceIntegrationData(layout, data.tabulate_function, data.type, reference_tabulator)
# end

function to_layout(data::MeshIntegrationData, layout)
    if isa(layout, GPULayout) && isa(data.layout, CPULayout)
        node_coordinates = integration_gpu_layout(data.node_coordinates)
        face_nodes = integration_gpu_layout(data.face_nodes)
        # reference_gradients = integration_gpu_layout(data.reference_gradients)
    elseif isa(layout, CPULayout) && isa(data.layout, GPULayout)
        # We assume that the mesh data are always vectors (dimension >= 2)
        coordinates_type = SVector{size(data.node_coordinates, 2), eltype(data.node_coordinates)}
        nodes_type = SVector{size(data.face_nodes, 2), eltype(data.face_nodes)}
        node_coordinates = integration_cpu_layout(data.node_coordinates, coordinates_type)
        face_nodes = integration_cpu_layout(data.face_nodes, nodes_type)
        # gradients_type = SVector{size(data.reference_gradients, 3), eltype(data.reference_gradients)}
        # reference_gradients = integration_cpu_layout(data.reference_gradients, gradients_type)
    else
        return data
    end
    # Weights already have the same contiguous layout on both sides, and all threads may load the same reference data point at a time.
    MeshIntegrationData(layout, node_coordinates, face_nodes, data.reference_gradients, data.reference_weights)
end


# TODO: expressions to be moved into compiler.jl when it is stable
# TODO: test the expressions against existing accessors

# TODO: remove this when we can evaluate code in GT. This is only used for external debugs
function _GT_debug_eval(expr)
    eval(expr)
end

# Di,Dj: 2 jacobian matrix dims; z: zero scalar value of the matrix;
# nlnodes: number of local nodes; lnode_to_node: local node to global node
# node_to_x: global node to coordinate; lnode_s: tabulated gradient



function jacobian_sum_function_expr(Di, Dj, z, nlnodes; is_simplex=false) # TODO: implement for simplex elements
    # copied from another file. maybe we want to change it later to have a faster jacobian computation
    result = Expr(:block)
    symbols = []
    for j in 1:Dj
        for i in 1:Di 
            v = Symbol("j$i$j")
            push!(symbols, v)
            push!(result.args, :($v=$z))
        end
    end
    
    body = []
    push!(body, :(node = lnode_to_node[lnode]))
    push!(body, :(x = node_to_x[node]))
    push!(body, :(s = lnode_s[lnode, point]))
    for i in 1:Di
        for j in 1:Dj
            v = Symbol("j$i$j")
            push!(body, :($v += x[$i] * s[$j]))
        end
    end

    loop = quote
        for lnode in 1:$nlnodes
            $(body...)
        end
    end

    push!(result.args, loop)
    push!(result.args, :(SMatrix{$Di, $Dj}($(symbols...))) )
    result
    
    # TODO: check acceleratedkernels. How does it port customized functions to GPUs?
    :((lnode_to_node, node_to_x, lnode_s, point) -> begin
        $result
    end
    )
end

# Here we only define the code generation for CPUs. TODO: for GPUs, we need to design data collectors from global memory to local/shared mem
function jacobian_expr(mesh, Q, mesh_expr, Q_expr, f_expr, q_expr) 
    # mesh, quadrature rule (FaceQuadrature), 4 symbols for mesh, quadrature, face (the element) id, quadrature point id
    data_expr = :(mesh_integration_data($mesh_expr, $Q_expr)) # TODO: add GPU convertion and hoist out the precomputation to CPU
    lnode_to_node_expr = :(face_nodes($data_expr)[$f_expr])
    node_to_x_expr = :(node_coordinates($data_expr))
    lnode_s_expr = :(reference_gradients($data_expr))

    z = zero(eltype(GT.node_coordinates(mesh)))
    z_element = zero(eltype(z))
    dom = GT.domain(Q)
    d = num_dims(dom)
    face_nodes = GT.face_nodes(mesh, d)
    Di = length(z)
    Dj = d

    nlnodes = length(first(face_nodes))
    jacobian_sum_function = jacobian_sum_function_expr(Di, Dj, z_element, nlnodes; is_simplex=false)
    # we place the generated jacobian function into a lambda expression
    # so that we can hoist the entire function out later for GPUs and remove redundancies if 2 lambdas are literally the same.
    :(($jacobian_sum_function)($lnode_to_node_expr, $node_to_x_expr, $lnode_s_expr, $q_expr))

end


function change_of_measure_expr(J, Di, Dj)
    if Di == Dj
        return :(abs(det($J)))
    elseif Di == 3 && Dj == 2
        return :(norm(cross($J[:,1], $J[:,2]))) # Note: Julia col-major
    else
        return :(sqrt(det(transpose($J)*$J)))
    end
    # potential optimization for 4x3, but we don't implement this.  
    # elseif M == 4 && N == 3
    #     # Cauchy-Binet: the squared volume is the sum of squared 3x3
    #     # minors. Reuse the six 2x2 minors of the first two columns.
    #     a12 = J[1,1]*J[2,2] - J[2,1]*J[1,2]
    #     a13 = J[1,1]*J[3,2] - J[3,1]*J[1,2]
    #     a14 = J[1,1]*J[4,2] - J[4,1]*J[1,2]
    #     a23 = J[2,1]*J[3,2] - J[3,1]*J[2,2]
    #     a24 = J[2,1]*J[4,2] - J[4,1]*J[2,2]
    #     a34 = J[3,1]*J[4,2] - J[4,1]*J[3,2]
    #     minors = SVector(
    #         J[1,3]*a23 - J[2,3]*a13 + J[3,3]*a12,
    #         J[1,3]*a24 - J[2,3]*a14 + J[4,3]*a12,
    #         J[1,3]*a34 - J[3,3]*a14 + J[4,3]*a13,
    #         J[2,3]*a34 - J[3,3]*a24 + J[4,3]*a23)
    #     return norm(minors)
    # end
end

function weight_expr(mesh, Q, mesh_expr, Q_expr, f_expr, q_expr) # mesh, quadrature rule (FaceQuadrature), face (the element) id, quadrature point id
    data_expr = :(mesh_integration_data($mesh_expr, $Q_expr))
    rw_expr = :(reference_weights($data_expr)[$q_expr])

    dom = GT.domain(Q)
    if is_reference_domain(dom)
        return rw_expr
    end

    J_expr = jacobian_expr(mesh, Q, mesh_expr, Q_expr, f_expr, q_expr)
    z = zero(eltype(GT.node_coordinates(mesh)))
    d = num_dims(dom)
    Di = length(z)
    Dj = d

    change_of_measure = change_of_measure_expr(J_expr, Di, Dj)

    :($change_of_measure * $rw_expr)
end


@inline function adjugate(J::StaticMatrix{1,1,T}) where {T}
    SMatrix{1,1,T}(one(T))
end

@inline function adjugate(J::StaticMatrix{2,2,T}) where {T}
    SMatrix{2,2,T}(J[2,2], -J[2,1], -J[1,2], J[1,1])
end

@inline function adjugate(J::StaticMatrix{3,3,T}) where {T}
    c11 = J[2,2]*J[3,3] - J[2,3]*J[3,2]
    c12 = J[2,3]*J[3,1] - J[2,1]*J[3,3]
    c13 = J[2,1]*J[3,2] - J[2,2]*J[3,1]
    c21 = J[1,3]*J[3,2] - J[1,2]*J[3,3]
    c22 = J[1,1]*J[3,3] - J[1,3]*J[3,1]
    c23 = J[1,2]*J[3,1] - J[1,1]*J[3,2]
    c31 = J[1,2]*J[2,3] - J[1,3]*J[2,2]
    c32 = J[1,3]*J[2,1] - J[1,1]*J[2,3]
    c33 = J[1,1]*J[2,2] - J[1,2]*J[2,1]
    SMatrix{3,3,T}(c11, c12, c13, c21, c22, c23, c31, c32, c33)
end

@inline function adjugate(J::StaticMatrix{N,N,T}) where {N,T}
    indices = SVector{N}(1:N)
    entries = ntuple(Val(N*N)) do k
        i = mod1(k, N)
        j = cld(k, N)
        rows = StaticArrays.deleteat(indices, j)
        cols = StaticArrays.deleteat(indices, i)
        cofactor = det(J[rows, cols])
        isodd(i + j) ? -cofactor : cofactor
    end
    SMatrix{N,N,T}(entries)
end

# For real full-column-rank geometry Jacobians, J^+ = (J^T J)^(-1) J^T.
# Here J^+ denotes the Moore-Penrose pseudoinverse and J^T the transpose.
# These explicit formulas assume full rank and do not truncate singular values.
@inline function jacobian_pinv(J::StaticMatrix{2,1,T}) where T
    # Curve: J = a, s = a^T a, and J^+ = a^T * (1/s).
    a1 = J[1,1]
    a2 = J[2,1]
    r = inv(a1*a1 + a2*a2)
    SMatrix{1,2}(a1*r, a2*r)
end

@inline function jacobian_pinv(J::StaticMatrix{3,1,T}) where T
    # Curve: J = a, s = a^T a, and J^+ = a^T * (1/s).
    a1 = J[1,1]
    a2 = J[2,1]
    a3 = J[3,1]
    r = inv(a1*a1 + a2*a2 + a3*a3)
    SMatrix{1,3}(a1*r, a2*r, a3*r)
end

@inline function jacobian_pinv(J::StaticMatrix{3,2,T}) where T
    # Surface: J = [a b], n = a x b, s = n^T n = det(J^T J), q = n * (1/s).
    # J^+ has rows (b x q)^T and (q x a)^T; x denotes the cross product.
    # Scaling the shared normal first costs 24 multiplies, 11 adds/subtracts,
    # and one division in total. The integration measure is sqrt(s).
    a1 = J[1,1]
    a2 = J[2,1]
    a3 = J[3,1]
    b1 = J[1,2]
    b2 = J[2,2]
    b3 = J[3,2]
    n1 = a2*b3 - a3*b2
    n2 = a3*b1 - a1*b3
    n3 = a1*b2 - a2*b1
    r = inv(n1*n1 + n2*n2 + n3*n3)
    q1 = n1*r
    q2 = n2*r
    q3 = n3*r
    u1 = b2*q3 - b3*q2
    u2 = b3*q1 - b1*q3
    u3 = b1*q2 - b2*q1
    v1 = q2*a3 - q3*a2
    v2 = q3*a1 - q1*a3
    v3 = q1*a2 - q2*a1
    SMatrix{2,3}(u1, v1, u2, v2, u3, v3)
end

@inline jacobian_pinv(J::StaticMatrix) = pinv(J)

# TODO: unroll modif, determinant, transpose, adjugate ... later
function modif_expr(func::typeof(ForwardDiff.gradient), v_expr, invJ_expr)
    :(transpose($invJ_expr) * $v_expr)
end

function modif_expr(func::typeof(ForwardDiff.jacobian), v_expr, invJ_expr)
    :($v_expr*$invJ_expr)
end

function modif_expr(func::typeof(value), v_expr, invJ_expr)
    v_expr
end

# tabulate function (can be gradient, value, etc.), space, quadrature rule, face id, quadrature point id, dof id
function shape_function_expr(func, mesh, V, Q, mesh_expr, V_expr, Q_expr, f_expr, q_expr, i_expr) 
    # mesh_data_expr = :(mesh_integration_data($mesh_expr, $Q_expr))
    space_data_expr = :(space_integration_data($func, $V_expr, $Q_expr))
    v_expr = :(reference_tabulator($space_data_expr)[$i_expr, $q_expr])

    if is_reference_domain(GT.domain(V)) || func === value
        return v_expr
    end

    J_expr = jacobian_expr(mesh, Q, mesh_expr, Q_expr, f_expr, q_expr)
    z = zero(eltype(GT.node_coordinates(mesh)))
    d = num_dims(GT.domain(V))
    Di = length(z)
    Dj = d

    # TODO: implement also for untabulated functions (without space_data_expr) 
    # TODO: (can we process external functions on GPUs?)

    invJ_expr = if Di == Dj # inv(J) = adj(J)/det(J) for square matrices
        :(adjugate($J_expr) * inv(det($J_expr)))
    else
        :(jacobian_pinv($J_expr))
    end

    modif_expr(func, v_expr, invJ_expr)

end

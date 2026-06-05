module GPUNewTests

using Test
using LinearAlgebra
using SparseArrays
using Adapt
using BenchmarkTools
import Atomix
import PartitionedArrays as PA
import GalerkinToolkit as GT
import KernelAbstractions as KA
using KernelAbstractions: @kernel, @index, @localmem
using StaticArrays

const HAS_CUDA = try
    @eval using CUDA
    true
catch
    false
end

const HAS_AMD = try
    @eval using AMDGPU
    true
catch
    false
end

function is_cuda_available()
    return HAS_CUDA && !isempty(CUDA.devices())
end

function is_rocm_available()
    return HAS_AMD && !isempty(AMDGPU.devices())
end

function select_backend()
    if is_cuda_available()
        dev = first(CUDA.devices())
        println("using CUDA backend: $(CUDA.name(dev))")
        println()
        return CUDABackend()
    elseif is_rocm_available()
        dev = first(AMDGPU.devices())
        println("using AMD backend: $(AMDGPU.HIP.name(dev))")
        println()
        return ROCBackend()
    else
        println("using CPU backend")
        println()
        return CPU()
    end
end

if is_cuda_available()
    macro alloc_shared_sta(T, dims)
        ex = :( CUDA.CuStaticSharedArray($T, $dims) )
        return esc(ex)
    end
    macro alloc_shared_dyn(T, dims)
        ex = :( CUDA.CuDynamicSharedArray($T, $dims) )
        return esc(ex)
    end
    macro call_kernel(name, threads, blocks, args...)
        func_name = name isa Symbol ? Symbol(name, :!) : name
        ex = :( @cuda threads=$threads blocks=$blocks $func_name($(args...)) )
        return esc(ex)
    end
elseif is_rocm_available()
    macro alloc_shared_sta(T, dims)
        ex = :( AMDGPU.@ROCStaticLocalArray($T, $dims) )
        return esc(ex)
    end
    macro alloc_shared_dyn(T, dims)
        ex = :( AMDGPU.@ROCDynamicLocalArray($T, $dims) )
        return esc(ex)
    end
    macro call_kernel(name, threads, blocks, args...)
        func_name = name isa Symbol ? Symbol(name, :!) : name
        ex = :( AMDGPU.@roc groupsize=$threads gridsize=$blocks $func_name($(args...)) )
        return esc(ex)
    end
    macro call_kernel_shmem(name, threads, blocks, shmem, args...)
        func_name = name isa Symbol ? Symbol(name, :!) : name
        ex = :( AMDGPU.@roc groupsize=$threads gridsize=$blocks shmem=$shmem $func_name($(args...)) )
        return esc(ex)
    end
else
    macro alloc_shared_sta(args...)
        error("No GPU backend available to compile @alloc_shared_sta")
    end
    macro alloc_shared_dyn(args...)
        error("No GPU backend available to compile @alloc_shared_dyn")
    end
    macro call_kernel(args...)
        error("No GPU backend available to compile @call_kernel")
    end
    macro call_kernel_shmem(args...)
        error("No GPU backend available to compile @call_kernel_shmem")
    end
end

const dev = select_backend()

function max_shared_memory_per_block()
    if is_cuda_available()
        return Int(CUDA.attribute(first(CUDA.devices()), CUDA.DEVICE_ATTRIBUTE_MAX_SHARED_MEMORY_PER_BLOCK))
    elseif is_rocm_available()
        return Int(AMDGPU.HIP.properties(first(AMDGPU.devices())).sharedMemPerBlock)
    else
        return typemax(Int)
    end
end

f(x) = 2*sin(sum(x))

# 0-form for analytical solution
function cpu_loop_1(dΩ_faces)
    s = 0.0
    for dΩ_face in dΩ_faces
        for dΩ_point in GT.each_point_new(dΩ_face)
            x = GT.coordinate(dΩ_point)
            dx = GT.weight(dΩ_point)
            s += f(x)*dx
        end
    end
    return s
end

# 0-form for discrete field
function cpu_loop_2(uh_faces)
    s = 0.0
    for uh_face in uh_faces
        for uh_point in GT.each_point_new(uh_face)
            ux = GT.field(GT.value,uh_point)
            dx = GT.weight(uh_point)
            s += ux*dx
        end
    end
    return s
end

# 0-form for discrete field but with a syntax that potentially
# reuses the jacobians.
function cpu_loop_3(uh_faces,dΩ_faces)
    s = 0.0
    for dΩ_face in dΩ_faces
        uh_face = uh_faces[dΩ_face]
        for dΩ_point in GT.each_point_new(dΩ_face)
            uh_point = uh_face[dΩ_point]
            ux = GT.field(GT.value,uh_point)
            dx = GT.weight(dΩ_point)
            s += ux*dx
        end
    end
    return s
end

# 1-form equivalent to matrix free
function cpu_loop_4!(b,bf,uh_faces)
    for uh_face in uh_faces
        dofs = GT.dofs(uh_face)
        n = GT.num_dofs(uh_face)
        fill!(bf,0)
        for uh_point in GT.each_point_new(uh_face)
            ux = GT.field(GT.gradient,uh_point)
            sx = GT.shape_functions(GT.gradient,uh_point)
            dx = GT.weight(uh_point)
            ux_dx = ux*dx
            for i in 1:n
                bf[i] += ux_dx⋅sx[i]
            end
        end
        for i in 1:n
            b[dofs[i]] += bf[i]
        end
    end
end

# 1-form equivalent to matrix free
# where the discrete field
# and trial functions are potentially different
function cpu_loop_5!(b,bf,uh_faces,V_faces,dΩ_faces)
    for dΩ_face in dΩ_faces
        V_face = V_faces[dΩ_face]
        uh_face = uh_faces[dΩ_face]
        dofs = GT.dofs(V_face)
        n = GT.num_dofs(V_face)
        fill!(bf,0)
        for dΩ_point in GT.each_point_new(dΩ_face)
            V_point = V_face[dΩ_point]
            uh_point = uh_face[dΩ_point]
            ux = GT.field(GT.gradient,uh_point)
            sx = GT.shape_functions(GT.gradient,V_point)
            dx = GT.weight(dΩ_point)
            ux_dx = ux*dx
            for i in 1:n
                bf[i] += ux_dx⋅sx[i]
            end
        end
        for i in 1:n
            b[dofs[i]] += bf[i]
        end
    end
end

# 2-form assembly in coo format
function cpu_loop_6_count(V_faces)
    num_nz = 0
    for V_face in V_faces
        n = GT.num_dofs(V_face)
        num_nz += n*n
    end
    num_nz
end

function cpu_loop_6_symbolic!(AI,AJ,V_faces)
    num_nz = 0
    for V_face in V_faces
        n = GT.num_dofs(V_face)
        dofs = GT.dofs(V_face)
        for j in 1:n
            gj = dofs[j]
            for i in 1:n
                gi = dofs[i]
                num_nz += 1
                AI[num_nz] = gi
                AJ[num_nz] = gj
            end
        end
    end
end

function cpu_loop_6_numeric!(AV,Af,V_faces)
    num_nz = 0
    for V_face in V_faces
        n = GT.num_dofs(V_face)
        fill!(Af,0)
        for V_point in GT.each_point_new(V_face)
            dx = GT.weight(V_point)
            sx = GT.shape_functions(GT.gradient,V_point)
            for j in 1:n
                sx_dx_j = sx[j]*dx
                for i in 1:n
                    Af[i,j] += sx[i]⋅sx_dx_j
                end
            end
        end
        for j in 1:n
            for i in 1:n
                num_nz += 1
                AV[num_nz] = Af[i,j]
            end
        end
    end
end

function cpu_loop_6_ltable!(ltable,V_faces)
    num_nz = 0
    face_id = 1
    for V_face in V_faces
        ltable[face_id] = num_nz
        n = GT.num_dofs(V_face)
        num_nz += n*n
        face_id += 1
    end
end

function cpu_loop_6_numeric_ltable!(AV,V_faces,ltable)
    face_id = 1
    for V_face in V_faces
        n = GT.num_dofs(V_face)
        offset = ltable[face_id]
        for V_point in GT.each_point_new(V_face)
            dx = GT.weight(V_point)
            sx = GT.shape_functions(GT.gradient,V_point)
            for j in 1:n
                sx_dx_j = sx[j]*dx
                for i in 1:n
                    AV[offset + (j-1)*n + i] += sx[i]⋅sx_dx_j
                end
            end
        end
        face_id += 1
    end
end

function main_cpu(params)
    (;face_nodes_layout,face_dofs_layout,n) = params

    # Start at CPU
    domain = (0,1,0,1)
    cells = (n,n)
    mesh = GT.cartesian_mesh(domain,cells)
    Ω = GT.interior(mesh)
    degree = 4
    dΩ = GT.quadrature(Ω,degree)

    k = 1
    V = GT.lagrange_space(Ω,k)
    u = GT.analytical_field(f,Ω)
    uh = GT.interpolate(u,V)

    println("num_dims = ", GT.num_dims(V))
    println("num_faces = ", GT.num_faces(mesh))
    println("num_free_dofs = ", GT.num_free_dofs(V))
    println("num_dirichlet_dofs = ", GT.num_dirichlet_dofs(V))
    println("num_dofs = ", GT.num_dofs(V))
    println("max_num_reference_dofs = ", GT.max_num_reference_dofs(V))

    tabulate = (GT.value,GT.gradient)
    dΩ_faces_cpu = GT.each_face_new(dΩ)
    V_faces_cpu = GT.each_face_new(V,dΩ;tabulate)
    uh_faces_cpu = GT.each_face_new(uh,dΩ;tabulate)

    # Change data layout
    dΩ_faces_cpu = GT.change_data_layout(dΩ_faces_cpu;face_nodes_layout)
    V_faces_cpu = GT.change_data_layout(V_faces_cpu;face_nodes_layout,face_dofs_layout)
    uh_faces_cpu = GT.change_data_layout(uh_faces_cpu;face_nodes_layout,face_dofs_layout)

    # This is not needed in practice, just to make sure that
    # we do not break anything when adapting.
    dΩ_faces_cpu = Adapt.adapt_structure(Array,dΩ_faces_cpu)
    V_faces_cpu = Adapt.adapt_structure(Array,V_faces_cpu)
    uh_faces_cpu = Adapt.adapt_structure(Array,uh_faces_cpu)

    @show r_cpu = cpu_loop_1(dΩ_faces_cpu)
    @show r_cpu = cpu_loop_2(uh_faces_cpu)
    @show r_cpu = cpu_loop_3(uh_faces_cpu,dΩ_faces_cpu)

    b = zeros(GT.num_free_dofs(V))
    nmax = GT.max_num_reference_dofs(V)
    bf = zeros(nmax)
    r_cpu = cpu_loop_4!(b,bf,uh_faces_cpu)
    @show norm(b)

    fill!(b,0)
    r_cpu = cpu_loop_5!(b,bf,uh_faces_cpu,V_faces_cpu,dΩ_faces_cpu)
    @show norm(b)

    num_nz = cpu_loop_6_count(V_faces_cpu)
    AI = zeros(Int32,num_nz)
    AJ = zeros(Int32,num_nz)
    AV = zeros(num_nz)
    Af = zeros(nmax,nmax)
    cpu_loop_6_symbolic!(AI,AJ,V_faces_cpu)
    cpu_loop_6_numeric!(AV,Af,V_faces_cpu)
    n_global = GT.num_dofs(V)
    A,Acache = PA.sparse_matrix(AI,AJ,AV,n_global,n_global;reuse=Val(true))
    x = GT.free_values(uh)
    b = A*x
    @show norm(b)

    # Loop using a previously built matrix
    fill!(AV,0)
    cpu_loop_6_numeric!(AV,Af,V_faces_cpu)
    PA.sparse_matrix!(A,AV,Acache)
    b = A*x
    @show norm(b)

    # Loop using a a lookup table
    fill!(AV,0)
    ltable = zeros(Int32,length(V_faces_cpu))
    cpu_loop_6_ltable!(ltable,V_faces_cpu)
    cpu_loop_6_numeric_ltable!(AV,V_faces_cpu,ltable)
    A,Acache = PA.sparse_matrix(AI,AJ,AV,n_global,n_global;reuse=Val(true))
    PA.sparse_matrix!(A,AV,Acache)
    b = A*x
    @show norm(b)
end

if is_cuda_available()
    function cuda_loop_1!(contributions,dΩ_faces)
        face_id = (blockIdx().x - 1) * blockDim().x + threadIdx().x
        if face_id > length(dΩ_faces)
            return nothing
        end
        dΩ_face = dΩ_faces[face_id]
        s = 0.0
        for dΩ_point in GT.each_point_new(dΩ_face)
            x = GT.coordinate(dΩ_point)
            dx = GT.weight(dΩ_point)
            s += f(x)*dx
        end
        contributions[face_id] = s
        return nothing
    end

    function cuda_loop_2!(contributions,uh_faces)
        face_id = (blockIdx().x - 1) * blockDim().x + threadIdx().x
        if face_id > length(uh_faces)
            return nothing
        end
        uh_face = uh_faces[face_id]
        s = 0.0
        for uh_point in GT.each_point_new(uh_face)
            ux = GT.field(GT.value,uh_point)
            dx = GT.weight(uh_point)
            s += ux*dx
        end
        contributions[face_id] = s
        return nothing
    end

    function cuda_loop_3!(contributions,uh_faces,dΩ_faces)
        face_id = (blockIdx().x - 1) * blockDim().x + threadIdx().x
        if face_id > length(uh_faces)
            return nothing
        end
        s = 0.0
        dΩ_face = dΩ_faces[face_id]
        uh_face = uh_faces[dΩ_face]
        for dΩ_point in GT.each_point_new(dΩ_face)
            uh_point = uh_face[dΩ_point]
            ux = GT.field(GT.value,uh_point)
            dx = GT.weight(dΩ_point)
            s += ux*dx
        end
        contributions[face_id] = s
        return nothing
    end

    function cuda_loop_4_atomic!(b,uh_faces)
        face_id = (blockIdx().x - 1) * blockDim().x + threadIdx().x
        if face_id > length(uh_faces)
            return nothing
        end
        uh_face = uh_faces[face_id]
        dofs = GT.dofs(uh_face)
        n = GT.num_dofs(uh_face)
        for uh_point in GT.each_point_new(uh_face)
            ux = GT.field(GT.gradient,uh_point)
            sx = GT.shape_functions(GT.gradient,uh_point)
            dx = GT.weight(uh_point)
            ux_dx = ux*dx
            for i in 1:n
                Atomix.@atomic b[dofs[i]] += ux_dx⋅sx[i]
            end
        end
        return nothing
    end

    function cuda_loop_4_global!(b,bf_global,uh_faces)
        face_id = (blockIdx().x - 1) * blockDim().x + threadIdx().x
        if face_id > length(uh_faces)
            return nothing
        end
        uh_face = uh_faces[face_id] # Eventually calls at_mesh_face
        dofs = GT.dofs(uh_face)
        n = GT.num_dofs(uh_face)
        bf = view(bf_global, :, face_id)
        fill!(bf,0)
        for uh_point in GT.each_point_new(uh_face)
            ux = GT.field(GT.gradient,uh_point)
            sx = GT.shape_functions(GT.gradient,uh_point)
            dx = GT.weight(uh_point)
            ux_dx = ux*dx
            for i in 1:n
                bf[i] += ux_dx⋅sx[i]
            end
        end
        for i in 1:n
            Atomix.@atomic b[dofs[i]] += bf[i]
        end
        return nothing
    end

    function cuda_loop_4_local!(b,::Val{max_dofs},uh_faces) where {max_dofs}
        face_id = (blockIdx().x - 1) * blockDim().x + threadIdx().x
        if face_id > length(uh_faces)
            return nothing
        end
        uh_face = uh_faces[face_id]
        dofs = GT.dofs(uh_face)
        n = GT.num_dofs(uh_face)
        bf = zeros(SVector{max_dofs, Float64})
        for uh_point in GT.each_point_new(uh_face)
            ux = GT.field(GT.gradient,uh_point)
            sx = GT.shape_functions(GT.gradient,uh_point)
            dx = GT.weight(uh_point)
            ux_dx = ux*dx
            bf = map(enumerate_static(bf)) do (i, bfi)
                bfi + (i <= n ? ux_dx⋅sx[i] : 0)
            end
        end
        for i in 1:n
            Atomix.@atomic b[dofs[i]] += bf[i]
        end
    end

    function cuda_loop_4_shared!(b,::Val{max_dofs},::Val{block_dim},uh_faces) where {max_dofs,block_dim}
        bf_shared = @alloc_shared_sta Float64 (max_dofs,block_dim)
        face_id = (blockIdx().x - 1) * blockDim().x + threadIdx().x
        if face_id > length(uh_faces)
            return nothing
        end
        uh_face = uh_faces[face_id]
        dofs = GT.dofs(uh_face)
        n = GT.num_dofs(uh_face)
        bf = view(bf_shared,:,threadIdx().x)
        fill!(bf, 0)
        for uh_point in GT.each_point_new(uh_face)
            ux = GT.field(GT.gradient,uh_point)
            sx = GT.shape_functions(GT.gradient,uh_point)
            dx = GT.weight(uh_point)
            ux_dx = ux*dx
            for i in 1:n
                bf[i] += ux_dx⋅sx[i]
            end
        end
        for i in 1:n
            Atomix.@atomic b[dofs[i]] += bf[i]
        end
    end

    function cuda_loop_5_atomic!(b,bf_global,uh_faces,V_faces,dΩ_faces)
        face_id = (blockIdx().x - 1) * blockDim().x + threadIdx().x
        if face_id > length(dΩ_faces)
            return nothing
        end
        dΩ_face = dΩ_faces[face_id]
        V_face = V_faces[dΩ_face]
        uh_face = uh_faces[dΩ_face]
        dofs = GT.dofs(V_face)
        n = GT.num_dofs(V_face)
        bf = view(bf_global, :, face_id)
        fill!(bf,0)
        for dΩ_point in GT.each_point_new(dΩ_face)
            V_point = V_face[dΩ_point]
            uh_point = uh_face[dΩ_point]
            ux = GT.field(GT.gradient,uh_point)
            sx = GT.shape_functions(GT.gradient,V_point)
            dx = GT.weight(dΩ_point)
            ux_dx = ux*dx
            for i in 1:n
                bf[i] += ux_dx⋅sx[i]
            end
        end
        for i in 1:n
            Atomix.@atomic b[dofs[i]] += bf[i]
        end
    end

    function cuda_loop_6_numeric_ltable!(AV,V_faces,ltable)
        face_id = (blockIdx().x - 1) * blockDim().x + threadIdx().x
        if face_id > length(V_faces)
            return nothing
        end
        V_face = V_faces[face_id]
        n = GT.num_dofs(V_face)
        offset = ltable[face_id]
        for V_point in GT.each_point_new(V_face)
            dx = GT.weight(V_point)
            sx = GT.shape_functions(GT.gradient,V_point)
            for j in 1:n
                sx_dx_j = sx[j]*dx
                for i in 1:n
                    AV[offset + (j-1)*n + i] += sx[i]⋅sx_dx_j
                end
            end
        end
    end

    function cuda_loop_6_numeric_ltable_local!(AV,V_faces,ltable,::Val{max_dofs}) where {max_dofs}
        face_id = (blockIdx().x - 1) * blockDim().x + threadIdx().x
        if face_id > length(V_faces)
            return nothing
        end
        A_local = zeros(MMatrix{max_dofs, max_dofs, Float64})
        V_face = V_faces[face_id]
        n = GT.num_dofs(V_face)
        offset = ltable[face_id]
        for V_point in GT.each_point_new(V_face)
            dx = GT.weight(V_point)
            sx = GT.shape_functions(GT.gradient,V_point)
            for j in 1:n
                sx_dx_j = sx[j]*dx
                for i in 1:n
                    A_local[i,j] += sx[i]⋅sx_dx_j
                end
            end
        end
        for j in 1:n
            for i in 1:n
                AV[offset + (j-1)*n + i] += A_local[i,j]
            end
        end
    end
elseif is_rocm_available()
    function hip_loop_1!(contributions,dΩ_faces)
        face_id = (workgroupIdx().x - 1) * workgroupDim().x + workitemIdx().x
        if face_id > length(dΩ_faces)
            return nothing
        end
        dΩ_face = dΩ_faces[face_id]
        s = 0.0
        for dΩ_point in GT.each_point_new(dΩ_face)
            x = GT.coordinate(dΩ_point)
            dx = GT.weight(dΩ_point)
            s += f(x)*dx
        end
        contributions[face_id] = s
        return nothing
    end

    function hip_loop_2!(contributions,uh_faces)
        face_id = (workgroupIdx().x - 1) * workgroupDim().x + workitemIdx().x
        if face_id > length(uh_faces)
            return nothing
        end
        uh_face = uh_faces[face_id]
        s = 0.0
        for uh_point in GT.each_point_new(uh_face)
            ux = GT.field(GT.value,uh_point)
            dx = GT.weight(uh_point)
            s += ux*dx
        end
        contributions[face_id] = s
        return nothing
    end

    function hip_loop_3!(contributions,uh_faces,dΩ_faces)
        face_id = (workgroupIdx().x - 1) * workgroupDim().x + workitemIdx().x
        if face_id > length(uh_faces)
            return nothing
        end
        s = 0.0
        dΩ_face = dΩ_faces[face_id]
        uh_face = uh_faces[dΩ_face]
        for dΩ_point in GT.each_point_new(dΩ_face)
            uh_point = uh_face[dΩ_point]
            ux = GT.field(GT.value,uh_point)
            dx = GT.weight(dΩ_point)
            s += ux*dx
        end
        contributions[face_id] = s
        return nothing
    end

    function hip_loop_4_atomic!(b,uh_faces)
        face_id = (workgroupIdx().x - 1) * workgroupDim().x + workitemIdx().x
        if face_id > length(uh_faces)
            return nothing
        end
        uh_face = uh_faces[face_id]
        dofs = GT.dofs(uh_face)
        n = GT.num_dofs(uh_face)
        for uh_point in GT.each_point_new(uh_face)
            ux = GT.field(GT.gradient,uh_point)
            sx = GT.shape_functions(GT.gradient,uh_point)
            dx = GT.weight(uh_point)
            ux_dx = ux*dx
            for i in 1:n
                Atomix.@atomic b[dofs[i]] += ux_dx⋅sx[i]
            end
        end
        return nothing
    end

    function hip_loop_4_global!(b,bf_global,uh_faces)
        face_id = (workgroupIdx().x - 1) * workgroupDim().x + workitemIdx().x
        if face_id > length(uh_faces)
            return nothing
        end
        uh_face = uh_faces[face_id]
        dofs = GT.dofs(uh_face)
        n = GT.num_dofs(uh_face)
        bf = view(bf_global, :, face_id)
        fill!(bf,0)
        for uh_point in GT.each_point_new(uh_face)
            ux = GT.field(GT.gradient,uh_point)
            sx = GT.shape_functions(GT.gradient,uh_point)
            dx = GT.weight(uh_point)
            ux_dx = ux*dx
            for i in 1:n
                bf[i] += ux_dx⋅sx[i]
            end
        end
        for i in 1:n
            Atomix.@atomic b[dofs[i]] += bf[i]
        end
        return nothing
    end

    function hip_loop_4_local!(b,::Val{max_dofs},uh_faces) where {max_dofs}
        face_id = (workgroupIdx().x - 1) * workgroupDim().x + workitemIdx().x
        if face_id > length(uh_faces)
            return nothing
        end
        uh_face = uh_faces[face_id]
        dofs = GT.dofs(uh_face)
        n = GT.num_dofs(uh_face)
        bf = zeros(SVector{max_dofs, Float64})
        for uh_point in GT.each_point_new(uh_face)
            ux = GT.field(GT.gradient,uh_point)
            sx = GT.shape_functions(GT.gradient,uh_point)
            dx = GT.weight(uh_point)
            ux_dx = ux*dx
            bf = map(enumerate_static(bf)) do (i, bfi)
                bfi + (i <= n ? ux_dx⋅sx[i] : 0)
            end
        end
        for i in 1:n
            Atomix.@atomic b[dofs[i]] += bf[i]
        end
    end

    function hip_loop_4_shared!(b,::Val{max_dofs},::Val{block_dim},uh_faces) where {max_dofs,block_dim}
        face_id = (workgroupIdx().x - 1) * workgroupDim().x + workitemIdx().x
        if face_id > length(uh_faces)
            return nothing
        end
        bf_shared = @alloc_shared_dyn Float64 (max_dofs,block_dim)
        uh_face = uh_faces[face_id]
        dofs = GT.dofs(uh_face)
        n = GT.num_dofs(uh_face)
        bf = view(bf_shared,:,workitemIdx().x)
        for uh_point in GT.each_point_new(uh_face)
            ux = GT.field(GT.gradient,uh_point)
            sx = GT.shape_functions(GT.gradient,uh_point)
            dx = GT.weight(uh_point)
            ux_dx = ux*dx
            for i in 1:n
                bf[i] += ux_dx⋅sx[i]
            end
        end
        for i in 1:n
            Atomix.@atomic b[dofs[i]] += bf[i]
        end
    end

    function hip_loop_5_atomic!(b,bf_global,uh_faces,V_faces,dΩ_faces)
        face_id = (workgroupIdx().x - 1) * workgroupDim().x + workitemIdx().x
        if face_id > length(dΩ_faces)
            return nothing
        end
        dΩ_face = dΩ_faces[face_id]
        V_face = V_faces[dΩ_face]
        uh_face = uh_faces[dΩ_face]
        dofs = GT.dofs(V_face)
        n = GT.num_dofs(V_face)
        bf = view(bf_global, :, face_id)
        fill!(bf,0)
        for dΩ_point in GT.each_point_new(dΩ_face)
            V_point = V_face[dΩ_point]
            uh_point = uh_face[dΩ_point]
            ux = GT.field(GT.gradient,uh_point)
            sx = GT.shape_functions(GT.gradient,V_point)
            dx = GT.weight(dΩ_point)
            ux_dx = ux*dx
            for i in 1:n
                bf[i] += ux_dx⋅sx[i]
            end
        end
        for i in 1:n
            Atomix.@atomic b[dofs[i]] += bf[i]
        end
    end

    function hip_loop_6_numeric_ltable!(AV,V_faces,ltable)
        face_id = (workgroupIdx().x - 1) * workgroupDim().x + workitemIdx().x
        if face_id > length(V_faces)
            return nothing
        end
        V_face = V_faces[face_id]
        n = GT.num_dofs(V_face)
        offset = ltable[face_id]
        for V_point in GT.each_point_new(V_face)
            dx = GT.weight(V_point)
            sx = GT.shape_functions(GT.gradient,V_point)
            for j in 1:n
                sx_dx_j = sx[j]*dx
                for i in 1:n
                    AV[offset + (j-1)*n + i] += sx[i]⋅sx_dx_j
                end
            end
        end
    end

    function hip_loop_6_numeric_ltable_local!(AV,V_faces,ltable,::Val{max_dofs}) where {max_dofs}
        face_id = (workgroupIdx().x - 1) * workgroupDim().x + workitemIdx().x
        if face_id > length(V_faces)
            return nothing
        end
        A_local = zeros(MMatrix{max_dofs, max_dofs, Float64})
        V_face = V_faces[face_id]
        n = GT.num_dofs(V_face)
        offset = ltable[face_id]
        for V_point in GT.each_point_new(V_face)
            dx = GT.weight(V_point)
            sx = GT.shape_functions(GT.gradient,V_point)
            for j in 1:n
                sx_dx_j = sx[j]*dx
                for i in 1:n
                    A_local[i,j] += sx[i]⋅sx_dx_j
                end
            end
        end
        for j in 1:n
            for i in 1:n
                AV[offset + (j-1)*n + i] += A_local[i,j]
            end
        end
    end
end

@kernel function gpu_loop_1!(contributions,dΩ_faces)
    face_id = @index(Global)
    if face_id <= length(dΩ_faces)
        dΩ_face = dΩ_faces[face_id]
        s = 0.0
        for dΩ_point in GT.each_point_new(dΩ_face)
            x = GT.coordinate(dΩ_point)
            dx = GT.weight(dΩ_point)
            s += f(x)*dx
        end
        contributions[face_id] = s
    end
end

@kernel function gpu_loop_2!(contributions,uh_faces)
    face_id = @index(Global)
    if face_id <= length(uh_faces)
        uh_face = uh_faces[face_id]
        s = 0.0
        for uh_point in GT.each_point_new(uh_face)
            ux = GT.field(GT.value,uh_point)
            dx = GT.weight(uh_point)
            s += ux*dx
        end
        contributions[face_id] = s
    end
end

@kernel function gpu_loop_3!(contributions,uh_faces,dΩ_faces)
    face_id = @index(Global)
    if face_id <= length(uh_faces)
        s = 0.0
        dΩ_face = dΩ_faces[face_id]
        uh_face = uh_faces[dΩ_face]
        for dΩ_point in GT.each_point_new(dΩ_face)
            uh_point = uh_face[dΩ_point]
            ux = GT.field(GT.value,uh_point)
            dx = GT.weight(dΩ_point)
            s += ux*dx
        end
        contributions[face_id] = s
    end
end

@kernel function gpu_loop_4_atomic!(b,uh_faces)
    face_id = @index(Global)
    if face_id <= length(uh_faces)
        uh_face = uh_faces[face_id]
        dofs = GT.dofs(uh_face)
        n = GT.num_dofs(uh_face)
        for uh_point in GT.each_point_new(uh_face)
            ux = GT.field(GT.gradient,uh_point)
            sx = GT.shape_functions(GT.gradient,uh_point)
            dx = GT.weight(uh_point)
            ux_dx = ux*dx
            for i in 1:n
                Atomix.@atomic b[dofs[i]] += ux_dx⋅sx[i]
            end
        end
    end
end

@kernel function gpu_loop_4_global!(b,bf_global,uh_faces)
    face_id = @index(Global)
    if face_id <= length(uh_faces)
        uh_face = uh_faces[face_id]
        dofs = GT.dofs(uh_face)
        n = GT.num_dofs(uh_face)
        bf = view(bf_global, :, face_id)
        fill!(bf,0)
        for uh_point in GT.each_point_new(uh_face)
            ux = GT.field(GT.gradient,uh_point)
            sx = GT.shape_functions(GT.gradient,uh_point)
            dx = GT.weight(uh_point)
            ux_dx = ux*dx
            for i in 1:n
                bf[i] += ux_dx⋅sx[i]
            end
        end
        for i in 1:n
            Atomix.@atomic b[dofs[i]] += bf[i]
        end
    end
end

@kernel function gpu_loop_4_local!(b,::Val{max_dofs},uh_faces) where {max_dofs}
    face_id = @index(Global)
    if face_id <= length(uh_faces)
        uh_face = uh_faces[face_id]
        dofs = GT.dofs(uh_face)
        n = GT.num_dofs(uh_face)
        bf = zeros(SVector{max_dofs, Float64})
        for uh_point in GT.each_point_new(uh_face)
            ux = GT.field(GT.gradient,uh_point)
            sx = GT.shape_functions(GT.gradient,uh_point)
            dx = GT.weight(uh_point)
            ux_dx = ux*dx
            bf = map(enumerate_static(bf)) do (i, bfi)
                bfi + (i <= n ? ux_dx⋅sx[i] : 0)
            end
        end
        for i in 1:n
            Atomix.@atomic b[dofs[i]] += bf[i]
        end
    end
end

@kernel function gpu_loop_4_shared!(b,::Val{max_dofs},::Val{block_dim},uh_faces) where {max_dofs,block_dim}
    bf_shared = @localmem Float64 (max_dofs,block_dim)
    face_id = @index(Global)
    if face_id <= length(uh_faces)
        uh_face = uh_faces[face_id]
        dofs = GT.dofs(uh_face)
        n = GT.num_dofs(uh_face)
        bf = view(bf_shared,:,@index(Local))
        fill!(bf, 0)
        for uh_point in GT.each_point_new(uh_face)
            ux = GT.field(GT.gradient,uh_point)
            sx = GT.shape_functions(GT.gradient,uh_point)
            dx = GT.weight(uh_point)
            ux_dx = ux*dx
            for i in 1:n
                bf[i] += ux_dx⋅sx[i]
            end
        end
        for i in 1:n
            Atomix.@atomic b[dofs[i]] += bf[i]
        end
    end
end

@kernel function gpu_loop_5_atomic!(b,bf_global,uh_faces,V_faces,dΩ_faces)
    face_id = @index(Global)
    if face_id <= length(dΩ_faces)
        dΩ_face = dΩ_faces[face_id]
        V_face = V_faces[dΩ_face]
        uh_face = uh_faces[dΩ_face]
        dofs = GT.dofs(V_face)
        n = GT.num_dofs(V_face)
        bf = view(bf_global, :, face_id)
        fill!(bf,0)
        for dΩ_point in GT.each_point_new(dΩ_face)
            V_point = V_face[dΩ_point]
            uh_point = uh_face[dΩ_point]
            ux = GT.field(GT.gradient,uh_point)
            sx = GT.shape_functions(GT.gradient,V_point)
            dx = GT.weight(dΩ_point)
            ux_dx = ux*dx
            for i in 1:n
                bf[i] += ux_dx⋅sx[i]
            end
        end
        for i in 1:n
            Atomix.@atomic b[dofs[i]] += bf[i]
        end
    end
end

@kernel function gpu_loop_6_numeric_ltable!(AV,V_faces,ltable)
    face_id = @index(Global)
    if face_id <= length(V_faces)
        V_face = V_faces[face_id]
        n = GT.num_dofs(V_face)
        offset = ltable[face_id]
        for V_point in GT.each_point_new(V_face)
            dx = GT.weight(V_point)
            sx = GT.shape_functions(GT.gradient,V_point)
            for j in 1:n
                sx_dx_j = sx[j]*dx
                for i in 1:n
                    AV[offset + (j-1)*n + i] += sx[i]⋅sx_dx_j
                end
            end
        end
    end
end

@kernel function gpu_loop_6_numeric_ltable_local!(AV,V_faces,ltable,::Val{max_dofs}) where {max_dofs}
    face_id = @index(Global)
    if face_id <= length(V_faces)
        A_local = zeros(MMatrix{max_dofs, max_dofs, Float64})
        V_face = V_faces[face_id]
        n = GT.num_dofs(V_face)
        offset = ltable[face_id]
        for V_point in GT.each_point_new(V_face)
            dx = GT.weight(V_point)
            sx = GT.shape_functions(GT.gradient,V_point)
            for j in 1:n
                sx_dx_j = sx[j]*dx
                for i in 1:n
                    A_local[i,j] += sx[i]⋅sx_dx_j
                end
            end
        end
        for j in 1:n
            for i in 1:n
                AV[offset + (j-1)*n + i] += A_local[i,j]
            end
        end
    end
end

function main_gpu(params)
    (;face_nodes_layout,face_dofs_layout,n) = params

    # Start at CPU
    domain = (0,1,0,1)
    cells = (n,n)
    mesh = GT.cartesian_mesh(domain,cells)
    Ω = GT.interior(mesh)
    degree = 4
    dΩ = GT.quadrature(Ω,degree)

    k = 1
    V = GT.lagrange_space(Ω,k)
    u = GT.analytical_field(f,Ω)
    uh = GT.interpolate(u,V)

    tabulate = (GT.value,GT.gradient)
    dΩ_faces_cpu = GT.each_face_new(dΩ)
    V_faces_cpu = GT.each_face_new(V,dΩ;tabulate)
    uh_faces_cpu = GT.each_face_new(uh,dΩ;tabulate)

    # Change data layout
    # This can also be done after moving memory to GPU (but not implemented yet)
    dΩ_faces_cpu = GT.change_data_layout(dΩ_faces_cpu;face_nodes_layout)
    V_faces_cpu = GT.change_data_layout(V_faces_cpu;face_nodes_layout,face_dofs_layout)
    uh_faces_cpu = GT.change_data_layout(uh_faces_cpu;face_nodes_layout,face_dofs_layout)

    # This is not needed in practice, just to make sure that
    # we do not break anything when adapting.
    dΩ_faces_gpu = adapt(dev, dΩ_faces_cpu)
    V_faces_gpu = adapt(dev, V_faces_cpu)
    uh_faces_gpu = adapt(dev, uh_faces_cpu)

    # Set the workspace location depending on the granularity
    threads_in_block = 256
    workspace_location = GT.GPUGlobalWorkspace
    # workspace_location = GT.GPUThreadWorkspace
    # workspace_location = GT.GPUSharedMemWorkspace{threads_in_block}
    dΩ_faces_gpu = GT.change_workspace_location(dΩ_faces_gpu;workspace_location)
    V_faces_gpu = GT.change_workspace_location(V_faces_gpu;workspace_location)
    uh_faces_gpu = GT.change_workspace_location(uh_faces_gpu;workspace_location)

    nfaces = length(dΩ_faces_gpu)
    nmax = GT.max_num_reference_dofs(V)
    contributions = KA.zeros(dev, Float64, nfaces)
    b_gpu = KA.zeros(dev, Float64, GT.num_free_dofs(V))
    bf_gpu = KA.zeros(dev, Float64, nmax, nfaces)

    # Launch kernel 1
    t1_gpu = @benchmark begin
        gpu_loop_1!($dev, $threads_in_block)($contributions, $dΩ_faces_gpu, ndrange=$nfaces)
        sum($contributions)
        KA.synchronize($dev)
    end
    @show sum(contributions)

    if is_cuda_available()
        blocks_in_grid = cld(nfaces, threads_in_block)
        t1_cuda = @benchmark begin
            @call_kernel cuda_loop_1 $threads_in_block $blocks_in_grid $contributions $dΩ_faces_gpu
            sum($contributions)
            CUDA.synchronize()
        end
        @show sum(contributions)
    elseif is_rocm_available()
        blocks_in_grid = cld(nfaces, threads_in_block)
        t1_hip = @benchmark begin
            @call_kernel hip_loop_1 $threads_in_block $blocks_in_grid $contributions $dΩ_faces_gpu
            sum($contributions)
            AMDGPU.synchronize()
        end
        @show sum(contributions)
    end

    println("Loop 1: KernelAbstractions throughput is ", nfaces / time(t1_gpu) * 1e9, " faces per second.")
    if is_cuda_available()
        println("Loop 1: CUDA throughput is ", nfaces / time(t1_cuda) * 1e9, " faces per second.")
        println("Loop 1: CUDA speedup is ", (nfaces / time(t1_cuda)) / (nfaces / time(t1_gpu)))
    elseif is_rocm_available()
        println("Loop 1: HIP throughput is ", nfaces / time(t1_hip) * 1e9, " faces per second.")
        println("Loop 1: HIP speedup is ", (nfaces / time(t1_hip)) / (nfaces / time(t1_gpu)))
    end

    contributions = KA.zeros(dev, Float64, nfaces)

    # Launch kernel 2
    t2_gpu = @benchmark begin
        gpu_loop_2!($dev, $threads_in_block)($contributions, $uh_faces_gpu, ndrange=$nfaces)
        sum($contributions)
        KA.synchronize($dev)
    end
    @show sum(contributions)
    if is_cuda_available()
        blocks_in_grid = cld(nfaces, threads_in_block)
        t2_cuda = @benchmark begin
            @call_kernel cuda_loop_2 $threads_in_block $blocks_in_grid $contributions $uh_faces_gpu
            sum($contributions)
            CUDA.synchronize()
        end
        @show sum(contributions)
    elseif is_rocm_available()
        blocks_in_grid = cld(nfaces, threads_in_block)
        t2_hip = @benchmark begin
            @call_kernel hip_loop_2 $threads_in_block $blocks_in_grid $contributions $uh_faces_gpu
            sum($contributions)
            AMDGPU.synchronize()
        end
        @show sum(contributions)
    end

    println("Loop 2: KernelAbstractions throughput is ", nfaces / time(t2_gpu) * 1e9, " faces per second.")
    if is_cuda_available()
        println("Loop 2: CUDA throughput is ", nfaces / time(t2_cuda) * 1e9, " faces per second.")
        println("Loop 2: CUDA speedup is ", (nfaces / time(t2_cuda) * 1e9) / (nfaces / time(t2_gpu) * 1e9))
    elseif is_rocm_available()
        println("Loop 2: HIP throughput is ", nfaces / time(t2_hip) * 1e9, " faces per second.")
        println("Loop 2: HIP speedup is ", (nfaces / time(t2_hip) * 1e9) / (nfaces / time(t2_gpu) * 1e9))
    end

    contributions = KA.zeros(dev, Float64, nfaces)

    # Launch kernel 3
    t3_gpu = @benchmark begin
        gpu_loop_3!($dev, $threads_in_block)($contributions, $uh_faces_gpu, $dΩ_faces_gpu, ndrange=$nfaces)
        sum($contributions)
        KA.synchronize($dev)
    end
    @show sum(contributions)
    if is_cuda_available()
        blocks_in_grid = cld(nfaces, threads_in_block)
        t3_cuda = @benchmark begin
            @call_kernel cuda_loop_3 $threads_in_block $blocks_in_grid $contributions $uh_faces_gpu $dΩ_faces_gpu
            sum($contributions)
            CUDA.synchronize()
        end
        @show sum(contributions)
    elseif is_rocm_available()
        blocks_in_grid = cld(nfaces, threads_in_block)
        t3_hip = @benchmark begin
            @call_kernel hip_loop_3 $threads_in_block $blocks_in_grid $contributions $uh_faces_gpu $dΩ_faces_gpu
            sum($contributions)
            AMDGPU.synchronize()
        end
        @show sum(contributions)
    end

    println("Loop 3: KernelAbstractions throughput is ", nfaces / time(t3_gpu) * 1e9, " faces per second.")
    if is_cuda_available()
        println("Loop 3: CUDA throughput is ", nfaces / time(t3_cuda) * 1e9, " faces per second.")
        println("Loop 3: CUDA speedup is ", (nfaces / time(t3_cuda) * 1e9) / (nfaces / time(t3_gpu) * 1e9))
    elseif is_rocm_available()
        println("Loop 3: HIP throughput is ", nfaces / time(t3_hip) * 1e9, " faces per second.")
        println("Loop 3: HIP speedup is ", (nfaces / time(t3_hip) * 1e9) / (nfaces / time(t3_gpu) * 1e9))
    end

    b_gpu = KA.zeros(dev, Float64, GT.num_free_dofs(V))

    # Launch kernel 4 atomic
    t4_atomic_gpu = @benchmark begin
        gpu_loop_4_atomic!($dev, $threads_in_block)($b_gpu, $uh_faces_gpu, ndrange=$nfaces)
        sqrt(sum($b_gpu.^2))
        KA.synchronize($dev)
    end setup=(fill!($b_gpu, 0.0))
    @show sqrt(sum(b_gpu.^2))
    if is_cuda_available()
        blocks_in_grid = cld(nfaces, threads_in_block)
        t4_atomic_cuda = @benchmark begin
            @call_kernel cuda_loop_4_atomic $threads_in_block $blocks_in_grid $b_gpu $uh_faces_gpu
            sqrt(sum($b_gpu.^2))
            CUDA.synchronize()
        end setup=(fill!($b_gpu, 0.0))
        @show sqrt(sum(b_gpu.^2))
    elseif is_rocm_available()
        blocks_in_grid = cld(nfaces, threads_in_block)
        t4_atomic_hip = @benchmark begin
            @call_kernel hip_loop_4_atomic $threads_in_block $blocks_in_grid $b_gpu $uh_faces_gpu
            sqrt(sum($b_gpu.^2))
            AMDGPU.synchronize()
        end setup=(fill!($b_gpu, 0.0))
        @show sqrt(sum(b_gpu.^2))
    end

    println("Loop 4 (atomic): KernelAbstractions throughput is ", nfaces / time(t4_atomic_gpu) * 1e9, " faces per second.")
    if is_cuda_available()
        println("Loop 4 (atomic): CUDA throughput is ", nfaces / time(t4_atomic_cuda) * 1e9, " faces per second.")
        println("Loop 4 (atomic): CUDA speedup is ", (nfaces / time(t4_atomic_cuda) * 1e9) / (nfaces / time(t4_atomic_gpu) * 1e9))
    elseif is_rocm_available()
        println("Loop 4 (atomic): HIP throughput is ", nfaces / time(t4_atomic_hip) * 1e9, " faces per second.")
        println("Loop 4 (atomic): HIP speedup is ", (nfaces / time(t4_atomic_hip) * 1e9) / (nfaces / time(t4_atomic_gpu) * 1e9))
    end

    b_gpu = KA.zeros(dev, Float64, GT.num_free_dofs(V))
    bf_gpu = KA.zeros(dev, Float64, nmax, nfaces)

    # Launch kernel 4 global
    t4_global_gpu = @benchmark begin
        gpu_loop_4_global!($dev, $threads_in_block)($b_gpu, $bf_gpu, $uh_faces_gpu, ndrange=$nfaces)
        sqrt(sum($b_gpu.^2))
        KA.synchronize($dev)
    end setup=(fill!($b_gpu, 0.0))
    @show sqrt(sum(b_gpu.^2))
    if is_cuda_available()
        blocks_in_grid = cld(nfaces, threads_in_block)
        t4_global_cuda = @benchmark begin
            @call_kernel cuda_loop_4_global $threads_in_block $blocks_in_grid $b_gpu $bf_gpu $uh_faces_gpu
            sqrt(sum($b_gpu.^2))
            CUDA.synchronize()
        end setup=(fill!($b_gpu, 0.0))
        @show sqrt(sum(b_gpu.^2))
    elseif is_rocm_available()
        blocks_in_grid = cld(nfaces, threads_in_block)
        t4_global_hip = @benchmark begin
            @call_kernel hip_loop_4_global $threads_in_block $blocks_in_grid $b_gpu $bf_gpu $uh_faces_gpu
            sqrt(sum($b_gpu.^2))
            AMDGPU.synchronize()
        end setup=(fill!($b_gpu, 0.0))
        @show sqrt(sum(b_gpu.^2))
    end

    println("Loop 4 (global): KernelAbstractions throughput is ", nfaces / time(t4_global_gpu) * 1e9, " faces per second.")
    if is_cuda_available()
        println("Loop 4 (global): CUDA throughput is ", nfaces / time(t4_global_cuda) * 1e9, " faces per second.")
        println("Loop 4 (global): CUDA speedup is ", (nfaces / time(t4_global_cuda) * 1e9) / (nfaces / time(t4_global_gpu) * 1e9))
    elseif is_rocm_available()
        println("Loop 4 (global): HIP throughput is ", nfaces / time(t4_global_hip) * 1e9, " faces per second.")
        println("Loop 4 (global): HIP speedup is ", (nfaces / time(t4_global_hip) * 1e9) / (nfaces / time(t4_global_gpu) * 1e9))
    end

    b_gpu = KA.zeros(dev, Float64, GT.num_free_dofs(V))

    # Launch kernel 4 local
    t4_local_gpu = @benchmark begin
        gpu_loop_4_local!($dev, $threads_in_block)($b_gpu, Val($nmax), $uh_faces_gpu, ndrange=$nfaces)
        sqrt(sum($b_gpu.^2))
        KA.synchronize($dev)
    end setup=(fill!($b_gpu, 0.0))
    @show sqrt(sum(b_gpu.^2))
    if is_cuda_available()
        blocks_in_grid = cld(nfaces, threads_in_block)
        t4_local_cuda = @benchmark begin
            @call_kernel cuda_loop_4_local $threads_in_block $blocks_in_grid $b_gpu Val($nmax) $uh_faces_gpu
            sqrt(sum($b_gpu.^2))
            CUDA.synchronize()
        end setup=(fill!($b_gpu, 0.0))
        @show sqrt(sum(b_gpu.^2))
    elseif is_rocm_available()
        blocks_in_grid = cld(nfaces, threads_in_block)
        t4_local_hip = @benchmark begin
            @call_kernel hip_loop_4_local $threads_in_block $blocks_in_grid $b_gpu Val($nmax) $uh_faces_gpu
            sqrt(sum($b_gpu.^2))
            AMDGPU.synchronize()
        end setup=(fill!($b_gpu, 0.0))
        @show sqrt(sum(b_gpu.^2))
    end

    println("Loop 4 (local): KernelAbstractions throughput is ", nfaces / time(t4_local_gpu) * 1e9, " faces per second.")
    if is_cuda_available()
        println("Loop 4 (local): CUDA throughput is ", nfaces / time(t4_local_cuda) * 1e9, " faces per second.")
        println("Loop 4 (local): CUDA speedup is ", (nfaces / time(t4_local_cuda) * 1e9) / (nfaces / time(t4_local_gpu) * 1e9))
    elseif is_rocm_available()
        println("Loop 4 (local): HIP throughput is ", nfaces / time(t4_local_hip) * 1e9, " faces per second.")
        println("Loop 4 (local): HIP speedup is ", (nfaces / time(t4_local_hip) * 1e9) / (nfaces / time(t4_local_gpu) * 1e9))
    end

    shmem_needed = sizeof(Float64) * nmax * threads_in_block
    shmem_max = max_shared_memory_per_block()
    if shmem_needed <= shmem_max
        b_gpu = KA.zeros(dev, Float64, GT.num_free_dofs(V))

        # Launch kernel 4 shared
        t4_shared_gpu = @benchmark begin
            gpu_loop_4_shared!($dev, $threads_in_block)($b_gpu, Val($nmax), Val($threads_in_block), $uh_faces_gpu, ndrange=$nfaces)
            sqrt(sum($b_gpu.^2))
            KA.synchronize($dev)
        end setup=(fill!($b_gpu, 0.0))
        @show sqrt(sum(b_gpu.^2))
        if is_cuda_available()
            blocks_in_grid = cld(nfaces, threads_in_block)
            t4_shared_cuda = @benchmark begin
                @call_kernel cuda_loop_4_shared $threads_in_block $blocks_in_grid $b_gpu Val($nmax) Val($threads_in_block) $uh_faces_gpu
                sqrt(sum($b_gpu.^2))
                CUDA.synchronize()
            end setup=(fill!($b_gpu, 0.0))
            @show sqrt(sum(b_gpu.^2))
        elseif is_rocm_available()
            blocks_in_grid = cld(nfaces, threads_in_block)
            shmem = sizeof(Float64) * nmax * threads_in_block
            t4_shared_hip = @benchmark begin
                @call_kernel_shmem hip_loop_4_shared $threads_in_block $blocks_in_grid $shmem $b_gpu Val($nmax) Val($threads_in_block) $uh_faces_gpu
                sqrt(sum($b_gpu.^2))
                AMDGPU.synchronize()
            end setup=(fill!($b_gpu, 0.0))
            @show sqrt(sum(b_gpu.^2))
        end
        println("Loop 4 (shared): KernelAbstractions throughput is ", nfaces / time(t4_shared_gpu) * 1e9, " faces per second.")
        if is_cuda_available()
            println("Loop 4 (shared): CUDA throughput is ", nfaces / time(t4_shared_cuda) * 1e9, " faces per second.")
            println("Loop 4 (shared): CUDA speedup is ", (nfaces / time(t4_shared_cuda) * 1e9) / (nfaces / time(t4_shared_gpu) * 1e9))
        elseif is_rocm_available()
            println("Loop 4 (shared): HIP throughput is ", nfaces / time(t4_shared_hip) * 1e9, " faces per second.")
            println("Loop 4 (shared): HIP speedup is ", (nfaces / time(t4_shared_hip) * 1e9) / (nfaces / time(t4_shared_gpu) * 1e9))
        end
    else
        println("Loop 4 (shared): skipped (requires $(shmem_needed) bytes, max is $(shmem_max) bytes per block)")
    end

    b_gpu = KA.zeros(dev, Float64, GT.num_free_dofs(V))
    bf_gpu = KA.zeros(dev, Float64, nmax, nfaces)

    # Launch kernel 5 atomic
    t5_atomic_gpu = @benchmark begin
        gpu_loop_5_atomic!($dev, $threads_in_block)($b_gpu, $bf_gpu, $uh_faces_gpu, $V_faces_gpu, $dΩ_faces_gpu, ndrange=$nfaces)
        sqrt(sum($b_gpu.^2))
        KA.synchronize($dev)
    end setup=(fill!($b_gpu, 0.0))
    @show sqrt(sum(b_gpu.^2))
    if is_cuda_available()
        blocks_in_grid = cld(nfaces, threads_in_block)
        t5_atomic_cuda = @benchmark begin
            @call_kernel cuda_loop_5_atomic $threads_in_block $blocks_in_grid $b_gpu $bf_gpu $uh_faces_gpu $V_faces_gpu $dΩ_faces_gpu
            sqrt(sum($b_gpu.^2))
            CUDA.synchronize()
        end setup=(fill!($b_gpu, 0.0))
        @show sqrt(sum(b_gpu.^2))
    elseif is_rocm_available()
        blocks_in_grid = cld(nfaces, threads_in_block)
        t5_atomic_hip = @benchmark begin
            @call_kernel hip_loop_5_atomic $threads_in_block $blocks_in_grid $b_gpu $bf_gpu $uh_faces_gpu $V_faces_gpu $dΩ_faces_gpu
            sqrt(sum($b_gpu.^2))
            AMDGPU.synchronize()
        end setup=(fill!($b_gpu, 0.0))
        @show sqrt(sum(b_gpu.^2))
    end
    println("Loop 5 (atomic): KernelAbstractions throughput is ", nfaces / time(t5_atomic_gpu) * 1e9, " faces per second.")
    if is_cuda_available()
        println("Loop 5 (atomic): CUDA throughput is ", nfaces / time(t5_atomic_cuda) * 1e9, " faces per second.")
        println("Loop 5 (atomic): CUDA speedup is ", (nfaces / time(t5_atomic_cuda) * 1e9) / (nfaces / time(t5_atomic_gpu) * 1e9))
    elseif is_rocm_available()
        println("Loop 5 (atomic): HIP throughput is ", nfaces / time(t5_atomic_hip) * 1e9, " faces per second.")
        println("Loop 5 (atomic): HIP speedup is ", (nfaces / time(t5_atomic_hip) * 1e9) / (nfaces / time(t5_atomic_gpu) * 1e9))
    end

    b = zeros(GT.num_free_dofs(V))
    num_nz = cpu_loop_6_count(V_faces_cpu)
    n_global = GT.num_dofs(V)
    x = GT.free_values(uh)
    AI = zeros(Int32,num_nz)
    AJ = zeros(Int32,num_nz)
    AV_gpu = KA.zeros(dev, Float64, num_nz)
    ltable_cpu = zeros(Int32,length(V_faces_cpu))
    ltable_gpu = KA.zeros(dev, Int32, length(V_faces_cpu))
    cpu_loop_6_ltable!(ltable_cpu,V_faces_cpu)
    KA.copy!(ltable_gpu, ltable_cpu)
    cpu_loop_6_symbolic!(AI,AJ,V_faces_cpu)

    # Launch kernel 6 numeric lookup table
    t6_numeric_ltable_gpu = @benchmark begin
        gpu_loop_6_numeric_ltable!($dev, $threads_in_block)($AV_gpu, $V_faces_gpu, $ltable_gpu, ndrange=$nfaces)
        KA.synchronize($dev)
    end setup=(fill!($AV_gpu, 0.0))
    AV_cpu = Array(AV_gpu)
    A,Acache = PA.sparse_matrix(AI,AJ,AV_cpu,n_global,n_global;reuse=Val(true))
    PA.sparse_matrix!(A,AV_cpu,Acache)
    b = A*x
    @show norm(b)

    if is_cuda_available()
        blocks_in_grid = cld(nfaces, threads_in_block)
        t6_numeric_ltable_cuda = @benchmark begin
            @call_kernel cuda_loop_6_numeric_ltable $threads_in_block $blocks_in_grid $AV_gpu $V_faces_gpu $ltable_gpu
            CUDA.synchronize()
        end setup=(fill!($AV_gpu, 0.0))
        AV_cpu = Array(AV_gpu)
        A,Acache = PA.sparse_matrix(AI,AJ,AV_cpu,n_global,n_global;reuse=Val(true))
        PA.sparse_matrix!(A,AV_cpu,Acache)
        b = A*x
        @show norm(b)
    elseif is_rocm_available()
        blocks_in_grid = cld(nfaces, threads_in_block)
        t6_numeric_ltable_hip = @benchmark begin
            @call_kernel hip_loop_6_numeric_ltable $threads_in_block $blocks_in_grid $AV_gpu $V_faces_gpu $ltable_gpu
            AMDGPU.synchronize()
        end setup=(fill!($AV_gpu, 0.0))
        AV_cpu = Array(AV_gpu)
        A,Acache = PA.sparse_matrix(AI,AJ,AV_cpu,n_global,n_global;reuse=Val(true))
        PA.sparse_matrix!(A,AV_cpu,Acache)
        b = A*x
        @show norm(b)
    end

    println("Loop 6 (numeric lookup table): KernelAbstractions throughput is ", nfaces / time(t6_numeric_ltable_gpu) * 1e9, " faces per second.")
    if is_cuda_available()
        println("Loop 6 (numeric lookup table): CUDA throughput is ", nfaces / time(t6_numeric_ltable_cuda) * 1e9, " faces per second.")
        println("Loop 6 (numeric lookup table): CUDA speedup is ", (nfaces / time(t6_numeric_ltable_cuda) * 1e9) / (nfaces / time(t6_numeric_ltable_gpu) * 1e9))
    elseif is_rocm_available()
        println("Loop 6 (numeric lookup table): HIP throughput is ", nfaces / time(t6_numeric_ltable_hip) * 1e9, " faces per second.")
        println("Loop 6 (numeric lookup table): HIP speedup is ", (nfaces / time(t6_numeric_ltable_hip) * 1e9) / (nfaces / time(t6_numeric_ltable_gpu) * 1e9))
    end

    b = zeros(GT.num_free_dofs(V))
    num_nz = cpu_loop_6_count(V_faces_cpu)
    n_global = GT.num_dofs(V)
    x = GT.free_values(uh)
    AI = zeros(Int32,num_nz)
    AJ = zeros(Int32,num_nz)
    AV_gpu = KA.zeros(dev, Float64, num_nz)
    ltable_cpu = zeros(Int32,length(V_faces_cpu))
    ltable_gpu = KA.zeros(dev, Int32, length(V_faces_cpu))
    cpu_loop_6_ltable!(ltable_cpu,V_faces_cpu)
    KA.copy!(ltable_gpu, ltable_cpu)
    cpu_loop_6_symbolic!(AI,AJ,V_faces_cpu)

    # Launch kernel 6 numeric lookup table local
    t6_numeric_ltable_local_gpu = @benchmark begin
        gpu_loop_6_numeric_ltable_local!($dev, $threads_in_block)($AV_gpu, $V_faces_gpu, $ltable_gpu, Val($nmax), ndrange=$nfaces)
        KA.synchronize($dev)
    end setup=(fill!($AV_gpu, 0.0))
    AV_cpu = Array(AV_gpu)
    A,Acache = PA.sparse_matrix(AI,AJ,AV_cpu,n_global,n_global;reuse=Val(true))
    PA.sparse_matrix!(A,AV_cpu,Acache)
    b = A*x
    @show norm(b)
    if is_cuda_available()
        blocks_in_grid = cld(nfaces, threads_in_block)
        t6_numeric_ltable_local_cuda = @benchmark begin
            @call_kernel cuda_loop_6_numeric_ltable_local $threads_in_block $blocks_in_grid $AV_gpu $V_faces_gpu $ltable_gpu Val($nmax)
            CUDA.synchronize()
        end setup=(fill!($AV_gpu, 0.0))
        AV_cpu = Array(AV_gpu)
        A,Acache = PA.sparse_matrix(AI,AJ,AV_cpu,n_global,n_global;reuse=Val(true))
        PA.sparse_matrix!(A,AV_cpu,Acache)
        b = A*x
        @show norm(b)
    elseif is_rocm_available()
        blocks_in_grid = cld(nfaces, threads_in_block)
        t6_numeric_ltable_local_hip = @benchmark begin
            @call_kernel hip_loop_6_numeric_ltable_local $threads_in_block $blocks_in_grid $AV_gpu $V_faces_gpu $ltable_gpu Val($nmax)
            AMDGPU.synchronize()
        end setup=(fill!($AV_gpu, 0.0))
        AV_cpu = Array(AV_gpu)
        A,Acache = PA.sparse_matrix(AI,AJ,AV_cpu,n_global,n_global;reuse=Val(true))
        PA.sparse_matrix!(A,AV_cpu,Acache)
        b = A*x
        @show norm(b)
    end

    println("Loop 6 (numeric lookup table local): KernelAbstractions throughput is ", nfaces / time(t6_numeric_ltable_local_gpu) * 1e9, " faces per second.")
    if is_cuda_available()
        println("Loop 6 (numeric lookup table local): CUDA throughput is ", nfaces / time(t6_numeric_ltable_local_cuda) * 1e9, " faces per second.")
        println("Loop 6 (numeric lookup table local): CUDA speedup is ", (nfaces / time(t6_numeric_ltable_local_cuda) * 1e9) / (nfaces / time(t6_numeric_ltable_local_gpu) * 1e9))
    elseif is_rocm_available()
        println("Loop 6 (numeric lookup table local): HIP throughput is ", nfaces / time(t6_numeric_ltable_local_hip) * 1e9, " faces per second.")
        println("Loop 6 (numeric lookup table local): HIP speedup is ", (nfaces / time(t6_numeric_ltable_local_hip) * 1e9) / (nfaces / time(t6_numeric_ltable_local_gpu) * 1e9))
    end

    b = zeros(GT.num_free_dofs(V))
    num_nz = cpu_loop_6_count(V_faces_cpu)
    n_global = GT.num_dofs(V)
    x = GT.free_values(uh)
    AI = zeros(Int32,num_nz)
    AJ = zeros(Int32,num_nz)
    AV_gpu = KA.zeros(dev, Float64, num_nz)
    ltable_cpu = zeros(Int32,length(V_faces_cpu))
    ltable_gpu = KA.zeros(dev, Int32, length(V_faces_cpu))
    cpu_loop_6_ltable!(ltable_cpu,V_faces_cpu)
    KA.copy!(ltable_gpu, ltable_cpu)
    cpu_loop_6_symbolic!(AI,AJ,V_faces_cpu)
end

for n in [10, 100, 1000, 10000]
    layouts = (GT.face_minor_array,GT.face_major_array)
    for face_dofs_layout in layouts
        for face_nodes_layout in layouts
            params = (;face_nodes_layout,face_dofs_layout,n)
            println("n = ", n)
            println("face_dofs_layout = ", face_dofs_layout)
            println("face_nodes_layout = ", face_nodes_layout)
            println()

            println("Running CPU version...")
            main_cpu(params)
            println()

            println("Running GPU version...")
            main_gpu(params)
            println()
        end
    end
    println()
end

end # module

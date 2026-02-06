module GPUNewTests

using CUDA
using Test
using LinearAlgebra
using SparseArrays
import Adapt
import PartitionedArrays as PA
import GalerkinToolkit as GT

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

function main_cpu(params)
    (;face_nodes_layout,face_dofs_layout) = params

    # Start at CPU
    domain = (0,1,0,1)
    cells = (4,4)
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
    dΩ_faces_cpu = GT.change_data_layout(dΩ_faces_cpu;face_nodes_layout)
    V_faces_cpu = GT.change_data_layout(V_faces_cpu;face_nodes_layout,face_dofs_layout)
    uh_faces_cpu = GT.change_data_layout(uh_faces_cpu;face_nodes_layout,face_dofs_layout)

    # This is not needed in practice, just to make sure that
    # we do not break anything when adapting.
    dΩ_faces_cpu = Adapt.adapt_structure(Array,dΩ_faces_cpu)
    V_faces_cpu = Adapt.adapt_structure(Array,V_faces_cpu)
    uh_faces_cpu = Adapt.adapt_structure(Array,uh_faces_cpu)

    #dΩ_face_gpu = CUDA.cu(dΩ_faces_cpu)
    #dΩ_face_gpu = GT.loop_options(dΩ_face_gpu;
    #    face_dofs_layout=:face_major, # :face_minor
    #    granularity=:face_per_thread, # :face_per_block
    #    shape_functions_location =:global_memory, # :shared_memory, :kernel_memory
    #   )

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
    fill(AV,0)
    cpu_loop_6_numeric!(AV,Af,V_faces_cpu)
    PA.sparse_matrix!(A,AV,Acache)
    b = A*x
    @show norm(b)

end

function kernel_1!(contributions,dΩ_faces_gpu)
    face_id = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if face_id > length(dΩ_faces_gpu)
        return nothing
    end
    dΩ_face = dΩ_faces_gpu[face_id]
    s = 0.0
    for dΩ_point in GT.each_point_new(dΩ_face)
        x = GT.coordinate(dΩ_point)
        dx = GT.weight(dΩ_point)
        s += f(x)*dx
    end
    contributions[face_id] = s
    return nothing
end

function kernel_2!(contributions,uh_faces_gpu)
    face_id = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if face_id > length(uh_faces_gpu)
        return nothing
    end
    uh_face = uh_faces_gpu[face_id]
    s = 0.0
    for uh_point in GT.each_point_new(uh_face)
        ux = GT.field(GT.value,uh_point)
        dx = GT.weight(uh_point)
        s += ux*dx
    end
    contributions[face_id] = s
    return nothing
end

function kernel_5!()
    #TODO
end

function main_gpu(params)
    (;face_nodes_layout,face_dofs_layout) = params

    # Start at CPU
    domain = (0,1,0,1)
    cells = (4,4)
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
    dΩ_faces_gpu = CUDA.cu(dΩ_faces_cpu)
    V_faces_gpu = CUDA.cu(V_faces_cpu)
    uh_faces_gpu = CUDA.cu(uh_faces_cpu)

    #TODO
    #granularity = Val(:face_per_thread) # Val(:face_per_block)
    #dΩ_faces_cpu = GT.change_loop_granularity(dΩ_faces_cpu,granularity)
    #V_faces_cpu = GT.change_loop_granularity(V_faces_cpu,granularity)
    #uh_faces_cpu = GT.change_loop_granularity(uh_faces_cpu,granularity)
    #
    #Maybe workspace location should not be independent and depend on granularity
    #workspace_location = Val(:global_memory) # Val(:shared_memory) Val(:thread_memory)
    #dΩ_faces_gpu = GT.change_workspace_location(dΩ_faces_gpu,workspace_location)
    #V_faces_gpu = GT.change_workspace_location(V_faces_gpu,workspace_location)
    #uh_faces_gpu = GT.change_workspace_location(uh_faces_gpu,workspace_location)

    nfaces = length(dΩ_faces_gpu)
    contributions = CUDA.zeros(Float64,nfaces)

    # Launch kernel 1
    threads_in_block = 256
    blocks_in_grid = ceil(Int, nfaces/256)
    @cuda threads=threads_in_block blocks=blocks_in_grid kernel_1!(contributions,dΩ_faces_gpu)
    @show r_gpu = sum(contributions)

    # Launch kernel 2
    threads_in_block = 256
    blocks_in_grid = ceil(Int, nfaces/256)
    @cuda threads=threads_in_block blocks=blocks_in_grid kernel_2!(contributions,uh_faces_gpu)
    @show r_gpu = sum(contributions)

end


layouts = (GT.face_minor_array,GT.face_major_array)
for face_dofs_layout in layouts
    for face_nodes_layout in layouts
        params = (;face_nodes_layout,face_dofs_layout)
        main_cpu(params)
        main_gpu(params)
    end
end




end # module

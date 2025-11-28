module GPUTests

using CUDA
using Test
import GalerkinToolkit as GT

# Goal integrate function f(x)
f(x) = 2*sin(sum(x))
# on a domain Ω defined by a computational mesh
# (a square in this example)

# Start at CPU
domain = (0,1,0,1)
cells = (4,4)
mesh = GT.cartesian_mesh(domain,cells)
Ω = GT.interior(mesh)
degree = 4
dΩ = GT.quadrature(Ω,degree)
dΩ_faces = GT.each_face(dΩ)

# This is still on CPU, but with a data
# layout more appropiate for GPUs.
dΩ_faces_cpu = GT.device_layout(dΩ_faces)

# Now, move data to GPU
dΩ_faces_gpu = CUDA.cu(dΩ_faces_cpu)

nfaces = length(dΩ_faces_gpu)
contributions = CUDA.zeros(Float64,nfaces)

function kernel!(contributions,dΩ_faces_gpu)
    face_id = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if face_id > length(dΩ_faces_gpu)
        return nothing
    end
    dΩ_face = dΩ_faces_gpu[face_id]
    s = 0.0
    for dΩ_point in GT.each_point(dΩ_face)
        x = GT.coordinate(dΩ_point)
        dx = GT.weight(dΩ_point)
        s += f(x)*dx
    end
    contributions[face_id] = s
    return nothing
end

# Launch kernel
threads_in_block = 256
blocks_in_grid = ceil(Int, nfaces/256)
@cuda threads=threads_in_block blocks=blocks_in_grid kernel!(contributions,dΩ_faces_gpu)

@show r_gpu = sum(contributions)

# For reference, this is how you can do this on the CPU
function cpu_loop(dΩ_faces)
    s = 0.0
    for dΩ_face in dΩ_faces
        for dΩ_point in GT.each_point(dΩ_face)
            x = GT.coordinate(dΩ_point)
            dx = GT.weight(dΩ_point)
            s += f(x)*dx
        end
    end
    return s
end
@show r_cpu = cpu_loop(dΩ_faces_cpu)
@test r_gpu≈r_cpu

@show r_cpu = cpu_loop(dΩ_faces)
@test r_gpu≈r_cpu

end # module

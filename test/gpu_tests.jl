module GPUTests

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
# layout more appropriate for GPUs.
dΩ_faces_cpu = GT.device_layout(dΩ_faces)

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
@show cpu_loop(dΩ_faces_cpu)

end # module

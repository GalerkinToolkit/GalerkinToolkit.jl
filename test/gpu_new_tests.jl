module GPUNewTests

using Test
using LinearAlgebra
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

k = 1
V = GT.lagrange_space(Ω,k)
u = GT.analytical_field(f,Ω)
uh = GT.interpolate(u,V)

dΩ_faces_cpu = GT.each_face_new(dΩ)

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
@show r_cpu = cpu_loop_1(dΩ_faces_cpu)

# 0-form for discrete field
uh_faces_cpu = GT.each_face_new(uh,dΩ;tabulate=(GT.value,))
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
@show r_cpu = cpu_loop_2(uh_faces_cpu)

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
@show r_cpu = cpu_loop_3(uh_faces_cpu,dΩ_faces_cpu)

# 1-form equivalent to matrix free
b = zeros(GT.num_free_dofs(V))
bf = zeros(GT.max_num_reference_dofs(V))
function cpu_loop_4!(b,bf,uh_faces)
    for uh_face in uh_faces
        dofs = GT.dofs(uh_face)
        n = GT.num_dofs(uh_face)
        fill!(bf,0)
        for uh_point in GT.each_point_new(uh_face)
            ux = GT.field(GT.value,uh_point)
            sx = GT.shape_functions(GT.value,uh_point)
            dx = GT.weight(uh_point)
            ux_dx = ux*dx
            for i in 1:n
                bf[i] += ux_dx*sx[i]
            end
        end
        b[dofs] += bf
    end
end
r_cpu = cpu_loop_4!(b,bf,uh_faces_cpu)
@show norm(b)

# 1-form equivalent to matrix free
# where the discrete field
# and trial functions are potentially different
fill!(b,0)
V_faces_cpu = GT.each_face_new(V,dΩ;tabulate=(GT.value,))
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
            ux = GT.field(GT.value,uh_point)
            sx = GT.shape_functions(GT.value,V_point)
            dx = GT.weight(dΩ_point)
            ux_dx = ux*dx
            for i in 1:n
                bf[i] += ux_dx*sx[i]
            end
        end
        b[dofs] += bf
    end
end
r_cpu = cpu_loop_5!(b,bf,uh_faces_cpu,V_faces_cpu,dΩ_faces_cpu)
@show norm(b)


end # module

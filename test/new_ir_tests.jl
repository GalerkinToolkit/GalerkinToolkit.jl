module NewIRTests

using Test
import GalerkinToolkit as GT
import PartitionedArrays as pa
using GalerkinToolkit: ∫
using LinearAlgebra
import ForwardDiff
using AbstractTrees

function hand_written_baseline(α,dΩ)
    v = GT.term(α, nothing).value
    face_point_dV = GT.weight_accessor(dΩ)
    face_point_J = GT.jacobian_accessor(dΩ)
    face_npoints = GT.num_points_accessor(dΩ)
    nfaces = GT.num_faces(GT.domain(dΩ))
    s = zero(v)
    for face in 1:nfaces
        point_dV = face_point_dV(face)
        point_J = face_point_J(face)
        npoints = face_npoints(face)
        for point in 1:npoints
            J = point_J(point)
            dV = point_dV(point,J)
            s += v*dV
        end
    end
    s
end

outdir = mkpath(joinpath(@__DIR__,"..","output"))

n = 10
domain = (0,2,0,2,0,2)
cells = (n,n,n)
mesh = GT.cartesian_mesh(domain,cells)
D = GT.num_dims(mesh)
Ω = GT.interior(mesh)
degree = 3
dΩ = GT.new_measure(Ω,degree)

α = GT.uniform_quantity(1.0)

@time v = hand_written_baseline(α,dΩ)
@time v = hand_written_baseline(α,dΩ)
@test v ≈ 8

#using InteractiveUtils
#@code_warntype hand_written_baseline(α,dΩ)

int = GT.contribution(x->α,dΩ)

expr = GT.generate(int,α,dΩ)
#display(expr)

f = GT.evaluate(expr)

domain_face_v = f(α,dΩ)

@time s = sum(domain_face_v, 1:GT.num_faces(Ω))
@time s = sum(domain_face_v, 1:GT.num_faces(Ω))

@test s ≈ 8


int = GT.contribution(x->α,dΩ)

expr = GT.generate(int,α,dΩ)

f = GT.evaluate(expr)

@time domain_face_v = f(α,dΩ)

domain_face = 2
sum(domain_face_v, 1:1)
@time s = sum(domain_face_v, 1:GT.num_faces(Ω))
@show s


expr_partial_param = GT.generate(int, α)
f_partial_param = GT.evaluate(expr_partial_param)

α = GT.uniform_quantity(2.0)
domain_face_v = f_partial_param(α)
sum(domain_face_v, 1:1)
@time s = sum(domain_face_v, 1:GT.num_faces(Ω))
@test s ≈ 16


@time domain_face_v = f(α,dΩ)
domain_face_v(1)
@time s = sum(domain_face_v, 1:GT.num_faces(Ω))
@show s
@test s ≈ 16


expr_noparam = GT.generate(int)
f_noparam = GT.evaluate(expr_noparam)
domain_face_v = f_noparam()
sum(domain_face_v, 1:1)
@time s = sum(domain_face_v, 1:GT.num_faces(Ω))
@test s ≈ 8

g = GT.uniform_quantity(sum)
int = GT.new_∫(GT.@qty(x->α * sum(x) + g(x) * (x[1] + 2.0 * x[2] + 3.0 * x[3])),dΩ)
expr_sum = GT.generate(int)
f_sum = GT.evaluate(expr_sum)
domain_face_v = f_sum()
sum(domain_face_v, 1:1)
@time s = sum(domain_face_v, 1:GT.num_faces(Ω))
@test s ≈ 208

domain = (0,1,0,1)
cells = (2,2)
mesh = GT.cartesian_mesh(domain,cells)
D = GT.num_dims(mesh)
Ω = GT.interior(mesh)
degree = 3
dΩ = GT.new_measure(Ω,degree)

@time domain_face_v = f(α,dΩ)
sum(domain_face_v, 1:1)
@time s = sum(domain_face_v, 1:GT.num_faces(Ω))
@show s
@test s ≈ 2

#face_point_J
#face_point_w 
#loop over faces, loop over pointsm and and up the output of face_point_w(face)(point,J)

#expr = :(sum(i->i*sum(j->i*j),a:sum(k->b*k,r)))
#expr = :(sum(k->(b+c)*k,r))
#display(expr)
#block = GT.statements(expr)
#display(block)
#
#expr = :(sum(k->(b+sum(i->(b+i+k),r))*k,r))
#display(expr)
#block = GT.statements(expr)
#display(block)
#
#expr = :(sum(i->i*sum(j->i*j),a:sum(k->b*k,r)))
#display(expr)
#block = GT.statements(expr)
#display(block)



#Ωref = GT.interior(mesh;is_reference_domain=true)
#ϕ = GT.physical_map(mesh,D)
#u = GT.analytical_field(x->sum(x),Ω)
#degree = 2
#
#
#
#dΩref = GT.new_measure(Ωref,degree)
#dΩref_old = GT.measure(Ωref,degree)
#
#
#int_old = ∫(dΩref_old) do q
#    1.0
#end
#sum_old = sum(int_old)
#
#int_new = () -> ∫(dΩref) do q
#    1.0
#end 
#f = GT.generate_sum(int_new)
#sum_new = f()
#
#@test sum_old ≈ sum_new
#
#
#dΩ = GT.new_measure(Ω,degree)
#dΩ_old = GT.measure(Ω,degree)
#int_old = ∫(dΩ_old) do q
#    1.0
#end
#sum_old = sum(int_old)
#
#sum_new = 4.0
#int_new = () -> ∫(dΩ) do q
#    1.0
#end 
#f = GT.generate_sum(int_new)
#sum_new = f()
#
#@test sum_old ≈ sum_new
#
#
#
#
#
#int_new = (a) -> ∫(dΩ) do q
#    a
#end 
#f = GT.generate_sum(int_new, 2.0)
#@test f(3.0) ≈ 12.0
#
#int_new = (a) -> ∫(dΩref) do q
#    a
#end 
#f = GT.generate_sum(int_new, 2.0)
#@test f(5.0) ≈ 320.


end # module

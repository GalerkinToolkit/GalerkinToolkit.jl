module NewIRTests

using Test
import GalerkinToolkit as GT
import PartitionedArrays as pa
using GalerkinToolkit: ∫
using LinearAlgebra
import ForwardDiff
using AbstractTrees
outdir = mkpath(joinpath(@__DIR__,"..","output"))


domain = (0,2,0,2)
cells = (8,8)
mesh = GT.cartesian_mesh(domain,cells)
D = GT.num_dims(mesh)
Ω = GT.interior(mesh)
degree = 3
dΩ = GT.new_measure(Ω,degree)

α = GT.uniform_quantity(1.0)

int = GT.contribution(x->α,dΩ)

expr = GT.generate(int,α,dΩ)

f = eval(expr)

@time domain_face_v = f(α,dΩ)

domain_face = 2
v = domain_face_v(domain_face)
@show v
domain_face_v(1)
@time s = sum(domain_face->domain_face_v(domain_face), 1:GT.num_faces(Ω))
@show s

α = GT.uniform_quantity(2.0)
@time domain_face_v = f(α,dΩ)
domain_face_v(1)
@time s = sum(domain_face->domain_face_v(domain_face), 1:GT.num_faces(Ω))
@show s

domain = (0,1,0,1)
cells = (2,2)
mesh = GT.cartesian_mesh(domain,cells)
D = GT.num_dims(mesh)
Ω = GT.interior(mesh)
degree = 3
dΩ = GT.new_measure(Ω,degree)

@time domain_face_v = f(α,dΩ)
domain_face_v(1)
@time s = sum(domain_face->domain_face_v(domain_face), 1:GT.num_faces(Ω))
@show s

#face_point_J
#face_point_w 
#loop over faces, loop over pointsm and and up the output of face_point_w(face)(point,J)


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

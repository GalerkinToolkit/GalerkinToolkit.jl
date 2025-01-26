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
Ωref = GT.interior(mesh;is_reference_domain=true)
ϕ = GT.physical_map(mesh,D)
u = GT.analytical_field(x->sum(x),Ω)
degree = 2
dΩref = GT.new_measure(Ωref,degree)
dΩref_old = GT.measure(Ωref,degree)


int_old = ∫(dΩref_old) do q
    1.0
end
sum_old = sum(int_old)

int_new = () -> ∫(dΩref) do q
    1.0
end 
f = GT.generate_sum(int_new)
sum_new = f()

@test sum_old ≈ sum_new


dΩ = GT.new_measure(Ω,degree)
dΩ_old = GT.measure(Ω,degree)
int_old = ∫(dΩ_old) do q
    1.0
end
sum_old = sum(int_old)

sum_new = 4.0
int_new = () -> ∫(dΩ) do q
    1.0
end 
f = GT.generate_sum(int_new)
sum_new = f()

@test sum_old ≈ sum_new





int_new = (a) -> ∫(dΩ) do q
    a
end 
f = GT.generate_sum(int_new, 2.0)
@test f(3.0) ≈ 12.0

int_new = (a) -> ∫(dΩref) do q
    a
end 
f = GT.generate_sum(int_new, 2.0)
@test f(5.0) ≈ 320.


end # module

module FieldsTests

import GalerkinToolkit as GT
using GalerkinToolkit: ∫, ×
using Test
import ForwardDiff
using LinearAlgebra
using BlockArrays
using WriteVTK
using StaticArrays
using PartitionedArrays

#using InteractiveUtils

domain = (0,1,0,1)
cells_per_dir = (4,4)
parts_per_dir = (2,2)
mesh = GT.cartesian_mesh(domain,cells_per_dir)

Ω = GT.interior(mesh)
Γ = GT.boundary(mesh)
Λ = GT.skeleton(mesh)
F = GT.domain(mesh,1)
dΩ = GT.quadrature(Ω,2)

order = 2

V = GT.lagrange_space(Γ,order)
uh = GT.rand_field(Float64,V)
plt = GT.plot(F)
GT.plot!(plt,uh,label="uh")
vtk_grid("plt",plt) |> close

V = GT.lagrange_space(Ω,order;dirichlet_boundary=Γ)

g = GT.analytical_field(Ω;piecewise=true) do x,name
    sum(x)
end

gx = GT.sample(g,dΩ)

g = GT.analytical_field(sum,Ω)
GT.interpolate(g,V)
GT.interpolate_free(g,V)
GT.interpolate_dirichlet(g,V)



VxV = V × V

v1xv2 = GT.zero_dirichlet_field(Float64,VxV)
v1,v2 = v1xv2

x = GT.free_values(v1xv2)
display(x)

x = GT.dirichlet_values(v1xv2)
display(x)


end # module

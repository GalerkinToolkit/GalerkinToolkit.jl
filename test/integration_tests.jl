module IntegrationTests

using Test
import GalerkinToolkit as gk
using GalerkinToolkit: ∫
using LinearAlgebra
import ForwardDiff

spx0 = gk.unit_simplex(0)
spx1 = gk.unit_simplex(1)
spx2 = gk.unit_simplex(2)
spx3 = gk.unit_simplex(3)

cube0 = gk.unit_n_cube(0)
cube1 = gk.unit_n_cube(1)
cube2 = gk.unit_n_cube(2)
cube3 = gk.unit_n_cube(3)

degree = 4
quad = gk.default_quadrature(spx0,degree)
quad = gk.default_quadrature(spx1,degree)
quad = gk.default_quadrature(spx2,degree)
quad = gk.default_quadrature(spx3,degree)

quad = gk.default_quadrature(cube0,degree)
@test sum(gk.weights(quad)) ≈ 1
quad = gk.default_quadrature(cube1,degree)
@test sum(gk.weights(quad)) ≈ 1
quad = gk.default_quadrature(cube2,degree)
@test sum(gk.weights(quad)) ≈ 1
quad = gk.default_quadrature(cube3,degree)
@test sum(gk.weights(quad)) ≈ 1

outdir = mkpath(joinpath(@__DIR__,"..","output"))

domain = (0,2,0,2)
cells = (8,8)
mesh = gk.cartesian_mesh(domain,cells)

Ω = gk.domain(mesh)
Ωref = gk.domain(mesh;is_reference_domain=true)
ϕ = gk.domain_map(Ωref,Ω)
u = gk.analytical_field(x->sum(x),Ω)

degree = 2
dΩref = gk.measure(Ωref,degree)
int = ∫(dΩref) do q
    x = ϕ(q)
    J = gk.call(ForwardDiff.jacobian,ϕ,q)
    dV = gk.call(J->abs(det(J)),J)
    gk.call(*,u(x),dV)
end

@test sum(int) ≈ 8

int = ∫(dΩref) do q
    J = gk.call(ForwardDiff.jacobian,ϕ,q)
    dV = gk.call(J->abs(det(J)),J)
    dV
end

@test sum(int) ≈ 4

end # module

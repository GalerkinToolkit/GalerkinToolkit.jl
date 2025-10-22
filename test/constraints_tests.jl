module ConstraintsTests

import GalerkinToolkit as GT
using LinearAlgebra
using Test
using ForwardDiff
import LinearSolve

const ∇ = ForwardDiff.gradient

function laplace_solve(V)
    Ω = GT.domain(V)
    order = GT.order(V)
    dΩ = GT.quadrature(Ω,2*order)
    a = (u,v)->GT.∫(x->∇(u,x)⋅∇(v,x),dΩ)
    l = 0
    uhd = GT.zero_dirichlet_field(Float64,V)
    p = GT.SciMLBase_LinearProblem(uhd,a,l)
    s = LinearSolve.solve(p)
    uh = GT.solution_field(uhd,s)
    int = GT.∫(x->uh(x),dΩ)
    sum(int)
end

domain = (0,1,0,1)
cells = (2,2)
mesh = GT.cartesian_mesh(domain,cells;periodic=(true,false))
Ω = GT.interior(mesh)
Γ = GT.boundary(mesh;group_names=["1-face-2"])

order = 1
V = GT.lagrange_space(Ω,order;dirichlet_boundary=Γ)
@test isa(V,GT.ConstrainedSpace)
C = GT.constraints(V)
@test isa(C,GT.PeriodicConstraints)

display(C)

face_C = GT.face_constraints(V)

face = 4
C = face_C[face]
display(C)

m,n = size(C)
a = zeros(m)
b = rand(n)
mul!(a,C,b)
@test a == b

a .= 1
b .= 0
mul!(b,transpose(C),a)
@test a == b

a2 = copy(a)
b2 = copy(b)
b2 .= 1
a2 .= 0
mul!(a2,C,b2)
@test a2 == b

laplace_solve(V)

V = GT.lagrange_space(Ω,order;dirichlet_boundary=Γ,periodic_scaling=0.5)

@test isa(V,GT.ConstrainedSpace)
C = GT.constraints(V)
@test isa(C,GT.PeriodicConstraints)

face_C = GT.face_constraints(V)

face = 4
C = face_C[face]
display(C)

m,n = size(C)
a = zeros(m)
b = rand(n)
mul!(a,C,b)

a .= 1
b .= 0
mul!(b,transpose(C),a)
@show b

a2 = copy(a)
b2 = copy(b)
b2 .= 1
a2 .= 0
mul!(a2,C,b2)
@test a2 == b

laplace_solve(V)

end # module

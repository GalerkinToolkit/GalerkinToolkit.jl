module ConstraintsTests

import GalerkinToolkit as GT
using LinearAlgebra
using Test

domain = (0,1,0,1)
cells = (2,2)
mesh = GT.cartesian_mesh(domain,cells;periodic=(true,false))
Ω = GT.interior(mesh)
Γ = GT.boundary(mesh;group_names=["1-face-2"])

order = 1
V0 = GT.lagrange_space(Ω,order;dirichlet_boundary=Γ)

C = GT.periodic_constraints(V0)

display(C)

V = GT.constrain(V0,C)

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

@show GT.periodic_dofs(V0)
C = GT.periodic_constraints(V0;scaling=0.5)


display(C)

V = GT.constrain(V0,C)

display(GT.face_dofs(V0))
display(GT.face_dofs(V))

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



end # module

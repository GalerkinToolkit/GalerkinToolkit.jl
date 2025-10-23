module ConstraintsTests

import GalerkinToolkit as GT
using LinearAlgebra
using Test
using ForwardDiff
import LinearSolve
using StaticArrays


# How to define these periodic conditions
# in a unit square mesh:
# top = 0.5*bottom
# right = 2*left

# Create mesh of the unit square
# We ask to allocate space for periodic info
domain = (0,1,0,1)
cells = (2,2)
periodic = Val(true)
mesh = GT.cartesian_mesh(domain,cells;periodic)

# Define the master boundary
# as the union of the bottom and left face
bottom_face = "1-face-1"
left_face = "1-face-3"
group_names = [bottom_face,left_face]
Γ_master = GT.boundary(mesh;group_names)

# Define a map from master to worker boundary
piecewise = Val(true)
coupling = GT.analytical_field(Γ_master;piecewise) do x,name
    if name === bottom_face
        SVector(x[1],1.0) # bottom -> top
    elseif name === left_face
        SVector(1.0,x[2]) # left -> right
    end
end

# Find couplings and
# update the periodic info in the mesh
tol = 1e-8
GT.find_periodic_nodes!(coupling;tol)

# NB. The particular computations above can also be done in one line as follows,
# but the syntax above is more general
# mesh = GT.cartesian_mesh(domain,cells;periodic=(true,true))

# Define the scaling
# top = 0.5*bottom
# right = 2*left
periodic_scaling = GT.analytical_field(Γ_master;piecewise) do x,name
    if name === bottom_face
        0.5
    elseif name === left_face
        2.0
    end
end

# Define interpolation space
Ω = GT.interior(mesh)
order = 1
V = GT.lagrange_space(Ω,order;periodic,periodic_scaling)
C = GT.constraints(V)

display(C)


domain = (0,1,0,1)
cells = (2,2)
mesh = GT.cartesian_mesh(domain,cells;periodic=(true,true))
display(collect(enumerate(GT.periodic_nodes(mesh))))
Ω = GT.interior(mesh)
order = 2
V = GT.lagrange_space(Ω,order;periodic=true)

domain = (0,1,0,1)
cells = (2,2)
mesh = GT.cartesian_mesh(domain,cells;periodic=Val(true))
Γ_master = GT.boundary(mesh;group_names=["1-face-1"])
u = GT.analytical_field(x->SVector(1-x[1],1.0),Γ_master)
GT.find_periodic_nodes!(u)

Ω = GT.interior(mesh)

order = 2
g = GT.analytical_field(x->2*x[1]-1,Γ_master)
V = GT.lagrange_space(Ω,order;periodic=true,periodic_scaling=g)

C = GT.constraints(V)

display(C)


function laplace_solve(V)
    ∇ = ForwardDiff.gradient
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


msh =  joinpath(@__DIR__,"..","assets","periodic.msh")
mesh = GT.mesh_from_msh(msh;periodic=Val(true))
@show GT.periodic_nodes(mesh)

order = 1
Ω = GT.interior(mesh)
V = GT.lagrange_space(Ω,order;periodic=true)
@test isa(V,GT.ConstrainedSpace)
C = GT.constraints(V)
@test isa(C,GT.PeriodicConstraints)

face_C = GT.face_constraints(V)

face = 4
C = face_C[face]
display(C)

laplace_solve(V)


domain = (0,1,0,1)
cells = (2,2)
mesh = GT.cartesian_mesh(domain,cells;periodic=(true,false))
Ω = GT.interior(mesh)
Γ = GT.boundary(mesh;group_names=["1-face-2"])

order = 1
V = GT.lagrange_space(Ω,order;dirichlet_boundary=Γ,periodic=true)
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

V = GT.lagrange_space(Ω,order;dirichlet_boundary=Γ,periodic=true,periodic_scaling=0.5)

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

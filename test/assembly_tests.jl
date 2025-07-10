module AssemblyTests

import GalerkinToolkit as GT
import PartitionedSolvers as PS
using GalerkinToolkit: ∫, ×
using Test
import ForwardDiff
using LinearAlgebra
using AbstractTrees
import WriteVTK

outdir = mkpath(joinpath(@__DIR__,"..","output"))

Ω = GT.unit_n_cube(Val(2))
mesh = GT.mesh(Ω)
Γdiri = GT.boundary(mesh)
@test GT.num_faces(Ω) == 1

order = 2
V = GT.lagrange_space(Ω,order;dirichlet_boundary=Γdiri)
GT.free_dofs(V)
V² = V × V
GT.allocate_vector(Float64,V²,Ω)
GT.allocate_matrix(Float64,V²,V²,Ω)

alloc = GT.allocate_vector(Float64,V,Ω;free_or_dirichlet=GT.FREE)
GT.contribute!(alloc,[1,2],[-1,1])
b = GT.compress(alloc)
@test length(b) == GT.num_free_dofs(V)

alloc = GT.allocate_vector(Float64,V,Ω;free_or_dirichlet=GT.DIRICHLET)
GT.contribute!(alloc,[1,2],[-1,1])
b = GT.compress(alloc)
@test length(b) == GT.num_dirichlet_dofs(V)

alloc = GT.allocate_matrix(Float64,V,V,Ω;free_or_dirichlet=(GT.FREE,GT.DIRICHLET))
GT.contribute!(alloc,zeros(2,2),[-1,1],[1,-1])
A = GT.compress(alloc)
@test size(A) == (GT.num_free_dofs(V),GT.num_dirichlet_dofs(V))

GT.num_free_dofs(V)
GT.num_dirichlet_dofs(V)
@test GT.num_free_dofs(V) == 1
@test GT.num_dirichlet_dofs(V) == 8

u = GT.analytical_field(sum,Ω)
uhd = GT.zero_dirichlet_field(Float64,V)
GT.interpolate_dirichlet!(u,uhd)
uh = GT.zero_field(Float64,V)
GT.interpolate_free!(u,uh)

# This returns the dofs on 1-faces, but in this context it is expected to return
# the dofs on faces on field component 1
GT.face_dofs(V,1)

#plt = GT.plot(Ω,refinement=20)
#GT.plot!(plt,u;label="u")
#GT.plot!(plt,uhd;label="uhd")
#GT.plot!(plt,uh;label="uh")
#WriteVTK.vtk_grid("check",plt) |> WriteVTK.close

display(GT.face_dofs(V))

degree = 2*order
dΩ = GT.measure(Ω,degree)

∇(u,q) = ForwardDiff.gradient(u,q)

a(u,v) = ∫( q->∇(u,q)⋅∇(v,q), dΩ)
l(v) = 0

p = GT.PartitionedSolvers_linear_problem(uhd,a,l)

A = PS.matrix(p)
display(A)
b = PS.rhs(p)
x = A\b
uh = GT.solution_field(uhd,x)

eh(q) = u(q) - uh(q)
∇eh(q) = ∇(u,q) - ∇(uh,q)

tol = 1.e-10

el2 = ∫( q->abs2(eh(q)), dΩ) |> sum |> sqrt
@test el2 < tol

eh1 = ∫( q->∇eh(q)⋅∇eh(q), dΩ) |> sum |> sqrt
@test el2 < tol

domain = (0,1,0,1)
cells = (4,4)
mesh = GT.cartesian_mesh(domain,cells)
GT.group_interior_faces!(mesh;physical_name="interior_faces")
GT.group_boundary_faces!(mesh;physical_name="boundary_faces")
Ω = GT.interior(mesh)
Γ = GT.boundary(mesh;group_names=["boundary_faces"])
Λ = GT.skeleton(mesh;group_names=["interior_faces"])

order = 2
degree = 2*order
dΩ = GT.measure(Ω,degree)
dΓ = GT.measure(Γ,degree)
dΛ = GT.measure(Λ,degree)

V = GT.lagrange_space(Ω,order;dirichlet_boundary=Γ)

face_point_x = GT.coordinate_accessor(dΩ)
face_point_J = GT.jacobian_accessor(dΩ)
face_point_dV = GT.weight_accessor(dΩ)
face_point_dof_s = GT.shape_function_accessor(GT.value,V,dΩ)
face_point_dof_∇s = GT.shape_function_accessor(ForwardDiff.gradient,V,dΩ)

face = 3
point = 4
dof = 2

point_x = face_point_x(face)
point_J = face_point_J(face)
point_dV = face_point_dV(face)
point_dof_s = face_point_dof_s(face)
point_dof_∇s = face_point_dof_∇s(face)

x = point_x(point)
J = point_J(point)
point_dV(point,J)

dof_s = point_dof_s(point)
dof_∇s = point_dof_∇s(point,J)

s = dof_s(dof)
∇s = dof_∇s(dof)
@show ∇s

T = Float64
uh = GT.zero_field(T,V)
ufun = x->0.0
u = GT.analytical_field(ufun,Ω)
GT.interpolate_dirichlet!(u,uh)
#face_dirichlet! = GT.dirichlet_accessor(uh,Ω)
#dirichlet! = face_dirichlet!(face)
#
face_dofs = GT.dofs_accessor(V,Ω)
#
n = maximum(map(GT.num_dofs,GT.reference_spaces(V)))
Auu = zeros(T,n,n)
bu = zeros(T,n)
#dirichlet!(Auu,bu)
#@show bu

dofs = face_dofs(face)

b_alloc = GT.allocate_vector(T,V,Ω)
GT.contribute!(b_alloc,bu,dofs)
b = GT.compress(b_alloc)

A_alloc = GT.allocate_matrix(T,V,V,Ω)
GT.contribute!(A_alloc,Auu,dofs,dofs)
A = GT.compress(A_alloc)

uhd = GT.zero_dirichlet_field(T,V)
GT.interpolate_dirichlet!(u,uhd)
f = x -> 0

face_npoints = GT.num_points_accessor(dΩ)

b_alloc = GT.allocate_vector(T,V,Ω)
A_alloc = GT.allocate_matrix(T,V,V,Ω)
for face in 1:GT.num_faces(Ω)
    point_x = face_point_x(face)
    point_J = face_point_J(face)
    point_dV = face_point_dV(face)
    point_dof_s = face_point_dof_s(face)
    point_dof_∇s = face_point_dof_∇s(face)
    npoints = face_npoints(face)
    dofs = face_dofs(face)
    #dirichlet! = face_dirichlet!(face)
    fill!(Auu,zero(eltype(Auu)))
    fill!(bu,zero(eltype(bu)))
    for point in 1:npoints
        x = point_x(point)
        J = point_J(point)
        dV = point_dV(point,J)
        dof_s = point_dof_s(point)
        dof_∇s = point_dof_∇s(point,J)
        for (i,dofi) in enumerate(dofs)
            v = dof_s(i)
            ∇v = dof_∇s(i)
            bu[i] += f(x)*v*dV
            for (j,dofj) in enumerate(dofs)
                ∇u = dof_∇s(j)
                Auu[i,j] += ∇v⋅∇u*dV
            end
        end
    end
    #dirichlet!(Auu,bu)
    GT.contribute!(b_alloc,bu,dofs)
    GT.contribute!(A_alloc,Auu,dofs,dofs)
end

b = GT.compress(b_alloc)
A = GT.compress(A_alloc)
x = A\b

uh = GT.solution_field(uhd,x)

face_point_val = GT.discrete_field_accessor(GT.value,uh,dΩ)

err = Ref(0.0)
for face in 1:GT.num_faces(Ω)
    npoints = face_npoints(face)
    point_x = face_point_x(face)
    point_J = face_point_J(face)
    point_dV = face_point_dV(face)
    point_val = face_point_val(face)
    for point in 1:npoints
        x = point_x(point)
        J = point_J(point)
        dV = point_dV(point,J)
        val = point_val(point,J)
        err[] += abs2(val-ufun(x))*dV
    end
end
tol = 1e-12
@test sqrt(err[]) < tol

face_point_val = GT.update(face_point_val;discrete_field=uh)

err = Ref(0.0)
for face in 1:GT.num_faces(Ω)
    npoints = face_npoints(face)
    point_x = face_point_x(face)
    point_J = face_point_J(face)
    point_dV = face_point_dV(face)
    point_val = face_point_val(face)
    for point in 1:npoints
        x = point_x(point)
        J = point_J(point)
        dV = point_dV(point,J)
        val = point_val(point,J)
        err[] += abs2(val-ufun(x))*dV
    end
end
tol = 1e-12
@test sqrt(err[]) < tol


# Poisson solve (API in physical domain)

domain = (0,1,0,1)
cells = (4,4)
mesh = GT.cartesian_mesh(domain,cells)
D = GT.num_dims(mesh)
Ω = GT.interior(mesh)
Λ = GT.skeleton(mesh)
Γdiri = GT.boundary(mesh)
order = 1
degree = 2*order
V = GT.lagrange_space(Ω,order;dirichlet_boundary=Γdiri)
dΩ = GT.measure(Ω,degree)
dΛ = GT.measure(Λ,degree)
uhd = GT.zero_dirichlet_field(Float64,V)
u = GT.analytical_field(sum,Ω)
GT.interpolate_dirichlet!(u,uhd)

n_Λ = GT.unit_normal(mesh,GT.num_dims(mesh)-1)
∇(u,q) = ForwardDiff.gradient(u,q)

jump(v,p) = v[2](p) - v[1](p)
jump(u,n,x) = u[2](x)*n[2](x) + u[1](x)*n[1](x)
mean(f,u,x) = 0.5*(f(u[1],x)+f(u[2],x))
γ = GT.uniform_quantity(1.0)
h_Λ = GT.face_diameter_field(Λ)

# The skeleton terms are not needed. They are added just
# to make sure that they are computed correctly.
a(u,v) = ∫( q->∇(u,q)⋅∇(v,q), dΩ) + ∫(x->(γ/h_Λ(x))*jump(v,n_Λ,x)⋅jump(u,n_Λ,x)-jump(v,n_Λ,x)⋅mean(∇,u,x)-mean(∇,v,x)⋅jump(u,n_Λ,x),dΛ)
l(v) = 0

p = GT.PartitionedSolvers_linear_problem(uhd,a,l)

x = PS.matrix(p)\PS.rhs(p)
uh = GT.solution_field(uhd,x)

eh(q) = u(q) - uh(q)
∇eh(q) = ∇(u,q) - ∇(uh,q)

tol = 1.0e-10
el2 = ∫( q->abs2(eh(q)), dΩ) |> sum |> sqrt
@test el2 < tol

eh1 = ∫( q->∇eh(q)⋅∇eh(q), dΩ) |> sum |> sqrt
@test eh1 < tol



domain = (0,1,0,1)
cells = (4,4)
mesh = GT.cartesian_mesh(domain,cells)
GT.group_interior_faces!(mesh;physical_name="interior_faces")
GT.group_boundary_faces!(mesh;physical_name="boundary_faces")
D = GT.num_dims(mesh)

Ω = GT.interior(mesh)
Ωref = GT.interior(mesh;is_reference_domain=true)
ϕ = GT.physical_map(mesh,D)

D = GT.num_dims(mesh)
Γdiri = GT.boundary(mesh;group_names=["1-face-1","1-face-3"])

Γref = GT.boundary(mesh;
                 is_reference_domain=true,
                 group_names=["1-face-2","1-face-4"])

Γ = GT.physical_domain(Γref)
Λ = GT.skeleton(mesh)

V = GT.lagrange_space(Ω,1;dirichlet_boundary=Γdiri)
@show GT.num_free_dofs(V)
@show GT.num_dirichlet_dofs(V)

degree = 2
dΩ = GT.measure(Ω,degree)
dΓ = GT.measure(Γ,degree)
dΛ = GT.measure(Λ,degree)



function a(u,v)
    ∫(dΩ) do q
        u(q)*v(q)
    end +
    ∫(dΓ) do q
        u(q)*v(q)
    end +
    ∫(dΛ) do p
        jump(u,p)*jump(v,p)
    end
end

A = GT.assemble_matrix(a,Float64,V,V)
@test size(A,1) == GT.num_free_dofs(V)

V² = V × V

function a((u1,u2),(v1,v2))
    ∫(dΩ) do q
        (v1(q)+v2(q))*(u1(q)+u2(q))
    end +
    ∫(dΛ) do q
        jump(u1,q)*jump(v2,q)
    end
end

A = GT.assemble_matrix(a,Float64,V²,V²)

function l(v)
    ∫(dΛ) do p
        jump(v,p)
    end
end

b = GT.assemble_vector(l,Float64,V)
@test sum(b)+1 ≈ 1

function l(v)
    ∫(dΩ) do q
        v(q)
    end +
    ∫(dΓ) do p
        v(p)
    end +
    ∫(dΛ) do p
        jump(v,p)
    end
end

b = GT.assemble_vector(l,Float64,V)

@test length(b) == GT.num_free_dofs(V)

function l((v1,v2))
    ∫(dΩ) do q
        v1(q)*v2(q)
    end +
    ∫(dΛ) do p
        jump(v1,p)*jump(v2,p)
    end
end

b = GT.assemble_vector(l,Float64,V²)


#V = GT.lagrange_space(Ωref,1;dirichlet_boundary=Γdiri)
#@show GT.num_free_dofs(V)
#@show GT.num_dirichlet_dofs(V)
#
#degree = 2
#dΩref = GT.measure(Ωref,degree)
#ϕ = GT.physical_map(mesh,D)
#
#dΓref = GT.measure(Γref,degree)
#α = GT.physical_map(mesh,D-1)
#β = GT.reference_map(mesh,D-1,D)
#
#Λref = GT.skeleton(mesh;
#                 is_reference_domain=true,
#                 group_names=["interior_faces"])
#
#Λ = GT.physical_domain(Λref)
#dΛref = GT.measure(Λref,degree)
#ϕ_Λref_Λ = GT.physical_map(mesh,D-1)
#ϕ_Λref_Ωref = GT.reference_map(mesh,D-1,D)
#
#function dV(J)
#    abs(det(J))
#end
#
#function dS(J)
#    Jt = transpose(J)
#    sqrt(det(Jt*J))
#end
#
#jump(u,ϕ,q) = u(ϕ(q)[2])-u(ϕ(q)[1])
#
#@show GT.is_boundary(Γref)
#
#function l(v)
#    ∫(dΩref) do q
#        J = ForwardDiff.jacobian(ϕ,q)
#        v(q)*dV(J)
#    end +
#    ∫(dΓref) do p
#        q = β(p)
#        J = ForwardDiff.jacobian(α,p)
#        v(q)*dS(J)
#    end +
#    ∫(dΛref) do p
#        J = ForwardDiff.jacobian(ϕ_Λref_Λ,p)
#        jump(v,ϕ_Λref_Ωref,p)*dS(J)
#    end
#end
#
#b = GT.assemble_vector(l,Float64,V)
#
#@test length(b) == GT.num_free_dofs(V)
#
#
#function l(v)
#    ∫(dΛref) do p
#        J = ForwardDiff.jacobian(ϕ_Λref_Λ,p)
#        jump(v,ϕ_Λref_Ωref,p)*dS(J)
#    end
#end
#
#b = GT.assemble_vector(l,Float64,V)
#@test sum(b)+1 ≈ 1
#
#V² = V × V
#
#function l((v1,v2))
#    ∫(dΩref) do q
#        J = ForwardDiff.jacobian(ϕ,q)
#        v1(q)*v2(q)*dV(J)
#    end +
#    ∫(dΛref) do p
#        J = ForwardDiff.jacobian(ϕ_Λref_Λ,p)
#        jump(v1,ϕ_Λref_Ωref,p)*jump(v2,ϕ_Λref_Ωref,p)*dS(J)
#    end
#end
#
#b = GT.assemble_vector(l,Float64,V²)
#
#function a(u,v)
#    ∫(dΩref) do q
#        J = ForwardDiff.jacobian(ϕ,q)
#        u(q)*v(q)*dV(J)
#    end +
#    ∫(dΓref) do p
#        q = β(p)
#        J = ForwardDiff.jacobian(α,p)
#        u(q)*v(q)*dS(J)
#    end +
#    ∫(dΛref) do p
#        J = ForwardDiff.jacobian(ϕ_Λref_Λ,p)
#        index = GT.generate_index(Λ)
#        jump(v,ϕ_Λref_Ωref,p)*jump(u,ϕ_Λref_Ωref,p)*dS(J)
#    end
#end
#
#A = GT.assemble_matrix(a,Float64,V,V)
#@test size(A,1) == GT.num_free_dofs(V)
#
#function a((u1,u2),(v1,v2))
#    ∫(dΩref) do q
#        J = ForwardDiff.jacobian(ϕ,q)
#        v1(q)*v2(q)*u1(q)*u2(q)*dV(J)
#    end +
#    ∫(dΛref) do p
#        J = ForwardDiff.jacobian(ϕ_Λref_Λ,p)
#        jump(v1,ϕ_Λref_Ωref,p)*jump(v2,ϕ_Λref_Ωref,p)*jump(u1,ϕ_Λref_Ωref,p)*jump(u2,ϕ_Λref_Ωref,p)*dS(J)
#    end
#end
#
#A = GT.assemble_matrix(a,Float64,V²,V²)
#
#x = similar(b,axes(A,2))
#fill!(x,0)
#uh = GT.solution_field(V²,x)
#uh1,uh2 = uh
#fill!(GT.free_values(uh2),1)
#@test x[end] == 1
#
#function dV(ϕ,q)
#    J = ForwardDiff.jacobian(ϕ,q)
#    abs(det(J))
#end
#
#a(u,v) = ∫( q->u(q)*v(q)*dV(ϕ,q), dΩref)
#
#f = GT.analytical_field(sum,Ω)
#
#l(v) = ∫( q->f(ϕ(q))*v(q)*dV(ϕ,q), dΩref)
#
#V = GT.lagrange_space(Ωref,1)#;dirichlet_boundary=Γdiri)
#
#p = GT.linear_problem(Float64,V,a,l)
#A = PS.matrix(p)
#b = PS.rhs(p)
#display(A)
#@show GT.num_free_dofs(V)
#@test size(A,1) == GT.num_free_dofs(V)
#@test length(b) == GT.num_free_dofs(V)
#x = PS.matrix(p)\PS.rhs(p)
#uh = GT.solution_field(V,x)
#
#function ∇(u,phi,q)
#   J = ForwardDiff.jacobian(phi,q)
#   g = ForwardDiff.gradient(u,q)
#   J\g
#end
#
#a(u,v) = ∫( q->∇(u,ϕ,q)⋅∇(v,ϕ,q)*dV(ϕ,q), dΩref)
#l(v) = 0
#
#p = GT.linear_problem(uh,a,l)
#display(PS.matrix(p))
#x = PS.matrix(p)\PS.rhs(p)
#
## Poisson solve (reference domain)
#
#domain = (0,1,0,1)
#cells = (4,4)
#mesh = GT.cartesian_mesh(domain,cells)
#GT.group_boundary_faces!(mesh;physical_name="boundary_faces")
#
#Ω = GT.interior(mesh)
#Ωref = GT.reference_domain(Ω)
#D = GT.num_dims(mesh)
#ϕ = GT.physical_map(mesh,D)
#
#Γdiri = GT.boundary(mesh;group_names=["boundary_faces"])
#
##V = GT.iso_parametric_space(Ωref;dirichlet_boundary=Γdiri)
#
#order = 3
#V = GT.lagrange_space(Ωref,order;dirichlet_boundary=Γdiri)
#
#u = GT.analytical_field(sum,Ω)
#uhd = GT.zero_dirichlet_field(Float64,V)
## TODO
##GT.interpolate_dirichlet!(q->u(ϕ(q)),uh)
#GT.interpolate_dirichlet!(u∘ϕ,uhd)
#
#function ∇(u,q)
#   J = ForwardDiff.jacobian(ϕ,q)
#   g = ForwardDiff.gradient(u,q)
#   J\g
#end
#
#function dV(q)
#    J = ForwardDiff.jacobian(ϕ,q)
#    abs(det(J))
#end
#
#degree = 2*order
#dΩref = GT.measure(Ωref,degree)
#
#a(u,v) = ∫( q->∇(u,q)⋅∇(v,q)*dV(q), dΩref)
#l(v) = 0
#
#p = GT.linear_problem(uhd,a,l)
#x = PS.matrix(p)\PS.rhs(p)
#uh = GT.solution_field(uhd,x)
#
## TODO
## Functions like this ones should
## work as AbstractQuantities?
#eh(q) = u(ϕ(q)) - uh(q)
#∇eh(q) = ForwardDiff.gradient(u,ϕ(q)) - ∇(uh,q)
#
#tol = 1.e-12
#el2 = ∫( q->abs2(eh(q))*dV(q), dΩref) |> sum |> sqrt
#@test el2 < tol
#
#eh1 = ∫( q->∇eh(q)⋅∇eh(q)*dV(q), dΩref) |> sum |> sqrt
#@test el2 < tol


# Nonlinear case

pl::Int = 3
flux(∇u) = norm(∇u)^(pl-2) * ∇u
dflux(∇du,∇u) = (pl-2)*norm(∇u)^(pl-4)*(∇u⋅∇du)*∇u+norm(∇u)^(pl-2)*∇du

f = GT.analytical_field(Ω) do x
    sum(x)
end

res(u) = v -> ∫( x-> ∇(v,x)⋅GT.call(flux,∇(u,x)) - f(x)*v(x) , dΩ)
jac(u) = (du,v) -> ∫( x-> ∇(v,x)⋅GT.call(dflux,∇(du,x),∇(u,x)) , dΩ)

uh = GT.rand_field(Float64,V)
p = GT.PartitionedSolvers_nonlinear_problem(uh,res,jac)

linsolve = PS.NLsolve_nlsolve_linsolve(PS.LinearAlgebra_lu,p)
s = PS.NLsolve_nlsolve(p;show_trace=true,method=:newton)
s = PS.solve(s)
uh = GT.solution_field(uh,s)

# 3d case

n = 2
domain = (0,1,0,1,0,1)
cells = (n,n,n)
mesh = GT.cartesian_mesh(domain,cells)
Ω = GT.interior(mesh)
k = 1
V = GT.lagrange_space(Ω,k)
dΩ = GT.measure(Ω,2*k)
∇3 = ForwardDiff.gradient
a(u,v) = GT.∫( x->∇3(u,x)⋅∇3(v,x), dΩ)
l(v) = 0
p = GT.PartitionedSolvers_linear_problem(Float64,V,a,l)
PS.matrix(p) |> display

#tini = 0.0
#tend = 10.0
#domain = (0,1,0,1)
#cells = (20,20)
#mesh = GT.cartesian_mesh(domain,cells)
#GT.group_boundary_faces!(mesh;physical_name="boundary_faces")
#Ω = GT.interior(mesh)
#Γ = GT.boundary(mesh;group_names=["boundary_faces"])
#g = x -> t -> exp(-2*t)*sinpi(t*x[1])*(x[2]^2-1)
#g0(t) = GT.analytical_field(x->g(x)(t),Ω)
#g1(t) = GT.analytical_field(x->ForwardDiff.derivative(g(x),t),Ω)
#α2(t) = GT.analytical_field(x->1+sin(t)*(x[1]^2+x[2]^2)/4, Ω)
#f2(t) = GT.analytical_field(x->sin(t)*sinpi(x[1])*sinpi(x[2]), Ω)
#order = 1
#dΩ = GT.measure(Ω,2*order)
#∇2 = ForwardDiff.gradient
#m(t,dtu,v) = ∫( x -> v(x)*dtu(x), dΩ)
#a(t,u,v) = ∫( x-> α2(t)(x)*∇2(v,x)⋅∇2(u,x), dΩ)
#l(t,v) = ∫( x->v(x)*f2(t)(x),dΩ)
#res(t,(u0,u1)) = v -> m(t,u1,v) + a(t,u0,v) - l(t,v)
#jac0(t,(u0,u1)) = (du,v) -> a(t,du,v)
#jac1(t,(u0,u1)) = (du,v) -> m(t,du,v)
#V = GT.lagrange_space(Ω,order;dirichlet_boundary=Γ)
#ut0 = GT.semi_discrete_field(Float64,V) do t,uh
#    GT.interpolate_dirichlet!(g0(t),uh)
#end
#ut1 = GT.semi_discrete_field(Float64,V) do t,uh
#    GT.interpolate_dirichlet!(g1(t),uh)
#end
#GT.interpolate_free!(g0(tini),ut0(tini))
##TODO
##this one is linear actually
#p = GT.nonlinear_ode((ut0,ut1),tini:tend,res,(jac0,jac1))
#x0 = GT.free_values(ut0(tini))
#x1 = GT.free_values(ut1(tini))
#PS.update(p,solution=(tini,x0,x1))



end # module

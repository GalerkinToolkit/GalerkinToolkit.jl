module HelpersTests

import GalerkinToolkit as GT
import ForwardDiff
using LinearAlgebra
using Test

domain = (0,1,0,1)
cells = (4,4)
mesh = GT.cartesian_mesh(domain,cells)
GT.label_interior_faces!(mesh;physical_name="interior_faces")
GT.label_boundary_faces!(mesh;physical_name="boundary_faces")
Ω = GT.interior(mesh)
Γ = GT.boundary(mesh;physical_names=["boundary_faces"])
Λ = GT.skeleton(mesh;physical_names=["interior_faces"])

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
dV = point_dV(point,J)

dof_s = point_dof_s(point)
dof_∇s = point_dof_∇s(point,J)

s = dof_s(dof)
∇s = dof_∇s(dof)
@show ∇s

T = Float64
uh = GT.zero_field(T,V)
ufun = sum
u = GT.analytical_field(ufun,Ω)
GT.interpolate_dirichlet!(u,uh)
face_dirichlet! = GT.dirichlet_accessor(uh,Ω)
dirichlet! = face_dirichlet!(face)

face_dofs = GT.dofs_accessor(V,Ω)

n = maximum(map(GT.num_dofs,GT.reference_fes(V)))
Auu = zeros(T,n,n)
bu = zeros(T,n)
dirichlet!(Auu,bu)
@show bu

dofs = face_dofs(face)

b_alloc = GT.allocate_vector(T,V,Ω)
GT.contribute!(b_alloc,bu,dofs)
b = GT.compress(b_alloc)

A_alloc = GT.allocate_matrix(T,V,V,Ω)
GT.contribute!(A_alloc,Auu,dofs,dofs)
A = GT.compress(A_alloc)

uhd = GT.dirichlet_field(T,V)
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
    dirichlet! = face_dirichlet!(face)
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
    dirichlet!(Auu,bu)
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

#display(A)


end # module

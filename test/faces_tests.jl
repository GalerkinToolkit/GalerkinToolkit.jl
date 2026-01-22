module FacesTests

using Test
import GalerkinToolkit as GT
using ForwardDiff
#using InteractiveUtils

domain = (0,1,0,1)
n = 4
cells = (n,n)
mesh = GT.cartesian_mesh(domain,cells)

@test length(GT.each_face_new(mesh,Val(1))) == GT.num_faces(mesh,Val(1))

for F in GT.each_face_new(mesh,Val(1))
    for A in GT.each_face_around_new(F,Val(2))
        a = GT.local_face(A,F)
        A2 = GT.parent_face(a)
        @test GT.id(A) == GT.id(A2)
        @test GT.num_dims(A) == GT.num_dims(A2)
        @test GT.mesh(A) === GT.mesh(A2)
        F2 = GT.global_face(a)
        @test GT.id(F) == GT.id(F2)
        @test GT.num_dims(F) == GT.num_dims(F2)
        @test GT.mesh(F) === GT.mesh(F2)
        P = GT.node_permutation(A,a)
        @assert GT.nodes(F) == GT.nodes(A)[ GT.nodes(a)[P]]
        P = GT.node_permutation(a)
        @assert GT.nodes(F) == GT.nodes(A)[ GT.nodes(a)[P]]
    end
end

for A in GT.each_face_new(mesh,Val(2))
    for a in GT.each_local_face(A,Val(1))
        A2 = GT.parent_face(a)
        @test GT.id(A) == GT.id(A2)
        @test GT.num_dims(A) == GT.num_dims(A2)
        @test GT.mesh(A) === GT.mesh(A2)
    end
end

Ω = GT.interior(mesh)
Γ = GT.boundary(mesh)
Λ = GT.skeleton(mesh)

k = 1
V = GT.lagrange_space(Ω,k)

dΩ = GT.measure(Ω,2*k)
dΓ = GT.measure(Γ,2*k)
dΛ = GT.measure(Λ,2*k)

D = GT.num_dims(mesh)
mesh_Dfaces = GT.each_face_new(mesh,D,dΩ)
mesh_Dfaces = GT.tabulate(GT.value,mesh_Dfaces)
mesh_Dfaces = GT.tabulate(ForwardDiff.gradient,mesh_Dfaces)


mesh_faces = GT.each_face_new(dΩ)
mesh_face = mesh_faces[2]
@show GT.nodes(mesh_face)
@show GT.node_coordinates(mesh_face)
@show GT.barycenter(mesh_face)
@show GT.diameter(mesh_face)
mesh_points = GT.each_point_new(mesh_face)
mesh_point = mesh_points[1]
@show GT.coordinate(GT.value,mesh_point)
@show GT.coordinate(ForwardDiff.jacobian,mesh_point)
@show GT.weight(mesh_point)

xx

end # module

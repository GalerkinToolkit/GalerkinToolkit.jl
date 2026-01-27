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
    @test GT.num_faces_around(F,Val(2)) in (1,2)
    for A in GT.each_face_around_new(F,Val(2))
        a = GT.local_face(A,F)
        P = GT.node_permutation(A,a)
        @assert GT.nodes(F) == GT.nodes(A)[ GT.nodes(a)[P]]
    end
end

for A in GT.each_face_new(mesh,Val(2))
    @test GT.num_local_faces(A,Val(1)) == 4
    for a in GT.each_local_face(A,Val(1))
        F = GT.global_face(A,a)
        P = GT.node_permutation(A,a)
        @assert GT.nodes(F) == GT.nodes(A)[ GT.nodes(a)[P]]
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

for mesh_face in GT.each_face_new(dΩ)
    GT.nodes(mesh_face)
    GT.node_coordinates(mesh_face)
    GT.num_nodes(mesh_face)
    for mesh_point in GT.each_point_new(mesh_face)
        GT.coordinate(mesh_point)
        GT.coordinate(GT.value,mesh_point)
        GT.coordinate(ForwardDiff.jacobian,mesh_point)
        GT.weight(mesh_point)
    end
end

mesh_dfaces = GT.each_face_new(mesh,D,dΛ)
mesh_dface = mesh_dfaces[2]
mesh_Dfaces = GT.each_face_around_new(mesh_dface)
mesh_Dface = mesh_Dfaces[2]
@show GT.nodes(mesh_Dface)
@show GT.node_coordinates(mesh_Dface)
mesh_points = GT.each_point_new(mesh_Dface)
mesh_point = mesh_points[1]
@show GT.coordinate(GT.value,mesh_point)
@show GT.coordinate(ForwardDiff.jacobian,mesh_point)
@show GT.weight(mesh_point)
@show GT.unit_normal(mesh_point)
@show GT.coordinate(mesh_point)
@show GT.jacobian(mesh_point)

for mesh_dface in mesh_dfaces
    for mesh_Dface in GT.each_face_around_new(mesh_dface)
        GT.nodes(mesh_Dface)
        GT.node_coordinates(mesh_Dface)
        GT.num_nodes(mesh_Dface)
        for mesh_point in GT.each_point_new(mesh_Dface)
            GT.coordinate(mesh_point)
            GT.coordinate(GT.value,mesh_point)
            GT.coordinate(ForwardDiff.jacobian,mesh_point)
            GT.weight(mesh_point)
            GT.unit_normal(mesh_point)
        end
    end
end

mesh_dfaces = GT.each_face_new(mesh,D,dΓ)
mesh_dface = mesh_dfaces[2]
mesh_Dface = mesh_dface
@show GT.nodes(mesh_Dface)
@show GT.node_coordinates(mesh_Dface)
mesh_points = GT.each_point_new(mesh_Dface)
mesh_point = mesh_points[1]
@show GT.coordinate(GT.value,mesh_point)
@show GT.coordinate(ForwardDiff.jacobian,mesh_point)
@show GT.weight(mesh_point)
@show GT.unit_normal(mesh_point)
@show GT.coordinate(mesh_point)
@show GT.jacobian(mesh_point)

for mesh_Dface in mesh_Dfaces
    GT.nodes(mesh_Dface)
    GT.node_coordinates(mesh_Dface)
    GT.num_nodes(mesh_Dface)
    for mesh_point in GT.each_point_new(mesh_Dface)
        GT.coordinate(mesh_point)
        GT.coordinate(GT.value,mesh_point)
        GT.coordinate(ForwardDiff.jacobian,mesh_point)
        GT.weight(mesh_point)
        GT.unit_normal(mesh_point)
    end
end


end # module

module Issue224

import GalerkinToolkit as GT
import PartitionedSolvers as PS
import StaticArrays
using LinearAlgebra
import ForwardDiff
using Test

#Node coordinates
S = StaticArrays.SVector{3,Float64}
node_coordinates = S[(0,0,0),(1,0,0),(0,1,0),(0,0,1)]

#Face nodes
face_nodes_0 = [[3],[2],[1],[4]]
face_nodes_1 = [[2,1],[3,1],[3,2],[4,1],[4,2],[4,3]]
face_nodes_2 = [[3,2,1],[4,2,1],[2,4,3],[4,1,3]]
face_nodes_3 = [[1,2,3,4]]
face_nodes = [
    face_nodes_0,
    face_nodes_1,
    face_nodes_2,
    face_nodes_3]

#Reference spaces
vertex = GT.unit_simplex(Val(0))
segment = GT.unit_simplex(Val(1))
triangle = GT.unit_simplex(Val(2))
tet = GT.unit_simplex(Val(3))
order = 1
vertex1 = GT.lagrange_space(vertex,order)
segment2 = GT.lagrange_space(segment,order)
triangle3 = GT.lagrange_space(triangle,order)
tet4 = GT.lagrange_space(tet,order)
reference_spaces_0 = (vertex1,)
reference_spaces_1 = (segment2,)
reference_spaces_2 = (triangle3,)
reference_spaces_3 = (tet4,)
reference_spaces = (
    reference_spaces_0,
    reference_spaces_1,
    reference_spaces_2,
    reference_spaces_3)

#Create mesh
mesh = GT.create_mesh(;
    node_coordinates,
    face_nodes,
    reference_spaces,
    is_cell_complex=Val(true),
   )


D = GT.num_dims(mesh)


k = 1
degree = 2 * k
Ω = GT.interior(mesh)
Γ = GT.boundary(mesh)
dΩ = GT.measure(Ω,degree)
dΓ = GT.measure(Γ,degree)
D = GT.num_dims(mesh)
n = GT.unit_normal(mesh,D-1)
u = GT.analytical_field(sum,Ω)
T = Float64
∇ = ForwardDiff.gradient

V = GT.lagrange_space(Ω,k)
a = (u,v) -> begin
    GT.∫( x -> v(x)*u(x)-v(x)*n(x)⋅∇(u,x)-n(x)⋅∇(v,x)*u(x), dΓ) +
    GT.∫( x -> ∇(u,x)⋅∇(v,x), dΩ)
end
l = v -> begin
    GT.∫( x -> v(x)*u(x)-n(x)⋅∇(v,x)*u(x), dΓ) +
    GT.∫( x -> v(x)*0, dΩ)
end
p = GT.PartitionedSolvers_linear_problem(T,V,a,l)
s = PS.solve(p)
uh = GT.solution_field(V,s)
int = GT.∫( x -> (u(x)-uh(x))^2, dΩ)
tol = 1.0e-10
@test sqrt(sum(int)) < tol

end # module

module Issue230

import Gmsh
import GalerkinToolkit as GT
import PartitionedSolvers as PS
using LinearAlgebra
import ForwardDiff
import StaticArrays

mesh = GT.with_gmsh() do gmsh
    lc = 1.0
    gmsh.model.geo.addPoint(0, 0, 0, lc, 1)
    gmsh.model.geo.addPoint(1, 0,  0, lc, 2)
    gmsh.model.geo.addPoint(1, 1, 0, lc, 3)
    gmsh.model.geo.addPoint(0, 1, 0, lc, 4)
    gmsh.model.geo.addLine(1, 2, 1)
    gmsh.model.geo.addLine(2, 3, 2)
    gmsh.model.geo.addLine(3, 4, 3)
    gmsh.model.geo.addLine(4, 1, 4)
    gmsh.model.geo.addCurveLoop([4, 1, 2, 3], 1)
    gmsh.model.geo.addPlaneSurface([1], 1)
    gmsh.model.geo.mesh.setAlgorithm(2, 1, 8)
    gmsh.model.geo.mesh.setRecombine(2, 1)
    gmsh.model.geo.extrude([(2, 1)], 0, 0, 1, [1.0/lc], [], true)
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(3)
    GT.mesh_from_gmsh(gmsh;complexify=true)
end

#import Makie
#import GLMakie
#shrink = 0.6
#fig = GT.makie_surfaces(mesh;dim=3,shrink,axis=(;aspect=:data))
#GT.makie_surfaces!(mesh;dim=2,shrink)
#GT.makie_edges!(mesh;dim=1,shrink)
#GT.makie_vertices!(mesh;dim=0)
#Makie.save("issue_230.png",fig)
#display(fig)

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
@assert sqrt(sum(int)) < 1.0e-9

end #module

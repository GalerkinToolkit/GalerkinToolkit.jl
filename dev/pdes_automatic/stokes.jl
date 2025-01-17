
# # Stokes
#
# !!! warning
#     TODO
#     * Problem statement
#     * Allow g1 and g1 to be defined on the boundary
#     * Build a pressure field with zero mean
#

import GalerkinToolkit as GT
import PartitionedSolvers as PS
import GLMakie as Makie
import ForwardDiff
import StaticArrays
import Tensors
using LinearAlgebra
import FileIO # hide

domain = (0,1,0,1)
cells = (20,20)
D = length(cells)
mesh = GT.cartesian_mesh(domain,cells)
Ω = GT.interior(mesh)
Γ1 = GT.boundary(mesh;physical_names=["1-face-2"])
Γ2 = GT.boundary(mesh;physical_names=["1-face-1","1-face-3","1-face-4"])
#TODO
#g1 = GT.analytical_field(x->SVector(1,0),Γ1)
#g2 = GT.analytical_field(x->SVector(0,0),Γ2)
#g = GT.piecewiese_field(g1,g2)
#Γ = GT.domain(g)
g1 = GT.analytical_field(x->StaticArrays.SVector(1,0),Ω)
g2 = GT.analytical_field(x->StaticArrays.SVector(0,0),Ω)
g = GT.piecewiese_field(g1,g2)
Γ = GT.piecewiese_domain(Γ1,Γ2)
order = 2
V = GT.lagrange_space(Ω,order;space_type=:Q,tensor_size=Val((D,)),dirichlet_boundary=Γ)
Q = GT.lagrange_space(Ω,order-1;space_type=:P,dirichlet_boundary=GT.last_dof())
VxQ = V × Q
u_field, p_field = 1,2
uhph_dirichlet = GT.dirichlet_field(Float64,VxQ)
GT.interpolate_dirichlet!(g,uhph_dirichlet,u_field)
dΩ = GT.measure(Ω,2*order)
∇ = ForwardDiff.jacobian
div(u,x) = tr(∇(u,x))
a((u,p),(v,q)) = GT.∫( x-> ∇(v,x)⋅∇(u,x) - div(v,x)*p(x) + q(x)*div(u,x), dΩ)
l((v,q)) = 0
p = GT.linear_problem(uhph_dirichlet,a,l)
s = PS.LinearAlgebra_lu(p)
s = PS.solve(s)
#TODO
#uh = GT.solution_field(uhph_dirichlet,s,u_field)
#ph = GT.solution_field(uhph_dirichlet,s,p_field;zeromean=true)
uh,ph = GT.solution_field(uhph_dirichlet,s)
Makie.plot(Ω,color=ph)
Makie.arrows!(uh;color=x->norm(uh(x)),lengthscale=0.1)
FileIO.save(joinpath(@__DIR__,"fig_pt_stokes.png"),Makie.current_figure()) # hide

# ![](fig_pt_stokes.png)


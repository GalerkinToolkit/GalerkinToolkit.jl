
# # Stokes lid-driven cavity
#
# ![](fig_stokes_1.png)
#
# ## Problem statement
#
# !!! warning
#     TODO
#
# ## Numerical scheme
#
# !!! warning
#     TODO
#
# ## Implementation
#

import GalerkinToolkit as GT
import LinearSolve
import GLMakie as Makie
import ForwardDiff
import StaticArrays
import Tensors
using LinearAlgebra
import FileIO # hide

#Mesh
domain = (0,1,0,1)
cells = (20,20)
D = length(cells)
mesh = GT.cartesian_mesh(domain,cells)

#Domains
Ω = GT.interior(mesh)
Γ = GT.boundary(mesh;group_names=["1-face-1","1-face-2","1-face-3","1-face-4"])
g = GT.analytical_field(Γ;piecewise=true) do x,name
    if name == "1-face-2"
        StaticArrays.SVector(1,0)
    else
        StaticArrays.SVector(0,0)
    end
end

#Velocity interpolation
order = 2
tensor_size=Val((D,))
dirichlet_boundary = Γ
V = GT.lagrange_space(Ω,order;tensor_size,dirichlet_boundary)

#Pressure interpolation
space_type=:P
dirichlet_boundary=GT.last_dof()
Q = GT.lagrange_space(Ω,order-1;space_type,dirichlet_boundary)

#Cartesian product space
VxQ = V × Q

#Dirichlet condition
uhph_dirichlet = GT.zero_dirichlet_field(Float64,VxQ)
uhd, = uhph_dirichlet
GT.interpolate_dirichlet!(g,uhd)

#Weak form
dΩ = GT.measure(Ω,2*order)
∇ = ForwardDiff.jacobian
div(u,x) = tr(∇(u,x))
a((u,p),(v,q)) = GT.∫( x-> ∇(v,x)⋅∇(u,x) - div(v,x)*p(x) + q(x)*div(u,x), dΩ)
l((v,q)) = 0

#Assemble linear problem and solve it
p = GT.SciMLBase_LinearProblem(uhph_dirichlet,a,l)
sol = LinearSolve.solve(p)

#Get solution fields
uh,ph = GT.solution_field(uhph_dirichlet,sol)

#Visualize solution
fig = Makie.Figure()
aspect = Makie.DataAspect()
ax = Makie.Axis(fig[1,1];aspect)
Makie.hidespines!(ax)
Makie.hidedecorations!(ax)
shading = Makie.NoShading
GT.makie_surfaces!(Ω;color=ph,shading)
GT.makie_arrows2d!(Ω,uh;color=x->norm(uh(x)),lengthscale=0.1)
FileIO.save(joinpath(@__DIR__,"fig_stokes_1.png"),Makie.current_figure()) # hide



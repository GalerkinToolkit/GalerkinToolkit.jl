# # Transient heat equation
# 
# ![](fig_transient_heat_equation_1.gif)
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

using LinearAlgebra
import GalerkinToolkit as GT
import DifferentialEquations
import ForwardDiff
import GLMakie as Makie
import FileIO # hide

#Parameters
mesh_size=0.02
R=0.15
T=2
N=100

#Generate mesh with GMSH Julia API
mesh = GT.with_gmsh() do gmsh
    R = 0.15
    dim = 2
    gmsh.option.set_number("General.Verbosity", 2)
    rect_tag = gmsh.model.occ.add_rectangle(0,0,0,1,1)
    circle_tag = gmsh.model.occ.add_circle(0.5,0.5,0,R)
    circle_curve_tag = gmsh.model.occ.add_curve_loop([circle_tag])
    circle_surf_tag = gmsh.model.occ.add_plane_surface([circle_curve_tag])
    gmsh.model.occ.cut([(dim,rect_tag)],[(dim,circle_surf_tag)]);
    gmsh.model.occ.synchronize()
    domain_tags = [1]
    outer_tags = [6,7,8,9]
    inner_tags = [5]
    gmsh.model.model.add_physical_group(dim,domain_tags,-1,"domain")
    gmsh.model.model.add_physical_group(dim-1,outer_tags,-1,"outer")
    gmsh.model.model.add_physical_group(dim-1,inner_tags,-1,"inner")
    gmsh.option.set_number("Mesh.MeshSizeMax",mesh_size)
    gmsh.model.mesh.generate(dim)
    GT.mesh_from_gmsh(gmsh)
end

#Domains
Ω = GT.interior(mesh;physical_names=["domain"])
Γ1 = GT.boundary(mesh;physical_names=["outer"])
Γ2 = GT.boundary(mesh;physical_names=["inner"])
Γ = GT.piecewise_domain(Γ1,Γ2)

#Interpolation
k = 1
V = GT.lagrange_space(Ω,k;dirichlet_boundary=Γ)
uh = GT.undef_field(Float64,V)

#Integration
degree = 2*GT.order(V)
dΩ = GT.measure(Ω,degree)

#Initial condition
u0 = GT.analytical_field(x->0.0,Ω)
GT.interpolate_free!(u0,uh)

#Time-dependent Dirichlet function
α = t -> sin(3*pi*t)
function dirichlet_dynamics!(t,uh,duh=nothing)
    g1 = GT.analytical_field(x->0.0,Ω)
    if uh !== nothing
        g2 = GT.analytical_field(x->α(t),Ω)
        g = GT.piecewise_field(g1,g2)
        GT.interpolate_dirichlet!(g,uh)
    end
    if duh !== nothing
        g2 = GT.analytical_field(x->ForwardDiff.derivative(α,t),Ω)
        g = GT.piecewise_field(g1,g2)
        GT.interpolate_dirichlet!(g,duh)
    end
end

#Definition of the ODE problem
C = 10
∇ = ForwardDiff.gradient
m = (u,v) -> GT.∫(x->C*v(x)*u(x),dΩ)
a = (u,v) -> -1*GT.∫(x->∇(u,x)⋅∇(v,x), dΩ)
r = (uh,t) -> v -> a(uh,v)
j = (uh,t) -> a
tspan = (0.0,T)
problem = GT.SciMLBase_ODEProblem(tspan,uh,m,r,j;dirichlet_dynamics!)

#Selection and setup of the ODE solver
dt = T/N
solver = DifferentialEquations.QNDF(autodiff=false);
initializealg=DifferentialEquations.NoInit()
adaptive=false
save_on=false
integrator = DifferentialEquations.init(
    problem,solver;initializealg,dt,adaptive,save_on)

#Setup Makie scene
axis = (aspect = Makie.DataAspect(),)
fig = Makie.Figure()
colorrange = (-1,1)
ax,sc = GT.makie_surface(fig[1,1],Ω;color=uh,axis,colorrange)
Makie.hidespines!(ax)
Makie.hidedecorations!(ax)
Makie.Colorbar(fig[1,2],sc)

#Record Makie scene while solving
fn = "fig_transient_heat_equation_1.gif"
file = joinpath(@__DIR__,fn)
Makie.record(fig,file,integrator) do integrator
    sc.color = GT.solution_field(integrator)
end
nothing # hide


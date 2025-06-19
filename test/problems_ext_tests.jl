module ProblemsExt

import GalerkinToolkit as GT
import LinearSolve
import NonlinearSolve
import DifferentialEquations
import ForwardDiff
import WriteVTK
using Test
using LinearAlgebra

# ODE problem

mesh_size=0.02
R=0.15
T=2.0
N=100

mesh = GT.with_gmsh() do gmsh
    R = 0.15
    dim = 2
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
    gmsh.option.setNumber("Mesh.MeshSizeMax",mesh_size)
    gmsh.model.mesh.generate(dim)
    GT.mesh_from_gmsh(gmsh)
end


Ω = GT.interior(mesh;physical_names=["domain"])
Γ1 = GT.boundary(mesh;physical_names=["outer"])
Γ2 = GT.boundary(mesh;physical_names=["inner"])
Γ = GT.piecewise_domain(Γ1,Γ2)
k = 1
V = GT.lagrange_space(Ω,k;dirichlet_boundary=Γ)

α = t -> sin(3*pi*t)

function dirichlet_dynamics!(t,uh,vh=nothing)
    g1 = GT.analytical_field(x->0.0,Ω)
    if uh !== nothing
        g2 = GT.analytical_field(x->α(t),Ω)
        g = GT.piecewise_field(g1,g2)
        GT.interpolate_dirichlet!(g,uh)
    end
    if vh !== nothing
        g2 = GT.analytical_field(x->ForwardDiff.derivative(α,t),Ω)
        g = GT.piecewise_field(g1,g2)
        GT.interpolate_dirichlet!(g,vh)
    end
end

uh = GT.zero_field(Float64,V)

u0 = GT.analytical_field(x->0.0,Ω)

degree = 2*GT.val_parameter(GT.order(V))
dΩ = GT.measure(Ω,degree)

C = 10
∇ = (u,q) -> ForwardDiff.gradient(u,q)
m = (u,v) -> GT.∫(x->C*v(x)*u(x),dΩ)
a = (u,v) -> -1*GT.∫(x->∇(u,x)⋅∇(v,x), dΩ)
r = (uh,t) -> v -> a(uh,v)
j = (uh,t) -> a
tspan = (0.0,T)
dt = T/N
GT.interpolate_free!(u0,uh)
@time prob = GT.SciMLBase_ODEProblem(tspan,uh,m,r,j;dirichlet_dynamics!)
@time timestepper = DifferentialEquations.Rodas5P(autodiff=false);
# TODO why init is so slow?
@time integrator = DifferentialEquations.init(
    prob, timestepper; initializealg=DifferentialEquations.NoInit(),
    dt,adaptive=false)

plt = GT.plot(Ω)
file = "ode"
WriteVTK.paraview_collection(file,plt) do pvd
    # TODO is this really lazy?
    # or are we storing all intermediate states inside integrator?
    for  (step, (x,t)) in enumerate(DifferentialEquations.intervals(integrator))
        @show (step,t)
        WriteVTK.vtk_grid("$(file)_$(step)",plt) do plt
            dirichlet_dynamics!(t,uh)
            GT.plot!(plt,GT.solution_field(uh,x);label="uh")
            pvd[step] = plt
        end
    end
end

# TODO We are solving a linear ode as a nonlinear one.
# Is it possible to solve linear ODEs as linear ODEs with DifferentialEquations?

## TODO I would like something like
#WriteVTK.paraview_collection(file,plt) do pvd
#    for  (step, i) in enumerate(integrator)
#        @show (step,t)
#        uh_at_t = GT.solution_field(i) #TODO
#        duh_at_t = GT.solution_field(i,Val(1)) #TODO
#        WriteVTK.vtk_grid("$(file)_$(step)",plt) do plt
#            GT.plot!(plt,uh_at_t;label="uh")
#            GT.plot!(plt,duh_at_t;label="duh")
#            pvd[step] = plt
#        end
#    end
#end

# Linear problem

n = 10
domain = (0,1,0,1)
cells = (n,n)
mesh = GT.cartesian_mesh(domain,cells)
D = GT.num_dims(mesh)
Ω = GT.interior(mesh)
Γ = GT.boundary(mesh)
order = 1
degree = 2*order
V = GT.lagrange_space(Ω,order;dirichlet_boundary=Γ)
dΩ = GT.measure(Ω,degree)
u = GT.analytical_field(sum,Ω)
uhd = GT.interpolate_dirichlet(u,V)
∇ = (u,q) -> ForwardDiff.gradient(u,q)
a = (u,v) -> GT.∫( q->∇(u,q)⋅∇(v,q), dΩ)
l = v -> 0

prob =  GT.SciMLBase_LinearProblem(uhd,a,l)
sol = LinearSolve.solve(prob)
uh = GT.solution_field(uhd,sol)

eh = q -> u(q) - uh(q)
∇eh =q -> ∇(u,q) - ∇(uh,q)
tol = 1.0e-10
el2 = GT.∫( q->abs2(eh(q)), dΩ) |> sum |> sqrt
eh1 = GT.∫( q->∇eh(q)⋅∇eh(q), dΩ) |> sum |> sqrt
@test el2 < tol
@test eh1 < tol

# Nonlinear problem

n = 10
domain = (0,1,0,1)
cells = (n,n)
mesh = GT.cartesian_mesh(domain,cells)
D = GT.num_dims(mesh)
Ω = GT.interior(mesh)
order = 1
degree = 2*order
V = GT.lagrange_space(Ω,order;dirichlet_boundary=Γ)
dΩ = GT.measure(Ω,degree)
uh = GT.rand_field(Float64,V)

∇ = (u,q) -> ForwardDiff.gradient(u,q)
const q = 3
flux(∇u) = norm(∇u)^(q-2) * ∇u
dflux(∇du,∇u) = (q-2)*norm(∇u)^(q-4)*(∇u⋅∇du)*∇u+norm(∇u)^(q-2)*∇du
res = u -> v -> GT.∫( x-> ∇(v,x)⋅GT.call(flux,∇(u,x)) - v(x), dΩ)
jac = u -> (du,v) -> GT.∫( x-> ∇(v,x)⋅GT.call(dflux,∇(du,x),∇(u,x)) , dΩ)

prob = GT.SciMLBase_NonlinearProblem(uh,res,jac)
sol = NonlinearSolve.solve(prob;show_trace=Val(true))
@test sol.retcode == NonlinearSolve.ReturnCode.Success
uh = GT.solution_field(uh,sol)
uhl2 = GT.∫( q->abs2(uh(q)), dΩ) |> sum |> sqrt
@test abs(uhl2-0.09133166701839236)<tol

end # module

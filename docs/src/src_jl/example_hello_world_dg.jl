# # Discontinuous Galerkin
#
# ![](fig_hello_world_dg_1.png)
#
# ## Problem statement
#
# We solve the same PDE as in the [Hello, World!](@ref) example, but this time using a discontinuous Galerkin scheme.
#
#  ## Numerical scheme
#
# We consider the symmetric interior penalty method [Arnold2002](@cite).
#

# ## Implementation
#  We solve the problem and visualize the solution. In this case, we draw the average of solution field
#  on the interior $2$-faces of the computational mesh. There faces are where the interior penalty is enforced.
#

# ## Implementation

import FileIO # hide
using LinearAlgebra
import GalerkinToolkit as GT
import PartitionedSolvers as PS
import ForwardDiff
import GLMakie as Makie
import LinearSolve

#Geometry
domain = (0,1,0,1,0,1)
cells = (4,4,4)
mesh = GT.cartesian_mesh(domain,cells)
D = GT.num_dims(mesh)
n = GT.unit_normal(mesh,D-1)
Ω = GT.interior(mesh)
Γd = GT.boundary(mesh)
Λ = GT.skeleton(mesh)
h_Λ = GT.face_diameter_field(Λ)
h_Γd = GT.face_diameter_field(Γd)

#Functions
const ∇ = ForwardDiff.gradient
g = GT.analytical_field(sum,Ω)
f = GT.analytical_field(x->0,Ω)
mean(f,u,x) = 0.5*(f(u[1],x)+f(u[2],x))
jump(u,n,x) = u[2](x)*n[2](x) + u[1](x)*n[1](x)

#Interpolation
k = 1
V = GT.lagrange_space(Ω,k;continuous=false)

#Integration
degree = 2*k
dΩ = GT.quadrature(Ω,degree)
dΛ = GT.quadrature(Λ,degree)
dΓd = GT.quadrature(Γd,degree)

#Weak form
γ = GT.uniform_quantity(k*(k+1)/10)
a = (u,v) -> begin
    #Laplace operator
    GT.∫(dΩ) do x
        ∇(u,x)⋅∇(v,x)
    end +
    #Interior penalty
    GT.∫(dΛ) do x
        (γ/h_Λ(x))*jump(v,n,x)⋅jump(u,n,x)-
        jump(v,n,x)⋅mean(∇,u,x)-
        mean(∇,v,x)⋅jump(u,n,x)
    end +
    #Nitsche term
    GT.∫(dΓd) do x
        (γ/h_Γd(x))*v(x)*u(x)-
        v(x)*n(x)⋅∇(u,x)-
        n(x)⋅∇(v,x)*u(x)
    end
end
l = v -> begin
    #RHS
    GT.∫(dΩ) do x
         v(x)*f(x)
    end +
    #Nietche term
    GT.∫(dΓd) do x
        (γ/h_Γd(x))*v(x)*g(x)-
        n(x)⋅∇(v,x)*g(x)
    end
end

#Linear problem
p = GT.SciMLBase_LinearProblem(Float64,V,a,l)
sol = LinearSolve.solve(p)
uh = GT.solution_field(V,sol)

#Error check
eh = x -> uh(x) - g(x)
el2 = GT.∫( x->abs2(eh(x)), dΩ) |> sum |> sqrt
@assert el2 < 1.0e-9

#Visualization
fig = Makie.Figure()
ax = Makie.Axis3(fig[1,1],aspect=:data)
Makie.hidespines!(ax)
Makie.hidedecorations!(ax)
GT.makie_surfaces!(Λ;color=x->mean(GT.value,uh,x))
FileIO.save(joinpath(@__DIR__,"fig_hello_world_dg_1.png"),Makie.current_figure()) # hide
nothing # hide


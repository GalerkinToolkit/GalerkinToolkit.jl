# # Laplace-Beltrami
#
# ![](fig_laplace_beltrami_1.png)

# ## Problem statement
#
#

# ## Implementation

import FileIO # hide
import GalerkinToolkit as GT
import GLMakie as Makie
import ForwardDiff
using LinearAlgebra
import LinearSolve


cells = (4,40)
mesh = GT.moebius_strip(cells;width=0.6)
Ω = GT.interior(mesh)
Γ = GT.boundary(mesh)
V = GT.lagrange_space(Ω,2;dirichlet_boundary=Γ)
∇ = ForwardDiff.gradient
f = GT.analytical_field(x->1,Ω)
dΩ = GT.measure(Ω,4)

a = (u,v) -> begin
    GT.∫(dΩ) do x
        ∇(u,x)⋅∇(v,x)
    end
end
l = (v) -> begin
    GT.∫(dΩ) do x
        v(x)*f(x)
    end
end

uhd = GT.zero_field(Float64,V)
p = GT.SciMLBase_LinearProblem(uhd,a,l)
sol = LinearSolve.solve(p)
uh = GT.solution_field(uhd,sol)


fig = Makie.Figure()
elevation = 0.24π
azimuth = -0.55π
aspect = :data
ax = Makie.Axis3(fig[1,1];aspect,elevation,azimuth)
Makie.hidespines!(ax)
Makie.hidedecorations!(ax)
shading=Makie.NoShading
GT.makie_surfaces!(Ω;color=uh,shading,refinement=2)
FileIO.save(joinpath(@__DIR__,"fig_laplace_beltrami_1.png"),Makie.current_figure()) # hide
nothing # hide


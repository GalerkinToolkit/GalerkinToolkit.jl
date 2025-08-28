# # Hello, World!
#
# ![](fig_hello_world_1.png)
#
#
# ## Problem statement
#
# In this example, we show how to solve the "Hello, world" PDE example:
# the Laplace equation on the unit hyper-cube $\Omega  =[0,1]^d$, $d=3$, with Dirichlet boundary conditions.
#
# ```math
# \left\lbrace
# \begin{aligned}
# -\Delta u = f  \ &\text{in} \ \Omega,\\
# u = g \ &\text{on}\ \partial\Omega,\\
# \end{aligned}
# \right.
# ```
# with $f=0$ and $g(x)=\text{sum}(x)$. In this case, we know that the solution is $u=g$ which allows us to check that we solve the
# problem correctly, by integration an error norm.

#  ## Numerical scheme
#
#  We use a conventional Galerkin finite element (FE) method with conforming Lagrangian FE spaces (see, e.g., [Johnson2009](@cite)).  The weak form equation we solve consists in finding $u_h\in V_g$ such that $a(u_h,v) = \ell(v)$ for all $v\in V_0$. To this end we build a space $V$ spanned by continuous and piece-wise Lagrangian basis functions. The auxiliary spaces $V_g$ and $V_0$ are the subsets of $V$ that fulfill the Dirichlet boundary condition $g$ and $0$ on $\partial\Omega$ respectively. The bilinear and linear forms are
# ```math
#   a(u,v) \doteq \int_{\Omega} \nabla v \cdot \nabla u \ {\rm d}\Omega, \quad b(v) \doteq \int_{\Omega} v\ f  \ {\rm  d}\Omega.
# ```
#
# This equation results in a system of linear algebraic equations that is solved using an external linear solver from [`LinearSolve.jl`](https://github.com/SciML/LinearSolve.jl).

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
cells = (10,10,10)
mesh = GT.cartesian_mesh(domain,cells)
Ω = GT.interior(mesh)
Γd = GT.boundary(mesh)

#Functions
∇ = ForwardDiff.gradient
g = GT.analytical_field(sum,Ω)
f = GT.analytical_field(x->0,Ω)

#Interpolation
k = 1
V = GT.lagrange_space(Ω,k;dirichlet_boundary=Γd)
T = Float64
uhd = GT.zero_dirichlet_field(T,V)
GT.interpolate_dirichlet!(g,uhd)

#Weak form
degree = 2*k
dΩ = GT.measure(Ω,degree)
a = (u,v) -> GT.∫( x->∇(u,x)⋅∇(v,x), dΩ)
l = v -> GT.∫( x->v(x)*f(x), dΩ)

#Linear problem
p = GT.SciMLBase_LinearProblem(uhd,a,l)
sol = LinearSolve.solve(p)
uh = GT.solution_field(uhd,sol)

#Error check
eh = x -> uh(x) - g(x)
el2 = GT.∫( x->abs2(eh(x)), dΩ) |> sum |> sqrt
@assert el2 < 1.0e-9

#Visualization
fig = Makie.Figure()
ax = Makie.Axis3(fig[1,1],aspect=:data)
Makie.hidespines!(ax)
Makie.hidedecorations!(ax)
GT.makie_surfaces!(Ω;color=uh)
FileIO.save(joinpath(@__DIR__,"fig_hello_world_1.png"),Makie.current_figure()) # hide
nothing # hide

# ## Explicit integration loops
#
# This other code version implements the integration loops manually instead of relying on the underlying automatic code generation.


#Manually written matrix assembly function
#Always use a function, never the global scope
function assemble_matrix!(A_alloc,Ad_alloc,V,dΩ)

    ∇ = ForwardDiff.gradient

    #Iterators to the quantities on the
    #integration points
    V_faces_0 = GT.foreach_face(V,dΩ)
    V_faces = GT.tabulate(∇,V_faces_0)

    #Temporaries
    n = GT.max_num_reference_dofs(V)
    T = Float64
    Auu = zeros(T,n,n)

    #Numerical integration loop
    for V_face in V_faces

        #Get quantities at current face
        dofs = GT.dofs(V_face)

        #Reset face matrix
        fill!(Auu,zero(T))

        #Loop over integration points
        for V_point in GT.foreach_point(V_face)

            #Get quantities at current integration point
            dV = GT.weight(V_point)
            dof_∇s = GT.shape_functions(∇,V_point)

            #Fill in face matrix
            for (i,dofi) in enumerate(dofs)
                ∇v = dof_∇s[i]
                for (j,dofj) in enumerate(dofs)
                    ∇u = dof_∇s[j]
                    Auu[i,j] += ∇v⋅∇u*dV
                end
            end
        end

        #Add face contribution to the
        #global allocations
        GT.contribute!(A_alloc,Auu,dofs,dofs)
        GT.contribute!(Ad_alloc,Auu,dofs,dofs)
    end
end

#Manually written integration of the error
#Always use a function, never the global scope
function integrate_l2_error(g,uh,dΩ)

    #Iterators to the quantities on the
    #integration points
    uh_faces_1 = GT.foreach_face(uh,dΩ)
    uh_faces_2 = GT.tabulate(GT.value,uh_faces_1)
    uh_faces = GT.compute(GT.coordinate,uh_faces_2)

    #Numerical integration loop
    s = 0.0
    for uh_face in uh_faces
        for uh_point in GT.foreach_point(uh_face)

            #Get quantities at current integration point
            x = GT.coordinate(uh_point)
            dV = GT.weight(uh_point)
            uhx = GT.field(GT.value,uh_point)

            #Add contribution
            s += abs2(uhx-g.definition(x))*dV
        end
    end

    #Compute the l2 norm
    el2 = sqrt(s)
end

#Allocate matrix for free columns
A_alloc = GT.allocate_matrix(T,V,V,Ω)

#Allocate matrix for dirichlet columns
free_or_dirichlet=(GT.FREE,GT.DIRICHLET)
Ad_alloc = GT.allocate_matrix(T,V,V,Ω;free_or_dirichlet)

#Fill allocations with the function we wrote above
assemble_matrix!(A_alloc,Ad_alloc,V,dΩ)

#Compress matrix into the final format
A = GT.compress(A_alloc)
Ad = GT.compress(Ad_alloc)

#Build the linear system
xd = GT.dirichlet_values(uhd)
b = -Ad*xd
p = LinearSolve.LinearProblem(A,b)

#Solve the problem
sol = LinearSolve.solve(p)
uh = GT.solution_field(uhd,sol)

#Integrate the error l2 norm
#with the function we wrote above
el2 = integrate_l2_error(g,uh,dΩ)
@assert el2 < 1.0e-9
nothing # hide


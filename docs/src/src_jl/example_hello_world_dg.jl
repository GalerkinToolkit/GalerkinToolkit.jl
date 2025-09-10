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
const g = GT.analytical_field(sum,Ω)
const f = GT.analytical_field(x->0,Ω)
mean(f,u,x) = 0.5*(f(u[1],x)+f(u[2],x))
jump(u,n,x) = u[2](x)*n[2](x) + u[1](x)*n[1](x)

#Interpolation
k = 1
const γ0 = k*(k+1)/10
γ = GT.uniform_quantity(γ0)
V = GT.lagrange_space(Ω,k;continuous=false)

#Integration
degree = 2*k
dΩ = GT.quadrature(Ω,degree)
dΛ = GT.quadrature(Λ,degree)
dΓd = GT.quadrature(Γd,degree)

#Weak form
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


# ## Explicit integration loops
#
# This other code version implements the integration loops manually instead of relying on the underlying automatic code generation.

function assemble_Ω!(A_alloc,b_alloc,V,dΩ)

    #Temporaries
    n = GT.max_num_reference_dofs(V)
    T = Float64
    Auu = zeros(T,n,n)
    bu = zeros(T,n)

    #Loop over bulk faces
    tabulate = (∇,GT.value)
    compute=(GT.coordinate,)
    for V_face in GT.each_face(V,dΩ;tabulate,compute)
        dofs = GT.dofs(V_face)
        ndofs = length(dofs)
        fill!(Auu,zero(T))
        fill!(bu,zero(T))
        for V_point in GT.each_point(V_face)
            dV = GT.weight(V_point)
            dof_∇s = GT.shape_functions(∇,V_point)
            dof_s = GT.shape_functions(GT.value,V_point)
            x = GT.coordinate(V_point)
            fx = f.definition(x)
            for i in 1:ndofs
                v = dof_s[i]
                bu[i] += v*fx*dV
            end
            for j in 1:ndofs
                ∇u = dof_∇s[j]
                for i in 1:ndofs
                    ∇v = dof_∇s[i]
                    Auu[i,j] += ∇v⋅∇u*dV
                end
            end
        end
        GT.contribute!(A_alloc,Auu,dofs,dofs)
        GT.contribute!(b_alloc,bu,dofs)
    end

end

function assemble_Λ!(A_alloc,V,dΛ)

    #Temporaries
    n = GT.max_num_reference_dofs(V)
    T = Float64
    Auu = zeros(T,n,n)

    #Loop over skeleton faces
    tabulate = (∇,GT.value)
    compute=(GT.unit_normal,)
    V_dfaces = GT.each_face(V,dΛ;tabulate,compute)
    dΛ_dfaces = GT.each_face(dΛ)
    for (dΛ_dface,V_dface) in zip(dΛ_dfaces,V_dfaces)
        h = GT.diameter(dΛ_dface)
        dΛ_points = GT.each_point(dΛ_dface)
        for V_Dface_i in GT.each_face_around(V_dface)
            dofs_i = GT.dofs(V_Dface_i)
            ndofs_i = length(dofs_i)
            V_points_i = GT.each_point(V_Dface_i)
            for V_Dface_j in GT.each_face_around(V_dface)
                V_points_j = GT.each_point(V_Dface_j)
                dofs_j = GT.dofs(V_Dface_j)
                ndofs_j = length(dofs_j)
                fill!(Auu,zero(T))
                for (dΛ_point,V_point_i,V_point_j) in zip(dΛ_points,V_points_i,V_points_j)
                    dS = GT.weight(dΛ_point)
                    s_i = GT.shape_functions(GT.value,V_point_i)
                    s_j = GT.shape_functions(GT.value,V_point_j)
                    ∇s_i = GT.shape_functions(∇,V_point_i)
                    ∇s_j = GT.shape_functions(∇,V_point_j)
                    n_i = GT.unit_normal(V_point_i)
                    n_j = GT.unit_normal(V_point_j)
                    for j in 1:ndofs_j
                        for i in 1:ndofs_i
                            jump_jump = (γ0/h)*(s_i[i]*n_i)⋅(s_j[j]*n_j)
                            jump_mean = (s_i[i]*n_i)⋅(0.5*∇s_j[j])
                            mean_jump = (0.5*∇s_i[i])⋅(s_j[j]*n_j)
                            Auu[i,j] += (jump_jump - jump_mean - mean_jump)*dS
                        end
                    end
                end
                GT.contribute!(A_alloc,Auu,dofs_i,dofs_j)
            end
        end
    end

end

function assemble_Γd!(A_alloc,b_alloc,V,dΓd)

    #Temporaries
    n = GT.max_num_reference_dofs(V)
    T = Float64
    Auu = zeros(T,n,n)
    bu = zeros(T,n)

    #Loop over Dirichlet boundary
    dΓd_dfaces = GT.each_face(dΓd)
    tabulate = (∇,GT.value)
    compute=(GT.unit_normal,)
    V_dfaces = GT.each_face(V,dΓd;tabulate,compute)
    for (V_dface,dΓd_dface) in zip(V_dfaces,dΓd_dfaces)
        dofs = GT.dofs(V_dface)
        ndofs = length(dofs)
        V_points = GT.each_point(V_dface)
        h = GT.diameter(dΓd_dface)
        dΓd_points = GT.each_point(dΓd_dface)
        fill!(Auu,zero(T))
        fill!(bu,zero(T))
        for (V_point,dΓd_point) in zip(V_points,dΓd_points)
            dS = GT.weight(dΓd_point)
            x = GT.coordinate(dΓd_point)
            dof_s = GT.shape_functions(GT.value,V_point)
            dof_∇s = GT.shape_functions(∇,V_point)
            n = GT.unit_normal(V_point)
            gx = g.definition(x)
            for i in 1:ndofs
                ∇v = dof_∇s[i]
                v = dof_s[i]
                bu[i] += ((γ0/h)*v*gx - n⋅∇v*gx)*dS
            end
            for j in 1:ndofs
                ∇u = dof_∇s[j]
                u = dof_s[j]
                for i in 1:ndofs
                    ∇v = dof_∇s[i]
                    v = dof_s[i]
                    Auu[i,j] += ((γ0/h)*v*u - v*n⋅∇u - n⋅∇v*u)*dS
                end
            end
        end
        GT.contribute!(A_alloc,Auu,dofs,dofs)
        GT.contribute!(b_alloc,bu,dofs)
    end

end

function integrate_l2_error(g,uh,dΩ)

    #Iterators to the quantities on the
    #integration points
    tabulate = (GT.value,)
    compute = (GT.coordinate,)
    uh_faces = GT.each_face(uh,dΩ;tabulate,compute)

    #Numerical integration loop
    s = 0.0
    for uh_face in uh_faces
        for uh_point in GT.each_point(uh_face)

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
T = Float64
A_alloc = GT.allocate_matrix(T,V,V,Ω,Λ,Γd)
b_alloc = GT.allocate_vector(T,V,Ω,Γd)

#Fill allocations with the function we wrote above
assemble_Ω!(A_alloc,b_alloc,V,dΩ)
assemble_Λ!(A_alloc,V,dΛ)
assemble_Γd!(A_alloc,b_alloc,V,dΓd)

#Compress matrix into the final format
A = GT.compress(A_alloc)
b = GT.compress(b_alloc)

#Build the linear system
p = LinearSolve.LinearProblem(A,b)

#Solve the problem
sol = LinearSolve.solve(p)
uh = GT.solution_field(V,sol)

#Integrate the error l2 norm
#with the function we wrote above
el2 = integrate_l2_error(g,uh,dΩ)
@assert el2 < 1.0e-9
nothing # hide





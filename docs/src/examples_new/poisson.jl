# # Hello, world!
#
# In this example, we show how to solve the "Hello, world" PDE example:
# the Poisson equation on the unit square with Dirichlet boundary conditions.
#
# ```math
# \left\lbrace
# \begin{aligned}
# -\Delta u = f  \ &\text{in} \ \Omega,\\
# u = g \ &\text{on}\ \partial\Omega,\\
# \end{aligned}
# \right.
# ```
# with $f=0$ and $g(x)=\text{sum}(x)$.
#
#

# ## Automatic assembly
#
# ### Code explained
#
# Load dependencies form Julia stdlib.

using LinearAlgebra

# Import other dependencies

import GalerkinToolkit as GT
import PartitionedSolvers as PS
import ForwardDiff
import GLMakie as Makie
import FileIO # hide

# Generate the computational mesh.

domain = (0,1,0,1)
cells = (10,10)
mesh = GT.cartesian_mesh(domain,cells)
nothing # hide

# Visualize the mesh.

Makie.plot(mesh,color=:pink,strokecolor=:blue)
FileIO.save(joinpath(@__DIR__,"fig_poisson_1.png"),Makie.current_figure()) # hide

# ![](fig_poisson_1.png)

# Define the Dirichlet boundary.

dirichlet_tag = "dirichlet"
GT.label_boundary_faces!(mesh;physical_name=dirichlet_tag)
nothing # hide

# Defile computational domains.

Ω = GT.interior(mesh)
Γd = GT.boundary(mesh;physical_names=[dirichlet_tag])
nothing # hide

# Define differential operators

const ∇ = ForwardDiff.gradient
Δ(f,x) = tr(ForwardDiff.jacobian(y->∇(f,y),x))

# Define manufactured fields.

g = GT.analytical_field_tmp(sum,Ω) # TODO remove _tmp
f = GT.analytical_field_tmp(x->-Δ(g.definition,x),Ω)
nothing # hide

# Define the interpolation space.

k = 1
V = GT.lagrange_space(Ω,k;dirichlet_boundary=Γd)
nothing # hide

# Interpolate Dirichlet values.

T = Float64
uhd = GT.dirichlet_field(T,V)
GT.interpolate_dirichlet!(g,uhd)
nothing # hide

# Visualize the Dirichlet field.

Makie.plot(Ω,color=uhd,strokecolor=:blue)
FileIO.save(joinpath(@__DIR__,"fig_poisson_2.png"),Makie.current_figure()) # hide

# ![](fig_poisson_2.png)

# Define numerical integration.

degree = 2*k
dΩ = GT.measure(Ω,degree)
nothing # hide

# Define weak form.

a = (u,v) -> GT.∫( x->∇(u,x)⋅∇(v,x), dΩ)
l = v -> GT.∫( x->v(x)*f(x), dΩ)
nothing # hide

# Assemble the problem using the automatic assembly loop generator

p = GT.linear_problem(uhd,a,l)
nothing # hide

# Solve the problem

s = PS.LinearAlgebra_lu(p)
s = PS.solve(s)
nothing # hide

# Build the FE solution.

uh = GT.solution_field(uhd,s)
nothing # hide

# Visualize the solution.

Makie.plot(Ω;color=uh,strokecolor=:black)
FileIO.save(joinpath(@__DIR__,"fig_poisson_3.png"),Makie.current_figure()) # hide

# ![](fig_poisson_3.png)

# Compute the L2 norm of the discretization error.

eh = x -> uh(x) - g(x)
el2 = GT.∫( x->abs2(eh(x)), dΩ) |> sum |> sqrt

# ### Final program.

module Automatic 

using LinearAlgebra
import GalerkinToolkit as GT
import PartitionedSolvers as PS
import ForwardDiff
import GLMakie as Makie

function main(;domain,cells)
    mesh = GT.cartesian_mesh(domain,cells)
    dirichlet_tag = "dirichlet"
    GT.label_boundary_faces!(mesh;physical_name=dirichlet_tag)
    Ω = GT.interior(mesh)
    Γd = GT.boundary(mesh;physical_names=[dirichlet_tag])
    ∇ = ForwardDiff.gradient
    Δ(f,x) = tr(ForwardDiff.jacobian(y->∇(f,y),x))
    g = GT.analytical_field_tmp(sum,Ω) # TODO remove _tmp
    f = GT.analytical_field_tmp(x->-Δ(g.definition,x),Ω)
    k = 1
    V = GT.lagrange_space(Ω,k;dirichlet_boundary=Γd)
    T = Float64
    uhd = GT.dirichlet_field(T,V)
    GT.interpolate_dirichlet!(g,uhd)
    degree = 2*k
    dΩ = GT.measure(Ω,degree)
    a = (u,v) -> GT.∫( x->∇(u,x)⋅∇(v,x), dΩ)
    l = v -> GT.∫( x->v(x)*f(x), dΩ)
    p = GT.linear_problem(uhd,a,l)
    s = PS.LinearAlgebra_lu(p)
    s = PS.solve(s)
    uh = GT.solution_field(uhd,s)
    eh = x -> uh(x) - g(x)
    el2 = GT.∫( x->abs2(eh(x)), dΩ) |> sum |> sqrt
end

end # module

# Run it for a 2d case

Automatic.main(domain=(0,1,0,1),cells=(10,10))

# Run it for a 3d case

Automatic.main(domain=(0,1,0,1,0,1),cells=(10,10,10))

# ## Hand-written assembly
#
# ### Code explained

# Create accessor functions for low-level integration quantities.

face_point_x = GT.coordinate_accessor(dΩ)
face_point_J = GT.jacobian_accessor(dΩ)
face_point_dV = GT.weight_accessor(dΩ)
face_npoints = GT.num_points_accessor(dΩ)
nothing # hide

# Create accessor functions for low-level interpolation quantities.

face_dofs = GT.dofs_accessor(V,Ω)
face_point_dof_s = GT.shape_function_accessor(GT.value,V,dΩ)
face_point_dof_∇s = GT.shape_function_accessor(ForwardDiff.gradient,V,dΩ)
face_dirichlet! = GT.dirichlet_accessor(uhd,Ω)
nothing # hide

# Allocate auxiliary element matrix and vector

n = maximum(map(GT.num_dofs,GT.reference_fes(V)))
Auu = zeros(T,n,n)
bu = zeros(T,n)

# Allocate space for the matrix and the vector

b_alloc = GT.allocate_vector(T,V,Ω)
A_alloc = GT.allocate_matrix(T,V,V,Ω)
nothing # hide

# Assembly loop

for face in 1:GT.num_faces(Ω)
    point_x = face_point_x(face)
    point_J = face_point_J(face)
    point_dV = face_point_dV(face)
    point_dof_s = face_point_dof_s(face)
    point_dof_∇s = face_point_dof_∇s(face)
    npoints = face_npoints(face)
    dofs = face_dofs(face)
    dirichlet! = face_dirichlet!(face)
    fill!(Auu,zero(eltype(Auu)))
    fill!(bu,zero(eltype(bu)))
    for point in 1:npoints
        x = point_x(point)
        J = point_J(point)
        dV = point_dV(point,J)
        dof_s = point_dof_s(point)
        dof_∇s = point_dof_∇s(point,J)
        for (i,dofi) in enumerate(dofs)
            v = dof_s(i)
            ∇v = dof_∇s(i)
            bu[i] += f.definition(x)*v*dV
            for (j,dofj) in enumerate(dofs)
                ∇u = dof_∇s(j)
                Auu[i,j] += ∇v⋅∇u*dV
            end
        end
    end
    dirichlet!(Auu,bu)
    GT.contribute!(b_alloc,bu,dofs)
    GT.contribute!(A_alloc,Auu,dofs,dofs)
end

# Compress the matrix and vector

b = GT.compress(b_alloc)
A = GT.compress(A_alloc)
nothing # hide

# Solve the system

free_values = A\b
uh = GT.solution_field(uhd,free_values)
nothing # hide

# Visualize the solution.

Makie.plot(Ω;color=uh,strokecolor=:black)
FileIO.save(joinpath(@__DIR__,"fig_poisson_4.png"),Makie.current_figure()) # hide

# ![](fig_poisson_4.png)

# Compute the L2 norm of the discretization error.

face_point_uhx = GT.discrete_field_accessor(GT.value,uh,dΩ)

err = Ref(0.0)
for face in 1:GT.num_faces(Ω)
    npoints = face_npoints(face)
    point_x = face_point_x(face)
    point_J = face_point_J(face)
    point_dV = face_point_dV(face)
    point_uhx = face_point_uhx(face)
    for point in 1:npoints
        x = point_x(point)
        J = point_J(point)
        dV = point_dV(point,J)
        uhx = point_uhx(point,J)
        err[] += abs2(uhx-g.definition(x))*dV
    end
end

sqrt(err[])

module Handwritten

function setup_example(;domain,cells)
    mesh = GT.cartesian_mesh(domain,cells)
    dirichlet_tag = "dirichlet"
    GT.label_boundary_faces!(mesh;physical_name=dirichlet_tag)
    Ω = GT.interior(mesh)
    Γd = GT.boundary(mesh;physical_names=[dirichlet_tag])
    ∇ = ForwardDiff.gradient
    Δ(f,x) = tr(ForwardDiff.jacobian(y->∇(f,y),x))
    g = GT.analytical_field_tmp(sum,Ω) # TODO remove _tmp
    f = GT.analytical_field_tmp(x->-Δ(g.definition,x),Ω)
    k = 1
    V = GT.lagrange_space(Ω,k;dirichlet_boundary=Γd)
    T = Float64
    uhd = GT.dirichlet_field(T,V)
    GT.interpolate_dirichlet!(g,uhd)
    degree = 2*k
    dΩ = GT.measure(Ω,degree)
    example = (;mesh,Ω,dΩ,V,uhd,T)
end

function setup_accessors(example)
    (;mesh,Ω,dΩ,V,uhd) = example
    face_point_x = GT.coordinate_accessor(dΩ)
    face_point_J = GT.jacobian_accessor(dΩ)
    face_point_dV = GT.weight_accessor(dΩ)
    face_npoints = GT.num_points_accessor(dΩ)
    face_dofs = GT.dofs_accessor(V,Ω)
    face_point_dof_s = GT.shape_function_accessor(GT.value,V,dΩ)
    face_point_dof_∇s = GT.shape_function_accessor(ForwardDiff.gradient,V,dΩ)
    face_dirichlet! = GT.dirichlet_accessor(uhd,Ω)
    accessors = (;
        face_point_x,face_point_J,face_point_J,face_point_dV,face_npoints,
        face_dofs,face_point_dof_s,face_point_dof_∇s,face_dirichlet!)
end

function assemble_problem(accesors,example)
    n = maximum(map(GT.num_dofs,GT.reference_fes(V)))
    Auu = zeros(T,n,n)
    bu = zeros(T,n)
    b_alloc = GT.allocate_vector(T,V,Ω)
    A_alloc = GT.allocate_matrix(T,V,V,Ω)
    for face in 1:GT.num_faces(Ω)
        point_x = face_point_x(face)
        point_J = face_point_J(face)
        point_dV = face_point_dV(face)
        point_dof_s = face_point_dof_s(face)
        point_dof_∇s = face_point_dof_∇s(face)
        npoints = face_npoints(face)
        dofs = face_dofs(face)
        dirichlet! = face_dirichlet!(face)
        fill!(Auu,zero(eltype(Auu)))
        fill!(bu,zero(eltype(bu)))
        for point in 1:npoints
            x = point_x(point)
            J = point_J(point)
            dV = point_dV(point,J)
            dof_s = point_dof_s(point)
            dof_∇s = point_dof_∇s(point,J)
            for (i,dofi) in enumerate(dofs)
                v = dof_s(i)
                ∇v = dof_∇s(i)
                bu[i] += f.definition(x)*v*dV
                for (j,dofj) in enumerate(dofs)
                    ∇u = dof_∇s(j)
                    Auu[i,j] += ∇v⋅∇u*dV
                end
            end
        end
        dirichlet!(Auu,bu)
        GT.contribute!(b_alloc,bu,dofs)
        GT.contribute!(A_alloc,Auu,dofs,dofs)
    end
    b = GT.compress(b_alloc)
    A = GT.compress(A_alloc)
    x = similar(b,axes(A,2))
    p = PS.linear_problem(x,A,b)
end

function 



end # module


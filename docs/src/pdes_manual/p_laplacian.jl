# # p-Laplacian
#
# ## Problem statement
#
# Solve the following p-Laplacian equation on the unit square,
#
# ```math
# \left\lbrace
# \begin{aligned}
# -\nabla \cdot \left( |\nabla u|^{q-2} \ \nabla u \right) = f\ &\text{in}\ \Omega,\\
# u = g \ &\text{on}\ \partial\Omega,\\
# \end{aligned}
# \right.
# ```
# with $f=1$ and $g=0$ and $q=3$.
#
# Solve it with a piece-wise bi-linear Lagrange interpolation, and visualize the result.
#

# ## Implementation
#
# Load dependencies form Julia stdlib.

using LinearAlgebra

# Import other dependencies

import GalerkinToolkit as GT
import PartitionedSolvers as PS
import ForwardDiff
import GLMakie as Makie

# Setup the objects defining this example

function setup_example(;domain,cells,file)
    mesh = GT.cartesian_mesh(domain,cells)
    dirichlet_tag = "dirichlet"
    GT.label_boundary_faces!(mesh;physical_name=dirichlet_tag)
    Ω = GT.interior(mesh)
    Γd = GT.boundary(mesh;physical_names=[dirichlet_tag])
    g = GT.analytical_field(x->0,Ω)
    f = GT.analytical_field(x->1,Ω)
    k = 1
    V = GT.lagrange_space(Ω,k;dirichlet_boundary=Γd)
    T = Float64
    uh = GT.rand_field(T,V)
    GT.interpolate_dirichlet!(g,uh)
    degree = 2*k
    dΩ = GT.measure(Ω,degree)
    example = (;mesh,Ω,dΩ,V,uh,T,f,g,file)
    state = (;example)
end
nothing # hide

# Create accessor functions for low-level integration quantities.

function setup_integration_accessors(state)
    (;dΩ) = state.example
    face_point_x = GT.coordinate_accessor(dΩ)
    face_point_J = GT.jacobian_accessor(dΩ)
    face_point_dV = GT.weight_accessor(dΩ)
    face_npoints = GT.num_points_accessor(dΩ)
    integration = (;
        face_point_x,face_point_J, face_point_dV,face_npoints)
    (;integration,state...)
end
nothing # hide

# Create accessor functions for low-level interpolation quantities.

function setup_interpolation_accessors(state)
    (;dΩ,Ω,V,uh) = state.example
    face_dofs = GT.dofs_accessor(V,Ω)
    face_point_dof_s = GT.shape_function_accessor(GT.value,V,dΩ)
    face_point_dof_∇s = GT.shape_function_accessor(ForwardDiff.gradient,V,dΩ)
    face_point_∇uh = GT.discrete_field_accessor(ForwardDiff.gradient,uh,Ω)
    interpolation = (;
        face_dofs,face_point_dof_s,face_point_dof_∇s,face_point_∇uh)
    (;interpolation,state...)
end
nothing # hide

# Define the algebraic non-linear problem object

function setup_nonlinear_problem(state)

    (;example) = state
    (;V,Ω,T,f) = example

    #Allocate auxiliary face matrix and vector
    n = maximum(map(GT.num_dofs,GT.reference_fes(V)))
    Auu = zeros(T,n,n)
    bu = zeros(T,n)

    #Allocate space for the global matrix and vector
    b_alloc = GT.allocate_vector(T,V,Ω)
    A_alloc = GT.allocate_matrix(T,V,V)

    #Fill in residual and jacobian according to the initial state
    fill_residual_and_jacobian!(state,allocs) # Defined later

    #Compress matrix and vector into the final format
    b,b_cache = GT.compress(b_alloc;reuse=Val(true))
    A,A_cache = GT.compress(A_alloc;reuse=Val(true))

    allocs = (;A_alloc,b_alloc,Auu,bu)

    workspace = nothing
    problem = PS.nonlinear_problem(x0,r0,j0,workspace) do p

        # Get the current solution vector
        x = PS.solution(p)

        # Build the current solution field
        uh = state.example.uh
        GT.solution_field!(uh,x)

        # Fill in residual and Jacobian
        fill_residual_and_jacobian!(state,allocs)

        # In-place compression of matrix and vector
        GT.compress!(b,b_cache,b_alloc)
        GT.compress!(A,A_cache,A_alloc)

        # Update the nonlinear problem object
        # Here, we computed the residual and Jacobian simultaneously,
        # but this API also allows to compute them separately.
        if PS.residual(p) !== nothing
            p = PS.update(p,b)
        end
        if PS.jacobian(p) !== nothing
            p = PS.update(p,A)
        end

        p
    end

    (;problem,state...)
end
nothing # hide

# Assembly loop

function fill_residual_and_jacobian!(state,allocs)
    (;A_alloc,b_alloc,Auu,bu) = allocs
    (;example) = state
    (;Ω,) = example

    # Reset allocations
    GT.reset!(b_alloc)
    GT.reset!(A_alloc)

    #Loop over the faces of the domain
    for face in 1:GT.num_faces(Ω)

        #Compute face tensors
        dofs = face_tensors!(Auu,bu,face,state) # Defined later

        #Add face contribution to global allocation
        GT.contribute!(b_alloc,bu,dofs)
        GT.contribute!(A_alloc,Auu,dofs,dofs)
    end

end
nothing # hide

# Compute local Jacobian and residual

function face_tensors!(Auu,bu,face,state)

    (;example,integration,interpolation) = state
    (;f) = example

    # Define flux and its derivative
    q = 3
    flux(∇u) = norm(∇u)^(q-2) * ∇u
    dflux(∇du,∇u) = (q-2)*norm(∇u)^(q-4)*(∇u⋅∇du)*∇u+norm(∇u)^(q-2)*∇du

    #Get quantities at current face
    npoints = integration.face_npoints(face)
    point_x = integration.face_point_x(face)
    point_J = integration.face_point_J(face)
    point_dV = integration.face_point_dV(face)
    point_dof_s = interpolation.face_point_dof_s(face)
    point_dof_∇s = interpolation.face_point_dof_∇s(face)
    point_∇uh = interpolation.face_point_∇uh(face)
    dofs = interpolation.face_dofs(face)

    #Reset face matrix and vector
    fill!(Auu,zero(eltype(Auu)))
    fill!(bu,zero(eltype(bu)))

    #Loop over integration points
    for point in 1:npoints

        #Get quantities at current integration point
        x = point_x(point)
        fx = f.definition(x)
        J = point_J(point)
        dV = point_dV(point,J)
        dof_s = point_dof_s(point)
        dof_∇s = point_dof_∇s(point,J)
        point_∇uh = point_∇uh(point,J)

        #Fill in face matrix and vector
        for (i,dofi) in enumerate(dofs)
            v = dof_s(i)
            ∇v = dof_∇s(i)
            bu[i] += flux(∇uh)⋅∇v - fx*v*dV
            for (j,dofj) in enumerate(dofs)
                ∇du = dof_∇s(j)
                Auu[i,j] += dflux(∇du,∇uh)⋅∇v*dV
            end
        end
    end

    #Return dof ids for this face
    dofs
end
nothing # hide

# Solve and visualize results

function solve_and_visualize!(state)
    (;problem,example) = state
    (;uh,file) = example

    # Setup solver
    solver = PS.newton_raphson(problem,verbose=true)

    # Get the lazy solver history
    solver_history = PS.history(solver)

    # Visualize
    color = Makie.Observable(uh)
    fig = Makie.plot(Ω;color,strokecolor=:black)
    fn = joinpath(@__DIR__,"fig_pt_plaplacian.gif")
    Makie.record(fig,file,solver_history;framerate=2) do s
        GT.solution_field!(uh,s)
        color[] = uh
    end
end

# Final program.

function main(;kwargs...)
    state1 = setup_example(;kwargs...)
    state2 = setup_integration_accessors(state1)
    state3 = setup_interpolation_accessors(state2)
    state4 = setup_nonlinear_problem(state3)
    solve_and_visualize!(state4)
end
nothing # hide

# Run it for a 2d case

file = joinpath(@__DIR__,"p_laplacian_2d_manual.gif")
Program.main(domain=(0,1,0,1),cells=(10,10),file)

# ![](p_laplacian_2d_manual.gif)




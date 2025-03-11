# # Hello, World! (manual assembly)
#
# ## Problem statement
# We solve the same problem as in the [Hello, World!](@ref) example, but in this case we explicitly write the numerical integration loops.
#
# ## Implementation

# Load dependencies form Julia stdlib.

using LinearAlgebra

# Import other dependencies

import GalerkinToolkit as GT
import PartitionedSolvers as PS
import ForwardDiff

# Setup the objects defining this example

function setup_example(;domain,cells)
    mesh = GT.cartesian_mesh(domain,cells)
    dirichlet_tag = "dirichlet"
    GT.label_boundary_faces!(mesh;physical_name=dirichlet_tag)
    Ω = GT.interior(mesh)
    Γd = GT.boundary(mesh;physical_names=[dirichlet_tag])
    ∇ = ForwardDiff.gradient
    Δ(f,x) = tr(ForwardDiff.jacobian(y->∇(f,y),x))
    g = GT.analytical_field(sum,Ω)
    f = GT.analytical_field(x->-Δ(g.definition,x),Ω)
    k = 1
    V = GT.lagrange_space(Ω,k;dirichlet_boundary=Γd)
    T = Float64
    uhd = GT.zero_dirichlet_field(T,V)
    GT.interpolate_dirichlet!(g,uhd)
    degree = 2*k
    dΩ = GT.measure(Ω,degree)
    example = (;mesh,Ω,dΩ,V,uhd,f,g)
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
    (;dΩ,Ω,V,uhd) = state.example
    face_dofs = GT.dofs_accessor(V,Ω)
    face_point_dof_s = GT.shape_function_accessor(GT.value,V,dΩ)
    face_point_dof_∇s = GT.shape_function_accessor(ForwardDiff.gradient,V,dΩ)
    face_dirichlet! = GT.dirichlet_accessor(uhd,Ω)
    interpolation = (;
        face_dofs,face_point_dof_s,face_point_dof_∇s,face_dirichlet!)
    (;interpolation,state...)
end
nothing # hide

# Compute local matrix and vector

function face_tensors!(Auu,bu,face,state)

    (;example,integration,interpolation) = state
    (;f) = example

    #Get quantities at current face
    npoints = integration.face_npoints(face)
    point_x = integration.face_point_x(face)
    point_J = integration.face_point_J(face)
    point_dV = integration.face_point_dV(face)
    point_dof_s = interpolation.face_point_dof_s(face)
    point_dof_∇s = interpolation.face_point_dof_∇s(face)
    dofs = interpolation.face_dofs(face)
    dirichlet! = interpolation.face_dirichlet!(face)

    #Reset face matrix and vector
    fill!(Auu,zero(eltype(Auu)))
    fill!(bu,zero(eltype(bu)))

    #Loop over integration points
    for point in 1:npoints

        #Get quantities at current integration point
        x = point_x(point)
        J = point_J(point)
        dV = point_dV(point,J)
        dof_s = point_dof_s(point)
        dof_∇s = point_dof_∇s(point,J)

        #Fill in face matrix and vector
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

    #Apply Dirichlet conditions on face vector
    dirichlet!(Auu,bu)

    #Return dof ids for this face
    dofs
end
nothing # hide

# Assemble linear problem

function assemble_problem(state)

    (;example) = state
    (;V,Ω,f) = example

    #Allocate auxiliary face matrix and vector
    n = GT.max_num_reference_dofs(V)
    T = Float64
    Auu = zeros(T,n,n)
    bu = zeros(T,n)

    #Allocate space for the global matrix and vector
    b_alloc = GT.allocate_vector(T,V,Ω)
    A_alloc = GT.allocate_matrix(T,V,V,Ω)

    #Loop over the faces of the domain
    for face in 1:GT.num_faces(Ω)

        #Compute face tensors
        dofs = face_tensors!(Auu,bu,face,state)

        #Add face contribution to global allocation
        GT.contribute!(b_alloc,bu,dofs)
        GT.contribute!(A_alloc,Auu,dofs,dofs)
    end

    #Compress matrix and vector into the final format
    b = GT.compress(b_alloc)
    A = GT.compress(A_alloc)

    #Allocate space for the solution
    sol = similar(b,axes(A,2))

    #Create linear problem object
    problem = PS.linear_problem(sol,A,b)
    (;problem,state...)
end
nothing # hide

# Solve the solution problem and build the solution field

function solve_problem(state)
    (;problem,example) = state
    (;uhd) = example
    solver = PS.LinearAlgebra_lu(problem)
    solver = PS.solve(solver)
    uh = GT.solution_field(uhd,solver)
    (;uh,state...)
end
nothing # hide

# Setup accessor functions to integrate error norms

function setup_postpro_accessors(state)
    (;uh,example) = state
    (;dΩ) = example
    #TODO add also h1 norm
    face_point_uhx = GT.discrete_field_accessor(GT.value,uh,dΩ)
    postpro = (;face_point_uhx)
    (;postpro,state...)
end
nothing # hide

# Compute the error norms

function integrate_error_norms(state)
    (;postpro,example,integration) = state
    (;Ω,dΩ,g) = example
    #TODO add also h1 norm
    T = Float64
    el2 = zero(T)
    for face in 1:GT.num_faces(Ω)
        npoints = integration.face_npoints(face)
        point_x = integration.face_point_x(face)
        point_J = integration.face_point_J(face)
        point_dV = integration.face_point_dV(face)
        point_uhx = postpro.face_point_uhx(face)
        for point in 1:npoints
            x = point_x(point)
            J = point_J(point)
            dV = point_dV(point,J)
            uhx = point_uhx(point,J)
            el2 += abs2(uhx-g.definition(x))*dV
        end
    end
    el2 = sqrt(el2)
    norms = (;el2)
end
nothing # hide

# ## Final Program

function main(;kwargs...)
    state1 = setup_example(;kwargs...)
    state2 = setup_integration_accessors(state1)
    state3 = setup_interpolation_accessors(state2)
    state4 = assemble_problem(state3)
    state5 = solve_problem(state4)
    state6 = setup_postpro_accessors(state5)
    norms = integrate_error_norms(state6)
end
nothing # hide

# Run it for a 2d case

main(domain=(0,1,0,1),cells=(10,10))

# Run it for a 3d case

main(domain=(0,1,0,1,0,1),cells=(10,10,10))




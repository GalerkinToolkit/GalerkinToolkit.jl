using LinearAlgebra
import GalerkinToolkit as GT
import PartitionedSolvers as PS
import ForwardDiff
# import GLMakie as Makie



mean(u) = 0.5*(u[1]+u[2])
jump(u,n) = u[2]*n[2] + u[1]*n[1]

function poisson_dg_assembly(n, k, iscube)
    outputs = Dict()
    if iscube
        domain=(0, 1, 0, 1, 0, 1)
        cells=(n, n, n)
        mesh = GT.cartesian_mesh(domain, cells) # TODO: use tetrahedron
    else
        error("loading from file not implemented!")
    end
    D = GT.num_dims(mesh)
    degree = 2*k
    γ = degree*(degree+1)
    γ = γ/10.0

    time = @elapsed begin
        dirichlet_tag = "dirichlet"
        GT.label_boundary_faces!(mesh;physical_name=dirichlet_tag)
        Ω = GT.interior(mesh)
        Γd = GT.boundary(mesh;physical_names=[dirichlet_tag])
        dΩ = GT.measure(Ω,degree)


        conformity = :L2
        GT.label_interior_faces!(mesh;physical_name="__INTERIOR_FACES__")
        Λ = GT.skeleton(mesh;physical_names=["__INTERIOR_FACES__"])
        dΛ = GT.measure(Λ,degree)
        n_Λ = GT.unit_normal(mesh,D-1)
        h_Λ = GT.face_diameter_field(Λ)

        n_Γd = GT.unit_normal(mesh,D-1)
        h_Γd = GT.face_diameter_field(Γd)
        dΓd = GT.measure(Γd,degree)
        V = GT.lagrange_space(Ω,k;conformity)



        # analytical field by k
        if k == 1
            u = GT.analytical_field(x -> 1 + x[1] + 2 * x[2],Ω)
            f = x -> 0.0
            # f = GT.analytical_field(x -> 0.0,Ω) # TODO: use constant
        elseif k == 2
            u = GT.analytical_field(x -> 1 + x[1] * x[1] + 2 * x[2] * x[2],Ω)
            # f = GT.analytical_field(x -> -6.0,Ω)
            f = x -> -6.0
        elseif k == 3
            u = GT.analytical_field(x -> 1 + x[1]*x[1]*x[1] + 2*x[2]*x[2]*x[2],Ω)
            f = GT.analytical_field(x -> -6*x[1] - 12*x[2],Ω)
        else
            error("k invalid")
        end

        ∇ = ForwardDiff.gradient
  

        V = GT.lagrange_space(Ω,k;dirichlet_boundary=Γd)
        T = Float64
        uhd = GT.dirichlet_field(T,V)
        GT.interpolate_dirichlet!(u,uhd)
        
        
        a = (u,v) -> GT.∫( x->∇(u,x)⋅∇(v,x), dΩ) +
                GT.∫( x -> (γ/h_Γd(x))*v(x)*u(x)-v(x)*n_Γd(x)⋅∇(u,x)-n_Γd(x)⋅∇(v,x)*u(x), dΓd) + 
                GT.∫( x -> (γ/h_Λ(x))*jump(v(x),n_Λ(x))⋅jump(u(x),n_Λ(x))-jump(v(x),n_Λ(x))⋅mean(∇(u,x))-mean(∇(v,x))⋅jump(u(x),n_Λ(x)), dΛ)
        l = v -> GT.∫( x->v(x)*f(x), dΩ) + 
            GT.∫( x->(γ/h_Γd(x))*v(x)*u(x)-n_Γd(x)⋅∇(v,x)*u(x), dΓd)
    end
    outputs["precomputing [s]"] = time


    # p = GT.linear_problem(Float64,V,a,l)

    time = @elapsed begin
        U = V
        A,b,cache = GT.assemble_matrix_and_vector(a,l,U,V,Float64;reuse=Val(true))
        x = similar(b,axes(A,2))
        p = PS.linear_problem(x,A,b)
    end
    outputs["assembly [s]"] = time

    time = @elapsed begin
        U = V
        A,b = GT.assemble_matrix_and_vector!(a,l,U,V,A,b,cache)
        p = PS.linear_problem(x,A,b)
    end
    outputs["assembly! [s]"] = time


    time = @elapsed begin
        s = PS.LinearAlgebra_lu(p)
        s = PS.solve(s)
    end

    outputs["solve [s]"] = time

    uh = GT.solution_field(uhd,s)
    eh = x -> uh(x) - u(x)
    el2 = GT.∫( x->abs2(eh(x)), dΩ) |> sum |> sqrt


    outputs["el2"] = el2
    # outputs["num_cells"] = num_cells(model)
    # outputs["num_dofs"] = length(b)
    # outputs["num_nz"] = nnz(A)
    outputs["n"] = n
    outputs["k"] = k
    outputs["experiment"] = "poisson_dg_assembly"
    outputs["lib"] = "GT"
    outputs["iscube"] = iscube

    return outputs
end
using LinearAlgebra
import GalerkinToolkit as GT
import PartitionedSolvers as PS
import ForwardDiff
# import GLMakie as Makie



function poisson_cg_assembly(n, k, iscube)
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

    time = @elapsed begin
        dirichlet_tag = "dirichlet"
        GT.label_boundary_faces!(mesh;physical_name=dirichlet_tag)
        Ω = GT.interior(mesh)
        Γd = GT.boundary(mesh;physical_names=[dirichlet_tag])
        dΩ = GT.measure(Ω,degree)

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
            u = GT.analytical_field(x -> -6*x[1] - 12*x[2],Ω)
        else
            error("k invalid")
        end

        ∇ = ForwardDiff.gradient
        V = GT.lagrange_space(Ω,k;dirichlet_boundary=Γd)
        T = Float64
        uhd = GT.dirichlet_field(T,V)
        GT.interpolate_dirichlet!(u,uhd)
        
        
        a = (u,v) -> GT.∫( x->∇(u,x)⋅∇(v,x), dΩ)
        l = v -> GT.∫( x->v(x)*f(x), dΩ)
    end
    outputs["precomputing [s]"] = time

    time = @elapsed begin
        U = GT.space(uhd)
        V = GT.space(uhd)
        A,b,cache = GT.assemble_matrix_and_vector_with_dirichlet(a,l,U,V,GT.dirichlet_values(uhd);reuse=Val(true))
        x = similar(b,axes(A,2))
        p = PS.linear_problem(x,A,b)
    end
    outputs["assembly [s]"] = time

    time = @elapsed begin
        U = GT.space(uhd)
        V = GT.space(uhd)
        A,b = GT.assemble_matrix_and_vector_with_dirichlet!(a,l,U,V,GT.dirichlet_values(uhd),A,b,cache)
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
    outputs["experiment"] = "poisson_cg_assembly"
    outputs["lib"] = "GT"
    outputs["iscube"] = iscube

    return outputs
end
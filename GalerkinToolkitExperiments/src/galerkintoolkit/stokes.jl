using LinearAlgebra
import GalerkinToolkit as GT
import PartitionedSolvers as PS
import ForwardDiff
using StaticArrays
# import GLMakie as Makie

laplacian(u,x) = tr(ForwardDiff.jacobian(y->ForwardDiff.gradient(u,y),x))
# TODO hide this
laplacian(u::GT.AbstractQuantity,x::GT.AbstractQuantity) = GT.call(laplacian,u,x)
Δ = laplacian

function stokes_assembly(n, k, iscube)
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

            
        if k == 2
            u = GT.analytical_field(x -> [(x[1] * x[1] + 2 * x[2] * x[2]), (-x[2] * x[2]), 0.0], Ω) # TODO: array type?
            v = GT.analytical_field(x -> x[1] + 3 * x[2], Γd)
            f = GT.analytical_field(x -> [-6.0+1.0, 2.0+3.0, 0.0], Ω)
            # f = x -> [-6.0+1.0, 2.0+3.0, 0.0] # TODO: static array
            g = GT.analytical_field(x -> 2 * x[1] - 2 * x[2], Γd)
        else
            error("k invalid")
        end

        ∇ = ForwardDiff.gradient
        V = GT.lagrange_space(Ω,k;shape=(3,), dirichlet_boundary=Γd) # TODO: examples
        T = Float64
        uhd = GT.dirichlet_field(Float64,V)
        GT.interpolate_dirichlet!(u,uhd)

        Q = GT.lagrange_space(Γd,k-1;conformity=:L2)
        VxQ = V × Q
        uhd2 = GT.dirichlet_field(Float64,Q)
        GT.interpolate_dirichlet!(v,uhd2)
        

        # a = (ufl.inner(ufl.grad(du), ufl.grad(dv)) - ufl.div(dv)*dp + dq*ufl.div(du))*ufl.dx
        # L = (ufl.inner(f, dv) + g*dq )*ufl.dx
        div(f, x) = tr(ForwardDiff.jacobian(f, x))

        function a((u, p),(v, q))
            GT.∫( x-> div(u, x)⋅div(v, x), dΩ) # - div(v, x) * p(x) + div(u, x) * q(x)
        end

        function l((v, q))
            GT.∫( x-> v(x)⋅f(x), dΩ) #  + g(x) * q(x)
        end
    end
    outputs["precomputing [s]"] = time

    time = @elapsed begin
        A,b,cache = GT.assemble_matrix_and_vector(a,l,VxQ,VxQ,Float64;reuse=Val(true))
        x = similar(b,axes(A,2))
        p = PS.linear_problem(x,A,b)
    end
    outputs["assembly [s]"] = time

    time = @elapsed begin
        A,b = GT.assemble_matrix_and_vector!(a,l,VxQ,VxQ,A,b,cache)
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
    outputs["experiment"] = "stokes_assembly"
    outputs["lib"] = "GT"
    outputs["iscube"] = iscube

    return outputs
end
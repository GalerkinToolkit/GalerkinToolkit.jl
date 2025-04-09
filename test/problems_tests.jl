module PoblemsTests

import GalerkinToolkit as GT
using Test
using LinearAlgebra
import ForwardDiff
import PartitionedSolvers as PS

function problems_tests_main()
    ∇ = ForwardDiff.gradient
    domain = (0,1,0,1)
    cells = (2,2)
    T = Float64
    mesh = GT.cartesian_mesh(domain,cells)

    k = 1
    degree = 2*k
    Ω = GT.interior(mesh)
    Γd = GT.boundary(mesh)
    Λ = GT.skeleton(mesh)
    dΩ = GT.measure(Ω,degree)
    dΓd = GT.measure(Γd,degree)
    dΛ = GT.measure(Λ,degree)

    n = GT.unit_normal(mesh,1)
    face_point_v = GT.sample(n,dΓd)
    face_point_v = GT.sample(n[1],dΛ)
    face_point_v = GT.sample(n[2],dΛ)

    f = GT.analytical_field(Ω) do x
        sum(x)
    end

    face_point_v = GT.sample(f,dΩ)
    display(face_point_v)
    face_point_v,cache = GT.sample(f,dΩ;reuse=Val(true))
    @time GT.update_sample!(face_point_v,cache)

    int = GT.∫(x-> GT.call(sum, x),dΩ)
    c = GT.assemble_scalar(int)
    @test c ≈ 1.0 

    face_to_c = GT.face_contribution(int,Ω)
    @test sum(face_to_c) ≈ 1


    #TODO: tabulate jacobian
    V = GT.lagrange_space(Ω,k)
    a = (u, v) -> GT.∫(x-> ∇(u,x) ⋅ ∇(v, x),dΩ)
    A = GT.assemble_matrix(a,T,V,V)

    l = v -> GT.∫(x->v(x),dΩ)
    a = (u, v) -> GT.∫(x->u(x) ⋅ v(x),dΩ) 

    b = GT.assemble_vector(l,T,V)
    @test sum(b) ≈ 1

    A = GT.assemble_matrix(a,T,V,V)
    @test sum(A) ≈ 1
    


    b,bcache = GT.assemble_vector(l,T,V;reuse=Val(true))
    @test sum(b) ≈ 1
    fill!(b, 0.0)
    @time GT.update_vector!(b,bcache)
    @time GT.update_vector!(b,bcache)
    @time GT.update_vector!(b,bcache)

    @test sum(b) ≈ 1


    A,Acache = GT.assemble_matrix(a,T,V,V;reuse=Val(true))
    @test sum(A) ≈ 1
    fill!(A, 0.0)
    @time GT.update_matrix!(A,Acache)
    @time GT.update_matrix!(A,Acache)

    @test sum(A) ≈ 1

    uhd = GT.zero_dirichlet_field(T, V)
    p = GT.linear_problem(uhd,a,l)
    

    Δ(f,x) = tr(ForwardDiff.jacobian(y->∇(f,y),x))
    g = GT.analytical_field(sum,Ω)
    f = GT.analytical_field(x->-Δ(g.definition,x),Ω)
    V = GT.lagrange_space(Ω,k;dirichlet_boundary=Γd)
    uhd = GT.zero_dirichlet_field(T,V)
    GT.interpolate_dirichlet!(g,uhd)
    a = (u,v) -> GT.∫( x->∇(u,x)⋅∇(v,x), dΩ)
    l = v -> GT.∫( x->v(x)*f(x), dΩ)
    p = GT.linear_problem(uhd,a,l)
    s = PS.LinearAlgebra_lu(p)
    s = PS.solve(s)
    uh = GT.solution_field(uhd,s)
    eh = x -> uh(x) - g(x)
    el2 = GT.∫( x->abs2(eh(x)), dΩ) |> sum |> sqrt

    @test el2 + 1.0 ≈ 1.0

end

problems_tests_main()


end # module

module PoblemsTests

import GalerkinToolkit as GT
using Test
using LinearAlgebra


function main()

    domain = (0,1,0,1)
    cells = (2,2)
    mesh = GT.cartesian_mesh(domain,cells)

    k = 1

    Ω = GT.interior(mesh)
    dΩ = GT.measure(Ω,2*k)

    V = GT.lagrange_space(Ω,k)

    l = v -> GT.∫(x->v(x),dΩ)
    a = (u, v) -> GT.∫(x->u(x) ⋅ v(x),dΩ)

    T = Float64
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

end

main()
main()


end # module

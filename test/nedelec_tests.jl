module NedelecTests

import GalerkinToolkit as GT
import ForwardDiff
using Test
using LinearAlgebra
using StaticArrays
using SparseArrays

curl(J::SMatrix{3,3}) = SVector(J[3,2]-J[2,3],J[1,3]-J[3,1],J[2,1]-J[1,2])
curl(J::SMatrix{2,2}) = J[2,1]-J[1,2]

# Reference element

for D in (2,3)
    geo = GT.unit_simplex(Val(D))
    fe = GT.nedelec_space(geo,0)
    nedges = D == 2 ? 3 : 6
    @test GT.num_dofs(fe) == nedges

    # Shape functions are dual to the edge moments.
    shapes = GT.shape_functions(fe)
    sigma = GT.dual_basis(fe)
    @test [sigma[i](shapes[j]) for i in 1:nedges, j in 1:nedges] ≈ Matrix(I,nedges,nedges) atol=1e-12
    # The moments are linear (not conjugate-linear) in complex fields.
    @test [sigma[i](x->im*shapes[j](x)) for i in 1:nedges, j in 1:nedges] ≈ im*Matrix(I,nedges,nedges) atol=1e-12

    # Tangential component of shape j along edge i is δ_ij at every point of the edge.
    cell_mesh = GT.mesh(geo)
    node_x = GT.node_coordinates(cell_mesh)
    edge_nodes = GT.face_nodes(cell_mesh,1)
    for i in 1:nedges, j in 1:nedges, s in (0.0,0.25,0.7,1.0)
        a,b = edge_nodes[i]
        t = node_x[b]-node_x[a]
        @test dot(shapes[j](node_x[a]+s*t),t) ≈ (i==j) atol=1e-12
    end

    # Lowest-order functions have a constant curl; the curl of the Whitney function
    # of edge (a,b) is 2∇λa×∇λb.
    xs = [SVector{D}(ntuple(k->0.1*k,D)),SVector{D}(ntuple(k->0.05+0.02k,D)),SVector{D}(ntuple(k->k==1 ? 0.6 : 0.1,D))]
    grad_λ = vcat([-ones(SVector{D})],[SVector{D}(ntuple(k->Float64(k==m),D)) for m in 1:D])
    for j in 1:nedges
        a,b = edge_nodes[j]
        if D == 3
            expected = 2*cross(grad_λ[a],grad_λ[b])
        else
            expected = 2*(grad_λ[a][1]*grad_λ[b][2]-grad_λ[a][2]*grad_λ[b][1])
        end
        for x in xs
            @test curl(ForwardDiff.jacobian(shapes[j],x)) ≈ expected atol=1e-12
        end
    end

    # The tabulator agrees with the shape functions.
    tab = GT.tabulator(fe)(GT.value,xs)
    @test all(tab[k,j] ≈ shapes[j](xs[k]) for k in 1:length(xs), j in 1:nedges)
    tab = GT.tabulator(fe)(ForwardDiff.jacobian,xs)
    @test all(tab[k,j] ≈ ForwardDiff.jacobian(shapes[j],xs[k]) for k in 1:length(xs), j in 1:nedges)

    # DOF ownership: one DOF per edge.
    @test GT.face_own_dofs(fe,1) == [[i] for i in 1:nedges]
    @test all(isempty,GT.face_own_dofs(fe,0))
    @test only(GT.face_dofs(fe,D)) == 1:nedges
    @test all(isempty,GT.face_own_dofs(fe,D))
    @test all(perms->length(perms)==2 && all(==([1]),perms),GT.face_own_dof_permutations(fe,1))
    if D == 3
        face_dofs = GT.face_dofs(fe,2)
        @test length(face_dofs) == 4
        @test all(dofs->length(dofs)==3,face_dofs)
        @test all(edge->count(dofs->edge in dofs,face_dofs)==2,1:nedges)
        face_nodes = GT.face_nodes(cell_mesh,2)
        for (dofs,nodes) in zip(face_dofs,face_nodes)
            @test all(edge->all(in(nodes),edge_nodes[edge]),dofs)
        end
    end
end

@test_throws ArgumentError GT.nedelec_space(GT.unit_simplex(Val(3)),-1)
@test_throws ArgumentError GT.nedelec_space(GT.unit_simplex(Val(3)),1.5)
@test_throws ArgumentError GT.nedelec_space(GT.unit_n_cube(Val(3)),0)
@test_throws ArgumentError GT.nedelec_space(GT.interior(GT.cartesian_mesh((0,1,0,1),(2,2))),0)

# Independent reference assembly: closed-form Whitney matrices from barycentric
# gradients, with edges oriented from the lower to the higher global node ID.

function whitney_matrices(mesh,D,dof_to_edge)
    node_x = GT.node_coordinates(mesh)
    cell_nodes = GT.face_nodes(mesh,D)
    edge_to_dof = Dict(edge=>dof for (dof,edge) in enumerate(dof_to_edge))
    n = length(dof_to_edge)
    M = zeros(n,n)
    K = zeros(n,n)
    G = zeros(n,n) # ∫ ∇φe : ∇φf, with ∇φ of edge (a,b) equal to ∇λb ∇λaᵀ - ∇λa ∇λbᵀ
    local_edges = D == 3 ? [(1,2),(1,3),(1,4),(2,3),(2,4),(3,4)] : [(1,2),(1,3),(2,3)]
    c = D == 3 ? 20 : 12 # ∫λiλj = |K|(1+δij)/((D+1)(D+2))
    for nodes in cell_nodes
        x = node_x[nodes]
        B = hcat([x[k+1]-x[1] for k in 1:D]...)
        vol = abs(det(B))/factorial(D)
        Binv = inv(B) # rows are ∇λ_2,...,∇λ_{D+1}
        ∇λ = vcat([-vec(sum(Binv,dims=1))],[Binv[k,:] for k in 1:D])
        m(i,j) = vol*(1+(i==j))/c
        g(i,j) = dot(∇λ[i],∇λ[j])
        rot(i,j) = D == 3 ? 2*cross(∇λ[i],∇λ[j]) : 2*(∇λ[i][1]*∇λ[j][2]-∇λ[i][2]*∇λ[j][1])
        edges = map(local_edges) do (i,j)
            nodes[i] < nodes[j] ? (i,j) : (j,i)
        end
        dofs = map(edges) do (i,j)
            edge_to_dof[(nodes[i],nodes[j])]
        end
        for (e,(a,b)) in enumerate(edges), (f,(p,q)) in enumerate(edges)
            M[dofs[e],dofs[f]] += m(a,p)*g(b,q) - m(a,q)*g(b,p) - m(b,p)*g(a,q) + m(b,q)*g(a,p)
            K[dofs[e],dofs[f]] += vol*dot(rot(a,b),rot(p,q))
            Je = ∇λ[b]*∇λ[a]' - ∇λ[a]*∇λ[b]'
            Jf = ∇λ[q]*∇λ[p]' - ∇λ[p]*∇λ[q]'
            G[dofs[e],dofs[f]] += vol*sum(Je .* Jf)
        end
    end
    M,K,G
end

frobenius(A,B) = sum(A .* B)

# GT DOF -> global edge as a sorted node pair; each DOF must map to one edge.
function dof_to_edge_map(V,mesh,D)
    fe = only(GT.reference_spaces(V))
    edge_nodes = GT.face_nodes(GT.mesh(GT.domain(fe)),1)
    cell_nodes = GT.face_nodes(mesh,D)
    face_dofs = GT.face_dofs(V)
    ndofs = GT.num_free_dofs(V)
    dof_to_edge = Vector{Tuple{Int,Int}}(undef,ndofs)
    assigned = falses(ndofs)
    consistent = true
    for cell in 1:length(cell_nodes)
        nodes = cell_nodes[cell]
        for (ledge,(a,b)) in enumerate(edge_nodes)
            dof = face_dofs[cell][ledge]
            edge = minmax(Int(nodes[a]),Int(nodes[b]))
            if assigned[dof]
                consistent &= dof_to_edge[dof] == edge
            else
                dof_to_edge[dof] = edge
                assigned[dof] = true
            end
        end
    end
    @test consistent
    @test all(assigned)
    @test allunique(dof_to_edge)
    dof_to_edge
end

# Dirichlet DOF -> global edge (sorted node pair).
function dirichlet_edges(V,mesh,D)
    fe = only(GT.reference_spaces(V))
    edge_nodes = GT.face_nodes(GT.mesh(GT.domain(fe)),1)
    cell_nodes = GT.face_nodes(mesh,D)
    edges = Vector{Tuple{Int,Int}}(undef,GT.num_dirichlet_dofs(V))
    for (cell,dofs) in enumerate(GT.face_dofs(V)), (ledge,(a,b)) in enumerate(edge_nodes)
        dof = dofs[ledge]
        dof < 0 || continue
        nodes = cell_nodes[cell]
        edges[-dof] = minmax(Int(nodes[a]),Int(nodes[b]))
    end
    edges
end

# Unstructured mesh of the unit square or cube with element size h. With Mesh.MeshSizeMax
# alone the mesh was the same for every h tried (0.3 to 0.45 in 2D, 0.4 to 2 in 3D), finer
# than all of them, so the minimum is set too and the mesh follows h. In 3D, h = 0.25 and
# 0.2 give tetrahedra in all 24 vertex orders; in 2D, the triangles came in 2 of the 6
# vertex orders at every h tried (0.1 to 0.45).
function gmsh_simplex_mesh(D,h)
    GT.with_gmsh() do gmsh
        gmsh.option.setNumber("General.Verbosity",0)
        gmsh.option.setNumber("Mesh.MeshSizeMin",h)
        gmsh.option.setNumber("Mesh.MeshSizeMax",h)
        if D == 3
            gmsh.model.occ.addBox(0,0,0,1,1,1)
        else
            gmsh.model.occ.addRectangle(0,0,0,1,1)
        end
        gmsh.model.occ.synchronize()
        gmsh.model.mesh.generate(D)
        GT.mesh_from_gmsh(gmsh)
    end
end

meshes = [
    (2,GT.cartesian_mesh((0,1,0,1),(3,2);simplexify=true)),
    (2,gmsh_simplex_mesh(2,0.12)),
    (3,GT.cartesian_mesh((0,2,0,1,0,0.5),(2,2,2);simplexify=true)),
    (3,gmsh_simplex_mesh(3,0.2)),
]

for (D,mesh) in meshes
    Ω = GT.interior(mesh)
    V = GT.nedelec_space(Ω,0)
    nedges = GT.num_faces(GT.topology(mesh),1)
    @test GT.num_free_dofs(V) == nedges
    @test GT.num_dirichlet_dofs(V) == 0
    orientation = GT.nedelec_cell_orientation(V)
    @test length(orientation) == GT.num_faces(mesh,D)
    @test all(o -> 1 <= o <= factorial(D+1),orientation)
    dof_to_edge = dof_to_edge_map(V,mesh,D)
    M_ref,K_ref,G_ref = whitney_matrices(mesh,D,dof_to_edge)

    dΩ = GT.measure(Ω,2)
    m(u,v) = GT.∫(x->u(x)⋅v(x),dΩ)
    k(u,v) = GT.∫(x->GT.call(curl,ForwardDiff.jacobian(u,x))⋅GT.call(curl,ForwardDiff.jacobian(v,x)),dΩ)
    g(u,v) = GT.∫(x->GT.call(frobenius,ForwardDiff.jacobian(u,x),ForwardDiff.jacobian(v,x)),dΩ)
    csc = (;matrix=GT.element_csc_assembly())
    for (execution,assembly_method) in ((GT.SerialExecution(),nothing),(GT.SerialExecution(),csc),(GT.CPUThreadedExecution(),csc))
        kwargs = assembly_method === nothing ? (;execution) : (;execution,assembly_method)
        M = GT.assemble_matrix(m,Float64,V,V;kwargs...)
        K = GT.assemble_matrix(k,Float64,V,V;kwargs...)
        @test norm(Matrix(M)-M_ref,Inf) <= 1e-12*norm(M_ref,Inf)
        @test norm(Matrix(K)-K_ref,Inf) <= 1e-12*norm(K_ref,Inf)
    end
    Gm = GT.assemble_matrix(g,Float64,V,V)
    @test norm(Matrix(Gm)-G_ref,Inf) <= 1e-12*norm(G_ref,Inf)

    # Tangential Dirichlet: the closure edges of the boundary faces.
    Γ = GT.boundary(mesh)
    VD = GT.nedelec_space(Ω,0;dirichlet_boundary=Γ)
    boundary_edges = Set{Tuple{Int,Int}}()
    for nodes in GT.face_nodes(mesh,D-1)[GT.faces(Γ)], i in 1:length(nodes), j in i+1:length(nodes)
        push!(boundary_edges,minmax(Int(nodes[i]),Int(nodes[j])))
    end
    @test GT.num_dirichlet_dofs(VD) == length(boundary_edges)
    @test GT.num_free_dofs(VD) == nedges - length(boundary_edges)
    MD = GT.assemble_matrix(m,Float64,VD,VD)
    @test size(MD) == (GT.num_free_dofs(VD),GT.num_free_dofs(VD))
    @test isposdef(Symmetric(Matrix(MD)))

    # A selected part of the boundary: exactly the edges in the closure of its faces become Dirichlet DOFs.
    group = D == 3 ? "2-face-1" : "1-face-1"
    if haskey(GT.group_faces(mesh,D-1),group)
        Γ1 = GT.boundary(mesh;group_names=[group])
        V1 = GT.nedelec_space(Ω,0;dirichlet_boundary=Γ1)
        expected = Set{Tuple{Int,Int}}()
        for nodes in GT.face_nodes(mesh,D-1)[GT.faces(Γ1)], i in 1:length(nodes), j in i+1:length(nodes)
            push!(expected,minmax(Int(nodes[i]),Int(nodes[j])))
        end
        @test !isempty(expected)
        @test Set(dirichlet_edges(V1,mesh,D)) == expected
        @test GT.num_dirichlet_dofs(V1) == length(expected)
    end
end

# Exactness: a field in NED0 (constant + rotation) is interpolated exactly, so the
# discrete field equals it at every quadrature point, with the same curl.

for (D,mesh) in meshes
    Ω = GT.interior(mesh)
    V = GT.nedelec_space(Ω,0)
    dΩ = GT.measure(Ω,2)
    if D == 3
        a = SVector(0.3,-1.2,0.7); b = SVector(1.1,0.4,-0.5)
        ufun = x -> a + cross(b,x)
        curl_u = 2*b
        jac_u = SMatrix{3,3}(0,b[3],-b[2],-b[3],0,b[1],b[2],-b[1],0)
    else
        a = SVector(0.3,-1.2); b = 0.8
        ufun = x -> a + b*SVector(-x[2],x[1])
        curl_u = 2*b
        jac_u = SMatrix{2,2}(0,b,-b,0)
    end
    u = GT.analytical_field(ufun,Ω)
    uh = GT.zero_field(Float64,V)
    GT.interpolate!(u,uh)
    @test GT.assemble_scalar(GT.∫(x->norm(u(x)-uh(x))^2,dΩ)) < 1e-24
    curl_error = c -> norm(c-curl_u)^2
    @test GT.assemble_scalar(GT.∫(x->GT.call(curl_error,GT.call(curl,ForwardDiff.jacobian(uh,x))),dΩ)) < 1e-22
    uh2 = GT.interpolate(u,V)
    @test GT.free_values(uh2) ≈ GT.free_values(uh)

    # The same field through the space/field accessors (values and Jacobians at every quadrature point).
    uh_dΩ = GT.field_accessor(uh,dΩ;tabulate=(GT.value,ForwardDiff.jacobian),compute=(GT.coordinate,))
    worst_value,worst_jac = 0.0,0.0
    for face in GT.each_face(uh_dΩ), point in GT.each_point(face)
        x = GT.coordinate(point)
        worst_value = max(worst_value,norm(GT.field(GT.value,point)-ufun(x)))
        worst_jac = max(worst_jac,norm(GT.field(ForwardDiff.jacobian,point)-jac_u))
    end
    @test worst_value < 1e-12
    @test worst_jac < 1e-11

    # Untabulated shape functions evaluated at physical points (cell centroids).
    s_value = GT.shape_function_accessor(GT.value,V)
    s_jac = GT.shape_function_accessor(ForwardDiff.jacobian,V)
    cell_nodes = GT.face_nodes(mesh,D)
    node_x = GT.node_coordinates(mesh)
    for cell in (1,length(cell_nodes))
        xc = sum(node_x[cell_nodes[cell]])/(D+1)
        dofs = GT.face_dofs(V)[cell]
        vals = GT.free_values(uh)[dofs]
        @test sum(vals[i]*s_value(cell)(i)(xc) for i in eachindex(dofs)) ≈ ufun(xc)
        @test sum(vals[i]*s_jac(cell)(i)(xc) for i in eachindex(dofs)) ≈ jac_u atol=1e-10
    end

    # Complex fields: the moments are complex-linear, so u + i u is reproduced exactly.
    uc = GT.analytical_field(x->(1+2im)*ufun(x),Ω)
    uhc = GT.zero_field(ComplexF64,V)
    GT.interpolate!(uc,uhc)
    @test GT.free_values(uhc) ≈ (1+2im)*GT.free_values(uh)
    @test GT.assemble_scalar(GT.∫(x->norm(uc(x)-uhc(x))^2,dΩ)) < 1e-24

    # The same field through vector assembly: M x = b with b_i = ∫u⋅φi.
    m(u,v) = GT.∫(x->u(x)⋅v(x),dΩ)
    l(v) = GT.∫(x->u(x)⋅v(x),dΩ)
    M = GT.assemble_matrix(m,Float64,V,V)
    rhs = GT.assemble_vector(l,Float64,V)
    @test M\rhs ≈ GT.free_values(uh)

    # Dirichlet values come from the same edge moments.
    VD = GT.nedelec_space(Ω,0;dirichlet_boundary=GT.boundary(mesh))
    uhD = GT.zero_field(Float64,VD)
    GT.interpolate_dirichlet!(u,uhD)
    GT.interpolate_free!(u,uhD)
    @test GT.assemble_scalar(GT.∫(x->norm(u(x)-uhD(x))^2,dΩ)) < 1e-24
end

# Convergence: curl curl u + u = f with tangential Dirichlet data, manufactured solution.

function solve_maxwell(D,n;k=0)
    domain = D == 3 ? (0,1,0,1,0,1) : (0,1,0,1)
    mesh = GT.cartesian_mesh(domain,ntuple(_->n,D);simplexify=true)
    Ω = GT.interior(mesh)
    Γ = GT.boundary(mesh)
    V = GT.nedelec_space(Ω,k;dirichlet_boundary=Γ)
    dΩ = GT.measure(Ω,2k+4)
    if D == 3
        ufun = x -> SVector(sinpi(x[2])*sinpi(x[3]),sinpi(x[1])*sinpi(x[3]),sinpi(x[1])*sinpi(x[2]))
        curlfun = x -> curl(ForwardDiff.jacobian(ufun,x))
        ffun = x -> (2*π^2+1)*ufun(x)
    else
        ufun = x -> SVector(sinpi(x[2]),sinpi(x[1]))
        curlfun = x -> curl(ForwardDiff.jacobian(ufun,x))
        ffun = x -> (π^2+1)*ufun(x)
    end
    u = GT.analytical_field(ufun,Ω)
    f = GT.analytical_field(ffun,Ω)
    curl_u = GT.analytical_field(curlfun,Ω)
    a(u,v) = GT.∫(x->GT.call(curl,ForwardDiff.jacobian(u,x))⋅GT.call(curl,ForwardDiff.jacobian(v,x)) + u(x)⋅v(x),dΩ)
    l(v) = GT.∫(x->f(x)⋅v(x),dΩ)
    A = GT.assemble_matrix(a,Float64,V,V)
    b = GT.assemble_vector(l,Float64,V)
    uh = GT.zero_field(Float64,V)
    GT.free_values(uh) .= A\b
    el2 = sqrt(GT.assemble_scalar(GT.∫(x->norm(u(x)-uh(x))^2,dΩ)))
    ecurl = sqrt(GT.assemble_scalar(GT.∫(x->norm(curl_u(x)-GT.call(curl,ForwardDiff.jacobian(uh,x)))^2,dΩ)))
    el2,ecurl
end

for D in (2,3)
    ns = D == 3 ? (2,4,8) : (4,8,16)
    errors = [solve_maxwell(D,n) for n in ns]
    rates_l2 = [log2(errors[i][1]/errors[i+1][1]) for i in 1:length(ns)-1]
    rates_curl = [log2(errors[i][2]/errors[i+1][2]) for i in 1:length(ns)-1]
    @test errors[end][1] < errors[1][1]
    @test errors[end][2] < errors[1][2]
    @test rates_l2[end] > 0.85
    @test rates_curl[end] > 0.85
end

# Piecewise Dirichlet data: each boundary piece receives the moments of its own field.
let
    mesh = GT.cartesian_mesh((0,1,0,1),(3,3);simplexify=true)
    Ω = GT.interior(mesh)
    pieces = [GT.boundary(mesh;group_names=[name]) for name in ("1-face-1","1-face-2")]
    V = GT.nedelec_space(Ω,0;dirichlet_boundary=GT.piecewise_domain(pieces...))
    cs = (SVector(1.0,1.0),SVector(3.0,-2.0))
    u = GT.piecewise_field(map(c->GT.analytical_field(x->c,Ω),cs)...)
    uh = GT.zero_field(Float64,V)
    GT.interpolate_dirichlet!(u,uh)
    locations = GT.dirichlet_dof_location(V)
    @test sort(unique(locations)) == [1,2]
    node_x = GT.node_coordinates(mesh)
    edges = dirichlet_edges(V,mesh,2)
    expected = map(edges,locations) do (na,nb),location
        dot(cs[location],node_x[nb]-node_x[na])
    end
    @test GT.dirichlet_values(uh) ≈ expected

    # An ordinary (non-piecewise) field fills every Dirichlet piece.
    c = SVector(2.0,3.0)
    uh_all = GT.zero_field(Float64,V)
    GT.dirichlet_values(uh_all) .= 123
    GT.interpolate_dirichlet!(GT.analytical_field(x->c,Ω),uh_all)
    @test GT.dirichlet_values(uh_all) ≈ [dot(c,node_x[nb]-node_x[na]) for (na,nb) in edges]
end

# 3D pieces that share edges: a shared edge belongs to the later piece. The expected
# locations come from the face groups directly, not from dirichlet_dof_location.
let
    mesh = GT.cartesian_mesh((0,1,0,1,0,1),(2,2,2);simplexify=true)
    Ω = GT.interior(mesh)
    names = ("2-face-1","2-face-3")
    pieces = [GT.boundary(mesh;group_names=[name]) for name in names]
    piece_edges = map(pieces) do Γ
        edges = Set{Tuple{Int,Int}}()
        for nodes in GT.face_nodes(mesh,2)[GT.faces(Γ)], i in 1:3, j in i+1:3
            push!(edges,minmax(Int(nodes[i]),Int(nodes[j])))
        end
        edges
    end
    shared = intersect(piece_edges...)
    @test !isempty(shared)
    V = GT.nedelec_space(Ω,0;dirichlet_boundary=GT.piecewise_domain(pieces...))
    edges = dirichlet_edges(V,mesh,3)
    @test Set(edges) == union(piece_edges...)
    expected_location = [edge in piece_edges[2] ? 2 : 1 for edge in edges]
    @test GT.dirichlet_dof_location(V) == expected_location
    cs = (SVector(1.0,-1.0,0.5),SVector(-2.0,0.25,3.0))
    uh = GT.zero_field(Float64,V)
    GT.interpolate_dirichlet!(GT.piecewise_field(map(c->GT.analytical_field(x->c,Ω),cs)...),uh)
    node_x = GT.node_coordinates(mesh)
    @test GT.dirichlet_values(uh) ≈ [dot(cs[l],node_x[nb]-node_x[na]) for ((na,nb),l) in zip(edges,expected_location)]
    c = SVector(0.5,1.5,-1.0)
    GT.interpolate_dirichlet!(GT.analytical_field(x->c,Ω),uh)
    @test GT.dirichlet_values(uh) ≈ [dot(c,node_x[nb]-node_x[na]) for (na,nb) in edges]
end

# The compiler's untabulated (general-operator) helper and the untabulated accessor
# against the closed-form Whitney functions λa∇λb - λb∇λa, edge (a,b) with ascending
# global node IDs.
let
    mesh = GT.cartesian_mesh((0,2,0,1,0,1),(2,2,2);simplexify=true)
    Ω = GT.interior(mesh)
    V = GT.nedelec_space(Ω,0)
    fe = only(GT.reference_spaces(V))
    geometry = only(GT.reference_spaces(mesh,3))
    node_x = GT.node_coordinates(mesh)
    edge_nodes = GT.face_nodes(GT.mesh(GT.domain(fe)),1)
    s_value = GT.shape_function_accessor(GT.value,V)
    s_jac = GT.shape_function_accessor(ForwardDiff.jacobian,V)
    for cell in (1,7,GT.num_faces(Ω))
        nodes = GT.face_nodes(mesh,3)[cell]
        x = node_x[nodes]
        B = hcat([x[k+1]-x[1] for k in 1:3]...)
        Binv = inv(B)
        ∇λ = vcat([-vec(sum(Binv,dims=1))],[Binv[k,:] for k in 1:3])
        λ(y) = (l = Binv*(y-x[1]); vcat(1-sum(l),l))
        xc = sum(x)/4 + SVector(0.01,0.0,-0.01)
        P = Int(GT.nedelec_cell_orientation(V)[cell])
        for dof in 1:6
            a,b = edge_nodes[dof]
            a,b = nodes[a] < nodes[b] ? (a,b) : (b,a)
            whitney(y) = (l = λ(y); l[a]*∇λ[b] - l[b]*∇λ[a])
            s = GT.physical_covariant_shape_function(GT.nedelec_oriented_shape_function(fe,dof,P),
                GT.shape_functions(geometry),node_x,nodes,zero(SVector{3,Float64}))
            @test s(xc) ≈ whitney(xc)
            @test ForwardDiff.jacobian(s,xc) ≈ ∇λ[b]*∇λ[a]' - ∇λ[a]*∇λ[b]'
            @test s_value(cell)(dof)(xc) ≈ whitney(xc)
            @test s_jac(cell)(dof)(xc) ≈ ∇λ[b]*∇λ[a]' - ∇λ[a]*∇λ[b]'
        end
    end
end

# Higher orders.

nedelec_dim(D,k) = D == 2 ? (k+1)*(k+3) : div((k+1)*(k+3)*(k+4),2)
nedelec_interior_dim(D,k) = D == 2 ? k*(k+1) : div((k-1)*k*(k+1),2)

# Reference element of order k: dimension, duality, ownership, and the orientation
# matrices. For every vertex permutation P, the oriented function of a slot must be the
# covariant push-forward, by the affine map R of the reference simplex with
# R(v̂ᵢ) = v̂_{P[i]}, of the reference function on the entity that R maps to the slot's
# entity (same own index). This is recomputed here from the entity vertex sets.
for (D,ks) in ((2,0:3),(3,0:2)), k in ks
    geo = GT.unit_simplex(Val(D))
    fe = GT.nedelec_space(geo,k)
    n = GT.num_dofs(fe)
    @test n == nedelec_dim(D,k)
    shapes = GT.shape_functions(fe)
    sigma = GT.dual_basis(fe)
    @test [sigma[i](shapes[j]) for i in 1:n, j in 1:n] ≈ Matrix(I,n,n) atol=1e-10
    @test [sigma[i](x->im*shapes[j](x)) for i in 1:n, j in 1:n] ≈ im*Matrix(I,n,n) atol=1e-10
    own = [GT.face_own_dofs(fe,d) for d in 0:D]
    @test all(isempty,own[1])
    @test all(dofs->length(dofs)==k+1,own[2])
    if D == 3
        @test all(dofs->length(dofs)==k*(k+1),own[3])
        @test all(dofs->length(dofs)==3*(k+1)+k*(k+1),GT.face_dofs(fe,2))
    end
    @test length(only(own[D+1])) == nedelec_interior_dim(D,k)
    @test sort(reduce(vcat,[reduce(vcat,dofs;init=Int32[]) for dofs in own])) == 1:n
    @test only(GT.face_dofs(fe,D)) == 1:n
    cell_mesh = GT.mesh(geo)
    vertices = [SVector{D,Float64}(x) for x in GT.node_coordinates(cell_mesh)]
    entity_nodes = [map(nodes->sort(Int.(collect(nodes))),GT.face_nodes(cell_mesh,d)) for d in 0:D]
    perms = fe.workspace.orientation_permutations
    @test length(perms) == factorial(D+1) && allunique(perms)
    ξs = [SVector{D}(ntuple(i->0.13*i,D)),SVector{D}(ntuple(i->0.3-0.07i,D))]
    worst = 0.0
    for (P,perm) in enumerate(perms)
        inverse = invperm(perm)
        b = vertices[perm[1]]
        Ainv = inv(hcat([vertices[perm[i+1]]-b for i in 1:D]...))
        for d in 1:D, (entity,dofs) in enumerate(own[d+1]), (i,slot) in enumerate(dofs)
            mapped = sort([inverse[v] for v in entity_nodes[d+1][entity]])
            source = own[d+1][findfirst(==(mapped),entity_nodes[d+1])][i]
            ψ = GT.nedelec_oriented_shape_function(fe,slot,P)
            for ξ in ξs
                expected = transpose(Ainv)*shapes[source](Ainv*(ξ-b))
                worst = max(worst,norm(ψ(ξ)-expected))
            end
        end
    end
    @test worst < 1e-9
    if k == 0
        C = fe.workspace.orientation_matrices
        @test all(P -> isdiag(C[:,:,P]) && all(c -> abs(abs(c)-1) < 1e-12,diag(C[:,:,P])),1:length(perms))
    end
end

function sheared_simplex_mesh(D,n)
    mesh = GT.cartesian_mesh(D == 3 ? (0,1,0,1,0,1) : (0,1,0,1),ntuple(_->n,D);simplexify=true)
    GT.replace_node_coordinates(mesh,
        [SVector{D}(ntuple(i->1.2x[i]+0.1(sum(x)-x[i]),D)) for x in GT.node_coordinates(mesh)])
end

# Tangential continuity of a discrete field with arbitrary DOF values, evaluated from
# both cells at points of shared faces (edges in 2D); normal components may jump.
function tangential_jumps(V,mesh,D;maxfaces=20)
    node_x = GT.node_coordinates(mesh)
    face_cells = GT.face_incidence(GT.topology(mesh),D-1,D)
    face_nodes = GT.face_nodes(mesh,D-1)
    cell_dofs = GT.face_dofs(V)
    vals = [sin(1.7i)+0.3 for i in 1:GT.num_free_dofs(V)]
    s = GT.shape_function_accessor(GT.value,V)
    uh(cell,x) = sum(vals[dof]*s(cell)(i)(x) for (i,dof) in enumerate(cell_dofs[cell]))
    weights = D == 3 ? ((1/3,1/3,1/3),(0.6,0.3,0.1),(0.1,0.2,0.7)) : ((0.5,0.5),(0.2,0.8),(0.9,0.1))
    tangential,total,nchecked = 0.0,0.0,0
    for (face,cells) in enumerate(face_cells)
        length(cells) == 2 || continue
        nchecked == maxfaces && break
        nchecked += 1
        xs = node_x[face_nodes[face]]
        tangents = [xs[i+1]-xs[1] for i in 1:D-1]
        for w in weights
            x = sum(w[i]*xs[i] for i in 1:D)
            jump = uh(cells[1],x)-uh(cells[2],x)
            tangential = max(tangential,maximum(t->abs(dot(t,jump)),tangents))
            total = max(total,norm(jump))
        end
    end
    tangential,total,nchecked
end

meshes_k = [(2,sheared_simplex_mesh(2,2)),(2,gmsh_simplex_mesh(2,0.1)),
            (3,sheared_simplex_mesh(3,2)),(3,gmsh_simplex_mesh(3,0.25))]

for (D,mesh) in meshes_k, k in (1,2)
    Ω = GT.interior(mesh)
    V = GT.nedelec_space(Ω,k)
    topo = GT.topology(mesh)
    nedges,nfaces,ncells = GT.num_faces(topo,1),GT.num_faces(topo,2),GT.num_faces(topo,D)
    expected = nedges*(k+1) + (D == 3 ? nfaces*k*(k+1) : 0) + ncells*nedelec_interior_dim(D,k)
    @test GT.num_free_dofs(V) == expected

    tangential,total,nchecked = tangential_jumps(V,mesh,D)
    @test nchecked > 0
    @test tangential < 1e-10*max(1,total)
    @test total > 1e-3 # normal components do jump: the check is not vacuous

    # A field of NED_k, including its top-degree part S_{k+1} (q(x)(-y,x) in 2D,
    # x × (q(x) e₃) in 3D, q homogeneous of degree k), is reproduced exactly.
    a = SVector{D}(ntuple(i->0.3i-0.5,D))
    Bm = SMatrix{D,D}(ntuple(i->sin(1.3i),D*D))
    quadratic = x -> SVector{D}(ntuple(c->x[c]*x[mod1(c+1,D)]-0.5*x[1]^2,D))
    top = D == 2 ? (x -> x[1]^k*SVector(-x[2],x[1])) : (x -> cross(x,SVector(0.0,0.0,x[1]^(k-1)*x[2])))
    ufun = k == 1 ? (x -> a + Bm*x + top(x)) : (x -> a + Bm*x + quadratic(x) + 0.7*top(x))
    u = GT.analytical_field(ufun,Ω)
    curl_u = GT.analytical_field(x->curl(ForwardDiff.jacobian(ufun,x)),Ω)
    uh = GT.zero_field(Float64,V)
    GT.interpolate!(u,uh)
    dΩ = GT.measure(Ω,2k+2)
    @test GT.assemble_scalar(GT.∫(x->norm(u(x)-uh(x))^2,dΩ)) < 1e-22
    @test GT.assemble_scalar(GT.∫(x->norm(curl_u(x)-GT.call(curl,ForwardDiff.jacobian(uh,x)))^2,dΩ)) < 1e-18
    uh_dΩ = GT.field_accessor(uh,dΩ;tabulate=(GT.value,ForwardDiff.jacobian),compute=(GT.coordinate,))
    worst_value,worst_jac = 0.0,0.0
    for face in GT.each_face(uh_dΩ), point in GT.each_point(face)
        x = GT.coordinate(point)
        worst_value = max(worst_value,norm(GT.field(GT.value,point)-ufun(x)))
        worst_jac = max(worst_jac,norm(GT.field(ForwardDiff.jacobian,point)-ForwardDiff.jacobian(ufun,x)))
    end
    @test worst_value < 1e-10
    @test worst_jac < 1e-9
    m(u,v) = GT.∫(x->u(x)⋅v(x),dΩ)
    l(v) = GT.∫(x->u(x)⋅v(x),dΩ)
    M = GT.assemble_matrix(m,Float64,V,V)
    @test M\GT.assemble_vector(l,Float64,V) ≈ GT.free_values(uh)
    csc = (;matrix=GT.element_csc_assembly())
    @test GT.assemble_matrix(m,Float64,V,V;execution=GT.CPUThreadedExecution(),assembly_method=csc) ≈ M

    # Tangential Dirichlet: the DOFs of the closure of the boundary faces.
    Γ = GT.boundary(mesh)
    VD = GT.nedelec_space(Ω,k;dirichlet_boundary=Γ)
    boundary_faces = GT.face_nodes(mesh,D-1)[GT.faces(Γ)]
    boundary_edges = Set{Tuple{Int,Int}}()
    for nodes in boundary_faces, i in 1:D, j in i+1:D
        push!(boundary_edges,minmax(Int(nodes[i]),Int(nodes[j])))
    end
    expected_diri = length(boundary_edges)*(k+1) + (D == 3 ? length(boundary_faces)*k*(k+1) : 0)
    @test GT.num_dirichlet_dofs(VD) == expected_diri
    uhD = GT.zero_field(Float64,VD)
    GT.interpolate!(u,uhD)
    @test GT.assemble_scalar(GT.∫(x->norm(u(x)-uhD(x))^2,dΩ)) < 1e-22
end

# Discrete Maxwell eigenvalues with tangential Dirichlet walls (unit square and cube):
# the kernel of the curl has the dimension of the gradients of the Lagrange space of
# order k+1, and the first nonzero eigenvalues approach π²(m²+n²(+p²)).
function maxwell_eigenvalues(D,n,k)
    mesh = GT.cartesian_mesh(D == 3 ? (0,1,0,1,0,1) : (0,1,0,1),ntuple(_->n,D);simplexify=true)
    Ω = GT.interior(mesh)
    Γ = GT.boundary(mesh)
    V = GT.nedelec_space(Ω,k;dirichlet_boundary=Γ)
    dΩ = GT.measure(Ω,2k+2)
    M = GT.assemble_matrix((u,v)->GT.∫(x->u(x)⋅v(x),dΩ),Float64,V,V)
    K = GT.assemble_matrix((u,v)->GT.∫(x->GT.call(curl,ForwardDiff.jacobian(u,x))⋅GT.call(curl,ForwardDiff.jacobian(v,x)),dΩ),Float64,V,V)
    λ = eigvals(Symmetric(Matrix(K)),Symmetric(Matrix(M)))
    nkernel = GT.num_free_dofs(GT.lagrange_space(Ω,k+1;dirichlet_boundary=Γ))
    λ,nkernel
end

for (D,n,k,exact) in ((2,4,1,π^2*[1,1,2,4,4]),(2,4,2,π^2*[1,1,2,4,4]),(3,2,1,π^2*[2,2,2,3,3]))
    λ,nkernel = maxwell_eigenvalues(D,n,k)
    zero_tol = 1e-8*maximum(λ)
    @test count(<(zero_tol),λ) == nkernel
    physical = filter(>=(zero_tol),λ)[1:length(exact)]
    @test maximum(abs.(physical-exact)./exact) < (D == 3 ? 0.1 : 0.02)
end

# Convergence at rate k+1 for curl curl u + u = f with tangential Dirichlet data.
for (D,k,ns) in ((2,1,(4,8,16)),(2,2,(4,8,16)),(3,1,(2,4,8)))
    errors = [solve_maxwell(D,n;k) for n in ns]
    rates_l2 = [log2(errors[i][1]/errors[i+1][1]) for i in 1:length(ns)-1]
    rates_curl = [log2(errors[i][2]/errors[i+1][2]) for i in 1:length(ns)-1]
    @test rates_l2[end] > k+1-0.15
    @test rates_curl[end] > k+1-0.15
end

# Piecewise Dirichlet data at order 1: each Dirichlet DOF gets the value of its own
# piece's field (compared with interpolating that field alone).
let
    mesh = GT.cartesian_mesh((0,1,0,1,0,1),(2,2,2);simplexify=true)
    Ω = GT.interior(mesh)
    pieces = [GT.boundary(mesh;group_names=[name]) for name in ("2-face-1","2-face-3")]
    V = GT.nedelec_space(Ω,1;dirichlet_boundary=GT.piecewise_domain(pieces...))
    fields = (x->SVector(1.0+x[2],-x[1],0.5),x->SVector(-2.0,0.25*x[3],3.0+x[1]))
    uh = GT.zero_field(Float64,V)
    GT.interpolate_dirichlet!(GT.piecewise_field(map(f->GT.analytical_field(f,Ω),fields)...),uh)
    alone = map(fields) do f
        uhf = GT.zero_field(Float64,V)
        GT.interpolate_dirichlet!(GT.analytical_field(f,Ω),uhf)
        GT.dirichlet_values(uhf)
    end
    locations = GT.dirichlet_dof_location(V)
    @test sort(unique(locations)) == [1,2]
    @test GT.dirichlet_values(uh) ≈ [alone[l][i] for (i,l) in enumerate(locations)]
end

# Traces on faces of dimension D-1. Fields of the space are reproduced exactly on
# the boundary, and face matrices and vectors agree with integrals of the interpolant.
for D in (2,3), k in (0,1,2)
    mesh = GT.cartesian_mesh(D == 3 ? (0,1,0,1,0,1) : (0,1,0,1),ntuple(_->2,D);simplexify=true)
    Ω = GT.interior(mesh)
    V = GT.nedelec_space(Ω,k)
    a = SVector{D}(ntuple(i->0.3i-0.5,D))
    B = SMatrix{D,D}(ntuple(i->sin(1.3i),D*D))
    # x ↦ x×q(x) (3D) or q(x)(-x2,x1) (2D), with q homogeneous of degree k-1 or k: in NED_k.
    rot = D == 3 ? (x->cross(x,SVector(0.0,0.0,x[1]^max(k-1,0)*x[2]))) : (x->x[1]^k*SVector(-x[2],x[1]))
    quad = x->SVector{D}(ntuple(c->x[c]*x[mod1(c+1,D)]-0.5*x[1]^2,D))
    f = k == 0 ? (D == 3 ? (x->SVector(0.3,-1.2,0.7)+cross(SVector(1.1,0.4,-0.5),x)) : (x->SVector(0.3,-1.2)+0.8*SVector(-x[2],x[1]))) :
        k == 1 ? (x->a+B*x+rot(x)) : (x->a+B*x+quad(x)+0.7*rot(x))
    u = GT.analytical_field(f,Ω)
    uh = GT.interpolate(u,V)
    x = GT.free_values(uh)
    dΓ = GT.measure(GT.boundary(mesh),2k+2)
    uu = GT.assemble_scalar(GT.∫(q->u(q)⋅u(q),dΓ))
    @test GT.assemble_scalar(GT.∫(q->norm(u(q)-uh(q))^2,dΓ)) < 1e-20*uu
    MΓ = GT.assemble_matrix((u,v)->GT.∫(q->u(q)⋅v(q),dΓ),Float64,V,V)
    @test x'*MΓ*x ≈ uu rtol=1e-12
    bΓ = GT.assemble_vector(v->GT.∫(q->u(q)⋅v(q),dΓ),Float64,V)
    @test bΓ'*x ≈ uu rtol=1e-12
end

# Curl on boundary faces: for u = a + b×x, curl u = 2b, so ∫|curl u|² over the boundary
# of the unit cube is 4|b|²·6.
let
    mesh = GT.cartesian_mesh((0,1,0,1,0,1),(2,2,2);simplexify=true)
    Ω = GT.interior(mesh)
    dΓ = GT.measure(GT.boundary(mesh),2)
    b = SVector(1.1,0.4,-0.5)
    expected = 4*dot(b,b)*6
    for k in (0,1)
        V = GT.nedelec_space(Ω,k)
        uh = GT.interpolate(GT.analytical_field(x->SVector(0.3,-1.2,0.7)+cross(b,x),Ω),V)
        cc = GT.assemble_scalar(GT.∫(q->GT.call(curl,ForwardDiff.jacobian(uh,q))⋅GT.call(curl,ForwardDiff.jacobian(uh,q)),dΓ))
        @test cc ≈ expected rtol=1e-12
        K = GT.assemble_matrix((u,v)->GT.∫(q->GT.call(curl,ForwardDiff.jacobian(u,q))⋅GT.call(curl,ForwardDiff.jacobian(v,q)),dΓ),Float64,V,V)
        x = GT.free_values(uh)
        @test x'*K*x ≈ expected rtol=1e-12
    end
end

# Interior faces. A face group inside the mesh is integrated from one incident cell,
# chosen by faces_around (default: the first incident cell), or as the boundary of a
# subdomain. The field below is (c,1) below the plane x_D = 1/2 and (c,2) above it: it
# lies in the space (its tangential trace on the plane is continuous), and its normal
# component on the plane depends on the side.
for D in (2,3)
    mesh = GT.cartesian_mesh(D == 3 ? (0,1,0,1,0,1) : (0,1,0,1),ntuple(_->2,D);simplexify=true)
    xn = GT.node_coordinates(mesh)
    fnodes = GT.face_nodes(mesh,D-1)
    cnodes = GT.face_nodes(mesh,D)
    height(nodes) = sum(i->xn[i][D],nodes)/length(nodes)
    GT.group_faces(mesh,D-1)["mid"] = findall(nodes->all(i->abs(xn[i][D]-0.5)<1e-12,nodes),fnodes)
    GT.group_faces(mesh,D)["low"] = findall(nodes->height(nodes)<0.5,cnodes)
    Ω = GT.interior(mesh)
    face_to_cells = GT.face_incidence(GT.topology(mesh),D-1,D)
    Γ = GT.boundary(mesh;group_names=["mid"])
    mid = GT.faces(Γ)
    below = [findfirst(c->height(cnodes[c])<0.5,face_to_cells[f]) for f in mid]
    above = [findfirst(c->height(cnodes[c])>0.5,face_to_cells[f]) for f in mid]
    @test all(f->length(face_to_cells[f])==2,mid)
    sides = (;below=GT.boundary(mesh;group_names=["mid"],faces_around=below),
              above=GT.boundary(mesh;group_names=["mid"],faces_around=above),
              low=GT.boundary(GT.interior(mesh;group_names=["low"]);group_names=["mid"]),
              dflt=Γ,
              first=GT.boundary(mesh;group_names=["mid"],faces_around=fill(1,length(mid))),
              mixed=GT.boundary(mesh;group_names=["mid"],
                  faces_around=[isodd(i) ? below[i] : above[i] for i in eachindex(mid)]))
    @test collect(GT.faces_around(sides.low)) == below
    c = D == 3 ? SVector(0.3,0.2) : SVector(0.3)
    e = SVector{D}(ntuple(i->Float64(i==D),D))
    t = SVector{D}(ntuple(i->i<D ? 1.0+i : 0.0,D))
    T = SMatrix{D,D}(ntuple(i->(1.0+0.1i)+0.3im*sin(i),D*D))
    for k in (0,1)
        V = GT.nedelec_space(Ω,k)
        uh = GT.interpolate(GT.analytical_field(x->SVector(c...,x[D]<0.5 ? 1.0 : 2.0),Ω),V)
        x = GT.free_values(uh)
        values = map(sides) do Γs
            dΓs = GT.measure(Γs,2k+2)
            (normal=GT.assemble_scalar(GT.∫(q->uh(q)[D]^2,dΓs)),
             mass=x'*GT.assemble_matrix((u,v)->GT.∫(q->u(q)⋅v(q),dΓs),Float64,V,V)*x,
             load_n=GT.assemble_vector(v->GT.∫(q->v(q)⋅GT.uniform_quantity(e),dΓs),Float64,V),
             load_t=GT.assemble_vector(v->GT.∫(q->v(q)⋅GT.uniform_quantity(t),dΓs),Float64,V),
             tensor=transpose(x)*GT.assemble_matrix((u,v)->GT.∫(q->im*(v(q)⋅(GT.uniform_quantity(T)*u(q))),dΓs),ComplexF64,V,V)*x)
        end
        @test values.below.normal ≈ 1 rtol=1e-12
        @test values.above.normal ≈ 4 rtol=1e-12
        @test values.low.normal ≈ 1 rtol=1e-12
        @test values.below.mass ≈ dot(c,c)+1 rtol=1e-12
        @test values.above.mass ≈ dot(c,c)+4 rtol=1e-12
        @test values.below.load_n'*x ≈ 1 rtol=1e-12
        @test values.above.load_n'*x ≈ 2 rtol=1e-12
        @test values.low.load_n ≈ values.below.load_n
        # A tangential load is the same from either side.
        @test values.below.load_t ≈ values.above.load_t rtol=1e-12
        @test values.below.load_t'*x ≈ dot(t[1:D-1],c) rtol=1e-12
        ub = SVector(c...,1.0); ua = SVector(c...,2.0)
        @test values.below.tensor ≈ im*dot(ub,T*ub) rtol=1e-12
        @test values.above.tensor ≈ im*dot(ua,T*ua) rtol=1e-12
        # The default side is the first incident cell of each face. The faces on the
        # plane have equal measure, so the normal integral counts the faces whose first
        # incident cell is below (1) or above (4).
        nbelow = count(==(1),below)
        @test values.dflt.normal ≈ values.first.normal rtol=1e-12
        @test values.dflt.normal ≈ (nbelow+4*(length(mid)-nbelow))/length(mid) rtol=1e-12
        @test values.dflt.mass ≈ values.first.mass rtol=1e-12
        @test values.dflt.load_n ≈ values.first.load_n
        @test values.dflt.load_t ≈ values.below.load_t rtol=1e-12
        @test values.dflt.tensor ≈ values.first.tensor rtol=1e-12
        # On these meshes the first incident cell of every face on the plane is the one below.
        # The mixed side takes the cell below on odd faces and the cell above on even ones.
        n = length(mid); nb = count(isodd,1:n); na = n-nb
        @test 0 < na && 0 < nb
        @test values.mixed.normal ≈ (nb+4na)/n rtol=1e-12
        @test values.mixed.mass ≈ (nb*(dot(c,c)+1)+na*(dot(c,c)+4))/n rtol=1e-12
        @test values.mixed.load_n'*x ≈ (nb+2na)/n rtol=1e-12
        @test values.mixed.load_t ≈ values.below.load_t rtol=1e-12
        @test values.mixed.tensor ≈ im*(nb*dot(ub,T*ub)+na*dot(ua,T*ua))/n rtol=1e-12
    end
end

# On the skeleton, the tangential jump of a field of the space vanishes and the normal
# jump does not.
let
    mesh = GT.cartesian_mesh((0,1,0,1,0,1),(2,2,2);simplexify=true)
    Ω = GT.interior(mesh)
    dΛ = GT.measure(GT.skeleton(mesh),4)
    n = GT.unit_normal(mesh,2)
    V = GT.nedelec_space(Ω,1)
    uh = GT.interpolate(GT.analytical_field(x->SVector(sin(x[1]+x[3]),cos(2x[2]-x[3]),x[1]*x[3]+exp(x[3])),Ω),V)
    tangential = GT.assemble_scalar(GT.∫(q->norm(cross(n[1](q),uh[2](q)-uh[1](q)))^2,dΛ))
    normal = GT.assemble_scalar(GT.∫(q->(n[1](q)⋅(uh[2](q)-uh[1](q)))^2,dΛ))
    @test tangential < 1e-24
    @test normal > 1e-4
    A = GT.assemble_matrix((u,v)->GT.∫(q->(u[2](q)-u[1](q))⋅(v[2](q)-v[1](q)),dΛ),Float64,V,V)
    x = GT.free_values(uh)
    @test x'*A*x ≈ normal rtol=1e-10
end

# Unsupported uses raise errors instead of silently using the Lagrange map.
let
    mesh = GT.cartesian_mesh((0,1,0,1),(2,2);simplexify=true)
    Ω = GT.interior(mesh)
    V = GT.nedelec_space(Ω,0)
    dΩ = GT.measure(Ω,2)
    g(u,v) = GT.∫(x->ForwardDiff.gradient(u,x)⋅ForwardDiff.gradient(v,x),dΩ)
    @test_throws ArgumentError GT.assemble_matrix(g,Float64,V,V)
    # Only cells and faces of dimension D-1 are supported (here, edges of tetrahedra); the message
    # identifies the check.
    mesh3 = GT.cartesian_mesh((0,1,0,1,0,1),(1,1,1);simplexify=true)
    V3 = GT.nedelec_space(GT.interior(mesh3),0)
    dE = GT.measure(GT.domain(mesh3,1),2)
    @test_throws r"faces of dimension D-1 only" GT.assemble_matrix((u,v)->GT.∫(x->u(x)⋅v(x),dE),Float64,V3,V3)
    uh3 = GT.zero_field(Float64,V3)
    @test_throws r"faces of dimension D-1 only" GT.assemble_scalar(GT.∫(x->uh3(x)⋅uh3(x),dE))
    # The same check guards the accessors used outside assembly.
    @test_throws r"faces of dimension D-1 only" GT.space_accessor(V3,dE;tabulate=(GT.value,))
    @test_throws r"faces of dimension D-1 only" GT.field_accessor(uh3,dE)
    @test_throws r"faces of dimension D-1 only" GT.shape_function_accessor_reference(GT.value,V3,dE)
    @test_throws r"faces of dimension D-1 only" GT.shape_function_accessor(GT.value,V3,dE)
    # Curved (quadratic) triangles are rejected.
    curved = GT.mesh_from_msh(joinpath(@__DIR__,"..","assets","higher_order_2D.msh"))
    @test_throws ArgumentError GT.nedelec_space(GT.interior(curved),0)
end

end # module

module PeriodicTests

import GalerkinToolkit as GT
import ForwardDiff
using Test
using LinearAlgebra
using StaticArrays

curl(J::SMatrix{3,3}) = SVector(J[3,2]-J[2,3],J[1,3]-J[3,1],J[2,1]-J[1,2])
curl(J::SMatrix{2,2}) = J[2,1]-J[1,2]

# Unit square (D=2) or cube (D=3) meshed by gmsh with element size h, periodic in the directions
# listed in `periodic`. The sides are in the physical groups "x0", "x1", "y0", "y1" (, "z0", "z1").
function gmsh_periodic_box(D,h;periodic=1:D)
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
        function side(dir,t)
            e = 1e-6
            lo = [-e,-e,-e]
            hi = [1+e,1+e,D == 3 ? 1+e : e]
            lo[dir] = t-e
            hi[dir] = t+e
            [tag for (_,tag) in gmsh.model.getEntitiesInBoundingBox(lo...,hi...,D-1)]
        end
        for dir in periodic
            T = [1.0,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1]
            T[4*dir] = 1
            gmsh.model.mesh.setPeriodic(D-1,side(dir,1),side(dir,0),T)
        end
        for dir in 1:D, t in (0,1)
            group = gmsh.model.addPhysicalGroup(D-1,side(dir,t))
            gmsh.model.setPhysicalName(D-1,group,string("xyz"[dir],t))
        end
        gmsh.model.mesh.generate(D)
        GT.mesh_from_gmsh(gmsh)
    end
end

# A point reduced to the periodic cell: coordinates in the periodic directions are taken modulo 1, and
# rounded so that points that are periodic images of each other compare equal.
function reduced(x,periodic)
    SVector(ntuple(length(x)) do i
        r = round(x[i];digits=8)
        (i in periodic ? mod(r,1.0) : r) + 0.0
    end)
end

centroid(xs) = sum(xs)/length(xs)

# The number of distinct reduced points among `xs`.
num_reduced(xs,periodic) = length(unique(map(x->reduced(x,periodic),xs)))

# Periodic nodes of the mesh: every node and its owner are periodic images of each other.
let
    mesh = gmsh_periodic_box(2,0.25)
    owner = GT.periodic_node_owners(mesh)
    x = GT.node_coordinates(mesh)
    @test owner[owner] == owner
    @test count(i->owner[i] != i,eachindex(owner)) > 0
    @test all(i->reduced(x[i],1:2) == reduced(x[owner[i]],1:2),eachindex(owner))
    # The corner nodes are one node.
    corners = findall(xi->all(c->abs(c)<1e-9 || abs(c-1)<1e-9,xi),x)
    @test length(corners) == 4
    @test length(unique(owner[corners])) == 1

    # A cycle in the identification is an error.
    cmesh = GT.mesh(;
        node_coordinates = GT.node_coordinates(mesh),
        face_nodes = GT.face_nodes(mesh),
        face_reference_id = GT.face_reference_id(mesh),
        reference_spaces = GT.reference_spaces(mesh),
        periodic_nodes = Int32[1,2] => Int32[2,1])
    @test_throws ArgumentError GT.periodic_node_owners(cmesh)

    # A node listed twice must have one master.
    pmesh(pairs) = GT.mesh(;
        node_coordinates = GT.node_coordinates(mesh),
        face_nodes = GT.face_nodes(mesh),
        face_reference_id = GT.face_reference_id(mesh),
        reference_spaces = GT.reference_spaces(mesh),
        periodic_nodes = pairs)
    @test_throws ArgumentError GT.periodic_node_owners(pmesh(Int32[3,3] => Int32[1,2]))
    @test GT.periodic_node_owners(pmesh(Int32[3,3] => Int32[1,1]))[3] == 1
    # A node listed as its own master and with another master, in either order.
    @test_throws ArgumentError GT.periodic_node_owners(pmesh(Int32[3,3] => Int32[3,1]))
    @test_throws ArgumentError GT.periodic_node_owners(pmesh(Int32[3,3] => Int32[1,3]))
    @test GT.periodic_node_owners(pmesh(Int32[3,3] => Int32[3,3]))[3] == 3
end

# DOF identification of Lagrange spaces: two local DOFs share a periodic DOF exactly when their nodes are
# periodic images of each other.
function check_lagrange(mesh,D,periodic,order)
    Ω = GT.interior(mesh)
    V = GT.lagrange_space(Ω,order;periodic=Val(true))
    V0 = GT.lagrange_space(Ω,order)
    x0 = GT.node_coordinates(V0)
    dofs = GT.face_dofs(V)
    dofs0 = GT.face_dofs(V0)
    dof_to_point = Dict{Int,SVector{D,Float64}}()
    consistent = true
    for cell in 1:length(dofs)
        for (dof,dof0) in zip(dofs[cell],dofs0[cell])
            point = reduced(x0[dof0],periodic)
            consistent &= get!(dof_to_point,dof,point) == point
        end
    end
    @test consistent
    @test length(dof_to_point) == GT.num_free_dofs(V)
    @test length(unique(values(dof_to_point))) == GT.num_free_dofs(V)
    @test GT.num_free_dofs(V) == num_reduced(x0,periodic)
    @test GT.num_free_dofs(V) < GT.num_free_dofs(V0)
    @test GT.num_dirichlet_dofs(V) == 0
    V,V0
end

# DOF identification of Nédélec spaces: the DOFs owned by a d-face are shared with the d-faces that are its
# periodic images, and with no other face.
function check_nedelec(mesh,D,periodic,order)
    Ω = GT.interior(mesh)
    V = GT.nedelec_space(Ω,order;periodic=Val(true))
    V0 = GT.nedelec_space(Ω,order)
    topo = GT.topology(mesh)
    node_x = GT.node_coordinates(mesh)
    dofs = GT.face_dofs(V)
    reffes = GT.reference_spaces(V)
    cell_to_rid = GT.face_reference_id(V)
    dof_to_entity = Dict{Int,Tuple{Int,SVector{D,Float64}}}()
    entity_to_dofs = Dict{Tuple{Int,SVector{D,Float64}},Set{Int}}()
    entity_to_num_own = Dict{Tuple{Int,SVector{D,Float64}},Int}()
    consistent = true
    for d in 1:D
        cell_to_dfaces = GT.face_incidence(topo,D,d)
        dface_to_nodes = GT.face_nodes(mesh,d)
        for cell in 1:length(dofs)
            reffe = reffes[cell_to_rid[cell]]
            for (ldface,dface) in enumerate(cell_to_dfaces[cell])
                own = GT.face_own_dofs(reffe,d)[ldface]
                isempty(own) && continue
                entity = (d,reduced(centroid(node_x[dface_to_nodes[dface]]),periodic))
                entity_to_num_own[entity] = length(own)
                for ldof in own
                    dof = dofs[cell][ldof]
                    consistent &= get!(dof_to_entity,dof,entity) == entity
                    push!(get!(()->Set{Int}(),entity_to_dofs,entity),dof)
                end
            end
        end
    end
    @test consistent
    @test length(dof_to_entity) == GT.num_free_dofs(V)
    @test all(entity->length(entity_to_dofs[entity]) == entity_to_num_own[entity],keys(entity_to_dofs))
    # The count of reduced edges, faces and cells times their DOFs.
    expected = 0
    for d in 1:D
        nown = length(GT.face_own_dofs(reffes[1],d)[1])
        dface_to_nodes = GT.face_nodes(mesh,d)
        expected += nown*num_reduced([centroid(node_x[nodes]) for nodes in dface_to_nodes],periodic)
    end
    @test GT.num_free_dofs(V) == expected
    @test GT.num_free_dofs(V) < GT.num_free_dofs(V0)
    V,V0
end

# A periodic function is interpolated to the same field by the periodic and the non-periodic spaces. A wrong
# permutation of the DOFs of a merged face would move values to other nodes or moments.
function interpolation_gap(V,V0,ufun,degree)
    Ω = GT.domain(V)
    u = GT.analytical_field(ufun,Ω)
    uh = GT.zero_field(Float64,V)
    uh0 = GT.zero_field(Float64,V0)
    GT.interpolate!(u,uh)
    GT.interpolate!(u,uh0)
    dΩ = GT.measure(Ω,degree)
    GT.assemble_scalar(GT.∫(x->norm(uh(x)-uh0(x))^2,dΩ)), GT.assemble_scalar(GT.∫(x->norm(uh0(x))^2,dΩ))
end

function periodic_scalar(D)
    if D == 3
        x -> sinpi(2x[1])*cospi(2x[2])*cospi(2x[3]) + 0.5*sinpi(2(x[1]+x[3])) + 0.3
    else
        x -> sinpi(2x[1])*cospi(2x[2]) + 0.5*sinpi(2(x[1]+x[2])) + 0.3
    end
end

function periodic_vector(D)
    if D == 3
        x -> SVector(sinpi(2x[2])+0.5,cospi(2(x[1]+x[3])),sinpi(2x[1])*cospi(2x[2])*sinpi(2x[3])-0.2)
    else
        x -> SVector(sinpi(2x[1])*cospi(2x[2])+0.5,cospi(2(x[1]+x[2])))
    end
end

# Euler characteristic of the identified mesh, from the DOF counts of P1, N0 and N1: vertices, edges and (in 3D)
# faces, with χ = V - E + F (- C). A box that is periodic in any direction has χ = 0, a box with none has χ = 1.
function euler_characteristic(mesh,D)
    Ω = GT.interior(mesh)
    nv = GT.num_free_dofs(GT.lagrange_space(Ω,1;periodic=Val(true)))
    ne = GT.num_free_dofs(GT.nedelec_space(Ω,0;periodic=Val(true)))
    ncells = GT.num_faces(GT.topology(mesh),D)
    if D == 3
        nf = (GT.num_free_dofs(GT.nedelec_space(Ω,1;periodic=Val(true))) - 2*ne) ÷ 2
        nv - ne + nf - ncells
    else
        nv - ne + ncells
    end
end

configs = [(2,0.25,(1,)),(2,0.2,(1,2)),(3,0.35,(1,2)),(3,0.35,(1,2,3))]
for (D,h,periodic) in configs
    mesh = gmsh_periodic_box(D,h;periodic)
    for order in (D == 3 ? (1,2) : (1,2,3))
        V,V0 = check_lagrange(mesh,D,periodic,order)
        gap,norm0 = interpolation_gap(V,V0,periodic_scalar(D),2order+2)
        @test norm0 > 0.1
        @test gap < 1e-20
    end
    for order in (D == 3 ? (0,1) : (0,1,2))
        V,V0 = check_nedelec(mesh,D,periodic,order)
        gap,norm0 = interpolation_gap(V,V0,periodic_vector(D),2order+2)
        @test norm0 > 0.1
        # The face and cell moments place their quadrature points by the vertex order of the cell, which follows
        # the owners in the periodic space and the nodes in the other one; the interpolants then differ at the
        # size of the quadrature error. A wrong DOF permutation gives a gap of the size of norm0.
        @test gap < 1e-10*norm0
    end
    @test euler_characteristic(mesh,D) == 0
end
let
    mesh = gmsh_periodic_box(3,0.35;periodic=())
    @test euler_characteristic(mesh,3) == 1
    # Without periodic nodes, a periodic space equals the standard one.
    Ω = GT.interior(mesh)
    @test GT.num_free_dofs(GT.lagrange_space(Ω,2;periodic=Val(true))) == GT.num_free_dofs(GT.lagrange_space(Ω,2))
    @test GT.num_free_dofs(GT.nedelec_space(Ω,1;periodic=Val(true))) == GT.num_free_dofs(GT.nedelec_space(Ω,1))
end

# Tangential continuity of a Nédélec field with arbitrary DOF values across the periodic faces: the field on a face
# at x = 1 (resp. y = 1, z = 1) equals in tangential direction the field on its image face at x = 0.
function periodic_tangential_jumps(V,mesh,D,dir)
    node_x = GT.node_coordinates(mesh)
    topo = GT.topology(mesh)
    face_cells = GT.face_incidence(topo,D-1,D)
    face_nodes = GT.face_nodes(mesh,D-1)
    cell_dofs = GT.face_dofs(V)
    vals = [sin(1.7i)+0.3 for i in 1:GT.num_free_dofs(V)]
    s = GT.shape_function_accessor(GT.value,V)
    uh(cell,x) = sum(vals[dof]*s(cell)(i)(x) for (i,dof) in enumerate(cell_dofs[cell]))
    on_side(face,t) = all(node->abs(node_x[node][dir]-t) < 1e-9,face_nodes[face])
    shift = SVector(ntuple(i->Float64(i==dir),D))
    side0 = Dict(round.(centroid(node_x[face_nodes[face]]);digits=8) => face
                 for face in 1:length(face_nodes) if on_side(face,0))
    weights = D == 3 ? ((1/3,1/3,1/3),(0.6,0.3,0.1),(0.1,0.2,0.7)) : ((0.5,0.5),(0.2,0.8),(0.9,0.1))
    tangential,total,nchecked = 0.0,0.0,0
    for face in 1:length(face_nodes)
        on_side(face,1) || continue
        xs = node_x[face_nodes[face]]
        face0 = side0[round.(centroid(xs)-shift;digits=8)]
        cell1 = only(face_cells[face])
        cell0 = only(face_cells[face0])
        tangents = [xs[i+1]-xs[1] for i in 1:D-1]
        for w in weights
            x = sum(w[i]*xs[i] for i in 1:D)
            jump = uh(cell1,x)-uh(cell0,x-shift)
            tangential = max(tangential,maximum(t->abs(dot(t,jump)),tangents))
            total = max(total,norm(jump))
        end
        nchecked += 1
    end
    tangential,total,nchecked
end

for (D,h,periodic) in [(2,0.2,(1,2)),(3,0.35,(1,2))]
    mesh = gmsh_periodic_box(D,h;periodic)
    for order in (0,1,2)
        V = GT.nedelec_space(GT.interior(mesh),order;periodic=Val(true))
        for dir in periodic
            tangential,total,nchecked = periodic_tangential_jumps(V,mesh,D,dir)
            @test nchecked > 0
            @test tangential < 1e-10*max(1,total)
            @test total > 1e-3 # normal components do jump: the check is not vacuous
        end
    end
end

# Convergence of curl curl u + u = f with a periodic manufactured solution. In 2D the box is periodic in x and y;
# in 3D it is periodic in x and y with u × n = 0 on z = 0 and z = 1.
function solve_periodic_maxwell(D,h,k)
    mesh = gmsh_periodic_box(D,h;periodic=(1,2))
    Ω = GT.interior(mesh)
    if D == 3
        φ = x -> sinpi(2x[1])*sinpi(2x[2])*sinpi(x[3])
        ufun = x -> SVector(cospi(2x[2])*sinpi(x[3]),cospi(2x[1])*sinpi(x[3]),sinpi(2x[1])*sinpi(2x[2])) +
            ForwardDiff.gradient(φ,x)
        curlfun = x -> curl(ForwardDiff.jacobian(ufun,x))
        ffun = x -> curl(ForwardDiff.jacobian(curlfun,x)) + ufun(x)
        V = GT.nedelec_space(Ω,k;periodic=Val(true),dirichlet_boundary=GT.boundary(mesh;group_names=["z0","z1"]))
    else
        ufun = x -> SVector(sinpi(2x[1])*cospi(2x[2])+0.5,cospi(2(x[1]+x[2])))
        curlfun = x -> curl(ForwardDiff.jacobian(ufun,x))
        ffun = x -> begin
            g = ForwardDiff.gradient(curlfun,x)
            SVector(g[2],-g[1]) + ufun(x)
        end
        V = GT.nedelec_space(Ω,k;periodic=Val(true))
    end
    dΩ = GT.measure(Ω,2k+4)
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
    heff = GT.num_faces(GT.topology(mesh),D)^(-1/D)
    heff,el2,ecurl,GT.num_dirichlet_dofs(V)
end

# Convergence of -Δu + u = f with a periodic manufactured solution in a box periodic in every direction.
function solve_periodic_poisson(D,h,order)
    mesh = gmsh_periodic_box(D,h)
    Ω = GT.interior(mesh)
    V = GT.lagrange_space(Ω,order;periodic=Val(true))
    if D == 3
        ufun = x -> sinpi(2x[1])*cospi(2x[2])*cospi(2x[3]) + 0.3
        ffun = x -> (12π^2+1)*(ufun(x)-0.3) + 0.3
    else
        ufun = x -> sinpi(2x[1])*cospi(2x[2]) + 0.3
        ffun = x -> (8π^2+1)*(ufun(x)-0.3) + 0.3
    end
    dΩ = GT.measure(Ω,2order+4)
    u = GT.analytical_field(ufun,Ω)
    f = GT.analytical_field(ffun,Ω)
    a(u,v) = GT.∫(x->ForwardDiff.gradient(u,x)⋅ForwardDiff.gradient(v,x) + u(x)*v(x),dΩ)
    l(v) = GT.∫(x->f(x)*v(x),dΩ)
    A = GT.assemble_matrix(a,Float64,V,V)
    b = GT.assemble_vector(l,Float64,V)
    uh = GT.zero_field(Float64,V)
    GT.free_values(uh) .= A\b
    el2 = sqrt(GT.assemble_scalar(GT.∫(x->norm(u(x)-uh(x))^2,dΩ)))
    eh1 = sqrt(GT.assemble_scalar(GT.∫(x->norm(ForwardDiff.gradient(u,x)-ForwardDiff.gradient(uh,x))^2,dΩ)))
    heff = GT.num_faces(GT.topology(mesh),D)^(-1/D)
    heff,el2,eh1
end

rate(r1,r2) = log(r1[2]/r2[2])/log(r1[1]/r2[1]), log(r1[3]/r2[3])/log(r1[1]/r2[1])

for k in (0,1,2)
    results = [solve_periodic_maxwell(2,h,k) for h in (0.2,0.1,0.05)]
    rates = [rate(results[i],results[i+1]) for i in 1:2]
    @test results[end][2] < results[1][2]
    # The smallest rate is in the test expression, so a failure prints it.
    @test minimum(r->r[1],rates) > k+1-0.4
    @test minimum(r->r[2],rates) > k+1-0.4
end

for k in (0,1)
    results = [solve_periodic_maxwell(3,h,k) for h in (0.35,0.25,0.18)]
    rates = [rate(results[i],results[i+1]) for i in 1:2]
    @test all(r->r[4] > 0,results)
    @test results[end][2] < results[1][2]
    @test rates[end][1] > k+1-0.4
    @test rates[end][2] > k+1-0.4
end

for (D,order,hs) in [(2,1,(0.1,0.05,0.025)),(2,2,(0.2,0.1,0.05)),(3,1,(0.3,0.2,0.14))]
    results = [solve_periodic_poisson(D,h,order) for h in hs]
    # Rates over the whole range: the unstructured 3D meshes are coarse for this solution, and the rate between two
    # of them varies (for P1 in L2, 1.75 and 1.58 between the meshes here, and 1.83 from h = 0.14 to 0.1).
    r = rate(results[1],results[end])
    @test r[1] > order+1-0.5
    @test r[2] > order-0.4
end

# Face-interior DOFs with nontrivial permutations (3D P4: three DOFs per face).
let
    mesh = gmsh_periodic_box(3,0.35;periodic=(1,2,3))
    V,V0 = check_lagrange(mesh,3,(1,2,3),4)
    gap,norm0 = interpolation_gap(V,V0,periodic_scalar(3),10)
    @test norm0 > 0.1
    @test gap < 1e-20
end

# Vector-valued Lagrange spaces.
let
    mesh = gmsh_periodic_box(2,0.2;periodic=(1,2))
    Ω = GT.interior(mesh)
    V = GT.lagrange_space(Ω,2;tensor_size=Val((2,)),periodic=Val(true))
    V0 = GT.lagrange_space(Ω,2;tensor_size=Val((2,)))
    @test GT.num_free_dofs(V) == 2*GT.num_free_dofs(GT.lagrange_space(Ω,2;periodic=Val(true)))
    gap,norm0 = interpolation_gap(V,V0,periodic_vector(2),6)
    @test norm0 > 0.1
    @test gap < 1e-20
end

# The mesh with periodic nodes added: in each direction of `periodic`, a node on the side x_dir = 1 is identified
# with the node at the same point of the side x_dir = 0 (one master per node, from the first direction that
# applies; periodic_node_owners follows the chains at edges and corners).
function with_periodic_nodes(mesh,periodic)
    x = GT.node_coordinates(mesh)
    key(p) = Tuple(round.(p;digits=10) .+ 0.0)
    point_to_node = Dict(key(p)=>Int32(i) for (i,p) in enumerate(x))
    slaves,masters = Int32[],Int32[]
    for (i,p) in enumerate(x)
        for dir in periodic
            if abs(p[dir]-1) < 1e-12
                push!(slaves,i)
                push!(masters,point_to_node[key(Base.setindex(p,0.0,dir))])
                break
            end
        end
    end
    GT.mesh(;
        node_coordinates = x,
        face_nodes = GT.face_nodes(mesh),
        face_reference_id = GT.face_reference_id(mesh),
        reference_spaces = GT.reference_spaces(mesh),
        periodic_nodes = slaves => masters)
end

# Quadrilateral and hexahedral meshes.
for (D,cells,orders) in ((2,(3,3),(1,2,3)),(3,(3,3,3),(1,2)))
    mesh = with_periodic_nodes(GT.cartesian_mesh(ntuple(i->isodd(i) ? 0 : 1,2D),cells),1:D)
    for order in orders
        V,V0 = check_lagrange(mesh,D,1:D,order)
        gap,norm0 = interpolation_gap(V,V0,periodic_scalar(D),2order+2)
        @test norm0 > 0.1
        @test gap < 1e-20
    end
end

# Convex quadrilaterals that are not parallelograms: the interior nodes of a 3 x 3 mesh are moved, so no cell is
# affine and the bound on the periods comes from the convex quadrilaterals.
let
    mesh = GT.cartesian_mesh((0,1,0,1),(3,3))
    moves = Dict((1,1)=>SVector(0.04,0.03),(2,1)=>SVector(-0.02,0.05),(1,2)=>SVector(0.01,-0.04),
                 (2,2)=>SVector(-0.05,-0.01))
    x = map(GT.node_coordinates(mesh)) do p
        k = (round(Int,3p[1]),round(Int,3p[2]))
        p + get(moves,k,zero(p))
    end
    mesh = with_periodic_nodes(GT.replace_node_coordinates(mesh,x),1:2)
    cell_x = [x[nodes] for nodes in GT.face_nodes(mesh,2)]
    @test all(c->norm(c[1]+c[4]-c[2]-c[3]) > 0.005,cell_x)
    for order in (1,2)
        check_lagrange(mesh,2,1:2,order)
    end
end

# Faces with the same vertex owners that are not periodic images of each other. In the unit square, periodic in x
# and y, the edges from (0,1/2) to (1/2,1) and from (1/2,0) to (1,1/2) have the same vertex owners, (0,1/2) and
# (1/2,0), but matching their ends takes two different periods. They are distinct edges. The mesh keeps both and
# fills the three polygons they cut out with fans of triangles.
let
    points = SVector{2,Float64}[]
    point_to_node = Dict{Tuple{Float64,Float64},Int32}()
    node(p) = get!(point_to_node,Float64.(p)) do
        push!(points,SVector{2,Float64}(p))
        Int32(length(points))
    end
    fans = [((1/6,5/6),[(0,1/2),(1/2,1),(1/4,1),(0,1),(0,3/4)]),
            ((5/6,1/6),[(1/2,0),(3/4,0),(1,0),(1,1/4),(1,1/2)]),
            ((1/2,1/2),[(0,0),(1/4,0),(1/2,0),(1,1/2),(1,3/4),(1,1),(3/4,1),(1/2,1),(0,1/2),(0,1/4)])]
    triangles = Vector{Int32}[]
    for (c,loop) in fans, i in eachindex(loop)
        push!(triangles,Int32[node(c),node(loop[i]),node(loop[mod1(i+1,length(loop))])])
    end
    slaves,masters = Int32[],Int32[]
    for (i,p) in enumerate(points)
        q = abs(p[1]-1) < 1e-12 ? (0.0,p[2]) : abs(p[2]-1) < 1e-12 ? (p[1],0.0) : nothing
        q === nothing && continue
        push!(slaves,i)
        push!(masters,point_to_node[q])
    end
    order = 1
    reference_spaces = map(d->(GT.lagrange_space(GT.unit_simplex(Val(d)),order),),(0,1,2))
    build(x,pairs=slaves=>masters) = GT.complexify(GT.mesh(;
        node_coordinates = x,
        face_nodes = [Vector{Int32}[],Vector{Int32}[],triangles],
        face_reference_id = [Int8[],Int8[],fill(Int8(1),length(triangles))],
        reference_spaces,
        periodic_nodes = pairs))
    mesh = build(points)
    @test GT.num_faces(mesh,2) == 20
    check_lagrange(mesh,2,(1,2),1)
    check_lagrange(mesh,2,(1,2),2)
    check_nedelec(mesh,2,(1,2),0)
    check_nedelec(mesh,2,(1,2),1)
    # P2, N0 and N1 DOFs: 10 vertices and 30 edges once identified (20 triangles, 2 DOFs per cell for N1).
    counts(m) = (GT.num_free_dofs(GT.lagrange_space(GT.interior(m),2;periodic=Val(true))),
                 GT.num_free_dofs(GT.nedelec_space(GT.interior(m),0;periodic=Val(true))),
                 GT.num_free_dofs(GT.nedelec_space(GT.interior(m),1;periodic=Val(true))))
    @test counts(mesh) == (40,30,100)
    # The translation test is relative to the size and position of the mesh: scaled, moved, and both.
    @test counts(build(1e-8*points)) == (40,30,100)
    @test counts(build([p+SVector(1e6,-1e6) for p in points])) == (40,30,100)
    @test counts(build([1e-3*p+SVector(1e6,-1e6) for p in points])) == (40,30,100)
    # Cells of size about 1e-8 at coordinates about 1e6 (spacing of the doubles about 1e-10) are not resolved
    # well enough to tell images from distinct faces: an error, not a guess.
    @test_throws r"do not resolve" counts(build([1e-8*p+SVector(1e6,-1e6) for p in points]))
    @test_throws r"do not resolve" counts(build([2.0^-26*p+SVector(1e6,-1e6) for p in points]))
    # Oblique periods: the mesh mapped by (u,v) -> (u-v,1e-4*v), with the left side slave to the right and the
    # top slave to the bottom. The periods are a = (1,0) and b = (-1,1e-4); the offsets of the nodes from their
    # owners are -a, b and b-a, and the two distinct edges above differ by a+b = (0,1e-4), a period much shorter
    # than any offset.
    oslaves,omasters = Int32[],Int32[]
    for (i,p) in enumerate(points)
        q = abs(p[1]) < 1e-12 ? (1.0,p[2]) : abs(p[2]-1) < 1e-12 ? (p[1],0.0) : nothing
        q === nothing && continue
        push!(oslaves,i)
        push!(omasters,point_to_node[q])
    end
    omesh = build([SVector(p[1]-p[2],1e-4*p[2]) for p in points],oslaves=>omasters)
    owner = GT.periodic_node_owners(omesh)
    ox = GT.node_coordinates(omesh)
    offsets = unique(round.(ox[i]-ox[owner[i]];digits=12) .+ 0.0 for i in eachindex(ox) if owner[i] != i)
    @test Set(offsets) == Set([SVector(-1.0,0.0),SVector(-1.0,1e-4),SVector(-2.0,1e-4)])
    @test counts(omesh) == (40,30,100)
    # The same oblique mesh as a surface z = g(u-v) in 3D: the bound on the periods holds for full-dimensional
    # domains only, so the domain is rejected.
    g(X) = (f = X-floor(X); 4f*(1-f))
    smesh = build([SVector(p[1]-p[2],1e-4*p[2],g(p[1]-p[2])) for p in points],oslaves=>omasters)
    @test_throws r"same dimension" counts(smesh)
    # The same oblique mesh with quadratic geometry, (u,v) -> (u-v, 1e-4*v + g(u-v)), exact on each triangle: no
    # cell is affine, so there is no bound on the periods and the mesh is rejected. With straight quadratic cells,
    # (u,v) -> (u-v, 1e-4*v), it is not.
    function build2(map_uv)
        ξs = GT.node_coordinates(GT.lagrange_space(GT.unit_simplex(Val(2)),2))
        uv_to_node = Dict{Tuple{Float64,Float64},Int32}()
        uvs = SVector{2,Float64}[]
        cells2 = map(triangles) do tri
            a,b,c = points[tri[1]],points[tri[2]],points[tri[3]]
            map(ξs) do ξ
                uv = a + ξ[1]*(b-a) + ξ[2]*(c-a)
                key = Tuple(round.(uv;digits=12) .+ 0.0)
                get!(uv_to_node,key) do
                    push!(uvs,uv)
                    Int32(length(uvs))
                end
            end
        end
        s2,m2 = Int32[],Int32[]
        for (i,uv) in enumerate(uvs)
            q = abs(uv[1]) < 1e-12 ? (1.0,uv[2]) : abs(uv[2]-1) < 1e-12 ? (uv[1],0.0) : nothing
            q === nothing && continue
            push!(s2,i)
            push!(m2,uv_to_node[Tuple(round.(q;digits=12) .+ 0.0)])
        end
        GT.complexify(GT.mesh(;
            node_coordinates = map(map_uv,uvs),
            face_nodes = [Vector{Int32}[],Vector{Int32}[],cells2],
            face_reference_id = [Int8[],Int8[],fill(Int8(1),length(cells2))],
            reference_spaces = map(d->(GT.lagrange_space(GT.unit_simplex(Val(d)),2),),(0,1,2)),
            periodic_nodes = s2 => m2))
    end
    P2count(m) = GT.num_free_dofs(GT.lagrange_space(GT.interior(m),2;periodic=Val(true)))
    @test_throws r"affine" P2count(build2(uv->SVector(uv[1]-uv[2],1e-4*uv[2]+g(uv[1]-uv[2]))))
    @test P2count(build2(uv->SVector(uv[1]-uv[2],1e-4*uv[2]))) == 40
    # A periodic node that is not the translate of its owner up to rounding is reported.
    moved = copy(points)
    moved[point_to_node[(1.0,0.25)]] += SVector(1e-6,0.0)
    @test_throws r"not translates" counts(build(moved))
end

# Thin affine cells: a 4 x 4 grid of the unit square, each square split along its lower-left to upper-right
# diagonal, mapped by (u,v) -> (u+v,1e-9*v) and periodic in u and v. The cells have inradius about 6e-11, far above
# the tolerance, and must not be taken as degenerate: 16 vertices, 48 edges and 32 triangles once identified.
let
    n = 4
    node(i,j) = Int32(1 + i + (n+1)*j)
    uvs = [SVector(i/n,j/n) for j in 0:n for i in 0:n]
    triangles = Vector{Int32}[]
    for j in 0:(n-1), i in 0:(n-1)
        push!(triangles,[node(i,j),node(i+1,j),node(i+1,j+1)])
        push!(triangles,[node(i,j),node(i+1,j+1),node(i,j+1)])
    end
    slaves = Int32[node(i,j) for j in 0:n for i in 0:n if i == n || j == n]
    masters = Int32[i == n ? node(0,j) : node(i,0) for j in 0:n for i in 0:n if i == n || j == n]
    mesh = GT.complexify(GT.mesh(;
        node_coordinates = [SVector(uv[1]+uv[2],1e-9*uv[2]) for uv in uvs],
        face_nodes = [Vector{Int32}[],Vector{Int32}[],triangles],
        face_reference_id = [Int8[],Int8[],fill(Int8(1),length(triangles))],
        reference_spaces = map(d->(GT.lagrange_space(GT.unit_simplex(Val(d)),1),),(0,1,2)),
        periodic_nodes = slaves => masters))
    Ω = GT.interior(mesh)
    @test GT.num_free_dofs(GT.lagrange_space(Ω,2;periodic=Val(true))) == 16+48
    @test GT.num_free_dofs(GT.nedelec_space(Ω,0;periodic=Val(true))) == 48
    @test GT.num_free_dofs(GT.nedelec_space(Ω,1;periodic=Val(true))) == 2*48+2*32
end

# A domain with two periodic images of a cell is rejected; a domain with one of them is not. The interval [0,2]
# with the nodes at x >= 1 identified with the nodes at x - 1: the cells in [1,2] are images of those in [0,1].
let
    mesh = GT.cartesian_mesh((0,2),(4,))
    x = GT.node_coordinates(mesh)
    point_to_node = Dict(round(p[1];digits=10)+0.0=>Int32(i) for (i,p) in enumerate(x))
    slaves = Int32[i for (i,p) in enumerate(x) if p[1] > 1-1e-12]
    masters = Int32[point_to_node[round(x[i][1]-1;digits=10)+0.0] for i in slaves]
    pmesh = GT.mesh(;
        node_coordinates = x,
        face_nodes = GT.face_nodes(mesh),
        face_reference_id = GT.face_reference_id(mesh),
        reference_spaces = GT.reference_spaces(mesh),
        periodic_nodes = slaves => masters)
    @test_throws ArgumentError GT.lagrange_space(GT.interior(pmesh),1;periodic=Val(true))
    cell_x = [sum(i->x[i][1],nodes)/length(nodes) for nodes in GT.face_nodes(pmesh,1)]
    GT.group_faces(pmesh,1)["left"] = Int32.(findall(<(1),cell_x))
    V = GT.lagrange_space(GT.interior(pmesh;group_names=["left"]),1;periodic=Val(true))
    @test GT.num_free_dofs(V) == 2
    # A domain of points (dimension below that of the coordinates) is rejected; a domain without cells has no DOFs.
    @test_throws r"same dimension" GT.lagrange_space(GT.boundary(pmesh),1;periodic=Val(true))
    GT.group_faces(pmesh,1)["none"] = Int32[]
    @test GT.num_free_dofs(GT.lagrange_space(GT.interior(pmesh;group_names=["none"]),1;periodic=Val(true))) == 0
end

# A DOF shared by the two sides of a periodic boundary is a Dirichlet DOF if either side is on the Dirichlet
# boundary.
let
    mesh = gmsh_periodic_box(2,0.2;periodic=(1,))
    Ω = GT.interior(mesh)
    node_x = GT.node_coordinates(mesh)
    nvertices = count(x->abs(x[1]) < 1e-9,node_x)
    nedges = count(nodes->all(node->abs(node_x[node][1]) < 1e-9,nodes),GT.face_nodes(mesh,1))
    for name in ("x0","x1")
        Γ = GT.boundary(mesh;group_names=[name])
        V1 = GT.lagrange_space(Ω,1;periodic=Val(true),dirichlet_boundary=Γ)
        V2 = GT.lagrange_space(Ω,2;periodic=Val(true),dirichlet_boundary=Γ)
        N0 = GT.nedelec_space(Ω,0;periodic=Val(true),dirichlet_boundary=Γ)
        @test GT.num_dirichlet_dofs(V1) == nvertices
        @test GT.num_dirichlet_dofs(V2) == nvertices+nedges
        @test GT.num_dirichlet_dofs(N0) == nedges
        @test GT.num_dofs(V1) == GT.num_free_dofs(GT.lagrange_space(Ω,1;periodic=Val(true)))
        @test GT.num_dofs(N0) == GT.num_free_dofs(GT.nedelec_space(Ω,0;periodic=Val(true)))
    end
end

# Meshes that are too coarse in a periodic direction, and discontinuous spaces, are rejected.
let
    # One cell across a period: a cell would be glued to itself.
    mesh = gmsh_periodic_box(2,2.0)
    @test_throws ArgumentError GT.lagrange_space(GT.interior(mesh),1;periodic=Val(true))
    @test_throws ArgumentError GT.nedelec_space(GT.interior(mesh),0;periodic=Val(true))
    # Two cells across a period: distinct edges have the same vertex owners, but they are not translates of each
    # other, so they keep their own DOFs.
    mesh = gmsh_periodic_box(2,0.5)
    @test count(x->abs(x[2]) < 1e-9,GT.node_coordinates(mesh)) == 3
    check_lagrange(mesh,2,(1,2),1)
    check_lagrange(mesh,2,(1,2),2)
    check_nedelec(mesh,2,(1,2),0)
    mesh = gmsh_periodic_box(2,0.25)
    @test_throws ArgumentError GT.lagrange_space(GT.interior(mesh),1;continuous=Val(false),periodic=Val(true))
end

end # module

module DomainTests

import GalerkinToolkit as GT
import GalerkinToolkit
import PartitionedArrays as pa
using Test
using WriteVTK
using LinearAlgebra

outdir = mkpath(joinpath(@__DIR__,"..","output"))

domain = (0,1,0,1)
cells = (4,4)
mesh = GT.cartesian_mesh(domain,cells)
GT.label_interior_faces!(mesh;physical_name="interior_faces")
GT.label_boundary_faces!(mesh;physical_name="boundary_faces")

Ω = GT.interior(mesh)
Ωref = GT.interior(mesh;is_reference_domain=true)

D = GT.num_dims(mesh)
Γref = GT.boundary(mesh;
                 is_reference_domain=true,
                 physical_names=["boundary_faces"])

Γ = GT.physical_domain(Γref)

Λref = GT.skeleton(mesh;
                 is_reference_domain=true,
                 physical_names=["interior_faces"])

Λ = GT.physical_domain(Λref)

@test Ω == Ω
@test Ω != Ωref

q1 = GT.face_quantity([ f for f in 1:GT.num_faces(Ω)],Ω)
q2 = GT.face_quantity([ 10*f for f in 1:GT.num_faces(Ω)],Ω)
q3 = GT.face_quantity([ 100*f for f in 1:GT.num_faces(Γ)],Γ)
q4 = GT.face_quantity([ 1000*f for f in 1:GT.num_faces(Λ)],Λ)

@test GT.is_boundary(Ω) == false
@test GT.is_boundary(Γ) == true
@test GT.is_boundary(Λ) == false

q = q1 + q2
index = GT.generate_index(Ω)
t = GT.term(q,index)
@test t.dim == D
storage = GT.index_storage(index)
expr = quote
    $(GT.unpack_index_storage(index,:storage))
    $(GT.dummy_face_index(index,D)) = 3
    $(GT.topological_sort(t.expr,())[1])
end
display(expr)
r = eval(expr)
@test r == 33

q = q1 + q3
index = GT.generate_index(Γ)
faces = GT.get_symbol!(index,GT.faces(Γ),"faces")
t = GT.term(q,index)
@test t.dim == D-1
storage = GT.index_storage(index)
expr = quote
    $(GT.unpack_index_storage(index,:storage))
    $(GT.dummy_face_index(index,D-1)) = $faces[2]
    $(GT.topological_sort(t.expr,())[1])
end
display(expr)
r = eval(expr)
@test r == 201

q = (q1 + q4)[2]
index = GT.generate_index(Λ)
faces = GT.get_symbol!(index,GT.faces(Λ),"faces")
t = GT.term(q,index)
@test t.dim == D-1
storage = GT.index_storage(index)
expr = quote
    $(GT.unpack_index_storage(index,:storage))
    $(GT.dummy_face_index(index,D-1)) = $faces[2]
    $(GT.topological_sort(t.expr,())[1])
end
display(expr)
r = eval(expr)
@test r == 2002

q = (q1 + q4)[2]
form_arity = 1
index = GT.generate_index(Λ,form_arity)
faces = GT.get_symbol!(index,GT.faces(Λ),"faces")
t = GT.term(q,index)
@test t.dim == D-1
storage = GT.index_storage(index)
expr = quote
    $(GT.unpack_index_storage(index,:storage))
    $(GT.dummy_face_index(index,D-1)) = $faces[2]
    $(GT.face_around_index(index,1)) = 2
    $(GT.topological_sort(t.expr,())[1])
end
display(expr)
r = eval(expr)
@test r == 2002

q = (q1 + q4)[2]
form_arity = 1
index = GT.generate_index(Λ,form_arity)
faces = GT.get_symbol!(index,GT.faces(Λ),"faces")
t = GT.term(q,index)
@test t.dim == D-1
storage = GT.index_storage(index)
expr = quote
    $(GT.unpack_index_storage(index,:storage))
    $(GT.dummy_face_index(index,D-1)) = $faces[2]
    $(GT.face_around_index(index,1)) = 1
    $(GT.topological_sort(t.expr,())[1])
end
display(expr)
r = eval(expr)
@test r == 0

xx





u2 = GT.face_constant_field(facedata,Ω)
indices = Dict(Ω=>GT.index(face=:face))
dict = Dict{Any,Symbol}()
expr = GT.term(u2)(indices,dict)
storage = GT.storage(dict)
expr = quote
    $(GT.unpack_storage(dict,:storage))
    face = $nfaces
    $expr(0)
end
display(expr)
@test eval(expr) == facedata[end]
indices = Dict(Ω=>GT.index(face=GT.FaceList(:faces)))
dict = Dict{Any,Symbol}()
expr = GT.term(u2)(indices,dict)
storage = GT.storage(dict)
expr = quote
    $(GT.unpack_storage(dict,:storage))
    faces = (4,5,6)
    $expr[2](0)
end
display(expr)
@test eval(expr) == facedata[5]

xxx

ϕ = GT.domain_map(Ωref,Ω)

ϕinv = GT.inverse_map(ϕ)

uref = u∘ϕ
u2 = uref∘ϕinv

vtk_grid(joinpath(outdir,"omega"),Ω,;plot_params=(;refinement=4)) do plt
    GT.plot!(plt,u;label="u")
    GT.plot!(plt,u2;label="u2")
end

vtk_grid(joinpath(outdir,"omega_ref"),Ωref;plot_params=(;refinement=4)) do plt
    GT.plot!(plt,uref;label="u")
end

D = GT.num_dims(mesh)
Γref = GT.boundary(mesh;
                 is_reference_domain=true,
                 physical_names=["boundary_faces"])

ϕ = GT.domain_map(Γref,Ωref)
g = uref∘ϕ
@test GT.domain(g) === GT.domain(ϕ)
@test GT.domain(g) === Γref
Γ = GT.physical_domain(Γref)

n = GT.unit_normal(Γref,Ω)
n2 = GT.unit_normal(Γ,Ω)
#h = GT.face_diameter_field(Γ)

vtk_grid(joinpath(outdir,"gamma_ref"),Γref) do plt
    GT.plot!(plt,g;label="u")
    GT.plot!(plt,n;label="n")
    GT.plot!(plt,q->n2(ϕ(q));label="n2") # TODO
    # GT.plot!(plt,h;label="h") TODO
    GT.plot!(plt;label="u2") do q
        x = ϕ(q)
        uref(x)
    end
end

vtk_grid(joinpath(outdir,"gamma"),Γ) do plt
    GT.plot!(plt,u;label="u")
    GT.plot!(plt,n;label="n")
    GT.plot!(plt,n2;label="n2")
    # TODO
    #GT.plot!(plt,h;label="h")
end

Λref = GT.skeleton(mesh;
                 is_reference_domain=true,
                 physical_names=["interior_faces"])

ϕ = GT.domain_map(Λref,Ωref)

Λ = GT.physical_domain(Λref)

# TODO
n = GT.unit_normal(Λref,Ω)
n2 = GT.unit_normal(Λ,Ω)
#h = GT.face_diameter_field(Λ)

jump(u) = u[+] - u[-]

vtk_grid(joinpath(outdir,"lambda_ref"),Λref) do plt
    GT.plot!(plt,n[+];label="n1")
    GT.plot!(plt,n[-];label="n2")
    GT.plot!(plt;label="jump_u2") do q
        jump((uref∘ϕ)(q))
    end
    #GT.plot!(plt,h;label="h")
    GT.plot!(plt;label="jump_u") do q
        jump(uref(ϕ(q)))
    end
end

jump2(u) = q -> jump(u(q))

vtk_grid(joinpath(outdir,"lambda"),Λ) do plt
    GT.plot!(plt,n2[+];label="n1")
    GT.plot!(plt,n2[-];label="n2")
    GT.plot!(plt;label="jump_u") do q
        jump(u(q))
    end
    GT.plot!(plt,jump2(u);label="jump_u2")
    #GT.plot!(plt,h;label="h")
end

# Parallel

domain = (0,1,0,1)
cells_per_dir = (4,4)
parts_per_dir = (2,2)
np = prod(parts_per_dir)
parts = pa.DebugArray(LinearIndices((np,)))
# TODO make this strategy the default one
partition_strategy = GT.partition_strategy(graph_nodes=:cells,graph_edges=:nodes,ghost_layers=0)
mesh = GT.cartesian_mesh(domain,cells_per_dir;parts_per_dir,parts,partition_strategy)

GT.physical_faces(mesh,1)
GT.node_coordinates(mesh)
GT.face_nodes(mesh,2)
GT.face_nodes(mesh)
GT.periodic_nodes(mesh)
GT.outwards_normals(mesh)


# TODO
#GT.label_interior_faces!(mesh;physical_name="interior_faces")
# TODO There is a bug when using a ghost layer
GT.label_boundary_faces!(mesh;physical_name="boundary_faces")

Ω = GT.interior(mesh)
Ωref = GT.interior(mesh;is_reference_domain=true)

@test Ω == Ω
@test Ω != Ωref

u = GT.analytical_field(x->sum(x),Ω)
ϕ = GT.domain_map(Ωref,Ω)
uref = u∘ϕ

pa.partition(Ω)

GT.faces(Ω)

vtk_grid(joinpath(outdir,"p_omega"),Ω;plot_params=(;refinement=4)) do plt
    GT.plot!(plt,u;label="u")
    GT.plot!(plt,q->u(q);label="u2")
end

D = GT.num_dims(mesh)
Γref = GT.boundary(mesh;
                 is_reference_domain=true,
                 physical_names=["boundary_faces"])

ϕ = GT.domain_map(Γref,Ωref)
g = uref∘ϕ

vtk_grid(joinpath(outdir,"p_gamma_ref"),Γref) do plt
    GT.plot!(plt,g;label="u")
    GT.plot!(plt;label="u2") do q
        x = ϕ(q)
        uref(x)
    end
end

end # module

module DomainTests

import GalerkinToolkit as gk
import PartitionedArrays as pa
using Test

outdir = mkpath(joinpath(@__DIR__,"..","output"))

domain = (0,1,0,1)
cells = (4,4)
mesh = gk.cartesian_mesh(domain,cells)
gk.label_interior_faces!(mesh;physical_name="interior_faces")
gk.label_boundary_faces!(mesh;physical_name="boundary_faces")

Ω = gk.interior(mesh)
Ωref = gk.interior(mesh;is_reference_domain=true)

@test Ω == Ω
@test Ω != Ωref

u = gk.analytical_field(x->sum(x),Ω)
ϕ = gk.domain_map(Ωref,Ω)
ϕinv = gk.inverse_map(ϕ)

uref = u∘ϕ
u2 = uref∘ϕinv

gk.vtk_plot(joinpath(outdir,"omega"),Ω;refinement=4) do plt
    gk.plot!(plt,u;label="u")
    gk.plot!(plt,u2;label="u2")
end


gk.vtk_plot(joinpath(outdir,"omega_ref"),Ωref;refinement=4) do plt
    gk.plot!(plt,uref;label="u")
end

D = gk.num_dims(mesh)
Γref = gk.boundary(mesh;
                 is_reference_domain=true,
                 physical_names=["boundary_faces"])

ϕ = gk.domain_map(Γref,Ωref)
g = uref∘ϕ
@test gk.domain(g) === gk.domain(ϕ)
@test gk.domain(g) === Γref
Γ = gk.physical_domain(Γref)

n = gk.unit_normal(Γref,Ω)
n2 = gk.unit_normal(Γ,Ω)
h = gk.face_diameter_field(Γ)

gk.vtk_plot(joinpath(outdir,"gamma_ref"),Γref) do plt
    gk.plot!(plt,g;label="u")
    gk.plot!(plt,n;label="n")
    gk.plot!(plt,q->n2(ϕ(q));label="n2")
    gk.plot!(plt,h;label="h")
    gk.plot!(plt;label="u2") do q
        x = ϕ(q)
        uref(x)
    end
end

gk.vtk_plot(joinpath(outdir,"gamma"),Γ) do plt
    gk.plot!(plt,u;label="u")
    gk.plot!(plt,n;label="n")
    gk.plot!(plt,n2;label="n2")
    gk.plot!(plt,h;label="h")
end

Λref = gk.skeleton(mesh;
                 is_reference_domain=true,
                 physical_names=["interior_faces"])

ϕ = gk.domain_map(Λref,Ωref)

Λ = gk.physical_domain(Λref)

n = gk.unit_normal(Λref,Ω)

n2 = gk.unit_normal(Λ,Ω)
h = gk.face_diameter_field(Λ)

jump(u,q) = u[2](q)-u[1](q)
jump(u,ϕ,q) = u(ϕ[2](q))[2]-u(ϕ[1](q))[1]

gk.vtk_plot(joinpath(outdir,"lambda_ref"),Λref) do plt
    gk.plot!(plt,n[1];label="n1")
    gk.plot!(plt,n[2];label="n2")
    gk.plot!(plt;label="jump_u2") do q
        jump(uref∘ϕ,q)
    end
    gk.plot!(plt,h;label="h")
    gk.plot!(plt;label="jump_u") do q
        jump(uref,ϕ,q)
    end
end

jump(u) = u[2]-u[1]
jump(u,x) = u(x)[2]-u(x)[1]
gk.vtk_plot(joinpath(outdir,"lambda"),Λ) do plt
    gk.plot!(plt,n2[1];label="n1")
    gk.plot!(plt,n2[2];label="n2")
    gk.plot!(plt;label="jump_u") do q
        jump(u,q)
    end
    gk.plot!(plt;label="jump_u2") do q
        jump(u(q))
    end
    gk.plot!(plt,h;label="h")
end

# Parallel

domain = (0,1,0,1)
cells_per_dir = (4,4)
parts_per_dir = (2,2)
np = prod(parts_per_dir)
parts = pa.DebugArray(LinearIndices((np,)))
# TODO make this strategy the default one
partition_strategy = gk.partition_strategy(graph_nodes=:cells,graph_edges=:nodes,ghost_layers=0)
mesh = gk.cartesian_mesh(domain,cells_per_dir;parts_per_dir,parts,partition_strategy)

gk.physical_faces(mesh,1)
gk.node_coordinates(mesh)
gk.face_nodes(mesh,2)
gk.face_nodes(mesh)
gk.periodic_nodes(mesh)
gk.outwards_normals(mesh)


# TODO
#gk.label_interior_faces!(mesh;physical_name="interior_faces")
# TODO There is a bug when using a ghost layer
gk.label_boundary_faces!(mesh;physical_name="boundary_faces")

Ω = gk.interior(mesh)
Ωref = gk.interior(mesh;is_reference_domain=true)

@test Ω == Ω
@test Ω != Ωref

u = gk.analytical_field(x->sum(x),Ω)
ϕ = gk.domain_map(Ωref,Ω)
uref = u∘ϕ

pa.partition(Ω)

gk.faces(Ω)

gk.vtk_plot(joinpath(outdir,"p_omega"),Ω;refinement=4) do plt
    gk.plot!(plt,u;label="u")
    gk.plot!(plt,q->u(q);label="u2")
end

D = gk.num_dims(mesh)
Γref = gk.boundary(mesh;
                 is_reference_domain=true,
                 physical_names=["boundary_faces"])

ϕ = gk.domain_map(Γref,Ωref)
g = uref∘ϕ

gk.vtk_plot(joinpath(outdir,"p_gamma_ref"),Γref) do plt
    gk.plot!(plt,g;label="u")
    gk.plot!(plt;label="u2") do q
        x = ϕ(q)
        uref(x)
    end
end

end # module

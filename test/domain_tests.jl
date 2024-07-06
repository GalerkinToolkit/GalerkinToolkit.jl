module DomainTests

import GalerkinToolkit as GT
import PartitionedArrays as pa
using Test

outdir = mkpath(joinpath(@__DIR__,"..","output"))

domain = (0,1,0,1)
cells = (4,4)
mesh = GT.cartesian_mesh(domain,cells)
GT.label_interior_faces!(mesh;physical_name="interior_faces")
GT.label_boundary_faces!(mesh;physical_name="boundary_faces")

Ω = GT.interior(mesh)
Ωref = GT.interior(mesh;is_reference_domain=true)

@test Ω == Ω
@test Ω != Ωref

u = GT.analytical_field(sum,Ω)

ϕ = GT.domain_map(Ωref,Ω)

ϕinv = GT.inverse_map(ϕ)

uref = u∘ϕ
u2 = uref∘ϕinv

GT.vtk_plot(joinpath(outdir,"omega"),Ω;refinement=4) do plt
    GT.plot!(plt,u;label="u")
    GT.plot!(plt,u2;label="u2")
end

GT.vtk_plot(joinpath(outdir,"omega_ref"),Ωref;refinement=4) do plt
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

GT.vtk_plot(joinpath(outdir,"gamma_ref"),Γref) do plt
    GT.plot!(plt,g;label="u")
    GT.plot!(plt,n;label="n")
    GT.plot!(plt,q->n2(ϕ(q));label="n2") # TODO
    # GT.plot!(plt,h;label="h") TODO
    GT.plot!(plt;label="u2") do q
        x = ϕ(q)
        uref(x)
    end
end

GT.vtk_plot(joinpath(outdir,"gamma"),Γ) do plt
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

GT.vtk_plot(joinpath(outdir,"lambda_ref"),Λref) do plt
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

GT.vtk_plot(joinpath(outdir,"lambda"),Λ) do plt
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

GT.vtk_plot(joinpath(outdir,"p_omega"),Ω;refinement=4) do plt
    GT.plot!(plt,u;label="u")
    GT.plot!(plt,q->u(q);label="u2")
end

D = GT.num_dims(mesh)
Γref = GT.boundary(mesh;
                 is_reference_domain=true,
                 physical_names=["boundary_faces"])

ϕ = GT.domain_map(Γref,Ωref)
g = uref∘ϕ

GT.vtk_plot(joinpath(outdir,"p_gamma_ref"),Γref) do plt
    GT.plot!(plt,g;label="u")
    GT.plot!(plt;label="u2") do q
        x = ϕ(q)
        uref(x)
    end
end

end # module

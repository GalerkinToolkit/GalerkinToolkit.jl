module PoissonBDDC

import GalerkinToolkit as glk
import ForwardDiff
using StaticArrays
using LinearAlgebra
using SparseArrays
using WriteVTK
using PartitionedArrays
using Metis

const CORNER = 0
const EDGE = 1
const FACE = 2

function main(params_in)

    # Process params
    params_default = default_params()
    params = add_default_params(params_in,params_default)

    # Setup
    pmesh = params[:pmesh]
    bddc_type = params[:bddc_type]
    u = params[:u]
    ranks = linear_indices(pmesh)
    np = length(ranks)

    # Setup dirichlet_bcs
    dirichlet_tags = params[:dirichlet_tags]
    dirichlet_bcs = map(pmesh) do mesh
        node_to_tag = zeros(glk.num_nodes(mesh))
        tag_to_name = dirichlet_tags
        glk.classify_mesh_nodes!(node_to_tag,mesh,tag_to_name)
        free_and_dirichlet_nodes = glk.partition_from_mask(i->i==0,node_to_tag)
        node_to_x = glk.node_coordinates(mesh)
        dirichlet_nodes = last(free_and_dirichlet_nodes)
        x_dirichlet = view(node_to_x,dirichlet_nodes)
        u_dirichlet = u.(x_dirichlet)
        (;free_and_dirichlet_nodes,u_dirichlet,node_to_tag)
    end

    coarse_dofs = setup_coarse_dofs(pmesh,dirichlet_bcs,bddc_type)

    display(coarse_dofs)

    # Visualize mesh
    example_path = params[:example_path]
    map(pmesh,ranks) do mesh,rank
        pvtk_grid(example_path,glk.vtk_args(mesh)...;part=rank,nparts=np) do vtk
            glk.vtk_physical_groups!(vtk,mesh)
            vtk["piece",VTKCellData()] = fill(rank,sum(glk.num_faces(mesh)))
            vtk["owner",VTKPointData()] = local_to_owner(glk.local_nodes(mesh))
            vtk["interface",VTKPointData()] = map(colors->Int(length(colors)!=1),glk.local_node_colors(mesh))
        end
    end

end

function add_default_params(params_in,params_default)
    UmD = setdiff(keys(params_in),keys(params_default))
    for k in UmD
        @warn "Parameter with key :$k is unused"
    end
    UiD = intersect(keys(params_default),keys(params_in))
    DmU = setdiff(keys(params_default),keys(params_in))
    a = [ k=>params_in[k] for k in UiD]
    b = [ k=>params_default[k] for k in DmU]
    Dict(vcat(a,b))
end

function default_params()
    mesh = glk.cartesian_mesh((0,1,0,1),(10,10),complexify=false)
    np = 4
    ranks = DebugArray(LinearIndices((np,)))
    pmesh = glk.partition_mesh(Metis.partition,ranks,mesh,via=:cells)
    outdir = mkpath(joinpath(@__DIR__,"..","output"))
    params = Dict{Symbol,Any}()
    params[:pmesh] = glk.cartesian_mesh((0,1,0,1),(10,10))
    params[:u] = (x) -> sum(x)
    params[:f] = (x) -> 0.0
    params[:g] = (x) -> 0.0
    params[:dirichlet_tags] = ["boundary"]
    params[:neumann_tags] = String[]
    params[:integration_degree] = 2
    params[:solve] = \
    params[:example_path] = joinpath(outdir,"poisson_parallel")
    params[:bddc_type] = (CORNER,EDGE,FACE)
    params
end

function setup_coarse_dofs(pmesh,dirichlet_bcs,types)
    ranks = linear_indices(pmesh)
    lnode_to_colors = map(glk.local_node_colors,pmesh)
    node_partition = map(glk.local_nodes,pmesh)
    cldof_to_ldofs, cldof_to_owner, nown = map(ranks,lnode_to_colors,dirichlet_bcs) do rank,lnode_to_colors,dirichlet_bcs
        free_and_dirichlet_lnodes = dirichlet_bcs.free_and_dirichlet_nodes
        nldofs = length(first(free_and_dirichlet_lnodes))
        lnode_to_ldof = glk.permutation(free_and_dirichlet_lnodes)
        nlnodes = length(lnode_to_colors)
        interface_node_to_lnode = findall(1:nlnodes) do lnode
            colors = lnode_to_colors[lnode]
            ldof = lnode_to_ldof[lnode]
            length(colors)!=1 && ldof <= nldofs
        end
        interface_node_to_colors = view(lnode_to_colors,interface_node_to_lnode)
        # assumes sorted colors
        clnode_to_colors = unique(interface_node_to_colors)
        clnode_to_owner = map(maximum,clnode_to_colors)
        nclnodes = length(clnode_to_colors)
        interface_node_to_clnode = indexin(interface_node_to_colors,clnode_to_colors)
        clnode_to_interface_nodes = glk.inverse_index_map(interface_node_to_clnode,nclnodes)
        f = interface_node->lnode_to_ldof[interface_node_to_lnode[interface_node]]
        clnode_to_interface_nodes.data .= f.(clnode_to_interface_nodes.data)
        clnode_to_ldofs = clnode_to_interface_nodes
        clnode_to_type = fill(EDGE,nclnodes)
        clnode_to_type[ map(i->length(i) == 2,clnode_to_colors) ] .= FACE
        clnode_to_type[ map(i->length(i) == 1,clnode_to_ldofs) ] .= CORNER
        mask = map(i->i in types,clnode_to_type)
        my_cldof_to_ldofs = JaggedArray(clnode_to_ldofs[mask])
        my_cldof_to_owner = clnode_to_owner[mask]
        my_nown = count(owner->owner==rank,my_cldof_to_owner)
        my_cldof_to_ldofs, my_cldof_to_owner, my_nown
    end |> tuple_of_arrays
    ncdofs = sum(nown)
    cdof_partition = variable_partition(nown,ncdofs)
    v = PVector{Vector{Int}}(undef,node_partition)
    map(ranks,cdof_partition,cldof_to_owner,cldof_to_ldofs,local_values(v),dirichlet_bcs) do rank, codofs, cldof_to_owner, cldof_to_ldofs, lnode_to_v, dirichlet_bcs
        ldof_to_lnode = first(dirichlet_bcs.free_and_dirichlet_nodes)
        codof_to_cdof = own_to_global(codofs)
        codof = 0
        for (cldof,owner) in enumerate(cldof_to_owner)
            if owner != rank
                continue
            end
            codof += 1
            cdof = codof_to_cdof[codof]
            ldofs = cldof_to_ldofs[cldof]
            lnodes = view(ldof_to_lnode,ldofs)
            lnode_to_v[lnodes] .= cdof
        end
    end
    consistent!(v) |> wait
    cldof_to_cdof = map(ranks,cldof_to_ldofs,local_values(v),dirichlet_bcs) do rank, cldof_to_ldofs, lnode_to_v,dirichlet_bcs
        ldof_to_lnode = first(dirichlet_bcs.free_and_dirichlet_nodes)
        ncldofs = length(cldof_to_ldofs)
        my_cldof_to_cdof = zeros(Int,ncldofs)
        for cldof in 1:ncldofs
            ldofs = cldof_to_ldofs[cldof]
            cdof = lnode_to_v[ldof_to_lnode[first(ldofs)]]
            my_cldof_to_cdof[cldof] = cdof
        end
        my_cldof_to_cdof
    end
    rank_to_cdofs = gather(cldof_to_cdof,destination=MAIN)
    map(cldof_to_ldofs,rank_to_cdofs) do cldof_to_ldofs,rank_to_cdofs
        (;cldof_to_ldofs, rank_to_cdofs, num_cdofs = ncdofs)
    end
end


end # module

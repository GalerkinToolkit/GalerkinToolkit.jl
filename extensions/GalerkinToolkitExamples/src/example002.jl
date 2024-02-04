module Example002

import GalerkinToolkit as gk
import GalerkinToolkitExamples: Example001
import ForwardDiff
using StaticArrays
using LinearAlgebra
using SparseArrays
using WriteVTK
using PartitionedArrays

function main(params_in)

    # Dict to collect results
    results = Dict{Symbol,Any}()

    # Process params
    params_default = default_params()
    params = add_default_params(params_in,params_default)

    # Setup main data structures
    state = setup(params)
    add_basic_info(results,params,state)

    # Assemble system and solve it
    A,b = assemble_system(state)
    x = solve_system(A,b,params)

    # Post process
    uh = setup_uh(x,state)
    integrate_error_norms(results,uh,state)
    export_results(uh,params,state)

    results
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
    np = 4
    parts = DebugArray(LinearIndices((np,)))
    domain = (0,1,0,1)
    cells_per_dir = (10,10)
    parts_per_dir = (2,2)
    ghost_layers = 0
    mesh = gk.cartesian_mesh(domain,cells_per_dir,parts_per_dir;parts,ghost_layers)
    outdir = mkpath(joinpath(@__DIR__,"..","output"))
    params = Example001.default_params()
    params[:mesh] = mesh
    params[:example_path] = joinpath(outdir,"example002")
    params
end

function setup(params)
    mesh = params[:mesh]
    local_states = map(partition(mesh)) do mesh
        local_params = copy(params)
        local_params[:mesh] = mesh
        Example001.setup(local_params)
    end

    node_partition = gk.node_partition(mesh)
    global_node_to_mask = pfill(false,node_partition)
    function fillmask!(node_to_mask,state)
        free_nodes = first(state.dirichlet_bcs.free_and_dirichlet_nodes)
        node_to_mask[free_nodes] .= true
    end
    map(fillmask!,partition(global_node_to_mask),local_states)
    global_dof_to_node, global_node_to_dof = find_local_indices(global_node_to_mask)
    dof_partition = partition(axes(global_dof_to_node,1))
    state = (;local_states,global_dof_to_node, global_node_to_dof)
    state
end

function add_basic_info(results,params,state)
    mesh = params[:mesh]
    nfree = length(state.global_dof_to_node)
    nnodes = gk.num_nodes(mesh)
    ncells = gk.num_faces(mesh,gk.num_dims(mesh))
    results[:nfree] = nfree
    results[:nnodes] = nnodes
    results[:ncells] = ncells
    results
end

function assemble_system(state)
    dof_partition = partition(axes(state.global_dof_to_node,1))
    b = pzeros(Float64,dof_partition)
    I,J,V,free_node_to_b = map(Example001.assemble_system_loop,state.local_states) |> tuple_of_arrays
    function map_dofs(state,node_to_dof,dofs,I,J,free_node_to_b,dof_to_b)
        free_node_to_node = first(state.dirichlet_bcs.free_and_dirichlet_nodes)
        dof_to_global_dof = local_to_global(dofs)
        function free_node_to_global_dof(free_node)
            # TODO this free_node -> dof should not be needed of the order free dofs conveniently
            node = free_node_to_node[free_node]
            dof = node_to_dof[node]
            global_dof = dof_to_global_dof[dof]
            global_dof
        end
        I .= free_node_to_global_dof.(I)
        J .= free_node_to_global_dof.(J)
        free_node_to_dof = view(node_to_dof,free_node_to_node)
        # TODO this one should not be needed of the order free dofs conveniently
        dof_to_b[free_node_to_dof] = free_node_to_b
        nothing
    end
    map(map_dofs,state.local_states,partition(state.global_node_to_dof),dof_partition,I,J,free_node_to_b,partition(b))
    A = psparse(I,J,V,dof_partition,dof_partition;assemble=false) |> fetch
    # TODO
    A = assemble(A) |> fetch
    b = assemble(b,partition(axes(A,1))) |> fetch
    A,b
end

function solve_system(A,b,params)
    solver = params[:solver]
    x = similar(b,axes(A,2))
    setup = solver.setup(x,A,b)
    solver.solve!(x,setup,b)
    solver.finalize!(setup)
    x
end

function setup_uh(x,state)
    dofs = axes(state.global_dof_to_node,1)
    local_states = state.local_states

    # TODO this similar can be avoided when using
    # a sub-assembled system
    global_dof_to_x = similar(x,dofs)
    global_dof_to_x .= x
    consistent!(global_dof_to_x) |> wait

    # TODO this map can be avoided
    # when the local numeration of free nodes coincides
    # with local dofs.
    function map_dofs(dof_to_x,node_to_dof,state)
        free_node_to_node = first(state.dirichlet_bcs.free_and_dirichlet_nodes)
        free_node_to_dof = view(node_to_dof,free_node_to_node)
        view(dof_to_x,free_node_to_dof)
    end
    free_node_to_x = map(map_dofs,partition(global_dof_to_x),partition(state.global_node_to_dof),state.local_states)

    uh = map(Example001.setup_uh,free_node_to_x,local_states)
    uh
end

function integrate_error_norms(results,uh,state)
    local_states = state.local_states
    eh1², el2² = map(Example001.integrate_error_norms_loop,uh,local_states) |> tuple_of_arrays
    eh1 = sqrt(sum(eh1²))
    el2 = sqrt(sum(el2²))
    results[:eh1] = eh1
    results[:el2] = el2
    results
end

function export_results(uh,params,state)
    example_path = params[:example_path]
    mesh = params[:mesh]
    ranks = linear_indices(gk.node_partition(mesh))
    np = length(ranks)
    map(partition(mesh),ranks,uh,gk.index_partition(mesh)) do mesh,rank,uh,ids
        pvtk_grid(example_path,gk.vtk_args(mesh)...;part=rank,nparts=np) do vtk
            gk.vtk_physical_faces!(vtk,mesh)
            vtk["piece",VTKCellData()] = fill(rank,sum(gk.num_faces(mesh)))
            vtk["owner",VTKPointData()] = local_to_owner(gk.node_indices(ids))
            vtk["uh",VTKPointData()] = uh
        end
    end
end

end # module

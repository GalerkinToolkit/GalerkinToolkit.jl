module CUDACompilerDataTests

using Test
using Adapt
using CUDA
using StaticArrays
import ForwardDiff
import GalerkinToolkit as GT

# Read every prepared array from parallel device code, including packed static
# vectors and matrices. Outputs have exactly the original shapes and element types.
@inline function read_array!(output,input,i)
    if i <= length(output)
        @inbounds output[i] = input[i]
    end
    nothing
end

function read_prepared!(outputs,mesh,space,domain,dofs)
    i = (blockIdx().x-1)*blockDim().x+threadIdx().x
    read_array!(outputs[1],GT.node_coordinates(mesh),i)
    read_array!(outputs[2],GT.face_nodes(mesh),i)
    read_array!(outputs[3],GT.reference_gradients(mesh),i)
    read_array!(outputs[4],GT.reference_weights(mesh),i)
    read_array!(outputs[5],GT.reference_tabulator(space),i)
    read_array!(outputs[6],GT.domain_faces(domain),i)
    read_array!(outputs[7],GT.face_dofs(dofs),i)
    nothing
end

# Exercise domain-local -> global cell -> node/DOF addressing after adaptation.
function read_selected_cells!(coordinates,ids,mesh,domain,dofs)
    i = (blockIdx().x-1)*blockDim().x+threadIdx().x
    if i <= length(coordinates)
        @inbounds begin
            cell = GT.domain_faces(domain)[i]
            node = GT.face_nodes(mesh)[cell][1]
            coordinates[i] = GT.node_coordinates(mesh)[node]
            ids[i] = GT.face_dofs(dofs)[1,cell]
        end
    end
    nothing
end

function check_prepared(host)
    mesh,space,domain,dofs = host
    arrays = (mesh.node_coordinates,mesh.face_nodes,mesh.reference_gradients,
        mesh.reference_weights,space.reference_tabulator,domain.domain_faces,dofs.face_dofs)
    device = adapt(CuArray,host)
    device_mesh,device_space,device_domain,device_dofs = device
    device_arrays = (device_mesh.node_coordinates,device_mesh.face_nodes,
        device_mesh.reference_gradients,device_mesh.reference_weights,
        device_space.reference_tabulator,device_domain.domain_faces,device_dofs.face_dofs)
    @test device_mesh.layout isa GT.ElementMajorLayout
    @test device_space.tabulate_function === space.tabulate_function
    @test isbitstype(typeof(CUDA.cudaconvert(device)))
    @test device_mesh.node_coordinates isa CuVector{eltype(mesh.node_coordinates)}
    @test device_mesh.face_nodes isa CuVector{eltype(mesh.face_nodes)}
    for (original,uploaded) in zip(arrays,device_arrays)
        @test uploaded isa CuArray
        @test eltype(uploaded) === eltype(original)
        @test size(uploaded) == size(original)
        @test Array(uploaded) == original
    end

    restored_mesh,restored_space,restored_domain,restored_dofs = adapt(Array,device)
    restored = (restored_mesh.node_coordinates,restored_mesh.face_nodes,
        restored_mesh.reference_gradients,restored_mesh.reference_weights,
        restored_space.reference_tabulator,restored_domain.domain_faces,restored_dofs.face_dofs)
    @test restored_mesh.layout isa GT.ElementMajorLayout
    @test restored_space.tabulate_function === space.tabulate_function
    for (original,roundtrip) in zip(arrays,restored)
        @test typeof(roundtrip) === typeof(original)
        @test roundtrip == original
    end

    outputs = map(a->similar(a),device_arrays)
    count = maximum(length,arrays)
    for _ in 1:2
        if count != 0
            CUDA.@sync @cuda threads=128 blocks=cld(count,128) read_prepared!(outputs,device...)
        end
        for (original,result) in zip(arrays,outputs)
            @test Array(result) == original
        end
    end

    coordinates = similar(device_mesh.node_coordinates,length(domain.domain_faces))
    ids = similar(device_dofs.face_dofs,length(domain.domain_faces))
    if !isempty(domain.domain_faces)
        CUDA.@sync @cuda threads=128 blocks=cld(length(ids),128) read_selected_cells!(
            coordinates,ids,device_mesh,device_domain,device_dofs)
    end
    @test Array(coordinates) == [mesh.node_coordinates[mesh.face_nodes[cell][1]] for cell in domain.domain_faces]
    @test Array(ids) == [dofs.face_dofs[1,cell] for cell in domain.domain_faces]
    nothing
end

@testset "Element-major prepared data adaptation" begin
    for T in (Float32,Float64), Ti in (Int32,Int64), nelements in (0,3)
        nnodes = nelements == 0 ? 0 : 6
        coordinates = [SVector{3,T}(n,2n,3n) for n in 1:nnodes]
        nodes = [SVector{4,Ti}(1,2,3,4) for _ in 1:nelements]
        gradients = [SVector{3,T}(n,p,n+p) for n in 1:4,p in 1:3]
        weights = T[0.25,0.5,0.25]
        mesh = GT.MeshIntegrationData(GT.ElementMajorLayout(),coordinates,nodes,gradients,weights)
        domain = GT.DomainIntegrationData(nelements == 0 ? Ti[] : Ti[3,1])
        dofs = GT.DOFData(reshape(Ti[isodd(i) ? i : -i for i in 1:4nelements],4,nelements))
        for (f,table) in ((GT.value,T[n+p for n in 1:4,p in 1:3]),
                          (ForwardDiff.gradient,gradients))
            space = GT.SpaceIntegrationData(f,table)
            check_prepared((mesh,space,domain,dofs))
        end
    end
end

@testset "Prepared single-field CG data on CUDA" begin
    for dim in 1:3,k in 1:3
        box = ntuple(i->isodd(i) ? 0 : 1,2dim)
        mesh = GT.cartesian_mesh(box,ntuple(_->2,dim);simplexify=false)
        omega = GT.interior(mesh)
        # Keep both positive/free and negative/Dirichlet DOF IDs in the transfer.
        V = GT.lagrange_space(omega,k;dirichlet_boundary=GT.boundary(mesh))
        GT.group_faces(mesh,dim)["gpu prepared subset"] = Int32[GT.num_faces(omega),1]
        subset = GT.interior(mesh;group_names=["gpu prepared subset"])
        quadrature = GT.quadrature(subset,2k)
        mesh_data = GT.mesh_integration_data(mesh,quadrature)
        domain_data = GT.domain_integration_data(subset)
        dof_data = GT.DOFData(V)
        @test any(<(0),dof_data.face_dofs)
        @test any(>(0),dof_data.face_dofs)
        for f in (GT.value,ForwardDiff.gradient)
            space_data = GT.space_integration_data(f,V,quadrature)
            check_prepared((mesh_data,space_data,domain_data,dof_data))
        end
    end
end

end # module

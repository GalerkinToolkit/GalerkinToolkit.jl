module CUDAPoissonNumericTests

using Test
using Adapt
using CUDA
using StaticArrays
using LinearAlgebra
import ForwardDiff
import GalerkinToolkit as GT

function poisson_numeric_paths!(jacobians,inverses,weights,gradients,coordinates,
        mesh,space,geometry_values,domain,jacobian)
    index = (blockIdx().x-1)*blockDim().x+threadIdx().x
    npoints = size(jacobians,1)
    if index <= length(jacobians)
        point = mod1(index,npoints)
        element = cld(index,npoints)
        @inbounds begin
            cell = GT.domain_faces(domain)[element]
            nodes = GT.face_nodes(mesh)[cell]
            J = jacobian(nodes,GT.node_coordinates(mesh),GT.reference_gradients(mesh),point)
            determinant = det(J)
            inverse = GT.adjugate(J)*inv(determinant)
            jacobians[point,element] = J
            inverses[point,element] = inverse
            weights[point,element] = abs(determinant)*GT.reference_weights(mesh)[point]
            coordinates[point,element] = GT.physical_point_value(GT.node_coordinates(mesh),
                nodes,GT.reference_tabulator(geometry_values),point)
            for dof in axes(gradients,1)
                gradients[dof,point,element] = transpose(inverse)*GT.reference_tabulator(space)[dof,point]
            end
        end
    end
    nothing
end

@testset "Poisson numeric paths on CUDA" begin
    warps = (
        x->SVector(x[1]+0.1x[1]^2),
        x->SVector(1.5x[1]+0.3x[2]+0.1x[1]*x[2],
                   -0.2x[1]+1.2x[2]+0.05x[1]*x[2]),
        x->SVector(1.5x[1]+0.2x[2]+0.1x[2]*x[3],
                   -0.3x[1]+1.2x[2]+0.05x[1]*x[3],
                   0.1x[1]+0.2x[2]+1.7x[3]+0.04x[1]*x[2]),
    )
    for dim in 1:3,k in 1:3,orientation in (-1,1)
        base = GT.cartesian_mesh(ntuple(i->isodd(i) ? 0 : 1,2dim),
            ntuple(_->2,dim);simplexify=false)
        nodes = map(GT.node_coordinates(base)) do x
            y = warps[dim](x)
            setindex(y,orientation*y[1],1)
        end
        mesh = GT.replace_node_coordinates(base,nodes)
        omega = GT.interior(mesh)
        V = GT.lagrange_space(omega,k)
        GT.group_faces(mesh,dim)["poisson numeric subset"] = Int32[GT.num_faces(omega),1]
        domain = GT.interior(mesh;group_names=["poisson numeric subset"])
        Q = GT.quadrature(domain,2k)
        mesh_data = GT.mesh_integration_data(mesh,Q)
        space_data = GT.space_integration_data(ForwardDiff.gradient,V,Q)
        geometry_values = GT.space_integration_data(GT.value,GT.mesh_space(mesh,Val(dim)),Q)
        domain_data = GT.domain_integration_data(domain)
        prepared = adapt(CuArray,(mesh_data,space_data,geometry_values,domain_data))
        npoints = length(mesh_data.reference_weights)
        nelements = GT.num_faces(domain)
        ndofs = size(space_data.reference_tabulator,1)
        jacobian = GT.JacobianSumCallable{dim,dim,0.0,length(first(mesh_data.face_nodes)),false}()

        jacobians = CuArray{SMatrix{dim,dim,Float64,dim*dim}}(undef,npoints,nelements)
        inverses = similar(jacobians)
        weights = CUDA.zeros(Float64,npoints,nelements)
        gradients = CuArray{SVector{dim,Float64}}(undef,ndofs,npoints,nelements)
        coordinates = CuArray{SVector{dim,Float64}}(undef,npoints,nelements)
        CUDA.@sync @cuda threads=128 blocks=cld(length(jacobians),128) poisson_numeric_paths!(
            jacobians,inverses,weights,gradients,coordinates,prepared...,jacobian)

        expected_J = Array{eltype(jacobians)}(undef,npoints,nelements)
        expected_inv = similar(expected_J)
        expected_weights = zeros(npoints,nelements)
        expected_gradients = Array{eltype(gradients)}(undef,ndofs,npoints,nelements)
        expected_coordinates = Array{eltype(coordinates)}(undef,npoints,nelements)
        reference_points = GT.coordinates(only(GT.reference_quadratures(Q)))
        J_accessor = GT.jacobian_accessor(Q)
        weight_accessor = GT.weight_accessor(Q)
        gradient_accessor = GT.shape_function_accessor(ForwardDiff.gradient,V,Q)
        coordinate_accessor = GT.physical_map_accessor(GT.value,mesh,Val(dim))
        for (element,cell) in enumerate(GT.faces(domain))
            # Quadrature accessors take domain-local IDs; the untabulated mesh
            # coordinate accessor below takes a global mesh cell ID.
            J = J_accessor(element)
            w = weight_accessor(element)
            g = gradient_accessor(element)
            x = coordinate_accessor(cell)
            for point in 1:npoints
                expected_J[point,element] = J(point)
                expected_inv[point,element] = inv(J(point))
                expected_weights[point,element] = w(point)
                expected_coordinates[point,element] = x(reference_points[point])
                values = g(point)
                for dof in 1:ndofs
                    expected_gradients[dof,point,element] = values(dof)
                end
            end
        end
        @test Array(jacobians) ≈ expected_J rtol=2e-12 atol=1e-12
        @test Array(inverses) ≈ expected_inv rtol=2e-12 atol=1e-12
        @test Array(weights) ≈ expected_weights rtol=2e-12 atol=1e-12
        @test Array(gradients) ≈ expected_gradients rtol=2e-12 atol=1e-12
        @test Array(coordinates) ≈ expected_coordinates rtol=2e-12 atol=1e-12
    end
end

end # module

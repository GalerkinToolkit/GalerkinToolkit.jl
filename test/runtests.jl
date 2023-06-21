module GalerkinToolkitTests

import GalerkinToolkit as glk
using StaticArrays
using WriteVTK
using SparseArrays
using LinearAlgebra

function isoparametric_poisson(mesh)

    # Manufactured solution
    u(x) = sum(x)
    f(x) = 0.0

    # Dirichlet values
    node_to_tag = zeros(glk.num_nodes(mesh))
    tag_to_name = ["boundary"]
    glk.classify_mesh_nodes!(node_to_tag,mesh,tag_to_name)
    free_and_dirichlet_nodes = glk.partition_from_mask(i->i==0,node_to_tag)
    node_to_x = glk.node_coordinates(mesh)
    dirichlet_nodes = last(free_and_dirichlet_nodes)
    x_dirichlet = view(node_to_x,dirichlet_nodes)
    u_dirichlet = u.(x_dirichlet)

    # Reference cell
    d = glk.num_dims(mesh)
    ref_cells = glk.reference_faces(mesh,d)
    face_to_rid = glk.face_reference_id(mesh,d)

    # Integration rule
    degree_per_dir = ntuple(i->4,Val(d))
    integration_rules = map(ref_cells) do ref_cell
        glk.quadrature(glk.geometry(ref_cell),degree_per_dir)
    end
    w = map(glk.weights,integration_rules)
    q = map(glk.coordinates,integration_rules)

    # Shape functions
    ab = map(integration_rules,ref_cells) do integration_rule,ref_cell
        q = glk.coordinates(integration_rule)
        interpolation = glk.interpolation(ref_cell)
        shape_functions = glk.shape_functions(interpolation)
        m = length(q)
        n = glk.num_nodes(interpolation)
        a = glk.gradient!(shape_functions)(zeros(eltype(q),m,n),q)
        b = glk.value!(shape_functions)(zeros(m,n),q)
        (a,b)
    end
    ∇s = map(first,ab)
    s = map(last,ab)

    # Count coo entries
    n_coo = 0
    face_to_nodes = glk.face_nodes(mesh,d)
    ncells = glk.num_faces(mesh,d)
    node_to_free_node = glk.permutation(free_and_dirichlet_nodes)
    nfree = length(first(free_and_dirichlet_nodes))
    for cell in 1:ncells
        nodes = face_to_nodes[cell]
        for nj in nodes
            if node_to_free_node[nj] > nfree
                continue
            end
            for ni in nodes
                if node_to_free_node[ni] > nfree
                    continue
                end
                n_coo += 1
            end
        end
    end

    # Allocate coo values
    I_coo = zeros(Int32,n_coo)
    J_coo = zeros(Int32,n_coo)
    V_coo = zeros(Float64,n_coo)
    b = zeros(nfree)

    # Allocate auxiliary buffers
    ∇x = map(i->similar(i,size(i,2)),∇s)
    Aes = map(i->zeros(size(i,2),size(i,2)),∇s)
    fes = map(i->zeros(size(i,2)),∇s)
    ues = map(i->zeros(size(i,2)),∇s)

    Tx = SVector{d,Float64}
    TJ = typeof(zero(Tx)*zero(Tx)')

    # Fill coo values
    i_coo = 0
    for cell in 1:ncells
        rid = face_to_rid[cell]
        nodes = face_to_nodes[cell]
        Ae = Aes[rid]
        ue = ues[rid]
        fe = fes[rid]
        nl = length(nodes)
        ∇se = ∇s[rid]
        ∇xe = ∇x[rid]
        se = s[rid]
        we = w[rid]
        nq = length(we)
        fill!(Ae,zero(eltype(Ae)))
        fill!(fe,zero(eltype(Ae)))
        fill!(ue,zero(eltype(Ae)))
        for k in 1:nl
            nk = nodes[k]
            gk = node_to_free_node[nk]
            if  gk <= nfree
                continue
            end
            ue[k] = u_dirichlet[gk-nfree]
        end
        for iq in 1:nq
            Jt = zero(TJ) 
            xint = zero(Tx)
            for k in 1:nl
                x = node_to_x[nodes[k]]
                ∇sqx = ∇se[iq,k]
                sqx = se[iq,k]
                Jt += ∇sqx*x'
                xint += sqx*x
            end
            detJt = det(Jt)
            invJt = inv(Jt)
            dV = abs(detJt)*we[iq]
            for k in 1:nl
                ∇xe[k] = invJt*∇se[iq,k]
            end
            for j in 1:nl
                for i in 1:nl
                    Ae[i,j] += ∇xe[i]⋅∇xe[j]*dV
                end
            end
            fx = f(xint)
            for k in 1:nl
                fe[k] +=  fx*se[iq,k]*dV
            end
        end
        for i in 1:nl
            for j in 1:nl
                fe[i] -= Ae[i,j]*ue[j]
            end
        end
        for i in 1:nl
            ni = nodes[i]
            gi = node_to_free_node[ni]
            if gi > nfree
                continue
            end
            b[gi] += fe[i]
        end
        for j in 1:nl
            nj = nodes[j]
            gj = node_to_free_node[nj]
            if  gj > nfree
                continue
            end
            for i in 1:nl
                ni = nodes[i]
                gi = node_to_free_node[ni]
                if gi > nfree
                    continue
                end
                i_coo += 1
                I_coo[i_coo] = gi
                J_coo[i_coo] = gj
                V_coo[i_coo] = Ae[i,j]
            end
        end
    end

    A = sparse(I_coo,J_coo,V_coo,nfree,nfree)
    u_free = A\b
    nnodes = glk.num_nodes(mesh)
    uh = zeros(nnodes)
    uh[first(free_and_dirichlet_nodes)] = u_free
    uh[last(free_and_dirichlet_nodes)] = u_dirichlet

    vtk_grid("isoparametric_poisson",glk.vtk_args(mesh)...) do vtk
        glk.vtk_physical_groups!(vtk,mesh)
        vtk["tag"] = node_to_tag
        vtk["uh"] = uh
    end

end

vtk_mesh_cell = nodes -> WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_VERTEX,nodes)
vertex =(;vtk_mesh_cell)
vtk_mesh_cell = nodes -> WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_LINE,nodes)
segment = (;vtk_mesh_cell)
vtk_mesh_cell = nodes -> WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_QUAD,nodes[[1,2,4,3]])
geometry = (;num_dims = Val(2), is_n_cube=true, is_axis_aligned=true, bounding_box=SVector{2,Float64}[(0,0),(1,1)])
order_per_dir = (1,1)
interpolation = glk.lagrange_interpolation(geometry,order_per_dir)
quad = (;geometry,vtk_mesh_cell,interpolation)
node_coordinates = SVector{2,Float64}[(0,0),(1,0),(2,0),(0,1),(1,1),(2,1),(0,2),(1,2),(2,2)]
face_nodes = [
   Vector{Int}[],
   [[1,2],[2,3],[1,4],[4,7],[7,8],[8,9],[3,6],[6,9]],
   [[1,2,4,5],[2,3,5,6],[4,5,7,8],[5,6,8,9]]
  ]
face_reference_id = [Int[],Int[1,1,1,1,1,1,1,1],[1,1,1,1]]
reference_faces = ([vertex],[segment],[quad])
physical_groups = [
  [],
  ["face_1"=>[1,2],"face_2"=>[3,4],"boundary"=>[1,2,3,4,5,6,7,8]],
  ["domain"=>[1,2,3,4]]]
mesh = (;
    num_dims=Val(2),node_coordinates,
    face_nodes,face_reference_id,
    reference_faces,physical_groups)

d = 2
vtk_grid("mesh2",glk.vtk_args(mesh,d)...) do vtk
    glk.vtk_physical_groups!(vtk,mesh,d)
end

d = 1
vtk_grid("mesh1",glk.vtk_args(mesh,d)...) do vtk
    glk.vtk_physical_groups!(vtk,mesh,d)
end

vtk_grid("mesh",glk.vtk_args(mesh)...) do vtk
    glk.vtk_physical_groups!(vtk,mesh)
end


isoparametric_poisson(mesh)




end # module

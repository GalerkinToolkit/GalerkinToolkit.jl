module GalerkinToolkitTests

import GalerkinToolkit as GK
using StaticArrays
using WriteVTK
using SparseArrays

function simple_poisson(mesh)

    node_to_tag = zeros(GK.num_nodes(mesh))
    tag_to_name = ["boundary"]
    GK.classify_mesh_nodes!(node_to_tag,mesh,tag_to_name)

    free_and_dirichlet_nodes = GK.partition_from_mask(i->i==0,node_to_tag)
    face_to_nodes = GK.face_nodes(mesh)

    d = GK.num_dims(mesh)
    ref_cells = GK.reference_faces(mesh,d)

    degree_per_dir = ntuple(i->4,Val(d))
    integration_rules = map(ref_cells) do ref_cell
        GK.quadrature(GK.geometry(ref_cell),degree_per_dir)
    end

    #∇s = map(integration_rules,ref_cells) do integration_rule,ref_cell
    #    f! = ref_cell |> GK.interpolation |> GK.shape_functions |> GK.gradient
    #    q = GK.coordinates(integration_rule)
    #    nlnodes = ref_cell |> GK.interpolation |> GK.node_coordinates |> length
    #    f!(,)
    #end

    vtk_grid("simple_poisson",GK.vtk_args(mesh)...) do vtk
        GK.vtk_physical_groups!(vtk,mesh)
        vtk["tag"] = node_to_tag
    end

end

vtk_mesh_cell = nodes -> WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_VERTEX,nodes)
vertex =(;vtk_mesh_cell)
vtk_mesh_cell = nodes -> WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_LINE,nodes)
segment = (;vtk_mesh_cell)
vtk_mesh_cell = nodes -> WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_QUAD,nodes[[1,2,4,3]])
geometry = (;num_dims = Val(2), is_n_cube=true, is_axis_aligned=true, bounding_box=SVector{2,Float64}[(0,0),(1,1)])
quad = (;geometry,vtk_mesh_cell)
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
vtk_grid("mesh2",GK.vtk_args(mesh,d)...) do vtk
    GK.vtk_physical_groups!(vtk,mesh,d)
end

d = 1
vtk_grid("mesh1",GK.vtk_args(mesh,d)...) do vtk
    GK.vtk_physical_groups!(vtk,mesh,d)
end

vtk_grid("mesh",GK.vtk_args(mesh)...) do vtk
    GK.vtk_physical_groups!(vtk,mesh)
end


simple_poisson(mesh)




end # module

module TestWriteVtk

using GalerkinToolkit
using Meshes
using Test
using WriteVTK

hex = Hexahedron(Point.([(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)]))

for d in 0:dimension(hex)
    fn = "hex_$d"
    vtk_grid(fn,vtk_args(hex,d)...) do vtk end
end

tet = Tetrahedron(Point.([(0,0,0),(1,0,0),(0,1,0),(0,0,1)]))

for d in 0:dimension(tet)
    fn = "tet_$d"
    vtk_grid(fn,vtk_args(tet,d)...) do vtk end
end

groups_0 = Dict(1=>physical_group([1,2,3],"myvertices"))
groups_1 = Dict(1=>physical_group([2,6,5],"myedges"))
groups_2 = Dict(1=>physical_group([2,1],"myfaces"))
groups_3 = Dict(1=>physical_group([1],"myvolume"))
groups = [groups_0,groups_1,groups_2,groups_3]

for d in 0:dimension(tet)
    fn = "tet_with_groups_$d"
    vtk_grid(fn,vtk_args(tet,d)...) do vtk
        vtk_physical_groups!(vtk,d,tet,groups[d+1])
    end
end

end # module

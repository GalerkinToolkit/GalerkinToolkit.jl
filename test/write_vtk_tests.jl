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

groups_0 = ["myvertices"=>[1,2,3]]
groups_1 = ["myedges"=>[2,6,5]]
groups_2 = ["myfaces"=>[2,1]]
groups_3 = ["myvolume"=>[1]]
groups = [groups_0,groups_1,groups_2,groups_3]

for d in 0:dimension(tet)
    fn = "tet_with_groups_$d"
    vtk_grid(fn,vtk_args(tet,d)...) do vtk
        vtk_physical_groups!(vtk,tet,d,physical_groups=groups[d+1])
    end
end

vtk_grid("tet_with_groups",vtk_args(tet)...) do vtk
    vtk_physical_groups!(vtk,tet,physical_groups=groups)
end

end # module

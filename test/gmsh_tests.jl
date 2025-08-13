module GmshTests

import GalerkinToolkit as GT

outdir = mkpath(joinpath(@__DIR__,"..","output"))
msh =  joinpath(@__DIR__,"..","assets","quad.msh")
mesh = GT.mesh_from_msh(msh)

mesh = GT.mesh_from_msh(msh,complexify=false)
mesh = GT.complexify(mesh)

msh =  joinpath(@__DIR__,"..","assets","demo.msh")
mesh = GT.mesh_from_msh(msh)


# 2D & 3D, cube & simplex 
for dim in 2:3
    for recombine in [true, false]
        mesh = GT.with_gmsh() do gmsh
            lc = 1.0
            gmsh.model.geo.addPoint(0, 0, 0, lc, 1)
            gmsh.model.geo.addPoint(1, 0,  0, lc, 2)
            gmsh.model.geo.addPoint(1, 1, 0, lc, 3)
            gmsh.model.geo.addPoint(0, 1, 0, lc, 4)
            gmsh.model.geo.addLine(1, 2, 1)
            gmsh.model.geo.addLine(2, 3, 2)
            gmsh.model.geo.addLine(3, 4, 3)
            gmsh.model.geo.addLine(4, 1, 4)
            gmsh.model.geo.addCurveLoop([4, 1, 2, 3], 1)
            gmsh.model.geo.addPlaneSurface([1], 1)
            
            gmsh.model.geo.mesh.setAlgorithm(2, 1, 8)
            if recombine
                gmsh.model.geo.mesh.setRecombine(2, 1)
            end
            if dim == 3
                gmsh.model.geo.extrude([(2, 1)], 0, 0, 1, [1/lc], [], recombine) # recombine to a single cube
            end
            gmsh.model.geo.synchronize()

            dimtag = gmsh.model.getEntities()
            entities = Dict([(i => []) for i in 0:3])
            for (dim, tag) in dimtag
                push!(entities[dim], tag)
            end
            
            gmsh.model.addPhysicalGroup(0, entities[0], 5)
            gmsh.model.addPhysicalGroup(1, entities[1], 6)
            gmsh.model.addPhysicalGroup(2, entities[2], 7)
            gmsh.model.setPhysicalName(0, 5, "boundary")
            gmsh.model.setPhysicalName(1, 6, "boundary")
            

            if dim == 3
                gmsh.model.addPhysicalGroup(3, entities[3], 8)
                gmsh.model.setPhysicalName(2, 7, "boundary")
                gmsh.model.setPhysicalName(3, 8, "volume")
            else
                gmsh.model.setPhysicalName(2, 7, "volume")
            end
            gmsh.model.mesh.generate(dim)
            GT.mesh_from_gmsh(gmsh)
        end
    end 
end



end # module

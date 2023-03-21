module P4estTests

using GalerkinToolkit
using Meshes
using StaticArrays
using Test

coarse_mesh = Quadrangle(Point.([(0,0),(2,0),(2,2),(0,2)]))

initial_level = 1
forest = forest_from_mesh(coarse_mesh,initial_level)

@show typeof(forest)

x = zeros(SVector{2,Float64},4)
for (itree,tree) in enumerate(forest)
    for leaf in tree
        anchor(leaf)
        level(leaf)
        node_coordinates!(x,forest,itree,leaf)
    end
end

recursive = true
ileaf = Ref(0)
refine!(forest,recursive) do itree,leaf
    ileaf[] += 1
    @show anchor(leaf)
    node_coordinates!(x,forest,itree,leaf)
    level(leaf) < 2 ? 1 : 0
end

@test length(first(forest)) == 16

recursive = false
coarsen!(forest,recursive) do itree,leafs
    @show typeof(leafs)
    ileaf[] += 1
    @test length(leafs) == 4
    for leaf in leafs
        node_coordinates!(x,forest,itree,leaf)
        @show x
    end
    true
end

refine!(forest,true) do itree,leaf
    anchor(leaf) == SVector(0,0) && level(leaf) < 4
end
balance!(forest)

partition!(forest)

allow_for_coarsening = true
partition!(forest,allow_for_coarsening) do itree,leaf
    level(leaf)
end

ghost_leafs = find_ghost_leafs(forest)
@test length(ghost_leafs) == 0

order = 1
forest_nodes = generate_nodes(forest,order)
forest_nodes = generate_nodes(forest,order,ghost_leafs)

leaf_to_nodes = leaf_nodes(forest_nodes)
display(leaf_to_nodes)

leaf_to_constraints = leaf_constraints(forest_nodes)

display(leaf_to_constraints)


end # module

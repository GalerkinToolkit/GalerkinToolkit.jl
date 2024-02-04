module Example002Tests

import GalerkinToolkit as gk
using GalerkinToolkitExamples: Example002
using Test
using PartitionedArrays

tol = 1.0e-10

params = Dict{Symbol,Any}()
results = Example002.main(params)
@test results[:eh1] < tol
@test results[:el2] < tol

domain = (0,10,0,10)
cells_per_dir = (20,20)
parts_per_dir = (2,2)
np = prod(parts_per_dir)
parts = DebugArray(LinearIndices((np,)))
ghost_layers = 0
mesh = gk.cartesian_mesh(domain,cells_per_dir,parts_per_dir;parts,ghost_layers)

params = Dict{Symbol,Any}()
results = Example002.main(params)
@test results[:eh1] < tol
@test results[:el2] < tol

end # module

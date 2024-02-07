module Example002Tests

import GalerkinToolkit as gk
using GalerkinToolkitExamples: Example002
using Test
using PartitionedArrays
using PetscCall

tol = 1.0e-10

params = Dict{Symbol,Any}()
results = Example002.main(params)
@test results[:eh1] < tol
@test results[:el2] < tol

params = Dict{Symbol,Any}()
domain = (0,10,0,10)
cells_per_dir = (20,20)
parts_per_dir = (2,2)
np = prod(parts_per_dir)
parts = DebugArray(LinearIndices((np,)))
ghost_layers = 0
mesh = gk.cartesian_mesh(domain,cells_per_dir,parts_per_dir;parts,ghost_layers)
params[:mesh] = mesh
results = Example002.main(params)
@test results[:eh1] < tol
@test results[:el2] < tol

options = "-pc_type gamg -ksp_type cg -ksp_error_if_not_converged true -ksp_converged_reason -ksp_rtol 1.0e-6"
PetscCall.init(args=split(options))

params = Dict{Symbol,Any}()
domain = (0,10,0,10,0,10)
cells_per_dir = (40,40,40)
parts_per_dir = (2,2,2)
np = prod(parts_per_dir)
parts = DebugArray(LinearIndices((np,)))
ghost_layers = 0
mesh = gk.cartesian_mesh(domain,cells_per_dir,parts_per_dir;parts,ghost_layers)
params[:mesh] = mesh
params[:export_vtu] = false
params[:solver] = Example002.ksp_solver()
results = Example002.main(params)
results = Example002.main(params)

end # module

module PoissonParallelTests

using Test
import GalerkinToolkit as gt
using PartitionedArrays
using Metis

include("../examples/poisson_bddc.jl")

using .PoissonBDDC: CORNER,FACE

with_debug() do distribute

    tol = 1.0e-7
    params = Dict{Symbol,Any}()
    mesh = gt.cartesian_mesh((0,2,0,2),(4,4),complexify=false)
    np = 4
    ranks = distribute(LinearIndices((np,)))
    ranks = LinearIndices((np,))
    pmesh = gt.partition_mesh(Metis.partition,ranks,mesh,via=:cells)
    params[:pmesh] = pmesh
    params[:dirichlet_tags] = ["1-face-1","1-face-3","1-face-4"]
    params[:neumann_tags] = ["1-face-2"]
    params[:u] = (x) -> sum(x)
    params[:f] = (x) -> 0.0
    params[:g] = (x) -> 1.0
    params[:bddc_type] = (CORNER,FACE,)
    @time results = PoissonBDDC.main(params)
    @test results[:eh1] < tol
    @test results[:el2] < tol

end


end # module

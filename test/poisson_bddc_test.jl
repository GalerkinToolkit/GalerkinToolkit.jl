module PoissonParallelTests

using Test
import GalerkinToolkit as gt
using PartitionedArrays
using Metis

include("../examples/poisson_bddc.jl")

with_debug() do distribute

    #tol = 1.0e-7
    #np = 8
    #load = 20
    #n = ceil(Int,load*(np^(1/3)))

    params = Dict{Symbol,Any}()
    n = 10
    np = 4 #MPI.Comm_size(MPI.COMM_WORLD)
    tol = 1.0e-6
    mesh = gt.cartesian_mesh((0,2,0,2,0,2),(n,n,n),complexify=false)
    ranks = distribute(LinearIndices((np,)))
    pmesh = gt.partition_mesh(Metis.partition,ranks,mesh,via=:cells)
    params[:pmesh] = pmesh
    #params[:dirichlet_tags] = ["1-face-1","1-face-3","1-face-4"]
    #params[:neumann_tags] = ["1-face-2"]
    params[:dirichlet_tags] = ["boundary"]
    params[:u] = (x) -> sum(x)
    params[:f] = (x) -> 0.0
    params[:g] = (x) -> 1.0
    @time results = PoissonBDDC.main(params)
    @test results[:eh1] < tol
    @test results[:el2] < tol

end


end # module

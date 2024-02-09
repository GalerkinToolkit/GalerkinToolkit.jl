module Example004Tests

using MPI

include("run_mpi_driver.jl")

file = joinpath(@__DIR__,"example004_driver_np_4.jl")
run_mpi_driver(file;procs=4)

end # module


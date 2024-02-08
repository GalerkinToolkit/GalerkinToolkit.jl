module Example002Tests

using MPI

include("run_mpi_driver.jl")

file = joinpath(@__DIR__,"example002_driver_np_4.jl")
run_mpi_driver(file;procs=4)

end # module


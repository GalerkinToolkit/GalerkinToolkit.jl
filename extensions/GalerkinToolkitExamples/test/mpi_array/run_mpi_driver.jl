using MPI
using Test
function run_mpi_driver(file;procs)
  repodir = normpath(joinpath(@__DIR__,"..",".."))
  mpiexec() do cmd
    if MPI.MPI_LIBRARY == "OpenMPI" || (isdefined(MPI, :OpenMPI) && MPI.MPI_LIBRARY == MPI.OpenMPI)
      run(`$cmd -n $procs --oversubscribe $(Base.julia_cmd()) --project=$repodir $(file)`)
    else
      run(`$cmd -n $procs $(Base.julia_cmd()) --project=$repodir $(file)`)
    end
    # This line will be reached if and only if the command launched by `run` runs without errors.
    # Then, if we arrive here, the test has succeeded.
    @test true  
  end
end

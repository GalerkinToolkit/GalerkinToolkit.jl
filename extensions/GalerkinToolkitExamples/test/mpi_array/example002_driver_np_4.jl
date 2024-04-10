module TMP

using PartitionedArrays

include(joinpath("..","example002_defs.jl"))
with_mpi(example002_tests_np_4)
with_mpi(example002_advanced_tests_np_4)

end # module

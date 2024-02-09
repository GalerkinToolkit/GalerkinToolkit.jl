module TMP

using PartitionedArrays

include(joinpath("..","example004_defs.jl"))
with_mpi(example004_tests_np_4)

end # module

# Element CSC CUDA correctness checks

Run these commands from the repository root on a machine with a functional CUDA
GPU:

```sh
julia --project=test/cuda test/cuda/instantiate.jl
julia --project=test/cuda test/cuda/runtests.jl
```

The setup script develops this checkout into the dedicated CUDA test environment
and installs its dependencies. The test script uses the environment selected by
Julia's `--project` argument. It can also be invoked with absolute paths from a
different working directory. Select a GPU with `CUDA_VISIBLE_DEVICES` if needed.

These checks cover single-field CG matrix storage: adaptation, both lookup
layouts, parallel atomic scatter, reset, sparse output, and reordered volume
subsets. Scalar indexing of device arrays is disabled. Missing GPU hardware fails
the suite. The ordinary CPU suite does not load CUDA. FEM integration and its
compiler data are prepared on the CPU; the complete integration kernel is a
later step.

Layout and placement are independent:

```julia
using Adapt, CUDA
import GalerkinToolkit as GT
device = adapt(CuArray, GT.to_layout(host, GT.EntryMajorLayout()))
argument = GT.element_csc_domain(device, 1)
A = GT.compress(device)
```

The entry-major table has physical storage `[element, local_entry]` and retains
logical indexing `[local_entry, element]` through a `Transpose` wrapper. Select
domain arguments after upload so that they share the parent's values buffer.
Adaptation preserves floating-point and integer element types.

The internal atomic contribution primitive accepts only numeric buffers and
prepared indices and returns `nothing`:

```julia
GT.contribute_element_csc_block_atomic!(nzval, positions, values, element, base, (nrows, ncols))
```

Its caller must skip disabled blocks (`base == -1`). All contribution methods
return `nothing`. Both paths are inlined for device code; GPU kernels can call
the buffer-only primitive directly. The tests exercise both paths, including
constrained entries and padded rectangular blocks. Timing and compiler
experiments are kept outside this correctness suite.

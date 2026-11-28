# CUDA assembly and prepared-data correctness checks

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
the suite. The ordinary CPU suite does not load CUDA.

Select `GT.CUDAExecution(threads=128)` and element CSC storage to compile the
variational form for one GPU.

```julia
using CUDA, ForwardDiff, LinearAlgebra
import GalerkinToolkit as GT
Ω = GT.interior(mesh)
k = 2
V = GT.lagrange_space(Ω, k)
dΩ = GT.measure(Ω, 2k)
a = (u,v) -> GT.∫(x -> dot(ForwardDiff.gradient(u,x), ForwardDiff.gradient(v,x)), dΩ)
A, cache = GT.assemble_matrix(a, Float64, V, V;
    execution=GT.CUDAExecution(), reuse=Val(true))
GT.update_matrix!(A, cache; parameters=())
```

The generator uses the same form optimization and worker-range lowering as CPU
assembly. It prepares geometry and reference data on the host, uploads numeric
inputs for each contribution during setup, and constructs statically sized
scratch within each GPU thread. The GPU parent allocation is uploaded once in
element-major layout. All domain arguments share its numerical values.

On update, supply the same number of parameters as during preparation. Scalars
and explicitly adapted GPU arrays and structures are supported. The backend
checks parameter compatibility before resetting the matrix. Scratch dimensions
are determined during preparation. Prepared inputs are a snapshot. Rebuild the cache
after changing CPU-side geometry, spaces, quadrature, or table contents. Matrix
values and prepared arrays remain in the cache across updates. Reset, contribution
launches, and compression/copies use the current task's CUDA stream, with one
stream synchronization at the public assembly/update boundary. CPU worker-local
allocation and its workspace factory are unchanged.

See `execution_tests.jl` for generated Poisson and mass-matrix checks. The suite
covers CG orders 1–3 on tensor-product and simplex meshes, overlapping domains,
parameter updates, and output reuse. The rectangular variable-coefficient case
uses an independent quadrature reference. Complex matrices (`ComplexF32`, `ComplexF64`)
are checked against CPU assembly; their entries are scattered as two real atomic adds.

Prepared mesh, space, domain, and DOF data can also be transferred directly:

```julia
using Adapt, CUDA, ForwardDiff
import GalerkinToolkit as GT
prepared = (
    GT.mesh_integration_data(mesh, Q),
    GT.space_integration_data(ForwardDiff.gradient, V, Q),
    GT.domain_integration_data(GT.domain(Q)),
    GT.DOFData(V),
)
device_prepared = adapt(CuArray, prepared)
host_prepared = adapt(Array, device_prepared)
```

Adaptation preserves element-major storage, array dimensions, scalar/index types,
and packed static vectors. It performs no permutation or scalar-type conversion.
Reference tabulations and weights use ordinary device arrays; they are not placed
in CUDA constant memory. Preparing the arrays and tabulating basis functions still
happens on the CPU. The tests read all prepared buffers in parallel GPU kernels
and check reordered domain-to-cell/node/DOF addressing for scalar CG spaces of
orders 1–3 in dimensions 1–3. Generated matrix assembly currently targets
single-field CG volume forms; vector assembly and other integration domains
remain separate work.

Poisson numeric checks also evaluate coordinates, Jacobians, determinant/adjugate
inverses, physical quadrature weights, and inverse-transpose basis gradients on
the GPU. They compare against existing CPU accessors on sheared and non-affine
cells, with reversed orientation and selected domains. Quadrature accessors use
domain-local IDs; prepared geometry uses global mesh IDs after domain lookup.

Constant-memory placement needs a per-table storage choice: the Float64 gradient
table for a cubic hexahedron has 64 DOFs and 64 quadrature points, occupying
96 KiB, beyond CUDA's 64 KiB constant-memory capacity. Smaller tables may benefit
when threads in a warp read the same entry. `CUDA.Const` marks a device array
read-only; it does not move that array into constant memory. See the
[CUDA.jl device-array API](https://cuda.juliagpu.org/stable/api/kernel/) and
[NVIDIA's constant-memory guidance](https://docs.nvidia.com/cuda/cuda-c-best-practices-guide/#constant-memory).

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

```@meta
CurrentModule = GalerkinToolkit
```

# Home

Welcome to the documentation page of GalerkinToolkit!


## What

GalerkinToolkit provides a collection of tools to build computer codes to solve partial differential equations (PDEs)
using finite element (FE) methods on different computing platforms, from laptops to modern GPU-based supercomputers.
GalerkingToolkit is fully implemented in the [Julia programming language](https://julialang.org/) and provides (or it will soon provide) tools for:

- Reading and partitioning computational meshes generated from external mesh generators ([`Gmsh`](https://gmsh.info/) at this moment).
- Defining discrete interpolation spaces on computational meshes including continuous and discontinuous interpolations.
- Integrating 0-, 1-, and 2-forms on triangulated manifolds for a variety of problems, including scalar- and vector-valued equations, and single- and multi-field systems.
- Discretizing PDEs into systems of linear-, nonlinear-, and differential-algebraic equations, including matrix-free approaches.
- Representing such algebraic problems in a way that can be readily solved using external tools such as [`PartitionedSolvers.jl`](https://github.com/PartitionedArrays/PartitionedArrays.jl), [`PetscCall.jl`](https://github.com/PartitionedArrays/PetscCall.jl), [`LinearSolve.jl`](https://github.com/SciML/LinearSolve.jl), [`NonLinearSolve.jl`](https://github.com/SciML/NonlinearSolve.jl), and [`DifferentialEquations.jl`](https://github.com/SciML/DifferentialEquations.jl).
- Automatic differentiation for non-linear and parametric problems.
- Visualizing results with [`Makie.jl`](https://github.com/MakieOrg/Makie.jl) and [`Paraview`](https://www.paraview.org/).

## Why

GalerkinToolkit is definitively not the first FE software project out there, but it has some unique design features:

- It combines the vision of frameworks like [FEniCS](https://fenicsproject.org/) and libraries like [Deal-ii](https://www.dealii.org/) in a single package. It provides a high-level API backed with automatic julia-to-julia code generation, and also easy access to the low-level building blocks of FE codes.
- It is designed to blend with the Julia package ecosystem, and reuses as much as possible from it. For instance, GalerkinToolkit does not define differential operators, it simply uses the operators already defined in [`ForwardDiff.jl`](https://github.com/JuliaDiff/ForwardDiff.jl). It leverages [`PartitionedArrays.jl`](https://github.com/PartitionedArrays/PartitionedArrays.jl) for distributed linear algebra data structures. You can use [`StaticArrays.jl`](https://github.com/JuliaArrays/StaticArrays.jl) and [`Tensors.jl`](https://github.com/Ferrite-FEM/Tensors.jl), when working with vector- and tensor- valued PDEs. It also makes straight-forward to use external solvers to solve the algebraic systems you get after discretizing a PDE. For visualization, it defines recipes for `Makie.jl` and defines helper functions to save results with `WriteVTK.jl` in vtk format.
- It is based on a new form compiler, the GalerkinToolkit form compiler (GTFC), that fixes the rigidity of domain-specific languages like UFL by introducing an alternative multi-level intermediate representation (MLIR) and leveraging the meta-programming features of the Julia programming language. 



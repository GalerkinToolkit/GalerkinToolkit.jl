# Introduction

## Features

**GalerkinToolkit** is a high-performance finite element toolbox built in the [Julia programming language](https://julialang.org/). Its mission is to provide general-purpose building blocks for finite element (FE) methods that can be flexibly combined to solve a wide variety of partial differential equations (PDEs), using diverse numerical schemes, and across a range of computing platforms -- from laptops to supercomputers.

The toolkit envisions a unified framework that addresses the needs of numerical analysts, domain scientists, and high-performance computing experts. To this end, it offers a rich API with multiple levels of abstraction:

- **High-level API**: Enables users to solve PDEs using mathematical abstractions to define PDEs in weak form, while hiding low-level implementation details.
- **Low-level API**: Grants direct access to underlying numerical quantities, allowing custom implementations of numerical schemes or advanced code optimizations not available through the high-level interface.

GalerkinToolkit provides tools for:

- Reading and partitioning computational meshes generated with external mesh generators.
- Defining discrete interpolation spaces on triangulated manifolds, supporting both continuous and discontinuous interpolations.
- Integrating functions, linear forms, and bilinear forms for scalar- and vector-valued equations, including single- and multi-field systems.
- Supporting both simplex and hypercube cell geometries, with high-order interpolation capabilities.
- Discretizing a wide range of PDEs into systems of linear, nonlinear, or differential-algebraic equations.
- Representing algebraic problems in formats compatible with external solvers.
- Visualization and post-process of results.

Even though not currently available, these other features are being implemented, or are planned to be implemented in the near future:

- Automatic differentiation for non-linear PDEs and gradient-based optimization.
- Matrix-free bilinear forms.
- H(div) and H(curl) interpolation spaces.
- Distributed assembly.
- Single- and multi-GPU support.


## Novelties

GalerkinToolkit is certainly not the first FEM software project, but it introduces several distinctive design novelties:

- **Unified high- and low-level APIs** --  It combines the vision of frameworks like [FEniCS](https://fenicsproject.org/) and libraries like [Deal-ii](https://www.dealii.org/) in a single package and in a single programming language, Julia. This addresses the two-language problem of previous FE projects that consider a Python front-end for easiness of use and a C/C++ backend for performance. With GalerkinToolkit, you can use a concise high-level syntax or directly implement integration loops using low-level building blocks, depending on your needs.

- **Deep integration with the Julia ecosystem**  --
  Rather than reinventing the wheel, GalerkinToolkit reuses existing Julia packages in numerous situations. For example:
  - Computational meshes generated with [`Gmsh`](https://gmsh.info/).
  - Differential operators from [`ForwardDiff.jl`](https://github.com/JuliaDiff/ForwardDiff.jl)
  - Vector/tensor types via [`StaticArrays.jl`](https://github.com/JuliaArrays/StaticArrays.jl) and [`Tensors.jl`](https://github.com/Ferrite-FEM/Tensors.jl)
  - External solvers for algebraic systems from [`PartitionedSolvers.jl`](https://github.com/PartitionedArrays/PartitionedArrays.jl), [`PetscCall.jl`](https://github.com/PartitionedArrays/PetscCall.jl), [`LinearSolve.jl`](https://github.com/SciML/LinearSolve.jl), [`NonLinearSolve.jl`](https://github.com/SciML/NonlinearSolve.jl), and [`DifferentialEquations.jl`](https://github.com/SciML/DifferentialEquations.jl).
  - Visualization with [`Makie.jl`](https://github.com/MakieOrg/Makie.jl) and [`WriteVTK.jl`](https://github.com/JuliaVTK/WriteVTK.jl)

- **A new form compiler: GTFC**  --
  The GalerkinToolkit Form Compiler (GTFC) uses Julia’s metaprogramming capabilities to generate efficient code from weak form definitions written in pure Julia. Unlike other compilers like FFC or TSFC, GTFC:
  - Does not rely on an external DSL like UFL. It considers a sub-set of the Julia programming language as alternative.
  - Allows user-defined types in the weak form.
  - Supports advanced use cases, such as coupling surface and volume arguments in multi-field weak forms both for continuous and discontinuous interpolations.

- **Gridap, reimagined**  --
  GalerkinToolkit began as a full reimplementation of the core ideas behind [Gridap](https://github.com/gridap/Gridap.jl). While Gridap is based on lazily mapped arrays, GalerkinToolkit centers on form compilation and functions designed to easily deal with quantities at integration points. These provide both a loop-free high-level API and manual control of integration loops.


## This manual

In next sections of the manual, we provide explanations on the different parts of the library. At the moment available manual sections are:

- Documentation [For developers](@ref).
















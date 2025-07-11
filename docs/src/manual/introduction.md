# Introduction

**GalerkinToolkit** is a high-performance finite element toolbox fully implemented in the [Julia programming language](https://julialang.org/). Its mission is to provide general-purpose building blocks of finite element (FE) methods to solve a wide variety of partial differential equations (PDEs), using diverse numerical schemes, and across a range of computing platforms, from laptops to supercomputers.

The toolkit envisions a unified framework that addresses the needs of numerical analysts, domain scientists, and high-performance computing experts; offering a rich API with multiple levels of abstraction:

- **High-level API**: Enables users to solve PDEs conveniently using mathematical abstractions instead of implementing low-level code by hand.
- **Low-level API**: Grants direct access to underlying numerical quantities, allowing custom implementations of numerical schemes or advanced code optimizations not available through the high-level interface.

## Features

GalerkinToolkit provides tools for:

- Reading and partitioning computational meshes generated with external mesh generators.
- Defining discrete interpolation spaces on triangulated manifolds, supporting both continuous and discontinuous interpolations.
- Integrating functions, linear, and bilinear forms for scalar- and vector-valued equations, including single- and multi-field systems.
- Supporting both simplex and hypercube face geometries, with high-order interpolation capabilities.
- Discretizing a wide range of PDEs into systems of linear, nonlinear, and differential-algebraic equations.
- Representing algebraic problems in formats compatible with external solvers.
- Visualization and post-process of results.

Even though not currently available, these other features are being implemented, or are planned to be implemented in the near future:

- Automatic differentiation for non-linear PDEs and gradient-based optimization methods.
- Matrix-free bilinear forms.
- H(div) and H(curl) interpolation spaces.
- Distributed assembly.
- Single- and multi-GPU support.


## Novelties

GalerkinToolkit is certainly not the first FEM software project, but it introduces it introduces a novel design:

- **Unified high- and low-level APIs**--  It combines the vision of frameworks like [FEniCS](https://fenicsproject.org/) and libraries like [Deal-ii](https://www.dealii.org/) in a single package and in a single programming language. It addresses the two-language problem of previous FE projects that consider a Python front-end for easiness of use and a C/C++ backend for performance. With GalerkinToolkit, you can use a concise high-level syntax or directly implement integration loops using low-level building blocks, depending on your needs.

- **Deep integration with the Julia ecosystem**--
  GalerkinToolkit reuses existing Julia packages in numerous situations. For example:
  - Computational meshes generated with [`Gmsh`](https://gmsh.info/).
  - Differential operators from [`ForwardDiff.jl`](https://github.com/JuliaDiff/ForwardDiff.jl)
  - Vector/tensor types via [`StaticArrays.jl`](https://github.com/JuliaArrays/StaticArrays.jl) and [`Tensors.jl`](https://github.com/Ferrite-FEM/Tensors.jl)
  - External solvers for algebraic systems from [`PartitionedSolvers.jl`](https://github.com/PartitionedArrays/PartitionedArrays.jl), [`PetscCall.jl`](https://github.com/PartitionedArrays/PetscCall.jl), [`LinearSolve.jl`](https://github.com/SciML/LinearSolve.jl), [`NonLinearSolve.jl`](https://github.com/SciML/NonlinearSolve.jl), and [`DifferentialEquations.jl`](https://github.com/SciML/DifferentialEquations.jl).
  - Visualization with [`Makie.jl`](https://github.com/MakieOrg/Makie.jl) and [`WriteVTK.jl`](https://github.com/JuliaVTK/WriteVTK.jl)

- **The GalerkinToolkit Form Compiler (GTFC)**--
  GTFC uses Julia’s metaprogramming capabilities to generate efficient code from high-level mathematical abstractions. Unlike other compilers like FFC or TSFC:
  - It does not rely on an external DSL like UFL. It considers a sub-set of the Julia programming language as alternative.
  - In consequence, it allows user-defined types in the weak form.
  - Moreover, it supports advanced use cases, such as coupling surface and volume arguments in multi-field weak forms both for continuous and discontinuous interpolations.

- **A major overhaul of Gridap**--
  GalerkinToolkit began as a full reimplementation of the core ideas behind [Gridap](https://github.com/gridap/Gridap.jl). While Gridap is based on lazily mapped arrays, GalerkinToolkit centers on form compilation and a low-level API designed to easily deal with quantities at integration points. These provide both a loop-free high-level API and manual control of integration loops.


## Installation

GalerkinToolkit is a registered package in the [official Julia package registry](https://github.com/JuliaRegistries/General). As such, you can install GalerkinToolkit easily using the [Julia package manager](https://docs.julialang.org/en/v1/stdlib/Pkg/).

Open Julia and type `]` to enter package mode.

```
julia> ]
```

Then, type
```
pkg> add GalerkinToolkit
```

This installs GalerkinToolkit and all its dependencies.

Press `ctrl+C` to go back to standard mode. Now, you can type

```
julia> import GalerkinToolkit as GT
```
to start using GalerkinToolkit.

You will need to have Julia installed in your system. See the official Julia installation instructions [here](https://julialang.org/install/).


## This manual

It provides detailed explanations so that users build an understanding of the library. We start detailing the geometrical foundations of the library: computational meshes and domains. Then, we cover the construction of interpolation spaces and the numerical integration capabilities. Finally, we explain the post-processing layer including the visualization of results.



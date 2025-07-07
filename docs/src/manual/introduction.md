# Introduction

## Overview

GalerkinToolkit is a high-performance finite element toolbox implemented in the [Julia programming language](https://julialang.org/). Its mission is to provide general-purpose building blocks of finite element (FE) methods that can be combined in multiple ways to solve a wide range of partial differential equations (PDEs), using different numerical schemes, and on different computing systems, from laptops to supercomputers. Its vision is to provide a unified framework that can address the needs of numerical analysts, domain scientists, and high-performance computing experts. To this end, GalerkinToolkit provides different levels of abstractions with a clear separation of concerns. A high-level API allows one to define complex numerical schemes using mathematical abstractions, conveniently hiding most of the implementation details. A low-level API is also available providing direct access to the underlying numerical quantities and making possible to implement schemes not available via the high-level API. In the current version, GalerkinToolkit provides tools for:

- Reading and partitioning computational meshes generated with [`Gmsh`](https://gmsh.info/).
- Defining discrete interpolation spaces on computational meshes including continuous and discontinuous interpolations.
- Integrating functionals, linear-, and bilinear-forms on triangulated manifolds for a variety of problems, including scalar- and vector-valued equations, and single- and multi-field systems.
- Discretizing PDEs into systems of linear-, nonlinear-, and differential-algebraic equations.
- Representing such algebraic problems in a way that can be readily solved using external tools such as [`PartitionedSolvers.jl`](https://github.com/PartitionedArrays/PartitionedArrays.jl), [`PetscCall.jl`](https://github.com/PartitionedArrays/PetscCall.jl), [`LinearSolve.jl`](https://github.com/SciML/LinearSolve.jl), [`NonLinearSolve.jl`](https://github.com/SciML/NonlinearSolve.jl), and [`DifferentialEquations.jl`](https://github.com/SciML/DifferentialEquations.jl).
- Visualizing results with [`Makie.jl`](https://github.com/MakieOrg/Makie.jl) and [`WriteVTK.jl`](https://github.com/JuliaVTK/WriteVTK.jl).










## Pre-requisites

You need to be fluent in Julia before using GalerkinToolkit. You can learn Julia using the learning materials in [julialang.org](https://julialang.org/) or the lecture notes in [https://www.francescverdugo.com/XM_40017/dev/](https://www.francescverdugo.com/XM_40017/dev/).

It is also required to be familiar with the key concepts in the FEM. The basics are explained in the [`Tutorials`](@ref) section. For more in depth introduction, you can use the following books:
  - C. Johnson [Johnson2009](@cite)
  - J. Whiteley [Whiteley2017](@cite)
  - S.C. Brenner and L. R. Scott [Brenner2007](@cite)


## Installation

You first need to have Julia installed in your system. See the official Julia installation instructions [here](https://julialang.org/install/).


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








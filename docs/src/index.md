```@meta
CurrentModule = GalerkinToolkit
```

# Home

Welcome to the documentation page of GalerkinToolkit!

!!! warning
    GalerkinToolkit is already useful in many situations, but keep it mind that it is under active development.
!!! warning
    This documentation page is under construction. Expect incomplete parts. You can contribute to fill them in. See the [Contributing](@ref) section for details.

## Contents

- **[Home](@ref)** -- This page.
- **[Examples](@ref)** -- They provide a quick overview of the main functionality of the library.
- **[Tutorials](@ref)** -- Learning material to get familiar with the building blocks of the library and the basics of the FEM. They are useful even if you are FEM expert as they walk you through the library step by step.
- **[Manual](@ref)** -- The user and developer guide. It gives the detailed explanations on how to use and extend the library.
- **[API](@ref)** -- All the docstrings are here.

!!! note
    This documentation is organized based on the [Diátaxis model](https://diataxis.fr/).

## What

GalerkinToolkit provides a collection of tools to build computer codes to solve partial differential equations (PDEs).  It provides the building blocks of finite element methods (FEMs) and supports different computing platforms, from laptops to modern GPU-based supercomputers.
GalerkingToolkit is fully implemented in the [Julia programming language](https://julialang.org/) and provides (or it will soon provide) tools for:

- Reading and partitioning computational meshes generated from external mesh generators ([`Gmsh`](https://gmsh.info/) at this moment).
- Defining discrete interpolation spaces on computational meshes including continuous and discontinuous interpolations.
- Integrating 0, 1, and 2-forms on triangulated manifolds for a variety of problems, including scalar- and vector-valued equations, and single- and multi-field systems.
- Discretizing PDEs into systems of linear-, nonlinear-, and differential-algebraic equations, including matrix-free operators.
- Representing such algebraic problems in a way that can be readily solved using external tools such as [`PartitionedSolvers.jl`](https://github.com/PartitionedArrays/PartitionedArrays.jl), [`PetscCall.jl`](https://github.com/PartitionedArrays/PetscCall.jl), [`LinearSolve.jl`](https://github.com/SciML/LinearSolve.jl), [`NonLinearSolve.jl`](https://github.com/SciML/NonlinearSolve.jl), and [`DifferentialEquations.jl`](https://github.com/SciML/DifferentialEquations.jl).
- Automatic differentiation for non-linear and parametric problems.
- Visualizing results with [`Makie.jl`](https://github.com/MakieOrg/Makie.jl) and [`WriteVTK.jl`](https://github.com/JuliaVTK/WriteVTK.jl).

From these items, the ones currently missing are:

- Automatic differentiation (expected by 2026).
- Matrix-free 2-forms (expected by 2026).
- Distributed computing (only partially available at the moment, expected by 2026).
- Single- and multi-GPU support (ongoing PhD project, expected by 2027).

## Why

GalerkinToolkit is definitively not the first FEM software project out there, but it has some unique design features:

- It combines the vision of frameworks like [FEniCS](https://fenicsproject.org/) and libraries like [Deal-ii](https://www.dealii.org/) in a single package. It provides a high-level API backed with automatic code generation, and also easy access to the low-level building blocks of FE codes. See in the [Examples](@ref) section examples using the high-level API and automatic-code generation and also examples implementing the integration loops "by hand" using the low-level building blocks.
- It is designed to blend with the Julia package ecosystem, and reuses as much as possible from it. For instance, GalerkinToolkit does not define differential operators, it simply uses the operators already defined in [`ForwardDiff.jl`](https://github.com/JuliaDiff/ForwardDiff.jl). You can use [`StaticArrays.jl`](https://github.com/JuliaArrays/StaticArrays.jl) and [`Tensors.jl`](https://github.com/Ferrite-FEM/Tensors.jl), when working with vector- and tensor- valued PDEs. It also makes straight-forward to use external solvers to solve the algebraic systems you get after discretizing a PDE. For visualization, it defines recipes for [`Makie.jl`](https://github.com/MakieOrg/Makie.jl) and defines helper functions to save results with [`WriteVTK.jl`](https://github.com/JuliaVTK/WriteVTK.jl) in [`vtk` format](https://vtk.org/).
- It is based on a new form compiler, the GalerkinToolkit form compiler (GTFC) that leverages the meta-programming features of the Julia programming language. The GTFC is similar in spirit to other form compilers like FFC and TSFC, but it has its own features. It is not based on the UFL domain specific language, but in a subset of the Julia programming language. This allows one to easily include user-defined types and operations in the PDE weak form, and simplifies the compiler implementation as it reuses the Julia compiler and the LLVM when possible. In addition, GTFC allows the definition of highly sophisticated weak forms combining arguments and integrals defined on domains of different number of space dimensions. For instance, you can have multi-field weak forms coupling arguments defined in a volume and on a surface, while allowing DG-like penalty terms for the argument defined on the volume. The GTFC also allows one to parametrize the weak form and its generated code with respect to discrete fields and other quantities for non-linear and time-dependent problems.
- GalerkinToolkit started as a complete re-write of the core of the [Gridap](https://github.com/gridap/Gridap.jl) project. If you liked Gridap, you will provably like GalerkinToolkit. The fundamental difference between GalerkinToolkit and Gridap is that GalerkinToolkit is not based on "lazily mapping arrays" like Gridap, but on a compiler of variational forms and accessor functions to get quantities at the integration points. The former is used to automatically generate integration loops and the later are used to implement integration loops "by-hand".

## How to start

Look into the [Examples](@ref) section for an overview of how to solve different types of PDEs with GalerkinToolkit. Learn how to use the library (and the basics of FEM) with the tutorials in the [Tutorials](@ref) section. In the [Manual](@ref) section you will find explanations about the software design and guides for users and developers.

You need to be fluent in Julia before using GalerkinToolkit. You can learn Julia using the learning materials in [julialang.org](https://julialang.org/) or the lecture notes in [https://www.francescverdugo.com/XM_40017/dev/](https://www.francescverdugo.com/XM_40017/dev/).


## Help and discussion

- You can open a new discussion to ask questions [here](https://github.com/GalerkinToolkit/GalerkinToolkit.jl/discussions).
- If you have found a bug, open an issue [here](https://github.com/GalerkinToolkit/GalerkinToolkit.jl/issues). Do not forget to include a (minimal) reproducer.

## How to cite

See the [`CITATION.cff`](https://github.com/GalerkinToolkit/GalerkinToolkit.jl/blob/main/CITATION.cff) file.

## Contributing

This package is under active development and there are several ways to contribute:

- by enhancing the documentation (e.g., fixing typos, enhancing doc strings, adding examples).
- by addressing one of the [issues waiting for help](https://github.com/GalerkinToolkit/GalerkinToolkit.jl/labels/help%20wanted).
- by adding more tests to increase the code coverage.
- by extending the current functionality. In this case, open a discussion [here](https://github.com/GalerkinToolkit/GalerkinToolkit.jl/discussions) to coordinate with the package maintainers before proposing significant changes.

Discuss with the package authors before working on any non-trivial contribution.

## Acknowledgments

Since July 2024, this package is being developed with support from the [Netherlands eScience Center](https://www.esciencecenter.nl/) under grant ID [NLESC.SS.2023.008](https://research-software-directory.org/projects/hp2sim).



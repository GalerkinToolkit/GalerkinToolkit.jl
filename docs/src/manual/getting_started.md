# Getting started

## Background knowledge

You need to be fluent in Julia before using GalerkinToolkit. You can learn Julia using the learning materials in [julialang.org](https://julialang.org/) or the lecture notes in [https://www.francescverdugo.com/XM_40017/dev/](https://www.francescverdugo.com/XM_40017/dev/).

It is also required to be familiar with the key concepts in the FEM. The basics are explained in the [`Tutorials`](@ref) section. For more in depth introduction, you can use the following books:
  - [Finite Element Methods: A Practical Guide](https://link.springer.com/book/10.1007/978-3-319-49971-0) by J. Whiteley.
  - Numerical Solution of Partial Differential Equations by the Finite Element Method by C. Johnson.
  - The Mathematical Theory of Finite Element Methods by S.C. Brenner and L. R. Scott.


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

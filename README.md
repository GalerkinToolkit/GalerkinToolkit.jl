# GalerkinToolkit

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://GalerkinToolkit.github.io/GalerkinToolkit.jl/stable/)
[![Dev Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://GalerkinToolkit.github.io/GalerkinToolkit.jl/dev/)
[![Build Status](https://github.com/GalerkinToolkit/GalerkinToolkit.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/GalerkinToolkit/GalerkinToolkit.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Test workflow status](https://github.com/GalerkinToolkit/GalerkinToolkit.jl/actions/workflows/Test.yml/badge.svg?branch=main)](https://github.com/GalerkinToolkit/GalerkinToolkit.jl/actions/workflows/Test.yml?query=branch%3Amain)
[![Lint workflow Status](https://github.com/GalerkinToolkit/GalerkinToolkit.jl/actions/workflows/Lint.yml/badge.svg?branch=main)](https://github.com/GalerkinToolkit/GalerkinToolkit.jl/actions/workflows/Lint.yml?query=branch%3Amain)
[![Docs workflow Status](https://github.com/GalerkinToolkit/GalerkinToolkit.jl/actions/workflows/Docs.yml/badge.svg?branch=main)](https://github.com/GalerkinToolkit/GalerkinToolkit.jl/actions/workflows/Docs.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/GalerkinToolkit/GalerkinToolkit.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/GalerkinToolkit/GalerkinToolkit.jl)
[![DOI](https://zenodo.org/badge/497260571.svg)](https://doi.org/10.5281/zenodo.13938389)

[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](CODE_OF_CONDUCT.md)
<!-- [![All Contributors](https://img.shields.io/github/all-contributors/GalerkinToolkit/GalerkinToolkit.jl?labelColor=5e1ec7&color=c0ffee&style=flat-square)](#contributors) -->
[![BestieTemplate](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/JuliaBesties/BestieTemplate.jl/main/docs/src/assets/badge.json)](https://github.com/JuliaBesties/BestieTemplate.jl)

## What

This package aims at providing a fully-fledged finite-element toolbox in pure Julia, with support for different computing systems from laptops to supercomputers and GPUs.

NB. This package is work in progress; a proof-of-concept API is already available (for CPUs). The package is not production ready at this point. Planned performance and documentation improvements are needed.

## Why

This package follows a new approach to implement finite-element methods based on the lessons learned in the [Gridap](https://github.com/gridap/Gridap.jl) project.

## Acknowledgments

Since July 2024, this package is being developed with support from the [Netherlands eScience Center](https://www.esciencecenter.nl/) under grant ID [NLESC.SS.2023.008](https://research-software-directory.org/projects/hp2sim).

## Documentation

- [**STABLE**](https://GalerkinToolkit.github.io/GalerkinToolkit.jl/stable) &mdash; **Documentation for the most recently tagged version.**
- [**LATEST**](https://GalerkinToolkit.github.io/GalerkinToolkit.jl/dev) &mdash; *Documentation for the in-development version.*

## How to Cite

If you use GalerkinToolkit.jl in your work, please cite using the reference given in [CITATION.cff](https://github.com/GalerkinToolkit/GalerkinToolkit.jl/blob/main/CITATION.cff).

## Help and discussion

- You can open a new discussion to ask questions [here](https://github.com/GalerkinToolkit/GalerkinToolkit.jl/discussions).
- If you have found a bug, open an issue [here](https://github.com/GalerkinToolkit/GalerkinToolkit.jl/issues). Do not forget to include a (minimal) reproducer.

## Contributing

If you want to make contributions of any kind, please first that a look into our [contributing guide directly on GitHub](docs/src/90-contributing.md) or the [contributing page on the website](https://GalerkinToolkit.github.io/GalerkinToolkit.jl/dev/90-contributing/)

## Examples

### "Hello world" Laplace PDE

The following code solves a Laplace PDE with Dirichlet boundary conditions.

```julia
import GalerkinToolkit as GT
import ForwardDiff
using LinearAlgebra
mesh = GT.mesh_from_gmsh("assets/demo.msh")
GT.label_boundary_faces!(mesh;physical_name="boundary")
Ω = GT.interior(mesh)
Γd = GT.boundary(mesh;physical_names=["boundary"])
k = 2
V = GT.lagrange_space(Ω,k;dirichlet_boundary=Γd)
uhd = GT.dirichlet_field(Float64,V)
u = GT.analytical_field(sum,Ω)
GT.interpolate_dirichlet!(u,uhd)
dΩ = GT.measure(Ω,2*k)
gradient(u) = x->ForwardDiff.gradient(u,x)
∇(u,x) = GT.call(gradient,u)(x)
a(u,v) = GT.∫( x->∇(u,x)⋅∇(v,x), dΩ)
l(v) = 0
x,A,b = GT.linear_problem(uhd,a,l)
x .= A\b
uh = GT.solution_field(uhd,x)
GT.vtk_plot("results",Ω) do plt
    GT.plot!(plt,u;label="u")
    GT.plot!(plt,uh;label="uh")
end
```

### Multi-field example

This code solves the same boundary value problem, but using an auxiliary field of Lagrange
multipliers to impose Dirichlet boundary conditions.

```julia
import GalerkinToolkit as GT
import ForwardDiff
using LinearAlgebra
mesh = GT.mesh_from_gmsh("assets/demo.msh")
GT.label_boundary_faces!(mesh;physical_name="boundary")
Ω = GT.interior(mesh)
Γd = GT.boundary(mesh;physical_names=["boundary"])
k = 2
V = GT.lagrange_space(Ω,k)
Q = GT.lagrange_space(Γd,k-1; conformity=:L2)
VxQ = V × Q
dΓd = GT.measure(Γd,2*k)
gradient(u) = x->ForwardDiff.gradient(u,x)
∇(u,x) = GT.call(gradient,u)(x)
a((u,p),(v,q)) =
    GT.∫( x->∇(u,x)⋅∇(v,x), dΩ) +
    GT.∫(x->
        (u(x)+p(x))*(v(x)+q(x))
        -u(x)*v(x)-p(x)*q(x), dΓd)
l((v,q)) = GT.∫(x->u(x)*q(x), dΓd)
x,A,b = GT.linear_problem(Float64,VxQ,a,l)
x .= A\b
uh,qh = GT.solution_field(VxQ,x)
GT.vtk_plot("results",Ω) do plt
    GT.plot!(plt,u;label="u")
    GT.plot!(plt,uh;label="uh")
end
```

# GalerkinToolkit

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://fverdugo.github.io/GalerkinToolkit.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://fverdugo.github.io/GalerkinToolkit.jl/dev/)
[![Build Status](https://github.com/fverdugo/GalerkinToolkit.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/fverdugo/GalerkinToolkit.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/fverdugo/GalerkinToolkit.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/fverdugo/GalerkinToolkit.jl)

## What

This package aims at providing a fully-fledged finite-element toolbox in pure Julia, with support for different computing systems from laptops to supercomputers and GPUs.  

NB. This package is work in progress; a proof-of-concept API is already available (for CPUs). The package is not production ready at this point. Planned performance and documentation improvements are needed.

## Why

This package follows a new approach to implement finite-element methods based on the lessons learned in the [Gridap](https://github.com/gridap/Gridap.jl) project.

## Acknowledgments

Since July 2024, this package is being developed with support from the [Netherlands eScience Center](https://www.esciencecenter.nl/) under grant ID [NLESC.SS.2023.008](https://research-software-directory.org/projects/hp2sim).

## Documentation

- [**STABLE**](https://fverdugo.github.io/GalerkinToolkit.jl/stable) &mdash; **Documentation for the most recently tagged version.**
- [**LATEST**](https://fverdugo.github.io/GalerkinToolkit.jl/dev) &mdash; *Documentation for the in-development version.*

## Help and discussion

- You can open a new discussion to ask questions [here](https://github.com/fverdugo/GalerkinToolkit.jl/discussions).
- If you have found a bug, open an issue [here](https://github.com/fverdugo/GalerkinToolkit.jl/issues). Do not forget to include a (minimal) reproducer.

## Contributing

This package is under active development and there are several ways to contribute:

- by enhancing the documentation (e.g., fixing typos, enhancing doc strings, adding examples).
- by addressing one of the [issues waiting for help](https://github.com/fverdugo/GalerkinToolkit.jl/labels/help%20wanted).
- by adding more tests to increase the code coverage.
- by extending the current functionality. In this case, open a discussion [here](https://github.com/fverdugo/GalerkinToolkit.jl/discussions) to coordinate with the package maintainers before proposing significant changes.

Discuss with the package authors before working on any non-trivial contribution.

## Examples

### "Hello world" Laplace PDE

The following code solves a Laplace PDE with Dirichlet boundary conditions.

```julia
import GalerkinToolkit as gk
import ForwardDiff
using LinearAlgebra
mesh = gk.mesh_from_gmsh("assets/demo.msh")
gk.label_boundary_faces!(mesh;physical_name="boundary")
Ω = gk.interior(mesh)
Γd = gk.boundary(mesh;physical_names=["boundary"])
k = 2
V = gk.lagrange_space(Ω,k;dirichlet_boundary=Γd)
uh = gk.zero_field(Float64,V)
u = gk.analytical_field(sum,Ω)
gk.interpolate_dirichlet!(u,uh)
dΩ = gk.measure(Ω,2*k)
gradient(u) = x->ForwardDiff.gradient(u,x)
∇(u,x) = gk.call(gradient,u)(x)
a(u,v) = gk.∫( x->∇(u,x)⋅∇(v,x), dΩ)
l(v) = 0
x,A,b = gk.linear_problem(uh,a,l)
x .= A\b
gk.vtk_plot("results",Ω) do plt
    gk.plot!(plt,u;label="u")
    gk.plot!(plt,uh;label="uh")
end
```

### Multi-field example

This code solves the same boundary value problem, but using an auxiliary field of Lagrange
multipliers to impose Dirichlet boundary conditions.

```julia
import GalerkinToolkit as gk
import ForwardDiff
using LinearAlgebra
mesh = gk.mesh_from_gmsh("assets/demo.msh")
gk.label_boundary_faces!(mesh;physical_name="boundary")
Ω = gk.interior(mesh)
Γd = gk.boundary(mesh;physical_names=["boundary"])
k = 2
V = gk.lagrange_space(Ω,k)
Q = gk.lagrange_space(Γd,k-1; conformity=:L2)
VxQ = V × Q
dΓd = gk.measure(Γd,2*k)
a((u,p),(v,q)) =
    gk.∫( x->∇(u,x)⋅∇(v,x), dΩ) +
    gk.∫(x->
        (u(x)+p(x))*(v(x)+q(x))
        -u(x)*v(x)-p(x)*q(x), dΓd)
l((v,q)) = gk.∫(x->u(x)*q(x), dΓd)
uh_qh = gk.zero_field(Float64,VxQ)
x,A,b = gk.linear_problem(uh_qh,a,l)
x .= A\b
uh,qh = uh_qh
gk.vtk_plot("results",Ω) do plt
    gk.plot!(plt,u;label="u")
    gk.plot!(plt,uh;label="uh")
end
```


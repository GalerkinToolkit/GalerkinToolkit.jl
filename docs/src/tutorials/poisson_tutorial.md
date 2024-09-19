```@meta
EditURL = "../../poisson_tutorial.jl"
```

In this tutorial, we will learn

   -  How to solve a simple PDE in Julia with Gridap
   -  How to load a discrete model (aka a FE mesh) from a file
   -  How to build a conforming Lagrangian FE space
   -  How to define the different terms in a weak form
   -  How to impose Dirichlet and Neumann boundary conditions
   -  How to visualize results


## Problem statement

In this first tutorial, we provide an overview of a complete simulation pipeline in Gridap: from the construction of the FE mesh to the visualization of the computed results. To this end, we consider a simple model problem: the Poisson equation.
 We want to solve the Poisson equation on the 3D domain depicted in next figure with Dirichlet and Neumann boundary conditions. Dirichlet boundary conditions are applied on $\Gamma_{\rm D}$, being the outer sides of the prism (marked in red). Non-homogeneous Neumann conditions are applied to the internal boundaries $\Gamma_{\rm G}$, $\Gamma_{\rm Y}$, and $\Gamma_{\rm B}$ (marked in green, yellow and blue respectively). And homogeneous Neumann boundary conditions are applied in $\Gamma_{\rm W}$, the remaining portion of the boundary (marked in white).

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*


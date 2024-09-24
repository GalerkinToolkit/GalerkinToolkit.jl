# # Poisson equation
#
# In this tutorial, we will learn
#
#    -  How to solve a simple PDE in Julia with GalerkinToolkit
#    -  How to load a discrete model (aka a FE mesh) from a file
#    -  How to build a conforming Lagrangian FE space
#    -  How to define the different terms in a weak form
#    -  How to impose Dirichlet and Neumann boundary conditions
#    -  How to visualize results
#
#
# ## Problem statement
#
# In this first tutorial, we provide an overview of a complete simulation pipeline in Gridap: from the construction of the FE mesh to the visualization of the computed results. To this end, we consider a simple model problem: the Poisson equation.
#  We want to solve the Poisson equation on the 3D domain depicted in next figure with Dirichlet and Neumann boundary conditions. Dirichlet boundary conditions are applied on $\Gamma_{\rm D}$, being the outer sides of the prism (marked in red). Non-homogeneous Neumann conditions are applied to the internal boundaries $\Gamma_{\rm G}$, $\Gamma_{\rm Y}$, and $\Gamma_{\rm B}$ (marked in green, yellow and blue respectively). And homogeneous Neumann boundary conditions are applied in $\Gamma_{\rm W}$, the remaining portion of the boundary (marked in white).
#
# ![](../assets/poisson/model-r1-2.png)
#
#  Formally, the problem to solve is: find the scalar field $u$ such that
#
# ```math
# \left\lbrace
# \begin{aligned}
# -\Delta u = f  \ &\text{in} \ \Omega,\\
# u = g \ &\text{on}\ \Gamma_{\rm D},\\
# \nabla u\cdot n = h \ &\text{on}\  \Gamma_{\rm N},\\
# \end{aligned}
# \right.
# ```
#  being $n$ the outwards unit normal vector to the Neumann boundary $\Gamma_{\rm N} \doteq \Gamma_{\rm G}\cup\Gamma_{\rm Y}\cup\Gamma_{\rm B}\cup\Gamma_{\rm W}$. In this example, we chose $f(x) = 1$, $g(x) = 2$, and $h(x)=3$ on $\Gamma_{\rm G}\cup\Gamma_{\rm Y}\cup\Gamma_{\rm B}$ and $h(x)=0$ on $\Gamma_{\rm W}$. The variable $x$ is the position vector $x=(x_1,x_2,x_3)$.
#
#  ## Numerical scheme
#
#  To solve this PDE, we use a conventional Galerkin finite element (FE) method with conforming Lagrangian FE spaces (see, e.g., [1] for specific details on this formulation). The weak form associated with this formulation is: find $u\in U_g$ such that $ a(u,v) = b(v) $ for all $v\in V_0$, where $U_g$ and $V_0$ are the subset of functions in $H^1(\Omega)$ that fulfill the Dirichlet boundary condition $g$ and $0$ respectively. The bilinear and linear forms for this problems are
# ```math
#   a(u,v) \doteq \int_{\Omega} \nabla v \cdot \nabla u \ {\rm d}\Omega, \quad b(v) \doteq \int_{\Omega} v\ f  \ {\rm  d}\Omega + \int_{\Gamma_{\rm N}} v\ h \ {\rm d}\Gamma_{\rm N}.
# ```
# The problem is solved numerically by approximating the spaces $U_g$ and $V_0$ by their discrete counterparts associated with a FE mesh of the computational domain $\Omega$. As we have anticipated, we consider standard conforming Lagrangian FE spaces for this purpose.
#
# The implementation of this numerical scheme in GalerkinToolkit is done in a user-friendly way thanks to the abstractions provided by the library. As it will be seen below, all the mathematical objects involved in the definition of the discrete weak problem have a correspondent representation in the code.
#
#  ## Setup
#
#  The step number 0 in order to solve the problem is to load the GalerkinToolkit library in the code. If you have configured your Julia environment properly, it is simply done with the line:

import GalerkinToolkit as GT
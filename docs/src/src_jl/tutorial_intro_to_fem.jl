# # Introduction to FEM

# In this tutorial, we will learn:
# - The gist of the finite element method (FEM).
# - How to solve a simple partial differential equation (PDE) with it.
# - How to express the key concept in code using GalerkinToolkit.
# - How to validate the code using the method manufactured solutions.

# ## Problem statement
#
# In this tutorial, we show how to solve a simple PDE with the FEM.
# To make this introduction really an introduction we consider what is considered the
# "hello, world" PDE: the Poisson equation. In particular, we want to solve the
# Poisson equation with Dirichlet boundary conditions a given
# domain $\Omega\subset\mathbb{R}^d$ with $d$ being the number of spatial
# dimensions ($d=2$ in this example).
# Precisely, our goal is to find a numerical approximation of the
# function $u:\Omega\rightarrow\mathbb{R}$ such that

# ```math
# \left\lbrace
# \begin{aligned}
# -\Delta u = f  \ &\text{in} \ \Omega,\\
# u = g \ &\text{on}\ \partial\Omega,\\
# \end{aligned}
# \right.
# ```
# where $\Delta u = \sum_{i=1}^d \partial^2 u / \partial x_i^2$ is the Laplace operator and $f$, $g$
# are two given functions $f,g:\Omega\rightarrow\mathbb{R}$. For simplicity
# we will define $\Omega$ as a disk of radius one. This is a simple geometry, but yet more complex than a 
# two-dimensional box. We this, we make want to illustrate that FEM can be used to solve PDEs on complex
# geometries beyond simple "boxes".
#
# ## The method of manufactured solutions
#
# In this example, we are going to select $f$ and $g$ in such a way $u$ is a known function.
# This will allow us to compare the numerical approximation computed with FEM against the theoretical
# exact solution $u$. This techniques is known as the method of manufactured solutions.
# Let us, "manufacture" $f$ and $g$ such that function $u(x)=(\sum_{i=1}^d x_i)^p$ is
# the solution of the PDE above. The scalar $p$ is a given integer $p>0$. It will be useful to see how
# the numerical solution will behave for different values of $p$.
#
# To manufacture function $f$ and $g$ we applying the PDE operators to the given function $u$. That is, $f$ needs to be computed as
# $f= -\Delta ((\sum_{i=1}^d x_i)^p)$ and $g$ is simply $g(x)=(\sum_{i=1}^d x_i)^p$. Applying the Laplace operator
# we get the closed-form expression for $f$, namely $f(x)= p(p-1)(\sum_{i=1}^d x_i)^{(p-2)}$.

# ## The finite element approximation
#
# The key idea of FEM is to transform a PDE into a system of linear algebraic equations of the form $Ax=b$,
# where $A$ is a matrix and $b$ is a vector. This reduces the problem of finding a function $u$ to finding
# vector $x$, which can be done on a computer using linear algebra libraries.
# To do this, FEM does not looks for the exact function $u$, but for functions that can be written as a 
# linear combination of a finite number of functions, namely
# ```math
# u^\mathrm{fem}(x)=\sum_{j=1}^N \alpha_i s_i(x),
# ```
# where $\alpha_i$ are the coefficients of the linear combination, $s_i$ are functions such that
# $s_i:\Omega\rightarrow\mathbb{R}$ and $N$ is an integer.
# The goal of FEM is to find suitable values for $\alpha_i$, $s_i(x)$ and $N$
# such that $u^\mathrm{fem}$ is a good approximation of the exact solution $u$, namely $u^\mathrm{fem}(x)\approx u(x)$ for points $x\in\Omega$.
# Spoiler alert: the more computational effort we put in building function $u^\mathrm{fem}$ the better 
# will be the approximation.
#
# ## Workflow
#
# Function $u^\mathrm{fem}$ is built with the following phases. First we define the auxiliary functions $s_i(x)$ and $N$.
# This step is often referred as the "numerical discretization". The next step is building a system of linear
# algebraic equations $Ax=b$, which is referred to as the "FEM assembly". Then one solves for the vector $x$ in
# what is called the "solver" or "solution" step.
# At this points, the coefficients $\alpha_i$ can be computed using both vector $x$ and the Dirichlet boundary conditions of the PDE.
# The next step is to find the coefficients $\alpha$. The final step is typically some post-process of
# function $u^\mathrm{fem}$. For instance, visualize it, compute some quantity of interest, etc. In summary,
# these are the key phases in a FEM computation:
# - Discretization
# - Assembly
# - Solution
# - Post-process
# 
#
# ## Discretization
#
# This phase stats by building a "triangulation" of the domain $\Omega$ in which the PDE is defined.
# A triangulation is a set of simpler domains $T_k\subset\mathbb{R}^d$, whose union is an approximation of $\Omega$,
# namely $\cup_{k=1}^M T_k\approx\Omega$. Each domain $T_k$ is often called an "element" or a "cell" or a "face",
# and they are typically simple polyhedra such as triangles, tetrahedral, hexahedra, etc. The integer $M$ denotes
# the number of elements here. The triangulation is also often called a computational mesh or a computational grid.
#
# Let's build a mesh for our domain $\Omega$ using code. First, let us load all packages that we will use in this tutorial:
#
using LinearAlgebra
import GalerkinToolkit as GT
import PartitionedSolvers as PS
import ForwardDiff
import GLMakie as Makie
import FileIO # hide

# A common practice in GalerkinToolkit is to use `using` for packages in the Julia standard library, and using `import` otherwise. 
# This makes clear from which package each function comes from, while assuming that developers already know
# functions in the standard library.
#
# The following cells build a mesh object using the external mesh generation GMSH.
# Note that the variable `mesh_size` controls how fine are the cells (smaller is finer). We start with a 
# coarse mesh to make visualization easier. In this tutorial we are not going to comment in detail all
# code lines. We will discuss only the parts relevant in this high level introduction.

mesh_size = 0.1
R = 1 #Radius
mesh = GT.with_gmsh() do
    dim = 2
    circle_tag = gmsh.model.occ.add_circle(0,0,0,R)
    circle_curve_tag = gmsh.model.occ.add_curve_loop([circle_tag])
    circle_surf_tag = gmsh.model.occ.add_plane_surface([circle_curve_tag])
    gmsh.model.occ.synchronize()
    gmsh.model.model.add_physical_group(dim,[circle_surf_tag],-1,"Omega")
    gmsh.option.setNumber("Mesh.MeshSizeMax",mesh_size)
    gmsh.model.mesh.generate(dim)
    GT.mesh_from_gmsh_module()
end

# We returned an object representing the triangulation. There are ways of accessing the low level information
# in this mesh, but we are not going to discuss them in this tutorial. Here, we only need to know how to
# visualize the mesh to know what we are doing. The mesh can be visualized both using Paraview and Makie.
# We use Makie in this tutorial.

Makie.plot(mesh;color=:pink)

# This gives a picture of the domain $\Omega$.
# If we want to see also the mesh cells, we can add a color to their edges like this:

Makie.plot(mesh;color=:pink,strokecolor=:blue)



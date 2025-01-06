module GridapRunner

# using GridapPaper
# using DrWatson

using Gridap
using Gridap.Io
using Gridap.FESpaces
using Gridap.Algebra
using Gridap.Geometry
using Gridap.ReferenceFEs
using SparseArrays
# using IterativeSolvers: cg, cg!
# using Preconditioners: AMGPreconditioner, SmoothedAggregation
# import GridapPETSc
# import GridapPETSc: PetscScalar, PetscInt, PETScSolver, PETScVector, PETScMatrix, PETSC, Mat, Vec, @check_error_code
# using MPI
# using SparseMatricesCSR


u_p1(x) = 1 + x[1] + 2*x[2]
u_p2(x) = 1 + x[1]*x[1] + 2*x[2]*x[2]
u_p3(x) = 1 + x[1]*x[1]*x[1] + 2*x[2]*x[2]*x[2]
f_p3(x) = -6*x[1] - 12*x[2]

function poisson_sol(k)
    if k == 1
      u = u_p1
      f = 0
    elseif k == 2
      u = u_p2
      f = -6
    elseif k == 3
      u = u_p3
      f = f_p3
    else
      error("k $k not valid")
    end
    u,f
  
end

function poisson_cg_assembly(n, k, iscube)
  outputs = Dict()

  outputs["assembly [s]"] = @elapsed begin
    u,f = poisson_sol(k)
    if iscube
      pmin = Point(0,0,0); pmax = Point(1,1,1); partition = (n,n,n)
      model_hex = CartesianDiscreteModel(pmin,pmax,partition)
      model = simplexify(model_hex)
    else
      bsonfile = datadir("experiments","piece_n$n.bson")
      model = from_bson_file(DiscreteModel,bsonfile)
    end
    reffe = ReferenceFE(lagrangian,Float64,k)
    V = TestFESpace(model,reffe,dirichlet_tags="boundary")
    U = TrialFESpace(V,u)
    Ω = Triangulation(model)
    degree = max((k-1)+(k-1),k+max(0,k-2))
    dΩ = Measure(Ω,Quadrature(strang,degree))
    # TODO Implement something like this
    # op = AffineFEOperator(a,l,assem,U,V)
    # AffineFEOperator!(op,a,l,assem,U,V)
    du = get_trial_fe_basis(U)
    dv = get_fe_basis(V)
    a = ∫( ∇(dv)⋅∇(du) )dΩ
    L = ∫( dv*f )dΩ
    uhd = zero(U)
    data = collect_cell_matrix_and_vector(U,V,a,L,uhd)
    Tm = SparseMatrixCSC{Float64,Int32}
    Tv = Vector{Float64}
    assem = SparseMatrixAssembler(Tm,Tv,U,V)
    A, b = assemble_matrix_and_vector(assem,data)
  end

  outputs["assembly! [s]"] = @elapsed assemble_matrix_and_vector!(A,b,assem,data)

  # Residual
  uh = interpolate(u,U)
  x = get_free_dof_values(uh)
  r = b - A*x
  rl2 = norm(r)
  bl2 = norm(b)

  outputs["rl2"] = rl2
  outputs["bl2"] = bl2
  outputs["rl2/bl2"] = rl2/bl2
  outputs["num_cells"] = num_cells(model)
  outputs["num_dofs"] = length(b)
  outputs["num_nz"] = nnz(A)
  outputs["n"] = n
  outputs["k"] = k
  outputs["experiment"] = "poisson_cg_assembly"
  outputs["lib"] = "gridap"
  outputs["iscube"] = iscube

  return outputs
end


end # module

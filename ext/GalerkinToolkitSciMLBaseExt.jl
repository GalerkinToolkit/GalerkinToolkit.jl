module GalerkinToolkitSciMLBaseExt

import GalerkinToolkit as GT
import SciMLBase
import PartitionedSolvers as PS
using LinearAlgebra

function GT.SciMLBase_LinearProblem(args...)
    p = GT.PartitionedSolvers_linear_problem(args...)
    A = PS.matrix(p)
    b = PS.rhs(p)
    x = PS.solution(p)
    SciMLBase.LinearProblem(A,b;u0=x)
end

function GT.solution_field(U::GT.AbstractSpace,sol::SciMLBase.LinearSolution)
    x = sol.u
    GT.solution_field(U,x)
end

function GT.solution_field(U::GT.DiscreteField,sol::SciMLBase.LinearSolution)
    x = sol.u
    GT.solution_field(U,x)
end

function GT.SciMLBase_NonlinearProblem(uh::GT.DiscreteField,r,j;
        matrix_strategy = GT.monolithic_matrix_assembly_strategy(),
        vector_strategy = GT.monolithic_vector_assembly_strategy(),
    )

    U = GT.space(uh)
    V = U
    x = vector_strategy.values_to_solution(GT.free_values(uh))
    T = eltype(x)
    parameters = (uh,)
    b,residual_cache = GT.assemble_vector(r(uh),T,V;parameters,reuse=Val(true),vector_strategy)
    A,jacobian_cache = GT.assemble_matrix(j(uh),T,U,V;parameters,reuse=Val(true),matrix_strategy)

    function f(dx,x,p)
        GT.solution_field!(uh,x)
        GT.update_vector!(dx,residual_cache;parameters)
        dx
    end

    function jac(J,x,p)
        GT.solution_field!(uh,x)
        GT.update_matrix!(J,jacobian_cache;parameters)
        J
    end

    jac_prototype = A
    nlfun = SciMLBase.NonlinearFunction(f;jac,jac_prototype)
    SciMLBase.NonlinearProblem(nlfun,x)
end

function GT.solution_field(U::GT.AbstractSpace,sol::SciMLBase.NonlinearSolution)
    x = sol.u
    GT.solution_field(U,x)
end

function GT.solution_field(U::GT.DiscreteField,sol::SciMLBase.NonlinearSolution)
    x = sol.u
    GT.solution_field(U,x)
end

function GT.SciMLBase_ODEProblem(interval,uh::GT.DiscreteField,m,r,j;
    dirichlet_dynamics! =nothing,
    matrix_strategy = GT.monolithic_matrix_assembly_strategy(),
    vector_strategy = GT.monolithic_vector_assembly_strategy(),
    )

    t = first(interval)
    U = GT.space(uh)
    V = U
    x = vector_strategy.values_to_solution(GT.free_values(uh))
    T = eltype(x)
    parameters = map(GT.parameter,(uh,t))
    b,residual_cache = GT.assemble_vector(r(parameters...),T,V;parameters,reuse=Val(true),vector_strategy)
    A,jacobian_cache = GT.assemble_matrix(j(parameters...),T,U,V;parameters,reuse=Val(true),matrix_strategy)
    M,Md = GT.assemble_matrix_with_free_and_dirichlet_columns(m,T,U,V;matrix_strategy)
    vh = GT.zero_field(T,U)
    vxd = vector_strategy.values_to_solution(GT.dirichlet_values(vh))

    function f(dx,x,p,t)
        if dirichlet_dynamics! !== nothing
            dirichlet_dynamics!(t,uh,vh)
        end
        GT.solution_field!(uh,x)
        vector_strategy.values_to_solution!(vxd,GT.dirichlet_values(vh))
        GT.update_vector!(dx,residual_cache;parameters)
        mul!(dx,Md,vxd,-1,1)
        dx
    end

    function jac(J,x,p,t)
        if dirichlet_dynamics! !== nothing
            dirichlet_dynamics!(t,uh,nothing)
        end
        GT.solution_field!(uh,x)
        GT.update_matrix!(J,jacobian_cache;parameters)
        J
    end

    jac_prototype = A
    mass_matrix = M
    nlfun = SciMLBase.ODEFunction{true,SciMLBase.NoSpecialize}(f;mass_matrix,jac,jac_prototype)
    workspace = (;uh,vh,dirichlet_dynamics!)
    p = workspace
    SciMLBase.ODEProblem(nlfun,x,interval,p)
end


end # module

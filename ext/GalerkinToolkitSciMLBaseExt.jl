module GalerkinToolkitSciMLBaseExt

import GalerkinToolkit as GT
import SciMLBase
import PartitionedSolvers as PS

function GT.SciMLBase_LinearProblem(args...)
    p = GT.linear_problem(args...)
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

function GT.SciMLBase_NonlinearProblem(uh::GT.DiscreteField,r,j)

    U = GT.space(uh)
    V = U
    x = GT.free_values(uh)
    T = eltype(x)
    parameters = (uh,)
    b,residual_cache = GT.assemble_vector(r(uh),T,V;parameters)
    A,jacobian_cache = GT.assemble_matrix(j(uh),T,U,V;parameters)

    function f(dx,x,p)
        duh = GT.solution_field(uh,x)
        GT.update_vector!(dx,residual_cache;parameters=(duh,))
        dx
    end

    function jac(J,x,p)
        duh = GT.solution_field(uh,x)
        GT.update_matrix!(A,jacobian_cache;parameters=(duh,))
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


end # module

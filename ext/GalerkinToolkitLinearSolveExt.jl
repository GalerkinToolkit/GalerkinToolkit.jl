module GalerkinToolkitLinearSolveExt

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

end # module

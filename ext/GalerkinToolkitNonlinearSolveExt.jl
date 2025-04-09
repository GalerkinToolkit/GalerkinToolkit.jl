module GalerkinToolkitNonlinearSolveExt

import GalerkinToolkit as GT
import SciMLBase

function GT.SciMLBase_NonlinearProblem(uh::DiscreteField,r,j)

    U = GT.space(uh)
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

end # module

module GalerkinToolkitOrdinaryDiffEqCoreExt

import GalerkinToolkit as GT
import OrdinaryDiffEqCore

function GT.solution_field(integrator::OrdinaryDiffEqCore.ODEIntegrator;derivative=Val(0))
    workspace = integrator.p
    (;uh,vh,dirichlet_dynamics!) = workspace
    t = integrator.t
    d = GT.val_parameter(derivative)
    if d == 0
        x = integrator(t)
        dirichlet_dynamics!(t,uh,nothing)
        GT.solution_field!(uh,x)
        uh
    elseif d == 1
        x = integrator(t,Val{1})
        dirichlet_dynamics!(t,nothing,vh)
        GT.solution_field!(vh,x)
        vh
    else
        # This would require to get also higher order
        # derivatives of the Dirichlet function.
        error("Case not implemented.")
    end
end

end # module

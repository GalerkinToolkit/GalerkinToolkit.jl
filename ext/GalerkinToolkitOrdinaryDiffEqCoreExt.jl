module GalerkinToolkitOrdinaryDiffEqCoreExt

import GalerkinToolkit as GT
import OrdinaryDiffEqCore

function GT.solution_field(integrator::OrdinaryDiffEqCore.ODEIntegrator;derivative=Val(0))
    workspace = integrator.p
    (;uhd,vhd,dirichlet_dynamics!) = workspace
    t = integrator.t
    dirichlet_dynamics!(t,uhd,nothing)
    d = GT.val_parameter(derivative)
    if d == 0
        x = integrator(t)
        GT.solution_field(uhd,x)
    elseif d == 1
        x = integrator(t,Val{1})
        GT.solution_field(vhd,x)
    else
        # This would require to get also higher order
        # derivatives of the Dirichlet function.
        error("Case not implemented.")
    end
end

end # module

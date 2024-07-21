module test_functions

import GalerkinToolkit as gk
using GalerkinToolkitExamples: Example003, Example001


tol = 1.0e-10
params = Dict{Symbol,Any}()


#params[:mesh] = gk.cartesian_mesh((0,10,0,10),(2,2))
params[:mesh] = gk.cartesian_mesh((0,3,0,2,0,1),(20,20,20))
# (80,80,80) is about 18 seconds for jacobian with cpu_v1 500k cells | flux
# (60,60,60) is about 7 seconds for jacobian with cpu_v1 | flux

# First test original correct solution to compare the rest to
params[:export_vtu] = false
ns_to_s = 1e+9
params[:float_type] = Dict(:Float => Float32, :Int => Int32)
params[:p] = 2
params[:jacobian_implementation] = :cpu_v1
params[:autodiff] = :flux


# From here you pull out the time measurement and the number of Newton iterations for convergence and store it in a df.

# Test the kernel generic of the cpu extension

timer = Dict("Initial_setup" => [],
                 "Jacobian" => [],
                 "Jacobian_transfer" => []
                )
params[:timer_dict] = timer

println("Jacobian extension ", params[:jacobian_implementation])
t = @elapsed results  = Example003.main(params)
#t = @elapsed results, x  = Example003.main(params)

#println(results[:timer_dict])
#println(t)
#println("iterations ",results[:iterations])
# println("Number of cells: ", results2[:ncells])

println(t)
end # module

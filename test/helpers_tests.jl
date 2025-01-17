module HelpersTests

import GalerkinToolkit as GT
#using InteractiveUtils

options = GT.options()
@show GT.real_type(options)
@show GT.int_type(options)

#@code_warntype GT.real_type(options)

end # module

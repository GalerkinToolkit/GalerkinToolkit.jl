module PassesTests

import GalerkinToolkit as GT
using Test
# using LinearAlgebra
# import ForwardDiff
# import PartitionedSolvers as PS

NORMAL_PASS_COUNTS = Dict()
function count_normal_pass_call(expr, id)
    NORMAL_PASS_COUNTS[id] = get(NORMAL_PASS_COUNTS, id, 0) + 1
    return expr
end
# ast = quote
#     local s2 = 0
#     for i in 1:10
#         for j in 1:20
#             for k in 1:30
#                 a2 = i * j
#                 b2 = i * 2
#                 c2 = j * 3
#                 d2 = i * k
#                 s2 += (a2 + b2) * (c2 + d2)
#             end
#         end
#     end
#     display(s2)
# end


original_expr = quote
        local s = 0
        for i in 1:10
            for j in 1:20
                for k in 1:30
                    a3 = count_normal_pass_call(i * j, 1)
                    b3 = count_normal_pass_call(i * 2, 2)
                    c3 = count_normal_pass_call(j * 3, 3)
                    d3 = count_normal_pass_call(i * k, 4)
                    s += (a3 + b3) * (c3 + d3)
                end
            end
        end
        s
    end

optimized_expr = GT.normal_licm(original_expr)
s1 = eval(original_expr)
origin_counts = NORMAL_PASS_COUNTS
NORMAL_PASS_COUNTS = Dict()
s2 = eval(optimized_expr)
optimized_counts = NORMAL_PASS_COUNTS
@test s1 == s2
@test optimized_counts == Dict(1 => 200, 2 => 10, 3=>200, 4 =>6000)



# 61041750
# eval(ast)
# ast2 = GT.normal_licm(ast)
# eval(ast)
# eval(ast2)
end # module
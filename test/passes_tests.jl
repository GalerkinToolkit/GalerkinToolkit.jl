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

# TODO: implement test cases for array_cse
# alloc_zeros = GT.alloc_zeros

# expr = quote
#     N = 20
#     Ni = 10
#     v2 = alloc_zeros("v2", Float64, N)
#     v = 0.0
#     for i in 1:Ni
#         v1 = f(i)
#         for k in 1:N
#             v2[k] = g(v1, s[k])
#         end
#         for j in 1:N
#             v3 = g(v1, s[j])
#             for k in 1:N
#                 v += v3 * v2[k]
#             end
#         end
#     end
# end

# expr = quote
#     N = 20
#     Ni = 10
#     v4 = alloc_zeros("v4", Float64, N)
#     v = 0.0
#     for i in 1:Ni
#         v1 = f1(i)
#         for k in 1:N
#             v2 = v1(s[k])
#             v3 = f3(v2)
#             v4[k] = f4(v3, v2)
#         end
#         for j in 1:N
#             v5 = v1(s[j])
#             v6 = f3(v5)
#             for k in 1:N
#                 v += v6 * v4[k]
#             end
#         end
#     end
# end

# expr = quote
#     N = 30
#     n1 = 10
#     n2 = 20
#     same_range_cache = alloc_zeros("same_range_cache", Float64, n1)
#     v = 0.0
#     for i in 1:n1
#         same_range_cache[i] = g(1.0, s[i])
#     end
#     for i in 1:n2
#         v += h(i)
#     end
#     for i in 1:n1
#         v += g(1.0, s[i])
#     end
# end


# expr = quote
#     N = 30
#     n1 = 10
#     n2 = 20
#     different_range_cache = alloc_zeros("different_range_cache", Float64, n2)
#     v = 0.0
#     for i in 1:n1
#         v += h(i)
#     end
#     for i in 1:n2
#         different_range_cache[i] = g(1.0, s[i])
#     end
#     for i in 1:n1
#         v += g(1.0, s[i])
#     end
# end
# expr = MacroTools.striplines(expr)
# new_expr, _ = GT.ast_array_cse(expr)
# new_expr = (new_expr) |> GT.ast_topological_sort |> GT.ast_remove_dead_code
# new_expr


end # module
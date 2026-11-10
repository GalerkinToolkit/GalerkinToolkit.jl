module PassesTests
# TODO: We can also improve this in a new work 
# we can place the array cse into licm: 
# weak LICM keeping arrays
# -> array cse supporting multi-dim arrays (hashing with placeholders for loop var)
# -> partial scalar replacement using the LCA to find out the minimum dimension needed for each array
# for the placeholders we need to sort them in the DFS order in the AST. It is like creating a function creating new dummy variables.
# 

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

optimized_expr = GT.ast_normal_licm(original_expr)
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

function array_cse_contains(node, expected)
    isequal(node,expected) || (node isa Expr && any(child -> array_cse_contains(child,expected),node.args))
end

function check_array_cse(expr)
    original = Base.remove_linenums!(deepcopy(expr))
    optimized,_ = GT.ast_array_cse(original)
    evaluate = code -> Base.invokelatest(Core.eval(GT,:(() -> $code)))
    expected = evaluate(original)
    @test evaluate(optimized) == expected
    cleaned = GT.ast_remove_dead_code(GT.ast_topological_sort(optimized))
    @test evaluate(cleaned) == expected
    optimized
end

@testset "Array CSE dependency lifetimes" begin
    # A completed fill may be reused while its fixed quadrature index is in scope.
    valid = check_array_cse(quote
        cache = zeros(Float64,27)
        result = zeros(Float64,27)
        for q in 1:2
            for i in 1:27
                cache[i] = atan(i,q)
            end
            for j in 1:27
                value = atan(j,q)
                result[j] += value
            end
        end
        result
    end)
    @test array_cse_contains(valid,:(value = cache[j]))

    # Reusing the name q in a later loop must not revive the previous table.
    expired = check_array_cse(quote
        cache = zeros(Float64,27)
        result = zeros(Float64,27)
        for q in 1:2
            for i in 1:27
                cache[i] = atan(i,q)
            end
        end
        for q in 1:2
            for j in 1:27
                value = atan(j,q)
                result[j] += value
            end
        end
        result
    end)
    @test array_cse_contains(expired,:(value = atan(j,q)))

    # The same lifetime rule applies to scalar caches introduced by this pass.
    check_array_cse(quote
        result = zeros(Float64,27)
        for q in 1:2
            for i in 1:27
                first_value = atan(i,q)
                result[i] += first_value
            end
        end
        for q in 1:2
            for j in 1:27
                second_value = atan(j,q)
                result[j] += second_value
            end
        end
        result
    end)

    # Tables indexed by every dependency remain valid after their filling loops.
    complete = check_array_cse(quote
        cache = zeros(Float64,27)
        result = zeros(Float64,27)
        for i in 1:27
            cache[i] = sin(i)
        end
        for j in 1:27
            value = sin(j)
            result[j] = value
        end
        result
    end)
    @test array_cse_contains(complete,:(value = cache[j]))
end

@testset "Array CSE dimension order" begin
    for (ni,nj) in ((27,27),(13,17))
        optimized = check_array_cse(quote
            cache = zeros(Float64,$nj,$ni)
            result = zeros(Float64,$ni,$nj)
            for i in 1:$ni
                for j in 1:$nj
                    cache[j,i] = atan(i,j)
                end
            end
            for p in 1:$ni
                for q in 1:$nj
                    value = atan(p,q)
                    result[p,q] = value
                end
            end
            result
        end)
        @test array_cse_contains(optimized,:(value = cache[q,p]))
    end
end

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


# =============================================================================
# End of existing tests. AI-generated initialization code motion tests follow.
# =============================================================================

module InitializationCodeMotionTests

using Test
import GalerkinToolkit as GT

@testset "Assembly preparation lifetimes" begin
    fixed_calls = Ref(0)
    refresh_calls = Ref(0)
    prepare = values -> (fixed_calls[] += 1; 2 .* values)
    refresh = (parameter,values) -> (refresh_calls[] += 1; values .* parameter[])
    expr = :((values,) -> ((parameter,) -> (alloc -> begin
        make_fixed = $prepare
        fixed = make_fixed(values)
        dynamic = ($refresh)(parameter,first((fixed,)))
        T = Float64
        buffer = alloc_zeros("buffer",T,2)
        scratch = alloc_zeros("scratch",eltype(dynamic),2)
        fill!(buffer,zero(T))
        for i in 1:2
            scratch[i] = dynamic[i]
            buffer[i] += scratch[i]
        end
        copyto!(alloc,buffer)
        nothing
    end)))
    staged = GT.ast_initialization_code_motion(expr)
    params = GT.evaluate(staged,([1.0,2.0],))
    @test fixed_calls[] == 1
    @test refresh_calls[] == 0
    parameter = Ref(3.0)
    kernel = Base.invokelatest(params,parameter)
    @test refresh_calls[] == 1
    output = zeros(2)
    Base.invokelatest(kernel,output)
    @test output == [6,12]
    output .= -1
    Base.invokelatest(kernel,output)
    @test output == [6,12]
    @test fixed_calls[] == 1
    @test refresh_calls[] == 1

    parameter[] = 4
    updated = Base.invokelatest(params,parameter)
    @test refresh_calls[] == 2
    @test getfield(updated,:buffer) === getfield(kernel,:buffer)
    Base.invokelatest(updated,output)
    @test output == [8,16]

    other_params = GT.evaluate(staged,([3.0,5.0],))
    other = Base.invokelatest(other_params,parameter)
    @test fixed_calls[] == 2
    @test getfield(other,:buffer) !== getfield(kernel,:buffer)
    Base.invokelatest(other,output)
    @test output == [24,40]
    Base.invokelatest(kernel,output)
    @test output == [6,12]
end

@testset "Explicit buffer resets and mutable reads" begin
    expr = :(() -> (() -> (output -> begin
        buffer = zeros(Float64,1)
        alias = buffer
        fill!(alias,0.0)
        before = sum(alias)
        buffer[1] += 1
        output[1] = before+buffer[1]
    end)))
    params = GT.evaluate(GT.ast_initialization_code_motion(expr),())
    kernel = Base.invokelatest(params)
    buffer = getfield(kernel,:buffer)
    output = zeros(1)
    for _ in 1:3
        Base.invokelatest(kernel,output)
        @test output == [1]
        @test getfield(kernel,:buffer) === buffer
    end

    reads = :((input,) -> (() -> (output -> begin
        alias = input
        total = sum(alias)
        for i in 1:1
            alias[i] += 1
        end
        output[1] = total
    end)))
    params = GT.evaluate(GT.ast_initialization_code_motion(reads),([1.0],))
    kernel = Base.invokelatest(params)
    Base.invokelatest(kernel,output)
    @test output == [1]
    Base.invokelatest(kernel,output)
    @test output == [2]

    calls = Ref(0)
    unknown = () -> (calls[] += 1)
    preparation_calls = Ref(0)
    prepare = input -> (preparation_calls[] += 1; input .+ 10)
    effect = :((input,) -> (() -> (output -> begin
        ($unknown)()
        fixed = ($prepare)(input)
        output[1] = first(fixed)
    end)))
    params = GT.evaluate(GT.ast_initialization_code_motion(effect),([3.0],))
    @test calls[] == 0
    @test preparation_calls[] == 1
    kernel = Base.invokelatest(params)
    Base.invokelatest(kernel,output)
    Base.invokelatest(kernel,output)
    @test calls[] == 2
    @test preparation_calls[] == 1
    @test output == [13]

    # Calls used as statements can mutate their arguments without a ! suffix.
    change = input -> (input[1] += 1; nothing)
    mutation = :((input,) -> (() -> (output -> begin
        before = sum(input)
        ($change)(input)
        after = sum(input)
        output[1] = before+after
    end)))
    params = GT.evaluate(GT.ast_initialization_code_motion(mutation),([1.0],))
    kernel = Base.invokelatest(params)
    Base.invokelatest(kernel,output)
    @test output == [3]
    Base.invokelatest(kernel,output)
    @test output == [5]

    # A compiler pass may assign a template mutator's returned value.
    assigned_effect = :((input,) -> (() -> (output -> begin
        update = setindex!
        before = sum(input)
        changed = update(input,4.0,1)
        after = sum(changed)
        output[1] = before+after
    end)))
    params = GT.evaluate(GT.ast_initialization_code_motion(assigned_effect),([1.0],))
    kernel = Base.invokelatest(params)
    Base.invokelatest(kernel,output)
    @test output == [5]
    Base.invokelatest(kernel,output)
    @test output == [8]

    copies = Ref(0)
    fresh_copy = input -> (copies[] += 1; copy(input))
    fresh = :((input,) -> (() -> (output -> begin
        values = ($fresh_copy)(input)
        values[1] += 1
        output[1] = values[1]
    end)))
    params = GT.evaluate(GT.ast_initialization_code_motion(fresh),([1.0],))
    @test copies[] == 0
    kernel = Base.invokelatest(params)
    Base.invokelatest(kernel,output)
    @test output == [2]
    Base.invokelatest(kernel,output)
    @test output == [2]
    @test copies[] == 2
end

@testset "Intermediate buffers are overwritten without implicit resets" begin
    expr = :((input,) -> ((scale,) -> (output -> begin
        buffer = alloc_zeros("scratch",Float64,3)
        for i in 1:2
            buffer[i] = input[i]*scale
        end
        copyto!(output,1,buffer,1,2)
    end)))
    params = GT.evaluate(GT.ast_initialization_code_motion(expr),([1.0,2.0],))
    kernel = Base.invokelatest(params,2.0)
    buffer = getfield(kernel,:buffer)
    fill!(buffer,NaN)
    buffer[3] = 17
    output = zeros(2)
    Base.invokelatest(kernel,output)
    @test output == [2,4]
    @test buffer[3] == 17
    buffer[1] = NaN
    buffer[2] = NaN
    updated = Base.invokelatest(params,3.0)
    @test getfield(updated,:buffer) === buffer
    Base.invokelatest(updated,output)
    @test output == [3,6]
    @test buffer[3] == 17
end

@testset "Assembly initialization guards" begin
    calls = Ref(0)
    prepare = input -> (calls[] += 1; [first(input)])
    expr = :((input,) -> (() -> (output -> begin
        if length(input) == 0
            return
        end
        sample = ($prepare)(input)
        output[1] = first(sample)
    end)))
    staged = GT.ast_initialization_code_motion(expr)
    params = GT.evaluate(staged,(Float64[],))
    kernel = Base.invokelatest(params)
    output = [7.0]
    Base.invokelatest(kernel,output)
    @test calls[] == 0
    @test output == [7]
    params = GT.evaluate(staged,([3.0],))
    @test calls[] == 1
    Base.invokelatest(Base.invokelatest(params),output)
    @test output == [3]

    dynamic_guard = :((input,) -> ((parameter,) -> (output -> begin
        if length(parameter) == 0
            return
        end
        sample = ($prepare)(input)
        output[1] = first(sample)+first(parameter)
    end)))
    params = GT.evaluate(GT.ast_initialization_code_motion(dynamic_guard),([4.0],))
    @test calls[] == 1
    kernel = Base.invokelatest(params,Float64[])
    Base.invokelatest(kernel,output)
    @test calls[] == 1
    @test output == [3]
    kernel = Base.invokelatest(params,[2.0])
    @test calls[] == 2
    Base.invokelatest(kernel,output)
    @test output == [6]
end

@testset "Custom assigned calls are pure expressions" begin
    calls = Ref(0)
    prepare = input -> (calls[] += 1; input .+ 1)
    expr = :((config,input) -> ((operation,) -> (output -> begin
        fixed = config.prepare(input)
        updated = operation(fixed)
        copyto!(output,updated)
    end)))
    params = GT.evaluate(GT.ast_initialization_code_motion(expr),((;prepare),[1.0,2.0]))
    @test calls[] == 1
    output = zeros(2)
    kernel = Base.invokelatest(params,input -> 2 .* input)
    Base.invokelatest(kernel,output)
    @test output == [4,6]
    kernel = Base.invokelatest(params,input -> 3 .* input)
    Base.invokelatest(kernel,output)
    @test output == [6,9]
    @test calls[] == 1
end

end


module InitializationCodeMotionASTTests

using Test
import GalerkinToolkit as GT

# A second AST representation exercises the pass through the existing interfaces.
struct MotionNode
    head::Symbol
    children::Vector{Any}
end
struct MotionName
    name::Symbol
end

GT.ast_head(node::MotionNode) = node.head
GT.ast_children(node::MotionNode) = node.children
GT.ast_is_leaf(::MotionNode) = false
GT.ast_is_definition(node::MotionNode) = GT.ast_head(node) === :(=)
GT.ast_is_incremental(node::MotionNode) = GT.ast_head(node) === :(+=)
GT.ast_is_call(node::MotionNode) = GT.ast_head(node) === :call
GT.ast_is_block(node::MotionNode) = GT.ast_head(node) === :block
GT.ast_is_loop(node::MotionNode) = GT.ast_head(node) === :for
GT.ast_is_index(node::MotionNode) = GT.ast_head(node) === :ref
GT.ast_is_alloc(node::MotionNode) = GT.ast_is_call(node) &&
    GT.ast_leaf_value(first(GT.ast_children(node))) in (:zeros,:alloc_zeros)
GT.ast_lhs(node::MotionNode) = first(GT.ast_children(node))
GT.ast_rhs(node::MotionNode) = last(GT.ast_children(node))
GT.ast_replace_children(node::MotionNode,children...) = MotionNode(GT.ast_head(node),Any[children...])
GT.ast_head(node::MotionName) = node.name
GT.ast_leaf_value(node::MotionName) = node.name

motion_ast(node::Expr) = MotionNode(node.head,Any[motion_ast(child) for child in node.args])
motion_ast(node::Symbol) = MotionName(node)
motion_ast(node) = node
motion_expr(node::MotionNode) = Expr(GT.ast_head(node),map(motion_expr,GT.ast_children(node))...)
motion_expr(node::MotionName) = GT.ast_leaf_value(node)
motion_expr(node) = node

@testset "Initialization code motion through AST interfaces" begin
    expr = :((values,) -> ((scale,) -> (output -> begin
        if length(values) == 0
            return
        end
        scratch = alloc_zeros("scratch",Float64,2)
        fixed = sum(values)
        coefficient = scale*fixed
        fill!(scratch,coefficient)
        copyto!(output,scratch)
    end)))
    input = motion_ast(Base.remove_linenums!(deepcopy(expr)))
    moved = GT.ast_initialization_code_motion(input)
    @test moved isa MotionNode
    @test isequal(motion_expr(moved),GT.ast_initialization_code_motion(Base.remove_linenums!(deepcopy(expr))))
    params = GT.evaluate(motion_expr(moved),([1.0,2.0],))
    output = zeros(2)
    kernel = Base.invokelatest(params,2.0)
    Base.invokelatest(kernel,output)
    @test output == [6,6]
    updated = Base.invokelatest(params,3.0)
    @test getfield(updated,:scratch) === getfield(kernel,:scratch)
    Base.invokelatest(updated,output)
    @test output == [9,9]
    empty_params = GT.evaluate(motion_expr(moved),(Float64[],))
    Base.invokelatest(Base.invokelatest(empty_params,2.0),output)
    @test output == [9,9]
end

@testset "Allocation names may be shadowed by parameters" begin
    calls = Ref(0)
    make_values = () -> (calls[] += 1; [1.0])
    expr = :(() -> ((zeros,) -> (output -> begin
        scratch = zeros()
        scratch[1] += 1
        output[1] = scratch[1]
    end)))
    params = GT.evaluate(GT.ast_initialization_code_motion(expr),())
    kernel = Base.invokelatest(params,make_values)
    @test calls[] == 0
    output = zeros(1)
    Base.invokelatest(kernel,output)
    @test output == [2]
    Base.invokelatest(kernel,output)
    @test output == [2]
    @test calls[] == 2
end

end # module InitializationCodeMotionASTTests

# End of AI-generated initialization code motion tests.

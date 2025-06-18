using LinearAlgebra

function ast_head(t)
    t
end

function ast_children(t)
    []
end

function ast_leaf_value(t)
    t
end

function ast_string(t)
    @assert ast_is_leaf(t)
    string(t)
end

function ast_children(t::Expr)
    t.args
end

function ast_dependencies(t)
    if ast_is_loop(t)
        ast_children(ast_loop_body(t))
    else
        ast_children(t)
    end
end

function ast_head(t::Expr)
    t.head
end

function ast_is_definition(t)
    false
end

function ast_is_definition(t::Expr)
    t.head === :(=)
end

function ast_is_incremental(t)
    false
end
function ast_is_incremental(t::Expr)
    t.head === :(+=)
end

function ast_is_call(t::Expr)
    t.head == :call
end

function ast_is_call(t)
    false
end

function ast_is_iterable(t)
    (t isa Expr) && ((t.head === :vect)  || (t.head === :tuple)) 
end

function ast_lhs(t::Expr)
    @assert ast_is_definition(t) || ast_is_incremental(t)
    t.args[1]
end

function ast_rhs(t::Expr)
    @assert ast_is_definition(t) || ast_is_incremental(t)
    t.args[2]
end

function ast_loop_signature_upperbound(t::Expr)
    if ast_is_loop(t)
        ast_children(ast_children(ast_loop_signature(t))[2])[3]
    elseif ast_is_definition(t)
        ast_children(ast_children(t)[2])[3]
    else
        error("arg error!")
    end
end

function ast_range(a, b)
    Expr(:call, :(:), a, b)
end

function ast_is_loop(t)
    false
end

function ast_is_loop(t::Expr)
    t.head === :for
end

function ast_loop_signature(t::Expr)
    t.args[1]
end

function ast_loop_index(t)
    a = ast_loop_signature(t)
    ast_lhs(a)
end

function ast_loop_body(t::Expr)
    t.args[2]
end

function ast_is_block(t)
    false
end

function ast_is_block(t::Expr)
    t.head === :block
end

function ast_is_leaf(t)
    length(ast_children(t)) == 0
end

function ast_is_number(t)
    t isa Number
end

function ast_is_bool(t)
    t isa Bool
end

function ast_is_nothing(t)
    return t === nothing || t == :nothing
end


function ast_is_zero(t)
    return   ((t isa Number) && t == 0.0) || (ast_is_call(t) && ast_children(t)[1] == :zero)
end


function ast_is_index(t)
    return (!ast_is_leaf(t)) && ast_head(t) == :ref
end

function ast_call(callee, args...)
    :($callee($(args...) ))
end

function ast_expr(head, args...)
    Expr(head, args...)
end

function ast_definition(lhs, rhs)
    Expr(:(=), lhs, rhs)
end

function ast_leaf(t)
    t
end

function ast_block(t=[])
    Expr(:block, t...)
end

function ast_for(range, body)
    Expr(:for, range, body)
end

function ast_index(a, indices...)
    Expr(:ref, a, indices...)
end

function ast_replace_children(t, args...)
    Expr(t.head, args...)
end

function ast_block_append_statements!(block, stmts...)
    push!(block.args, stmts...)
    return block
end


function ast_accumulate_vars(ast)
    lhs_vars = Set()
    result = Set()

    function ast_accumulate_vars_impl!(node)
        if ast_is_leaf(node)
            return 
        elseif ast_is_definition(node) || ast_is_incremental(node)
            lhs = ast_lhs(node)
            if lhs in lhs_vars
                push!(result, lhs)
            else
                push!(lhs_vars, lhs)
            end
        else
            map(ast_accumulate_vars_impl!, ast_children(node))
        end
    end

    ast_accumulate_vars_impl!(ast)
    return result
end


function ast_scope_indices(ast)
    scope_indices = IdDict()
    function visit!(node,parent_indices)
        indices = parent_indices
        if ast_is_loop(node)
            index = ast_loop_index(node)
            indices = Any[parent_indices...,index]
            scope_indices[node] = indices
        end
        for child in ast_dependencies(node)
            visit!(child,indices)
        end
        return
    end
    indices = []
    scope_indices[ast] = indices
    visit!(ast,indices)
    scope_indices
end

function ast_lhs_scopes(ast)
    @assert ast_is_block(ast)
    lhs_scopes = IdDict()
    function visit!(node,parent_scopes)
        scopes = parent_scopes
        if ast_is_definition(node)
            lhs = ast_lhs(node)
            lhs_scopes[lhs] = scopes
        end
        if ast_is_loop(node)
            scopes = Any[parent_scopes...,node]
        end
        for child in ast_dependencies(node)
            visit!(child,scopes)
        end
        return
    end
    visit!(ast,[ast])
    lhs_scopes
end

function ast_lhs_rhs(ast)
    lhs_rhs = IdDict()
    function visit!(node)
        for child in ast_dependencies(node)
            visit!(child)
        end
        if ast_is_definition(node)
            lhs = ast_lhs(node)
            rhs = ast_rhs(node)
            lhs_rhs[lhs] = rhs
        end
    end
    visit!(ast)
    lhs_rhs
end

"""
    ast_lhs_indices(ast)

"""
function ast_lhs_indices(ast)
    scope_indices = ast_scope_indices(ast)
    lhs_scopes = ast_lhs_scopes(ast)
    lhs_rhs = ast_lhs_rhs(ast)
    lhs_indices = IdDict()
    function visit!(node,scope)
        if haskey(lhs_indices,node)
            return lhs_indices[node]
        end
        if haskey(lhs_rhs,node)
            rhs = lhs_rhs[node]
            scope = lhs_scopes[node][end]
            children = ast_children(rhs)
            nested_indices = map(child->visit!(child,scope),ast_children(rhs))
            indices = unique(reduce(vcat,nested_indices;init=[]))
        else
            indices = scope_indices[scope]
            indices = intersect([node],indices)
        end
        return indices
    end
    for lhs in keys(lhs_rhs)
        scope = lhs_scopes[lhs][end]
        lhs_indices[lhs] = visit!(lhs,scope)
    end
    lhs_indices
end

"""
    ast_hoist!(ast)

Hoist loop invariant definitions in `ast`. The new AST is overwritten in `ast`.

Hypotheses:

- The root of the AST is a code block.
- The code block contains only for-loops, function calls, array indexing, definitions, and increments. No while loops, if statements, lambdas, etc.
- Functions that appear in a rhs of a definition are pure functions.
- RHSs contain a function call at most.
"""
function ast_hoist!(ast)
    lhs_scopes = ast_lhs_scopes(ast)
    scope_indices = ast_scope_indices(ast)
    lhs_indices = ast_lhs_indices(ast)
    node_done = IdDict()
    function visit!(node)
        if haskey(node_done,node)
            return
        end
        for child in ast_dependencies(node)
            visit!(child)
        end
        if ! ast_is_definition(node)
            return
        end
        lhs = ast_lhs(node)
        scopes = lhs_scopes[lhs]
        parent_scopes = scopes[1:end-1]
        current_scope = scopes[end]
        for (i,scope) in enumerate(parent_scopes)
            if issubset(lhs_indices[lhs],scope_indices[scope])
                current_children = ast_dependencies(current_scope)
                current_i = findfirst(child->child===node,current_children)
                deleteat!(current_children,current_i)
                children = ast_dependencies(scope)
                i = findfirst(child->child===scopes[i+1],children)
                insert!(children,i,node)
                node_done[node] = true
                break
            end
        end
        return
    end
    visit!(ast)
    ast
end

function ast_clean_up!(ast)
    function visit!(node)
        deps = ast_dependencies(node)
        filter!(child->!isa(child,Core.LineNumberNode),deps)
        deps = map(visit!,deps)
    end
    visit!(ast)
end


function ast_array_aliasing(ast)
    
    # TODO: this is another hard optimization. maybe we need to use the same tabulator for test and trial
end


function ast_proto_block(ast, var_count = 0)
    # TODO: remove dead code, including mem allocs.
    block = ast_block()
    var_proto = Dict()
    loop_indices = Set()

    loop_var_proto = ast_leaf(1)

    function new_var_proto()
        var_count += 1
        name = Symbol("proto_$(var_count)")
        return ast_leaf(name)
    end

    function replace_loop_var(node)
        if ast_is_definition(node)
            lhs = new_var_proto()
            rhs = replace_loop_var(ast_rhs(node))
            var_proto[ast_lhs(node)] = lhs
            return ast_definition(lhs, rhs)
        elseif ast_is_leaf(node)
            if node in loop_indices
                return loop_var_proto
            elseif haskey(var_proto, node)
                return var_proto[node]
            else
                return node
            end
        else
            new_children = map(replace_loop_var, ast_children(node))
            return ast_replace_children(node, new_children...)
        end
    end
    
    function ast_proto_block_impl(node)
        if ast_is_block(node)
            map(ast_proto_block_impl, ast_children(node))
        elseif ast_is_loop(node)
            loop_index = ast_loop_index(node)
            push!(loop_indices, loop_index)
            ast_proto_block_impl(ast_loop_body(node))
            pop!(loop_indices, loop_index)
        elseif ast_is_definition(node)
            new_stmt = replace_loop_var(node)
            ast_block_append_statements!(block, new_stmt)
        end
        
    end

    ast_proto_block_impl(ast)

    return block, var_proto, var_count
end



# a: define, b: use
function lca_deps(a::Int, b::Int)
    function low_bit(a::Int)
        a & (-a) # trick to get the lowest bit of an integer.
    end

    lca = 0
    if a == b
        lca = a
    else
        diff = a ⊻ b
        lca = (low_bit(diff) - 1) & a
    end
    result = a & ~(lca)
    return result # bitmap a - lca
end

# Note: must be applied after flatten
function ast_tabulate(ast, var_count = 0, loop_var_maxlength = Dict())
    var_dependencies = Dict() # bitset
    var_alloc_shape = Dict()
    loop_index_range = Dict()
    var_alloc_index = Dict() # loop vars indexing alloc arrays
    loop_index_depth = Dict()
    dependencies = []

    binomial_tree_blocks = [ast_block()]
    accumulate_vars = ast_accumulate_vars(ast)

    # 1 get prototype (prototype init code block & Dict original code var => prototype var )
    proto_block, var_proto, var_count = ast_proto_block(ast, var_count)


    function get_deps_impl!(node, deps_list)
        if ast_is_leaf(node)
            if haskey(var_dependencies, node)
                push!(deps_list, (node, var_dependencies[node]))
            end
        else
            for i in ast_children(node)
                get_deps_impl!(i, deps_list)
            end
        end
    end

    function get_deps(node)
        deps_list = []
        get_deps_impl!(node, deps_list)
        return deps_list
    end

    function reduce_deps(vardeps_list)
        deps_list = map(last, vardeps_list)
        reduce(|, deps_list, init=0)
    end


    # update the dependencies between loops, and also compute lca
    function update_alloc_info(var, use_deps, define_deps)
        if haskey(loop_index_depth, var)
            return # skip loop vars
        end
        # get lca 
        lca_dependencies = lca_deps(define_deps, use_deps)


        # if is alloc, then update alloc shape and alloc index
        if lca_dependencies != 0
            max_depth = length(dependencies)
            alloc_deps = []
            alloc_shape = []

            for i in max_depth:-1:1
                if 2^(i-1) & lca_dependencies > 0
                    loop_idx = dependencies[i]
                    push!(alloc_deps, loop_idx)
                    push!(alloc_shape, loop_index_range[loop_idx])
                end
            end
            
            var_alloc_shape[var] = alloc_shape
            var_alloc_index[var] = alloc_deps
        end
    end


    function ast_tabulate_impl!(node, depth = 0)
        if ast_is_loop(node)
            loop_idx = ast_loop_index(node)
            loop_signature = ast_loop_signature(node)
            loop_range = ast_loop_signature_upperbound(loop_signature)
            if ast_is_number(loop_range)
                # pass 
            elseif haskey(loop_var_maxlength, loop_idx)
                # TODO: find the exact max dim size
                loop_range = loop_var_maxlength[loop_idx]
            elseif haskey(var_proto, loop_range)
                loop_range = var_proto[loop_range] # TODO: maybe incorrect.
            else
                error("loop range not found!")
            end
            
            loop_index_depth[loop_idx] = depth
            push!(dependencies, loop_idx)
            loop_index_range[loop_idx] = loop_range
            var_dependencies[loop_idx] = 2^(depth)

            # grow binomial tree 
            size_binomial_tree = length(binomial_tree_blocks)
            new_blocks = [ast_block() for _ in 1:size_binomial_tree]
            push!(binomial_tree_blocks, new_blocks...)

            ast_tabulate_impl!(ast_loop_body(node), depth + 1)
            
            # collapse binomial tree, remove the last dep and push the node into its parent as a loop
            for i in 1:size_binomial_tree
                tabulate_loop_body = binomial_tree_blocks[i + size_binomial_tree]
                if length(ast_children(tabulate_loop_body)) > 0
                    new_loop = ast_for(loop_signature, tabulate_loop_body)
                    ast_block_append_statements!(binomial_tree_blocks[i], new_loop)
                end
            end
            binomial_tree_blocks = binomial_tree_blocks[1:size_binomial_tree]
            pop!(dependencies)
            
        elseif ast_is_block(node) 
            map(x -> ast_tabulate_impl!(x, depth), ast_children(node))
        elseif ast_is_definition(node) && ((ast_lhs(node) in accumulate_vars) || !ast_is_leaf(ast_lhs(node))) # push to the last binomial tree node
            var_dependencies[ast_lhs(node)] = (2^depth)-1
            ast_block_append_statements!(binomial_tree_blocks[end], node)
        elseif ast_is_definition(node)
            children_deps = get_deps(ast_rhs(node))
            node_deps = reduce_deps(children_deps)
            for (child_var, child_deps) in children_deps
                update_alloc_info(child_var, node_deps, child_deps)
            end
            var_dependencies[ast_lhs(node)] = node_deps
            ast_block_append_statements!(binomial_tree_blocks[node_deps + 1], node)
        else # do not hoist. these are fill!, contribute!, +=, and init loop vars. push to the last binomial tree node #TODO: do we need to hoist the entire sum loop? is it possible?
            ast_block_append_statements!(binomial_tree_blocks[end], node)
        end
    end

    function ast_tabulate_alloc_block(var_shape, var_proto)
        block = ast_block()
        var_shape_sorted = sort(collect(var_shape), by=x -> string(first(x)))
        zeros_ast = ast_leaf(:zeros)
        typeof_ast = ast_leaf(:typeof)
        for (var, shape) in var_shape_sorted
            proto = var_proto[var]
            proto_ast = ast_call(typeof_ast, proto)
            rhs_ast = ast_call(zeros_ast, proto_ast, shape...)
            stmt_ast = ast_definition(var, rhs_ast)
            ast_block_append_statements!(block, stmt_ast)
        end
        return block
    end

    function ast_tabulate_replace_tabulated_symbols(node)
        if ast_is_leaf(node)
            if haskey(var_alloc_index, node)
                ast_index(node, var_alloc_index[node]...)
            else
                node
            end
        else 
            new_children = map(ast_tabulate_replace_tabulated_symbols, ast_children(node))
            ast_replace_children(node, new_children...)
        end
    end
    
    # 2 tabulate the code, without mem alloc
    ast_tabulate_impl!(ast)
    tabulated_block = binomial_tree_blocks[1]

    # 3 allocate memory for tabulated arrays
    alloc_block = ast_tabulate_alloc_block(var_alloc_shape, var_proto)
    
    # 4 replace the tabulated symbols
    loop = ast_tabulate_replace_tabulated_symbols(tabulated_block)

    result = proto_block

    ast_block_append_statements!(result, ast_children(alloc_block)...)
    @assert ast_is_block(loop)
    ast_block_append_statements!(result, ast_children(loop)...)


    return result, var_count
end 

# struct LoopTree 
#     expr_split_by_loops::Vector
#     loop_signature
#     loop_order::Dict
#     children::Dict
#     explicit_inner_loop_signatures::Vector
# end

# function ast_tabulate_old(ast)
#     # TODO: remove empty loops
#     # TODO: variable aliasing. 
#     var_dependencies = Dict()
#     var_alloc_shape = Dict()
#     var_alloc_index = Dict() # loop vars indexing alloc arrays
#     loop_index_signature = Dict()
#     loop_index_depth = Dict()
    

#     dependency_tree = LoopTree([ast_block()], nothing, Dict(), Dict(), []) # this is a dynamic binomial tree. 
#     # children is a dict from loop signature to the loop body block
#     # loop order is the topo order of loop signatures

#     accumulate_vars = ast_accumulate_vars(ast)

#     function insert_stmt(dependencies, node)
#         pos = dependency_tree
#         for dep in dependencies
#             signature = loop_index_signature[dep]
#             if !haskey(pos.children, signature)
#                 pos.children[signature] = LoopTree([ast_block()], signature, Dict(), Dict(), [])
#             end
#             pos = pos.children[signature]
#         end
#         ast_block_append_statements!(pos.expr_split_by_loops[end], node)
#     end

#     function get_deps_impl!(node, deps_list)
#         if ast_is_leaf(node)
#             if haskey(var_dependencies, node)
#                 push!(deps_list, (node, var_dependencies[node]))
#             end
#         else
#             for i in ast_children(node)
#                 get_deps_impl!(i, deps_list)
#             end
#         end
#     end

#     function get_deps(node)
#         deps_list = []
#         get_deps_impl!(node, deps_list)
#         return deps_list
#     end

#     function reduce_deps(vardeps_list)
#         deps_list = map(last, vardeps_list)
#         if length(deps_list) == 0
#             return ()
#         end

#         dependencies = union(deps_list...)
#         dependencies_with_depth = [(loop_index_depth[dep], dep) for dep in dependencies]
#         sort!(dependencies_with_depth, by=first)
#         map(last, dependencies_with_depth)
#     end

#     # update the dependencies between loops, and also compute lca
#     function update_loop_order(var, use_deps, define_deps)
#         if haskey(loop_index_depth, var)
#             return # skip loop vars
#         end
#         lca_depth = 0
#         lca = dependency_tree
#         max_lca_depth = min(length(use_deps), length(define_deps))
#         while lca_depth < max_lca_depth && use_deps[lca_depth+1] == define_deps[lca_depth+1]
#             signature = loop_index_signature[use_deps[lca_depth+1]]
#             lca = lca.children[signature]
#             lca_depth += 1
#         end
#         if lca_depth < length(define_deps) # need to allocate a memory, insert loop order dependency
#             use_signature = loop_index_signature[use_deps[lca_depth+1]]
#             def_signature = loop_index_signature[define_deps[lca_depth+1]]
#             if !haskey(lca.loop_order, use_signature)
#                 lca.loop_order[use_signature] = Set()
#             end
#             push!(lca.loop_order[use_signature], def_signature)
#             alloc_deps = reverse(define_deps[lca_depth+1:end])
#             alloc_shape = map(x -> ast_loop_signature_upperbound(loop_index_signature[x]), alloc_deps) # TODO: find the exact maximum dim size
#             var_alloc_shape[var] = alloc_shape
#             var_alloc_index[var] = alloc_deps
#         end
#     end

#     function insert_explicit_signature_if_exists(dependencies, inner_loop_signature)
#         pos = dependency_tree
#         for i in dependencies
#             signature = loop_index_signature[i]
#             if !haskey(pos.children, signature)
#                 return
#             end
#             pos = pos.children[signature]
#         end
#         push!(pos.expr_split_by_loops, ast_block())
#         push!(pos.explicit_inner_loop_signatures, inner_loop_signature)
#     end

#     function ast_tabulate_impl(node, dependencies = [])
#         if ast_is_loop(node)
#             loop_idx = ast_loop_index(node)
#             loop_signature = ast_loop_signature(node)
#             loop_index_signature[loop_idx] = loop_signature
#             depth = length(dependencies)
#             loop_index_depth[loop_idx] = depth
#             new_dependencies = [dependencies..., loop_idx]
#             var_dependencies[loop_idx] = (loop_idx,)
            
#             ast_tabulate_impl(ast_loop_body(node), new_dependencies)
#             insert_explicit_signature_if_exists(dependencies, loop_signature)
#         elseif ast_is_block(node) 
#             map(x -> ast_tabulate_impl(x, dependencies), ast_children(node))
#         elseif ast_is_definition(node) && (ast_lhs(node) in accumulate_vars)
#             insert_stmt(dependencies, node)
#             var_dependencies[ast_lhs(node)] = dependencies
#         elseif ast_is_definition(node)
#             children_deps = get_deps(ast_rhs(node))
#             node_deps = reduce_deps(children_deps)
#             for (child_var, child_deps) in children_deps
#                 update_loop_order(child_var, node_deps, child_deps)
#             end
#             var_dependencies[ast_lhs(node)] = node_deps
#             insert_stmt(node_deps, node)
#         else # do not hoist. these are fill!, contribute!, +=, and init loop vars #TODO: do we need to hoist the entire sum loop? is it possible?
#             insert_stmt(dependencies, node)
#         end
#     end


#     # order:  block 1 -> inplicit loop 1 -> explicit loop 1 -> .... block n -> inplicit loop n -> explicit loop n -> block n + 1 -> other implicit loops
#     # order of implicit loops are obtained by topo sort. when we reconstruct an explicit loop, we also include implicit loops it depends on
#     # if a loop occurs before any assign statement or other loops, then it is treated as an implicit loop. With this strategy we do not generate empty loops
#     # an inner loop before any assign statement should not be a accumulation loop
#     # if a explicit loop signature appears twice, then loop fusion is applied and we just need to process the second loop
#     function reconstruct_loop(loop_node::LoopTree)
#         explicit_signatures = loop_node.explicit_inner_loop_signatures
#         loop_order = loop_node.loop_order
#         explicit_loop_count = Dict()
#         for i in explicit_signatures
#             if haskey(explicit_loop_count, i)
#                 explicit_loop_count[i] += 1
#             else
#                 explicit_loop_count[i] = 1
#             end
#         end

#         processed_loops = Set()
#         result_block = ast_block()

#         function reconstruct_dep_loops(signature, deps_only=false)
#             if signature in processed_loops
#                 return
#             end

#             if haskey(loop_order, signature)
#                 implicit_loop_signatures = loop_order[signature] |> collect
#                 implicit_loop_signatures_sorted = sort(implicit_loop_signatures, by=string)
#                 for s in implicit_loop_signatures_sorted
#                     reconstruct_dep_loops(s)
#                 end
#             end

#             if deps_only == false && haskey(loop_node.children, signature)
#                 loop = reconstruct_loop(loop_node.children[signature])
#                 ast_block_append_statements!(result_block, loop)
#                 push!(processed_loops, signature)
#             end
#         end


#         for i in 1:length(explicit_signatures)
#             # block i 
#             ast_block_append_statements!(result_block, ast_children(loop_node.expr_split_by_loops[i])...)
            
#             signature = explicit_signatures[i]
            
#             explicit_loop_count[signature] -= 1
#             if explicit_loop_count[signature] == 0
#                 # inplicit i 
#                 reconstruct_dep_loops(signature, true)

#                 # explicit i
#                 child = loop_node.children[signature]
#                 child_loop = reconstruct_loop(child)
#                 ast_block_append_statements!(result_block, child_loop) # append a loop
#                 push!(processed_loops, signature)
#             end
            
#         end
#         ast_block_append_statements!(result_block, ast_children(loop_node.expr_split_by_loops[end])...)

#         # other implicit loops
#         last_children = [k for (k, v) in loop_node.children]
#         for signature in last_children
#             reconstruct_dep_loops(signature)
#         end

#         result = if loop_node.loop_signature === nothing
#             result_block
#         else 
#             ast_for(loop_node.loop_signature, result_block)
#         end
#     end


#     function ast_tabulate_alloc_block(var_shape, var_proto)
#         block = ast_block()
#         var_shape_sorted = sort(collect(var_shape), by=x -> string(first(x)))
#         zeros_ast = ast_leaf(:zeros)
#         typeof_ast = ast_leaf(:typeof)
#         for (var, shape) in var_shape_sorted
#             proto = var_proto[var]
#             proto_ast = ast_call(typeof_ast, proto)
#             rhs_ast = ast_call(zeros_ast, proto_ast, shape...)
#             stmt_ast = ast_definition(var, rhs_ast)
#             ast_block_append_statements!(block, stmt_ast)
#         end
#         return block
#     end

#     function ast_tabulate_replace_tabulated_symbols!(loop_node::LoopTree)
#         for (k, child) in loop_node.children
#             ast_tabulate_replace_tabulated_symbols!(child)
#         end 

#         function replace_tabulated_symbols_expr(t)
#             if ast_is_leaf(t)
#                 if haskey(var_alloc_index, t)
#                     ast_index(t, var_alloc_index[t]...)
#                 else 
#                     t
#                 end
#             else
#                 new_children = map(replace_tabulated_symbols_expr, ast_children(t))
#                 ast_replace_children(t, new_children...)
#             end
#         end
#         for i in 1:length(loop_node.expr_split_by_loops)
#             block = loop_node.expr_split_by_loops[i]
#             new_block = replace_tabulated_symbols_expr(block)
#             loop_node.expr_split_by_loops[i] = new_block
#         end
#     end

#     @assert ast_is_block(ast)

#     # 1. tabulate without mem alloc
#     ast_tabulate_impl(ast)
    

#     # 2.a get prototype
#     proto_block, var_proto = ast_proto_block(ast)
#     display(proto_block)

#     # 2.b allocate memory for tabulated arrays
#     alloc_block = ast_tabulate_alloc_block(var_alloc_shape, var_proto)
#     display(alloc_block)
    
#     # TODO: 2.c replace the tabulated symbols
#     ast_tabulate_replace_tabulated_symbols!(dependency_tree)
    

#     # 3. reconstruct the loop
#     loop = reconstruct_loop(dependency_tree)
#     return loop
# end





function ast_flatten(ast, var_count_init = 0)
    # TODO: split flatten and topo sort?
    # in a compiler workflow, we need to allocate new variables with incremental id
    # also, we need to make the result of a workflow deterministic. 
    # so we have a var count passed from the workflow
    var_count = var_count_init
    function new_var() # TODO: decouple
        var_count += 1
        name = Symbol("var_$(var_count)")
        return ast_leaf(name)
    end


    # TODO: type: -1 -> block statement (probably array ops), 0-> always generate symbol, 1-> generate after 1 step
    function ast_flatten_impl!(node, block, expr_var, genereate_var_mode = 0)
        if ast_is_leaf(node)
            if haskey(expr_var, node)
                return expr_var[node]
            else
                return node
            end
        elseif ast_is_block(node) 
            for i in ast_children(node)
                ast_flatten_impl!(i, block, expr_var, -1)
            end
            return nothing
        elseif ast_is_loop(node)
            args = ast_children(node)
            expr = ast_for(args[1], ast_flatten_impl_block(args[2]))
            # TODO: do we need to flatten loop ranges?
            ast_block_append_statements!(block, expr)
            return expr
        elseif ast_is_definition(node) # TODO: this is not correct for general cases. when a = c and b = c appear in 2 loop levels, we cannot refer that
            lhs, rhs = ast_children(node)
            new_rhs = ast_flatten_impl!(rhs, block, expr_var, 1)
            if ast_is_leaf(new_rhs) && ast_is_leaf(lhs)
                expr_var[lhs] = new_rhs
                return new_rhs
            else
                expr = ast_replace_children(node, lhs, new_rhs)
                ast_block_append_statements!(block, expr)
                return lhs
            end
        elseif ast_is_incremental(node) # statements, assuming that we do not need to optimize the lhs as everything will be indexed by the loop vars directly.
            # the rhs is a depth-1 expression, or otherwise flatten it
            lhs, rhs = ast_children(node)
            new_rhs = ast_flatten_impl!(rhs, block, expr_var)
            expr = ast_replace_children(node, lhs, new_rhs)
            ast_block_append_statements!(block, expr)
            return lhs
        else # function calls or indexing
            args = map(ast_children(node)) do arg 
                ast_flatten_impl!(arg, block, expr_var)
            end
            new_expr = ast_replace_children(node, args...)
            if genereate_var_mode > 0
                return new_expr
            elseif genereate_var_mode == -1
                ast_block_append_statements!(block, new_expr)
                return new_expr
            elseif haskey(expr_var, new_expr)
                return expr_var[new_expr]
            else
                symbol = new_var()
                expr_var[new_expr] = symbol
                expr = ast_definition(symbol, new_expr)
                ast_block_append_statements!(block, expr)
                return symbol
            end
        end
    end

    function ast_flatten_impl_block(node)
        block = ast_block()
        for i in ast_children(node)
            ast_flatten_impl!(i, block, Dict(), -1)
        end
        return block
    end
    # the input should be a block
    @assert ast_is_block(ast)
    global_block = ast_flatten_impl_block(ast)

    global_block, var_count
end




# assuming the blocks are all integers and masks are all integer comparisons (==)
function ast_constant_folding(ast)
    # TODO: add rule: a[?] += 0 => remove stmt 
    lhs_constants = Dict()
    single_arg_functions = Set(map(ast_leaf, [+, :+, -, :-, LinearAlgebra.tr, :(LinearAlgebra.tr), sum, :sum, ])) # is there an abstract level to represent these leaf nodes? 
    multiple_arg_functions_any = Set(map(ast_leaf, [*, :*, :⋅, :dot, ⋅, ]))
    zero_calls = Set(map(ast_leaf, [zero, :zero]))

    zero_vars = Dict()
    accumulate_vars = ast_accumulate_vars(ast)
    explicit_iterables = Dict()

    function is_zero(t)
        haskey(zero_vars, t) || ast_is_zero(t)
    end

    function get_zero_proto(t)
        if haskey(zero_vars, t)
            ast_children(zero_vars[t])[2]
        elseif ast_is_zero(t) && ast_is_call(t)
            ast_children(t)[2]
        else
            t
        end
    end

    function ast_constant_folding_impl(node, is_lhs=false)
        if (!is_lhs) && haskey(lhs_constants, node)
            return lhs_constants[node]
        elseif ast_is_leaf(node)
            return node
        else
            
            head = ast_head(node)
            args_0 = ast_children(node)
            args_1 = if ast_is_definition(node) || ast_is_incremental(node)
                (ast_constant_folding_impl(args_0[1], true), map(x -> ast_constant_folding_impl(x, false), args_0[2:end])...)
            else 
                map(x -> ast_constant_folding_impl(x, false), args_0)
            end
            args = ast_is_block(args_1) ? filter(x -> x !== nothing, args_1) : args_1

            if ast_is_definition(node) && ast_is_call(args[2]) && (ast_children(args[2])[1] in zero_calls) && ast_is_leaf(ast_lhs(node)) && !(ast_lhs(node) in accumulate_vars)
                zero_vars[ast_lhs(node)] = args[2]
            end

            if ast_is_definition(node) && ast_is_iterable(args[2]) && ast_is_leaf(ast_lhs(node)) && !(ast_lhs(node) in accumulate_vars)
                explicit_iterables[ast_lhs(node)] = args[2]
            end

            if ast_is_definition(node) && ast_is_leaf(args[1]) && ast_is_leaf(args[2])
                lhs_constants[args[1]] = args[2]
                return nothing
            elseif ast_is_incremental(node) && is_zero(args[2])
                return nothing
            elseif ast_is_index(node) && ast_is_number(args[2]) # assuming that the explicit array is 1d\
                index = ast_leaf_value(args[2])
                if (ast_is_iterable(args[1]))
                    return ast_children(args[1])[index]
                elseif haskey(explicit_iterables, args[1])
                    return ast_children(explicit_iterables[args[1]])[index]
                end
            elseif ast_is_call(node)
                if args[1] == ast_leaf(:(==))
                    lhs, rhs = args[2:3]
                    if (ast_is_number(lhs) && ast_is_number(rhs)) || (ast_is_nothing(lhs) && ast_is_nothing(rhs))
                        return lhs == rhs
                    end
                elseif args[1] == ast_leaf(:ifelse) 
                    if args[2] == true
                        return args[3]
                    elseif args[2] == false
                        return args[4]
                    end
                elseif length(args) == 2 && (args[1] in single_arg_functions) && is_zero(args[2])
                    return ast_call(ast_leaf(:zero),  ast_call(args[1], get_zero_proto(args[2]) ) )
                elseif length(args) > 2 && (args[1] in multiple_arg_functions_any) && any(is_zero, args[2:end])
                    new_args = map(args[2:end]) do arg
                        if is_zero(arg)
                            get_zero_proto(arg) 
                        else
                            arg
                        end
                    end
                    return ast_call(ast_leaf(:zero), ast_call(args[1], new_args...))
                elseif args[1] == ast_leaf(:+) || args[1] == ast_leaf(+) # a + b [+...]
                    new_args = filter(x -> !is_zero(x), args[2:end])
                    if length(new_args) == 1
                        return new_args[1]
                    else
                        return ast_call(args[1], new_args...)
                    end
                elseif args[1] == ast_leaf(:-) || args[1] == ast_leaf(-) # a - b
                    if is_zero(args[2])
                        if is_zero(args[3])
                            return ast_call( ast_leaf(:zero), ast_call(ast_leaf(:-), get_zero_proto(args[2]), get_zero_proto(args[3])) )
                        else
                            return ast_call(ast_leaf(:-), args[3])
                        end
                    elseif is_zero(args[3])
                        return args[2]
                    end
                end
            elseif head == ast_leaf(:(&&))
                if all(x -> ast_is_bool(x), args)
                    return all(ast_leaf_value, args)
                end
            end
            return ast_expr(head, args...) 
        end
    
    end
    # fold +, -, *, sum, dot, tr or other linearalg operations
    # multiple arg functions = 0 with any arg = 0: dot, *
    # multiple args special case: +, -,
    # single arg: +, -, tr, sum
    ast_constant_folding_impl(ast, false)
end

# loop unroll
# lambdas to loops
# constant folding 
# flatten 
# tabulate 
# dead code elimination
# hoist
# topo sort 


# assuming that the loop range of unrolled loops are all constant integers
# unroll_loops: a set of unroll loop variables
function ast_loop_unroll(ast, unroll_loops)

    unroll_indices = []
    loop_var_unrolled_val = Dict()
    var_unrolled = Dict()

    function replace_var_unroll(node)
        if ast_is_leaf(node)
            if haskey(loop_var_unrolled_val, node)
                loop_var_unrolled_val[node]
            elseif haskey(var_unrolled, node)
                return var_unrolled[node]
            else
                node
            end
        else
            if ast_is_definition(node) && ast_is_leaf(ast_lhs(node))
                lhs = ast_lhs(node)
                suffix = map(x -> "_$x", unroll_indices)
                new_lhs_str = string(ast_string(lhs), suffix...)
                new_lhs = Symbol(new_lhs_str)
                var_unrolled[lhs] = new_lhs
                new_rhs = replace_var_unroll(ast_rhs(node))
                ast_definition(new_lhs, new_rhs)
            else
                children = map(replace_var_unroll, ast_children(node))
                ast_replace_children(node, children...)
            end
        end
    end

    # unroll_indices: a list of unroll index values. this is used to generate the suffix
    # loop_var_unrolled_val: a dict unroll index to the unrolled value. this is used to inline the loop var.
    function ast_loop_unroll_impl(node)
        if ast_is_block(node)
            children = map(ast_loop_unroll_impl, ast_children(node))
            ast_block(children)
        elseif ast_is_loop(node)
            loop_var = ast_loop_index(node)
            loop_signature = ast_loop_signature(node)
            body = ast_loop_body(node)
            if loop_var in unroll_loops
                result = ast_block()
                range = ast_loop_signature_upperbound(node)
                @assert ast_is_number(range)
                range = ast_leaf_value(range)
                push!(unroll_indices, 1)
                for i in 1:range 
                    unroll_indices[end] = i
                    loop_var_unrolled_val[loop_var] = i
                    code = ast_loop_unroll_impl(body)
                    ast_block_append_statements!(result, ast_children(code)...)
                end
                pop!(unroll_indices)
                return result
            else
                new_body = ast_loop_unroll_impl(body)
                range = ast_children(loop_signature)[2]
                new_range = ast_loop_unroll_impl(range)
                new_signature = ast_definition(loop_var, new_range)
                ast_for(new_signature, new_body)
            end
        else
            new_expr = replace_var_unroll(node)
        end
    end

    ast_loop_unroll_impl(ast)
end

# TODO: a[1] -> a_1
function ast_array_unroll(ast)
    array_def = Dict() # the last n elements represent the array shape
    array_unroll_indices = Dict()
    alloc_funcs = Set(map(ast_leaf, [:alloc_zeros, :zeros])) # TODO: decouple

    function identify_array_unrolls(node)
        if ast_is_index(node) && ast_is_leaf(ast_children(node)[1])
            # usage of arrays
            children = ast_children(node)
            var = children[1]
            indexing = children[2:end]
            unroll_mask = map(x -> ast_is_number(x) ? ast_leaf_value(x) : 0, indexing)
            if haskey(array_unroll_indices, var)
                array_unroll_indices[var] = map((x, y) -> (x == 0 || y == 0) ? 0 : max(x, y), array_unroll_indices[var], unroll_mask)
            else
                array_unroll_indices[var] = unroll_mask
            end
        elseif ast_is_definition(node) && ast_is_call(ast_rhs(node)) && (ast_children(ast_rhs(node))[1] in alloc_funcs)
            # definition of arrays. assuming that the last n args represent the array shape
            lhs = ast_lhs(node)
            array_def[lhs] = ast_rhs(node) 
        else
            map(identify_array_unrolls, ast_children(node))
        end
    end

    function ast_array_unroll_impl(node)
        if ast_is_index(node) && ast_is_leaf(ast_children(node)[1]) && haskey(array_def, ast_children(node)[1])
            # usage
            var = ast_children(node)[1]

            unroll_indices_mask = array_unroll_indices[var]
            len = length(unroll_indices_mask)
            old_indices = ast_children(node)[2:end]
            new_indices = [old_indices[i] for i in 1:len if unroll_indices_mask[i] == 0]
            unrolled_indices = [old_indices[i] for i in 1:len if unroll_indices_mask[i] != 0]
            
            var_str = ast_string(var)
            suffix_str = ["_$(i)" for i in unrolled_indices]
            new_var_str = string(var_str, suffix_str...)
            new_var = ast_leaf(Symbol(new_var_str))
            if length(new_indices) > 0
                return ast_index(new_var, new_indices...)
            else
                return new_var
            end

        elseif ast_is_definition(node) && haskey(array_def, ast_lhs(node)) && haskey(array_unroll_indices, ast_lhs(node))
            # definition of arrays. assuming that the last n args represent the array shape
            lhs = ast_lhs(node)
            rhs = array_def[lhs]
            unroll_indices = array_unroll_indices[lhs]
            len = length(unroll_indices)
            old_indices = ast_children(rhs)[end-len+1:end]
            alloc_indices = [old_indices[i] for i in 1:len if unroll_indices[i] == 0]
            unrolled_max_indices = filter(x ->x > 0, unroll_indices)
            new_rhs_children = [ast_children(rhs)[1:end-len]..., alloc_indices...]
            new_rhs = ast_replace_children(rhs, new_rhs_children...)
            
            if any(x -> x > 0, unroll_indices)
                return nothing
            elseif any(x -> x > 0, unroll_indices)
                result = ast_block()
                ranges = [1:i for i in unrolled_max_indices]
                for i in CartesianIndices((ranges...,))
                    suffices = ["_$(i[j])" for j in 1:length(unrolled_max_indices)]
                    lhs_str = ast_string(lhs)
                    new_lhs_str = string(lhs_str, reverse(suffices)...)
                    new_lhs = ast_leaf(Symbol(new_lhs_str))
                    stmt = ast_definition(new_lhs, new_rhs)
                    ast_block_append_statements!(result, stmt)
                end
                
                return result
            else
                return node
            end             
        elseif ast_is_leaf(node)
            return node
        else
            children = map(ast_array_unroll_impl, ast_children(node))
            ast_replace_children(node, children...)
        end
    end

    identify_array_unrolls(ast)
    ast_array_unroll_impl(ast)
end


# TODO: loop fusion, check correctness
function ast_loop_fusion(ast)
    function ast_loop_fusion_impl(node)
        if ast_is_block(node)
            children = ast_children(node)
            signature_count = Dict()
            signature_block = Dict()
            for child in children
                if ast_is_loop(child)
                    signature = ast_loop_signature(child)
                    signature_count[signature] =  get!(signature_count, signature, 0) + 1
                end
            end
            children = map(children) do child
                if ast_is_loop(child)
                    signature = ast_loop_signature(child)
                    signature_count[signature] -= 1
                    if !haskey(signature_block, signature)
                        signature_block[signature] = ast_block()
                    end
                    block = signature_block[signature]
                    ast_block_append_statements!(block, ast_children(ast_loop_body(child)))
                    if signature_count[signature] == 0
                        ast_for(sigature, signature_block[signature])
                    else
                        nothing
                    end
                else
                    child
                end
            end 
            children = filter(x -> x !== nothing, children)
            children = map(ast_loop_fusion_impl, children)
            ast_block(children)
        elseif ast_is_loop(node)
            signature = ast_loop_signature(node)
            body = ast_loop_fusion_impl(ast_loop_body(node))
            ast_for(signature, body)
        else
            node
        end
    
    end

    ast_loop_fusion_impl(ast)
end

# remove dead code
function ast_remove_dead_code(ast)
    # 1. nothing 
    # 2. define/alloc but not used
    # 3. else?

    used_vars = Dict()
    removed_statements = Set(map(ast_leaf, [nothing, :nothing]))

    function ast_remove_dead_code_impl(node)
        if ast_is_leaf(node)
            used_vars[node] = true
            return node
        elseif ast_is_definition(node) && ast_is_leaf(ast_lhs(node))
            lhs = ast_lhs(node)
            if haskey(used_vars, lhs) && used_vars[lhs] == true
                used_vars[lhs] = false
                ast_remove_dead_code_impl(ast_rhs(node))
                return node
            else
                return nothing
            end
        elseif ast_is_loop(node)
            ast_remove_dead_code_impl(ast_loop_signature_upperbound(node))
            signature = ast_loop_signature(node)
            body = ast_remove_dead_code_impl(ast_loop_body(node))
            return ast_for(signature, body)
        elseif ast_is_block(node)
            children = reverse(ast_children(node))
            new_children = reverse(map(ast_remove_dead_code_impl, children))
            new_children_filtered = filter(x -> !(x in removed_statements), new_children)
            return ast_block(new_children_filtered)
        else
            children = map(ast_remove_dead_code_impl, ast_children(node))
            return node
        end
    end

    ast_remove_dead_code_impl(ast)
end


# TODO: topological sort. input ast is already flattened
function ast_topological_sort(ast)
    # TODO: do we need to alias arrays?
    expr_var = Dict()
    var_depth = [Set()]
    accumulate_vars = ast_accumulate_vars(ast)
    function ast_topological_sort_impl(node, depth = 1)
        if haskey(expr_var, node)
            return expr_var[node]
        elseif ast_is_leaf(node)
            return node
        elseif ast_is_loop(node)
            push!(var_depth, Set())
            signature = ast_loop_signature(node)
            signature_children = map(x -> ast_topological_sort_impl(x, depth), ast_children(signature))
            new_signature = ast_replace_children(signature, signature_children...)
            body = ast_topological_sort_impl(ast_loop_body(node), depth + 1)
            
            for i in var_depth[end]
                delete!(expr_var, i)
            end
            pop!(var_depth)
            return ast_for(new_signature, body)
        else
            new_children = map(ast_children(node)) do child
                new_child = ast_topological_sort_impl(child, depth)
                haskey(expr_var, new_child) ? expr_var[new_child] : new_child
            end 
            new_children_filtered = ast_is_block(node) ? filter(x -> x !== nothing, new_children) : new_children
            
            new_node = ast_replace_children(node, new_children_filtered...)
            if ast_is_definition(new_node) && ast_is_leaf(ast_lhs(new_node)) && !(ast_lhs(new_node) in accumulate_vars)
                lhs = ast_lhs(new_node)
                rhs = ast_rhs(new_node)
                if ast_is_leaf(rhs)
                    expr_var[lhs] = rhs
                    push!(var_depth[depth], lhs)
                    return nothing
                else
                    expr_var[rhs] = lhs
                    push!(var_depth[depth], rhs)
                    return new_node
                end
            end
            new_node 
        end
    end
    
    ast_topological_sort_impl(ast)
end


# expr0: perfect loop 
# expr1: flatten + topo sort +
# expr2: hoist/tabulate +
# expr3: loop unroll  +
# expr4: constant folding + 
# expr : dead code elimination +

# TODO: is it possible to add information in comments?
# pass orders 1:   loop unroll  -> loop fusion -> lambdas to loops -> constant folding -> flatten -> tabulate -> topo sort (with array aliasing) -> cleanup
# pass orders 2:   lambdas to loops -> flatten -> topo sort -> loop unroll  -> loop fusion -> constant folding ->  tabulate -> loop fusion -> cleanup
# constant folding -> flatten, because we need to fold zero operations like (zero(T1) + a = a)
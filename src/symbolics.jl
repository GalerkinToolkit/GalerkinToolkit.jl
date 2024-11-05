
#macro term(expr)
#    transform(a) = a # :($(Expr(:quote,a)))
#    function transform(expr::Expr)
#        if  expr.head === :call
#            transform_call(expr)
#        elseif  expr.head === :ref
#            transform_ref(expr)
#        else
#            transform_default(expr)
#        end
#    end
#    function transform_call(expr::Expr)
#        args = map(transform,expr.args)
#        quote
#            Expr(:call,$(args...))
#        end
#    end
#    function transform_ref(expr::Expr)
#        args = map(transform,expr.args)
#        quote
#            Expr(:ref,$(args...))
#        end
#    end
#    function transform_default(expr::Expr)
#        head = expr.head
#        args = map(transform,expr.args)
#        Expr(head,args...)
#    end
#    transform(expr) |> esc
#end

macro term(expr)
    vars = Set{Symbol}()
    function findvars!(a)
        nothing
    end
    function findvars!(expr::Expr)
        if  expr.head === :(=)
            var = expr.args[1]
            push!(vars,var)
        else
            map(findvars!,expr.args)
        end
        nothing
    end
    transform(a) = a
    function transform(a::Symbol)
        if a in vars
            a
        else
            :($(Expr(:quote,a)))
        end
    end
    function transform(expr::Expr)
        if  expr.head === :call
            transform_call(expr)
        elseif  expr.head === :ref
            transform_ref(expr)
        elseif  expr.head === :$
            transform_interpolation(expr)
        elseif  expr.head === :(=)
            transform_default(expr)
        elseif  expr.head === :(->)
            transform_lambda(expr)
        elseif  expr.head === :do
            transform_do(expr)
        elseif  expr.head === :block
            transform_default(expr)
        elseif  expr.head === :tuple
            transform_tuple(expr)
        else
            error("Expr with head=$(expr.head) not supported in macro @term")
        end
    end
    function transform_call(expr::Expr)
        args = map(transform,expr.args)
        quote
            Expr(:call,$(args...))
        end
    end
    function transform_ref(expr::Expr)
        args = map(transform,expr.args)
        quote
            Expr(:ref,$(args...))
        end
    end
    function transform_lambda(expr::Expr)
        args = map(transform,expr.args)
        quote
            Expr(:(->),$(args...))
        end
    end
    function transform_do(expr::Expr)
        expr2 = Expr(:call,expr.args[1].args[1],expr.args[2],expr.args[1].args[2])
        transform(expr2)
    end
    function transform_interpolation(expr::Expr)
        expr.args[1]
    end
    function transform_tuple(expr::Expr)
        args = map(transform,expr.args)
        quote
            Expr(:tuple,$(args...))
        end
    end
    function transform_default(expr::Expr)
        head = expr.head
        args = map(transform,expr.args)
        Expr(head,args...)
    end
    findvars!(expr)
    transform(expr) |> esc
end

function topological_sort(expr,deps)
    temporary = gensym()
    expr_L = [Expr(:block) for _ in 0:length(deps)]
    marks = Dict{UInt,Symbol}()
    marks_deps = Dict{UInt,Int}()
    function visit(expr_n::Symbol)
        id_n = hash(expr_n)
        marks[id_n] = expr_n
        i = findfirst(e->expr_n===e,deps)
        j = i === nothing ? 0 : Int(i)
        marks_deps[id_n] = j
        expr_n , j
    end
    function visit(expr_n)
        id_n = hash(expr_n)
        if haskey(marks,id_n)
            if marks[id_n] !== temporary
                return marks[id_n], marks_deps[id_n]
            else
                error("Graph has cycles! This is not possible.")
            end
        end
        marks[id_n] = temporary
        var = gensym()
        if isa(expr_n,Expr)
            args = expr_n.args
            r = map(visit,args)
            args_var = map(first,r)
            j = maximum(map(last,r))
            expr_n_new = Expr(expr_n.head,args_var...)
            assignment = :($var = $expr_n_new)
        else
            j = 0
            assignment = :($var = $expr_n)
        end
        marks[id_n] = var
        marks_deps[id_n] = j
        push!(expr_L[j+1].args,assignment)
        var, j
    end
    visit(expr)
    expr_L
end

function concat_blocks(exprs)
    Expr(:block, vcat(map(x -> x.args, exprs)...)...)
end

function get_export_symbols(expr_src, exprs_dst)
    vars = Set{Symbol}()
    result = Set{Symbol}()
    function findvars!(a)
        nothing
    end
    function findvars!(expr::Expr)
        if  expr.head === :(=)
            var = expr.args[1]
            push!(vars, var)
        else
            map(findvars!, expr.args)
        end
        nothing
    end

    function findexports!(a)
        if a in vars
            push!(result, a)
        end
        nothing
    end
    function findexports!(expr::Expr)
        map(findexports!, expr.args)
        nothing
    end

    findvars!(expr_src)
    map(findexports!, exprs_dst)
    result
end

### TODO: rename it
function rewrite_export_symbols(expr_src, exprs_dst, len_loop_src, loop_idx_src, prototype_idx_src=1)
    export_symbols = get_export_symbols(expr_src, exprs_dst)
    newname_export_symbols = marks = Dict{Symbol,Symbol}()
    for s in export_symbols
        newname_export_symbols[s] = gensym()
    end


    # TODO: allocate memory on stack?
    # set the first elem

    expr_alloc = (quote 
        $v = Vector{typeof($k)}(undef, $len_loop_src)
        $v[$prototype_idx_src] = $k
    end for (k, v) in newname_export_symbols)
    
    
    expr_precomputation = quote
        $loop_idx_src = $prototype_idx_src
        $(expr_src)
        $(concat_blocks(expr_alloc))
    end

    function rewrite(a)
        if a in export_symbols
            :($(newname_export_symbols[a])[$loop_idx_src])
        else
            a 
        end
    end

    function rewrite(expr::Expr)
        Expr(expr.head, map(rewrite, expr.args)...)
    end
    
    # TODO: rewrite expr_src & expr_dst
    
    # return mem alloc, src loop, dst loop
    (expr_precomputation, rewrite(expr_src), map(rewrite, exprs_dst))
end

function topological_sort_power(expr,deps)
    # topological sort with 2^(length(deps)) dep blocks
    temporary = gensym()
    expr_L = [Expr(:block) for _ in 1:2^length(deps)]
    marks = Dict{UInt,Symbol}()
    marks_deps = Dict{UInt,Int}()
    function visit(expr_n::Symbol)
        id_n = hash(expr_n)
        marks[id_n] = expr_n
        i = findfirst(e->expr_n===e,deps)
        j = i === nothing ? 0 : 2^(Int(i)-1)
        marks_deps[id_n] = j
        expr_n , j
    end
    function visit(expr_n)
        id_n = hash(expr_n)
        if haskey(marks,id_n)
            if marks[id_n] !== temporary
                return marks[id_n], marks_deps[id_n]
            else
                error("Graph has cycles! This is not possible.")
            end
        end
        marks[id_n] = temporary
        var = gensym()
        if isa(expr_n,Expr)
            args = expr_n.args
            r = map(visit,args)
            args_var = map(first,r)
            j = reduce(|, map(last,r))
            expr_n_new = Expr(expr_n.head,args_var...)
            assignment = :($var = $expr_n_new)
        else
            j = 0
            assignment = :($var = $expr_n)
        end
        marks[id_n] = var
        marks_deps[id_n] = j
        push!(expr_L[j+1].args,assignment)
        var, j
    end
    visit(expr)
    expr_L
end

function simplify(expr)
    expr
end

const gradient = ForwardDiff.gradient
const jacobian = ForwardDiff.jacobian
const hessian = ForwardDiff.hessian
function simplify(expr::Expr)
    # r_invmap needs to be processed before r_tabulator
    # We use call since @rule is not able to identify function calls on variable slots (phi in this case)
    r_invmap = @slots a b c phi @rule call( a ∘ inverse_map_impl(phi,b) , call(phi,c) ) --> call(a,c)

    # $(\nabla (f\circ\varphi^{-1}))(\varphi(x))$ -> $(J(\varphi(x)))^{-T} (\nabla f)(x)  $
    # matrix inversion and jacobian may introduce floating point error
    r_invmap_gradient = @slots a b c phi @rule gradient(a ∘ inverse_map_impl(phi,b), call(phi,c))  -->  inv(jacobian(phi, c))' * gradient(a, c)

    # we assume to have the same type of reference face (b[e] == h[i])
    r_tabulator = @slots a b c d e f g h i @rule call(face_function(a,b,c,d,e),reference_value(f,h,i)[g]) --> face_function_value(map(reference_tabulator,a,f),b,c,d,e,g)
    r_jacobian_tabulator = @slots a b c d e f g h i @rule jacobian(face_function(a,b,c,d,e),reference_value(f,h,i)[g]) --> jacobian_face_function_value(map(gradient_reference_tabulator,a,f),b,c,d,e,g)
    

    r_shape_gradient_tabulator = @slots rid_to_fs face_to_rid face dof rid_to_coords face_to_rid2 sface point @rule gradient(face_shape_function(rid_to_fs,face_to_rid,face,dof), reference_value(rid_to_coords, face_to_rid2, sface)[point]) --> 
        face_shape_function_value(map(gradient_reference_tabulator, rid_to_fs, rid_to_coords), face_to_rid, face, dof, face_to_rid2, sface, point)

    r_shape_tabulator = @slots rid_to_fs face_to_rid face dof rid_to_coords face_to_rid2 sface point @rule call(face_shape_function(rid_to_fs,face_to_rid,face,dof), reference_value(rid_to_coords, face_to_rid2, sface)[point]) --> 
        face_shape_function_value(map(reference_tabulator, rid_to_fs, rid_to_coords), face_to_rid, face, dof, face_to_rid2, sface, point)

    
    rules = [
             r_invmap_gradient,
             r_invmap,
             r_tabulator,
             r_jacobian_tabulator,
             r_shape_tabulator,
             r_shape_gradient_tabulator,
            ]

    expr2 = nothing
    for rule in rules
        expr2 = rule(expr)
        if expr2 !== nothing
            break
        end
    end
    if expr2 === nothing
        args = map(simplify,expr.args)
        Expr(expr.head,args...)
    else
        simplify(expr2)
    end
end

function unpack_storage(dict,state)
    expr = Expr(:block)
    for k in Base.values(dict) |> collect |> sort
        push!(expr.args,:($k = $state.$k))
    end
    expr
end



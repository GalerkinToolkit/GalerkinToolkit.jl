
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
        elseif  expr.head === :.
            transform_dot(expr)
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
    function transform_dot(expr::Expr)
        expr
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
    function setup_expr_default(expr_n)
        args = expr_n.args
        r = map(visit,args)
        args_var = map(first,r)
        j = maximum(map(last,r))
        expr_n_new = Expr(expr_n.head,args_var...)
        expr_n_new, j
    end
    function setup_expr_call_or_ref(expr_n)
        args = expr_n.args[2:end]
        r = map(visit,args)
        args_var = map(first,r)
        j = maximum(map(last,r))
        expr_n_new = Expr(expr_n.head,expr_n.args[1],args_var...)
        expr_n_new, j
    end
    function setup_expr_lambda(expr_n)
        body = expr_n.args[2]
        body_sorted = topological_sort(body,())[1]
        j = length(deps)
        expr_n_new = Expr(expr_n.head,expr_n.args[1],body_sorted)
        expr_n_new, j
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
            if expr_n.head === :call || expr_n.head === :ref
                expr_n_new, j = setup_expr_call_or_ref(expr_n)
            elseif expr_n.head === :(->)
                expr_n_new, j = setup_expr_lambda(expr_n)
            else
                expr_n_new, j = setup_expr_default(expr_n)
            end
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

#const gradient = ForwardDiff.gradient
#const jacobian = ForwardDiff.jacobian
#const hessian = ForwardDiff.hessian
#function simplify(expr::Expr)
#    # r_invmap needs to be processed before r_tabulator
#    # We use call since @rule is not able to identify function calls on variable slots (phi in this case)
#    r_invmap = @slots a b c phi @rule call( a ∘ inverse_map_impl(phi,b) , call(phi,c) ) --> call(a,c)
#
#    # $(\nabla (f\circ\varphi^{-1}))(\varphi(x))$ -> $(J(\varphi(x)))^{-T} (\nabla f)(x)  $
#    # matrix inversion and jacobian may introduce floating point error
#    r_invmap_gradient = @slots a b c phi @rule gradient(a ∘ inverse_map_impl(phi,b), call(phi,c))  -->  inv(jacobian(phi, c))' * gradient(a, c)
#
#    # we assume to have the same type of reference face (b[e] == h[i])
#    r_tabulator = @slots a b c d e f g h i @rule call(face_function(a,b,c,d,e),reference_value(f,h,i)[g]) --> face_function_value(map(reference_tabulator,a,f),b,c,d,e,g)
#    r_jacobian_tabulator = @slots a b c d e f g h i @rule jacobian(face_function(a,b,c,d,e),reference_value(f,h,i)[g]) --> jacobian_face_function_value(map(gradient_reference_tabulator,a,f),b,c,d,e,g)
#    
#
#    r_shape_gradient_tabulator = @slots rid_to_fs face_to_rid face dof rid_to_coords face_to_rid2 sface point @rule gradient(face_shape_function(rid_to_fs,face_to_rid,face,dof), reference_value(rid_to_coords, face_to_rid2, sface)[point]) --> 
#        face_shape_function_value(map(gradient_reference_tabulator, rid_to_fs, rid_to_coords), face_to_rid, face, dof, face_to_rid2, sface, point)
#
#    r_shape_tabulator = @slots rid_to_fs face_to_rid face dof rid_to_coords face_to_rid2 sface point @rule call(face_shape_function(rid_to_fs,face_to_rid,face,dof), reference_value(rid_to_coords, face_to_rid2, sface)[point]) --> 
#        face_shape_function_value(map(reference_tabulator, rid_to_fs, rid_to_coords), face_to_rid, face, dof, face_to_rid2, sface, point)
#
#    
#    rules = [
#             r_invmap_gradient,
#             r_invmap,
#             r_tabulator,
#             r_jacobian_tabulator,
#             r_shape_tabulator,
#             r_shape_gradient_tabulator,
#            ]
#
#    expr2 = nothing
#    for rule in rules
#        expr2 = rule(expr)
#        if expr2 !== nothing
#            break
#        end
#    end
#    if expr2 === nothing
#        args = map(simplify,expr.args)
#        Expr(expr.head,args...)
#    else
#        simplify(expr2)
#    end
#end

function simplify(expr)
    expr
end

function simplify(expr::Expr)
    expr = rewrite(expr)
    expr = inline_lambdas(expr)
end

function rewrite(expr)
    expr
end

function rewrite(expr::Expr)
    #r_map_getindex_enumerate = @slots f r i @rule getindex(map(f,enumerate(r)),i) --> call(f,(i,r[i]))
    r_map_getindex_1 = @slots f r i @rule getindex(map(f,r),i) --> call(f,r[i])
    r_map_getindex_2 = @slots f r s i @rule getindex(map(f,r,s),i) --> call(f,r[i],s[i])
    rules = [
             #r_map_getindex_enumerate,
             r_map_getindex_1,
             r_map_getindex_2,
            ]

    expr2 = nothing
    for rule in rules
        expr2 = rule(expr)
        if expr2 !== nothing
            break
        end
    end
    if expr2 === nothing
        args = map(rewrite,expr.args)
        Expr(expr.head,args...)
    else
        rewrite(expr2)
    end
end

function inline_lambdas(expr)
    MacroTools.postwalk(expr) do ex
        if !isa(ex,Expr)
            return ex
        end
        if ex.head !== :call
            return ex
        end
        if ex.args[1] !== :call
            return ex
        end
        lambda = ex.args[2]
        if !isa(lambda,Expr)
            return ex
        end
        if lambda.head !== :(->)
            return ex
        end
        signature, body = lambda.args
        if isa(signature,Symbol)
            x = ex.args[3]
            return substitute(body,signature=>x)
        end
        if !isa(signature,Expr)
            return ex
        end
        if signature.head === :tuple
            actual_x = ex.args[3:end]
            dummy_x = signature.args
            for i in 1:length(dummy_x)
                if isa(dummy_x[i],Symbol)
                    body = substitute(body,dummy_x[i]=>actual_x[i])
                end
            end
            return body
        end
        ex
    end
end

function substitute(expr,old_new)
    old,new = old_new
    MacroTools.postwalk(ex -> ex === old ? new : ex ,expr)
end

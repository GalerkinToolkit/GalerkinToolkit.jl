
macro term(expr)
    transform(a) = a # :($(Expr(:quote,a)))
    function transform(expr::Expr)
        if  expr.head === :call
            transform_call(expr)
        elseif  expr.head === :ref
            transform_ref(expr)
        else
            transform_default(expr)
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
    function transform_default(expr::Expr)
        head = expr.head
        args = map(transform,expr.args)
        Expr(head,args...)
    end
    transform(expr) |> esc
end

function topological_sort(expr,deps)
    temporary = gensym()
    expr_L = [Expr(:block) for _ in 0:length(deps)]
    marks = Dict{UInt,Symbol}()
    marks_deps = Dict{UInt,Int}()
    function visit(expr_n::Symbol)
        id_n = objectid(expr_n)
        marks[id_n] = expr_n
        i = findfirst(e->expr_n===e,deps)
        j = i === nothing ? 0 : Int(i)
        marks_deps[id_n] = j
        expr_n , j
    end
    function visit(expr_n)
        id_n = objectid(expr_n)
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

function simplify(expr)
    expr
end

function simplify(expr::Expr)
    r001 = @slots a b c d e f g @rule face_function(a,b,c,d,e)(reference_value(f,b,e)[g]) --> face_function_value(reference_tabulators(a,f),b,c,d,e,g)
    expr2 = r001(expr)
    if expr2 === nothing
        args = map(simplify,expr.args)
        Expr(expr.head,args...)
    else
        expr2
    end
end

function unpack_state(dict,state)
    expr = Expr(:block)
    for k in Base.values(dict) |> collect |> sort
        push!(expr.args,:($k = $state.$k))
    end
    expr
end



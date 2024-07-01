
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

function topological_sort(expr)
    temporary = gensym()
    expr_L = Expr(:block)
    marks = Dict{UInt,Symbol}()
    function visit(expr_n)
        id_n = objectid(expr_n)
        if haskey(marks,id_n)
            if marks[id_n] !== temporary
                return marks[id_n]
            else
                error("Graph has cycles! This is not possible.")
            end
        end
        marks[id_n] = temporary
        var = gensym()
        if isa(expr_n,Expr)
            args = expr_n.args
            args_var = map(visit,args)
            expr_n_new = Expr(expr_n.head,args_var...)
            assignment = :($var = $expr_n_new)
        else
            assignment = :($var = $expr_n)
        end
        marks[id_n] = var
        push!(expr_L.args,assignment)
        var
    end
    visit(expr)
    expr_L
end

module TMP

using ForwardDiff

struct Var end
struct Call end
struct Lambda end
struct Expansion end
struct Index end

struct Term{A}
    head::A
    args::Vector{Any}
end

variable() = Term(Var(),Any[])

call(f,args...) = Term(Call(),Any[f,args...])

lambda(var,term) = Term(Lambda(),Any[var,term])

expand(term,var,range) = Term(Expansion(),Any[term,var,range])

# We need the range to be able to "undo" an expand
# Really? Think about this
index(term,i,range) = Term(Index(),Any[term,i,range])

function substitute(term::Term{Var},var::Term{Var},val)
    term === var ? val : term
end

function substitute(term::Term,var::Term{Var},val)
    args = Any[ substitute(t,var,val) for t in term.args ]
    Term(term.head,args)
end

function substitute(term,var::Term{Var},val)
    term
end

function materialize(term)
    term
end

function materialize(term::Term{Var})
    error("A variable cannot be materalized")
end

# TODO don't use args
# use getters instead

function materialize(term::Term{Call})
    #TODO This is not efficient
    # Convert to DAG
    # compute topological ordering
    # Generate code
    # call code
    args = map(materialize,term.args[2:end])
    f = materialize(term.args[1])
    f(args...)
end

function materialize(term::Term{Lambda})
    #TODO This is not efficient
    # Convert to DAG
    # compute topological ordering
    # Generate code
    x -> begin
        substitute(term.args[2],term.args[1],x) |> materialize
    end
end

# TODO in-place version
function materialize(term::Term{Expansion})
    #TODO This is not efficient
    # Convert to DAG
    # compute topological ordering
    # Generate code
    l = lambda(term.args[2],term.args[1])
    f = materialize(l)
    r = materialize(term.args[3])
    [f(i) for i in r]
end

function materialize(term::Term{Index})
    #TODO This is not efficient
    # Convert to DAG
    # compute topological ordering
    # Generate code
    a = materialize(term.args[1])
    i = materialize(term.args[2])
    r = materialize(term.args[3])
    @boundscheck @assert i in r
    a[i]
end

x = variable()
y = 3
t = call(+,x,y)
dump(t)

t = substitute(t,x,50.)
dump(t)

@show materialize(t)

t = call(+,x,y)

t = lambda(x,t)

f = materialize(t)
@show f(50.)

t = call(+,x,y)
t = expand(t,x,1:4)
@show materialize(t)

t = index(t,3,1:4)
@show materialize(t)

nldofs = 3
npoints = 4
ncells = 5

sfuns = [ (x) -> sum(x.*i) for i in 1:nldofs ]
nodes = [(0,0),(1,0),(0,1)]
qp = [(1,1),(2,2),(3,3),(4,4)]
coords = [ [ (rand(),rand()) for _ in 1:nldofs] for _ in 1:ncells]

i = variable()
q = variable()
j = variable()
c = variable()
k = variable()
x = variable()

s = call(index(sfuns,i,1:nldofs),q)
t = lambda(q,expand(s,i,1:nldofs))
f = materialize(t)
@show f((1,2))
@show f((2,4))

f = variable()
σ = lambda(f,call(f,call(getindex,nodes,k)))

xck = call(getindex,call(getindex,coords,c),k)
sk = index(expand(s,i,1:nldofs),k,1:nldofs)
ϕ = call(sum,lambda(k,call(*,xck,sk)),1:nldofs)

c = 100
t = index(ϕ,c,1:ncells)
t = lambda(q,t)
f = materialize(t)
@show f((3,2))


aaaaa


invϕ = call(call(invmap,lambda(q,ϕ)),x)

a = variable()

u = call(lambda(q,lambda(a,s)),invϕ)



#∇ϕ = call(call(ForwardDiff.jacobian,lambda(q,ϕ)),q)


dump(σ)
dump(ϕ)




end # module

module TMP

using ForwardDiff
using StaticArrays
using LinearAlgebra

struct Var end
struct Call end
struct Lambda end
struct Expand end

struct Term{A}
    head::A
    args::Vector{Any}
end

variable() = Term(Var(),Any[])

call(f,args...) = Term(Call(),Any[f,args...])

lambda(var,term) = Term(Lambda(),Any[var,term])

expand(term,var,range) = Term(Expand(),Any[term,var,range])

function substitute(term::Term{Var},vars,vals)
    for (var,val) in zip(vars,vals)
        term === var && return val
    end
    term
end

function substitute(term::Term,vars,vals)
    args = Any[ substitute(t,vars,vals) for t in term.args ]
    Term(term.head,args)
end

function substitute(term,vars,vals)
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
    (x...) -> begin
        substitute(term.args[2],term.args[1],x) |> materialize
    end
end

# TODO in-place version
function materialize(term::Term{Expand})
    #TODO This is not efficient
    # Convert to DAG
    # compute topological ordering
    # Generate code
    @assert length(term.args[2]) == length(term.args[3])
    l = lambda(term.args[2],term.args[1])
    f = materialize(l)
    r = map(materialize,term.args[3])
    length(r) == 1 && return [f(i) for i in r[1]]
    length(r) == 2 && return [f(i,j) for i in r[1], j in r[2]]
    error("expand up to 2 dimensions")
end

x = variable()
y = 3
t = call(+,x,y)
dump(t)

t = substitute(t,[x],[50.])
dump(t)

@show materialize(t)

t = call(+,x,y)

t = lambda([x],t)

f = materialize(t)
@show f(50.)

t = call(+,x,y)
t = expand(t,[x],[1:4])
@show materialize(t)

i = variable()
j = variable()
t = call(*,i,j)
t = expand(t,[i,j],[1:4,1:6])
@show materialize(t) |> size

nldofs = 3
npoints = 4
ncells = 5

nodes = SVector{2,Float64}[(0,0),(1,0),(0,1)]
points = SVector{2,Float64}[(0,0),(1,0),(0,1),(1,1)]
cell_coords = [ rand(SVector{2,Float64},nldofs) for i in 1:ncells]
monomials = [ (x) -> prod(x.^nodes[i]) for i in 1:nldofs ]
q = SVector(2,3)
monomials[1](q)
ForwardDiff.gradient(monomials[2],q)

k = variable()
mk = call(getindex,monomials,k)
m = expand(mk,[k],[1:nldofs])

materialize(m)[1](q)
ForwardDiff.gradient(materialize(m)[1],q)

k = variable()
f = variable()
σ = expand(lambda([f],call(f,call(getindex,nodes,k))),[k],[1:nldofs])

i = variable()
j = variable()
mi = call(getindex,m,i)
σj = call(getindex,σ,j)
Cij = call(σj,mi)

C = expand(Cij,[i,j],[1:nldofs,1:nldofs])
@show materialize(C)

D = call(\,C,I)

k = variable()
q = variable()
mqk = call(call(getindex,m,k),q)
mq = expand(mqk,[k],[1:nldofs])

sq = call(*,D,mq)
k = variable()
sqk = call(getindex,sq,k)
sk = lambda([q],sqk)
s = expand(sk,[k],[1:nldofs])

e = variable()
k = variable()
q = variable()
xek = call(getindex,call(getindex,cell_coords,e),k)
sk = call(getindex,s,k)
skq = call(sk,q)
ϕeq = call(sum,lambda([k],call(*,xek,skq)),1:nldofs)
ϕe = lambda([q],ϕeq)
ϕ = expand(ϕe,[e],[1:ncells])

display( materialize(ϕ)[1](SVector(0,0)))

e = variable()
q = variable()
qx = call(getindex,points,q)
ϕe = call(getindex,ϕ,e)
∇ϕeq = call(ForwardDiff.jacobian,ϕe,qx)
invJeq = call(inv,∇ϕeq)
detJeq = call(det,∇ϕeq)

j = variable()
i = variable()
ujq = call(call(getindex,s,j),qx)
viq = call(call(getindex,s,i),qx)
∇ujq = call(dot,ujq,invJeq)
∇viq = call(dot,viq,invJeq)

Keqij = call(*,call(⋅,∇ujq,∇viq),detJeq)
Keij = call(sum,lambda([q],Keqij),1:npoints)

Ke = expand(Keij,[i,j],[1:nldofs,1:nldofs])
t = lambda([e],Ke)
f = materialize(t)
f(1)

xxx
#call(sum,lambda([k],,1)


aaa



sfuns = [ (x) -> sum(x*i) for i in 1:nldofs ]


qp = SVector{2,Float64}[(1,1),(2,2),(3,3),(4,4)]
coords = [ [ SVector{2,Float64}(rand(),rand()) for _ in 1:nldofs] for _ in 1:ncells]

i = variable()
q = variable()
j = variable()
c = variable()
k = variable()
x = variable()

s = call(call(getindex,sfuns,i),q)
t = lambda(q,expand(s,i,1:nldofs))
f = materialize(t)
@show f(SVector(1,2))
@show f(SVector(2,4))

f = variable()
σ = lambda(f,call(f,call(getindex,nodes,k)))

xck = call(getindex,call(getindex,coords,c),k)
sk = call(getindex,expand(s,i,1:nldofs),k)
ϕ = call(sum,lambda(k,call(*,xck,sk)),1:nldofs)

t = call(getindex,expand(ϕ,c,1:ncells),4)
t = lambda(q,t)
f = materialize(t)
@show f(SVector(3,2))


function invmap end

invϕ = call(call(invmap,lambda(q,ϕ)),x)

#a = variable()
#u = call(lambda(q,call(getindex,expand(s,i,1:nldofs),a),invϕ)
#b = variable()
#v = call(lambda(q,call(getindex,expand(s,i,1:nldofs),v),invϕ)



#∇ϕ = call(call(ForwardDiff.jacobian,lambda(q,ϕ)),q)


dump(σ)
dump(ϕ)




end # module

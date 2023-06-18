module TMP

using ForwardDiff
using StaticArrays
using LinearAlgebra
using SparseArrays

struct Var end
struct Call end
struct Lambda end
struct Expand end
struct Assemble end

struct Term{A}
    head::A
    args::Vector{Any}
end

variable() = Term(Var(),Any[])

call(f,args...) = Term(Call(),Any[f,args...])

lambda(var,term) = Term(Lambda(),Any[var,term])

expand(term,var,range) = Term(Expand(),Any[term,var,range])

assemble(mat,rows,cols,var,range) = Term(Assemble(),Any[mat,rows,cols,var,range])

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

function materialize!(A,term::Term{Assemble})
    #TODO This is not efficient
    # This the full assembly code needs to be generated
    # We cannot call invokelatest at each cell
    tmat = lambda(term.args[4],term.args[1])
    trows = lambda(term.args[4],term.args[2])
    tcols = lambda(term.args[4],term.args[3])
    fmat = materialize(tmat)
    frows = materialize(trows)
    fcols = materialize(tcols)
    r = map(materialize,term.args[5])
    @assert length(r) == 1 # In the future it can be r>1, e.g. for HDG
    LinearAlgebra.fillstored!(A,zero(eltype(A)))
    for e in r[1]
        mat = fmat(e) # TODO in-place
        rows = frows(e) # TODO in-place
        cols = fcols(e) # TODO in-place
        A[rows,cols] .+= mat
    end
    A
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
ngdofs = 7
cell_dofs = [ rand(1:ngdofs,nldofs) for i in 1:ncells]
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
∇sqk = call(ForwardDiff.gradient,sk,q)
s = expand(sk,[k],[1:nldofs])
∇sk = lambda([q],∇sqk)
∇s = expand(∇sk,[k],[1:nldofs])

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

t = lambda([e,q],∇ϕeq)
f = materialize(t)
@show f(4,2)

j = variable()
i = variable()
∇ujq_ref = call(call(getindex,∇s,j),qx)
∇viq_ref = call(call(getindex,∇s,i),qx)
∇uejq = call(*,call(transpose,∇ujq_ref),invJeq)
∇veiq = call(*,call(transpose,∇viq_ref),invJeq)

t = lambda([e,q,i],∇veiq)
f = materialize(t)
@show f(5,2,1)

Keqij = call(*,call(⋅,∇uejq,∇veiq),detJeq)
Keij = call(sum,lambda([q],Keqij),1:npoints)

Ke = expand(Keij,[i,j],[1:nldofs,1:nldofs])
t = lambda([e],Ke)
f = materialize(t)
display(f(5))

dofse = call(getindex,cell_dofs,e)
A = spzeros(ngdofs,ngdofs)
t = assemble(Ke,dofse,dofse,[e],[1:ncells])
dump(t)
materialize!(A,t)
display(A)





end # module


function q_monomial_basis(orders::Tuple)
  cis = CartesianIndices(orders.+1)
  exps = [ Tuple(cis[i]).-1  for i in 1:length(cis) ]
  map(Monomial,exps)
end

function p_monomial_basis(orders::Tuple)
  k = first(orders)
  @boundscheck @assert all(i->i==k,orders) msg
  m = q_monomial_basis(orders)
  filter(mi->sum(mi.exps)<=k,m)
end

function p̃_monomial_basis(orders::Tuple)
  k = first(orders)
  @boundscheck @assert all(i->i==k,orders)
  m = p_monomial_basis(orders)
  filter(i->sum(i.exps)>(k-1),m)
end

function s̃_basis(orders::Tuple{Any,Any})
  k = first(orders)
  @boundscheck @assert all(i->i==k,orders)
  e1 = SVector(1,0)
  e2 = SVector(0,1)
  function f(α)
    m1 = Monomial(α-1,k-α+1)
    m2 = Monomial(α,k-α)
    linear_combination([-e1,e2],[m1,m2])
  end
  [ f(α) for α in 1:k ]
end

function s̃_basis(orders::Tuple{Any,Any,Any})
  k = first(orders)
  @boundscheck @assert all(i->i==k,orders)
  e1 = SVector(1,0,0)
  e2 = SVector(0,1,0)
  e3 = SVector(0,0,1)
  function a(α,β)
    m1 = Monomial(α-1,k-α-β-2,β-1)
    m2 = Monomial(α,k-α-β-1,β-1)
    linear_combination([-e1,e2],[m1,m2])
  end
  function b(α,β)
    m1 = Monomial(k-α-β+1,β-1,α)
    m3 = Monomial(k-α-β+2,β-1,α-1)
    linear_combination([-e1,e3],[m1,m3])
  end
  function c(α)
    m2 = Monomial(0,α-1,k-α+1)
    m3 = Monomial(0,α,k-α)
    linear_combination([-e2,e3],[m2,m3])
  end
  m = [ c(α) for α in 1:k]
  for β in 1:k
    for α in 1:(k+1-β)
      push!(m,a(α,β))
      push!(m,b(α,β))
    end
  end
  m
end

function q_equispaced_nodes(orders::Tuple)
  monomials = q_monomial_basis(orders)
  map(monomials) do monomial
    SVector(monomial.exps./orders)
  end
end

function p_equispaced_nodes(orders::Tuple)
  monomials = p_monomial_basis(orders)
  map(monomials) do monomial
    SVector(monomial.exps./orders)
  end
end





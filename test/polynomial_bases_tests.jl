module PolynomialBasesTests
using GalerkinToolkit

nodes = p_equispaced_nodes((2,2))
nodes = q_equispaced_nodes((2,4))

basis = p_monomial_basis((1,1))
vector_basis = cartesian_product(basis,basis,basis)
vector_basis = vector_valued_basis(basis,Val(3))



x = (1,2)

r = evaluate(vector_basis,x)
display(r)








end

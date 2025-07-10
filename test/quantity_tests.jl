module QuantityTests

import GalerkinToolkit as GT
import GalerkinToolkit
import PartitionedArrays as pa
using Test
using WriteVTK
using LinearAlgebra
using StaticArrays
using ForwardDiff
using AbstractTrees

outdir = mkpath(joinpath(@__DIR__,"..","output"))

domain = (0,1,0,1)
cells = (4,4)
mesh = GT.cartesian_mesh(domain,cells)

D = GT.num_dims(mesh)
fs = GT.reference_spaces(mesh,D-1)
GT.label_interior_faces!(mesh;physical_name="interior_faces")
GT.label_boundary_faces!(mesh;physical_name="boundary_faces")

Ω = GT.interior(mesh)
Ωref = GT.interior(mesh;is_reference_domain=true)

D = GT.num_dims(mesh)
Γref = GT.boundary(mesh;
                 is_reference_domain=true,
                 group_names=["boundary_faces"])

Γ = GT.physical_domain(Γref)

Λref = GT.skeleton(mesh;
                 is_reference_domain=true,
                 group_names=["interior_faces"])

Λ = GT.physical_domain(Λref)

@test Ω == Ω
@test Ω != Ωref

q1 = GT.face_quantity([ f for f in 1:GT.num_faces(Ω)],Ω)
q2 = GT.face_quantity([ 10*f for f in 1:GT.num_faces(Ω)],Ω)
q3 = GT.face_quantity([ 100*f for f in 1:GT.num_faces(Γ)],Γ)
q4 = GT.face_quantity([ 1000*f for f in 1:GT.num_faces(Λ)],Λ)

@test GT.is_boundary(Ω) == false
@test GT.is_boundary(Γ) == true
@test GT.is_boundary(Λ) == false

q = q1 + q2
index = GT.generate_index(Ω)
faces = GT.get_symbol!(index,GT.faces(Ω),"faces")
t = GT.term(q,index)
print_tree(t)
expr = GT.expression(t)
@test GT.free_dims(t) == [D]
storage = GT.index_storage(index)
expr = quote
    $(GT.unpack_index_storage(index,:storage))
    $(GT.face_index(index,D)) = $faces[3]
    $(GT.topological_sort(expr,())[1])
end
display(expr)
r = eval(expr)
@test r == 33

q = q1 + q3
index = GT.generate_index(Γ)
faces = GT.get_symbol!(index,GT.faces(Γ),"faces")
t = GT.term(q,index)
print_tree(t)
expr = GT.expression(t)
@test GT.free_dims(t) == [D-1]
storage = GT.index_storage(index)
expr = quote
    $(GT.unpack_index_storage(index,:storage))
    $(GT.face_index(index,D-1)) = $faces[2]
    $(GT.topological_sort(expr,())[1])
end
display(expr)
r = eval(expr)
@test r == 201


q = (q1 + q4)[2]
index = GT.generate_index(Λ)
faces = GT.get_symbol!(index,GT.faces(Λ),"faces")
t = GT.term(q,index)
print_tree(t)
@test GT.free_dims(t) == [D-1]
expr = GT.expression(t)
storage = GT.index_storage(index)
expr = quote
    $(GT.unpack_index_storage(index,:storage))
    $(GT.face_index(index,D-1)) = $faces[2]
    $(GT.topological_sort(GT.simplify(expr),())[1])
end
display(expr)
r = eval(expr)
@test r == 2002

#q = (q1 + q4)[2]
#index = GT.generate_index(Λ)
#faces = GT.get_symbol!(index,GT.faces(Λ),"faces")
#t = GT.term(q,index)
#print_tree(t)
#@test GT.free_dims(t) == [D-1]
#expr = GT.expression(t)
#storage = GT.index_storage(index)
#expr = quote
#    $(GT.unpack_index_storage(index,:storage))
#    $(GT.face_index(index,D-1)) = $faces[2]
#    $(GT.topological_sort(GT.simplify(expr),())[1])
#    #$(GT.topological_sort(t.expr,())[1])
#end
#display(expr)
#r = eval(expr)
#@test r == 2002

#q = (q1 + q4)[2]
#form_arity = 1
#index = GT.generate_index(Λ,form_arity)
#faces = GT.get_symbol!(index,GT.faces(Λ),"faces")
#t = GT.term(q,index)
#print_tree(t)
#@test GT.free_dims(t) == [D-1]
#expr = GT.expression(t)
#storage = GT.index_storage(index)
#expr = quote
#    $(GT.unpack_index_storage(index,:storage))
#    $(GT.face_index(index,D-1)) = $faces[2]
#    $(GT.face_around_index(index,1)) = 1
#    #$(GT.topological_sort(GT.simplify(expr),())[1])
#    $(GT.topological_sort(expr,())[1])
#end
#display(expr)
#r = eval(expr)
#@test r == 0

rid_to_dof_to_s = [[x -> dof*sum(x) for dof in 1:5]]
dof = gensym("dummy-dof")
s2 = GT.shape_function_quantity(rid_to_dof_to_s,Ω;reference=true,dof)
x2 = GT.point_quantity([SVector{2,Float64}[[0,0],[1,1]]],Ω;reference=true)

q = s2(x2)
index = GT.generate_index(Ω)
faces = GT.get_symbol!(index,GT.faces(Ω),"faces")
t = GT.term(q,index)
print_tree(t)
expr = GT.expression(t)
@test GT.free_dims(t) == [D]
storage = GT.index_storage(index)
expr = quote
    $(GT.unpack_index_storage(index,:storage))
    $(GT.face_index(index,D)) = $faces[3]
    $dof = 5
    $(GT.point_index(index)) = 2
    $(GT.topological_sort(expr,())[1])
end
display(expr)
r = eval(expr)
@test r == 10

q = ForwardDiff.gradient(s2,x2)
index = GT.generate_index(Ω)
faces = GT.get_symbol!(index,GT.faces(Ω),"faces")
t = GT.term(q,index)
print_tree(t)
expr = GT.expression(t)
@test GT.free_dims(t) == [D]
storage = GT.index_storage(index)
expr = quote
    $(GT.unpack_index_storage(index,:storage))
    $(GT.face_index(index,D)) = $faces[3]
    $dof = 5
    $(GT.point_index(index)) = 2
    $(GT.topological_sort(expr,())[1])
end
display(expr)
r = eval(expr)
@test r == [5,5]

g_to_v = sum.(GT.node_coordinates(mesh))
u2 = GT.quantity() do index
    g_to_value = GT.get_symbol!(index,g_to_v,"g_to_value")
    face_to_dof_to_g = GT.get_symbol!(index,GT.face_nodes(mesh,D),"face_to_dof_to_g")
    face = GT.face_index(index,D)
    expr = :($g_to_value[$face_to_dof_to_g[$face][$dof]])
    p = zero(eltype(g_to_v))
    coeffs = GT.expr_term([D],expr,p,index)
    funs = GT.term(s2,index)
    expr = :(length($face_to_dof_to_g[$face]))
    ndofs = GT.expr_term([D],expr,0,index)
    GT.discrete_function_term(coeffs,funs,dof,ndofs)
end

q = ForwardDiff.gradient(u2,x2)

index = GT.generate_index(Ω)
faces = GT.get_symbol!(index,GT.faces(Ω),"faces")
t = GT.term(q,index)
print_tree(t)
expr = GT.expression(t)
@test GT.free_dims(t) == [D]
storage = GT.index_storage(index)
expr = quote
    $(GT.unpack_index_storage(index,:storage))
    $(GT.face_index(index,D)) = $faces[3]
    $(GT.point_index(index)) = 2
    $(GT.topological_sort(expr,())[1])
end
display(expr)
r = eval(expr)
@test r == [8.25, 8.25]

u2 = GT.physical_map(mesh,D)

q = u2(x2)
index = GT.generate_index(Ω)
faces = GT.get_symbol!(index,GT.faces(Ω),"faces")
t = GT.term(q,index)
print_tree(t)
expr = GT.expression(t)
@test GT.free_dims(t) == [D]
storage = GT.index_storage(index)
expr = quote
    $(GT.unpack_index_storage(index,:storage))
    $(GT.face_index(index,D)) = $faces[3]
    $(GT.point_index(index)) = 2
    $(GT.topological_sort(expr,())[1])
end
display(expr)
r = eval(expr)
@test r == [0.75, 0.25]

q = ForwardDiff.jacobian(u2,x2)
index = GT.generate_index(Ω)
faces = GT.get_symbol!(index,GT.faces(Ω),"faces")
t = GT.term(q,index)
print_tree(t)
expr = GT.expression(t)
@test GT.free_dims(t) == [D]
storage = GT.index_storage(index)
expr = quote
    $(GT.unpack_index_storage(index,:storage))
    $(GT.face_index(index,D)) = $faces[3]
    $(GT.point_index(index)) = 2
    $(GT.topological_sort(expr,())[1])
end
display(expr)
r = eval(expr)
@test r == [0.25 0.0; 0.0 0.25]

u2 = GT.physical_map(mesh,D)
u2inv = GT.inverse_physical_map(mesh,D)
q = (s2∘u2inv)(u2(x2))
index = GT.generate_index(Ω)
faces = GT.get_symbol!(index,GT.faces(Ω),"faces")
t = GT.term(q,index)
print_tree(t)
expr = GT.expression(t)
@test GT.free_dims(t) == [D]
storage = GT.index_storage(index)
expr = quote
    $(GT.unpack_index_storage(index,:storage))
    $(GT.face_index(index,D)) = $faces[3]
    $(GT.point_index(index)) = 2
    $(GT.topological_sort(expr,())[1])
end
display(expr)
r = eval(expr)
@test r == 10

x1 = GT.point_quantity([SVector{1,Float64}[[0],[1]]],Γ;reference=true)
u1 = GT.physical_map(mesh,D-1)
u2inv = GT.inverse_physical_map(mesh,D)
q = (s2∘u2inv)(u1(x1))
index = GT.generate_index(Γ)
faces = GT.get_symbol!(index,GT.faces(Γ),"faces")
t = GT.term(q,index)
print_tree(t)
expr = GT.expression(t)
@test GT.free_dims(t) == [D-1]
storage = GT.index_storage(index)
expr = quote
    $(GT.unpack_index_storage(index,:storage))
    $(GT.face_index(index,D-1)) = $faces[3]
    $dof = 5
    $(GT.point_index(index)) = 2
    $(GT.topological_sort(expr,())[1])
end
display(expr)
r = eval(expr)
@test r == 5

x1 = GT.point_quantity([SVector{1,Float64}[[0],[1]]],Γ;reference=true)
phi = GT.physical_map(mesh,D-1)
phiD = GT.physical_map(mesh,D)
phiDinv = GT.inverse_physical_map(mesh,D)
q = (phiD∘phiDinv)(phi(x1))
index = GT.generate_index(Γ)
faces = GT.get_symbol!(index,GT.faces(Γ),"faces")
t = GT.term(q,index)
print_tree(t)
expr = GT.expression(t)
@test GT.free_dims(t) == [D-1]
storage = GT.index_storage(index)
expr = quote
    $(GT.unpack_index_storage(index,:storage))
    $(GT.face_index(index,D-1)) = $faces[3]
    $dof = 5
    $(GT.point_index(index)) = 2
    $(GT.topological_sort(expr,())[1])
end
display(expr)
r = eval(expr)
@test [0.5, 0.0] == r

x1 = GT.point_quantity([SVector{1,Float64}[[0],[1]]],Γ;reference=true)
phi = GT.physical_map(mesh,D-1)
n1 = GT.unit_normal(mesh,D-1)
q = n1(phi(x1))
index = GT.generate_index(Γ)
faces = GT.get_symbol!(index,GT.faces(Γ),"faces")
t = GT.term(q,index)
print_tree(t)
expr = GT.expression(t)
@test GT.free_dims(t) == [D-1]
storage = GT.index_storage(index)
expr = quote
    $(GT.unpack_index_storage(index,:storage))
    $(GT.face_index(index,D-1)) = $faces[3]
    $dof = 5
    $(GT.point_index(index)) = 2
    $(GT.topological_sort(expr,())[1])
end
display(expr)
r = eval(expr)
@test r == [0.0, -1.0]

x1 = GT.point_quantity([SVector{1,Float64}[[0],[1]]],Λ;reference=true)
phi = GT.physical_map(mesh,D-1)
n1 = GT.unit_normal(mesh,D-1)
q = n1(phi(x1))[1]
index = GT.generate_index(Λ)
faces = GT.get_symbol!(index,GT.faces(Λ),"faces")
t = GT.term(q,index)
print_tree(t)
expr = GT.expression(t)
@test GT.free_dims(t) == [D-1]
storage = GT.index_storage(index)
expr = quote
    $(GT.unpack_index_storage(index,:storage))
    $(GT.face_index(index,D-1)) = $faces[3]
    $dof = 5
    $(GT.point_index(index)) = 2
    $(GT.topological_sort(expr,())[1])
end
display(expr)
r = eval(expr)
@test r == [0.0, 1.0]

rid_to_dof_to_s = [[x -> dof*sum(x) for dof in 1:5]]
form_arity = 1
axis = 1
field = 1
s2 = GT.form_argument(axis,field,rid_to_dof_to_s,Ω;reference=true)
x2 = GT.point_quantity([SVector{2,Float64}[[0,0],[1,1]]],Ω;reference=true)
q = s2(x2)
index = GT.generate_index(Ω,form_arity)
faces = GT.get_symbol!(index,GT.faces(Ω),"faces")
t = GT.term(q,index)
print_tree(t)
expr = GT.expression(t)
@test GT.free_dims(t) == [D]
storage = GT.index_storage(index)
expr = quote
    $(GT.unpack_index_storage(index,:storage))
    $(GT.face_index(index,D)) = $faces[3]
    $(GT.dof_index(index,axis)) = 5
    $(GT.field_index(index,axis)) = 1
    $(GT.point_index(index)) = 2
    $(GT.topological_sort(expr,())[1])
end
display(expr)
r = eval(expr)
@test r == 10

rid_to_dof_to_s = [[x -> dof*sum(x) for dof in 1:5]]
form_arity = 1
axis = 1
field = 1
s2 = GT.form_argument(axis,field,rid_to_dof_to_s,Ω;reference=true)
x2 = GT.point_quantity([SVector{2,Float64}[[0,0],[1,1]]],Ω;reference=true)
q = s2(x2)
index = GT.generate_index(Ω,form_arity)
faces = GT.get_symbol!(index,GT.faces(Ω),"faces")
t = GT.term(q,index)
print_tree(t)
expr = GT.expression(t)
@test GT.free_dims(t) == [D]
storage = GT.index_storage(index)
expr = quote
    $(GT.unpack_index_storage(index,:storage))
    $(GT.face_index(index,D)) = $faces[3]
    $(GT.dof_index(index,axis)) = 5
    $(GT.field_index(index,axis)) = 2
    $(GT.point_index(index)) = 2
    $(GT.topological_sort(expr,())[1])
end
display(expr)
r = eval(expr)
@test r == 0

rid_to_dof_to_s = [[x -> dof*sum(x) for dof in 1:5]]
form_arity = 1
axis = 1
field = 1
s2 = GT.form_argument(axis,field,rid_to_dof_to_s,Ω;reference=true)
x1 = GT.point_quantity([SVector{1,Float64}[[0],[1]]],Γ;reference=true)
q = s2(x1)
index = GT.generate_index(Γ,form_arity)
faces = GT.get_symbol!(index,GT.faces(Γ),"faces")
t = GT.term(q,index)
print_tree(t)
expr = GT.expression(t)
@test GT.free_dims(t) == [D-1]
storage = GT.index_storage(index)
expr = quote
    $(GT.unpack_index_storage(index,:storage))
    $(GT.face_index(index,D-1)) = $faces[3]
    $(GT.dof_index(index,axis)) = 5
    $(GT.field_index(index,axis)) = 1
    $(GT.point_index(index)) = 2
    $(GT.topological_sort(expr,())[1])
end
display(expr)
r = eval(expr)
@test r == 5

rid_to_dof_to_s = [[x -> dof*sum(x) for dof in 1:5]]
form_arity = 1
axis = 1
field = 1
s2 = GT.form_argument(axis,field,rid_to_dof_to_s,Ω;reference=true)
x1 = GT.point_quantity([SVector{1,Float64}[[0],[1]]],Λ;reference=true)
q = s2(x1)[1]
index = GT.generate_index(Λ,form_arity)
faces = GT.get_symbol!(index,GT.faces(Λ),"faces")
t = GT.term(q,index)
print_tree(t)
expr = GT.expression(t)
@test GT.free_dims(t) == [D-1]
storage = GT.index_storage(index)
expr = quote
    $(GT.unpack_index_storage(index,:storage))
    $(GT.face_index(index,D-1)) = $faces[3]
    $(GT.dof_index(index,axis)) = 5
    $(GT.field_index(index,axis)) = 1
    $(GT.face_around_index(index,axis)) = 1
    $(GT.point_index(index)) = 2
    $(GT.topological_sort(expr,())[1])
end
display(expr)
r = eval(expr)
@test r == 5

rid_to_dof_to_s = [[x -> dof*sum(x) for dof in 1:5]]
form_arity = 1
axis = 1
field = 1
s2 = GT.form_argument(axis,field,rid_to_dof_to_s,Ω;reference=true)
x1 = GT.point_quantity([SVector{1,Float64}[[0],[1]]],Λ;reference=true)
q = s2(x1)[2]
index = GT.generate_index(Λ,form_arity)
faces = GT.get_symbol!(index,GT.faces(Λ),"faces")
t = GT.term(q,index)
print_tree(t)
expr = GT.expression(t)
@test GT.free_dims(t) == [D-1]
storage = GT.index_storage(index)
expr = quote
    $(GT.unpack_index_storage(index,:storage))
    $(GT.face_index(index,D-1)) = $faces[3]
    $(GT.dof_index(index,axis)) = 5
    $(GT.field_index(index,axis)) = 1
    $(GT.face_around_index(index,axis)) = 1
    $(GT.point_index(index)) = 2
    $(GT.topological_sort(expr,())[1])
end
display(expr)
r = eval(expr)
@test r == 0

rid_to_dof_to_s = [[x -> dof*sum(x) for dof in 1:5]]
form_arity = 1
axis = 1
field = 1
s2 = GT.form_argument(axis,field,rid_to_dof_to_s,Ω;reference=true)
phiDinv = GT.inverse_physical_map(mesh,D)
phid = GT.physical_map(mesh,D-1)
r2 = s2 ∘ phiDinv
x1 = GT.point_quantity([SVector{1,Float64}[[0],[1]]],Γ;reference=true)
q = r2(phid(x1))
index = GT.generate_index(Γ,form_arity)
faces = GT.get_symbol!(index,GT.faces(Γ),"faces")
t = GT.term(q,index)
print_tree(t)
expr = GT.expression(t)
@test GT.free_dims(t) == [D-1]
storage = GT.index_storage(index)
expr = quote
    $(GT.unpack_index_storage(index,:storage))
    $(GT.face_index(index,D-1)) = $faces[3]
    $(GT.dof_index(index,axis)) = 5
    $(GT.field_index(index,axis)) = 1
    $(GT.point_index(index)) = 2
    $(GT.topological_sort(expr,())[1])
end
display(expr)
r = eval(expr)
@test r == 5

rid_to_dof_to_s = [[x -> dof*sum(x) for dof in 1:5]]
form_arity = 1
axis = 1
field = 1
s2 = GT.form_argument(axis,field,rid_to_dof_to_s,Ω;reference=true)
phiDinv = GT.inverse_physical_map(mesh,D)
phid = GT.physical_map(mesh,D-1)
r2 = s2 ∘ phiDinv
x1 = GT.point_quantity([SVector{1,Float64}[[0],[1]]],Λ;reference=true)
q = r2(phid(x1))[1]
index = GT.generate_index(Λ,form_arity)
faces = GT.get_symbol!(index,GT.faces(Λ),"faces")
t = GT.term(q,index)
print_tree(t)
expr = GT.expression(t)
@test GT.free_dims(t) == [D-1]
storage = GT.index_storage(index)
expr = quote
    $(GT.unpack_index_storage(index,:storage))
    $(GT.face_index(index,D-1)) = $faces[3]
    $(GT.dof_index(index,axis)) = 5
    $(GT.field_index(index,axis)) = 1
    $(GT.point_index(index)) = 2
    $(GT.face_around_index(index,axis)) = 1
    $(GT.topological_sort(expr,())[1])
end
display(expr)
r = eval(expr)
@test r == 10

rid_to_dof_to_s = [[x -> dof*sum(x) for dof in 1:5]]
form_arity = 1
axis = 1
field = 1
s2 = GT.form_argument(axis,field,rid_to_dof_to_s,Ω;reference=true)
phiDinv = GT.inverse_physical_map(mesh,D)
phid = GT.physical_map(mesh,D-1)
r2 = s2 ∘ phiDinv
x1 = GT.point_quantity([SVector{1,Float64}[[0],[1]]],Λ;reference=true)
q = r2(phid(x1))[2]
index = GT.generate_index(Λ,form_arity)
faces = GT.get_symbol!(index,GT.faces(Λ),"faces")
t = GT.term(q,index)
print_tree(t)
expr = GT.expression(t)
@test GT.free_dims(t) == [D-1]
storage = GT.index_storage(index)
expr = quote
    $(GT.unpack_index_storage(index,:storage))
    $(GT.face_index(index,D-1)) = $faces[3]
    $(GT.dof_index(index,axis)) = 5
    $(GT.field_index(index,axis)) = 1
    $(GT.point_index(index)) = 2
    $(GT.face_around_index(index,axis)) = 1
    $(GT.topological_sort(expr,())[1])
end
display(expr)
r = eval(expr)
@test r == 0

end #module

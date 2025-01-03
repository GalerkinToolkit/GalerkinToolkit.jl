module DomainTests

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
fs = GT.reference_faces(mesh,D-1)
GT.label_interior_faces!(mesh;physical_name="interior_faces")
GT.label_boundary_faces!(mesh;physical_name="boundary_faces")

Ω = GT.interior(mesh)
Ωref = GT.interior(mesh;is_reference_domain=true)

D = GT.num_dims(mesh)
Γref = GT.boundary(mesh;
                 is_reference_domain=true,
                 physical_names=["boundary_faces"])

Γ = GT.physical_domain(Γref)

Λref = GT.skeleton(mesh;
                 is_reference_domain=true,
                 physical_names=["interior_faces"])

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

#u1 = GT.analytical_field(sum,Ω)
#u2 = GT.face_map(mesh,D)
#x1 = GT.point_quantity([SVector{2,Float64}[[0,0],[1,1]]],Ω;reference=true)
#x2 = GT.point_quantity(fill(SVector{2,Float64}[[4,5],[6,7]],GT.num_faces(Ω)),Ω)
#
#q = u1(x2)
#index = GT.generate_index(Ω)
#faces = GT.get_symbol!(index,GT.faces(Ω),"faces")
#t = GT.term(q,index)
#print_tree(t)
#@test GT.num_dims(t) == D
#expr = GT.expression(t)
#storage = GT.index_storage(index)
#expr = quote
#    $(GT.unpack_index_storage(index,:storage))
#    $(GT.face_index(index,D)) = $faces[3]
#    $(GT.point_index(index)) = 2
#    $(GT.topological_sort(GT.simplify(expr),())[1])
#    #$(GT.topological_sort(t.expr,())[1])
#end
#display(expr)
#r = eval(expr)
#@test r == 13
#
#q = u1(x1)
#index = GT.generate_index(Ω)
#faces = GT.get_symbol!(index,GT.faces(Ω),"faces")
#t = GT.term(q,index)
#print_tree(t)
#@test GT.num_dims(t) == D
#expr = GT.expression(t)
#storage = GT.index_storage(index)
#expr = quote
#    $(GT.unpack_index_storage(index,:storage))
#    $(GT.face_index(index,D)) = $faces[3]
#    $(GT.point_index(index)) = 2
#    $(GT.topological_sort(GT.simplify(expr),())[1])
#    #$(GT.topological_sort(t.expr,())[1])
#end
#display(expr)
#r = eval(expr)
#@test r == [1.0, 1.0]
#
#q = ForwardDiff.gradient(u1,x2)
#index = GT.generate_index(Ω)
#faces = GT.get_symbol!(index,GT.faces(Ω),"faces")
#t = GT.term(q,index)
#print_tree(t)
#@test GT.num_dims(t) == D
#expr = GT.expression(t)
#storage = GT.index_storage(index)
#expr = quote
#    $(GT.unpack_index_storage(index,:storage))
#    $(GT.face_index(index,D)) = $faces[3]
#    $(GT.point_index(index)) = 2
#    $(GT.topological_sort(GT.simplify(expr),())[1])
#    #$(GT.topological_sort(t.expr,())[1])
#end
#display(expr)
#r = eval(expr)
#@test r == [1.0, 1.0]
#
#xxxx
#
#q = u2(x1)
#index = GT.generate_index(Ω)
#faces = GT.get_symbol!(index,GT.faces(Ω),"faces")
#t = GT.term(q,index)
#@test t.dim == 2
#storage = GT.index_storage(index)
#expr = quote
#    $(GT.unpack_index_storage(index,:storage))
#    $(GT.face_index(index,D)) = $faces[3]
#    $(GT.point_index(index)) = 2
#    $(GT.topological_sort(GT.simplify(t.expr),())[1])
#    #$(GT.topological_sort(t.expr,())[1])
#end
#display(expr)
#r = eval(expr)
#@test r == [0.75, 0.25]
#
#q = ForwardDiff.jacobian(u2,x1)
#index = GT.generate_index(Ω)
#faces = GT.get_symbol!(index,GT.faces(Ω),"faces")
#t = GT.term(q,index)
#@test t.dim == 2
#storage = GT.index_storage(index)
#expr = quote
#    $(GT.unpack_index_storage(index,:storage))
#    $(GT.face_index(index,D)) = $faces[3]
#    $(GT.point_index(index)) = 2
#    $(GT.topological_sort(GT.simplify(t.expr),())[1])
#    #$(GT.topological_sort(t.expr,())[1])
#end
#display(expr)
#r = eval(expr)
#@test r == [0.25 0.0; 0.0 0.25]
#
#u3 = GT.inverse_face_map(mesh,D)
#q = u3(x1)
#index = GT.generate_index(Ω)
#faces = GT.get_symbol!(index,GT.faces(Ω),"faces")
#t = GT.term(q,index)
#@test t.dim == 2
#storage = GT.index_storage(index)
#expr = quote
#    $(GT.unpack_index_storage(index,:storage))
#    $(GT.face_index(index,D)) = $faces[3]
#    $(GT.point_index(index)) = 2
#    $(GT.topological_sort(GT.simplify(t.expr),())[1])
#    #$(GT.topological_sort(t.expr,())[1])
#end
#display(expr)
#r = eval(expr)
#@test r == [2.0, 4.0]
#
#u4 = GT.inverse_face_map(mesh,D-1)
#q = u4(x1)
#index = GT.generate_index(Γ)
#faces = GT.get_symbol!(index,GT.faces(Γ),"faces")
#t = GT.term(q,index)
#@test t.dim == 1
#storage = GT.index_storage(index)
#expr = quote
#    $(GT.unpack_index_storage(index,:storage))
#    $(GT.face_index(index,D-1)) = $faces[3]
#    $(GT.point_index(index)) = 2
#    $(GT.topological_sort(GT.simplify(t.expr),())[1])
#    #$(GT.topological_sort(t.expr,())[1])
#end
#display(expr)
#r = eval(expr)
#@test r == [3.0]
#
#u5 = u1∘GT.inverse_face_map(mesh,D)
#q = u5(x1)
#index = GT.generate_index(Ω)
#faces = GT.get_symbol!(index,GT.faces(Ω),"faces")
#t = GT.term(q,index)
#@test t.dim == 2
#storage = GT.index_storage(index)
#expr = quote
#    $(GT.unpack_index_storage(index,:storage))
#    $(GT.face_index(index,D)) = $faces[3]
#    $(GT.point_index(index)) = 2
#    $(GT.topological_sort(GT.simplify(t.expr),())[1])
#    #$(GT.topological_sort(t.expr,())[1])
#end
#display(expr)
#r = eval(expr)
#@test r == 6.0
#
#u = GT.face_map(mesh,D-1,D)
#q = u(x1)
#index = GT.generate_index(Γ)
#faces = GT.get_symbol!(index,GT.faces(Γ),"faces")
#t = GT.term(q,index)
#@test t.dim == 1
#storage = GT.index_storage(index)
#expr = quote
#    $(GT.unpack_index_storage(index,:storage))
#    $(GT.face_index(index,D-1)) = $faces[3]
#    $(GT.point_index(index)) = 2
#    $(GT.topological_sort(GT.simplify(t.expr),())[1])
#    #$(GT.topological_sort(t.expr,())[1])
#end
#display(expr)
#r = eval(expr)
#@test [1.0, 0.0] == r
#
#u = GT.face_map(mesh,D-1,D)
#x3 = GT.point_quantity([SVector{1,Float64}[[0],[1]]],Λ;reference=true)
#q = u[2](x3)
#index = GT.generate_index(Λ)
#faces = GT.get_symbol!(index,GT.faces(Λ),"faces")
#t = GT.term(q,index)
#@test t.dim == 1
#storage = GT.index_storage(index)
#expr = quote
#    $(GT.unpack_index_storage(index,:storage))
#    $(GT.face_index(index,D-1)) = $faces[3]
#    $(GT.point_index(index)) = 2
#    #$(GT.topological_sort(GT.simplify(t.expr),())[1])
#    $(GT.topological_sort(t.expr,())[1])
#end
#display(expr)
#r = eval(expr)
#@test [1.0, 0.0] == r
#
#xxx
#
#
#
#
#
#
#
#
#u2 = GT.face_constant_field(facedata,Ω)
#indices = Dict(Ω=>GT.index(face=:face))
#dict = Dict{Any,Symbol}()
#expr = GT.term(u2)(indices,dict)
#storage = GT.storage(dict)
#expr = quote
#    $(GT.unpack_storage(dict,:storage))
#    face = $nfaces
#    $expr(0)
#end
#display(expr)
#@test eval(expr) == facedata[end]
#indices = Dict(Ω=>GT.index(face=GT.FaceList(:faces)))
#dict = Dict{Any,Symbol}()
#expr = GT.term(u2)(indices,dict)
#storage = GT.storage(dict)
#expr = quote
#    $(GT.unpack_storage(dict,:storage))
#    faces = (4,5,6)
#    $expr[2](0)
#end
#display(expr)
#@test eval(expr) == facedata[5]
#
#xxx
#
#ϕ = GT.domain_map(Ωref,Ω)
#
#ϕinv = GT.inverse_map(ϕ)
#
#uref = u∘ϕ
#u2 = uref∘ϕinv
#
#vtk_grid(joinpath(outdir,"omega"),Ω,;plot_params=(;refinement=4)) do plt
#    GT.plot!(plt,u;label="u")
#    GT.plot!(plt,u2;label="u2")
#end
#
#vtk_grid(joinpath(outdir,"omega_ref"),Ωref;plot_params=(;refinement=4)) do plt
#    GT.plot!(plt,uref;label="u")
#end
#
#D = GT.num_dims(mesh)
#Γref = GT.boundary(mesh;
#                 is_reference_domain=true,
#                 physical_names=["boundary_faces"])
#
#ϕ = GT.domain_map(Γref,Ωref)
#g = uref∘ϕ
#@test GT.domain(g) === GT.domain(ϕ)
#@test GT.domain(g) === Γref
#Γ = GT.physical_domain(Γref)
#
#n = GT.unit_normal(Γref,Ω)
#n2 = GT.unit_normal(Γ,Ω)
##h = GT.face_diameter_field(Γ)
#
#vtk_grid(joinpath(outdir,"gamma_ref"),Γref) do plt
#    GT.plot!(plt,g;label="u")
#    GT.plot!(plt,n;label="n")
#    GT.plot!(plt,q->n2(ϕ(q));label="n2") # TODO
#    # GT.plot!(plt,h;label="h") TODO
#    GT.plot!(plt;label="u2") do q
#        x = ϕ(q)
#        uref(x)
#    end
#end
#
#vtk_grid(joinpath(outdir,"gamma"),Γ) do plt
#    GT.plot!(plt,u;label="u")
#    GT.plot!(plt,n;label="n")
#    GT.plot!(plt,n2;label="n2")
#    # TODO
#    #GT.plot!(plt,h;label="h")
#end
#
#Λref = GT.skeleton(mesh;
#                 is_reference_domain=true,
#                 physical_names=["interior_faces"])
#
#ϕ = GT.domain_map(Λref,Ωref)
#
#Λ = GT.physical_domain(Λref)
#
## TODO
#n = GT.unit_normal(Λref,Ω)
#n2 = GT.unit_normal(Λ,Ω)
##h = GT.face_diameter_field(Λ)
#
#jump(u) = u[+] - u[-]
#
#vtk_grid(joinpath(outdir,"lambda_ref"),Λref) do plt
#    GT.plot!(plt,n[+];label="n1")
#    GT.plot!(plt,n[-];label="n2")
#    GT.plot!(plt;label="jump_u2") do q
#        jump((uref∘ϕ)(q))
#    end
#    #GT.plot!(plt,h;label="h")
#    GT.plot!(plt;label="jump_u") do q
#        jump(uref(ϕ(q)))
#    end
#end
#
#jump2(u) = q -> jump(u(q))
#
#vtk_grid(joinpath(outdir,"lambda"),Λ) do plt
#    GT.plot!(plt,n2[+];label="n1")
#    GT.plot!(plt,n2[-];label="n2")
#    GT.plot!(plt;label="jump_u") do q
#        jump(u(q))
#    end
#    GT.plot!(plt,jump2(u);label="jump_u2")
#    #GT.plot!(plt,h;label="h")
#end
#
## Parallel
#
#domain = (0,1,0,1)
#cells_per_dir = (4,4)
#parts_per_dir = (2,2)
#np = prod(parts_per_dir)
#parts = pa.DebugArray(LinearIndices((np,)))
## TODO make this strategy the default one
#partition_strategy = GT.partition_strategy(graph_nodes=:cells,graph_edges=:nodes,ghost_layers=0)
#mesh = GT.cartesian_mesh(domain,cells_per_dir;parts_per_dir,parts,partition_strategy)
#
#GT.physical_faces(mesh,1)
#GT.node_coordinates(mesh)
#GT.face_nodes(mesh,2)
#GT.face_nodes(mesh)
#GT.periodic_nodes(mesh)
#GT.outwards_normals(mesh)
#
#
## TODO
##GT.label_interior_faces!(mesh;physical_name="interior_faces")
## TODO There is a bug when using a ghost layer
#GT.label_boundary_faces!(mesh;physical_name="boundary_faces")
#
#Ω = GT.interior(mesh)
#Ωref = GT.interior(mesh;is_reference_domain=true)
#
#@test Ω == Ω
#@test Ω != Ωref
#
#u = GT.analytical_field(x->sum(x),Ω)
#ϕ = GT.domain_map(Ωref,Ω)
#uref = u∘ϕ
#
#pa.partition(Ω)
#
#GT.faces(Ω)
#
#vtk_grid(joinpath(outdir,"p_omega"),Ω;plot_params=(;refinement=4)) do plt
#    GT.plot!(plt,u;label="u")
#    GT.plot!(plt,q->u(q);label="u2")
#end
#
#D = GT.num_dims(mesh)
#Γref = GT.boundary(mesh;
#                 is_reference_domain=true,
#                 physical_names=["boundary_faces"])
#
#ϕ = GT.domain_map(Γref,Ωref)
#g = uref∘ϕ
#
#vtk_grid(joinpath(outdir,"p_gamma_ref"),Γref) do plt
#    GT.plot!(plt,g;label="u")
#    GT.plot!(plt;label="u2") do q
#        x = ϕ(q)
#        uref(x)
#    end
#end

end # module

module SymbolicsTests

using Test
using LinearAlgebra
import GalerkinToolkit as GT
import ForwardDiff


# function comp_terms(a, b)
#     if a == b || (isa(a, LineNumberNode) && isa(b, LineNumberNode))
#         return true
#     end
#     if a.head != b.head || size(a.args) != size(b.args)
#         return false
#     end
#     return all(map(comp_terms, a.args, b.args))
# end

function expr_has(expr, a)
    if expr == a
        return true
    end
    if isa(expr, LineNumberNode) || isa(expr, Symbol) || isa(expr, Function) || isa(expr, Number)
        return false
    end
    any(map(x -> expr_has(x, a), expr.args))
end

# manual expr blocks to check whether the rules work 
expr = quote
    call( a ∘ inverse_map_impl(phi,b) , call(phi,c) )
    gradient(a ∘ inverse_map_impl(phi,b), call(phi,c)) 
end
expr = GT.simplify(expr)
@test !expr_has(expr, :inverse_map_impl)


expr = quote
    call(face_function(a,b,c,d,e),reference_value(f,h,i)[g])
    jacobian(face_function(a,b,c,d,e),reference_value(f,h,i)[g]) 
end
expr = GT.simplify(expr)
@test !expr_has(expr, :face_function)


expr = quote
    call(face_shape_function(rid_to_fs,face_to_rid,face,dof), reference_value(rid_to_coords, face_to_rid2, sface)[point])
    gradient(face_shape_function(rid_to_fs,face_to_rid,face,dof), reference_value(rid_to_coords, face_to_rid2, sface)[point])
end
expr = GT.simplify(expr)
@test !expr_has(expr, :face_shape_function)


expr = quote
    gradient(face_shape_function(rid_to_fs,face_to_rid,face,dof) ∘ inverse_map_impl(face_function(a,b,c,d,e), bb), call(face_function(a,b,c,d,e), (reference_value(f, h, i))[point]))
    call( a ∘ inverse_map_impl(face_shape_function(rid_to_fs,face_to_rid,face,dof),b) , call(face_shape_function(rid_to_fs,face_to_rid,face,dof),reference_value(rid_to_coords, face_to_rid2, sface)[point]) )
    gradient( a ∘ inverse_map_impl(face_function(a,b,c,d,e),b) , call(face_function(a,b,c,d,e),reference_value(f,h,i)[g]) )
    :(call(var"##constant_quantity_value#378", var"##constant_quantity_value#327", call(var"##constant_quantity_value#378", call(var"##constant_quantity_value#365", gradient(face_function_free_and_dirichlet(var"##refid_to_funs#330", var"##face_to_refid#329", var"##face_to_dofs#331", var"##free_vals#332", var"##diri_vals#333", (var"##sface_to_tfaces#328"[sface])[1]) ∘ inverse_map_impl(face_function(var"##refid_to_funs#335", var"##face_to_refid#338", var"##face_to_nodes#336", var"##node_to_coords#337", var"##sface_to_face#334"[sface]), var"##constant_quantity_value#339"), call(face_function(var"##refid_to_funs#335", var"##face_to_refid#338", var"##face_to_nodes#336", var"##node_to_coords#337", var"##sface_to_face#334"[sface]), (reference_value(var"##refid_to_coords#342", var"##face_to_refid#338", var"##sface_to_face#334"[sface]))[point])), gradient(face_shape_function(var"##refid_to_funs#350", var"##face_to_refid#349", (var"##sface_to_tfaces#348"[sface])[1], idof) ∘ inverse_map_impl(face_function(var"##refid_to_funs#335", var"##face_to_refid#338", var"##face_to_nodes#336", var"##node_to_coords#337", var"##sface_to_face#334"[sface]), var"##constant_quantity_value#339"), call(face_function(var"##refid_to_funs#335", var"##face_to_refid#338", var"##face_to_nodes#336", var"##node_to_coords#337", var"##sface_to_face#334"[sface]), (reference_value(var"##refid_to_coords#342", var"##face_to_refid#338", var"##sface_to_face#334"[sface]))[point]))), call(var"##constant_quantity_value#378", (reference_value(var"##refid_to_ws#368", var"##face_to_refid#338", var"##sface_to_face#334"[sface]))[point], call(var"##constant_quantity_value#377", jacobian(face_function(var"##refid_to_funs#335", var"##face_to_refid#338", var"##face_to_nodes#336", var"##node_to_coords#337", var"##sface_to_face#334"[sface]), (reference_value(var"##refid_to_coords#376", var"##face_to_refid#338", var"##sface_to_face#334"[sface]))[point]))))))
end
expr = GT.simplify(expr)
@test !expr_has(expr, :face_shape_function) && !expr_has(expr, :inverse_map_impl) && !expr_has(expr, :face_function)



# real equations to check whether the rewrited term gets the same math result
# TODO: check it in `integration_tests.jl`? or we can pass an argument to `sum(::Integral)` and get 2 result with/without `simplify`

domain = (0,1,0,1)
cells = (4,4)
mesh = GT.cartesian_mesh(domain,cells)
GT.label_boundary_faces!(mesh;physical_name="boundary_faces")
GT.label_interior_faces!(mesh;physical_name="interior_faces")

Ω = GT.interior(mesh)
Ωref = GT.interior(mesh;is_reference_domain=true)
ϕ = GT.domain_map(Ωref,Ω)
u = GT.analytical_field(x->prod(x),Ω)

degree = 2
dΩref = GT.measure(Ωref,degree)
dΩ = GT.measure(Ω,degree)

int = GT.∫(dΩ) do q
    u(q)
end

s = sum(int)
@test s ≈ 0.25


end # module

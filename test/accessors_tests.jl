module AccessorsTests

using Test
import GalerkinToolkit as GT
using ForwardDiff
#using InteractiveUtils

domain = (0,1,0,1)
n = 4
cells = (n,n)
mesh = GT.cartesian_mesh(domain,cells)

Ω = GT.interior(mesh)
Γ = GT.boundary(mesh)
Λ = GT.skeleton(mesh)

k = 1
V = GT.lagrange_space(Ω,k)

dΩ = GT.measure(Ω,2*k)
dΓ = GT.measure(Γ,2*k)
dΛ = GT.measure(Λ,2*k)

mesh_dΩ = GT.mesh_accessor(mesh,dΩ)
mesh_dΩ = GT.tabulate(GT.value,mesh_dΩ)
mesh_dΩ = GT.tabulate(ForwardDiff.gradient,mesh_dΩ)
mesh_faces = GT.foreach_face(mesh_dΩ)
mesh_face = mesh_faces[2]
@show GT.nodes(mesh_face)
@show GT.node_coordinates(mesh_face)
mesh_points = GT.foreach_point(mesh_face)
mesh_point = mesh_points[1]
@show GT.coordinate(GT.value,mesh_point)
@show GT.coordinate(ForwardDiff.jacobian,mesh_point)
@show GT.weight(mesh_point)

for mesh_face in GT.foreach_face(mesh_dΩ)
    GT.nodes(mesh_face)
    GT.node_coordinates(mesh_face)
    GT.num_nodes(mesh_face)
    for mesh_point in GT.foreach_point(mesh_face)
        GT.coordinate(mesh_point)
        GT.coordinate(GT.value,mesh_point)
        GT.coordinate(ForwardDiff.jacobian,mesh_point)
        GT.weight(mesh_point)
    end
end

mesh_dΓ = GT.mesh_accessor(mesh,dΓ)
mesh_dΓ = GT.tabulate(GT.value,mesh_dΓ)
mesh_dΓ = GT.tabulate(ForwardDiff.gradient,mesh_dΓ)
mesh_faces = GT.foreach_face(mesh_dΓ)
mesh_face = mesh_faces[2]
@show GT.nodes(mesh_face)
@show GT.node_coordinates(mesh_face)
mesh_points = GT.foreach_point(mesh_face)
mesh_point = mesh_points[1]
@show GT.coordinate(GT.value,mesh_point)
@show GT.coordinate(ForwardDiff.jacobian,mesh_point)
@show GT.weight(mesh_point)

for mesh_face in GT.foreach_face(mesh_dΓ)
    GT.nodes(mesh_face)
    GT.node_coordinates(mesh_face)
    GT.num_nodes(mesh_face)
    for mesh_point in GT.foreach_point(mesh_face)
        GT.coordinate(mesh_point)
        GT.coordinate(GT.value,mesh_point)
        GT.coordinate(ForwardDiff.jacobian,mesh_point)
        GT.weight(mesh_point)
    end
end

mesh_dΛ = GT.mesh_accessor(mesh,dΛ)
mesh_dΛ = GT.tabulate(GT.value,mesh_dΛ)
mesh_dΛ = GT.tabulate(ForwardDiff.gradient,mesh_dΛ)
mesh_dΛ = GT.compute(GT.unit_normal,mesh_dΛ)
mesh_dfaces = GT.foreach_face(mesh_dΛ)
mesh_dface = mesh_dfaces[2]
mesh_Dfaces = GT.foreach_face_around(mesh_dface)
mesh_Dface = mesh_Dfaces[2]
@show GT.nodes(mesh_Dface)
@show GT.node_coordinates(mesh_Dface)
mesh_points = GT.foreach_point(mesh_Dface)
mesh_point = mesh_points[1]
@show GT.coordinate(GT.value,mesh_point)
@show GT.coordinate(ForwardDiff.jacobian,mesh_point)
@show GT.weight(mesh_point)
@show GT.unit_normal(mesh_point)

for mesh_dface in GT.foreach_face(mesh_dΛ)
    for mesh_Dface in GT.foreach_face_around(mesh_dface)
        GT.nodes(mesh_Dface)
        GT.node_coordinates(mesh_Dface)
        GT.num_nodes(mesh_Dface)
        for mesh_point in GT.foreach_point(mesh_Dface)
            GT.coordinate(mesh_point)
            GT.coordinate(GT.value,mesh_point)
            GT.coordinate(ForwardDiff.jacobian,mesh_point)
            GT.weight(mesh_point)
            GT.unit_normal(mesh_point)
        end
    end
end

dΩ_faces = GT.foreach_face(dΩ)
dΩ_face = dΩ_faces[2]
dΩ_points = GT.foreach_point(dΩ_face)
dΩ_point = dΩ_points[2]
@show GT.coordinate(dΩ_point)
@show GT.weight(dΩ_point)

dΓ_faces = GT.foreach_face(dΓ)
dΓ_face = dΓ_faces[2]
dΓ_points = GT.foreach_point(dΓ_face)
dΓ_point = dΓ_points[2]
@show GT.coordinate(dΓ_point)
@show GT.weight(dΓ_point)

dΛ_dfaces = GT.foreach_face(dΛ)
dΛ_dface = dΛ_dfaces[2]
dΛ_points = GT.foreach_point(dΛ_dface)
dΛ_point = dΛ_points[2]
@show GT.coordinate(dΛ_point)
@show GT.weight(dΛ_point)

V_dΩ = GT.space_accessor(V,dΩ)
V_dΩ = GT.tabulate(GT.value,V_dΩ)
V_dΩ = GT.tabulate(ForwardDiff.gradient,V_dΩ)
V_faces = GT.foreach_face(V_dΩ)
V_face = V_faces[1]
V_points = GT.foreach_point(V_face)
V_point = V_points[1]
@show GT.shape_functions(GT.value,V_point)

uh = GT.zero_field(Float64,V)
uh_dΩ = GT.field_accessor(uh,dΩ)
uh_dΩ = GT.tabulate(GT.value,uh_dΩ)
uh_faces = GT.foreach_face(uh_dΩ)
uh_face = uh_faces[1]
uh_points = GT.foreach_point(uh_face)
uh_point = uh_points[1]
@show GT.field(GT.value,uh_point)

f_l_x = GT.node_coordinate_accessor(mesh,2)
@show f_l_x(2)(1)

f_refspace = GT.reference_space_accessor(mesh,2)
refspace = f_refspace(2)
@show refspace
f_nrefnodes = GT.reference_space_accessor(GT.num_nodes,mesh,2)
nrefnodes = f_nrefnodes(2)
@show nrefnodes

topo = GT.topology(mesh)
f_reftopo = GT.reference_topology_accessor(topo,2)
reftopo = f_reftopo(2)
@show reftopo
f_v = GT.reference_topology_accessor(GT.vertex_permutations,topo,1)
v = f_v(2)
@show v

# OK
#@code_warntype GT.shape_function_accessor_reference(GT.value,V,dΩ)

#@code_warntype GT.reference_spaces(V)
#@code_warntype map(GT.num_dofs,GT.reference_spaces(V))
#@code_warntype GT.num_dofs_accessor_interior(V,Ω)

#@code_warntype GT.shape_function_accessor_physical(GT.value,V,dΩ)

Ωface_point_x = GT.coordinate_accessor(dΩ)
Γface_point_x = GT.coordinate_accessor(dΓ)
Λface_point_x = GT.coordinate_accessor(dΛ)

@show Ωface_point_x(2)(4)
@show Γface_point_x(3)(2)
@show Λface_point_x(4)(2)

Ωface_point_w = GT.weight_accessor(dΩ)
Γface_point_w = GT.weight_accessor(dΓ)
Λface_point_w = GT.weight_accessor(dΛ)

@show Ωface_point_w(2)(4)
@show Γface_point_w(3)(2)
@show Λface_point_w(4)(2)

Ωface_point_J = GT.jacobian_accessor(dΩ)
Γface_point_J = GT.jacobian_accessor(dΓ)
Λface_point_J = GT.jacobian_accessor(dΛ)

@show Ωface_point_J(2)(4)
@show Γface_point_J(3)(2)
@show Λface_point_J(4)(2)

Ωface_point_J = GT.jacobian_accessor(dΩ,2)
Γface_point_J = GT.jacobian_accessor(dΓ,2)
Λface_point_J = GT.jacobian_accessor(dΛ,2)

@show Ωface_point_J(2)(4)
@show Γface_point_J(3)(2)
@show Λface_point_J(4,1)(2)

Ωface_point_J = GT.jacobian_accessor(dΩ,Val(2))
Γface_point_J = GT.jacobian_accessor(dΓ,Val(2))
Λface_point_J = GT.jacobian_accessor(dΛ,Val(2))

@show Ωface_point_J(2)(4)
@show Γface_point_J(3)(2)
@show Λface_point_J(4,1)(2)

face_phi = GT.physical_map_accessor(GT.value,mesh,2)
@show face_phi(2)([1,1])

face_phi = GT.physical_map_accessor(GT.value,mesh,1)
@show face_phi(2)([1])

face_phi = GT.physical_map_accessor(ForwardDiff.jacobian,mesh,1)
@show face_phi(2)([1])

face_point_phi = GT.physical_map_accessor(ForwardDiff.jacobian,dΓ,2)
@show face_point_phi(5)(1)

face_point_phi = GT.physical_map_accessor(ForwardDiff.jacobian,dΓ,1)
@show face_point_phi(5)(1)

face_npoints = GT.num_points_accessor(dΩ)
@show face_npoints(3)

face_dof_s = GT.shape_function_accessor(GT.value,V)
@show face_dof_s(3)(1)([1,1])

J = Ωface_point_J(3)(1)
face_point_dof_s = GT.shape_function_accessor(ForwardDiff.gradient,V,dΩ)
@show face_point_dof_s(3)(1,J)(4)

J = Γface_point_J(5)(1)
face_point_dof_s = GT.shape_function_accessor(ForwardDiff.gradient,V,dΓ)
@show face_point_dof_s(5)(1,J)(4)

face_point_dof_s = GT.shape_function_accessor(ForwardDiff.gradient,V,dΛ)
J = Λface_point_J(5,1)(1)
@show face_point_dof_s(5,1)(1,J)(4)
J = Λface_point_J(5,2)(1)
@show face_point_dof_s(5,2)(1,J)(4)

face_point_dof_s = GT.shape_function_accessor(GT.value,V,dΩ)
@show face_point_dof_s(3)(1)(4)

face_point_dof_s = GT.shape_function_accessor(GT.value,V,dΓ)
@show face_point_dof_s(5)(1)(4)

face_point_dof_s = GT.shape_function_accessor(GT.value,V,dΛ)
@show face_point_dof_s(5,1)(1)(4)
@show face_point_dof_s(5,2)(1)(4)

face_point_dof_s = GT.shape_function_accessor(ForwardDiff.gradient,V,dΩ)
@show face_point_dof_s(3)(1)(4)

face_point_dof_s = GT.shape_function_accessor(ForwardDiff.gradient,V,dΓ)
@show face_point_dof_s(5)(1)(4)

face_point_dof_s = GT.shape_function_accessor(ForwardDiff.gradient,V,dΛ)
@show face_point_dof_s(5,1)(1)(4)
@show face_point_dof_s(5,2)(1)(4)

#face_dof_s = GT.form_argument_accessor(GT.value,V)
#@show face_dof_s(3)(1)([1,1])
#
#face_point_dof_s = GT.form_argument_accessor(GT.value,V,dΩ)
#@show face_point_dof_s(3)(1)(4)
#
#face_point_dof_s = GT.form_argument_accessor(GT.value,V,dΓ)
#@show face_point_dof_s(5)(1)(4)
#
#face_point_dof_s = GT.form_argument_accessor(GT.value,V,dΛ)
#@show face_point_dof_s(5,1)(1)(4)
#@show face_point_dof_s(5,2)(1)(4)
#@show face_point_dof_s(5,1)(1)(4,1)
#@show face_point_dof_s(5,2)(1)(4,1)
#@show face_point_dof_s(5,1)(1)(4,2)
#@show face_point_dof_s(5,2)(1)(4,2)
#
#face_point_dof_s = GT.form_argument_accessor(ForwardDiff.gradient,V,dΩ)
#@show face_point_dof_s(3)(1)(4)
#
#face_point_dof_s = GT.form_argument_accessor(ForwardDiff.gradient,V,dΓ)
#@show face_point_dof_s(5)(1)(4)
#
#face_point_dof_s = GT.form_argument_accessor(ForwardDiff.gradient,V,dΛ)
#@show face_point_dof_s(5,1)(1)(4,1)
#@show face_point_dof_s(5,2)(1)(4,1)
#@show face_point_dof_s(5,1)(1)(4,2)
#@show face_point_dof_s(5,2)(1)(4,2)

uh = GT.rand_field(Float64,V)

J = Ωface_point_J(3)(1)
face_point_dof_s = GT.discrete_field_accessor(ForwardDiff.gradient,uh,dΩ)
@show face_point_dof_s(3)(1,J)

J = Γface_point_J(5)(1)
face_point_dof_s = GT.discrete_field_accessor(ForwardDiff.gradient,uh,dΓ)
@show face_point_dof_s(5)(1,J)

face_point_dof_s = GT.discrete_field_accessor(ForwardDiff.gradient,uh,dΛ)
J = Λface_point_J(5,1)(1)
@show face_point_dof_s(5,1)(1,J)
J = Λface_point_J(5,2)(1)
@show face_point_dof_s(5,2)(1,J)

face_point_dof_s = GT.discrete_field_accessor(GT.value,uh,dΩ)
@show face_point_dof_s(3)(1)

face_point_dof_s = GT.discrete_field_accessor(GT.value,uh,dΓ)
@show face_point_dof_s(5)(1)

face_point_dof_s = GT.discrete_field_accessor(GT.value,uh,dΛ)
@show face_point_dof_s(5,1)(1)
@show face_point_dof_s(5,2)(1)

face_point_dof_s = GT.discrete_field_accessor(ForwardDiff.gradient,uh,dΩ)
@show face_point_dof_s(3)(1)

face_point_dof_s = GT.discrete_field_accessor(ForwardDiff.gradient,uh,dΓ)
@show face_point_dof_s(5)(1)

face_point_dof_s = GT.discrete_field_accessor(ForwardDiff.gradient,uh,dΛ)
@show face_point_dof_s(5,1)(1)
@show face_point_dof_s(5,2)(1)

face_point_dof_s = GT.discrete_field_accessor(GT.value,uh,dΩ)
face_point_dof_s = GT.update(face_point_dof_s,discrete_field=uh)
@show face_point_dof_s(3)(1)

#face_diri! =  GT.dirichlet_accessor(uh,Ω)
#face_diri! =  GT.update(face_diri!,discrete_field=uh)
#A = ones(4,4)
#b = zeros(4)
#face_diri!(1)(A,b)
#@show b

face_point_n = GT.unit_normal_accessor(dΓ)
@show face_point_n(6)(1)

face_point_n = GT.unit_normal_accessor(dΛ)
@show face_point_n(6,1)(1)
@show face_point_n(6,2)(1)

diams = GT.face_diameter(Γ)
@show diams

diams = GT.face_diameter(Λ)
@show diams

end # module


function physical_map_accessor(f,mesh,D)
end

function shape_function_accessor_reference(f,space::AbstractSpace)
end

function shape_function_accessor_modifier(f::ForwardDiff.value,space::AbstractSpace)
    function face_dof_modifier(face)
        function dof_modifier(dof)
            function modifier(v,J)
                v
            end
        end
    end
end

function shape_function_accessor_modifier(f::ForwardDiff.gradient,space::AbstractSpace)
    function face_dof_modifier(face)
        function dof_modifier(dof)
            function modifier(v,J)
                transpose(J)\v
            end
        end
    end
end

function shape_function_accessor_modifier(f::ForwardDiff.jacobian,space::AbstractSpace)
    function face_dof_modifier(face)
        function dof_modifier(dof)
            function modifier(v,J)
                v/J
            end
        end
    end
end

function shape_function_accessor_physical(f,space::AbstractSpace)
    domain = GT.domain(space)
    D = num_dims(domain)
    mesh = GT.mesh(domain)
    face_dof_s_ref = shape_function_accessor_reference(value,space)
    face_dof_modif = shape_function_accessor_modifier(value,space)
    face_phi = physical_map_accessor(mesh,Val(D))
    face_Dphi = physical_map_accessor(ForwardDiff.jacobian,mesh,Val(D))
    x0 = zero(SVector{D,Float64})
    function face_dof_s_phys(face)
        phi = face_phi(face)
        Dphi = face_Dphi(face)
        invphi = inv_map(phi,x0)
        dof_modif = face_dof_modif(face)
        function dof_s_phys(dof)
            modif = dof_modif(dof)
            s_ref = dof_s_ref(dof)
            function s_phys(x)
                v = s_ref(x)
                J = Dphi(x)
                modif(v,J)
            end
            x->f(s_phys,x)
        end
    end
end

function shape_function_accessor_reference_interior(f,space::AbstractSpace,measure::AbstractQuadrature)
end

function shape_function_accessor_reference_boundary(f,space::AbstractSpace,measure::AbstractQuadrature)
end

function shape_function_accessor_reference_skeleton(f,space::AbstractSpace,measure::AbstractQuadrature)
end

function shape_function_accessor_physical_interior(f,space::AbstractSpace,measure::AbstractQuadrature)
    Dface_dof_x_s = shape_function_accessor_physical(f,space)
    face_point_x = coordinate_accessor(measure)
    face_Dface = faces(domain(measure))
    function face_point_dof_s(face)
        point_x = face_point_x(face)
        Dface = face_to_Dface[face]
        dof_x_s = Dface_dof_x_s(Dface)
        function point_dof_s(point,J=nothing)
            x = point_x(point)
            function dof_s(dof)
                dof_x_s(dof)(x)
            end
        end
    end
end

function shape_function_accessor_physical_boundary(f,space::AbstractSpace,measure::AbstractQuadrature)
end

function shape_function_accessor_physical_skeleton(f,space::AbstractSpace,measure::AbstractQuadrature)
end

for T in (value,ForwardDiff.gradient,ForwardDiff.jacobian)
    @eval begin

        function shape_function_accessor_physical_interior(f::$T,space::AbstractSpace,measure::AbstractQuadrature)
            face_point_dof_v = shape_function_accessor_reference(f,space,measure)
            face_dof_modif = shape_function_accessor_modifier(f,space)
            face_to_Dface = faces(domain(measure))
            Dface_to_rid = face_reference_id(space)
            J = 
            v = GT.prototype(face_point_dof_v)
            modif = GT.prototype(face_dof_modif)
            prototype = modif(v,J)
            rid_to_dof_s = map(reference_spaces(space)) do reffe
                zeros(T,num_dofs(reffe))
            end
            function face_point_dof_s(face)
                point_dof_v = face_point_dof_v(face)
                dof_modif = face_dof_modif(face)
                Dface = face_to_Dface[face]
                rid = Dface_to_rid[Dface]
                dof_s = rid_to_dof_s[rid]
                function point_dof_s(point,J=nothing)
                    dof_v = point_dof_v(point)
                    ndofs = length(dof_s)
                    for dof in 1:ndofs
                        v = dof_v(dof)
                        modif = dof_modif(dof)
                        dof_s[dof] = modif(v,J)
                    end
                    function dof_f(dof)
                        dof_s[dof]
                    end
                end
            end
            accessor()
        end

        function shape_function_accessor_physical_boundary(f::$T,space::AbstractSpace,measure::AbstractQuadrature)
        end

        function shape_function_accessor_physical_skeleton(f::$T,space::AbstractSpace,measure::AbstractQuadrature)
        end

    end
end



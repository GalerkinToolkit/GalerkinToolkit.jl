

# Constrained Space

function constrain(a::AbstractSpace,constraints)
    constraints_accessor = setup_constraints_accessor(a,constraints)
    ConstrainedSpace(a,constraints,constraints_accessor)
end

struct ConstrainedSpace{A,B,C} <: AbstractSpace
    parent::A
    constraints::B
    constraints_accessor::C
end

function face_dofs(a::ConstrainedSpace)
    face_dofs(a.constraints_accessor)
end

function face_constraints(a::ConstrainedSpace)
    each_face(a.constraints_accessor)
end

function face_constraints(a::AbstractSpace)
    each_face(identity_constraints_accessor(face_dofs(a)))
end

function constraints(a::ConstrainedSpace)
    a.constraints
end

function constraints(a::AbstractSpace)
    identity_constraints(dofs(a),dofs(a))
end

function interpolate_impl!(
    f::AbstractField,
    u::DiscreteField,
    space::ConstrainedSpace,
    free_or_diri::FreeOrDirichlet)
    parent_space = u.parent

    xf = free_values(u)
    xd = dirichlet_values(u)
    xp = u.workspace
    u_parent = discrete_field(parent_space,xp,xd)
    interpolate_impl!(f,u_parent,parent_space,free_or_diri)
    C = constrains(space)
    restrict_values!(xf,C,xp)
    u
end

function setup_discrete_field_workspace(space::ConstrainedSpace,fv,dv)
    parent_space = space.parent
    pv = similar(fv,free_dofs(parent_space))
    pv
end

# Types of constraints

# Identity (do nothing)

struct IdentityConstraints{A,B} <: AbstractType
    rows::A
    cols::B
end

function prolongate_values!(b,C::IdentityConstraints,a)
    if b !== a
        copyto!(b,a)
    end
    b
end

function restrict_values!(b,C::IdentityConstraints,a)
    if b !== a
        copyto!(b,a)
    end
    b
end

function identity_constraints_accessor(face_dofs)
    face = nothing
    IdentityConstraintsAccessor(face_dofs,face)
end

struct IdentityConstraintsAccessor{A,B} <: NewAbstractAccessor
    face_dofs::A
    face::B
end

num_faces(a::IdentityConstraintsAccessor) = length(a.face_dofs)

function at_face(a::IdentityConstraintsAccessor,face)
    IdentityConstraintsAccessor(a.face_dofs,face)
end

function num_dofs(a::IdentityConstraintsAccessor)
    face = a.face
    dofs = a.face_dofs[face]
    length(dofs)
end

function num_parent_dofs(a::IdentityConstraintsAccessor)
    num_dofs(a)
end

function constrain_shape_functions!(b,C::IdentityConstraintsAccessor,a)
    if b !== a
        copyto!(b,a)
    end
    b
end

# Periodic

struct PeriodicConstraints{A,B,C,D} <: AbstractType
    parent_dof_free_or_periodic_dof::A
    free_dof_parent_dof::B
    periodic_dof_parent_dof::C
    periodic_dof_scaling::D
end

#function prolongate_values!(parent_dof_value,C::PeriodicConstraints,free_dof_value)
#end

function restrict_values!(free_dof_value,C::PeriodicConstraints,parent_dof_value)
    free_dof_parent_dof = C.free_dof_parent_dof
    for free_dof in 1:length(free_dof_parent_dof)
        parent_dof = free_dof_parent_dof[free_dof]
        free_dof_value[free_dof] = parent_dof_value[parent_dof]
    end
    free_dof_value
end

function setup_constraints_accessor(parent_space,C::PeriodicConstraints)
    face_parent_dofs = GT.face_dofs(parent_space)
    parent_dof_free_or_periodic_dof = C.parent_dof_free_or_periodic_dof
    periodic_dof_parent_dof = C.periodic_dof_parent_dof
    face_parent_dofs = jagged_array(face_parent_dofs)
    function dof_map(parent_dof)
        # Do nothing for Dirichlet
        if parent_dof < 0
            return parent_dof
        end
        free_or_periodic_dof = parent_dof_free_or_periodic_dof[parent_dof]
        if free_or_periodic_dof < 0
            periodic_dof = - free_or_periodic_dof
            parent_dof = periodic_dof_parent_dof[periodic_dof]
        else
            parent_dof = free_or_periodic_dof
        end
        return parent_dof
    end
    face_dofs_data = map(dof_map,face_parent_dofs.data)
    face_dofs = jagged_array(face_dofs_data,face_parent_dofs.ptrs)
    if C.periodic_dof_scaling === nothing
        identity_constraints_accessor(face_dofs)
    else
        PeriodicConstraintsAccessor(C,face_parent_dofs,face_dofs)
    end
end

struct PeriodicConstraintsAccessor{A,B,C,D} <: NewAbstractAccessor
    constraints::A
    face_dofs::B
    face_parent_dofs::C
    face::D
end

function num_faces(a::PeriodicConstraintsAccessor)
    length(a.face_dofs)
end

function at_face(a::PeriodicConstraintsAccessor,face)
    PeriodicConstraintsAccessor(
                                a.constraints,
                                a.face_dofs,
                                a.face_parent_dofs,
                                face
                               )
end

function face_dofs(C::PeriodicConstraintsAccessor)
    C.face_dofs
end

function num_dofs(W::PeriodicConstraintsAccessor)
    face = W.face
    dofs = W.face_dofs[face]
    length(dofs)
end

function num_parent_dofs(W::PeriodicConstraintsAccessor)
    face = W.face
    parent_dofs = W.face_parent_dofs[face]
    length(parent_dofs)
end

#function prolongate_values!(b,W::PeriodicConstraintsAccessor,a)
#end

# This is like multiplying with the transpose of W
function constrain_shape_functions!(ldof_s,W::PeriodicConstraintsAccessor,parent_ldof_s)
    C = W.C
    periodic_dof_scaling = C.periodic_dof_scaling
    face = W.face
    parent_ldof_parent_dof = W.face_parent_dofs[face]
    parent_dof_free_or_periodic_dof = C.parent_dof_free_or_periodic_dof
    n_parent_ldofs = length(parent_ldof_parent_dof)
    for parent_ldof in 1:n_parent_ldofs
        parent_dof = parent_ldof_parent_dof[parent_ldof]
        ldof = parent_ldof
        free_or_periodic_dof = parent_dof_free_or_periodic_dof[parent_dof]
        if free_dof_parent_dof < 0
            periodic_dof = - free_or_periodic_dof
            scaling = periodic_dof_scaling[periodic_dof]
            ldof_s[ldof] = scaling*parent_ldof_s[parent_ldof]
        else
            ldof_s[ldof] = parent_ldof_s[parent_ldof]
        end
    end
    ldof_s
end


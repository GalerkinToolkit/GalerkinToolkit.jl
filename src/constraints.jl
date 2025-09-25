


function identity_constraint()
    IdentityConstraint()
end

struct IdentityConstraint <: AbstractType end

function constraints(a::AbstractSpace)
    identity_constraint()
end

function face_constraints(a::AbstractSpace)
    nfaces = length(face_dofs(a))
    Fill(identity_constraint(),nfaces)
end

function constrain_shape_functions!(C::IdentityConstraint,ldof_s,parent_ldof_s,face=nothing)
    parent_ldof_s
end

# Constrained Space

function constrain(a::AbstractSpace,constraints)
    b = ConstrainedSpace(a,constraints)
    workspace = setup_workspace(b)
    ConstrainedSpace(a,constraints,workspace)
end

struct ConstrainedSpace{A,B,W} <: AbstractSpace
    parent::A
    constraints::B
    workspace::C
end

function face_dofs(a::ConstrainedSpace)
    a.workspace.face_dofs
end

function face_constraints(a::ConstrainedSpace)
    FaceConstraints(a.workspace,nothing)
end

function constraints(a::ConstrainedSpace)
    a.constraints
end

struct FaceConstraints{A,B} <: AbstractType
    workspace::A
    face::B
end

function Base.getindex(a::FaceConstraints,face)
    FaceConstraints(a.workspace,face)
end

function constrain_shape_functions!(C::FaceConstraints,ldof_s,parent_ldof_s)
    constrain_shape_functions!(C.workspace,ldof_s,parent_ldof_s,C.face)
end

function setup_workspace(a::ConstrainedSpace)
    constrained_space_workspace(a.parent,a.constraints)
end

struct PeriodicConstraints{A,B,C,D} <: AbstractType
    parent_dof_free_or_periodic_dof::A
    free_dof_parent_dof::B
    periodic_dof_parent_dof::C
    periodic_dof_scaling::D
end

struct PeriodicConstraintsWorkspace{A,B,C} <: AbstractType
    constraints::A
    face_dofs::B
    face_parent_dofs::C
end

function restrict_parent_values!(C::PeriodicConstraints,free_dof_value,parent_dof_value)
    free_dof_parent_dof = C.free_dof_parent_dof
    for free_dof in 1:length(free_dof_parent_dof)
        parent_dof = free_dof_parent_dof[free_dof]
        free_dof_value[free_dof] = parent_dof_value[parent_dof]
    end
    free_dof_value
end

function constrained_space_workspace(parent_space,C::PeriodicConstraints)
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
    PeriodicConstraintsWorkspace(C,face_parent_dofs,face_dofs)
end

function face_dofs(C::PeriodicConstraintsWorkspace)
    C.face_dofs
end

function constrain_shape_functions!(W::PeriodicConstraintsWorkspace,ldof_s,parent_ldof_s,face)
    C = W.C
    periodic_dof_scaling = C.periodic_dof_scaling
    if periodic_dof_scaling === nothing
        return parent_ldof_s
    end
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
    view(ldof_s,1:n_parent_ldofs)
end


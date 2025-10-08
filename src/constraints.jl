

# Constrained Space

function constrain(a::AbstractSpace,constraints)
    workspace = setup_constraints_workspace(a,constraints)
    ConstrainedSpace(a,constraints,workspace)
end

struct ConstrainedSpace{A,B,C} <: AbstractSpace
    parent::A
    constraints::B
    workspace::C
end

function face_dofs(a::ConstrainedSpace)
    face_dofs(a.workspace)
end

function face_constraints(a::ConstrainedSpace)
    each_face(FaceConstraintsAccessor(a.constraints,a.workspace))
end

struct FaceConstraintsAccessor{A,B} <: NewAbstractAccessor
    constraints::A
    workspace::B
end

function num_faces(a::FaceConstraintsAccessor)
    num_faces(a.workspace)
end

function at_face(a::FaceConstraintsAccessor,face)
    constraints_at_face(a.constraints,a.workspace,face)
end

function face_constraints(a::AbstractSpace)
    C = constraints(a)
    workspace = setup_constraints_workspace(a,C)
    each_face(FaceConstraintsAccessor(C,workspace))
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
    free_values!(xf,C,xp)
    u
end

function setup_discrete_field_workspace(space::ConstrainedSpace,fv,dv)
    parent_space = space.parent
    pv = similar(fv,free_dofs(parent_space))
    pv
end

# Types of constraints

# Identity (do nothing)

function identity_constraints(rows,cols)
    IdentityConstraints(Int,rows,cols)
end

struct IdentityConstraints{T,A,B} <: AbstractMatrix{T}
    eltype::Type{T}
    rows::A
    cols::B
end

Base.size(a::IdentityConstraints) = map(length,axes(a))
Base.axes(a::IdentityConstraints) = (a.rows,a.cols)

function LinearAlgebra.mul!(b,C::IdentityConstraints,a)
    if b !== a
        copyto!(b,a)
    end
    b
end

function LinearAlgebra.mul!(
    b,
    C::LinearAlgebra.Transpose{T,<:IdentityConstraints} where T,
    a)
    mul!(b,C.parent,a)
end

function free_values!(b,C::IdentityConstraints,a)
    mul!(b,C,a)
end

function setup_constraints_workspace(a,constraints::IdentityConstraints)
    IdentityConstraintsWorkspace(face_dofs(a))
end

struct IdentityConstraintsWorkspace{A} <: AbstractType
    face_dofs::A
end

function face_dofs(a::IdentityConstraintsWorkspace)
    a.face_dofs
end

function num_faces(a::IdentityConstraintsWorkspace)
    length(a.face_dofs)
end

function constraints_at_face(constraints::IdentityConstraints,workspace,face)
    face_dofs = workspace
    dofs = workspace.face_dofs[face]
    n = length(dofs)
    rows = Base.OneTo(n)
    cols = rows
    identity_constraints(rows,cols)
end

# Periodic

function create_periodic_constraints(
        parent_dof_free_or_periodic_dof,
        free_dof_parent_dof,
        periodic_dof_parent_dof,
        periodic_dof_scaling,
    )

    if periodic_dof_scaling === nothing
        T = Int
    else
        T = eltype(periodic_dof_scaling)
    end
    PeriodicConstraints(
                        T,
                        parent_dof_free_or_periodic_dof,
                        free_dof_parent_dof,
                        periodic_dof_parent_dof,
                        periodic_dof_scaling,
                       )
end

struct PeriodicConstraints{T,A,B,C,D} <: AbstractMatrix{T}
    eltype::Type{T}
    parent_dof_free_or_periodic_dof::A
    free_dof_parent_dof::B
    periodic_dof_parent_dof::C
    periodic_dof_scaling::D
end

function Base.size(a::PeriodicConstraints)
    nparent = length(a.parent_dof_free_or_periodic_dof)
    nfree = length(a.free_dof_parent_dof)
    (nparent,nfree)
end

function LinearAlgebra.mul!(parent_dof_value,C::PeriodicConstraints,free_dof_value)
    error("not implemented")
end

function free_values!(free_dof_value,C::PeriodicConstraints,parent_dof_value)
    free_dof_parent_dof = C.free_dof_parent_dof
    for free_dof in 1:length(free_dof_parent_dof)
        parent_dof = free_dof_parent_dof[free_dof]
        free_dof_value[free_dof] = parent_dof_value[parent_dof]
    end
    free_dof_value
end

function setup_constraints_workspace(parent_space,C::PeriodicConstraints)
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
    face = nothing
    T = eltype(C)
    PeriodicConstraintsWorkspace(face_parent_dofs,face_dofs)
end

struct PeriodicConstraintsWorkspace{A,B} <: AbstractType
    face_dofs::A
    face_parent_dofs::B
end

function face_dofs(C::PeriodicConstraintsWorkspace)
    C.face_dofs
end

function num_faces(C::PeriodicConstraintsWorkspace)
    length(C.face_dofs)
end

function constraints_at_face(C::PeriodicConstraints,workspace,face)
    PeriodicConstraintsAtFace(eltype(C),C,workspace,face)
end

struct PeriodicConstraintsAtFace{T,A,B} <: AbstractMatrix{T}
    eltype::Type{T}
    constraints::A
    workspace::B
    face::Int
end

function Base.size(C::PeriodicConstraintsAtFace)
    face = W.face
    dofs = W.face_dofs[face]
    parent_dofs = W.face_parent_dofs[face]
    (length(parent_dofs),length(dofs))
end

function LinearAlgebra.mul!(
    ldof_s,
    at::LinearAlgebra.Transpose{T,<:PeriodicConstraintsAtFace} where T,
    parent_ldof_s)
    a = at.parent
    C = a.constraints
    W = a.workspace
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


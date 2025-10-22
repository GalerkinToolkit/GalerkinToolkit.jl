

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

workspace(space::ConstrainedSpace) = space.workspace

domain(space::ConstrainedSpace) = domain(space.parent)
order(space::ConstrainedSpace) = order(space.parent)
reference_spaces(space::ConstrainedSpace) = reference_spaces(space.parent)
face_reference_id(space::ConstrainedSpace) = face_reference_id(space.parent)
continuous(space::ConstrainedSpace) = continuous(space.parent)
dirichlet_boundary(space::ConstrainedSpace) = dirichlet_boundary(space.parent)

free_dofs(space::ConstrainedSpace) = axes(space.constraints,2)
dirichlet_dofs(space::ConstrainedSpace) = dirichlet_dofs(space.parent)

function max_num_face_dofs(space)
    max_num_reference_dofs(space)
end

function max_num_face_dofs(space::ConstrainedSpace)
    maximum(length,face_dofs(space))
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
    identity_constraints(dofs(a,FREE),dofs(a,FREE))
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

function Base.getindex(a::IdentityConstraints,i::Integer,j::Integer)
    T = a.eltype
    T(i==j)
end

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

function periodic_constraints(parent_space::AbstractSpace;periodic_scaling=nothing)
    periodic_constraints_impl(parent_space,periodic_scaling)
end

function periodic_constraints_impl(parent_space,periodic_scaling::Nothing)
    periodic_constraints_impl(parent_space)
end

function periodic_constraints_impl(parent_space,periodic_scaling::DiscreteField)
    C = periodic_constraints_impl(parent_space)
    parent_dof_value = GT.free_values(periodic_scaling)
    periodic_dof_scaling = parent_dof_value[C.periodic_dof_parent_dof]
    T = eltype(parent_dof_value)
    PeriodicConstraints(
                        T,
                        C.parent_dof_free_or_periodic_dof,
                        C.free_dof_parent_dof,
                        C.periodic_dof_parent_dof,
                        periodic_dof_scaling,
                       )
end

function periodic_constraints_impl(parent_space,periodic_scaling::AbstractField)
    uh = interpolate(periodic_scaling,parent_space)
    periodic_constraints_impl(parent_space,uh)
end

function periodic_constraints_impl(parent_space,periodic_scaling::Number)
    C = periodic_constraints_impl(parent_space)
    periodic_dof_scaling = Fill(periodic_scaling,length(C.periodic_dof_parent_dof))
    T = eltype(periodic_dof_scaling)
    PeriodicConstraints(
                        T,
                        C.parent_dof_free_or_periodic_dof,
                        C.free_dof_parent_dof,
                        C.periodic_dof_parent_dof,
                        periodic_dof_scaling,
                       )
end

function periodic_constraints_impl(parent_space)
    periodic_dofs = GT.periodic_dofs(parent_space)
    msg = "Periodic dof points to Dirichlet dof. Do not define as Dirichlet the master boundary of a periodic coupling."
    @boundscheck @assert all(dof->dof>0,periodic_dofs) msg
    n_parent_dofs = GT.num_free_dofs(parent_space)
    free_and_periodic_dofs = GT.partition_from_mask(i->periodic_dofs[i] == i,1:n_parent_dofs)
    free_dof_parent_dof = first(free_and_periodic_dofs)
    periodic_dof_parent_dof = periodic_dofs[last(free_and_periodic_dofs)]
    msg0 = "Cyclic periodic constraints not allowed"
    @boundscheck @assert all(i->periodic_dofs[i]==i,periodic_dof_parent_dof) msg0
    dof_permutation = GT.permutation(free_and_periodic_dofs)
    n_free_dofs = length(free_dof_parent_dof)
    f = dof -> begin
        dof2 = dof_permutation[dof]
        T = typeof(dof2)
        if dof2 > n_free_dofs
            return T(n_free_dofs-dof2)
        end
        dof2
    end
    parent_dof_free_or_periodic_dof = f.(1:n_parent_dofs)
    T = Int
    periodic_dof_scaling = nothing
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

function Base.getindex(C::PeriodicConstraints,parent_dof::Integer,j::Integer)
    T = C.eltype
    parent_dof_free_or_periodic_dof = C.parent_dof_free_or_periodic_dof
    free_dof_parent_dof = C.free_dof_parent_dof
    periodic_dof_parent_dof = C.periodic_dof_parent_dof
    periodic_dof_scaling = C.periodic_dof_scaling
    free_or_periodic_dof = parent_dof_free_or_periodic_dof[parent_dof]
    if free_or_periodic_dof > 0
        free_dof = free_or_periodic_dof
        parent_dof_owner = free_dof_parent_dof[free_dof]
        periodic_scaling = T(1)
    else
        periodic_dof = -free_or_periodic_dof
        parent_dof_owner = periodic_dof_parent_dof[periodic_dof]
        if periodic_dof_scaling !== nothing
            periodic_scaling = periodic_dof_scaling[periodic_dof]
        else
            periodic_scaling = T(1)
        end
    end
    free_dof_owner = parent_dof_free_or_periodic_dof[parent_dof_owner]
    if free_dof_owner == j
        periodic_scaling
    else
        zero(periodic_scaling)
    end
end

#function LinearAlgebra.mul!(parent_dof_value,C::PeriodicConstraints,free_dof_value)
#    error("not implemented")
#end

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
        T = eltype(parent_dof_free_or_periodic_dof)
        if free_or_periodic_dof < 0
            periodic_dof = - free_or_periodic_dof
            parent_dof_owner = periodic_dof_parent_dof[periodic_dof]
            free_dof = parent_dof_free_or_periodic_dof[parent_dof_owner]
        else
            free_dof = T(free_or_periodic_dof)
        end
        return free_dof
    end
    face_dofs_data = map(dof_map,face_parent_dofs.data)
    face_dofs = jagged_array(face_dofs_data,face_parent_dofs.ptrs)
    face = nothing
    T = eltype(C)
    PeriodicConstraintsWorkspace(face_dofs,face_parent_dofs)
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
    if C.periodic_dof_scaling !== nothing
        PeriodicConstraintsAtFace(eltype(C),C,workspace,face)
    else
        dofs = workspace.face_dofs[face]
        n = length(dofs)
        rows = Base.OneTo(n)
        cols = rows
        identity_constraints(rows,cols)
    end
end

struct PeriodicConstraintsAtFace{T,A,B} <: AbstractMatrix{T}
    eltype::Type{T}
    constraints::A
    workspace::B
    face::Int
end

function Base.size(C::PeriodicConstraintsAtFace)
    face = C.face
    W = C.workspace
    dofs = W.face_dofs[face]
    parent_dofs = W.face_parent_dofs[face]
    (length(parent_dofs),length(dofs))
end

function Base.getindex(C::PeriodicConstraintsAtFace,parent_ldof::Integer,ldof::Integer)
    T = C.eltype
    W = C.workspace
    face = C.face
    dofs = W.face_dofs[face]
    dof = dofs[ldof]
    parent_dofs = W.face_parent_dofs[face]
    parent_dof = parent_dofs[parent_ldof]
    if dof < 0 && parent_dof < 0 && dof == parent_dof
        return T(1)
    end
    if dof < 0
        return T(0)
    end
    if parent_dof < 0
        return T(0)
    end
    C.constraints[parent_dof,dof]
end

function LinearAlgebra.mul!(
    ldof_s,
    at::LinearAlgebra.Transpose{T,<:PeriodicConstraintsAtFace} where T,
    parent_ldof_s)
    a = at.parent
    C = a.constraints
    W = a.workspace
    periodic_dof_scaling = C.periodic_dof_scaling
    face = a.face
    parent_ldof_parent_dof = W.face_parent_dofs[face]
    parent_dof_free_or_periodic_dof = C.parent_dof_free_or_periodic_dof
    n_parent_ldofs = length(parent_ldof_parent_dof)
    for parent_ldof in 1:n_parent_ldofs
        ldof = parent_ldof
        parent_dof = parent_ldof_parent_dof[parent_ldof]
        if parent_dof < 0
            ldof_s[ldof] = parent_ldof_s[parent_ldof]
            continue
        end
        free_or_periodic_dof = parent_dof_free_or_periodic_dof[parent_dof]
        if free_or_periodic_dof < 0
            periodic_dof = - free_or_periodic_dof
            periodic_scaling = periodic_dof_scaling[periodic_dof]
            ldof_s[ldof] = periodic_scaling*parent_ldof_s[parent_ldof]
        else
            ldof_s[ldof] = parent_ldof_s[parent_ldof]
        end
    end
    ldof_s
end


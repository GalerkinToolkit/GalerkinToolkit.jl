
struct GenericDofGlue{A,B,C}
    face_dofs::A
    face_constraints::B
    free_and_dirichlet::TwoPartPartition{C}
end

face_dofs(a::GenericDofGlue,d) = a.face_dofs[d+1]
face_constraints(a::GenericDofGlue,d) = a.face_constraints[d+1]
free_and_dirichlet(a::GenericDofGlue) = a.free_and_dirichlet

function dof_glue_from_mesh_nodes(mesh,node_to_id)
    @assert ! has_hanging_nodes(mesh)
    @assert ! periodic_nodes(mesh)
    face_dofs = face_nodes(mesh)
    face_constraints = map(i->Fill(I,length(i))face_dofs)
    free_and_dirichlet = partition_from_mask(i->i==0,node_to_id)
    GenericDofGlue(face_dofs,face_constraints,free_and_dirichlet)
end

function classify_mesh_nodes(mesh,names)
    T = Int32
    nnodes = num_nodes(mesh)
    node_to_id = fill(T(INVALID_ID),nnodes)
    has_physical_groups(mesh) || return node_to_id
    D = dimension(mesh)
    for d in D:-1:0
        face_to_nodes = face_nodes(mesh,d)
        for (name2,faces_in_group) in physical_groups(mesh,d)
            id = findfirst(i->i===name2,names)
            if id === nothing
                continue
            end
            for face in faces_in_group
                nodes = face_to_nodes[face]
                node_to_id[nodes] = id
            end
        end
    end
    node_to_id
end

struct DictView{A,B,T} <: AbstractVector{T}
    parent::A
    indices::B
    function DictView(parent::AbstractDict,indices)
        A = typeof(parent)
        B = typeof(indices)
        T = valtype(parent)
        new{A,B,T}(parent,indices)
    end
end

Base.size(a::DictView) = (length(a.indices),)
Base.IndexStyle(::Type{<:DictView}) = IndexLinear()
function Base.getindex(a::DictView,i::Int)
    a.parent[a.indices[i]]
end
function Base.setindex!(a::DictView,v,i::Int)
    a.parent[a.indices[i]] = v
    v
end



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

function classify_mesh_faces!(face_to_tag,mesh,tag_to_name,d)
    has_physical_groups(mesh) || return face_to_tag
    face_groups = physical_groups(mesh,d)
    for (tag,name) in enumerate(tag_to_name)
        for (name2,faces) in face_groups
            if name != name2
                continue
            end
            face_to_tag[faces] .= tag
        end
    end
    face_to_tag
end

function classify_mesh_faces(mesh,tag_to_name,d)
    nfaces = num_faces(mesh,d)
    face_to_tag = zeros(Int32,nfaces)
    classify_mesh_faces!(face_to_tag,mesh,tag_to_name,d)
end

function classify_mesh_nodes!(node_to_tag,mesh,tag_to_name,dmax)
    has_physical_groups(mesh) || return node_to_tag
    for d in dmax:-1:0
        face_to_nodes = face_nodes(mesh,d)
        face_groups = physical_groups(mesh,d)
        for (tag,name) in enumerate(tag_to_name)
            for (name2,faces) in face_groups
                if name != name2
                    continue
                end
                for face in faces
                    nodes = face_to_nodes[face]
                    node_to_tag[nodes] .= tag
                end
            end
        end
    end
    node_to_tag
end

function classify_mesh_nodes(mesh,tag_to_name,dmax)
    T = Int32
    nnodes = num_nodes(mesh)
    node_to_tag = fill(T(INVALID_ID),nnodes)
    classify_mesh_nodes!(node_to_tag,mesh,tag_to_name,dmax)
end

function restrict_face_dofs(nnodes,cell_to_nodes,cells_in_domain)
    node_to_touched = fill(false,nnodes)
    for nodes in view(cell_to_nodes,cells_in_domain)
        node_to_touched[nodes] .= true
    end
    dofs_and_nondofs = partition_from_mask(node_to_touched)
    domain_cell_to_dofs = JaggedArray{Int32,Int32}(view(cell_to_nodes,cells_in_domain))
    node_permutation = permutation(dofs_and_nondofs)
    domain_cell_to_dofs.data[:] = view(node_permutation,domain_cell_to_dofs.data)
    domain_cell_to_dofs, dofs_and_nondofs
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



abstract type AbstractDomain{A} <: GT.AbstractType end
domain(a::AbstractDomain) = a
mesh(a::AbstractDomain) = a.mesh
mesh_id(a::AbstractDomain) = a.mesh_id
physical_names(a::AbstractDomain) = a.physical_names
face_dim(a::AbstractDomain) = GT.val_parameter(a.face_dim)
# TODO two functions for the same
num_dims(a::AbstractDomain) = face_dim(a)
num_ambient_dims(a::AbstractDomain) = num_ambient_dims(mesh(a))
is_reference_domain(a::AbstractDomain) = a.is_reference_domain |> GT.val_parameter
face_around(a::AbstractDomain) = a.face_around

function interior(mesh;
    mesh_id = objectid(mesh),
    physical_names=GT.physical_names(mesh,num_dims(mesh)),
    is_reference_domain = Val(false))
    D = num_dims(mesh)
    domain(mesh;face_dim=D,face_around=1,mesh_id,physical_names,is_reference_domain)
end

function skeleton(mesh;
    mesh_id = objectid(mesh),
    physical_names=GT.physical_names(mesh,num_dims(mesh)-1),
    is_reference_domain = Val(false))
    D = num_dims(mesh)
    domain(mesh;face_dim=D-1,face_around=nothing,mesh_id,physical_names,is_reference_domain)
end

function boundary(mesh::Union{AbstractMesh,PMesh};
    face_around=1,
    mesh_id = objectid(mesh),
    physical_names=GT.physical_names(mesh,num_dims(mesh)-1),
    is_reference_domain = Val(false))
    D = num_dims(mesh)
    domain(mesh;face_dim=D-1,face_around,mesh_id,physical_names,is_reference_domain)
end

function is_boundary(dom::AbstractDomain)
    dom.face_around !== nothing && (num_dims(dom) + 1) == num_dims(mesh(dom))
end

function domain(mesh;
    mesh_id = objectid(mesh),
    face_dim = Val(GT.num_dims(mesh)),
    physical_names=GT.physical_names(mesh,face_dim),
    is_reference_domain = Val(false),
    face_around=nothing,
    cache=nothing,
    )

    if val_parameter(is_reference_domain)
        ReferenceDomain(
                        mesh,
                        mesh_id,
                        physical_names,
                        Val(val_parameter(face_dim)),
                        face_around,
                        cache,
                       )
    else
        PhysicalDomain(
                        mesh,
                        mesh_id,
                        physical_names,
                        Val(val_parameter(face_dim)),
                        face_around,
                        cache,
                       )
    end |> setup_domain
end

function setup_domain(domain::AbstractDomain)
    if domain.cache !== nothing
        return domain
    end
    faces = GT.faces(domain)
    inverse_faces = GT.inverse_faces(domain)
    mesh = domain |> GT.mesh
    d = domain |> GT.face_dim
    refid_to_funs = map(GT.shape_functions,GT.reference_faces(mesh,d))
    cache = domain_cache(;faces,inverse_faces,refid_to_funs)
    replace_cache(domain,cache)
end

function setup_domain(domain::AbstractDomain{<:PMesh})
    domain
end

function domain_cache(;kwargs...)
    (;kwargs...)
end

struct PhysicalDomain{A,B,C,D,E,F} <: AbstractDomain{A}
    mesh::A
    mesh_id::B
    physical_names::C
    face_dim::Val{D}
    face_around::E
    cache::F
end
is_reference_domain(a::PhysicalDomain) = false

struct ReferenceDomain{A,B,C,D,E,F} <: AbstractDomain{A}
    mesh::A
    mesh_id::B
    physical_names::C
    face_dim::Val{D}
    face_around::E
    cache::F
end
is_reference_domain(a::ReferenceDomain) = true

function replace_mesh(domain::AbstractDomain,mesh)
    face_dim = GT.face_dim(domain)
    mesh_id = GT.mesh_id(domain)
    physical_names = GT.physical_names(domain)
    is_reference_domain = GT.is_reference_domain(domain)
    face_around = GT.face_around(domain)
    GT.domain(mesh;face_dim,mesh_id,physical_names,is_reference_domain,face_around)
end

function replace_face_around(domain::AbstractDomain,face_around)
    mesh = GT.mesh(domain)
    face_dim = GT.face_dim(domain)
    mesh_id = GT.mesh_id(domain)
    physical_names = GT.physical_names(domain)
    is_reference_domain = GT.is_reference_domain(domain)
    cache = domain.cache
    GT.domain(mesh;face_dim,mesh_id,physical_names,is_reference_domain,face_around,cache)
end

function replace_cache(domain::AbstractDomain,cache)
    mesh = GT.mesh(domain)
    face_dim = GT.face_dim(domain)
    mesh_id = GT.mesh_id(domain)
    physical_names = GT.physical_names(domain)
    is_reference_domain = GT.is_reference_domain(domain)
    face_around = GT.face_around(domain)
    GT.domain(mesh;face_dim,mesh_id,physical_names,is_reference_domain,face_around,cache)
end

function PartitionedArrays.partition(domain::AbstractDomain{<:PMesh})
    pmesh = GT.mesh(domain)
    map(pmesh.mesh_partition) do mesh
        replace_mesh(domain,mesh)
    end
end

function Base.:(==)(a::AbstractDomain,b::AbstractDomain)
    flag = true
    # TODO check also that one mesh is not a sequential one and the other a parallel one
    flag = flag && (GT.mesh_id(a) == GT.mesh_id(b))
    flag = flag && (GT.physical_names(a) == GT.physical_names(b))
    flag = flag && (GT.face_dim(a) == GT.face_dim(b))
    flag = flag && (GT.is_reference_domain(a) == GT.is_reference_domain(b))
    flag
end

function reference_domain(domain::PhysicalDomain)
    mesh = GT.mesh(domain)
    face_dim = GT.face_dim(domain)
    mesh_id = GT.mesh_id(domain)
    physical_names = GT.physical_names(domain)
    is_reference_domain = Val(true)
    face_around = GT.face_around(domain)
    cache = domain.cache
    GT.domain(mesh;face_dim,mesh_id,physical_names,is_reference_domain,face_around,cache)
end

function reference_domain(domain::ReferenceDomain)
    domain
end

function physical_domain(domain::ReferenceDomain)
    mesh = GT.mesh(domain)
    face_dim = GT.face_dim(domain)
    mesh_id = GT.mesh_id(domain)
    physical_names = GT.physical_names(domain)
    face_around = GT.face_around(domain)
    is_reference_domain = Val(false)
    cache = domain.cache
    GT.domain(mesh;face_dim,mesh_id,physical_names,is_reference_domain,face_around,cache)
end

function physical_domain(domain::PhysicalDomain)
    domain
end

function faces(domain::AbstractDomain)
    if domain.cache !== nothing
        return domain.cache.faces
    end
    mesh = domain |> GT.mesh
    D = GT.face_dim(domain)
    Dface_to_tag = zeros(Int,GT.num_faces(mesh,D))
    tag_to_name = GT.physical_names(domain)
    GT.classify_mesh_faces!(Dface_to_tag,mesh,D,tag_to_name)
    physical_Dfaces = findall(i->i!=0,Dface_to_tag)
    physical_Dfaces
end

function faces(domain::AbstractDomain{<:PMesh})
    map(GT.faces,partition(domain))
end

function num_faces(domain::AbstractDomain)
    length(faces(domain))
end

function num_faces(domain::AbstractDomain{<:PMesh})
    map(GT.num_faces,partition(domain))
end

function inverse_faces(domain::AbstractDomain)
    if domain.cache !== nothing
        return domain.cache.inverse_faces
    end
    d = num_dims(domain)
    ndfaces = num_faces(mesh(domain),d)
    dface_to_face = zeros(Int32,ndfaces)
    face_to_dface = faces(domain)
    dface_to_face[face_to_dface] = 1:length(face_to_dface)
    dface_to_face
end

function inverse_faces(domain::AbstractDomain{<:PMesh})
    map(GT.inverse_faces,partition(domain))
end

# TODO glues not needed anymore.
abstract type AbstractDomainGlue{A} <: GT.AbstractType end
#mesh(a::AbstractDomainGlue) = a.mesh
#domain(a::AbstractDomainGlue) = a.domain
#codomain(a::AbstractDomainGlue) = a.codomain
#
#function PartitionedArrays.partition(a::AbstractDomainGlue{<:PMesh})
#    if hasproperty(a,:face_around)
#        map(partition(GT.domain(a)),partition(GT.codomain(a))) do dom,cod
#            GT.domain_glue(dom,cod;a.face_around)
#        end
#    else
#        map(GT.domain_glue,partition(GT.domain(a)),partition(GT.codomain(a)))
#    end
#end
#
#function domain_glue(domain::AbstractDomain,codomain::AbstractDomain;strict=true)
#    msg = "Trying to combine domains on different meshes"
#    @assert GT.mesh_id(domain) == GT.mesh_id(codomain) msg
#    mesh = GT.mesh(domain)
#    Ddom = GT.face_dim(domain)
#    Dcod = GT.face_dim(codomain)
#    face_around = GT.face_around(domain)
#    if Ddom == Dcod
#        InteriorGlue(mesh,domain,codomain)
#    elseif Ddom < Dcod
#        if face_around === nothing
#            CoboundaryGlue(mesh,domain,codomain)
#        else
#            BoundaryGlue(mesh,domain,codomain)
#        end
#    else
#        if strict
#            error("This case does not make sense")
#        else
#            return nothing
#        end
#    end
#end
#
struct InteriorGlue{A,B,C} <: AbstractDomainGlue{A}
    mesh::A
    domain::B
    codomain::C
end

struct BoundaryGlue{A,B,C} <: AbstractDomainGlue{A}
    mesh::A
    domain::B
    codomain::C
end

struct CoboundaryGlue{A,B,C} <: AbstractDomainGlue{A}
    mesh::A
    domain::B
    codomain::C
end
#
#function target_face(glue::InteriorGlue)
#    mesh = glue |> GT.domain |> GT.mesh
#    domain = glue |> GT.domain
#    codomain = glue |> GT.codomain
#    d = domain |> GT.face_dim
#    D = codomain |> GT.face_dim
#    @assert d == D
#    sface_to_dface = domain |> GT.faces
#    tface_to_Dface = codomain |> GT.faces
#    Dface_to_tface = zeros(Int32,GT.num_faces(mesh,D))
#    tface_to_tface = LinearIndices(tface_to_Dface)
#    Dface_to_tface[tface_to_Dface] = tface_to_tface
#    sface_to_tface = Dface_to_tface[sface_to_dface]
#    nsfaces = length(sface_to_tface)
#    ptrs = collect(Int32,1:(nsfaces+1))
#    face_around = 1
#    sface_to_tfaces = JaggedArray(sface_to_tface,ptrs)
#    sface_to_lfaces = JaggedArray(fill(Int32(1),nsfaces),ptrs)
#    sface_to_faces_around = JaggedArray(fill(Int32(face_around),nsfaces),ptrs)
#    sface_to_tfaces, sface_to_lfaces, sface_to_faces_around
#end
#
#function target_face(glue::CoboundaryGlue)
#    mesh = glue |> GT.domain |> GT.mesh
#    domain = glue |> GT.domain
#    codomain = glue |> GT.codomain
#    d = domain |> GT.face_dim
#    D = codomain |> GT.face_dim
#    @assert d < D
#    sface_to_dface = domain |> GT.faces
#    tface_to_Dface = codomain |> GT.faces
#    Dface_to_tface = zeros(Int32,GT.num_faces(mesh,D))
#    tface_to_tface = LinearIndices(tface_to_Dface)
#    Dface_to_tface[tface_to_Dface] = tface_to_tface
#    topo = GT.topology(mesh)
#    dface_to_lfaces = GT.face_local_faces(topo,d,D)
#    sface_to_lfaces = JaggedArray(view(dface_to_lfaces,sface_to_dface))
#    dface_to_Dfaces = GT.face_incidence(topo,d,D)
#    sface_to_Dfaces = JaggedArray(view(dface_to_Dfaces,sface_to_dface))
#    data = sface_to_Dfaces.data
#    f(Dface) = Dface_to_tface[Dface]
#    data .= f.(data)
#    sface_to_tfaces = sface_to_Dfaces
#    sface_to_faces_around_data = zeros(Int32,length(data))
#    ptrs = sface_to_tfaces.ptrs
#    nsfaces = length(sface_to_tfaces)
#    for sface in 1:nsfaces
#        pini = ptrs[sface]
#        pend = ptrs[sface+1]-1
#        for (ip,p) in enumerate(pini:pend)
#            sface_to_faces_around_data[p] = ip
#        end
#    end
#    sface_to_faces_around = JaggedArray(sface_to_faces_around_data,ptrs)
#    sface_to_tfaces, sface_to_lfaces, sface_to_faces_around
#end
#
#function target_face(glue::BoundaryGlue)
#    domain = replace_face_around(glue.domain,nothing)
#    glue2 = domain_glue(domain,glue.codomain)
#    sface_to_tfaces, sface_to_lfaces, sface_to_faces_around = target_face(glue2)
#    face_around = GT.face_around(glue.domain)
#    sface_to_tface = map(tfaces->tfaces[face_around],sface_to_tfaces)
#    sface_to_lface = map(tfaces->tfaces[face_around],sface_to_lfaces)
#    sface_to_face_around = map(tfaces->tfaces[face_around],sface_to_faces_around)
#    nsfaces = length(sface_to_tface)
#    ptrs = collect(Int32,1:(nsfaces+1))
#    sface_to_tfaces = JaggedArray(sface_to_tface,ptrs)
#    sface_to_lfaces = JaggedArray(sface_to_lface,ptrs)
#    sface_to_faces_around = JaggedArray(sface_to_face_around,ptrs)
#    sface_to_tfaces, sface_to_lfaces, sface_to_faces_around
#end

abstract type AbstractQuantity <: GT.AbstractType end
#mesh(a::AbstractQuantity) = a.mesh
term(a::AbstractQuantity) = a.term
#prototype(a::AbstractQuantity) = a.prototype
prototype(a) = a.prototype
#domain(a::AbstractQuantity) = a.domain
#function PartitionedArrays.partition(a::AbstractQuantity)
#    prototype = a |> GT.prototype
#    map(GT.term(a),partition(GT.domains(a))) do term,domains
#        GT.quantity(term,prototype,domain)
#    end
#end

function quantity(term)
    Quantity(term)
end

struct Quantity{A} <: AbstractQuantity
    term::A
end

function index(;
    domain = nothing,
    face=nothing,
    point=nothing,
    field=nothing,
    dof=nothing,
    face_around=nothing,
    prefix = gensym,
    dict=IdDict{Any,Symbol}(),
    )
    data = (;domain,face,point,field,dof,face_around,dict,prefix)
    Index(data)
end

struct Index{A} <: GT.AbstractType
    data::A
end

function generate_index(dom::AbstractDomain,form_arity=0;
    prefix=gensym,
    field=[prefix("field-in-axis-$i") for i in 1:form_arity],
    face_around = [prefix("face-around-in-axis-$a") for a in 1:form_arity],
    )
    d = num_dims(dom)
    D = num_dims(mesh(dom))
    face = [ prefix("$d-face-dummy") for d in 0:D ]
    point = prefix("point")
    dof = [prefix("dof-in-axis-$i") for i in 1:form_arity]
    face[d+1] = prefix("$d-face")
    index(;domain=dom,face,point,field,dof,face_around,prefix)
end

function num_dims(index::Index)
    num_dims(index.data.domain)
end

function mesh(index::Index)
    index.data.domain.mesh
end

function form_arity(index::Index)
    length(face_around_index(index))
end

function index_storage(index)
    (;( key=>val for (val,key) in index.data.dict )...)
end

function unpack_index_storage(index,state)
    expr = Expr(:block)
    for k in Base.values(index.data.dict) |> collect |> sort
        push!(expr.args,:($k = $state.$k))
    end
    expr
end

function get_symbol!(index,val,name="";prefix=index.data.prefix)
    get!(index.data.dict,val,prefix(name))
end

function face_index(index,d)
    index.data.face[d+1]
end

max_dim(index::Index) = length(index.data.face)-1

function point_index(index)
    index.data.point
end

function face_around_term(index,d,D)
    if index.data.domain.face_around === nothing
        return nothing
    end
    if d == num_dims(index) && (d+1 == D)
        constant_term(index.data.domain.face_around,index)
    else
        nothing
    end
end

function face_around_index(index,a)
    index.data.face_around[a]
end

function face_around_index(index)
    index.data.face_around
end

function field_index(index,a)
    index.data.field[a]
end

function dof_index(index,a)
    index.data.dof[a]
end

function constant_quantity(v)
    quantity() do index
        constant_term(v,index)
        #expr = get_symbol!(index,v,"constant_quantity_value")
        #dim = ANY_DIM
        #(;expr,dim,prototype=v)
    end
end

function face_quantity(data,mesh::AbstractMesh,d;reference=Val(false))
    dim = d
    quantity() do index
        face = face_index(index,dim)
        p = zero(eltype(data))
        if val_parameter(reference)
            #reference_face_term(dim,rid_to_value,index)
            face_to_rid = get_symbol!(index,face_reference_id(mesh,dim),"face_to_rid")
            rid_to_value = get_symbol!(index,data,"rid_to_value")
            expr = @term $rid_to_value[$face_to_refid[$face]]
            expr_term([dim],expr,p,index)
        else
            #sface_to_value = data
            #face_to_sface = inverse_faces(dom)
            #physical_face_term(dim,sface_to_value,face_to_sface,index)
            data_sym = get_symbol!(index,data,"face_to_value_$(dim)d")
            expr = @term $data_sym[$face]
            expr_term([dim],expr,p,index)
        end
        #(;expr,dim,prototype=p)
    end
end

function face_quantity(data,dom::AbstractDomain;reference=Val(false))
    dim = num_dims(dom)
    quantity() do index
        face = face_index(index,dim)
        p = zero(eltype(data))
        if val_parameter(reference)
            #reference_face_term(dim,rid_to_value,index)
            face_to_rid = get_symbol!(index,face_reference_id(mesh(dom),dim),"face_to_rid")
            rid_to_value = get_symbol!(index,data,"rid_to_value")
            expr = @term $rid_to_value[$face_to_refid[$face]]
            expr_term([dim],expr,p,index)
        else
            #sface_to_value = data
            #face_to_sface = inverse_faces(dom)
            #physical_face_term(dim,sface_to_value,face_to_sface,index)
            face_to_sface = inverse_faces(dom)
            data_sym = get_symbol!(index,data,"face_to_value_$(dim)d")
            face_to_sface_sym = get_symbol!(index,face_to_sface,"face_to_sface_$(dim)d")
            expr = @term $data_sym[$face_to_sface_sym[$face]]
            expr_term([dim],expr,p,index)
        end
        #(;expr,dim,prototype=p)
    end
end

function term(a::AbstractQuantity,index)
    a.term(index)
end

#function promote_dim(a)
#    a
#end
#
#function promote_dim(a,b)
#    if a == ANY_DIM
#        b
#    elseif b == ANY_DIM
#        a
#    else
#        a == b ? a : nothing
#    end
#end
#
#function promote_dim(a,args...)
#    dim = promote_dim(args...)
#    if dim === nothing
#        return nothing
#    end
#    promote_dim(a,dim)
#end

#function constant_quantity(v,domain::AbstractDomain{<:PMesh})
#    pmesh = GT.mesh(domain)
#    term = map(pmesh.mesh_partition) do _
#        index -> v
#    end
#    GT.quantity(term,v,domain)
#end

#function index(;
#    face=nothing,
#    local_face=nothing,
#    face_around=nothing,
#    point=nothing,
#    field_per_dim =nothing,
#    dof_per_dim=nothing,
#    face_around_per_dim=nothing,
#    state=gensym("state"),
#    dict=IdDict{Any,Symbol}()
#    )
#    Index(
#          face,
#          local_face,
#          face_around,
#          point,
#          field_per_dim,
#          dof_per_dim,
#          face_around_per_dim,
#          state,
#          dict,
#         )
#end

#struct Index{A,B,D,E,F,G,H} <: AbstractType
#    face::A
#    local_face::B
#    face_around::D
#    point::E
#    field_per_dim::F
#    dof_per_dim::G
#    face_around_per_dim::H
#    state::Symbol #TODO remove
#    dict::IdDict{Any,Symbol}
#end


#function storage(index::Index)
#    (;( key=>val for (val,key) in index.dict )...)
#end


#function replace_face(index::Index,face)
#    Index(
#          face,
#          index.local_face,
#          index.face_around,
#          index.point,
#          index.field_per_dim,
#          index.dof_per_dim,
#          index.face_around_per_dim,
#          index.state,
#          index.dict,
#         )
#end
#
#function replace_local_face(index::Index,local_face)
#    Index(
#          index.face,
#          local_face,
#          index.face_around,
#          index.point,
#          index.field_per_dim,
#          index.dof_per_dim,
#          index.face_around_per_dim,
#          index.state,
#          index.dict,
#         )
#end
#
#function replace_face_around(index::Index,face_around)
#    Index(
#          index.face,
#          index.local_face,
#          face_around,
#          index.point,
#          index.field_per_dim,
#          index.dof_per_dim,
#          index.face_around_per_dim,
#          index.state,
#          index.dict,
#         )
#end
#
#function replace_point(index::Index,point)
#    Index(
#          index.face,
#          index.local_face,
#          index.face_around,
#          point,
#          index.field_per_dim,
#          index.dof_per_dim,
#          index.face_around_per_dim,
#          index.state,
#          index.dict,
#         )
#end
#
#function replace_field_per_dim(index::Index,field_per_dim)
#    Index(
#          index.face,
#          index.local_face,
#          index.face_around,
#          index.point,
#          field_per_dim,
#          index.dof_per_dim,
#          index.face_around_per_dim,
#          index.state,
#          index.dict,
#         )
#end
#
#function replace_dof_per_dim(index::Index,dof_per_dim)
#    Index(
#          index.face,
#          index.local_face,
#          index.face_around,
#          index.point,
#          index.field_per_dim,
#          dof_per_dim,
#          index.face_around_per_dim,
#          index.state,
#          index.dict,
#         )
#end

function return_prototype(f,args...)
    f(args...)
end

function return_prototype(::typeof(getindex),a,b)
    zero(eltype(a))
end

function return_prototype(::typeof(getindex),a,b::Integer)
    if length(a) != 0
        first(a)
    else
        zero(eltype(a))
    end
end

function call(f,args...)
    f(args...)
end

function get_symbol!(index,val::typeof(call),name="";prefix=index.data.prefix)
    :(GalerkinToolkit.call)
end

function call(g,q::AbstractQuantity)
    quantity() do index
        unary_call_term(g,term(q,index))
        #g_sym = get_symbol!(index,g,"callee")
        #tq = term(q,index)
        #p = GT.return_prototype(g,prototype(tq))
        #expr = :($g_sym($(tq.expr)))
        #(;expr,dim=tq.dim,prototype=p)
    end
end

function call(g,a::AbstractQuantity,b::AbstractQuantity)
    quantity() do index
        ta = term(a,index)
        tb = term(b,index)
        t = binary_call_term(g,term(a,index),term(b,index))
        dims = union(free_dims(ta),free_dims(tb))
        if length(dims) <= 1
            return t
        end
        @assert length(dims) == 2
        dim = minimum(dims)
        dim2 = maximum(dims)
        cell_around = face_around_term(index,dim,dim2)
        if cell_around !== nothing
            boundary_term(dim,dim2,t,cell_around)
        else
            cell_around = gensym("cell_around")
            skeleton_term(dim,dim2,t,cell_around)
        end
        #g_sym = get_symbol!(index,g,"callee")
        #ta = term(a,index)
        #tb = term(b,index)
        #p = return_prototype(g,prototype(ta),prototype(tb))
        #expr = :($g_sym($(ta.expr),$(tb.expr)))
        #dim = promote_dim(ta.dim,tb.dim)
        #p2 = p
        #if dim == nothing
        #    dim = min(ta.dim,tb.dim)
        #    dim2 = max(ta.dim,tb.dim)
        #    dummy = face_index(index,dim2)
        #    face = face_index(index,dim)
        #    #faces = actual_faces_index(index,dim,dim2)
        #    topo = topology(mesh(index))
        #    face_to_cells = get_symbol!(index,face_incidence(topo,dim,dim2),"face_to_cells")
        #    axis_to_face_around = face_around_index(index)
        #    #dom = target_domain(index)
        #    #gluable = face_around(dom) !== nothing && num_dims(dom) == dim && dim2 == dim+1
        #    #if gluable
        #    #    face = face(index,dim)
        #    #    face_to_sface = get_symbol!(index,inverse_faces(dom))
        #    #    sface_to_cell_around = get_symbol!(index,face_around(dom))
        #    #    topo = topology(mesh(dom))
        #    #    face_to_cells = get_symbol!(index,face_incidence(topo,dim,dim2))
        #    #    expr = @term begin
        #    #        sface = $face_to_sface[$face]
        #    #        cell_around = $face_to_cells[$sface]


        #    #        fun = $dummy -> $expr
        #    #        fun()
        #    #    end
        #    cells = @term $face_to_cells[$face]
        #    cell_around = face_around_term(index,dim,dim2)
        #    cell_around_sym = get_symbol!(index,cell_around,"cell_around")
        #    if isa(cell_around,Int)
        #        #expr = @term GalerkinToolkit.call($dummy -> $expr,$cells[$cell_around_sym])
        #        actual = :($cells[$cell_around_sym])
        #        expr = substitute(expr,dummy=>actual)
        #        p2 = p
        #    elseif isa(cell_around,AbstractArray)
        #        face_to_sface = get_symbol!(inverse_faces(domain(index)),"face_to_sface")
        #        actual = :($cells[$cell_around_sym[$face_to_sface[$face]]])
        #        expr = substitute(expr,dummy=>actual)
        #        #expr = @term GalerkinToolkit.call($dummy -> $expr,$cells[$cell_around_sym[$face_to_sface[$face]]])
        #        p2 = p
        #    elseif cell_around !== nothing
        #        error()
        #    elseif length(axis_to_face_around) == 0
        #        expr = @term begin
        #            fun = $dummy -> $expr
        #            map(fun,$cells)
        #        end
        #        p2 = [p,]
        #    else
        #        i = gensym("i-dummy")
        #        exprs = map(axis_to_face_around) do face_around
        #            :(LinearAlgebra.I[$i,$face_around])
        #        end
        #        if length(exprs) > 1
        #            deltas = Expr(:call,:tuple,exprs...)
        #            mask = :(prod($deltas))
        #        else
        #            mask = exprs[1]
        #        end
        #        expr = @term begin
        #            fun = ($i,$dummy) -> $expr*$mask
        #            map(fun,1:length($cells),$cells)
        #        end
        #        p2 = [p,]
        #    end
        #end
        #(;expr,dim,prototype=p2)
    end
end

function call(g,args::AbstractQuantity...)
    quantity() do index
        multivar_call_term(g,map(a->term(a,index),args))
        #g_sym = get_symbol!(index,v,gensym("callee"))
        #ts = map(a->term(a,index),args)
        #p = GT.return_prototype(g,(map(GT.prototype,ts)...))
        #dims = map(t->ts.dim,ts)
        #dim = promote_dim(dims...)
        #@assert dim !== nothing
        #exprs = map(t->ts.expr,ts)
        ## We use call since @rule is not able to
        ## identify function calls on variable slots
        #expr = @term call($g_sym,$(exprs...))
        #(;expr,dim,prototype=p)
    end
end

#function call(g::AbstractQuantity,args::AbstractQuantity...)
#    fs = map(GT.term,args)
#    domains = union(domains(g),map(domains,args)...)
#    #msg = "All quantities need to be defined on the same domain"
#    #@assert all(dom->dom==domain,map(GT.domain,args)) msg
#    # TODO check everything except reference/physical domain?
#    # Maybe check reference/physical domain only when evaluating functions?
#    prototype = GT.return_prototype(GT.prototype(g),(map(GT.prototype,args)...))
#    g_term = GT.term(g)
#    GT.quantity(prototype,domains) do indices,dict
#        f_exprs = map(f->f(indices,dict),fs)
#        g_expr = g_term(indices,dict)
#        # We use call since @rule is not able to
#        # identify function calls on variable slots
#        @term call($g_expr,$(f_exprs...))
#        #:($g_expr($(f_exprs...)))
#    end
#end

#function call(g::AbstractQuantity,args::AbstractQuantity{<:PMesh}...)
#    pg = partition(g)
#    pargs = map(partition,args)
#    q = map(pg,pargs...) do myg,myargs...
#        GT.call(myg,myargs...)
#    end
#    domain = args |> first |> GT.domain
#    term = map(GT.term,q)
#    prototype = map(GT.prototype,q) |> PartitionedArrays.getany
#    GT.quantity(term,prototype,domain)
#end

#function call(g,args::AbstractQuantity...)
#    g_qty = GT.constant_quantity(g)
#    call(g_qty,args...)
#end

#function (f::AbstractQuantity)(x::AbstractQuantity)
#    domain = GT.domain(x)
#    codomain = GT.domain(f)
#    flag = physical_domain(domain) == physical_domain(codomain)
#    if flag
#        call(f,x)
#    else
#        align_and_call(f,x)
#    end
#end

function (f::AbstractQuantity)(x::AbstractQuantity)
    call(call,f,x)
end

#function align_and_call(f,x)
#    domain = GT.domain(x)
#    codomain = GT.domain(f)
#    glue = GT.domain_glue(domain,codomain)
#    align_and_call(f,x,glue)
#end
#
#function align_and_call(f,x,glue::InteriorGlue)
#    g = GT.align_field(f,GT.domain(glue))
#    call(g,x)
#end
#
#function align_and_call(f,x,glue::BoundaryGlue)
#    g = GT.align_field(f,GT.domain(glue))
#    call(g,x)
#end
#
#function align_and_call(f,x,glue::CoboundaryGlue)
#    g = GT.align_field(f,GT.domain(glue))
#    call(g,x)
#end

#function Base.getindex(q::AbstractQuantity,::typeof(+))
#    face_around = 1
#    term = GT.term(q)
#    GT.quantity(GT.prototype(q),GT.domain(q)) do index
#        index2 = GT.replace_face_around(index,face_around)
#        term(index2)
#    end
#end
#
#function Base.getindex(q::AbstractQuantity,::typeof(-))
#    face_around = 2
#    term = GT.term(q)
#    GT.quantity(GT.prototype(q),GT.domain(q)) do index
#        index2 = GT.replace_face_around(index,face_around)
#        term(index2)
#    end
#end

function analytical_field(f)
    constant_quantity(f)
end

function face_constant_field(data,dom::AbstractDomain)
    q = face_quantity(data,dom)
    call(a->(x->a),q)
end

function physical_map(mesh::AbstractMesh,d)
    quantity() do index
        dval = Val(val_parameter(d))
        physical_map_term(dval,index)
    end
end

function inverse_physical_map(mesh::AbstractMesh,d)
    x0 = zero(SVector{val_parameter(d),Float64})
    x = constant_quantity(x0)
    phi = physical_map(mesh,d)
    call(inv_map,phi,x)
end

#function face_map(mesh::AbstractMesh,d)
#    node_to_x = node_coordinates(mesh)
#    x0 = zero(eltype(node_to_x))
#    p = x -> x0
#    quantity() do index
#        drefid_to_refdface = GT.reference_faces(mesh,d)
#        dface_to_drefid = face_reference_id(mesh,d)
#        caller = value
#        dof = gensym("geom_dof_dummy")
#        funs = reference_shape_function_term(d,drefid_to_refdface,dface_to_drefid,caller,dof,index)
#        coeffs = face_dof_term(d,node_to_coord,dface_to_nodes,dof,index)
#        ndofs = reference_face_term(d,drefid_to_ndofs,dface_to_drefid,index)
#        discrete_function_term(coeff,funs,ndofs)
#        #node_to_coord = get_symbol!(index,node_to_x,"node_to_x")
#        #drefid_to_refdface = GT.reference_faces(mesh,d)
#        #drefid_to_dof_to_shape = get_symbol!(index,map(shape_functions,drefid_to_refdface),"refid_to_dof_to_shape")
#        #dface_to_nodes = get_symbol!(index,face_nodes(mesh,d),"face_to_nodes")
#        #dface_to_drefid = get_symbol!(index,face_reference_id(mesh,d),"face_to_refid")
#        #dof = gensym("dummy_dof_geom")
#        #x = gensym("dummy_x")
#        #dface = face_index(index,d)
#        #expr = @term begin
#        #    node = $dface_to_nodes[$dface][$dof]
#        #    coord = $node_to_coord[node]
#        #    drefid = $dface_to_drefid[$dface]
#        #    dof_to_shape = $drefid_to_dof_to_shape[drefid]
#        #    fun = $dof -> coord*$dof_to_shape[$dof]($x)
#        #    $x -> sum(fun,1:length($dof_to_shape))
#        #end
#        #(;expr,dim=d,prototype=p)
#    end
#end

function point_quantity(data,dom::AbstractDomain;
    reference=Val(false),
    face_reference_id = face_reference_id(mesh(dom),num_dims(dom)),
    )
    dim = num_dims(dom)
    #p = zero(eltype(eltype(data)))
    quantity() do index
        if val_parameter(reference)
            reference_point_term(dim,data,face_reference_id,index)
        else
            physical_point_term(dim,data,inverse_faces(dom),index)
        end
        #face = face_index(index,dim)
        #point = point_index(index)
        #expr = if val_parameter(reference)
        #    face_to_rid = get_symbol!(index,face_reference_id(mesh(dom),dim),"face_to_rid")
        #    rid_to_points = get_symbol!(index,data,"rid_to_points")
        #    @term begin
        #        rid = $face_to_rid[$face]
        #        $rid_to_points[rid][$point]
        #    end
        #else
        #    face_to_sface = get_symbol!(index,inverse_faces(dom),"face_to_sface")
        #    sface_to_points = get_symbol!(index,data,"sface_to_points")
        #    @term begin
        #        sface = $face_to_sface[$face]
        #        $sface_to_points[sface][$point]
        #    end
        #end
        #(;expr,dim,prototype=p)
    end
end

function shape_function_quantity(data,dom::AbstractDomain;
    reference=Val(false),
    face_reference_id = face_reference_id(mesh(dom),num_dims(dom)),
    dof = gensym("dummy-dof")
    )
    dim = num_dims(dom)
    #p = zero(eltype(eltype(data)))
    quantity() do index
        if val_parameter(reference)
            reference_shape_function_term(dim,data,face_reference_id,dof,index)
        else
            error("Not implemented, but possible to implement it.")
        end
    end
end

function reference_map(mesh::AbstractMesh,d,D)
    quantity() do index
        dval = Val(val_parameter(d))
        Dval = Val(val_parameter(D))
        t = reference_map_term(dval,Dval,index)
        cell_around = face_around_term(index,d,D)
        if cell_around !== nothing
            boundary_term(d,D,t,cell_around)
        else
            cell_around = gensym("cell_around")
            skeleton_term(d,D,t,cell_around)
        end
    end
end

#function face_map(mesh::AbstractMesh,d,D)
#    #x0 = zero(SVector{D,Float64})
#    #p = x -> x0
#    #topo = topology(mesh)
#    #a, b = GT.face_incidence_ext(topo,d,D)
#    #perms = GT.face_permutation_ids(topo,D,d)
#    #    Drefid_to_refDface = GT.reference_faces(mesh,D)
#    #    drefid_to_refdface = GT.reference_faces(mesh,d)
#    #coords = map(rfe->face_node_coordinates(rfe,d),Drefid_to_refDface)
#    #funs = map(shape_functions,drefid_to_refdface)
#    #ridsD = GT.face_reference_id(mesh,D)
#    #ridsd = GT.face_reference_id(mesh,d)
#    quantity() do index
#        drefid_to_refdface = GT.reference_faces(mesh,d)
#        dface_to_drefid = face_reference_id(mesh,d)
#        caller = value
#        dof = gensym("geom_dof_dummy")
#        funs = reference_shape_function_term(d,drefid_to_refdface,dface_to_drefid,caller,dof,index)
#        coeffs = reference_dof_term()
#
#
#
#
#        dface = face_index(index,d)
#        dim = d
#        dface_to_Dfaces = get_symbol!(index,a,"dface_to_Dfaces")
#        dface_to_ldfaces = get_symbol!(index,b,"dface_to_ldfaces")
#        Dface_to_ldface_to_perm = get_symbol!(index,perms,"Dface_to_ldface_to_perm")
#        Dface_to_Drefid = get_symbol!(index,ridsD,"Dface_to_Drefid")
#        dface_to_drefid = get_symbol!(index,ridsd,"dface_to_drefid")
#        Drefid_to_ldface_to_perm_to_coords = get_symbol!(index,coords,"Drefid_to_ldface_to_perm_to_coords")
#        drefid_to_dof_to_shape = get_symbol!(index,funs,"drefid_to_dof_to_shape")
#        dof = gensym("dummy_dof_geom")
#        x = gensym("dummy_x")
#        Dface_around = gensym("dummy_Dface_around")
#        expr = @term begin
#            Dface = $dface_to_Dfaces[$dface][$Dface_around]
#            ldface = $dface_to_ldfaces[$dface][$Dface_around]
#            Drefid = $Dface_to_Drefid[Dface]
#            perm = $Dface_to_ldface_to_perm[Dface][ldface]
#            coeff = $Drefid_to_ldface_to_perm_to_coords[Drefid][ldface][perm][$dof]
#            drefid = $dface_to_drefid[$dface]
#            dof_to_shape = $drefid_to_dof_to_shape[drefid]
#            fun = $dof -> coeff*$dof_to_shape[$dof]($x)
#            $x -> sum(fun,1:length($dof_to_shape))
#        end
#        Dface_around_target = face_around_term(index,d,D)
#        Dface_around_actual = get_symbol!(index,Dface_around_target,"Dface_around")
#        if d == D
#            expr = expr
#            p2 = p
#        elseif isa(Dface_around_target,Int)
#            expr = substitute(expr,Dface_around=>Dface_around_actual)
#            p2 = p
#        elseif isa(Dface_around_target,AbstractArray)
#            face_to_sface = get_symbol!(inverse_faces(domain(index)),"face_to_sface")
#            actual = :($Dface_around_actual[$face_to_sface[$dface]])
#            expr = substitute(expr,Dface_around=>actual)
#            p2 = [p]
#        else
#            expr = @term begin
#                map($Dface_around -> $expr,1:length($dface_to_Dfaces[$dface]))
#            end
#            p2 = [p]
#        end
#        (;expr,dim,prototype=p2)
#    end
#end

function face_node_coordinates(refDface,d)
    D = num_dims(refDface)
    if d == D
        [[node_coordinates(refDface)]]
    else
        boundary = refDface |> GT.geometry |> GT.boundary
        lface_to_nodes = GT.face_nodes(boundary,d)
        node_to_coords = GT.node_coordinates(boundary)
        lface_to_lrefid = GT.face_reference_id(boundary,d)
        lrefid_to_lrefface = GT.reference_faces(boundary,d)
        lrefid_to_perm_to_ids = map(GT.node_permutations,lrefid_to_lrefface)
        map(1:GT.num_faces(boundary,d)) do lface
            lrefid = lface_to_lrefid[lface]
            nodes = lface_to_nodes[lface]
            perm_to_ids = lrefid_to_perm_to_ids[lrefid]
            map(perm_to_ids) do ids
                coords = node_to_coords[nodes[ids]]
                coords
            end
        end
    end
end

function reference_map(refdface::AbstractMeshFace,refDface::AbstractMeshFace)
    d = num_dims(refdface)
    dof_to_f = shape_functions(refdface)
    boundary = refDface |> GT.geometry |> GT.boundary
    lface_to_nodes = GT.face_nodes(boundary,d)
    node_to_coords = GT.node_coordinates(boundary)
    lface_to_lrefid = GT.face_reference_id(boundary,d)
    lrefid_to_lrefface = GT.reference_faces(boundary,d)
    lrefid_to_perm_to_ids = map(GT.node_permutations,lrefid_to_lrefface)
    map(1:GT.num_faces(boundary,d)) do lface
        lrefid = lface_to_lrefid[lface]
        nodes = lface_to_nodes[lface]
        perm_to_ids = lrefid_to_perm_to_ids[lrefid]
        map(perm_to_ids) do ids
            dof_to_coeff = node_to_coords[nodes[ids]]
            ndofs = length(dof_to_coeff)
            x -> sum(dof->dof_to_coeff[dof]*dof_to_f[dof](x),1:ndofs)
        end
    end
end

#function face_local_map(refDface,refDface)
#end

function face_incidence_ext(topo,d,D)
    dface_to_Dface_around_to_Dface = GT.face_incidence(topo,d,D) |> JaggedArray
    Dface_to_ldface_to_dface = GT.face_incidence(topo,D,d)
    data = copy(dface_to_Dface_around_to_Dface.data)
    fill!(data,0)
    ptrs = dface_to_Dface_around_to_Dface.ptrs
    dface_to_Dface_around_to_ldface = JaggedArray(data,ptrs)
    ndfaces = num_faces(topo,d)
    for dface in 1:ndfaces
        Dface_around_to_Dface = dface_to_Dface_around_to_Dface[dface]
        for (Dface_around, Dface) in enumerate(Dface_around_to_Dface)
            ldface_to_dface = Dface_to_ldface_to_dface[Dface]
            for (ldface2,dface2) in enumerate(ldface_to_dface)
                if dface == dface2
                    dface_to_Dface_around_to_ldface[dface][Dface_around] = ldface2
                    break
                end
            end
        end
    end
    dface_to_Dface_around_to_Dface,dface_to_Dface_around_to_ldface
end

#function domain_map(domain::AbstractDomain,codomain::AbstractDomain)
#    d = num_dims(domain)
#    D = num_dims(codomain)
#    topo = topology(mesh(domain))
#    a, b = GT.face_incidence_ext(topo,d,D)
#    dface_to_Dfaces = get_symbol!(index,a,"dface_to_Dfaces")
#    dface_to_ldfaces = get_symbol!(index,b,"dface_to_ldfaces")
#    Dface_to_ldface_to_perm = get_symbol!(index,GT.face_permutation_ids(topo,D,d),"Dface_to_ldface_to_perm")
#    Dface_to_Drefid = get_symbol!(index,GT.face_reference_id(mesh,D),"Dface_to_Drefid")
#    dface_to_drefid = get_symbol!(index,GT.face_reference_id(mesh,d),"dface_to_drefid")
#    Drefid_to_refDface = GT.reference_faces(mesh,D)
#    drefid_to_refdface = GT.reference_faces(mesh,d)
#    Drefid_to_ldface_to_perm_to_coords = get_symbol!(index,map(face_node_coordinates,Drefid_to_refDface),"Drefid_to_ldface_to_perm_to_coords")
#    drefid_to_dof_to_shape = get_symbol!(index,map(shape_functions,drefid_to_refdface),"drefid_to_dof_to_shape")
#    dof = gensym("dummy_dof_geom")
#    x = gensym("dummy_x")
#    quantity(p) do index
#        dface = face_index(index,d)
#        dim = d
#        expr = @term begin
#            Dface = $dface_to_Dfaces[$dface][$Dface_around]
#            ldface = $dface_to_ldfaces[$dface][$Dface_around]
#            Drefid = $Dface_to_Drefid[Dface]
#            perm = $Dface_to_ldface_to_perm[Dface][ldface]
#            coeff = $Drefid_to_ldface_to_perm_to_coords[Drefid][ldface][perm][$dof]
#            drefid = $dface_to_drefid[$dface]
#            dof_to_shape = $drefid_to_dof_to_shape[drefid]
#            fun = $dof -> coeff*$dof_to_shape[$dof]($x)
#            $x -> sum(fun,1:length($dof_to_shape))
#        end
#        Dface_around_target = face_around_term(index,d,D)
#        Dface_around_actual = get_symbol!(index,Dface_around_target,"Dface_around")
#        if d == D
#            expr = expr
#        elseif isa(Dface_around_target,Int)
#            expr = substitute(expr,Dface_around=>Dface_around_actual)
#        elseif isa(Dface_around_target,AbstractArray)
#            face_to_sface = get_symbol!(inverse_faces(domain(index)),"face_to_sface")
#            actual = :($Dface_around_actual[$face_to_sface[$dface]])
#            expr = substitute(expr,Dface_around=>actual)
#        else
#            expr = @term begin
#                map($Dface_around -> $expr,$dface_to_Dfaces[$dface])
#            end
#        end
#        (;expr,dim)
#    end
#end
#

#function domain_map(domain::AbstractDomain,codomain::AbstractDomain)
#    glue = GT.domain_glue(domain,codomain)
#    domain_map(glue,domain,codomain)
#end

#function domain_map(domain::AbstractDomain{<:PMesh},codomain::AbstractDomain{<:PMesh})
#    q = map(GT.domain_map,partition(domain),partition(codomain))
#    prototype = map(GT.prototype,q) |> PartitionedArrays.getany
#    term = map(GT.term,q)
#    GT.quantity(term,prototype,domain)
#end

#function domain_map(glue::InteriorGlue,::ReferenceDomain,::ReferenceDomain)
#    domain = glue.domain
#    prototype = identity
#    term = indices,dict -> :identity
#    GT.quantity(term,prototype,domain)
#end
#
#function domain_map(glue::InteriorGlue,::PhysicalDomain,::PhysicalDomain)
#    domain = glue.domain
#    prototype = identity
#    term = indices,dict -> :identity
#    GT.quantity(term,prototype,domain)
#end

#function linear_combination(funs,coefs)
#    q -> begin
#        sum(1:length(coefs)) do i
#            x = coefs[i]
#            fun = funs[i]
#            coeff = fun(q)
#            coeff*x
#        end
#    end
#end

#function linear_combination(funs,coefs,ids)
#    q -> begin
#        sum(1:length(ids)) do i
#            x = coefs[ids[i]]
#            fun = funs[i]
#            coeff = fun(q)
#            coeff*x
#        end
#    end
#end

# TODO a better name
#function face_function(rid_to_fs,face_to_rid,face_to_dofs,dofs_to_value,face)
#    fs = reference_value(rid_to_fs,face_to_rid,face)
#    dofs = face_to_dofs[face]
#    values = view(dofs_to_value,dofs)
#    linear_combination(fs,values)
#end

#function reference_point(rid_to_x,face_to_rid,face,point)
#    reference_value(rid_to_x,face_to_rid,face)[point]
#end

#function reference_tabulators(rid_to_fs,rid_to_xs)
#    map((fs,xs)->map(x->map(f->f(x),fs),xs),rid_to_fs,rid_to_xs)
#end

#function reference_tabulator(fs,xs)
#    f = fs[1]
#    z = zero(eltype(xs))
#    T = typeof(f(z))
#    nx = length(xs)
#    nf = length(fs)
#    A = Matrix{T}(undef,nf,nx)
#    for j in 1:nx
#        for i in 1:nf
#            A[i,j] = fs[i](xs[j])
#        end
#    end
#    A
#end
#
#function face_function_value(rid_to_tab,face_to_rid,face_to_dofs,dofs_to_value,face,point)
#    tab = reference_value(rid_to_tab,face_to_rid,face)
#    dofs = face_to_dofs[face]
#    values = view(dofs_to_value,dofs)
#    n = length(values)
#    sum(i->tab[i,point]*values[i],1:n)
#end
#
#function jacobian_face_function_value(rid_to_tab,face_to_rid,face_to_dofs,dofs_to_value,face,point)
#    tab = reference_value(rid_to_tab,face_to_rid,face)
#    dofs = face_to_dofs[face]
#    values = view(dofs_to_value,dofs)
#    n = length(values)
#    sum(i-> values[i] * transpose(tab[i,point]),1:n) # TODO: optimize it by precomputing the transpose or unrolling the loop
#end
#
#
#function gradient_reference_tabulator(fs,xs)
#    # TODO: generalize this function
#    f = fs[1]
#    z = zero(eltype(xs))
#    T = typeof(ForwardDiff.gradient(f, z))
#    nx = length(xs)
#    nf = length(fs)
#    A = Matrix{T}(undef,nf,nx)
#    for j in 1:nx
#        for i in 1:nf
#            A[i,j] = ForwardDiff.gradient(fs[i], xs[j])
#        end
#    end
#    A
#end
#
#function face_shape_function_value(rid_to_tab, face_to_rid, face, dof, face_to_rid2, sface, point)
#    # @assert face_to_rid[face] == face_to_rid2[sface]
#    tab = reference_value(rid_to_tab,face_to_rid,face)
#    tab[dof, point]
#end
#
#function reference_value(rid_to_value,face_to_rid,face)
#    rid = face_to_rid[face]
#    value = rid_to_value[rid]
#    value
#end

#function domain_map(glue::InteriorGlue,::ReferenceDomain,::PhysicalDomain)
#    domain = glue.domain
#    mesh = domain |> GT.mesh
#    d = domain |> GT.face_dim
#    node_to_coords = GT.node_coordinates(mesh)
#    sface_to_face = domain |> GT.faces
#    face_to_nodes = GT.face_nodes(mesh,d)
#    face_to_refid = GT.face_reference_id(mesh,d)
#    if haskey(domain.cache,:refid_to_funs)
#        refid_to_funs = domain.cache.refid_to_funs
#    else
#        refid_to_refface = GT.reference_faces(mesh,d)
#        refid_to_funs = map(GT.shape_functions,refid_to_refface)
#    end
#    T = eltype(GT.node_coordinates(mesh))
#    x = zero(T)
#    prototype = y->x
#    GT.quantity(prototype,domain) do indices,dict
#        index = indices[domain]
#        sface = index.face
#        state = index.state
#        sface_to_face_sym = get!(dict,sface_to_face,gensym("sface_to_face"))
#        refid_to_funs_sym = get!(dict,refid_to_funs,gensym("refid_to_funs"))
#        face_to_nodes_sym = get!(dict,face_to_nodes,gensym("face_to_nodes"))
#        node_to_coords_sym = get!(dict,node_to_coords,gensym("node_to_coords"))
#        face_to_refid_sym = get!(dict,face_to_refid,gensym("face_to_refid"))
#        if isa(sface,FaceList)
#        else
#            @term begin
#                face_function(
#                              $refid_to_funs_sym,
#                              $face_to_refid_sym,
#                              $face_to_nodes_sym,
#                              $node_to_coords_sym,
#                              $sface_to_face_sym[$sface])
#            end
#        end
#    end
#end
#
#function domain_map(glue::InteriorGlue,::PhysicalDomain,::ReferenceDomain)
#    error("Physical to reference map not implemented yet")
#end
#
#function domain_map(glue::CoboundaryGlue,::PhysicalDomain,::PhysicalDomain)
#    error("Case not yet implemented")
#end
#
#function domain_map(glue::CoboundaryGlue,::ReferenceDomain,::ReferenceDomain)
#    domain = glue.domain
#    codomain = glue |> GT.codomain
#    mesh = codomain |> GT.mesh
#    D = codomain |> GT.face_dim
#    Drefid_to_refDface = GT.reference_faces(mesh,D)
#    refDface = first(Drefid_to_refDface)
#    boundary = refDface |> GT.geometry |> GT.boundary
#    node_to_coords = GT.node_coordinates(boundary)
#    T = eltype(node_to_coords)
#    x = zero(T)
#    prototype = y -> x
#    sface_to_tfaces, sface_to_lfaces, = glue |> GT.target_face
#    tface_to_Dface = codomain |> GT.faces
#    d = domain |> GT.face_dim
#    topo = mesh |> GT.topology
#    Dface_to_lface_to_perm = GT.face_permutation_ids(topo,D,d)
#    Dface_to_Drefid = GT.face_reference_id(mesh,D)
#    Drefid_to_refDface = GT.reference_faces(mesh,D)
#    Drefid_to_lface_to_perm_to_coords = map(Drefid_to_refDface) do refDface
#        boundary = refDface |> GT.geometry |> GT.boundary
#        lface_to_nodes = GT.face_nodes(boundary,d)
#        node_to_coords = GT.node_coordinates(boundary)
#        lface_to_lrefid = GT.face_reference_id(boundary,d)
#        lrefid_to_lrefface = GT.reference_faces(boundary,d)
#        lrefid_to_perm_to_ids = map(GT.node_permutations,lrefid_to_lrefface)
#        map(1:GT.num_faces(boundary,d)) do lface
#            lrefid = lface_to_lrefid[lface]
#            nodes = lface_to_nodes[lface]
#            perm_to_ids = lrefid_to_perm_to_ids[lrefid]
#            map(perm_to_ids) do ids
#                coords = node_to_coords[nodes[ids]]
#                coords
#            end
#        end
#    end
#    sface_to_dface = domain |> GT.faces
#    dface_to_drefid = GT.face_reference_id(mesh,d)
#    drefid_to_refdface = GT.reference_faces(mesh,d)
#    drefid_to_funs = map(GT.shape_functions,drefid_to_refdface)
#    GT.quantity(prototype,domain) do index
#        sface = index.face
#        dict = index.dict
#        sface_to_tfaces_sym = get!(dict,sface_to_tfaces,gensym("sface_to_tfaces"))
#        sface_to_lfaces_sym = get!(dict,sface_to_lfaces,gensym("sface_to_lfaces"))
#        tface_to_Dface_sym = get!(dict,tface_to_Dface,gensym("tface_to_Dface"))
#        sface_to_dface_sym = get!(dict,sface_to_dface,gensym("sface_to_dface"))
#        Dface_to_lface_to_perm_sym = get!(dict,Dface_to_lface_to_perm,gensym("Dface_to_lface_to_perm"))
#        Dface_to_Drefid_sym = get!(dict,Dface_to_Drefid,gensym("Dface_to_Drefid"))
#        dface_to_drefid_sym = get!(dict,dface_to_drefid,gensym("dface_to_drefid"))
#        Drefid_to_lface_to_perm_to_coords_sym = get!(dict,Drefid_to_lface_to_perm_to_coords,gensym("Drefid_to_lface_to_perm_to_coords"))
#        drefid_to_funs_sym = get!(dict,drefid_to_funs,gensym("drefid_to_funs"))
#        face_around = index.face_around
#        @assert face_around !== nothing
#        @term begin
#            boundary_face_function(
#                $sface,
#                $sface_to_tfaces_sym,
#                $sface_to_lfaces_sym,
#                $tface_to_Dface_sym,
#                $sface_to_dface_sym,
#                $Dface_to_lface_to_perm_sym,
#                $Dface_to_Drefid_sym,
#                $dface_to_drefid_sym,
#                $Drefid_to_lface_to_perm_to_coords_sym,
#                $drefid_to_funs_sym,
#                $face_around,
#           )
#        end
#    end
#end
#
#function domain_map(glue::CoboundaryGlue,::ReferenceDomain,::PhysicalDomain)
#    error("Case not yet implemented")
#end
#
#function domain_map(glue::CoboundaryGlue,::PhysicalDomain,::ReferenceDomain)
#    error("Case not yet implemented")
#end
#
#function domain_map(glue::BoundaryGlue,::PhysicalDomain,::PhysicalDomain)
#    error("Case not yet implemented")
#end
#
#function boundary_face_function(
#        sface,
#        sface_to_tfaces,
#        sface_to_lfaces,
#        tface_to_Dface,
#        sface_to_dface,
#        Dface_to_lface_to_perm,
#        Dface_to_Drefid,
#        dface_to_drefid,
#        Drefid_to_lface_to_perm_to_coords,
#        drefid_to_funs,
#        face_around,
#    )
#    tfaces = sface_to_tfaces[sface]
#    lfaces = sface_to_lfaces[sface]
#    tface = tfaces[face_around]
#    lface = lfaces[face_around]
#    Dface = tface_to_Dface[tface]
#    dface = sface_to_dface[sface]
#    perm = Dface_to_lface_to_perm[Dface][lface]
#    Drefid = Dface_to_Drefid[Dface]
#    drefid = dface_to_drefid[dface]
#    coords = Drefid_to_lface_to_perm_to_coords[Drefid][lface][perm]
#    funs = drefid_to_funs[drefid]
#    linear_combination(funs,coords)
#end
#
#function domain_map(glue::BoundaryGlue,::ReferenceDomain,::ReferenceDomain)
#    domain = glue.domain
#    codomain = glue |> GT.codomain
#    mesh = codomain |> GT.mesh
#    D = codomain |> GT.face_dim
#    Drefid_to_refDface = GT.reference_faces(mesh,D)
#    refDface = first(Drefid_to_refDface)
#    boundary = refDface |> GT.geometry |> GT.boundary
#    node_to_coords = GT.node_coordinates(boundary)
#    T = eltype(node_to_coords)
#    x = zero(T)
#    prototype = y -> x
#    sface_to_tfaces, sface_to_lfaces, = glue |> GT.target_face
#    tface_to_Dface = codomain |> GT.faces
#    d = domain |> GT.face_dim
#    topo = mesh |> GT.topology
#    Dface_to_lface_to_perm = GT.face_permutation_ids(topo,D,d)
#    Dface_to_Drefid = GT.face_reference_id(mesh,D)
#    Drefid_to_refDface = GT.reference_faces(mesh,D)
#    Drefid_to_lface_to_perm_to_coords = map(Drefid_to_refDface) do refDface
#        boundary = refDface |> GT.geometry |> GT.boundary
#        lface_to_nodes = GT.face_nodes(boundary,d)
#        node_to_coords = GT.node_coordinates(boundary)
#        lface_to_lrefid = GT.face_reference_id(boundary,d)
#        lrefid_to_lrefface = GT.reference_faces(boundary,d)
#        lrefid_to_perm_to_ids = map(GT.node_permutations,lrefid_to_lrefface)
#        map(1:GT.num_faces(boundary,d)) do lface
#            lrefid = lface_to_lrefid[lface]
#            nodes = lface_to_nodes[lface]
#            perm_to_ids = lrefid_to_perm_to_ids[lrefid]
#            map(perm_to_ids) do ids
#                coords = node_to_coords[nodes[ids]]
#                coords
#            end
#        end
#    end
#    sface_to_dface = domain |> GT.faces
#    dface_to_drefid = GT.face_reference_id(mesh,d)
#    drefid_to_refdface = GT.reference_faces(mesh,d)
#    drefid_to_funs = map(GT.shape_functions,drefid_to_refdface)
#    GT.quantity(prototype,domain) do index
#        sface = index.face
#        dict = index.dict
#        sface_to_tfaces_sym = get!(dict,sface_to_tfaces,gensym("sface_to_tfaces"))
#        sface_to_lfaces_sym = get!(dict,sface_to_lfaces,gensym("sface_to_lfaces"))
#        tface_to_Dface_sym = get!(dict,tface_to_Dface,gensym("tface_to_Dface"))
#        sface_to_dface_sym = get!(dict,sface_to_dface,gensym("sface_to_dface"))
#        Dface_to_lface_to_perm_sym = get!(dict,Dface_to_lface_to_perm,gensym("Dface_to_lface_to_perm"))
#        Dface_to_Drefid_sym = get!(dict,Dface_to_Drefid,gensym("Dface_to_Drefid"))
#        dface_to_drefid_sym = get!(dict,dface_to_drefid,gensym("dface_to_drefid"))
#        Drefid_to_lface_to_perm_to_coords_sym = get!(dict,Drefid_to_lface_to_perm_to_coords,gensym("Drefid_to_lface_to_perm_to_coords"))
#        drefid_to_funs_sym = get!(dict,drefid_to_funs,gensym("drefid_to_funs"))
#        face_around = 1
#        @term begin
#            boundary_face_function(
#                $sface,
#                $sface_to_tfaces_sym,
#                $sface_to_lfaces_sym,
#                $tface_to_Dface_sym,
#                $sface_to_dface_sym,
#                $Dface_to_lface_to_perm_sym,
#                $Dface_to_Drefid_sym,
#                $dface_to_drefid_sym,
#                $Drefid_to_lface_to_perm_to_coords_sym,
#                $drefid_to_funs_sym,
#                $face_around,
#           )
#        end
#    end
#end
#
#function domain_map(glue::BoundaryGlue,::ReferenceDomain,::PhysicalDomain)
#    error("Case not yet implemented")
#end
#
#function domain_map(glue::BoundaryGlue,::PhysicalDomain,::ReferenceDomain)
#    error("Case not yet implemented")
#end

#function align_field(a::AbstractQuantity,domain::AbstractDomain)
#    glue = GT.domain_glue(domain,GT.domain(a))
#    align_field(a,glue)
#end
#
#function align_field(a::AbstractQuantity{<:PMesh},domain::AbstractDomain{<:PMesh})
#    q = map(GT.align_field,partition(a),partition(domain))
#    prototype = map(GT.prototype,q) |> PartitionedArrays.getany
#    term = map(GT.term,q)
#    GT.quantity(term,prototype,domain)
#end
#
#function align_field(a::AbstractQuantity,glue::InteriorGlue)
#    domain = glue |> GT.domain
#    prototype = GT.prototype(a)
#    term_a = GT.term(a)
#    sface_to_tfaces, = GT.target_face(glue)
#    GT.quantity(prototype,domain) do index
#        sface = index.face
#        sface_to_tfaces_sym = get!(index.dict,sface_to_tfaces,gensym("sface_to_tfaces"))
#        tface = @term $sface_to_tfaces_sym[$sface][1]
#        index2 = replace_face(index,tface)
#        ai = term_a(index2)
#        ai
#    end
#end
#
#function align_field(a::AbstractQuantity,glue::CoboundaryGlue)
#    prototype = GT.prototype(a)
#    domain = glue |> GT.domain
#    term_a = GT.term(a)
#    sface_to_tfaces, sface_to_lfaces, = glue |> GT.target_face
#    GT.quantity(prototype,domain) do index
#        face_around = index.face_around
#        @assert face_around !== nothing
#        sface = index.face
#        dict = index.dict
#        sface_to_tfaces_sym = get!(dict,sface_to_tfaces,gensym("sface_to_tfaces"))
#        tface = @term $sface_to_tfaces_sym[$sface][$face_around]
#        index2 = replace_face(index,tface)
#        ai = term_a(index2)
#        ai
#    end
#end
#
#function align_field(a::AbstractQuantity,glue::BoundaryGlue)
#    prototype = GT.prototype(a)
#    domain = glue |> GT.domain
#    term_a = GT.term(a)
#    sface_to_tfaces, sface_to_lfaces, = glue |> GT.target_face
#    face_around = GT.face_around(glue.domain)
#    GT.quantity(prototype,domain) do index
#        sface = index.face
#        dict = index.dict
#        sface_to_tfaces_sym = get!(dict,sface_to_tfaces,gensym("sface_to_tfaces"))
#        tface = @term $sface_to_tfaces_sym[$sface][1]
#        index2 = replace_face(index,tface)
#        index3 = replace_face_around(index2,face_around)
#        ai = term_a(index3)
#        ai
#    end
#end

#function inverse_map(q::AbstractQuantity)
#    domain = q |> GT.domain
#    D = domain |> GT.num_dims
#    x0 = zero(SVector{D,Float64})
#    x = constant_quantity(x0)
#    args = (q,x)
#    prototype = GT.return_prototype(inv_map,(map(GT.prototype,args)...))
#    fs = map(GT.term,args)
#    GT.quantity(prototype,domain) do index
#        f_exprs = map(f->f(index),fs)
#        @term inv_map($(f_exprs[1]),$(f_exprs[2]))
#    end
#end

#function inverse_face_map(dom::AbstractDomain)
#    inverse_face_map(mesh(dom),num_dims(dom))
#end
#
#function inverse_face_map(mesh::AbstractMesh,d)
#    x0 = zero(SVector{val_parameter(d),Float64})
#    x = constant_quantity(x0)
#    phi = face_map(mesh,d)
#    call(inv_map,phi,x)
#end

function Base.:∘(a::AbstractQuantity,phi::AbstractQuantity)
    call(∘,a,phi)
end

function get_symbol!(index,val::typeof(Base.:∘),name="";prefix=index.data.prefix)
    :(Base.:∘)
end

#function Base.:∘(a::AbstractQuantity,phi::AbstractQuantity)
#    compose(a,phi)
#end
#
#function compose(a::AbstractQuantity,phi::AbstractQuantity)
#    glue = GT.domain_glue(GT.domain(phi),GT.domain(a))
#    compose(a,phi,glue)
#end

#function compose(a::AbstractQuantity{<:PMesh},phi::AbstractQuantity{<:PMesh})
#    q = map(GT.compose,partition(a),partition(phi))
#    term = map(GT.term,q)
#    prototype = map(GT.prototype,q) |> PartitionedArrays.getany
#    domain = GT.domain(phi)
#    GT.quantity(term,prototype,domain)
#end

#function compose(a::AbstractQuantity,phi::AbstractQuantity,glue::InteriorGlue)
#    @assert GT.domain(a) == GT.codomain(glue)
#    g = GT.prototype(a)
#    f = GT.prototype(phi)
#    prototype = x-> g(f(x))
#    domain = phi |> GT.domain
#    term_a = GT.term(a)
#    term_phi = GT.term(phi)
#    sface_to_tfaces, = GT.target_face(glue)
#    GT.quantity(prototype,domain) do index
#        sface = index.face
#        sface_to_tfaces_sym = get!(index.dict,sface_to_tfaces,gensym("sface_to_tfaces"))
#        tface = @term $sface_to_tfaces_sym[$sface][1]
#        index2 = replace_face(index,tface)
#        ai = term_a(index2)
#        phii = term_phi(index)
#        @term $ai∘$phii
#    end
#
#end
#
#function compose(a::AbstractQuantity,phi::AbstractQuantity,glue::CoboundaryGlue)
#    @assert GT.domain(a) == GT.codomain(glue)
#    g = GT.prototype(a)
#    f = GT.prototype(phi)
#    prototype = x-> g(f(x))
#    domain = phi |> GT.domain
#    term_a = GT.term(a)
#    term_phi = GT.term(phi)
#    sface_to_tfaces, sface_to_lfaces, = glue |> GT.target_face
#    face_around = GT.face_around(glue.domain)
#    GT.quantity(prototype,domain) do index
#        face_around = index.face_around
#        @assert face_around !== nothing
#        sface = index.face
#        sface_to_tfaces_sym = get!(index.dict,sface_to_tfaces,gensym("sface_to_tfaces"))
#        tface = @term $sface_to_tfaces_sym[$sface][$face_around]
#        index2 = replace_face(index,tface)
#        ai = term_a(index2)
#        phii = term_phi(index)
#        @term $ai∘$phii
#    end
#end
#
#function compose(a::AbstractQuantity,phi::AbstractQuantity,glue::BoundaryGlue)
#    @assert GT.domain(a) == GT.codomain(glue)
#    g = GT.prototype(a)
#    f = GT.prototype(phi)
#    prototype = x-> g(f(x))
#    domain = phi |> GT.domain
#    term_a = GT.term(a)
#    term_phi = GT.term(phi)
#    sface_to_tfaces, sface_to_lfaces, = glue |> GT.target_face
#    face_around = GT.face_around(glue.domain)
#    GT.quantity(prototype,domain) do index
#        sface = index.face
#        sface_to_tfaces_sym = get!(index.dict,sface_to_tfaces,gensym("sface_to_tfaces"))
#        tface = @term $sface_to_tfaces_sym[$sface][1]
#        index2 = replace_face(index,tface)
#        index3 = replace_face_around(index2,face_around)
#        ai = term_a(index3)
#        phii = term_phi(index)
#        @term $ai∘$phii
#    end
#end

#function unit_normal(domain::AbstractDomain)
#    error("not implemented yet")
#end
#
#function unit_normal(domain::AbstractDomain,codomain::AbstractDomain)
#    glue = GT.domain_glue(domain,codomain)
#    unit_normal(domain,codomain,glue)
#end
#
function map_unit_normal(φ_fun,ϕ_fun,n)
    q -> begin
        p = φ_fun(q)
        J = ForwardDiff.jacobian(ϕ_fun,p)
        Jt = transpose(J)
        pinvJt = transpose(inv(Jt*J)*Jt)
        v = pinvJt*n
        m = sqrt(inner(v,v))
        if m < eps()
            return zero(v)
        else
            return v/m
        end
    end
end

function unit_normal(mesh::AbstractMesh,d)
    D = num_dims(mesh)
    @assert d == D-1
    phiD = physical_map(mesh,D)
    phidinv = inverse_physical_map(mesh,d)
    phidD = reference_map(mesh,d,D)
    ctype_to_refface = GT.reference_faces(mesh,D)
    Drid_to_lface_to_n_data= map(ctype_to_refface) do refface
        boundary = refface |> GT.geometry |> GT.boundary
        boundary |> GT.outwards_normals # TODO also rename?
    end
    topo = topology(mesh)
    dface_to_Dfaces_data, dface_to_ldfaces_data = GT.face_incidence_ext(topo,d,D)
    n = quantity() do index
        dface = face_index(index,d)
        Dface = face_index(index,D)
        Dface_to_rid = get_symbol!(a.index,face_reference_id(mesh,D),"Dface_to_Drid")
        dface_to_Dfaces = get_symbol!(index(a),dface_to_Dfaces_data,"dface_to_Dfaces")
        dface_to_ldfaces = get_symbol!(index(a),dface_to_ldfaces_data,"dface_to_ldfaces")
        Drid_to_lface_to_n = get_symbol!(index(a),Drid_to_lface_to_n_data,"Drid_to_lface_to_n")
        expr = @term begin
            Drid = $Dface_to_Drid[$Dface]
            Dfaces = $dface_to_Dfaces[$dface]
            Dface_around = GalerkinToolkit.find_face_adound($Dface,Dfaces)
            ldface = $dface_to_ldfaces[$dface][Dface_around]
            Drid_to_lface_to_n[Drif][ldface]
        end
        p = Drid_to_lface_to_n_data |> eltype |> eltype |> zero
        expr_term([d,D],expr,p,index)
    end
    call(map_unit_normal,phidD∘phidinv,phiD,n)
end

#
#function unit_normal(domain::ReferenceDomain,codomain::PhysicalDomain,glue::BoundaryGlue)
#    Γref = domain
#    Ω = codomain
#    Ωref = GT.reference_domain(Ω)
#    D = GT.num_dims(Ω)
#    mesh = GT.mesh(Ω)
#    φ = domain_map(Γref,Ωref)
#    ϕ = GT.domain_map(Ωref,Ω)
#    sface_to_tfaces, sface_to_lfaces, = GT.target_face(glue)
#    tface_to_face = GT.faces(Ωref)
#    face_to_ctype = GT.face_reference_id(mesh,D)
#    ctype_to_refface = GT.reference_faces(mesh,D)
#    ctype_to_lface_to_n= map(ctype_to_refface) do refface
#        boundary = refface |> GT.geometry |> GT.boundary
#        boundary |> GT.outwards_normals # TODO also rename?
#    end
#    ϕ_term = GT.term(ϕ)
#    φ_term = GT.term(φ)
#    prototype = GT.prototype(φ)
#    GT.quantity(prototype,Γref) do index
#        sface = index.face
#        sface_to_tfaces_sym = get!(index.dict,sface_to_tfaces,gensym("sface_to_tfaces"))
#        sface_to_lfaces_sym = get!(index.dict,sface_to_lfaces,gensym("sface_to_lfaces"))
#        tface_to_face_sym = get!(index.dict,tface_to_face,gensym("tface_to_face"))
#        face_to_ctype_sym = get!(index.dict,face_to_ctype,gensym("face_to_ctype"))
#        ctype_to_lface_to_n_sym = get!(index.dict,ctype_to_lface_to_n,gensym("ctype_to_lface_to_n"))
#        tface = @term $sface_to_tfaces_sym[$sface][1]
#        index2 = replace_face(index,tface)
#        ϕ_fun = ϕ_term(index2)
#        φ_fun = φ_term(index)
#        @term begin
#            lface = $sface_to_lfaces_sym[$sface][1]
#            face = $tface_to_face_sym[$tface]
#            ctype = $face_to_ctype_sym[face]
#            lface_to_n = $ctype_to_lface_to_n_sym[ctype]
#            n = lface_to_n[lface]
#            boundary_unit_normal($φ_fun,$ϕ_fun,n)
#        end
#    end
#end
#
## TODO a lot of code duplication
#function unit_normal(domain::PhysicalDomain,codomain::PhysicalDomain,glue::BoundaryGlue)
#    Γ = domain
#    Γref = GT.physical_domain(Γ)
#    Ω = codomain
#    Ωref = GT.reference_domain(Ω)
#    ϕ = GT.domain_map(Ωref,Ω)
#    D = GT.num_dims(Ω)
#    mesh = GT.mesh(Ω)
#    ϕ = GT.domain_map(Ωref,Ω)
#    ϕinv = GT.inverse_map(ϕ)
#    sface_to_tfaces, sface_to_lfaces, = GT.target_face(glue)
#    tface_to_face = GT.faces(Ωref)
#    face_to_ctype = GT.face_reference_id(mesh,D)
#    ctype_to_refface = GT.reference_faces(mesh,D)
#    ctype_to_lface_to_n= map(ctype_to_refface) do refface
#        boundary = refface |> GT.geometry |> GT.boundary
#        boundary |> GT.outwards_normals # TODO also rename?
#    end
#    ϕinv_term = GT.term(ϕinv)
#    ϕ_term = GT.term(ϕ)
#    prototype = GT.prototype(ϕ)
#    GT.quantity(prototype,Γref) do index
#        sface = index.face
#        sface_to_tfaces_sym = get!(index.dict,sface_to_tfaces,gensym("sface_to_tfaces"))
#        sface_to_lfaces_sym = get!(index.dict,sface_to_lfaces,gensym("sface_to_lfaces"))
#        tface_to_face_sym = get!(index.dict,tface_to_face,gensym("tface_to_face"))
#        face_to_ctype_sym = get!(index.dict,face_to_ctype,gensym("face_to_ctype"))
#        ctype_to_lface_to_n_sym = get!(index.dict,ctype_to_lface_to_n,gensym("ctype_to_lface_to_n"))
#        tface = @term $sface_to_tfaces_sym[$sface][1]
#        index2 = replace_face(index,tface)
#        ϕinv_fun = ϕinv_term(index2)
#        ϕ_fun = ϕ_term(index2)
#        @term begin
#            lface = $sface_to_lfaces_sym[$sface][1]
#            face = $tface_to_face_sym[$tface]
#            ctype = $face_to_ctype_sym[face]
#            lface_to_n = $ctype_to_lface_to_n_sym[ctype]
#            n = lface_to_n[lface]
#            boundary_unit_normal($ϕinv_fun,$ϕ_fun,n)
#        end
#    end
#end
#
## TODO a lot of code duplication
#function unit_normal(domain::ReferenceDomain,codomain::PhysicalDomain,glue::CoboundaryGlue)
#    Γref = domain
#    Ω = codomain
#    Ωref = GT.reference_domain(Ω)
#    D = GT.num_dims(Ω)
#    mesh = GT.mesh(Ω)
#    φ = domain_map(Γref,Ωref)
#    ϕ = GT.domain_map(Ωref,Ω)
#    sface_to_tfaces, sface_to_lfaces, = GT.target_face(glue)
#    tface_to_face = GT.faces(Ωref)
#    face_to_ctype = GT.face_reference_id(mesh,D)
#    ctype_to_refface = GT.reference_faces(mesh,D)
#    ctype_to_lface_to_n= map(ctype_to_refface) do refface
#        boundary = refface |> GT.geometry |> GT.boundary
#        boundary |> GT.outwards_normals # TODO also rename?
#    end
#    ϕ_term = GT.term(ϕ)
#    φ_term = GT.term(φ)
#    prototype = GT.prototype(φ)
#    GT.quantity(prototype,Γref) do index
#        sface = index.face
#        face_around = index.face_around
#        @assert face_around !== nothing
#        sface_to_tfaces_sym = get!(index.dict,sface_to_tfaces,gensym("sface_to_tfaces"))
#        sface_to_lfaces_sym = get!(index.dict,sface_to_lfaces,gensym("sface_to_lfaces"))
#        tface_to_face_sym = get!(index.dict,tface_to_face,gensym("tface_to_face"))
#        face_to_ctype_sym = get!(index.dict,face_to_ctype,gensym("face_to_ctype"))
#        ctype_to_lface_to_n_sym = get!(index.dict,ctype_to_lface_to_n,gensym("ctype_to_lface_to_n"))
#        tface = @term $sface_to_tfaces_sym[$sface][$face_around]
#        index2 = replace_face(index,tface)
#        index3 = replace_face_around(index2,nothing)
#        ϕ_fun = ϕ_term(index3)
#        φ_fun = φ_term(index)
#        @term begin
#            lface = $sface_to_lfaces_sym[$sface][$face_around]
#            face = $tface_to_face_sym[$tface]
#            ctype = $face_to_ctype_sym[face]
#            lface_to_n = $ctype_to_lface_to_n_sym[ctype]
#            n = lface_to_n[lface]
#            boundary_unit_normal($φ_fun,$ϕ_fun,n)
#        end
#    end
#end
#
#function unit_normal(domain::PhysicalDomain,codomain::PhysicalDomain,glue::CoboundaryGlue)
#    Γ = domain
#    Γref = GT.physical_domain(Γ)
#    Ω = codomain
#    Ωref = GT.reference_domain(Ω)
#    ϕ = GT.domain_map(Ωref,Ω)
#    D = GT.num_dims(Ω)
#    mesh = GT.mesh(Ω)
#    ϕ = GT.domain_map(Ωref,Ω)
#    ϕinv = GT.inverse_map(ϕ)
#    sface_to_tfaces, sface_to_lfaces, = GT.target_face(glue)
#    tface_to_face = GT.faces(Ωref)
#    face_to_ctype = GT.face_reference_id(mesh,D)
#    ctype_to_refface = GT.reference_faces(mesh,D)
#    ctype_to_lface_to_n= map(ctype_to_refface) do refface
#        boundary = refface |> GT.geometry |> GT.boundary
#        boundary |> GT.outwards_normals # TODO also rename?
#    end
#    ϕinv_term = GT.term(ϕinv)
#    ϕ_term = GT.term(ϕ)
#    prototype = GT.prototype(ϕ)
#    GT.quantity(prototype,Γref) do index
#        sface = index.face
#        face_around = index.face_around
#        @assert face_around !== nothing
#        sface_to_tfaces_sym = get!(index.dict,sface_to_tfaces,gensym("sface_to_tfaces"))
#        sface_to_lfaces_sym = get!(index.dict,sface_to_lfaces,gensym("sface_to_lfaces"))
#        tface_to_face_sym = get!(index.dict,tface_to_face,gensym("tface_to_face"))
#        face_to_ctype_sym = get!(index.dict,face_to_ctype,gensym("face_to_ctype"))
#        ctype_to_lface_to_n_sym = get!(index.dict,ctype_to_lface_to_n,gensym("ctype_to_lface_to_n"))
#        tface = @term $sface_to_tfaces_sym[$sface][$face_around]
#        index2 = replace_face(index,tface)
#        index3 = replace_face_around(index2,nothing)
#        ϕinv_fun = ϕinv_term(index3)
#        ϕ_fun = ϕ_term(index3)
#        @term begin
#            lface = $sface_to_lfaces_sym[$sface][$face_around]
#            face = $tface_to_face_sym[$tface]
#            ctype = $face_to_ctype_sym[face]
#            lface_to_n = $ctype_to_lface_to_n_sym[ctype]
#            n = lface_to_n[lface]
#            boundary_unit_normal($ϕinv_fun,$ϕ_fun,n)
#        end
#    end
#end

function piecewiese_field(fields::AbstractQuantity...)
    PiecewiseField(fields)
end

struct PiecewiseField{A} <: AbstractType
    fields::A
end

function domain(u::PiecewiseField)
    domains = map(GT.domain,u.fields)
    PiecewiseDomain(domains)
end

function piecewiese_domain(domains::AbstractDomain...)
    PiecewiseDomain(domains)
end

struct PiecewiseDomain{A} <: AbstractType
    domains::A
end

# Operations

# Base

for op in (:+,:-,:sqrt,:abs,:abs2,:real,:imag,:conj,:transpose,:adjoint,:*,:/,:\,:^,:getindex)
    @eval begin
        function get_symbol!(index,val::typeof(Base.$op),name="";prefix=gensym)
            $( Expr(:quote,op) )
        end
    end
end

for op in (:+,:-,:sqrt,:abs,:abs2,:real,:imag,:conj,:transpose,:adjoint)
  @eval begin
      (Base.$op)(a::AbstractQuantity) = call(Base.$op,a)
  end
end

function Base.getindex(a::AbstractQuantity,b::AbstractQuantity)
    call(getindex,a,b)
end

function Base.getindex(a::AbstractQuantity,b::Number)
    Base.getindex(a,constant_quantity(b))
end

for op in (:+,:-,:*,:/,:\,:^)
  @eval begin
      (Base.$op)(a::AbstractQuantity,b::AbstractQuantity) = call(Base.$op,a,b)
      (Base.$op)(a::Number,b::AbstractQuantity) = call(Base.$op,GT.constant_quantity(a),b)
      (Base.$op)(a::AbstractQuantity,b::Number) = call(Base.$op,a,GT.constant_quantity(b))
  end
end

# LinearAlgebra

for op in (:inv,:det,:norm,:tr)
  @eval begin
      function get_symbol!(index,val::typeof(LinearAlgebra.$op),name="";prefix=gensym)
          $( Expr(:quote,op) )
      end
    (LinearAlgebra.$op)(a::AbstractQuantity) = call(LinearAlgebra.$op,a)
  end
end

for op in (:dot,:cross)
  @eval begin
      function get_symbol!(index,val::typeof(LinearAlgebra.$op),name="";prefix=gensym)
          $( Expr(:quote,op) )
      end
      (LinearAlgebra.$op)(a::AbstractQuantity,b::AbstractQuantity) = call(LinearAlgebra.$op,a,b)
      (LinearAlgebra.$op)(a::Number,b::AbstractQuantity) = call(LinearAlgebra.$op,GT.constant_quantity(a),b)
      (LinearAlgebra.$op)(a::AbstractQuantity,b::Number) = call(LinearAlgebra.$op,a,GT.constant_quantity(b))
  end
end

# ForwardDiff


#function call_function_symbol(g,args::AbstractQuantity...)
#    fs = map(GT.term,args)
#    domain = args |> first |> GT.domain
#    prototype = GT.return_prototype(g,(map(GT.prototype,args)...))
#    g_expr = nameof(g)
#    GT.quantity(prototype,domain) do index
#        f_exprs = map(f->f(index),fs)
#        :($g_expr($(f_exprs...)))
#    end
#end

#function call_function_symbol(g,args::AbstractQuantity{<:PMesh}...)
#    pargs = map(partition,args)
#
#    q = map(pargs...) do myargs...
#        call_function_symbol(g, myargs...)
#    end
#    domain = args |> first |> GT.domain
#    term = map(GT.term,q)
#    prototype = map(GT.prototype,q) |> PartitionedArrays.getany
#    GT.quantity(term,prototype,domain)
#end

for op in (:(ForwardDiff.gradient),:(ForwardDiff.jacobian),:(ForwardDiff.hessian))
  @eval begin
      function get_symbol!(index,val::typeof($op),name="";prefix=gensym)
          $( Expr(:quote,op) )
      end
      ($op)(a::AbstractQuantity,b::AbstractQuantity) = call($op,a,b)
  end
end


abstract type AbstractDomain{A} <: gk.AbstractType end
domain(a::AbstractDomain) = a
mesh(a::AbstractDomain) = a.mesh
mesh_id(a::AbstractDomain) = a.mesh_id
physical_names(a::AbstractDomain) = a.physical_names
face_dim(a::AbstractDomain) = gk.val_parameter(a.face_dim)
# TODO two functions for the same
num_dims(a::AbstractDomain) = face_dim(a)
is_reference_domain(a::AbstractDomain) = a.is_reference_domain |> gk.val_parameter
face_around(a::AbstractDomain) = a.face_around

function interior(mesh;
    mesh_id = objectid(mesh),
    physical_names=gk.physical_names(mesh,num_dims(mesh)),
    is_reference_domain = Val(false))
    D = num_dims(mesh)
    domain(mesh;face_dim=D,face_around=nothing,mesh_id,physical_names,is_reference_domain)
end

function skeleton(mesh;
    mesh_id = objectid(mesh),
    physical_names=gk.physical_names(mesh,num_dims(mesh)-1),
    is_reference_domain = Val(false))
    D = num_dims(mesh)
    domain(mesh;face_dim=D-1,face_around=nothing,mesh_id,physical_names,is_reference_domain)
end

function boundary(mesh::Union{AbstractFEMesh,PMesh};
    face_around=1,
    mesh_id = objectid(mesh),
    physical_names=gk.physical_names(mesh,num_dims(mesh)-1),
    is_reference_domain = Val(false))
    D = num_dims(mesh)
    domain(mesh;face_dim=D-1,face_around,mesh_id,physical_names,is_reference_domain)
end

function domain(mesh;
    mesh_id = objectid(mesh),
    face_dim = Val(gk.num_dims(mesh)),
    physical_names=gk.physical_names(mesh,face_dim),
    is_reference_domain = Val(false),
    face_around=nothing,
    )

    if val_parameter(is_reference_domain)
        ReferenceDomain(
                        mesh,
                        mesh_id,
                        physical_names,
                        Val(val_parameter(face_dim)),
                        face_around,
                       )
    else
        PhysicalDomain(
                        mesh,
                        mesh_id,
                        physical_names,
                        Val(val_parameter(face_dim)),
                        face_around,
                       )
    end

end

struct PhysicalDomain{A,B,C,D,E} <: AbstractDomain{A}
    mesh::A
    mesh_id::B
    physical_names::C
    face_dim::Val{D}
    face_around::E
end
is_reference_domain(a::PhysicalDomain) = false

struct ReferenceDomain{A,B,C,D,E} <: AbstractDomain{A}
    mesh::A
    mesh_id::B
    physical_names::C
    face_dim::Val{D}
    face_around::E
end
is_reference_domain(a::ReferenceDomain) = true

function replace_mesh(domain::AbstractDomain,mesh)
    face_dim = gk.face_dim(domain)
    mesh_id = gk.mesh_id(domain)
    physical_names = gk.physical_names(domain)
    is_reference_domain = gk.is_reference_domain(domain)
    face_around = gk.face_around(domain)
    gk.domain(mesh;face_dim,mesh_id,physical_names,is_reference_domain,face_around)
end

function replace_face_around(domain::AbstractDomain,face_around)
    mesh = gk.mesh(domain)
    face_dim = gk.face_dim(domain)
    mesh_id = gk.mesh_id(domain)
    physical_names = gk.physical_names(domain)
    is_reference_domain = gk.is_reference_domain(domain)
    gk.domain(mesh;face_dim,mesh_id,physical_names,is_reference_domain,face_around)
end

function PartitionedArrays.partition(domain::AbstractDomain{<:PMesh})
    pmesh = gk.mesh(domain)
    map(pmesh.mesh_partition) do mesh
        replace_mesh(domain,mesh)
    end
end

function Base.:(==)(a::AbstractDomain,b::AbstractDomain)
    flag = true
    # TODO check also that one mesh is not a sequential one and the other a parallel one
    flag = flag && (gk.mesh_id(a) == gk.mesh_id(b))
    flag = flag && (gk.physical_names(a) == gk.physical_names(b))
    flag = flag && (gk.face_dim(a) == gk.face_dim(b))
    flag = flag && (gk.is_reference_domain(a) == gk.is_reference_domain(b))
    flag
end

function reference_domain(domain::PhysicalDomain)
    mesh = gk.mesh(domain)
    face_dim = gk.face_dim(domain)
    mesh_id = gk.mesh_id(domain)
    physical_names = gk.physical_names(domain)
    is_reference_domain = Val(true)
    face_around = gk.face_around(domain)
    gk.domain(mesh;face_dim,mesh_id,physical_names,is_reference_domain,face_around)
end

function reference_domain(domain::ReferenceDomain)
    domain
end

function physical_domain(domain::ReferenceDomain)
    mesh = gk.mesh(domain)
    face_dim = gk.face_dim(domain)
    mesh_id = gk.mesh_id(domain)
    physical_names = gk.physical_names(domain)
    face_around = gk.face_around(domain)
    is_reference_domain = Val(false)
    gk.domain(mesh;face_dim,mesh_id,physical_names,is_reference_domain,face_around)
end

function physical_domain(domain::PhysicalDomain)
    domain
end

function faces(domain::AbstractDomain)
    mesh = domain |> gk.mesh
    D = gk.face_dim(domain)
    Dface_to_tag = zeros(Int,gk.num_faces(mesh,D))
    tag_to_name = gk.physical_names(domain)
    gk.classify_mesh_faces!(Dface_to_tag,mesh,D,tag_to_name)
    physical_Dfaces = findall(i->i!=0,Dface_to_tag)
    physical_Dfaces
end

function faces(domain::AbstractDomain{<:PMesh})
    map(gk.faces,partition(domain))
end

function num_faces(domain::AbstractDomain)
    length(faces(domain))
end

function num_faces(domain::AbstractDomain{<:PMesh})
    map(gk.num_faces,partition(domain))
end

abstract type AbstractDomainGlue{A} <: gk.AbstractType end
mesh(a::AbstractDomainGlue) = a.mesh
domain(a::AbstractDomainGlue) = a.domain
codomain(a::AbstractDomainGlue) = a.codomain

function PartitionedArrays.partition(a::AbstractDomainGlue{<:PMesh})
    if hasproperty(a,:face_around)
        map(partition(gk.domain(a)),partition(gk.codomain(a))) do dom,cod
            gk.domain_glue(dom,cod;a.face_around)
        end
    else
        map(gk.domain_glue,partition(gk.domain(a)),partition(gk.codomain(a)))
    end
end

function domain_glue(domain::AbstractDomain,codomain::AbstractDomain)
    msg = "Trying to combine domains on different meshes"
    @assert gk.mesh_id(domain) == gk.mesh_id(codomain) msg
    mesh = gk.mesh(domain)
    Ddom = gk.face_dim(domain)
    Dcod = gk.face_dim(codomain)
    face_around = gk.face_around(domain)
    if Ddom == Dcod
        InteriorGlue(mesh,domain,codomain)
    elseif Ddom < Dcod
        if face_around === nothing
            CoboundaryGlue(mesh,domain,codomain)
        else
            BoundaryGlue(mesh,domain,codomain,face_around)
        end
    else
        error("This case does not make sense")
    end
end

struct InteriorGlue{A,B,C} <: AbstractDomainGlue{A}
    mesh::A
    domain::B
    codomain::C
end

struct BoundaryGlue{A,B,C,D} <: AbstractDomainGlue{A}
    mesh::A
    domain::B
    codomain::C
    face_around::D
end

struct CoboundaryGlue{A,B,C} <: AbstractDomainGlue{A}
    mesh::A
    domain::B
    codomain::C
end

function target_face(glue::InteriorGlue)
    mesh = glue |> gk.domain |> gk.mesh
    domain = glue |> gk.domain
    codomain = glue |> gk.codomain
    d = domain |> gk.face_dim
    D = codomain |> gk.face_dim
    @assert d == D
    sface_to_dface = domain |> gk.faces
    tface_to_Dface = codomain |> gk.faces
    Dface_to_tface = zeros(Int32,gk.num_faces(mesh,D))
    tface_to_tface = LinearIndices(tface_to_Dface)
    Dface_to_tface[tface_to_Dface] = tface_to_tface
    sface_to_tface = Dface_to_tface[sface_to_dface]
end

function target_face(glue::CoboundaryGlue)
    mesh = glue |> gk.domain |> gk.mesh
    domain = glue |> gk.domain
    codomain = glue |> gk.codomain
    d = domain |> gk.face_dim
    D = codomain |> gk.face_dim
    @assert d < D
    sface_to_dface = domain |> gk.faces
    tface_to_Dface = codomain |> gk.faces
    Dface_to_tface = zeros(Int32,gk.num_faces(mesh,D))
    tface_to_tface = LinearIndices(tface_to_Dface)
    Dface_to_tface[tface_to_Dface] = tface_to_tface
    topo = gk.topology(mesh)
    dface_to_lfaces = gk.face_local_faces(topo,d,D)
    sface_to_lfaces = JaggedArray(view(dface_to_lfaces,sface_to_dface))
    dface_to_Dfaces = gk.face_incidence(topo,d,D)
    sface_to_Dfaces = JaggedArray(view(dface_to_Dfaces,sface_to_dface))
    data = sface_to_Dfaces.data
    f(Dface) = Dface_to_tface[Dface]
    data .= f.(data)
    sface_to_tfaces = sface_to_Dfaces
    sface_to_tfaces, sface_to_lfaces
end

function target_face(glue::BoundaryGlue)
    domain = replace_face_around(glue.domain,nothing)
    glue2 = domain_glue(domain,glue.codomain)
    sface_to_tfaces, sface_to_lfaces = target_face(glue2)
    face_around = glue.face_around
    sface_to_tface = map(tfaces->tfaces[face_around],sface_to_tfaces)
    sface_to_lface = map(tfaces->tfaces[face_around],sface_to_lfaces)
    sface_to_tface, sface_to_lface
end

abstract type AbstractQuantity{A} <: gk.AbstractType end
mesh(a::AbstractQuantity) = a.mesh
term(a::AbstractQuantity) = a.term
prototype(a::AbstractQuantity) = a.prototype
domain(a::AbstractQuantity) = a.domain
function PartitionedArrays.partition(a::AbstractQuantity)
    prototype = a |> gk.prototype
    map(gk.term(a),partition(gk.domain(a))) do term,domain
        gk.quantity(term,prototype,domain)
    end
end

function quantity(term,prototype,domain)
    mesh = gk.mesh(domain)
    Quantity(mesh,term,prototype,domain)
end

struct Quantity{A,B,C,D} <: AbstractQuantity{A}
    mesh::A
    term::B
    prototype::C
    domain::D
end

function constant_quantity(v,domain::AbstractDomain)
    gk.quantity(v,domain) do index
        v
    end
end

function constant_quantity(v,domain::AbstractDomain{<:PMesh})
    pmesh = gk.mesh(domain)
    term = map(pmesh.mesh_partition) do _
        index -> v
    end
    gk.quantity(term,v,domain)
end

function index(;
    face=nothing,
    local_face=nothing,
    face_around=nothing,
    point=nothing,
    field_per_dim =nothing,
    dof_per_dim=nothing,
    face_around_per_dim=nothing
    )
    Index(
          face,
          local_face,
          face_around,
          point,
          field_per_dim,
          dof_per_dim,
          face_around_per_dim
         )
end

struct Index{A,B,D,E,F,G,H}
    face::A
    local_face::B
    face_around::D
    point::E
    field_per_dim::F
    dof_per_dim::G
    face_around_per_dim::H
end

function replace_face(index::Index,face)
    Index(
          face,
          index.local_face,
          index.face_around,
          index.point,
          index.field_per_dim,
          index.dof_per_dim,
          index.face_around_per_dim
         )
end

function replace_local_face(index::Index,local_face)
    Index(
          index.face,
          local_face,
          index.face_around,
          index.point,
          index.field_per_dim,
          index.dof_per_dim,
          index.face_around_per_dim
         )
end

function replace_face_around(index::Index,face_around)
    Index(
          index.face,
          index.local_face,
          face_around,
          index.point,
          index.field_per_dim,
          index.dof_per_dim,
          index.face_around_per_dim
         )
end

function replace_point(index::Index,point)
    Index(
          index.face,
          index.local_face,
          index.face_around,
          point,
          index.field_per_dim,
          index.dof_per_dim,
          index.face_around_per_dim
         )
end

function replace_field_per_dim(index::Index,field_per_dim)
    Index(
          index.face,
          index.local_face,
          index.face_around,
          index.point,
          field_per_dim,
          index.dof_per_dim,
          index.face_around_per_dim
         )
end

function replace_dof_per_dim(index::Index,dof_per_dim)
    Index(
          index.face,
          index.local_face,
          index.face_around,
          index.point,
          index.field_per_dim,
          dof_per_dim,
          index.face_around_per_dim
         )
end

function return_prototype(f,args...)
    f(args...)
end

function call(f,args...)
    f(args...)
end

function call(g,args::AbstractQuantity...)
    fs = map(gk.term,args)
    domain = args |> first |> gk.domain
    #msg = "All quantities need to be defined on the same domain"
    #@assert all(dom->dom==domain,map(gk.domain,args)) msg
    # TODO check everything except reference/physical domain?
    # Maybe check reference/physical domain only when evaluating functions?
    prototype = gk.return_prototype(g,map(gk.prototype,args)...)
    gk.quantity(prototype,domain) do index
        g(map(f->f(index),fs)...)
    end
end

function call(g,args::AbstractQuantity{<:PMesh}...)
    pargs = map(partition,args)
    q = map(pargs...) do myargs...
        gk.call(g,myargs...)
    end
    domain = args |> first |> gk.domain
    term = map(gk.term,q)
    prototype = map(gk.prototype,q) |> PartitionedArrays.getany
    gk.quantity(term,prototype,domain)
end

function (f::AbstractQuantity)(x::AbstractQuantity)
    domain = gk.domain(x)
    codomain = gk.domain(f)
    flag = physical_domain(domain) == physical_domain(codomain)
    if flag
        call(call,f,x)
    else
        g = gk.align_field(f,domain)
        call(call,g,x)
    end
end

function analytical_field(f,dom::AbstractDomain)
    constant_quantity(f,dom)
end

function face_constant_field(data,dom::AbstractDomain)
    prototype = x->zero(eltype(data))
    gk.quantity(prototype,dom) do index
        face = index.face
        x->data[face]
    end
end

function domain_map(domain::AbstractDomain,codomain::AbstractDomain)
    glue = gk.domain_glue(domain,codomain)
    domain_map(glue,domain,codomain)
end

function domain_map(domain::AbstractDomain{<:PMesh},codomain::AbstractDomain{<:PMesh})
    q = map(gk.domain_map,partition(domain),partition(codomain))
    prototype = map(gk.prototype,q) |> PartitionedArrays.getany
    term = map(gk.term,q)
    gk.quantity(term,prototype,domain)
end

function domain_map(glue::InteriorGlue,::ReferenceDomain,::ReferenceDomain)
    domain = glue.domain
    prototype = identity
    term = identity
    gk.quantity(term,prototype,domain)
end

function domain_map(glue::InteriorGlue,::PhysicalDomain,::PhysicalDomain)
    domain = glue.domain
    prototype = identity
    term = identity
    gk.quantity(term,prototype,domain)
end

function domain_map(glue::InteriorGlue,::ReferenceDomain,::PhysicalDomain)
    domain = glue.domain
    mesh = domain |> gk.mesh
    d = domain |> gk.face_dim
    node_to_coords = gk.node_coordinates(mesh)
    sface_to_face = domain |> gk.faces
    face_to_nodes = gk.face_nodes(mesh,d)
    face_to_refid = gk.face_reference_id(mesh,d)
    refid_to_refface = gk.reference_faces(mesh,d)
    refid_to_funs = map(gk.shape_functions,refid_to_refface)
    T = eltype(gk.node_coordinates(mesh))
    x = zero(T)
    prototype = y->x
    gk.quantity(prototype,domain) do index
        sface = index.face
        face = sface_to_face[sface]
        refid = face_to_refid[face]
        funs = refid_to_funs[refid]
        nodes = face_to_nodes[face]
        coords = node_to_coords[nodes]
        q -> begin
            sum(1:length(coords)) do i
                x = coords[i]
                fun = funs[i]
                coeff = fun(q)
                coeff*x
            end
        end
    end
end

function domain_map(glue::InteriorGlue,::PhysicalDomain,::ReferenceDomain)
    error("Physical to reference map not implemented yet")
end

function domain_map(glue::CoboundaryGlue,::PhysicalDomain,::PhysicalDomain)
    error("Case not yet implemented")
end

function domain_map(glue::CoboundaryGlue,::ReferenceDomain,::ReferenceDomain)
    domain = glue.domain
    codomain = glue |> gk.codomain
    mesh = codomain |> gk.mesh
    D = codomain |> gk.face_dim
    Drefid_to_refDface = gk.reference_faces(mesh,D)
    refDface = first(Drefid_to_refDface)
    boundary = refDface |> gk.geometry |> gk.boundary
    node_to_coords = gk.node_coordinates(boundary)
    T = eltype(node_to_coords)
    x = zero(T)
    prototype = y -> [x,x]
    sface_to_tfaces, sface_to_lfaces = glue |> gk.target_face
    tface_to_Dface = codomain |> gk.faces
    d = domain |> gk.face_dim
    topo = mesh |> gk.topology
    Dface_to_lface_to_perm = gk.face_permutation_ids(topo,D,d)
    Dface_to_Drefid = gk.face_reference_id(mesh,D)
    Drefid_to_refDface = gk.reference_faces(mesh,D)
    Drefid_to_lface_to_perm_to_coords = map(Drefid_to_refDface) do refDface
        boundary = refDface |> gk.geometry |> gk.boundary
        lface_to_nodes = gk.face_nodes(boundary,d)
        node_to_coords = gk.node_coordinates(boundary)
        lface_to_lrefid = gk.face_reference_id(boundary,d)
        lrefid_to_lrefface = gk.reference_faces(boundary,d)
        lrefid_to_perm_to_ids = map(gk.node_permutations,lrefid_to_lrefface)
        map(1:gk.num_faces(boundary,d)) do lface
            lrefid = lface_to_lrefid[lface]
            nodes = lface_to_nodes[lface]
            perm_to_ids = lrefid_to_perm_to_ids[lrefid]
            map(perm_to_ids) do ids
                coords = node_to_coords[nodes[ids]]
                coords
            end
        end
    end
    sface_to_dface = domain |> gk.faces
    dface_to_drefid = gk.face_reference_id(mesh,d)
    drefid_to_refdface = gk.reference_faces(mesh,d)
    drefid_to_funs = map(gk.shape_functions,drefid_to_refdface)
    gk.quantity(prototype,domain) do index
        sface = index.face
        tfaces = sface_to_tfaces[sface]
        lfaces = sface_to_lfaces[sface]
        n_faces_around = length(lfaces)
        q -> begin
            map(1:n_faces_around) do face_around
                tface = tfaces[face_around]
                lface = lfaces[face_around]
                Dface = tface_to_Dface[tface]
                dface = sface_to_dface[sface]
                perm = Dface_to_lface_to_perm[Dface][lface]
                Drefid = Dface_to_Drefid[Dface]
                drefid = dface_to_drefid[dface]
                coords = Drefid_to_lface_to_perm_to_coords[Drefid][lface][perm]
                funs = drefid_to_funs[drefid]
                sum(1:length(coords)) do i
                    x = coords[i]
                    fun = funs[i]
                    coeff = fun(q)
                    coeff*x
                end
            end
        end
    end
end

function domain_map(glue::CoboundaryGlue,::ReferenceDomain,::PhysicalDomain)
    error("Case not yet implemented")
end

function domain_map(glue::CoboundaryGlue,::PhysicalDomain,::ReferenceDomain)
    error("Case not yet implemented")
end

function domain_map(glue::BoundaryGlue,::PhysicalDomain,::PhysicalDomain)
    error("Case not yet implemented")
end

function domain_map(glue::BoundaryGlue,::ReferenceDomain,::ReferenceDomain)
    domain = glue.domain
    codomain = glue |> gk.codomain
    mesh = codomain |> gk.mesh
    D = codomain |> gk.face_dim
    Drefid_to_refDface = gk.reference_faces(mesh,D)
    refDface = first(Drefid_to_refDface)
    boundary = refDface |> gk.geometry |> gk.boundary
    node_to_coords = gk.node_coordinates(boundary)
    T = eltype(node_to_coords)
    x = zero(T)
    prototype = y -> x
    sface_to_tface, sface_to_lfaces = glue |> gk.target_face
    tface_to_Dface = codomain |> gk.faces
    d = domain |> gk.face_dim
    topo = mesh |> gk.topology
    Dface_to_lface_to_perm = gk.face_permutation_ids(topo,D,d)
    Dface_to_Drefid = gk.face_reference_id(mesh,D)
    Drefid_to_refDface = gk.reference_faces(mesh,D)
    Drefid_to_lface_to_perm_to_coords = map(Drefid_to_refDface) do refDface
        boundary = refDface |> gk.geometry |> gk.boundary
        lface_to_nodes = gk.face_nodes(boundary,d)
        node_to_coords = gk.node_coordinates(boundary)
        lface_to_lrefid = gk.face_reference_id(boundary,d)
        lrefid_to_lrefface = gk.reference_faces(boundary,d)
        lrefid_to_perm_to_ids = map(gk.node_permutations,lrefid_to_lrefface)
        map(1:gk.num_faces(boundary,d)) do lface
            lrefid = lface_to_lrefid[lface]
            nodes = lface_to_nodes[lface]
            perm_to_ids = lrefid_to_perm_to_ids[lrefid]
            map(perm_to_ids) do ids
                coords = node_to_coords[nodes[ids]]
                coords
            end
        end
    end
    sface_to_dface = domain |> gk.faces
    dface_to_drefid = gk.face_reference_id(mesh,d)
    drefid_to_refdface = gk.reference_faces(mesh,d)
    drefid_to_funs = map(gk.shape_functions,drefid_to_refdface)
    gk.quantity(prototype,domain) do index
        sface = index.face
        tface = sface_to_tface[sface]
        lface = sface_to_lfaces[sface]
        Dface = tface_to_Dface[tface]
        dface = sface_to_dface[sface]
        perm = Dface_to_lface_to_perm[Dface][lface]
        Drefid = Dface_to_Drefid[Dface]
        drefid = dface_to_drefid[dface]
        coords = Drefid_to_lface_to_perm_to_coords[Drefid][lface][perm]
        funs = drefid_to_funs[drefid]
        q -> begin
            sum(1:length(coords)) do i
                x = coords[i]
                fun = funs[i]
                coeff = fun(q)
                coeff*x
            end
        end
    end
end

function domain_map(glue::BoundaryGlue,::ReferenceDomain,::PhysicalDomain)
    error("Case not yet implemented")
end

function domain_map(glue::BoundaryGlue,::PhysicalDomain,::ReferenceDomain)
    error("Case not yet implemented")
end

function align_field(a::AbstractQuantity,domain::AbstractDomain)
    glue = gk.domain_glue(domain,gk.domain(a))
    align_field(a,glue)
end

function align_field(a::AbstractQuantity{<:PMesh},domain::AbstractDomain{<:PMesh})
    q = map(gk.align_field,partition(a),partition(domain))
    prototype = map(gk.prototype,q) |> PartitionedArrays.getany
    term = map(gk.term,q)
    gk.quantity(term,prototype,domain)
end

function align_field(a::AbstractQuantity,glue::InteriorGlue)
    domain = glue |> gk.domain
    prototype = gk.prototype(a)
    term_a = gk.term(a)
    sface_to_tface = gk.target_face(glue)
    gk.quantity(prototype,domain) do index
        sface = index.face
        tface = sface_to_tface[sface]
        index2 = replace_face(index,tface)
        ai = term_a(index2)
        ai
    end
end

function align_field(a::AbstractQuantity,glue::CoboundaryGlue)
    pa = gk.prototype(a)
    prototype = x -> [pa(x),pa(x)]
    domain = glue |> gk.domain
    term_a = gk.term(a)
    sface_to_tfaces, sface_to_lfaces = glue |> gk.target_face
    gk.quantity(prototype,domain) do index
        sface = index.face
        tfaces = sface_to_tfaces[sface]
        lfaces = sface_to_lfaces[sface]
        n_faces_around = length(tfaces)
        x -> begin
            # TODO This should be a tuple
            map(1:n_faces_around) do face_around
                tface = sface_to_tfaces[sface][face_around]
                lface = sface_to_lfaces[sface][face_around]
                index2 = replace_face(index,tface)
                index3 = replace_face_around(index2,face_around)
                ai = term_a(index3)
                ai(x)
            end
        end
    end
end

function align_field(a::AbstractQuantity,glue::BoundaryGlue)
    prototype = gk.prototype(a)
    domain = glue |> gk.domain
    term_a = gk.term(a)
    sface_to_tface, sface_to_lface = glue |> gk.target_face
    face_around = glue.face_around
    gk.quantity(prototype,domain) do index
        sface = index.face
        tface = sface_to_tface[sface]
        lface = sface_to_lface[sface]
        index2 = replace_face(index,tface)
        index3 = replace_face_around(index2,face_around)
        ai = term_a(index3)
        ai
    end
end

function inverse_map(f)
    function invf(fx)
        x = zero(fx)
        tol = 1.0e-12
        niters = 100
        for _ in 1:niters
            J = ForwardDiff.jacobian(f,x)
            dx = J\(fx-f(x))
            x += dx
            if norm(dx) < tol
                return x
            end
        end
        error("Max iterations reached")
        x
    end
end

function return_prototype(f::typeof(inverse_map),g)
    g
end

function inverse_map(q::AbstractQuantity)

    gk.call(inverse_map,q)
end

function Base.:∘(a::AbstractQuantity,phi::AbstractQuantity)
    compose(a,phi)
end

function compose(a::AbstractQuantity,phi::AbstractQuantity)
    glue = gk.domain_glue(gk.domain(phi),gk.domain(a))
    compose(a,phi,glue)
end

function compose(a::AbstractQuantity{<:PMesh},phi::AbstractQuantity{<:PMesh})
    q = map(gk.compose,partition(a),partition(phi))
    term = map(gk.term,q)
    prototype = map(gk.prototype,q) |> PartitionedArrays.getany
    domain = gk.domain(phi)
    gk.quantity(term,prototype,domain)
end

function compose(a::AbstractQuantity,phi::AbstractQuantity,glue::InteriorGlue)
    @assert gk.domain(a) == gk.codomain(glue)
    g = gk.prototype(a)
    f = gk.prototype(phi)
    prototype = x-> g(f(x))
    domain = phi |> gk.domain
    term_a = gk.term(a)
    term_phi = gk.term(phi)
    sface_to_tface = gk.target_face(glue)
    gk.quantity(prototype,domain) do index
        sface = index.face
        tface = sface_to_tface[sface]
        index2 = replace_face(index,tface)
        ai = term_a(index2)
        phii = term_phi(index)
        x-> ai(phii(x))
    end

end

function compose(a::AbstractQuantity,phi::AbstractQuantity,glue::CoboundaryGlue)
    @assert gk.domain(a) == gk.codomain(glue)
    g = gk.prototype(a)
    f = gk.prototype(phi)
    prototype = x-> map(g,f(x))
    domain = phi |> gk.domain
    term_a = gk.term(a)
    term_phi = gk.term(phi)
    sface_to_tfaces, sface_to_lfaces = glue |> gk.target_face
    gk.quantity(prototype,domain) do index
        sface = index.face
        tfaces = sface_to_tfaces[sface]
        lfaces = sface_to_lfaces[sface]
        n_faces_around = length(tfaces)
        phii = term_phi(index)
        x -> begin
            ys = phii(x)
            # TODO This should be a tuple
            map(1:n_faces_around) do face_around
                tface = sface_to_tfaces[sface][face_around]
                lface = sface_to_lfaces[sface][face_around]
                index2 = replace_face(index,tface)
                index3 = replace_face_around(index2,face_around)
                ai = term_a(index3)
                y = ys[face_around]
                ai(y)
            end
        end
    end
end

function compose(a::AbstractQuantity,phi::AbstractQuantity,glue::BoundaryGlue)
    @assert gk.domain(a) == gk.codomain(glue)
    g = gk.prototype(a)
    f = gk.prototype(phi)
    prototype = x-> g(f(x))
    domain = phi |> gk.domain
    term_a = gk.term(a)
    term_phi = gk.term(phi)
    sface_to_tface, sface_to_lface = glue |> gk.target_face
    face_around = glue.face_around
    gk.quantity(prototype,domain) do index
        sface = index.face
        tface = sface_to_tface[sface]
        lface = sface_to_lface[sface]
        index2 = replace_face(index,tface)
        index3 = replace_face_around(index2,face_around)
        ai = term_a(index3)
        phii = term_phi(index)
        x -> ai(phii(x))
    end
end

function plot(domain::AbstractDomain;kwargs...)
    mesh = gk.mesh(domain)
    d = gk.face_dim(domain)
    domface_to_face = gk.faces(domain)
    vismesh = gk.visualization_mesh(mesh,d,domface_to_face;kwargs...)
    node_data = Dict{String,Any}()
    face_data = Dict{String,Any}()
    Plot(mesh,domain,vismesh,node_data,face_data)
end

function plot(domain::AbstractDomain{<:PMesh};kwargs...)
    mesh = gk.mesh(domain)
    args = map(partition(domain)) do mydom
        plt = plot(mydom;kwargs...)
        (plt.visualization_mesh, plt.node_data, plt.face_data)
    end |> tuple_of_arrays
    Plot(mesh,domain,args...)
end

struct Plot{A,B,C,D,E}
    mesh::A
    domain::B
    visualization_mesh::C
    node_data::D
    face_data::E
end
domain(plt::Plot) = plt.domain
visualization_mesh(plt::Plot) = plt.visualization_mesh
function PartitionedArrays.partition(plt::Plot{<:PMesh})
    map(Plot,
        partition(plt.mesh),partition(plt.domain),
        plt.visualization_mesh,plt.node_data,plt.face_data)
end

function reference_coordinates(plt::Plot)
    domain = gk.reference_domain(plt.domain)
    d = gk.face_dim(domain)
    domface_to_face = gk.faces(domain)
    mesh = gk.mesh(domain)
    vmesh, vglue = gk.visualization_mesh(plt)
    refid_to_snode_to_coords = vglue.reference_coordinates
    d = gk.num_dims(vmesh)
    face_to_refid = gk.face_reference_id(mesh,d)
    prototype = first(first(refid_to_snode_to_coords))
    gk.quantity(prototype,domain) do index
        domface = index.face[1]
        point = index.point
        face = domface_to_face[domface]
        refid = face_to_refid[face]
        refid_to_snode_to_coords[refid][point]
    end
end

function reference_coordinates(plt::Plot{<:PMesh})
    q = map(gk.reference_coordinates,partition(plt))
    term = map(gk.term,q)
    prototype = map(gk.prototype,q) |> PartitionedArrays.getany
    gk.quantity(term,prototype,plt.domain)
end

function coordinates(plt::Plot)
    domain = plt |> gk.domain
    gk.coordinates(plt,domain)
end

function coordinates(plt::Plot,::ReferenceDomain)
    gk.reference_coordinates(plt)
end

function coordinates(plt::Plot,::PhysicalDomain)
    domain_phys = plt |> gk.domain
    domain_ref = domain_phys |> reference_domain
    phi = gk.domain_map(domain_ref,domain_phys)
    q = gk.reference_coordinates(plt)
    phi(q)
end

function plot!(field,plt::Plot;label)
    plot!(plt,field;label)
end

function plot!(plt::Plot,field;label)
    q = gk.coordinates(plt)
    f_q = field(q)
    term = gk.term(f_q)
    T = typeof(gk.prototype(f_q))
    plot_impl!(plt,term,label,T)
end

function plot_impl!(plt,term,label,::Type{T}) where T
    vmesh,vglue = plt.visualization_mesh
    nnodes = gk.num_nodes(vmesh)
    data = zeros(T,nnodes)
    face_to_nodes = vglue.face_fine_nodes
    for face in 1:length(face_to_nodes)
        nodes = face_to_nodes[face]
        for point in 1:length(nodes)
            index = gk.index(;face,point)
            v = term(index)
            data[nodes[point]] = v
        end
    end
    plt.node_data[label] = data
    plt
end

function plot!(plt::Plot{<:PMesh},field;label)
    q = gk.coordinates(plt)
    f_q = field(q)
    term = gk.term(f_q)
    T = typeof(gk.prototype(f_q))
    map(partition(plt),term) do myplt, myterm
        plot_impl!(myplt,myterm,label,T)
    end
    plt
end

function vtk_plot(f,filename,args...;kwargs...)
    plt = gk.plot(args...;kwargs...)
    vtk_plot_impl(f,filename,plt)
end

function vtk_plot_impl(f,filename,plt::Plot)
    function translate(v)
        v
    end
    function translate(v::AbstractVector{<:SVector{2}})
        z = zero(eltype(eltype(v)))
        map(vi->SVector((vi...,z)),v)
    end
    vmesh,_ = plt.visualization_mesh
    d = gk.face_dim(plt.domain)
    r = f(plt)
    vtk_grid(filename,gk.vtk_args(vmesh,d)...) do vtk
        for (k,v) in plt.node_data
            vtk[k,WriteVTK.VTKPointData()] = translate(v)
        end
        for (k,v) in plt.face_data
            vtk[k,WriteVTK.VTKPointData()] = translate(v)
        end
        r
    end
end

function vtk_plot_impl(f,filename,pplt::Plot{<:PMesh})
    r = f(pplt)
    pmesh = pplt.mesh
    d = gk.face_dim(pplt.domain)
    parts = linear_indices(pmesh.mesh_partition)
    nparts = length(parts)
    map(partition(pplt),pmesh.face_partition[d+1],parts) do plt,myfaces,part
        vmesh,vglue = plt.visualization_mesh
        vcell_to_islocal =Int.(local_to_owner(myfaces) .== part)[vglue.parent_face]
        vcell_to_owner =local_to_owner(myfaces)[vglue.parent_face]
        pvtk_grid(filename,gk.vtk_args(vmesh,d)...;part,nparts) do vtk
            vtk["__PART__",WriteVTK.VTKCellData()] = fill(part,num_faces(vmesh,d))
            vtk["__LOCAL__",WriteVTK.VTKCellData()] = vcell_to_islocal
            vtk["__OWNER__",WriteVTK.VTKCellData()] = vcell_to_owner
            for (k,v) in plt.node_data
                vtk[k,WriteVTK.VTKPointData()] = v
            end
            for (k,v) in plt.face_data
                vtk[k,WriteVTK.VTKPointData()] = v
            end
        end
    end
    r
end

function unit_normal(domain::AbstractDomain)
    error("not implemented yet")
end

function unit_normal(domain::AbstractDomain,codomain::AbstractDomain)
    Γref = domain
    Ω = codomain
    Ωref = gk.reference_domain(Ω)
    D = gk.num_dims(Ω)
    mesh = gk.mesh(Ω)
    φ = domain_map(Γref,Ωref)
    glue = gk.domain_glue(Γref,Ωref)
    ϕ = gk.domain_map(Ωref,Ω)
    sface_to_tface, sface_to_lface = gk.target_face(glue)
    tface_to_face = gk.faces(Ωref)
    face_to_ctype = gk.face_reference_id(mesh,D)
    ctype_to_refface = gk.reference_faces(mesh,D)
    ctype_to_lface_to_n= map(ctype_to_refface) do refface
        boundary = refface |> gk.geometry |> gk.boundary
        boundary |> gk.outwards_normals # TODO also rename?
    end
    ϕ_term = gk.term(ϕ)
    φ_term = gk.term(φ)
    prototype = gk.prototype(ϕ)
    gk.quantity(prototype,Γref) do index
        sface = index.face
        tface = sface_to_tface[sface]
        lface = sface_to_lface[sface]
        face = tface_to_face[tface]
        ctype = face_to_ctype[face]
        lface_to_n = ctype_to_lface_to_n[ctype]
        n = lface_to_n[lface]
        index2 = replace_face(index,tface)
        ϕ_fun = ϕ_term(index2)
        φ_fun = φ_term(index)
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
end

function piecewiese_field(fields::AbstractQuantity...)
    PiecewiseField(fields)
end

struct PiecewiseField{A}
    fields::A
end

function domain(u::PiecewiseField)
    domains = map(gk.domain,u.fields)
    PiecewiseDomain(domains)
end

function piecewiese_domain(domains::AbstractDomain...)
    PiecewiseDomain(domains)
end

struct PiecewiseDomain{A}
    domains::A
end

# Operations

# Base

function Base.getindex(a::AbstractQuantity,i::Integer...)
    call(b->b[i...],a)
end

for op in (:+,:-,:sqrt,:abs,:abs2,:real,:imag,:conj,:transpose,:adjoint)
  @eval begin
    (Base.$op)(a::AbstractQuantity) = call(Base.$op,a)
  end
end

for op in (:+,:-,:*,:/,:\)
  @eval begin
      (Base.$op)(a::AbstractQuantity,b::AbstractQuantity) = call(Base.$op,a,b)
      (Base.$op)(a::Number,b::AbstractQuantity) = call(Base.$op,gk.constant_quantity(a,gk.domain(b)),b)
      (Base.$op)(a::AbstractQuantity,b::Number) = call(Base.$op,a,gk.constant_quantity(b,domain(a)))
  end
end

# LinearAlgebra

for op in (:inv,:det)
  @eval begin
    (LinearAlgebra.$op)(a::AbstractQuantity) = call(LinearAlgebra.$op,a)
  end
end

for op in (:dot,:cross)
  @eval begin
      (LinearAlgebra.$op)(a::AbstractQuantity,b::AbstractQuantity) = call(LinearAlgebra.$op,a,b)
      (LinearAlgebra.$op)(a::Number,b::AbstractQuantity) = call(LinearAlgebra.$op,gk.constant_quantity(a,gk.domain(b)),b)
      (LinearAlgebra.$op)(a::AbstractQuantity,b::Number) = call(LinearAlgebra.$op,a,gk.constant_quantity(b,domain(a)))
  end
end

# ForwardDiff

for op in (:gradient,:jacobian,:hessian)
  @eval begin
      (ForwardDiff.$op)(a::AbstractQuantity,b::AbstractQuantity) = call(ForwardDiff.$op,a,b)
      (ForwardDiff.$op)(a::Number,b::AbstractQuantity) = call(ForwardDiff.$op,gk.constant_quantity(a,gk.domain(b)),b)
      (ForwardDiff.$op)(a::AbstractQuantity,b::Number) = call(ForwardDiff.$op,a,gk.constant_quantity(b,domain(a)))
  end
end

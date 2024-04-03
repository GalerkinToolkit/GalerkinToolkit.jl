
abstract type AbstractDomain <: gk.AbstractType end
domain(a::AbstractDomain) = a
mesh(a::AbstractDomain) = a.mesh
mesh_id(a::AbstractDomain) = a.mesh_id
physical_names(a::AbstractDomain) = a.physical_names
face_dim(a::AbstractDomain) = gk.val_parameter(a.face_dim)
is_reference_domain(a::AbstractDomain) = a.is_reference_domain |> gk.val_parameter

function domain(mesh;
    mesh_id = objectid(mesh),
    face_dim = gk.num_dims(mesh),
    physical_names=gk.physical_names(mesh,face_dim),
    is_reference_domain = Val(false),
    )

    Domain(
           mesh,
           mesh_id,
           physical_names,
           face_dim,
           is_reference_domain,
          )
end

struct Domain{A,B,C,D,E} <: AbstractDomain
    mesh::A
    mesh_id::B
    physical_names::C
    face_dim::D
    is_reference_domain::E
end

abstract type AbstractDomainStyle <: gk.AbstractType end
is_reference_domain(a::AbstractDomainStyle) = a.is_reference_domain |> gk.val_parameter

struct GlobalDomain{A} <: AbstractDomainStyle
    is_reference_domain::Val{A}
end

struct LocalDomain{A} <: AbstractDomainStyle
    is_reference_domain::Val{A}
end

function domain_style(domain::AbstractDomain)
    flag = gk.is_reference_domain(domain)
    domain_style(domain,flag)
end

function domain_style(domain::AbstractDomain,is_reference_domain::Bool)
    GlobalDomain(Val(is_reference_domain))
end

function domain_style(domain::AbstractDomain,is_reference_domain::Tuple{Bool,Bool})
    LocalDomain(Val(is_reference_domain))
end

function Base.:(==)(a::AbstractDomain,b::AbstractDomain)
    flag = true
    flag = flag && (gk.mesh_id(a) == gk.mesh_id(b))
    flag = flag && (gk.physical_names(a) == gk.physical_names(b))
    flag = flag && (gk.face_dim(a) == gk.face_dim(b))
    flag = flag && (gk.is_reference_domain(a) == gk.is_reference_domain(b))
    flag
end

function reference_domain(domain::AbstractDomain)
    reference_domain(domain,domain |> gk.domain_style)
end

function reference_domain(domain::AbstractDomain,::GlobalDomain{false})
    Domain(
           domain |> gk.mesh,
           domain |> gk.mesh_id,
           domain |> gk.physical_names,
           domain |> gk.face_dim,
           Val(true),
          )
end

function reference_domain(domain::AbstractDomain,::GlobalDomain{true})
    domain
end

function faces(domain::AbstractDomain)
    faces(domain,gk.domain_style(domain))
end

function faces(domain::AbstractDomain,::GlobalDomain)
    D = gk.face_dim(domain)
    mesh = gk.mesh(domain)
    Dface_to_tag = zeros(Int,gk.num_faces(mesh,D))
    tag_to_name = gk.physical_names(domain)
    gk.classify_mesh_faces!(Dface_to_tag,mesh,D,tag_to_name)
    physical_Dfaces = findall(i->i!=0,Dface_to_tag)
    physical_Dfaces
end

function domain_glue(domain,codomain)
    msg = "Trying to combine domains on different meshes"
    @assert gk.mesh_id(domain) == gk.mesh_id(codomain) msg
    DomainGlue(domain,codomain)
end

struct DomainGlue{A,B}
    domain::A
    codomain::B
end
domain(a::DomainGlue) = a.domain
codomain(a::DomainGlue) = a.codomain

abstract type AbstractDomainGlueStyle <: gk.AbstractType end
domain_style(a::AbstractDomainGlueStyle) = a.domain_style
codomain_style(a::AbstractDomainGlueStyle) = a.codomain_style

struct InteriorGlue{A,B} <: AbstractDomainGlueStyle
    domain_style::GlobalDomain{A}
    codomain_style::GlobalDomain{B}
end

struct BoundaryGlue{A,B} <: AbstractDomainGlueStyle
    domain_style::LocalDomain{A}
    codomain_style::GlobalDomain{B}
end

struct CoboundaryGlue{A,B} <: AbstractDomainGlueStyle
    domain_style::GlobalDomain{A}
    codomain_style::GlobalDomain{B}
end

function domain_glue_style(glue::DomainGlue)
    d1 = glue |> gk.domain |> gk.face_dim |> gk.val_parameter
    d2 = glue |> gk.codomain |> gk.face_dim |> gk.val_parameter
    domain_glue_style(glue,d1,d2)
end

function domain_glue_style(glue,d1::Integer,d2::Integer)
    domain_style = glue |> gk.domain |> gk.domain_style
    codomain_style = glue |> gk.codomain |> gk.domain_style
    if d1 == d2
        InteriorGlue(domain_style,codomain_style)
    elseif d1 < d2
        CoboundaryGlue(domain_style,codomain_style)
    else
        error("This case does not make sense")
    end
end

function domain_glue_style(glue,d1::Tuple{Integer,Integer},d2::Integer)
    domain_style = glue |> gk.domain |> gk.domain_style
    codomain_style = glue |> gk.codomain |> gk.domain_style
    d1g,d1l = d1
    if d1l == d2
        BoundaryGlue(domain_style,codomain_style)
    elseif d1l < d2
        error("Case not supported yet")
    else
        error("This case does not make sense")
    end
end

function target_face(glue::DomainGlue)
    glue_style = gk.domain_glue_style(glue)
    target_face(glue,glue_style)
end

function target_index_term(glue::DomainGlue)
    glue_style = gk.domain_glue_style(glue)
    target_index_term(glue,glue_style)
end

function target_face(glue::DomainGlue,style::InteriorGlue)
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

function target_index_term(glue::DomainGlue,::InteriorGlue)
    sface_to_tface = glue |> gk.target_face
    index -> begin
        sface = index.face
        tface = sface_to_tface[sface]
        replace_face(index,tface)
    end
end

function target_face(glue::DomainGlue,style::CoboundaryGlue)
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

function target_index_term(glue::DomainGlue,::CoboundaryGlue)
    sface_to_tfaces, sface_to_lfaces = glue |> gk.target_face
    index -> begin
        face_around = index.face_around
        @assert face_around !== nothing
        sface = index.face
        tface = sface_to_tfaces[sface][face_around]
        lface = sface_to_lfaces[sface][face_around]
        index2 = replace_face(index,tface)
        replace_local_face(index2,lface)
    end
end

abstract type AbstractQuantity <: gk.AbstractType end
term(a::AbstractQuantity) = a.term
prototype(a::AbstractQuantity) = a.prototype
domain(a::AbstractQuantity) = a.domain

quantity(term,prototype,domain) = Quantity(term,prototype,domain)
struct Quantity{T,B,C} <: AbstractQuantity
    term::T
    prototype::B
    domain::C
end

function index(;
    face=nothing,
    local_face=nothing,
    face_around=nothing,
    point=nothing,
    field_per_dim=nothing,
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

struct Index{A,B,C,D,E,F,G}
    face::A
    local_face::B
    face_around::C
    point::D
    field_per_dim::E
    dof_per_dim::F
    face_around_per_dim::G
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

function return_prototype(f,args...)
    f(args...)
end

function call(f,args...)
    f(args...)
end

function call(g,args::AbstractQuantity...)
    fs = map(gk.term,args)
    domain = args |> first |> gk.domain
    msg = "All quantities need to be defined on the same domain"
    @assert all(dom->dom==domain,map(gk.domain,args)) msg
    prototype = gk.return_prototype(g,map(gk.prototype,args)...)
    gk.quantity(prototype,domain) do index
        g(map(f->f(index),fs)...)
    end
end

function (f::AbstractQuantity)(x::AbstractQuantity)
    call(call,f,x)
end

function analytical_field(f,dom)
    AnalyticalField(f,dom)
end

struct AnalyticalField{A,B} <: AbstractQuantity
    f::A
    domain::B
end
term(a::AnalyticalField) = index -> a.f
prototype(a::AnalyticalField) = a.f

function domain_map(domain,codomain;face_around=nothing)
    glue = gk.domain_glue(domain,codomain)
    glue_style = gk.domain_glue_style(glue)
    domain_map(glue,glue_style,face_around)
end

function domain_map(glue::DomainGlue,glue_style,face_around)
    @assert face_around === nothing
    DomainMap(glue,face_around)
end

function domain_map(glue::DomainGlue,::CoboundaryGlue,face_around::Integer)
    DomainMap(glue,face_around)
end

function domain_map(glue::DomainGlue,::CoboundaryGlue,::Nothing)
    sface_to_tfaces, = gk.target_face(glue)
    msg = "Not all faces are interior faces!"
    @assert all(tfaces->length(tfaces)==2,sface_to_tfaces) msg
    face_around_plus = 1
    face_around_minus = 2
    phi_plus = DomainMap(glue,face_around_plus)
    phi_minus = DomainMap(glue,face_around_minus)
    (phi_plus,phi_minus)
end

struct DomainMap{A,B} <: AbstractQuantity
    domain_glue::A
    face_around::B
end
domain_glue(a::DomainMap) = a.domain_glue
domain(a::DomainMap) = a |> gk.domain_glue |> gk.domain
codomain(a::DomainMap) = a |> gk.domain_glue |> gk.codomain

function prototype(a::DomainMap)
    glue_style = a |> gk.domain_glue |> gk.domain_glue_style
    prototype(a,glue_style)
end

function term(a::DomainMap)
    glue_style = a |> gk.domain_glue |> gk.domain_glue_style
    term(a,glue_style)
end

function Base.:∘(a::AbstractQuantity,phi::DomainMap)
    @assert gk.domain(a) == gk.codomain(phi)
    g = gk.prototype(a)
    f = gk.prototype(phi)
    prototype = x-> g(f(x))
    term_a = gk.term(a)
    term_phi = gk.term(phi)
    glue = domain_glue(phi)
    target_index_term = gk.target_index_term(glue)
    face_around = phi.face_around
    domain = phi |> gk.domain
    gk.quantity(prototype,domain) do index
        index2 = replace_face_around(index,face_around)
        index3 = target_index_term(index2)
        ai = term_a(index3)
        phii = term_phi(index)
        x-> ai(phii(x))
    end
end

function prototype(a::DomainMap,::InteriorGlue{true,true})
    identity
end

function term(a::DomainMap,::InteriorGlue{true,true})
    index -> identity
end

function prototype(a::DomainMap,::InteriorGlue{false,false})
    identity
end

function term(a::DomainMap,::InteriorGlue{false,false})
    index -> identity
end

function prototype(a::DomainMap,::InteriorGlue{true,false})
    glue = a |> gk.domain_glue
    domain = glue |> gk.domain
    mesh = domain |> gk.mesh
    T = eltype(gk.node_coordinates(mesh))
    x = zero(T)
    y->x
end
function term(a::DomainMap,::InteriorGlue{true,false})
    glue = a |> gk.domain_glue
    domain = glue |> gk.domain
    mesh = domain |> gk.mesh
    d = domain |> gk.face_dim
    node_to_coords = gk.node_coordinates(mesh)
    sface_to_face = domain |> gk.faces
    face_to_nodes = gk.face_nodes(mesh,d)
    face_to_refid = gk.face_reference_id(mesh,d)
    refid_to_refface = gk.reference_faces(mesh,d)
    refid_to_funs = map(gk.shape_functions,refid_to_refface)
    index -> begin
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

function prototype(a::DomainMap,::InteriorGlue{false,true})
    error("Physical to reference map not implemented yet")
end

function term(a::DomainMap,::InteriorGlue{false,true})
    error("Physical to reference map not implemented yet")
end

function prototype(phi::DomainMap,::CoboundaryGlue{true,true})
    glue = phi |> gk.domain_glue
    codomain = glue |> gk.codomain
    mesh = codomain |> gk.mesh
    D = codomain |> gk.face_dim
    Drefid_to_refDface = gk.reference_faces(mesh,D)
    refDface = first(Drefid_to_refDface)
    boundary = refDface |> gk.geometry |> gk.boundary
    node_to_coords = gk.node_coordinates(boundary)
    T = eltype(node_to_coords)
    x = zero(T)
    y -> x
end

function term(phi::DomainMap,::CoboundaryGlue{true,true})
    glue = phi |> gk.domain_glue
    domain = glue |> gk.domain
    codomain = glue |> gk.codomain
    sface_to_tfaces, sface_to_lfaces = glue |> gk.target_face
    tface_to_Dface = codomain |> gk.faces
    d = domain |> gk.face_dim
    D = codomain |> gk.face_dim
    mesh = domain |> gk.mesh
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
    face_around = phi.face_around
    @assert face_around !== nothing
    index -> begin
        sface = index.face
        tface = sface_to_tfaces[sface][face_around]
        lface = sface_to_lfaces[sface][face_around]
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

function prototype(a::DomainMap,::CoboundaryGlue{false,false})
    identity
end

function term(a::DomainMap,::CoboundaryGlue{false,false})
    index -> identity
end

function prototype(a::DomainMap,::CoboundaryGlue)
    error("Case not yet implemented")
end

function term(a::DomainMap,::CoboundaryGlue)
    error("Case not yet implemented")
end

function (phi::DomainMap)(x::AbstractQuantity)
    MappedPoint(phi,x)
end

struct MappedPoint{A,B} <: AbstractQuantity
    phi::A
    x::B
end
function prototype(y::MappedPoint)
    phi_p = gk.prototype(phi)
    x_p = gk.prototype(x)
    phi_p(x_p)
end
function term(y::MappedPoint)
    x = y.x
    phi = y.phi
    term_x = gk.term(x)
    term_phi = gk.term(phi)
    glue = gk.domain_glue(phi)
    term = index -> begin
        f = term_phi(index)
        f(term_x(index))
    end
end
function (a::AbstractQuantity)(y::MappedPoint)
    phi = y.phi
    x = y.x
    (a∘phi)(x)
end

function plot(domain::AbstractDomain;kwargs...)
    plot(domain,gk.domain_style(domain))
end

function plot(domain::AbstractDomain,::GlobalDomain;kwargs...)
    style = domain |> gk.domain_style
    d = gk.face_dim(domain)
    domface_to_face = gk.faces(domain)
    mesh = gk.mesh(domain)
    vismesh = gk.visualization_mesh(mesh,d,domface_to_face;kwargs...)
    node_data = Dict{String,Any}()
    face_data = Dict{String,Any}()
    Plot(domain,vismesh,node_data,face_data)
end

function plot(domain::AbstractDomain,style;kwargs...)
    error("case not implemented")
end

struct Plot{A,B,C,D}
    domain::A
    visualization_mesh::B
    node_data::C
    face_data::D
end
domain(plt::Plot) = plt.domain
visualization_mesh(plt::Plot) = plt.visualization_mesh

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

function coordinates(plt::Plot)
    style = plt |> gk.domain |> gk.domain_style
    gk.coordinates(plt,style)
end

function coordinates(plt::Plot,::GlobalDomain{true})
    gk.reference_coordinates(plt)
end

function coordinates(plt::Plot,::GlobalDomain{false})
    domain_phys = plt |> gk.domain
    domain_ref = domain_phys |> reference_domain
    phi = gk.domain_map(domain_ref,domain_phys)
    q = gk.reference_coordinates(plt)
    phi(q)
end

function plot!(plt::Plot,field;label)
    plot_impl!(field,plt;label)
    plt
end

function plot!(field,plt::Plot;label)
    plot_impl!(field,plt;label)
    plt
end

function plot_impl!(field,plt::Plot;label)
    q = gk.coordinates(plt)
    f_q = field(q)
    term = gk.term(f_q)
    T = typeof(gk.prototype(f_q))
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
    data, :node_data
end

function vtk_plot(f,filename,args...;kwargs...)
    plt = plot(args...;kwargs...)
    vmesh, = plt.visualization_mesh
    d = gk.face_dim(plt.domain)
    vtk_grid(filename,gk.vtk_args(vmesh,d)...) do vtk
        f(VtkPlot(plt,vtk))
    end
end

struct VtkPlot{A,B}
    plt::A
    vtk::B
end

function plot!(field,plt::VtkPlot;label)
    plot!(plt,field;label)
end

function plot!(plt::VtkPlot,field;label)
    data,data_type = plot_impl!(field,plt.plt;label)
    if data_type === :node_data
        plt.vtk[label,WriteVTK.VTKPointData()] = data
    elseif data_type === :face_data
        plt.vtk[label,WriteVTK.VTKCellData()] = data
    else
        error("Unreachable line reached")
    end
    plt
end


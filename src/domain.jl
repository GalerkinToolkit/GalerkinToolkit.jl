
abstract type AbstractDomain <: GalerkinToolkitDataType end
domain(a::AbstractDomain) = a
mesh(a::AbstractDomain) = a.mesh
mesh_id(a::AbstractDomain) = a.mesh_id
physical_names(a::AbstractDomain) = a.physical_names
face_dim(a::AbstractDomain) = gk.val_parameter(a.face_dim)
is_reference_domain(a::AbstractDomain) = a.is_reference_domain |> gk.val_parameter

function domain(mesh;
    mesh_id = objectid(mesh),
    physical_names=gk.physical_names(mesh,gk.num_dims(mesh)),
    face_dim = gk.num_dims(mesh),
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

function Base.:(==)(a::AbstractDomain,b::AbstractDomain)
    flag = true
    flar = flag && gk.mesh_id(a) == gk.mesh_id(b)
    flar = flag && gk.physical_names(a) == gk.physical_names(b)
    flar = flag && gk.face_dim(a) == gk.face_dim(b)
    flar = flag && gk.is_reference_domain(a) == gk.is_reference_domain(b)
    flag
end

abstract type AbstractDomainStyle end
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

abstract type AbstractDomainGlueStyle end
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
    doman = glue |> gk.domain
    codoman = glue |> gk.codomain
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
    doman = glue |> gk.domain
    codoman = glue |> gk.codomain
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

abstract type AbstractQuantity end
term(a::AbstractQuantity) = a.term
prototype(a::AbstractQuantity) = a.prototype

quantity(term,prototype) = Quantity(term,prototype)
struct Quantity{T,B} <: AbstractQuantity
    term::T
    prototype::B
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
    prototype = gk.return_prototype(g,map(gk.prototype,args)...)
    gk.quantity(prototype) do index
        g(map(f->f(index),fs)...)
    end
end

function (f::AbstractQuantity)(x::AbstractQuantity)
    call(call,f,x)
end

abstract type AbstractField <: AbstractQuantity end
domain(a::AbstractField) = a.domain

function analytical_field(f,dom)
    AnalyticalField(f,dom)
end

struct AnalyticalField{A,B} <: AbstractField
    f::A
    domain::B
end
term(a::AnalyticalField) = index -> a.f
prototype(a::AbstractField) = a.f

function domain_map(domain,codomain;face_around=nothing)
    glue = gk.domain_glue(domain,codomain)
    glue_style = gk.domain_glue_style(glue)
    domain_map(glue,glue_style,face_around)
end

function domain_map(glue::DomainGlue,glue_style,face_around)
    @assert face_around === nothing
    DomainMap(glue,face_around)
end

function domain_map(glue::DomainGlue,::CoboundaryGlue,::Integer)
    DomainMap(glue,face_around)
end

function domain_map(glue::DomainGlue,::CoboundaryGlue,::Nothing)
    face_around_plus = 1
    face_around_minus = 2
    phi_plus = DomainMap(glue,face_around_plus)
    phi_minus = DomainMap(glue,face_around_minus)
    (phi_plus,phi_minus)
end

struct DomainMap{A,B} <: AbstractField
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

function prototype(a::DomainMap,::InteriorGlue)
    identity
end

function term(a::DomainMap,::InteriorGlue{true,true})
    index -> identity
end

function term(a::DomainMap,::InteriorGlue{true,false})
    glue = a |> gk.domain_glue
    domain = glue |> gk.domain
    mesh = domain |> gk.mesh
    d = domain |> gk.face_dim
    node_to_coords = gk.node_coordinates(mesh)
    sface_to_face = domain |> gk.faces
    face_to_nodes = gk.face_nodes(mesh,d)
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

function term(a::DomainMap,::InteriorGlue{false,true})
    error("Physical to reference map not implemented yet")
end

function term(a::DomainMap,::InteriorGlue{false,false})
    index -> identity
end

function Base.:∘(a::AbstractQuantity,phi::DomainMap)
    g = gk.prototype(a)
    f = gk.prototype(phi)
    prototype = x-> g(f(x))
    term_a = gk.term(a)
    term_phi = gk.term(phi)
    glue = domain_glue(phi)
    target_index_term = gk.target_index_term(glue)
    face_around = phi.face_around
    gk.quantity(prototype) do index
        index2 = replace_face_around(index,face_around)
        index3 = target_index_term(index2)
        ai = term_a(index3)
        phii = term_phi(index)
        x-> ai(phii(x))
    end
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

function measure(dom::AbstractDomain,degree)
    measure(default_quadrature,dom,degree)
end
function measure(f,dom::AbstractDomain,degree)
    Measure(f,dom,degree)
end
struct Measure{A,B,C}
    quadrature_rule::A
    domain::B
    degree::C
end
domain(a::Measure) = a.domain

#TODO this is just a reference quadrature
function quadrature(measure::Measure)
    domain = gk.domain(measure)
    @assert length(gk.face_dim(domain)) == 1
    @assert all(gk.is_reference_domain(domain))
    @assert gk.is_permuted_domain == false
    mesh = gk.mesh(domain)
    d = gk.face_dim(domain)[end]
    drefid_refdface = gk.reference_faces(mesh,d)
    ndrefids = length(drefid_refdface)
    refid_to_quad = map(1:ndrefids) do drefid
        refdface = drefid_refdface[drefid]
        geo = gk.geometry(refdface)
        measure.quadrature_rule(geo,measure.degree)
    end
    domface_to_face, = gk.faces(domain)
    face_to_refid = gk.face_reference_id(mesh,d)
    prototype = first(refid_to_quad)
    gk.quantity(prototype) do index
        domface = index.face[1]
        face = domface_to_face[domface]
        refid = face_to_refid[face]
        refid_to_quad[refid]
    end
end

#TODO these are just the reference coordinates
function coordinates(measure::Measure)
    quadrature = gk.quadrature(measure)
    quadrature_term = gk.term(quadrature)
    prototype = gk.coordinates(gk.prototype(quadrature))
    gk.quantity(prototype) do index
        quad = quadrature_term(index)
        gk.coordinates(quad)[index.point]
    end
end

#TODO these are just the reference weights
function weights(measure::Measure)
    quadrature = gk.quadrature(measure)
    quadrature_term = gk.term(quadrature)
    prototype = gk.weights(gk.prototype(quadrature))
    gk.quantity(prototype) do index
        quad = quadrature_term(index)
        gk.weights(quad)[index.point]
    end
end

function num_points(measure::Measure)
    quadrature = gk.quadrature(measure)
    quadrature_term = gk.term(quadrature)
    prototype = 1
    gk.quantity(prototype) do index
        quad = quadrature_term(index)
        length(gk.weights(quad))
    end
end

function integrate(f,measure::Measure)
    q = gk.coordinates(measure)
    w = gk.weights(measure) |> gk.term
    np = gk.num_points(measure) |> gk.term
    f_q = f(q)
    prototype = gk.prototype(w)*gk.prototype(f_q)+gk.prototype(w)*gk.prototype(f_q)
    domain_contribution = gk.quantity(prototype) do index
        sum(1:np(index)) do point
            @assert index.point === nothing
            index2 = replace_point(index,point)
            w(index2)*f_q(index2)
        end
    end
    domain = gk.domain(measure)
    contribution = domain => domain_contribution
    contributions = (contribution,)
    gk.integral(contributions)
end
const ∫ = integrate

integral(contribs) = Integral(contribs)
struct Integral{A}
    contributions::A
end
contributions(a::Integral) = a.contributions

function plot(domain::AbstractDomain;kwargs...)
    @assert length(gk.face_dim(domain)) == 1
    d = gk.face_dim(domain)[1]
    domface_to_face = gk.faces(domain)
    mesh = gk.mesh(domain)
    vmesh, vglue = gk.visualization_mesh(mesh,d,domface_to_face;kwargs...)
    Plot(domain,(vmesh,vglue),Dict{String,Any}(),Dict{String,Any}())
end

struct Plot{A,B,C,D}
    domain::A
    visualization_mesh::B
    node_data::C
    face_data::D
end
domain(plt::Plot) = plt.domain
visualization_mesh(plt::Plot) = plt.visualization_mesh

function coordinates(plt::Plot)
    domain = plt.domain
    d = gk.face_dim(domain)[1]
    domface_to_face = gk.faces(domain)
    mesh = gk.mesh(domain)
    vmesh,vglue = gk.visualization_mesh(plt)
    refid_to_snode_to_coords = vglue.reference_coordinates
    d = gk.num_dims(vmesh)
    face_to_refid = gk.face_reference_id(mesh,d)
    prototype = first(first(refid_to_snode_to_coords))
    # TODO use a domain_map here
    if gk.val_parameter(gk.is_reference_domain(domain)[1])
        gk.quantity(prototype) do index
            domface = index.face[1]
            point = index.point
            face = domface_to_face[domface]
            refid = face_to_refid[face]
            refid_to_snode_to_coords[refid][point]
        end
    else
        node_to_x = gk.node_coordinates(mesh)
        face_to_nodes = gk.face_nodes(mesh,d)
        refid_to_refface = gk.reference_faces(mesh,d)
        refid_to_A = map(refid_to_refface,refid_to_snode_to_coords) do refface, x
            A = gk.tabulator(refface)(value,x)
            A
        end
        gk.quantity(prototype) do index
            domface = index.face[1]
            point = index.point
            face = domface_to_face[domface]
            refid = face_to_refid[face]
            A = refid_to_A[refid]
            nodes = face_to_nodes[face]
            x = view(node_to_x,nodes)
            (A*x)[point]
        end
    end
end

function plot!(plt::Plot,field;label)
    plot_impl!(plt,field;label)
    plt
end

function plot_impl!(plt::Plot,field::AbstractField;label)
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

function plot_impl!(plt::Plot,integral::Integral;label)
    error("Not implemented")
end

function vtk_plot(f,filename,args...;kwargs...)
    plt = plot(args...;kwargs...)
    vmesh, = plt.visualization_mesh
    vtk_grid(filename,gk.vtk_args(vmesh)...) do vtk
        f(VtkPlot(plt,vtk))
    end
end

struct VtkPlot{A,B}
    plt::A
    vtk::B
end

function plot!(plt::VtkPlot,field;label)
    data,data_type = plot_impl!(plt.plt,field;label)
    if data_type === :node_data
        plt.vtk[label,WriteVTK.VTKPointData()] = data
    elseif data_type === :face_data
        plt.vtk[label,WriteVTK.VTKCellData()] = data
    else
        error("Unreachable line reached")
    end
    plt
end



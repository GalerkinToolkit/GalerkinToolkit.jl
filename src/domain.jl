
abstract type AbstractQuantity end
term(a::AbstractQuantity) = a.term

quantity(term) = Quantity(term)
struct Quantity{T} <: AbstractQuantity
    term::T
end

function index(;
    face=nothing,
    field=nothing,
    dof=nothing,
    side=nothing,
    point=nothing)
    Index(face,field,dof,side,point)
end

struct Index{A,B,C,D,E}
    face::A
    field::B
    dof::C
    side::D
    point::E
end

function call(g,args::AbstractQuantity...)
    fs = map(gk.term,args)
    gk.quantity() do index
        g(map(f->f(index),fs)...)
    end
end

abstract type AbstractDomain <: GalerkinToolkitDataType end
domain(a::AbstractDomain) = a
mesh(a::AbstractDomain) = a.mesh
physical_names(a::AbstractDomain) = a.physical_names
face_dims(a::AbstractDomain) = a.face_dims
is_reference_domain(a::AbstractDomain) = a.is_reference_domain
is_permuted_domain(a::AbstractDomain) = a.is_permuted_domain

struct Domain{A,B,C,D,E} <: AbstractDomain
    mesh::A
    physical_names::B
    face_dims::C
    is_reference_domain::D
    is_permuted_domain::E
end

function domain(mesh::AbstractFEMesh;
    physical_names=gk.physical_names(mesh,gk.num_dims(mesh)),
    face_dims = (gk.num_dims(mesh),),
    is_reference_domain = (Val(false),),
    is_permuted_domain = Val(true),
    )

    Domain(
           mesh,
           physical_names,
           face_dims,
           is_reference_domain,
           is_permuted_domain,
          )
end

function face_ids(domain::AbstractDomain)
    @assert length(gk.face_dims(domain)) == 1
    D = face_dims[1]
    Dface_to_tag = zeros(Int,gk.num_faces(mesh,D))
    gk.classify_mesh_faces!(Dface_to_tag,mesh,D,tag_to_name)
    physical_Dfaces = findall(i->i!=0,Dface_to_tag)
    (physical_Dfaces,)
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

function quadrature(measure::Measure)
    domain = gk.domain(measure)
    @assert length(gk.face_dims(domain)) == 1
    @assert all(gk.is_reference_domain(domain))
    @assert gk.is_permuted_domain == false
    mesh = gk.mesh(domain)
    d = gk.face_dims(domain)[end]
    drefid_refdface = gk.reference_faces(mesh,d)
    ndrefids = length(drefid_refdface)
    refid_to_quad = map(1:ndrefids) do drefid
        refdface = drefid_refdface[drefid]
        geo = gk.geometry(refdface)
        measure.quadrature_rule(geo,measure.degree)
    end
    domface_to_face, = gk.face_ids(domain)
    face_to_refid = gk.face_reference_id(mesh,d)
    gk.quantity() do index
        domface = index.face[1]
        face = domface_to_face[domface]
        refid = face_to_refid[face]
        refid_to_quad[refid]
    end
end

function coordinates(measure::Measure)
    quadrature = gk.quadrature(measure)
    quadrature_term = gk.term(quadrature)
    gk.quantity() do index
        quad = quadrature_term(index)
        gk.coordinates(quad)[index.point]
    end
end

function weights(measure::Measure)
    quadrature = gk.quadrature(measure)
    quadrature_term = gk.term(quadrature)
    gk.quantity() do index
        quad = quadrature_term(index)
        gk.weights(quad)[index.point]
    end
end

function num_points(measure::Measure)
    quadrature = gk.quadrature(measure)
    quadrature_term = gk.term(quadrature)
    gk.quantity() do index
        quad = quadrature_term(index)
        length(gk.weights(quad))
    end
end

function integrate(f,measure::Measure)
    q = gk.coordinates(measure)
    w = gk.weights(measure) |> gk.term
    np = gk.num_points(measure) |> gk.term
    f_q = f(q)
    domain_contribution = gk.quantity() do index
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
    @assert length(gk.face_dims(domain)) == 1
    d = face_dims(domain)[1]
    domface_to_face, = gk.face_ids(domain)
    mesh = gk.mesh(domain)
    vmesh, vglue = visualization_mesh(mesh,dim,domface_to_face;kwargs...)
    Plot(domain,(vmesh,vgluevmesh,vglue),Dict{String,Any}(),Dict{String,Any}())
end

struct Plot{A,B,C,D}
    domain::A
    vizualization_mesh::B
    node_data::C
    face_data::D
end

function plot!(plt::Plot,field;label)
    plot_impl!(plt,field;label)
    plt
end

function plot_impl!(plt::Plot,field::AbstractField;label)
    q = gk.coordinates(plt)
    f_q = field(q)
    term = gk.term(f_q)
    T = gk.quantity_type(f_q)
    vmesh,vglue = plt.visualization_mesh
    nnodes = gk.num_nodes(vmesh)
    data = zeros(T,nnodes)
    face_to_nodes = vglue.face_fine_nodes
    for face in 1:gk.num_faces(plt.domain)
        nodes = face_to_nodes[face]
        for point in 1:length(nodes)
            index = gk.index(;face,point)
            data[nodes[point]] = term(index)
        end
    end
    plt.node_data[label] = data
    plt
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
    plot!(plt.plt,field;label)
    plt.vtk[label] = plt.plt.node_data[label]
    plt
end



module GalerkinToolkit

using WriteVTK
using FastGaussQuadrature
using StaticArrays
using LinearAlgebra
using ForwardDiff
using Gmsh
using PartitionedArrays

val_parameter(a) = a
val_parameter(::Val{a}) where a = a

num_dims(a) = val_parameter(a.num_dims)
node_coordinates(a) = a.node_coordinates
reference_faces(a) = a.reference_faces
face_nodes(a) = a.face_nodes
face_reference_id(a) = a.face_reference_id
vtk_mesh_cell(a) = a.vtk_mesh_cell
physical_groups(a) = a.physical_groups
geometry(a) = a.geometry
is_n_cube(a) = hasproperty(a,:is_n_cube) ? val_parameter(a.is_n_cube) : false
is_simplex(a) = hasproperty(a,:is_simplex) ? val_parameter(a.is_simplex) : false
is_axis_aligned(a) = a.is_axis_aligned
bounding_box(a) = a.bounding_box
coordinates(a) = a.coordinates
weights(a) = a.weights
interpolation(a) = a.interpolation
shape_functions(a) = a.shape_functions
gradient!(a) = a.gradient!
value!(a) = a.value!
order(a) = a.order
monomial_exponents(a) = a.monomial_exponents
lex_to_user_nodes(a) = a.lex_to_user_nodes

reference_faces(a,d) = reference_faces(a)[val_parameter(d)+1]
face_nodes(a,d) = face_nodes(a)[val_parameter(d)+1]
face_reference_id(a,d) = face_reference_id(a)[val_parameter(d)+1]
num_faces(a) = map(length,face_reference_id(a))
num_faces(a,d) = length(face_reference_id(a,d))
physical_groups(a,d) = physical_groups(a)[val_parameter(d)+1]
num_nodes(a) = length(node_coordinates(a))
num_ambient_dims(a) = length(eltype(node_coordinates(a)))
function face_offsets(a)
    D = num_dims(a)
    offsets = zeros(Int,D+1)
    for d in 1:D
        offsets[d+1] = offsets[d] + num_faces(a,d-1)
    end
    offsets
end

function is_unit_n_cube(geom)
    !(is_axis_aligned(geom) && is_n_cube(geom)) && return false
    my_bounding_box = bounding_box(geom)
    all(i->i==0,first(my_bounding_box)) && all(i->i==1,last(my_bounding_box))
end

function is_unit_simplex(geom)
    !(is_axis_aligned(geom) && is_simplex(geom)) && return false
    my_bounding_box = bounding_box(geom)
    all(i->i==0,first(my_bounding_box)) && all(i->i==1,last(my_bounding_box))
end

evaluate(f,x) = f(x)

function vtk_points(mesh)
    function barrirer(coords)
        nnodes = length(coords)
        points = zeros(3,nnodes)
        for node in 1:nnodes
            coord = coords[node]
            for i in 1:length(coord)
                points[i,node] = coord[i]
            end
        end
        points
    end
    coords = node_coordinates(mesh)
    barrirer(coords)
end

function vtk_cells(mesh,d)
    function barrirer(face_to_refid,face_to_nodes,refid_mesh_cell)
        cells = map(face_to_refid,face_to_nodes) do refid, nodes
            mesh_cell = refid_mesh_cell[refid]
            mesh_cell(nodes)
        end
        cells
    end
    face_to_nodes = face_nodes(mesh,d)
    face_to_refid = face_reference_id(mesh,d)
    refid_refface = reference_faces(mesh,d)
    refid_mesh_cell = map(vtk_mesh_cell,refid_refface)
    barrirer(face_to_refid,face_to_nodes,refid_mesh_cell)
end

"""
    args = vtk_args(mesh,d)

Return the arguments `args` to be passed in final position
to functions like `WriteVTK.vtk_grid`.
"""
function vtk_args(mesh,d)
    points = vtk_points(mesh)
    cells = vtk_cells(mesh,d)
    points, cells
end

function vtk_args(mesh)
    points = vtk_points(mesh)
    D = num_dims(mesh)
    allcells = [vtk_cells(mesh,d) for d in 0:D if num_faces(mesh,d) != 0]
    cells = reduce(vcat,allcells)
    points, cells
end

function vtk_physical_groups!(vtk,geo,d;physical_groups=physical_groups(geo,d))
    ndfaces = num_faces(geo,d)
    for group in physical_groups
        name,faces = group
        face_mask = zeros(Int,ndfaces)
        face_mask[faces] .= 1
        vtk[name,WriteVTK.VTKCellData()] = face_mask
    end
    vtk
end

function vtk_physical_groups!(vtk,geo;physical_groups=physical_groups(geo))
    nfaces = sum(num_faces(geo))
    offsets = face_offsets(geo)
    D = num_dims(geo)
    data = Dict{String,Vector{Int}}()
    for d in 0:D
        for group in physical_groups[d+1]
            name, = group
            if !haskey(data,name)
                face_mask = zeros(Int,nfaces)
                data[name] = face_mask
            end
        end
    end
    for d in 0:D
        for group in physical_groups[d+1]
            offset = offsets[d+1]
            name,faces = group
            face_mask = data[name]
            face_mask[faces.+offset] .= 1
        end
    end
    for (name,face_mask) in data
        vtk[name,WriteVTK.VTKCellData()] = face_mask
    end
    vtk
end

function classify_mesh_nodes!(node_to_tag,mesh,tag_to_name,dmax=num_dims(mesh))
    fill!(node_to_tag,zero(eltype(node_to_tag)))
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

partition_from_mask(a) = partition_from_mask(identity,a)

function partition_from_mask(f,node_to_mask)
    T = Vector{Int32}
    free_nodes = convert(T,findall(f,node_to_mask))
    dirichlet_nodes = convert(T,findall(i->!f(i),node_to_mask))
    nfree = length(free_nodes)
    ndiri = length(dirichlet_nodes)
    permutation = T(undef,nfree+ndiri)
    permutation[free_nodes] = 1:nfree
    permutation[dirichlet_nodes] = (1:ndiri) .+ nfree
    TwoPartPartition(free_nodes,dirichlet_nodes,permutation)
end

struct TwoPartPartition{A} <: AbstractVector{A}
    first::A
    last::A
    permutation::A
end

permutation(a::TwoPartPartition) = a.permutation
Base.size(a::TwoPartPartition) = (2,)
Base.IndexStyle(::Type{<:TwoPartPartition}) = IndexLinear()
function Base.getindex(a::TwoPartPartition,i::Int)
    @boundscheck @assert i in (1,2)
    if i == 1
        a.first
    else
        a.last
    end
end

function tensor_product_quadrature(
    degree_per_dir,
    limits_per_dir;
    weight_type=Float64,
    coordinate_type=SVector{length(limits_per_dir),weight_type},
    allocator=zeros)
    D = length(degree_per_dir)
    n_per_dir = map(d->ceil(Int,(d+1)/2),degree_per_dir)
    quad_per_dir = map(n_per_dir,limits_per_dir) do n,limits
        x,w = FastGaussQuadrature.gausslegendre(n)
        a,b = limits
        x .= (0.5*(b-a)) .*x .+ (0.5*(b+a))
        w .*= 0.5*(b-a)
        (;coordinates=x,weights=w)
    end
    coords_per_dir = map(coordinates,quad_per_dir)
    weights_per_dir = map(weights,quad_per_dir)
    m = prod(map(length,weights_per_dir))
    w = allocator(weight_type,m)
    x = allocator(coordinate_type,m)
    tensor_product!(identity,x,coords_per_dir)
    tensor_product!(prod,w,weights_per_dir)
    (;coordinates=x,weights=w)
end

function tensor_product!(f,result,values_per_dir)
    shape = Tuple(map(length,values_per_dir))
    cis = CartesianIndices(shape)
    lis = LinearIndices(cis)
    for ci in cis
        li = lis[ci]
        result[li] = f(map((q,i)->q[i],values_per_dir,Tuple(ci)))
    end
    result
end

function quadrature(geom,degree;kwargs...)
    if is_n_cube(geom)
        if is_axis_aligned(geom)
            my_bounding_box = bounding_box(geom)
            D = length(first(my_bounding_box))
            limits_per_dir = ntuple(i->(my_bounding_box[1][i],my_bounding_box[2][i]),Val(D))
            degree_per_dir = ntuple(i->degree,Val(D))
            tensor_product_quadrature(degree_per_dir,limits_per_dir;kwargs...)
        else
        error("Not implemented")
        end
    else
        error("Not implemented")
    end
end

function lagrange_monomial_exponents(f,order,D;kwargs...)
    function tensor_monomial_exponents(f,degree_per_dir; monomial_type=SVector{length(degree_per_dir),Int}, allocator=zeros)
        terms_per_dir = Tuple(map(d->d+1,degree_per_dir))
        D = length(terms_per_dir)
        cis = CartesianIndices(terms_per_dir)
        m = count(ci->f(Tuple(ci) .- 1,degree_per_dir),cis)
        li = 0
        result = zeros(monomial_type,m)
        for ci in cis
            t = Tuple(ci) .- 1
            if f(t,degree_per_dir)
                li += 1
                result[li] = t
            end
        end
        result
    end
    d = val_parameter(D)
    order_per_dir = ntuple(i->order,Val(d))
    tensor_monomial_exponents((e,os)->f(e,order),order_per_dir;kwargs...)
end

function lagrange_shape_functions(lagrange_monomial_exponents,node_coordinates)
    monomials = map(e->(x-> prod(x.^e)),lagrange_monomial_exponents)
    monomials_t = permutedims(monomials)
    A = evaluate.(monomials_t,node_coordinates)
    B = A\I
    function value!(r,x;scratch=similar(r))
        C = broadcast!(evaluate,scratch,monomials_t,x)
        mul!(r,C,B)
        r
    end
    function gradient!(r,x;scratch=similar(r,eltype(x)))
        C = broadcast!(ForwardDiff.gradient,scratch,monomials_t,x)
        mul!(r,C,B)
        r
    end
    (;value!,gradient!)
end

function lagrange_interpolation(
    geom,order;
    type=:default,
    lex_to_user_nodes = nothing,
    coordinate_type=SVector{num_dims(geom),Float64},
    kwargs...)
    f = if is_unit_n_cube(geom)
        if type === :default
            (e,o)->true
        elseif type === :Q
            (e,o)->true
        else
            error("Not implemented")
        end
    elseif is_unit_simplex(geom)
        if type === :default
            (e,o)->sum(e)<=o
        elseif type === :P
            (e,o)->sum(e)<=o
        else
            error("Not implemented")
        end
    else
        error("Not implemented")
    end
    d = num_dims(geom)
    monomial_exponents = lagrange_monomial_exponents(f,order,Val(d);kwargs...)
    T = eltype(coordinate_type)
    node_coordinates = map(monomial_exponents) do exponent
        c = map(exponent) do e
            if order != 0
                T(e/order)
            else
                T(e)
            end
        end
        coordinate_type(c)
    end
    nnodes = length(node_coordinates)
    lex_to_user_nodes = lex_to_user_nodes !== nothing ? lex_to_user_nodes : collect(1:nnodes)
    node_coordinates[lex_to_user_nodes] = node_coordinates
    shape_functions = lagrange_shape_functions(monomial_exponents,node_coordinates)
    (;shape_functions,node_coordinates,order,monomial_exponents,lex_to_user_nodes)
end

function lagrange_reference_face(geometry,args...;kwargs...)
    interpolation = lagrange_interpolation(geometry,args...;kwargs...)
    vtk_mesh_cell = vtk_mesh_cell_from_geometry(geometry,interpolation)
    (;geometry,interpolation,vtk_mesh_cell)
end

function vtk_mesh_cell_from_geometry(geom,interpolation)
    d = num_dims(geom)
    nnodes = num_nodes(interpolation)
    lex_to_user = lex_to_user_nodes(interpolation)
    if d == 0 && nnodes == 1
        cell_type = WriteVTK.VTKCellTypes.VTK_VERTEX
        vtk_to_lex = [1]
    elseif d == 1 && (is_simplex(geom) || is_n_cube(geom)) && nnodes == 2
        cell_type = WriteVTK.VTKCellTypes.VTK_LINE
        vtk_to_lex = [1,2]
    elseif d == 1 && (is_simplex(geom) || is_n_cube(geom)) && nnodes == 3
        cell_type = WriteVTK.VTKCellTypes.VTK_QUADRATIC_EDGE
        vtk_to_lex = [1,3,2]
    elseif d == 2 && is_n_cube(geom) && nnodes == 4
        cell_type = WriteVTK.VTKCellTypes.VTK_QUAD
        vtk_to_lex = [1,2,4,3]
    elseif d == 2 && is_simplex(geom) && nnodes == 3
        cell_type = WriteVTK.VTKCellTypes.VTK_TRIANGLE
        vtk_to_lex = [1,2,3]
    elseif d == 2 && is_simplex(geom) && nnodes == 6
        cell_type = WriteVTK.VTKCellTypes.VTK_QUADRATIC_TRIANGLE
        vtk_to_lex = [1,3,6,2,4,5]
    elseif d == 3 && is_n_cube(geom) && nnodes == 8
        cell_type = WriteVTK.VTKCellTypes.VTK_HEXAHEDRON
        vtk_to_lex = [1,2,4,3,5,6,8,7]
    elseif d == 3 && is_simplex(geom) && nnodes == 4
        cell_type = WriteVTK.VTKCellTypes.VTK_TETRA
        vtk_to_lex = [1,2,3,4]
    else
        error("Not implemented")
    end
    nodes -> WriteVTK.MeshCell(cell_type,nodes[lex_to_user][vtk_to_lex])
end

function unit_n_cube(D)
    d = val_parameter(D)
    num_dims = Val(d)
    is_n_cube=true
    is_axis_aligned=true
    bounding_box=SVector{d,Float64}[ntuple(i->0,Val(d)),ntuple(i->1,Val(d))]
    is_simplex = d in (0,1)
    geometry = (;num_dims,is_n_cube,is_simplex,is_axis_aligned,bounding_box)
    geometry
end

function unit_simplex(D)
    d = val_parameter(D)
    num_dims = Val(d)
    is_simplex=true
    is_axis_aligned=true
    is_n_cube = d in (0,1)
    bounding_box=SVector{d,Float64}[ntuple(i->0,Val(d)),ntuple(i->1,Val(d))]
    geometry = (;num_dims,is_n_cube,is_simplex,is_axis_aligned,bounding_box)
    geometry
end


const INVALID = 0

function default_gmsh_options()
    [
     "General.Terminal"=>1,
     "Mesh.SaveAll"=>1,
     "Mesh.MedImportGroupsOfNodes"=>1
    ]
end

function with_gmsh(f;options=default_gmsh_options())
    gmsh.initialize()
    for (k,v) in options
        gmsh.option.setNumber(k,v)
    end
    try
        return f()
    finally
        gmsh.finalize()
    end
end

function mesh_from_gmsh(file;renumber=true,kwargs...)
    @assert ispath(file) "File not found: $(file)"
    with_gmsh(;kwargs...) do
        gmsh.open(file)
        renumber && gmsh.model.mesh.renumberNodes()
        renumber && gmsh.model.mesh.renumberElements()
        mesh_from_gmsh_module()
    end
end

function mesh_from_gmsh_module()
    entities = gmsh.model.getEntities()
    nodeTags, coord, parametricCoord = gmsh.model.mesh.getNodes()

    # find num_dims
    ddim = -1
    for e in entities
        ddim = max(ddim,e[1])
    end
    if ddim == -1
        error("No entities in the msh file.")
    end
    D = ddim

    # find embedded_dimension
    dtouched = [false,false,false]
    for node in nodeTags
        if !(coord[(node-1)*3+1] + 1 ≈ 1)
            dtouched[1] = true
        end
        if !(coord[(node-1)*3+2] + 1 ≈ 1)
            dtouched[2] = true
        end
        if !(coord[(node-1)*3+3] + 1 ≈ 1)
            dtouched[3] = true
        end
    end
    if dtouched[3]
        adim = 3
    elseif dtouched[2]
        adim = 2
    elseif dtouched[1]
        adim = 1
    else
        adim = 0
    end

    # Setup node coords
    nmin = minimum(nodeTags)
    nmax = maximum(nodeTags)
    nnodes = length(nodeTags)
    if !(nmax == nnodes && nmin == 1)
        error("Only consecutive node tags allowed.")
    end
    my_node_to_coords = zeros(SVector{adim,Float64},nnodes)
    m = zero(MVector{adim,Float64})
    for node in nodeTags
        for j in 1:adim
            k = (node-1)*3 + j
            xj = coord[k]
            m[j] = xj
        end
        my_node_to_coords[node] = m
    end

    # Setup face nodes
    offsets = zeros(Int32,D+1)
    my_face_nodes = Vector{JaggedArray{Int32,Int32}}(undef,D+1)
    for d in 0:D
        elemTypes, elemTags, nodeTags = gmsh.model.mesh.getElements(d)
        ndfaces = 0
        for t in 1:length(elemTypes)
            ndfaces += length(elemTags[t])
        end
        if ndfaces != 0
            nmin::Int = minimum( minimum, elemTags )
            nmax::Int = maximum( maximum, elemTags )
            if !( (nmax-nmin+1) == ndfaces)
                error("Only consecutive elem tags allowed.")
            end
            offsets[d+1] = nmin-1
        end
        ptrs = zeros(Int32,ndfaces+1)
        dface = 0
        for t in 1:length(elemTypes)
            elementName, dim, order, numNodes, nodeCoord =
            gmsh.model.mesh.getElementProperties(elemTypes[t])
            for e in 1:length(elemTags[t])
                dface += 1
                ptrs[dface+1] = numNodes
            end
        end
        length_to_ptrs!(ptrs)
        ndata = ptrs[end]-1
        data = zeros(Int32,ndata)
        dface = 1
        for t in 1:length(elemTypes)
            p = ptrs[dface]-Int32(1)
            for (i,node) in enumerate(nodeTags[t])
                data[p+i] = node
            end
            dface += length(elemTags[t])
        end
        my_face_nodes[d+1] = JaggedArray(data,ptrs)
    end

    # Setup face_reference_id
    my_face_reference_id = Vector{Vector{Int32}}(undef,D+1)
    for d in 0:D
        elemTypes, elemTags, nodeTags = gmsh.model.mesh.getElements(d)
        ndfaces = length(my_face_nodes[d+1])
        dface_to_refid = zeros(Int8,ndfaces)
        refid = 0
        dface = 0
        for t in 1:length(elemTypes)
            refid += 1
            for e in 1:length(elemTags[t])
                dface += 1
                dface_to_refid[dface] = refid
            end
        end
        my_face_reference_id[d+1] = dface_to_refid
    end

    # Setup reference faces
    my_reference_faces = ()
    for d in D:-1:0
        elemTypes, elemTags, nodeTags = gmsh.model.mesh.getElements(d)
        refdfaces = ()
        for t in 1:length(elemTypes)
            refface = reference_face_from_gmsh_eltype(elemTypes[t])
            refdfaces = (refdfaces...,refface)
        end
        if refdfaces == ()
            refdfaces = reference_faces(first(first(my_reference_faces)),d)
        end
        my_reference_faces = (refdfaces,my_reference_faces...)
    end

    ## Setup periodic nodes
    #node_to_main_node = fill(Int32(INVALID),nnodes)
    #for (dim,tag) in entities
    #    tagMaster, nodeTags, nodeTagsMaster, = gmsh.model.mesh.getPeriodicNodes(dim,tag)
    #    for i in 1:length(nodeTags)
    #        node = nodeTags[i]
    #        main_node = nodeTagsMaster[i]
    #        node_to_main_node[node] = main_node
    #    end
    #end
    #my_free_and_periodic_nodes = partition_from_mask(i->i==INVALID,node_to_main_node)
    #my_periodic_nodes = last(my_free_and_periodic_nodes)
    #my_periodic_to_master = node_to_main_node[my_periodic_nodes]
    #my_periodic_to_coeff = ones(length(my_periodic_to_master))

    # Setup physical groups
    my_groups = [ Pair{String,Vector{Int32}}[] for d in 0:D]
    for d in 0:D
        offset = Int32(offsets[d+1])
        dimTags = gmsh.model.getPhysicalGroups(d)
        for (dim,tag) in dimTags
            @boundscheck @assert dim == d
            g_entities = gmsh.model.getEntitiesForPhysicalGroup(dim,tag)
            ndfaces_in_physical_group = 0
            for entity in g_entities
                elemTypes, elemTags, nodeTags = gmsh.model.mesh.getElements(dim,entity)
                for t in 1:length(elemTypes)
                    ndfaces_in_physical_group += length(elemTags[t])
                end
            end
            dfaces_in_physical_group = zeros(Int32,ndfaces_in_physical_group)
            ndfaces_in_physical_group = 0
            for entity in g_entities
                elemTypes, elemTags, nodeTags = gmsh.model.mesh.getElements(dim,entity)
                for t in 1:length(elemTypes)
                    for etag in elemTags[t]
                        ndfaces_in_physical_group += 1
                        dfaces_in_physical_group[ndfaces_in_physical_group] = Int32(etag)-offset
                    end
                end
            end
            groupname = gmsh.model.getPhysicalName(dim,tag)
            my_group = groupname => dfaces_in_physical_group
            push!(my_groups[d+1],my_group)
        end
    end

    mesh = (;
            num_dims=Val(D),
            node_coordinates = my_node_to_coords,
            face_nodes = my_face_nodes,
            face_reference_id = my_face_reference_id,
            reference_faces = my_reference_faces,
            physical_groups = my_groups)
    mesh
end

function reference_face_from_gmsh_eltype(eltype)
    if eltype == 1
        order = 1
        geom = unit_n_cube(Val(1))
        lex_to_gmsh = [1,2]
    elseif eltype == 2
        order = 1
        geom = unit_simplex(Val(2))
        lex_to_gmsh = [1,2,3]
    elseif eltype == 3
        order = 1
        geom = unit_n_cube(Val(2))
        lex_to_gmsh = [1,2,4,3]
    elseif eltype == 4
        order = 1
        geom = unit_simplex(Val(3))
        lex_to_gmsh = [1,2,3,4]
    elseif eltype == 5
        order = 1
        lex_to_gmsh = [1,2,4,3,5,6,8,7]
    elseif eltype == 15
        order = 1
        geom = unit_n_cube(Val(0))
        lex_to_gmsh = [1]
    elseif eltype == 8
        order = 2
        geom = unit_n_cube(Val(1))
        lex_to_gmsh = [1,3,2]
    elseif eltype == 9
        order = 2
        geom = unit_simplex(Val(2))
        lex_to_gmsh = [1,4,2,6,5,3]
    else
        en, = gmsh.model.mesh.getElementProperties(eltype)
        error("Unsupported element type. elemType: $eltype ($en)")
    end
    lagrange_reference_face(geom,order;lex_to_user_nodes=lex_to_gmsh)
end


end # module

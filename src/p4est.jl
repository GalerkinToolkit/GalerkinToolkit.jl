
pxest_topology(::P4est.LibP4est.p4est) = Val(4)
pxest_topology(::P4est.LibP4est.p8est) = Val(8)
pxest_topology(::Ptr{P4est.LibP4est.p4est}) = Val(4)
pxest_topology(::Ptr{P4est.LibP4est.p8est}) = Val(8)
pxest_topology(::P4est.p4est_tree_t) = Val(4)
pxest_topology(::P4est.p8est_tree_t) = Val(8)
pxest_topology(::P4est.p4est_connectivity_t) = Val(4)
pxest_topology(::P4est.p8est_connectivity_t) = Val(8)
pxest_topology(::Ptr{P4est.p4est_ghost_t}) = Val(4)
pxest_topology(::Ptr{P4est.p8est_ghost_t}) = Val(8)
pxest_topology(::P4est.p4est_lnodes_t) = Val(4)
pxest_topology(::P4est.p8est_lnodes_t) = Val(8)
pxest_quadrant_t(::Val{4}) = P4est.p4est_quadrant_t
pxest_quadrant_t(::Val{8}) = P4est.p8est_quadrant_t
pxest_tree_t(::Val{4}) = P4est.p4est_tree_t
pxest_tree_t(::Val{8}) = P4est.p8est_tree_t

struct Leaf4
    p4est_quadrant::P4est.p4est_quadrant_t
end

struct Leaf8
    p4est_quadrant::P4est.p8est_quadrant_t
end

const Leaf = Union{Leaf4,Leaf8}
Leaf(p4est_quadrant::P4est.p4est_quadrant_t) = Leaf4(p4est_quadrant)
Leaf(p4est_quadrant::P4est.p8est_quadrant_t) = Leaf8(p4est_quadrant)

function anchor(a::Leaf4)
    SVector(a.p4est_quadrant.x,a.p4est_quadrant.y)
end

function anchor(a::Leaf8)
    SVector(a.p4est_quadrant.x,a.p4est_quadrant.y,a.p4est_quadrants.z)
end

function level(a::Leaf)
    Int(a.p4est_quadrant.level)
end

struct Tree4 <: AbstractVector{Leaf4}
    p4est_tree::P4est.p4est_tree_t
end

struct Tree8 <: AbstractVector{Leaf8}
    p4est_tree::P4est.p8est_tree_t
end

const Tree = Union{Tree4,Tree8}
Tree(p4est_tree::P4est.p4est_tree_t) = Tree4(p4est_tree)
Tree(p4est_tree::P4est.p8est_tree_t) = Tree8(p4est_tree)

Base.IndexStyle(::Type{<:Tree}) = IndexLinear()

function Base.size(tree::Tree)
    quadrants =  tree.p4est_tree.quadrants
    nquadrants = Int(quadrants.elem_count)
    (nquadrants,)
end

function Base.getindex(tree::Tree,iquadrant::Int)
    T = pxest_topology(tree.p4est_tree)
    quadrants =  tree.p4est_tree.quadrants
    p4est_quadrant = unsafe_load(Ptr{pxest_quadrant_t(T)}(quadrants.array),iquadrant)
    Leaf(p4est_quadrant)
end

mutable struct Forest4{A} <: AbstractVector{Tree4}
    ground::A
    p4est_ptr::Ptr{P4est.LibP4est.p4est}
end

mutable struct Forest8{A} <: AbstractVector{Tree8}
    ground::A
    p4est_ptr::Ptr{P4est.LibP4est.p8est}
end

const Forest = Union{Forest4{A},Forest8{A}} where A
Forest(ground,p4est_ptr::Ptr{P4est.LibP4est.p4est}) = Forest4(ground,p4est_ptr)
Forest(ground,p4est_ptr::Ptr{P4est.LibP4est.p8est}) = Forest8(ground,p4est_ptr)

Base.IndexStyle(::Type{<:Forest}) = IndexLinear()

function Base.size(forest::Forest)
    p4est_ptr = forest.p4est_ptr
    p4est = unsafe_load(p4est_ptr)
    trees_ptr = p4est.trees
    trees = unsafe_load(trees_ptr)
    ntrees = Int(trees.elem_count)
    (ntrees,)
end

function Base.getindex(forest::Forest,itree::Int)
    p4est_ptr = forest.p4est_ptr
    p4est = unsafe_load(p4est_ptr)
    trees_ptr = p4est.trees
    trees = unsafe_load(trees_ptr)
    T = pxest_topology(p4est)
    p4est_tree = unsafe_load(Ptr{pxest_tree_t(T)}(trees.array),itree)
    Tree(p4est_tree)
end

dimension(a::Forest) = dimension(a.ground)
embedded_dimension(a::Forest) = embedded_dimension(a.ground)
max_level(a::Forest) = P4est.P4EST_MAXLEVEL
ground(a::Forest) = a.ground

function forest_from_mesh(coarse_mesh,initial_level=0)
    if !MPI.Initialized()
        MPI.Init()
    end
    D = dimension(coarse_mesh)
    node_to_coords = node_coordinates(coarse_mesh)
    cell_to_nodes = face_nodes(coarse_mesh,D)
    n_nodes = length(node_to_coords)
    n_cells = length(cell_to_nodes)
    refcells = reference_faces(coarse_mesh,D)
    @assert length(refcells) == 1
    refcell = first(refcells)
    T = pxest_find_topology(refcell)
    conn_ptr = pxest_init_connectivity(T,n_nodes,n_cells)
    perm = pxest_node_permutation(refcell)
    pxest_fill_connectivity!(T,conn_ptr,D,node_to_coords,cell_to_nodes,perm)
    pxest_connectivity_complete(T,conn_ptr)
    @boundscheck @assert Bool(pxest_connectivity_is_valid(T,conn_ptr))
    min_quadrants = 0
    fill_uniform = 1
    p4est_ptr = pxest_new_ext(
      T,
      MPI.COMM_WORLD,
      conn_ptr,
      min_quadrants,
      initial_level,
      fill_uniform,
      sizeof(Int),
      C_NULL,
      C_NULL)
    forest = Forest(coarse_mesh,p4est_ptr)
    finalizer(forest) do forest
        p4est_ptr = forest.p4est_ptr
        p4est = unsafe_load(p4est_ptr)
        T = pxest_topology(p4est)
        conn_ptr = p4est.connectivity
        pxest_connectivity_destroy(T,conn_ptr)
        pxest_destroy(T,forest.p4est_ptr)
        nothing
    end
end

"""
    perm = pxest_node_permutation(elem)
    nodeids_pxest_order = nodeids_elem_order[perm]
"""
function pxest_node_permutation end

function pxest_find_topology(refcell)
    D = dimension(refcell)
    if is_hypercube(refcell) && D==2
        Val(4)
    elseif is_hypercube(refcell) && D==3
        Val(8)
    else
        error("cell topology not implemented yet")
    end
end

function pxest_init_connectivity(::Val{4},n_nodes,n_cells)
    p4est_connectivity_new(
        p4est_topidx_t(n_nodes),
        p4est_topidx_t(n_cells),
        p4est_topidx_t(0),
        p4est_topidx_t(0))
end

function pxest_init_connectivity(::Val{8},n_nodes,n_cells)
    p8est_connectivity_new(
        p4est_topidx_t(n_nodes),
        p4est_topidx_t(n_cells),
        p4est_topidx_t(0),
        p4est_topidx_t(0),
        p4est_topidx_t(0),
        p4est_topidx_t(0))
end

function pxest_fill_connectivity!(T,conn_ptr,Dc,node_to_coords,cell_to_nodes,perm)
    n_nodes = length(node_to_coords)
    n_cells = length(cell_to_nodes)
    D = length(eltype(node_to_coords))
    conn = unsafe_load(conn_ptr)
    vertices = unsafe_wrap(Array, conn.vertices, n_nodes*3)
    current = 1
    for node in 1:n_nodes
        coords = node_to_coords[node]
        for j in 1:D
            vertices[current] = Cdouble(coords[j])
            current += 1
        end
        if D==2
            vertices[current] = Cdouble(0.0)
            current += 1
        end
    end
    tree_to_vertex = unsafe_wrap(Array, conn.tree_to_vertex, n_cells*(2^Dc))
    current = 1
    for cell in 1:n_cells
        nodes = cell_to_nodes[cell]
        for inode in perm
            node = nodes[inode]
            tree_to_vertex[current]=p4est_topidx_t(node-1)
            current += 1
        end
    end
    PXEST_FACES = pxest_num_faces(T)
    tree_to_tree = unsafe_wrap(Array, conn.tree_to_tree, conn.num_trees*PXEST_FACES )
    tree_to_face = unsafe_wrap(Array, conn.tree_to_face, conn.num_trees*PXEST_FACES )
    n_trees = conn.num_trees
    for tree in 1:n_trees
      for face in 1:PXEST_FACES
        tree_to_tree[PXEST_FACES * (tree-1) + face] = tree-1
        tree_to_face[PXEST_FACES * (tree-1) + face] = face-1
      end
    end
    conn_ptr
end

pxest_num_faces(::Val{4}) = 4
pxest_num_faces(::Val{8}) = 6
pxest_connectivity_complete(::Val{4},conn_ptr) = p4est_connectivity_complete(conn_ptr)
pxest_connectivity_complete(::Val{8},conn_ptr) = p8est_connectivity_complete(conn_ptr)
pxest_connectivity_is_valid(::Val{4},conn_ptr) = p4est_connectivity_is_valid(conn_ptr)
pxest_connectivity_is_valid(::Val{8},conn_ptr) = p8est_connectivity_is_valid(conn_ptr)
pxest_new_ext(::Val{4},args...) = p4est_new_ext(args...)
pxest_new_ext(::Val{8},args...) = p8est_new_ext(args...)
pxest_destroy(::Val{4},p4est_ptr) = P4est.p4est_destroy(p4est_ptr)
pxest_destroy(::Val{8},p4est_ptr) = P4est.p8est_destroy(p4est_ptr)
pxest_connectivity_destroy(::Val{4},conn_ptr) = P4est.p4est_connectivity_destroy(conn_ptr)
pxest_connectivity_destroy(::Val{8},conn_ptr) = P4est.p8est_connectivity_destroy(conn_ptr)

function node_coordinates!(node_to_coordinates,forest::Forest4,tree_id,leaf::Leaf4)
    p4est_ptr = forest.p4est_ptr
    p4est = unsafe_load(p4est_ptr)
    trees_ptr = p4est.trees
    trees = unsafe_load(trees_ptr)
    ntrees = Int(trees.elem_count)
    conn_ptr = p4est.connectivity
    conn = unsafe_load(conn_ptr)
    nvertices = Int(conn.num_vertices)
    tree_to_vertex_ptr = conn.tree_to_vertex
    tree_to_vertices = unsafe_wrap(Array,tree_to_vertex_ptr,(4,ntrees))
    vertices = view(tree_to_vertices,:,tree_id)
    vertex_coords_ptr = conn.vertices
    vertex_coords = unsafe_wrap(Array,vertex_coords_ptr,(3,nvertices,))
    node = 0
    q0 = anchor(leaf)
    l = level(leaf)
    L = max_level(forest)
    for jy in 0:1
        for jx in 0:1
            q = q0 + 2^(L-l)*SVector(jx,jy)
            s2 = q ./ 2^L
            s1 = 1. .- s2
            s = (s1,s2)
            ivertex = 0
            x = zero(SVector{2,Float64})
            for iy in 1:2
                yfactor = s[iy][2]
                for ix in 1:2
                    ivertex += 1
                    xfactor = yfactor*s[ix][1]
                    xc1 = vertex_coords[1,1+vertices[ivertex]]
                    xc2 = vertex_coords[2,1+vertices[ivertex]]
                    xc = SVector(xc1,xc2)
                    x += xfactor*xc
                end
            end
            node += 1
            node_to_coordinates[node] = x
        end
    end
end

function refine!(f,forest::Forest,recursive=0)
    function callback(p4est_ptr,which_tree,quadrant_ptr)
        quadrant = unsafe_load(quadrant_ptr)
        Cint(f(Int(which_tree)+1,Leaf(quadrant)))
    end
    p4est_ptr = forest.p4est_ptr
    T = pxest_topology(p4est_ptr)
    pxest_refine(callback,T,p4est_ptr,recursive)
    forest
end

function pxest_refine(callback,::Val{4},p4est_ptr,recursive)
    callback_ptr = @cfunction $callback Cint (Ptr{P4est.p4est_t},P4est.p4est_topidx_t,Ptr{P4est.p4est_quadrant_t})
    P4est.p4est_refine(p4est_ptr,recursive,callback_ptr,C_NULL)
end

function pxest_refine(callback,::Val{8},p4est_ptr,recursive)
    callback_ptr = @cfunction $callback Cint (Ptr{P4est.p8est_t},P4est.p4est_topidx_t,Ptr{P4est.p8est_quadrant_t})
    P4est.p8est_refine(p4est_ptr,recursive,callback_ptr,C_NULL)
end

function coarsen!(f,forest::Forest,recursive=0)
    function callback(p4est_ptr,which_tree,leafs_ptr_ptr)
        leafs_ptr = unsafe_load(leafs_ptr_ptr)
        leafs = unsafe_wrap(Array,leafs_ptr,4)
        Cint(f(Int(which_tree)+1,leafs))
    end
    p4est_ptr = forest.p4est_ptr
    T = pxest_topology(p4est_ptr)
    pxest_coarsen(callback,T,p4est_ptr,recursive)
    forest
end

function pxest_coarsen(callback,::Val{4},p4est_ptr::Ptr{P4est.p4est_t},recursive)
    callback_ptr = @cfunction $callback Cint (Ptr{P4est.p4est_t},P4est.p4est_topidx_t,Ptr{Ptr{Leaf4}})
    P4est.p4est_coarsen(p4est_ptr,recursive,callback_ptr,C_NULL)
end

function pxest_coarsen(callback,::Val{8},p4est_ptr,recursive)
    callback_ptr = @cfunction $callback Cint (Ptr{P4est.p8est_t},P4est.p4est_topidx_t,Ptr{Ptr{Leaf8}})
    P4est.p8est_coarsen(p4est_ptr,recursive,callback_ptr,C_NULL)
end

const CONNECT_CORNER = P4est.P4EST_CONNECT_CORNER
const CONNECT_FACE = P4est.P4EST_CONNECT_FACE
const CONNECT_FULL = P4est.P4EST_CONNECT_FULL

function balance!(forest::Forest,connection=CONNECT_CORNER)
    p4est_ptr = forest.p4est_ptr
    T = pxest_topology(p4est_ptr)
    pxest_balance(T,p4est_ptr,connection)
    forest
end

function pxest_balance(::Val{4},p4est_ptr::Ptr{P4est.p4est_t},connection)
    P4est.p4est_balance(p4est_ptr,connection,C_NULL)
end

function pxest_balance(::Val{8},p4est_ptr,connection)
    P4est.p8est_balance(p4est_ptr,connection,C_NULL)
end

function partition!(forest::Forest,allow_for_coarsening=false)
    p4est_ptr = forest.p4est_ptr
    T = pxest_topology(p4est_ptr)
    pxest_partition(T,p4est_ptr,allow_for_coarsening)
    forest
end

function partition!(f,forest::Forest,allow_for_coarsening)
    function callback(p4est_ptr,which_tree,quadrant_ptr)
        quadrant = unsafe_load(quadrant_ptr)
        Cint(f(Int(which_tree)+1,Leaf(quadrant)))
    end
    p4est_ptr = forest.p4est_ptr
    T = pxest_topology(p4est_ptr)
    pxest_partition(callback,T,p4est_ptr,allow_for_coarsening)
    forest
end

function pxest_partition(::Val{4},p4est_ptr::Ptr{P4est.p4est_t},allow_for_coarsening)
    p4est_partition(p4est_ptr, allow_for_coarsening, C_NULL)
end

function pxest_partition(::Val{8},p4est_ptr,allow_for_coarsening)
    p8est_partition(p4est_ptr, allow_for_coarsening, C_NULL)
end

function pxest_partition(callback,::Val{4},p4est_ptr::Ptr{P4est.p4est_t},allow_for_coarsening)
    callback_ptr = @cfunction $callback Cint (Ptr{P4est.p4est_t},P4est.p4est_topidx_t,Ptr{P4est.p4est_quadrant_t})
    p4est_partition(p4est_ptr,allow_for_coarsening,callback_ptr)
end

function pxest_partition(callback,::Val{8},p4est_ptr,allow_for_coarsening)
    callback_ptr = @cfunction $callback Cint (Ptr{P4est.p8est_t},P4est.p4est_topidx_t,Ptr{P4est.p8est_quadrant_t})
    p8est_partition(p4est_ptr,allow_for_coarsening,callback_ptr)
end

mutable struct GhostLeafs4 <: AbstractVector{Leaf4}
    p4est_ghost_ptr::Ptr{P4est.p4est_ghost_t}
end

mutable struct GhostLeafs8 <: AbstractVector{Leaf8}
    p4est_ghost_ptr::Ptr{P4est.p8est_ghost_t}
end

const GhostLeafs = Union{GhostLeafs4,GhostLeafs8}
GhostLeafs(p4est_ghost_ptr::Ptr{P4est.p4est_ghost_t}) = GhostLeafs4(p4est_ghost_ptr)
GhostLeafs(p4est_ghost_ptr::Ptr{P4est.p8est_ghost_t}) = GhostLeafs8(p4est_ghost_ptr)

Base.IndexStyle(::Type{<:GhostLeafs}) = IndexLinear()

function Base.size(ghost_leafs::GhostLeafs)
    p4est_ghost_ptr = ghost_leafs.p4est_ghost_ptr
    p4est_ghost = unsafe_load(p4est_ghost_ptr)
    ghosts = p4est_ghost.ghosts
    nghosts = Int(ghosts.elem_count)
    (nghosts,)
end

function Base.getindex(ghost_leafs::GhostLeafs,iquadrant::Int)
    p4est_ghost_ptr = ghost_leafs.p4est_ghost_ptr
    p4est_ghost = unsafe_load(p4est_ghost_ptr)
    ghosts = p4est_ghost.ghosts
    T = pxest_topology(p4est_ghost_ptr)
    p4est_quadrant = unsafe_load(Ptr{pxest_quadrant_t(T)}(ghosts.array),iquadrant)
    Leaf(p4est_quadrant)
end

function find_ghost_leafs(forest::Forest,connection=CONNECT_FULL)
    p4est_ptr = forest.p4est_ptr
    T = pxest_topology(p4est_ptr)
    p4est_ghost_ptr = pxest_ghost_new(T,p4est_ptr,connection)
    ghost = GhostLeafs(p4est_ghost_ptr)
    finalizer(ghost) do ghost
        pxest_ghost_destroy(T,ghost.p4est_ghost_ptr)
    end
end

function pxest_ghost_new(::Val{4},p4est_ptr,connection)
    P4est.p4est_ghost_new(p4est_ptr,connection)
end

function pxest_ghost_new(::Val{8},p4est_ptr,connection)
    P4est.p8est_ghost_new(p4est_ptr,connection)
end

function pxest_ghost_destroy(::Val{4},p4est_ghost_ptr)
    P4est.p4est_ghost_destroy(p4est_ghost_ptr)
end

function pxest_ghost_destroy(::Val{8},p4est_ghost_ptr)
    P4est.p4est_ghost_destroy(p4est_ghost_ptr)
end

function mesh_from_forest(forest::Forest,order=1,ghost_leafs::GhostLeafs=find_ghost_leafs(forest))
    p4est_ptr = forest.p4est_ptr
    T = pxest_topology(p4est_ptr)
    ghost_ptr = ghost_leafs.p4est_ghost_ptr
    lnodes_ptr = pxest_lnodes_new(T,p4est_ptr,ghost_ptr,order)
    lnodes = unsafe_load(lnodes_ptr)
    leaf_to_nodes = pxest_leaf_nodes(lnodes_ptr,order)
    leaf_to_code, code_to_constraints = pxest_leaf_constraints(lnodes_ptr,order)
    nfree = lnodes.num_local_nodes
    constraints = pxest_numerate_hanging_nodes!(nfree,leaf_to_nodes,leaf_to_code,code_to_constraints)
    node_coordinates = pxest_node_coordinates(T,size(constraints,1),leaf_to_nodes,lnodes,order,forest)
    leaf_reference_id = Fill(1,length(leaf_to_nodes))
    reference_leafs = (pxest_reference_leaf(T),)
    pxest_lnodes_destroy(T,lnodes_ptr)
    D = pxest_dimension(T)
    face_to_nodes = [JaggedArray{Int32,Int32}([Int32[]]) for d in 0:D]
    face_to_nodes[end] = leaf_to_nodes
    face_reference_id = [Fill(1,0) for d in 0:D]
    face_reference_id[end] = leaf_reference_id
    reference_faces = (ntuple(i->(),D)...,reference_leafs)
    mesh = GenericMesh(node_coordinates,face_to_nodes,face_reference_id,reference_faces)
    set_haning_node_constraints(mesh,constraints)
end

pxest_dimension(::Val{4}) = 2
pxest_dimension(::Val{8}) = 3

pxest_reference_leaf(::Val{4}) = P4estQuad{Float64}()

struct P4estQuad{T}
    node_coordinates::Vector{SVector{2,T}}
end
function P4estQuad{T}() where T
    P4estQuad(SVector{2,T}[(0,0),(1,0),(0,1),(1,1)])
end
dimension(a::P4estQuad) = 2
is_hypercube(a::P4estQuad) = true
node_coordinates(a::P4estQuad) = a.node_coordinates
function face_reference_id(a::P4estQuad,d)
    d==0 && return fill(Int8(1),4)
    d==1 && return fill(Int8(1),4)
    d==2 && return fill(Int8(1),1)
    throw(DomainError(d))
end
function reference_faces(::Type{<:P4estQuad},::Val{0})
    [Meshes.Point(SVector{0,Float64}())]
end
function reference_faces(::Type{<:P4estQuad},::Val{1})
    [Meshes.Segment(Meshes.Point(0),Meshes.Point(1))]
end
function reference_faces(::Type{<:P4estQuad},::Val{2})
    [P4estQuad()]
end
reference_faces(a::P4estQuad,::Val{d}) where d = ref_faces(typeof(a),Val(d))
function face_nodes(a::P4estQuad,d)
    d==0 && return GenericJaggedArray(Vector{Int32}[[1],[2],[3],[4]])
    d==1 && return GenericJaggedArray(Vector{Int32}[[1,3],[2,4],[1,2],[3,4]])
    d==2 && return GenericJaggedArray(Vector{Int32}[[1,2,3,4]])
    throw(DomainError(d))
end
const P4EST_VERTEX_PERMUTATION = [1,2,4,3]
function vtk_mesh_cell(a::P4estQuad)
    nodes -> WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_QUAD,nodes[P4EST_VERTEX_PERMUTATION])
end

function pxest_lnodes_new(::Val{4},p4est_ptr,ghost_ptr,order)
    P4est.p4est_lnodes_new(p4est_ptr,ghost_ptr,order)
end

function pxest_lnodes_new(::Val{8},p4est_ptr,ghost_ptr,order)
    P4est.p8est_lnodes_new(p4est_ptr,ghost_ptr,order)
end

function pxest_lnodes_destroy(::Val{4},lnodes_ptr)
    P4est.p4est_lnodes_destroy(lnodes_ptr)
end

function pxest_lnodes_destroy(::Val{8},lnodes_ptr)
    P4est.p8est_lnodes_destroy(lnodes_ptr)
end

function pxest_leaf_nodes(lnodes_ptr,order)
    lnodes = unsafe_load(lnodes_ptr)
    T = pxest_topology(lnodes)
    n_cell_nodes = pxest_num_leaf_nodes(T,order)
    n_cells = Int(lnodes.num_local_elements)
    cell_to_nodes_ptrs = fill(Int32(n_cell_nodes),n_cells+1)
    cell_to_nodes_ptrs[1] = 0
    length_to_ptrs!(cell_to_nodes_ptrs)
    cell_to_nodes_data = copy(unsafe_wrap(Array,lnodes.element_nodes,n_cells*n_cell_nodes))
    cell_to_nodes_data .+= 1
    cell_to_nodes = JaggedArray(cell_to_nodes_data,cell_to_nodes_ptrs)
    cell_to_nodes
end

function pxest_num_leaf_nodes(::Val{4},order)
    @assert order == 1
    4
end

function pxest_num_leaf_nodes(::Val{8},order)
    @assert order == 1
    8
end

function pxest_node_coordinates(::Val{t},n,leaf_to_nodes,lnodes,order,forest) where t
  T = Val(t)
  D = pxest_dimension(T)
  node_to_coords = zeros(SVector{D,Float64},n)
  ileaf = 0
  x1 = similar(node_to_coords,pxest_num_leaf_nodes(T,order))
  for (itree,tree) in enumerate(forest)
      for leaf in tree
          ileaf += 1
          nodes = leaf_to_nodes[ileaf]
          node_coordinates!(x1,forest,itree,leaf)
          for i in 1:length(x1)
              node_to_coords[nodes[i]] = x1[i]
          end
      end
  end
  node_to_coords
end

function pxest_leaf_constraints(lnodes_ptr,order)
    @assert order == 1
    lnodes = unsafe_load(lnodes_ptr)
    T = pxest_topology(lnodes)
    n_cells = Int(lnodes.num_local_elements)
    cache = pxest_setup_hanging_nodes_cache(T)
    face_code_ptr = lnodes.face_code
    face_code = collect(Int32,unsafe_wrap(Array,face_code_ptr,(n_cells,)))
    free_nodes = Int32[]
    hanging_nodes = Int[]
    master_nodes = JaggedArray{Int32,Int32}([Int32[]])
    master_coeffs = JaggedArray{Float64,Int32}([Int32[]])
    permutation = Int32[]
    A = typeof(HangingNodeConstraints(0,free_nodes,hanging_nodes,master_nodes,master_coeffs,permutation))
    code_constraints = Dict{Int32,A}()
    for cell in 1:n_cells
        code = face_code[cell]
        if !haskey(code_constraints,code)
            code_constraints[code] = pxest_setup_hanging_nodes(T,code,cache)
        end
    end
    (face_code, code_constraints)
end

pxest_corner_faces(::Val{4}) = collect(transpose(unsafe_wrap(Array,cglobal((:p4est_corner_faces,P4est.LibP4est.libp4est),Cint),(2,4))))
pxest_face_corners(::Val{4}) = collect(transpose(unsafe_wrap(Array,cglobal((:p4est_face_corners,P4est.LibP4est.libp4est),Cint),(2,4))))
pxest_corner_faces(::Val{8}) = collect(transpose(unsafe_wrap(Array,cglobal((:p8est_corner_faces,P4est.LibP4est.libp4est),Cint),(3,8))))
pxest_face_corners(::Val{8}) = collect(transpose(unsafe_wrap(Array,cglobal((:p8est_face_corners,P4est.LibP4est.libp4est),Cint),(4,6))))
pxest_corner_edges(::Val{8}) = collect(transpose(unsafe_wrap(Array,cglobal((:p8est_corner_edges,P4est.LibP4est.libp4est),Cint),(3,8))))
pxest_edge_corners(::Val{8}) = collect(transpose(unsafe_wrap(Array,cglobal((:p8est_edge_corners,P4est.LibP4est.libp4est),Cint),(2,12))))

function pxest_setup_hanging_nodes_cache(::Val{4})
    (;
     corner_faces = pxest_corner_faces(Val(4)),
     face_corners = pxest_face_corners(Val(4))
    )
end

function pxest_setup_hanging_nodes_cache(::Val{8})
    (;
     corner_faces = pxest_corner_faces(Val(8)),
     face_corners = pxest_face_corners(Val(8)),
     corner_edges = pxest_corner_edges(Val(8)),
     edge_corners = pxest_edge_corners(Val(8)),
    )
end

function pxest_setup_hanging_nodes(::Val{4},face_code,cache)
    p4est_corner_faces = cache.corner_faces
    p4est_face_corners = cache.face_corners
    p4est_half = div(P4est.P4EST_CHILDREN,2*one(P4est.P4EST_CHILDREN))
    ones = Cint(P4est.P4EST_CHILDREN - 1)
    c = Cint(face_code & ones)
    c1 = c+one(c)
    work = Cint(face_code >> P4est.P4EST_DIM)
    hanging_nodes = Int32[]
    master_nodes = Vector{Int32}[]
    master_coeffs = Vector{Float64}[]
    for i1 in one(P4est.P4EST_DIM):P4est.P4EST_DIM
        if work & Cint(1) != Cint(0)
            i = Cint(i1-one(i1))
            ef = p4est_corner_faces[c1,i1];
            ef1 = ef + one(ef)
            hanging = xor(c,xor(ones,(Cint(1) << i)))
            hanging1 = hanging + one(hanging)
            my_masters = Int32[0,0]
            my_coeffs = Float64[0,0]
            for j1 in one(p4est_half):p4est_half
                master = p4est_face_corners[ef1,j1]
                master1 = master + one(master)
                my_masters[j1] = master1
                my_coeffs[j1] = 0.5
            end
            push!(hanging_nodes,hanging1)
            push!(master_nodes,my_masters)
            push!(master_coeffs,my_coeffs)
        end
        work = work >> Cint(1)
    end
    permutation = fill(Int32(INVALID_ID),4)
    nhanging = length(hanging_nodes)
    permutation[hanging_nodes] = (1:nhanging) .+ (4-nhanging)
    free_nodes = collect(Int32,findall(i->i==Int32(INVALID_ID),permutation))
    permutation[free_nodes] = 1:length(free_nodes)
    master_nodes_jagged = JaggedArray(master_nodes)
    aux = JaggedArray(master_coeffs)
    master_coeffs_jagged = JaggedArray(aux.data,master_nodes_jagged.ptrs)
    HangingNodeConstraints(4,free_nodes,hanging_nodes,master_nodes_jagged,master_coeffs_jagged,permutation)
end

function pxest_setup_hanging_nodes(::Val{8},face_code,cache)
    p8est_corner_faces = cache.corner_faces
    p8est_face_corners = cache.face_corners
    p8est_half = div(P4est.P8EST_CHILDREN,2*one(P4est.P8EST_CHILDREN))
    ones = Cint(P4est.P8EST_CHILDREN - 1)
    c = Cint(face_code & ones)
    c1 = c+one(c)
    work = Cint(face_code >> P4est.P8EST_DIM)
    hanging_nodes = Int32[]
    master_nodes = Vector{Int32}[]
    master_coeffs = Vector{Float64}[]
    for i1 in one(P4est.P8EST_DIM):P4est.P8EST_DIM
        if work & Cint(1) != Cint(0)
            i = Cint(i1-one(i1))
            ef = p8est_corner_faces[c1,i1];
            ef1 = ef + one(ef)
            hanging = xor(c,xor(ones,(Cint(1) << i)))
            hanging1 = hanging + one(hanging)
            my_masters = Int32[0,0,0,0]
            my_coeffs = Float64[0,0]
            for j1 in one(p8est_half):p8est_half
                master = p8est_face_corners[ef1,j1]
                master1 = master + one(master)
                my_masters[j1] = master1
                my_coeffs[j1] = 0.25
            end
            push!(hanging_nodes,hanging1)
            push!(master_nodes,my_masters)
            push!(master_coeffs,my_coeffs)
        end
        work = work >> Cint(1)
    end
    p8est_corner_edges = cache.corner_edges
    p8est_edge_corners = cache.edge_corners
    for i1 in one(P4est.P8EST_DIM):P4est.P8EST_DIM
        if work & Cint(1) != Cint(0)
            i = Cint(i1-one(i1))
            ef = p8est_corner_edges[c1,i1];
            ef1 = ef + one(ef)
            hanging = xor(c,(Cint(1) << i))
            hanging1 = hanging + one(hanging)
            my_masters = Int32[0,0]
            my_coeffs = Float64[0,0]
            for k in 1:2
                master = p8est_edge_corners[ef1,k]
                master1 = master + one(master)
                my_masters[k] = master1
                my_coeffs[j1] = 0.5
            end
            push!(hanging_nodes,hanging1)
            push!(master_nodes,my_masters)
            push!(master_coeffs,my_coeffs)
        end
        work = work >> Cint(1)
    end
    permutation = fill(Int32(INVALID_ID),8)
    nhanging = length(hanging_nodes)
    permutation[hanging_nodes] = (1:nhanging) .+ (8-nhanging)
    free_nodes = collect(Int32,findall(i->i!=Int32(INVALID_ID),permutation))
    permutation[free_nodes] = 1:length(free_nodes)
    master_nodes_jagged = JaggedArray(master_nodes)
    aux = JaggedArray(master_coeffs)
    master_coeffs_jagged = JaggedArray(aux.data,master_nodes_jagged.ptrs)
    HangingNodeConstraints(8,free_nodes,hanging_nodes,master_nodes_jagged,master_coeffs_jagged,permutation)
end

function pxest_numerate_hanging_nodes!(n_nodes,leaf_to_nodes,leaf_to_code,code_to_constraints)
    n_pre_hanging = 0
    nleafs = length(leaf_to_nodes)
    for leaf in 1:nleafs
        code = leaf_to_code[leaf]
        constraints = code_to_constraints[code]
        n_pre_hanging += length(hanging_nodes(constraints))
    end
    pre_hanging_to_masters_ptrs = zeros(Int32,n_pre_hanging+1)
    pre_hanging = 0
    for leaf in 1:nleafs
        code = leaf_to_code[leaf]
        constraints = code_to_constraints[code]
        for masters in master_nodes(constraints)
            pre_hanging += 1
            pre_hanging_to_masters_ptrs[pre_hanging+1] = length(masters)
        end
    end
    length_to_ptrs!(pre_hanging_to_masters_ptrs)
    pre_hanging_to_masters_data = zeros(Int32,pre_hanging_to_masters_ptrs[end]-1)
    pre_hanging_to_masters = JaggedArray(pre_hanging_to_masters_data,pre_hanging_to_masters_ptrs)
    pre_hanging = 0
    for leaf in 1:nleafs
        code = leaf_to_code[leaf]
        constraints = code_to_constraints[code]
        for masters in master_nodes(constraints)
            pre_hanging += 1
            mymasters = pre_hanging_to_masters[pre_hanging]
            for (imaster,master) in enumerate(masters)
                mymasters[imaster] = master
            end
        end
    end
    master_to_pre_hanging_ptrs = zeros(Int32,n_nodes+1)
    for masters in pre_hanging_to_masters
        for master in masters
            master_to_pre_hanging_ptrs[1+master] += 1
        end
    end
    max_hanging_per_master = maximum(view(master_to_pre_hanging_ptrs,2:n_nodes))
    length_to_ptrs!(master_to_pre_hanging_ptrs)
    master_to_pre_hanging_data = zeros(Int32,master_to_pre_hanging_ptrs[end]-1)
    for (pre_hanging,masters) in enumerate(pre_hanging_to_masters)
        for master in masters
            p = master_to_pre_hanging_ptrs[master]
            master_to_pre_hanging_data[p] = pre_hanging
            master_to_pre_hanging_ptrs[master] += Int32(1)
        end
    end
    rewind_ptrs!(master_to_pre_hanging_ptrs)
    master_to_pre_hanging = JaggedArray(master_to_pre_hanging_data,master_to_pre_hanging_ptrs)
    pre_hanging_to_hanging = fill(Int32(INVALID_ID),n_pre_hanging)
    pre_hanging_set_1 = fill(Int32(INVALID_ID),max_hanging_per_master)
    pre_hanging_set_2 = fill(Int32(INVALID_ID),max_hanging_per_master)
    hanging = 0
    for pre_hanging in 1:n_pre_hanging
        if pre_hanging_to_hanging[pre_hanging] != Int32(INVALID_ID)
            continue
        end
        hanging += 1
        masters = pre_hanging_to_masters[pre_hanging]
        fill!(pre_hanging_set_1,Int32(INVALID_ID))
        fill!(pre_hanging_set_2,Int32(INVALID_ID))
        master = first(masters)
        pre_hangings_1 = master_to_pre_hanging[master]
        n_pre_hanging_set_1 = length(pre_hangings_1)
        copyto!(pre_hanging_set_1,pre_hangings_1)
        for master in masters
            pre_hangings_2 = master_to_pre_hanging[master]
            n_pre_hanging_set_2 = length(pre_hangings_2)
            copyto!(pre_hanging_set_2,pre_hangings_2)
            intersection!(pre_hanging_set_1,pre_hanging_set_2,n_pre_hanging_set_1,n_pre_hanging_set_2)
        end
        for pre_hanging_1 in pre_hanging_set_1
            if pre_hanging_1 != Int32(INVALID_ID)
                pre_hanging_to_hanging[pre_hanging_1] = hanging
            end
        end
    end
    n_hanging = hanging
    pre_hanging_to_hanging .+= n_nodes
    pre_hanging = 0
    for leaf in 1:nleafs
        code = leaf_to_code[leaf]
        constraints = code_to_constraints[code]
        mynodes = leaf_to_nodes[leaf]
        for myhaning in hanging_nodes(constraints)
            pre_hanging += 1
            mynodes[myhaning] = pre_hanging_to_hanging[pre_hanging]
        end
    end
    hanging_to_masters_ptrs = zeros(Int32,n_hanging+1)
    for pre_hanging in 1:n_pre_hanging
        hanging = pre_hanging_to_hanging[pre_hanging] - n_nodes
        masters = pre_hanging_to_masters[pre_hanging]
        hanging_to_masters_ptrs[hanging+1] = length(masters)
    end
    length_to_ptrs!(hanging_to_masters_ptrs)
    hanging_to_masters_data = zeros(Int32,hanging_to_masters_ptrs[end]-1)
    hanging_to_coeffs_data = zeros(Float64,hanging_to_masters_ptrs[end]-1)
    hanging_to_masters = JaggedArray(hanging_to_masters_data,hanging_to_masters_ptrs)
    hanging_to_coeffs = JaggedArray(hanging_to_coeffs_data,hanging_to_masters_ptrs)
    for pre_hanging in 1:n_pre_hanging
        masters = pre_hanging_to_masters[pre_hanging]
        hanging = pre_hanging_to_hanging[pre_hanging] - n_nodes
        hanging_to_masters[hanging] .= masters
        hanging_to_coeffs[hanging] .= 1 ./ length(masters)
    end
    my_hanging_nodes = (1:n_hanging) .+ n_nodes
    my_free_nodes = 1:n_nodes
    mypermutation = 1:(n_hanging+n_nodes)
    HangingNodeConstraints(n_nodes,my_free_nodes,my_hanging_nodes,hanging_to_masters,hanging_to_coeffs,mypermutation)
end

const INVALID_ID = 0

function intersection!(a,b,na,nb)
  function findeq!(i,a,b,nb)
    for j in 1:nb
      if a[i] == b[j]
        return
      end
    end
    a[i] = INVALID_ID
    return
  end
  for i in 1:na
    if a[i] == INVALID_ID
      continue
    end
    findeq!(i,a,b,nb)
  end
end

#pxest_dimension(p4est::P4est.LibP4est.p4est) = 2
#pxest_dimension(p4est::P4est.LibP4est.p8est) = 3
#pxest_dimension(p4est::P4est.p4est_tree_t) = 2
#pxest_dimension(p4est::P4est.p8est_tree_t) = 3
#pxest_dimension(p4est::P4est.p4est_quadrant_t) = 2
#pxest_dimension(p4est::P4est.p8est_quadrant_t) = 3
#pxest_tree_t(::Val{2}) = P4est.p4est_tree_t
#pxest_tree_t(::Val{3}) = P4est.p8est_tree_t
#pxest_quadrant_t(::Val{2}) = P4est.p4est_quadrant_t
#pxest_quadrant_t(::Val{3}) = P4est.p8est_quadrant_t
#
#
#
#
#
#
#function pxest_setup_lnodes(p4est_ptr::Ptr{P4est.LibP4est.p4est})
#    order = 1
#    ghost_ptr = P4est.p4est_ghost_new(p4est_ptr,P4est.P4EST_CONNECT_FULL)
#    lnodes_ptr = P4est.p4est_lnodes_new(p4est_ptr,ghost_ptr,order)
#    lnodes_ptr
#end
#
#function pxest_setup_lnodes(p4est_ptr::Ptr{P4est.LibP4est.p8est})
#    order = 1
#    ghost_ptr = P4est.p8est_ghost_new(p8est_ptr,P4est.P8EST_CONNECT_FULL)
#    lnodes_ptr = P4est.p8est_lnodes_new(p8est_ptr,ghost_ptr,order)
#    lnodes_ptr
#end
#
#function pxest_node_ids(p4est_ptr,lnodes_ptr=pxest_setup_lnodes(p4est_ptr))
#    p4est = unsafe_load(p4est_ptr)
#    lnodes = unsafe_load(lnodes_ptr)
#    D = pxest_dimension(p4est)
#    n_cell_nodes = 2^D
#    n_cells = Int(lnodes.num_local_elements)
#    cell_to_nodes_ptrs = fill(Int32(n_cell_nodes),n_cells+1)
#    cell_to_nodes_ptrs[1] = 0
#    length_to_ptrs!(cell_to_nodes_ptrs)
#    cell_to_nodes_data = copy(unsafe_wrap(Array,lnodes.element_nodes,n_cells*n_cell_nodes))
#    cell_to_nodes = JaggedArray(cell_to_nodes_data,cell_to_nodes_ptrs)
#    trees_ptr = p4est.trees
#    trees = unsafe_load(trees_ptr)
#    n_trees = Int(trees.elem_count)
#    face_code_ptr = lnodes.face_code
#    face_code = unsafe_wrap(Array,face_code_ptr,(n_cells,))
#    face_code_to_hanging_elem_nodes = Dict{Int32,Vector{Int32}}()
#    face_code_to_master_elem_nodes = Dict{Int32,JaggedArray{Int32,Int32}}()
#    n_pre_hanging = 0
#    cache = pxest_setup_hanging_nodes_cache(Val(D))
#    for cell in 1:n_cells
#        code = face_code[cell]
#        if !haskey(face_code_to_hanging_elem_nodes,code)
#            hanging_elem_nodes, master_elem_nodes = pxest_setup_hanging_nodes(Val(D),code,cache)
#            face_code_to_hanging_elem_nodes[code] = hanging_elem_nodes
#            face_code_to_master_elem_nodes[code] = master_elem_nodes
#        end
#        hanging_elem_nodes = face_code_to_hanging_elem_nodes[code]
#        n_pre_hanging += length(hanging_elem_nodes)
#    end
#    pre_hanging_to_masters_ptrs = zeros(Int32,n_pre_hanging+1)
#    pre_hanging = 0
#    for cell in 1:n_cells
#        code = face_code[cell]
#        master_elem_nodes = face_code_to_master_elem_nodes[code]
#        for my_elem_masters in master_elem_nodes
#            pre_hanging += 1
#            pre_hanging_to_masters_ptrs[pre_hanging+1] = length(my_elem_masters)
#        end
#    end
#    length_to_ptrs!(pre_hanging_to_masters_ptrs)
#    pre_hanging_to_masters_data = zeros(Int32,pre_hanging_to_masters_ptrs[end]-1)
#    pre_hanging_to_masters = JaggedArray(pre_hanging_to_masters_data,pre_hanging_to_masters_ptrs)
#    pre_hanging = 0
#    for cell in 1:n_cells
#        my_nodes = cell_to_nodes[cell]
#        code = face_code[cell]
#        master_elem_nodes = face_code_to_master_elem_nodes[code]
#        for my_elem_masters in master_elem_nodes
#            pre_hanging += 1
#            masters = pre_hanging_to_masters[pre_hanging]
#            for imaster in 1:length(masters)
#                masters[imaster] = my_nodes[my_elem_masters[imaster]]
#            end
#        end
#    end
#    n_nodes = lnodes.num_local_nodes
#    master_to_pre_hanging_ptrs = zeros(Int32,n_nodes+1)
#    for masters in pre_hanging_to_masters
#        for master in masters
#            master_to_pre_hanging_ptrs[1+master] += 1
#        end
#    end
#    max_hanging_per_master = maximum(view(master_to_pre_hanging_ptrs,2:n_nodes))
#    length_to_ptrs!(master_to_pre_hanging_ptrs)
#    master_to_pre_hanging_data = zeros(Int32,master_to_pre_hanging_ptrs[end]-1)
#    for (pre_hanging,masters) in enumerate(pre_hanging_to_masters)
#        for master in masters
#            p = master_to_pre_hanging_ptrs[master]
#            master_to_pre_hanging_data[p] = pre_hanging
#            master_to_pre_hanging_ptrs[master] += Int32(1)
#        end
#    end
#    rewind!(master_to_pre_hanging_ptrs)
#    master_to_pre_hanging = JaggedArray(master_to_pre_hanging_data,master_to_pre_hanging_ptrs)
#    n_pre_hanging = length(pre_hanging_to_masters)
#    pre_hanging_to_hanging = fill(Int32(INVALID_ID),n_pre_hanging)
#    pre_hanging_set_1 = fill(Int32(INVALID_ID),max_hanging_per_master)
#    pre_hanging_set_2 = fill(Int32(INVALID_ID),max_hanging_per_master)
#    haning = 0
#    for pre_hanging in 1:n_pre_hanging
#        if pre_hanging_to_hanging[pre_hanging] != Int32(INVALID_ID)
#            continue
#        end
#        haning += 1
#        masters = pre_hanging_to_masters[pre_hanging]
#        fill!(pre_hanging_set_1,Int32(INVALID_ID))
#        fill!(pre_hanging_set_2,Int32(INVALID_ID))
#        master = first(masters)
#        pre_hangings_1 = master_to_pre_hanging[master]
#        n_pre_hanging_set_1 = length(pre_hangings_1)
#        copyto!(pre_hanging_set_1,pre_hangings_1)
#        for master in masters
#            pre_hangings_2 = master_to_pre_hanging[master]
#            n_pre_hanging_set_2 = length(pre_hangings_2)
#            copyto!(pre_hanging_set_2,pre_hangings_2)
#            intersection!(pre_hanging_set_1,pre_hanging_set_2,n_pre_hanging_set_1,n_pre_hanging_set_2)
#        end
#        for pre_hanging_1 in pre_hanging_set_1
#            if pre_hanging_1 != Int32(INVALID_ID)
#                pre_hanging_to_hanging[pre_hanging_1] = haning
#            end
#        end
#    end
#    n_hanging = haning
#    pre_hanging_to_hanging .+= n_nodes
#    pre_hanging = 0
#    for cell in 1:n_cells
#        code = face_code[cell]
#        hanging_elem_nodes = face_code_to_hanging_elem_nodes[code]
#        my_nodes = cell_to_nodes[cell]
#        for my_haning in hanging_elem_nodes
#            pre_hanging += 1
#            my_nodes[my_haning] = pre_hanging_to_hanging[pre_hanging]
#        end
#    end
#    hanging_to_masters_ptrs = zeros(Int32,n_hanging+1)
#    for pre_hanging in 1:n_pre_hanging
#        haning = pre_hanging_to_hanging[pre_hanging]
#        masters = pre_hanging_to_masters[pre_hanging]
#        hanging_to_masters_ptrs[hanging+1] = length(masters)
#    end
#    length_to_ptrs!(hanging_to_masters_ptrs)
#    hanging_to_masters_data = zeros(Int32,hanging_to_masters_ptrs[end]-1)
#    hanging_to_coeffs_data = zeros(Int32,hanging_to_masters_ptrs[end]-1)
#    hanging_to_masters = JaggedArray(hanging_to_masters_data,hanging_to_masters_ptrs)
#    hanging_to_coeffs = JaggedArray(hanging_to_coeffs_data,hanging_to_masters_ptrs)
#    for pre_hanging in 1:n_pre_hanging
#        masters = pre_hanging_to_masters[pre_hanging]
#        haning = pre_hanging_to_hanging[pre_hanging]
#        hanging_to_masters[haning] .= masters
#        hanging_to_coeffs[haning] .= 1 ./ length(masters)
#    end
#    hanging_nodes = (1:n_hanging) .+ n_nodes
#    cell_to_nodes, hanging_nodes, hanging_to_masters, hanging_to_coeffs
#end
#
#    #for itree in 1:n_trees
#    #    tree = unsafe_load(Ptr{pxest_tree_t(Val(D))}(trees.array),itree)
#    #    quadrants =  tree.quadrants
#    #    n_quadrants = Int(quadrants.elem_count)
#    #    for iquadrant in 1:n_quadrants
#    #        cell += 1
#    #        n_pre_hanging += pxest_count_hanging(Val(D),face_code[cell])
#    #    end
#    #end
#
#
#function intersection!(a,b,na,nb)
#  function findeq!(i,a,b,nb)
#    for j in 1:nb
#      if a[i] == b[j]
#        return
#      end
#    end
#    a[i] = INVALID_ID
#    return
#  end
#  for i in 1:na
#    if a[i] == INVALID_ID
#      continue
#    end
#    findeq!(i,a,b,nb)
#  end
#end


# Point
dimension(a::Meshes.Point) = 0
dimension(a::Type{<:Meshes.Point}) = 0
embedded_dimension(a::Meshes.Point) = Meshes.embeddim(a)
is_simplex(a::Meshes.Point) = true
is_hypercube(a::Meshes.Point) = true
function face_reference_id(a::Meshes.Point,d)
  d!=0 && throw(DomainError(d))
  [Int8(1)]
end
function reference_faces(::Type{<:Meshes.Point},::Val{d}) where d
  d!=0 && throw(DomainError(d))
  [Meshes.Point(SVector{0,Float64}())]
end
reference_faces(a::Meshes.Point,::Val{d}) where d = reference_faces(typeof(a),Val(d))
function face_nodes(a::Meshes.Point,::Val{0})
  JaggedArray([[Int32(1)]])
end
node_coordinates(a::Meshes.Point) = [Meshes.coordinates(a)]
function vtk_mesh_cell(a::Meshes.Point)
  nodes -> WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_VERTEX,nodes)
end
vtk_cell_type(a::Meshes.Point) = WriteVTK.VTKCellTypes.VTK_VERTEX

# Segment
dimension(a::Meshes.Segment) = 1
dimension(a::Type{<:Meshes.Segment}) = 1
embedded_dimension(a::Meshes.Segment) = Meshes.embeddim(a)
is_simplex(a::Meshes.Segment) = true
is_hypercube(a::Meshes.Segment) = true
function face_reference_id(a::Meshes.Segment,d)
  d==0 && return [Int8(1),Int8(1)]
  d==1 && return [Int8(1)]
  throw(DomainError(d))
end
function reference_faces(::Type{<:Meshes.Segment},::Val{0})
  [Meshes.Point(SVector{0,Float64}())]
end
function reference_faces(::Type{<:Meshes.Segment},::Val{1})
  [Meshes.Segment(Meshes.Point(0),Meshes.Point(1))]
end
reference_faces(a::Meshes.Segment,::Val{d}) where d = reference_faces(typeof(a),Val(d))
function face_nodes(a::Meshes.Segment,d)
  d==0 && return JaggedArray([[Int32(1)],[Int32(2)]])
  d==1 && return JaggedArray([[Int32(1)],[Int32(2)]])
  throw(DomainError(d))
end
node_coordinates(a::Meshes.Segment) = collect(Meshes.coordinates.(Meshes.vertices(a)))
function vtk_mesh_cell(a::Meshes.Segment)
  nodes -> WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_LINE,nodes)
end
vtk_cell_type(a::Meshes.Segment) = WriteVTK.VTKCellTypes.VTK_LINE

# Quadrangle
dimension(a::Meshes.Quadrangle) = 2
dimension(a::Type{<:Meshes.Quadrangle}) = 2
embedded_dimension(a::Meshes.Quadrangle) = Meshes.embeddim(a)
is_hypercube(a::Meshes.Quadrangle) = true
function face_reference_id(a::Meshes.Quadrangle,d)
  d==0 && return fill(Int8(1),4)
  d==1 && return fill(Int8(1),4)
  d==2 && return fill(Int8(1),1)
  throw(DomainError(d))
end
function reference_faces(::Type{<:Meshes.Quadrangle},::Val{0})
  [Meshes.Point(SVector{0,Float64}())]
end
function reference_faces(::Type{<:Meshes.Quadrangle},::Val{1})
  [Meshes.Segment(Meshes.Point(0),Meshes.Point(1))]
end
function reference_faces(::Type{<:Meshes.Quadrangle},::Val{2})
  [Meshes.Quadrangle(Meshes.Point.([(0,0),(1,0),(1,1),(0,1)]))]
end
reference_faces(a::Meshes.Quadrangle,::Val{d}) where d = reference_faces(typeof(a),Val(d))
function face_nodes(a::Meshes.Quadrangle,d)
  d==0 && return JaggedArray(Vector{Int32}[[1],[2],[3],[4]])
  d==1 && return JaggedArray(Vector{Int32}[[1,2],[2,3],[3,4],[4,1]])
  d==2 && return JaggedArray(Vector{Int32}[[1,2,3,4]])
  throw(DomainError(d))
end
node_coordinates(a::Meshes.Quadrangle) = collect(Meshes.coordinates.(Meshes.vertices(a)))
function vtk_mesh_cell(a::Meshes.Quadrangle)
  nodes -> WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_QUAD,nodes)
end
vtk_cell_type(a::Meshes.Quadrangle) = WriteVTK.VTKCellTypes.VTK_QUAD
pxest_node_permutation(a::Meshes.Quadrangle) = [1,2,4,3]
function add_physical_groups_hypercube(a::Meshes.Quadrangle)
    add_trivial_physical_groups(a)
end

# Triangle
dimension(a::Meshes.Triangle) = 2
dimension(a::Type{<:Meshes.Triangle}) = 2
embedded_dimension(a::Meshes.Triangle) = Meshes.embeddim(a)
is_simplex(a::Meshes.Triangle) = true
function face_reference_id(a::Meshes.Triangle,d)
  d==0 && return fill(Int8(1),3)
  d==1 && return fill(Int8(1),3)
  d==2 && return fill(Int8(1),1)
  throw(DomainError(d))
end
function reference_faces(::Type{<:Meshes.Triangle},::Val{0})
  [Meshes.Point(SVector{0,Float64}())]
end
function reference_faces(::Type{<:Meshes.Triangle},::Val{1})
  [Meshes.Segment(Meshes.Point(0),Meshes.Point(1))]
end
function reference_faces(::Type{<:Meshes.Triangle},::Val{2})
  [Meshes.Triangle(Meshes.Point.([(0,0),(1,0),(0,1)]))]
end
reference_faces(a::Meshes.Triangle,::Val{d}) where d = reference_faces(typeof(a),Val(d))
function face_nodes(a::Meshes.Triangle,d)
  d==0 && return JaggedArray(Vector{Int32}[[1],[2],[3]])
  d==1 && return JaggedArray(Vector{Int32}[[1,2],[2,3],[3,1]])
  d==2 && return JaggedArray(Vector{Int32}[[1,2,3]])
  throw(DomainError(d))
end
node_coordinates(a::Meshes.Triangle) = collect(Meshes.coordinates.(Meshes.vertices(a)))
function vtk_mesh_cell(a::Meshes.Triangle)
  nodes -> WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_TRIANGLE,nodes)
end
vtk_cell_type(a::Meshes.Triangle) = WriteVTK.VTKCellTypes.VTK_TRIANGLE

# Tetrahedron
dimension(a::Meshes.Tetrahedron) = 3
dimension(a::Type{<:Meshes.Tetrahedron}) = 3
embedded_dimension(a::Meshes.Tetrahedron) = Meshes.embeddim(a)
is_simplex(a::Meshes.Tetrahedron) = true
node_coordinates(a::Meshes.Tetrahedron) = collect(Meshes.coordinates.(Meshes.vertices(a)))
function face_reference_id(a::Meshes.Tetrahedron,d)
  d==0 && return fill(Int8(1),4)
  d==1 && return fill(Int8(1),6)
  d==2 && return fill(Int8(1),4)
  d==3 && return fill(Int8(1),1)
  throw(DomainError(d))
end
function reference_faces(::Type{<:Meshes.Tetrahedron},::Val{0})
  [Meshes.Point(SVector{0,Float64}())]
end
function reference_faces(::Type{<:Meshes.Tetrahedron},::Val{1})
  [Meshes.Segment(Meshes.Point(0),Meshes.Point(1))]
end
function reference_faces(::Type{<:Meshes.Tetrahedron},::Val{2})
  [Meshes.Triangle(Meshes.Point.([(0,0),(1,0),(0,1)]))]
end
function reference_faces(::Type{<:Meshes.Tetrahedron},::Val{3})
  [Meshes.Tetrahedron(Meshes.Point.([(0,0,0),(1,0,0),(0,1,0),(0,0,1)]))]
end
reference_faces(a::Meshes.Tetrahedron,::Val{d}) where d = reference_faces(typeof(a),Val(d))
function face_nodes(a::Meshes.Tetrahedron,d)
  d==0 && return JaggedArray(Vector{Int32}[[1],[2],[3],[4]])
  d==1 && return JaggedArray(Vector{Int32}[[1,2],[2,3],[3,1],[1,4],[2,4],[3,4]])
  d==2 && return JaggedArray(Vector{Int32}[[1,3,2],[1,2,4],[2,3,4],[3,1,4]])
  d==3 && return JaggedArray(Vector{Int32}[[1,2,3,4]])
  throw(DomainError(d))
end
function vtk_mesh_cell(a::Meshes.Tetrahedron)
  nodes -> WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_TETRA,nodes)
end
vtk_cell_type(a::Meshes.Tetrahedron) = WriteVTK.VTKCellTypes.VTK_TERA

# Hexahedron
dimension(a::Meshes.Hexahedron) = 3
dimension(a::Type{<:Meshes.Hexahedron}) = 3
embedded_dimension(a::Meshes.Hexahedron) = Meshes.embeddim(a)
is_hypercube(a::Meshes.Hexahedron) = true
node_coordinates(a::Meshes.Hexahedron) = collect(Meshes.coordinates.(Meshes.vertices(a)))
function face_reference_id(a::Meshes.Hexahedron,d)
  d==0 && return fill(Int8(1),8)
  d==1 && return fill(Int8(1),12)
  d==2 && return fill(Int8(1),6)
  d==3 && return fill(Int8(1),1)
  throw(DomainError(d))
end
function reference_faces(::Type{<:Meshes.Hexahedron},::Val{0})
  [Meshes.Point(SVector{0,Float64}())]
end
function reference_faces(::Type{<:Meshes.Hexahedron},::Val{1})
  [Meshes.Segment(Meshes.Point(0),Meshes.Point(1))]
end
function reference_faces(::Type{<:Meshes.Hexahedron},::Val{2})
  [Meshes.Quadrangle(Meshes.Point.([(0,0),(1,0),(1,1),(0,1)]))]
end
function reference_faces(::Type{<:Meshes.Hexahedron},::Val{3})
  [Meshes.Hexahedron(Meshes.Point.([(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)]))]
end
reference_faces(a::Meshes.Hexahedron,::Val{d}) where d = reference_faces(typeof(a),Val(d))
function face_nodes(a::Meshes.Hexahedron,d)
  d==0 && return JaggedArray(Vector{Int32}[[1],[2],[3],[4],[5],[6],[7],[8]])
  d==1 && return JaggedArray(Vector{Int32}[[1,2],[2,3],[3,4],[4,1],[5,6],[6,7],[7,8],[8,5],[1,5],[2,6],[3,7],[4,8]])
  d==2 && return JaggedArray(Vector{Int32}[[1,4,3,2],[1,2,6,5],[2,3,7,6],[3,4,8,7],[4,1,5,8],[5,6,7,8]])
  d==3 && return JaggedArray(Vector{Int32}[[1,2,3,4,5,6,7,8]])
  throw(DomainError(d))
end
function vtk_mesh_cell(a::Meshes.Hexahedron)
  nodes -> WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_HEXAHEDRON,nodes)
end
vtk_cell_type(a::Meshes.Hexahedron) = WriteVTK.VTKCellTypes.VTK_HEXAHEDRON


function complexify(geom::Union{UnitNCube{0},UnitSimplex{0}})
    Tv = real_type(options(geom))
    Ti = int_type(options(geom))
    Tr = reference_int_type(options(geom))
    node_coordinates = [SVector{0,Tv}()]
    face_nodes = [Ti[1]]
    face_reference_id = [Tr[1]]
    space = lagrange_space(geom)
    reference_spaces = ((space,),)
    mesh(;
                node_coordinates,
                face_nodes,
                face_reference_id,
                reference_spaces)
end

function complexify(geom::Union{UnitNCube{1},UnitSimplex{1}})
    Tv = real_type(options(geom))
    Ti = int_type(options(geom))
    Tr = reference_int_type(options(geom))
    geom0 = if is_unit_n_cube(geom)
        unit_n_cube(Val(0);options=options(geom))
    else
        unit_simplex(Val(0);options=options(geom))
    end
    geom1 = geom
    space0 = lagrange_space(geom0)
    space1 = lagrange_space(geom1)
    node_coordinates = SVector{1,Tv}[(0,),(1,)]
    face_nodes = Vector{Vector{Ti}}[[[1],[2]],[[1,2]]]
    face_reference_id = Vector{Tr}[[1,1],[1]]
    reference_spaces = ((space0,),(space1,))
    outwards_normals = SVector{1,Tv}[(-1,),(1,)]
    mesh(;
                node_coordinates,
                face_nodes,
                face_reference_id,
                reference_spaces,
                outward_normals
               )
end

function complexify(geom::UnitSimplex{2})
    Tv = real_type(options(geom))
    Ti = int_type(options(geom))
    Tr = reference_int_type(options(geom))
    geom0 = unit_simplex(Val(0);options=options(geom))
    geom1 = unit_simplex(Val(1);options=options(geom))
    geom2 = geom
    space0 = lagrange_space(geom0)
    space1 = lagrange_space(geom1)
    space2 = lagrange_space(geom2)
    node_coordinates = SVector{2,Tv}[(0,0),(1,0),(0,1)]
    face_nodes = Vector{Vector{Ti}}[[[1],[2],[3]],[[1,2],[1,3],[2,3]],[[1,2,3]]]
    face_reference_id = Vector{Tr}[[1,1,1],[1,1,1],[1]]
    reference_spaces = ((space0,),(space1,),(space2,))
    n1 = sqrt(2)/2
    outwards_normals = SVector{2,Tv}[(0,-1),(-1,0),(n1,n1)]
    mesh(;
                node_coordinates,
                face_nodes,
                face_reference_id,
                reference_spaces,
                outward_normals
               )
end

function complexify(geom::UnitNCube{2})
    Tv = real_type(options(geom))
    Ti = int_type(options(geom))
    Tr = reference_int_type(options(geom))
    geom0 = unit_n_cube(Val(0);options=options(geom))
    geom1 = unit_n_cube(Val(1);options=options(geom))
    geom2 = geom
    space0 = lagrange_space(geom0)
    space1 = lagrange_space(geom1)
    space2 = lagrange_space(geom2)
    node_coordinates = SVector{2,Tv}[(0,0),(1,0),(0,1),(1,1)]
    face_nodes = Vector{Vector{Ti}}[[[1],[2],[3],[4]],[[1,2],[3,4],[1,3],[2,4]],[[1,2,3,4]]]
    face_reference_id = Vector{Tr}[[1,1,1,1],[1,1,1,1],[1]]
    outwards_normals = SVector{2,Tv}[(0,-1),(0,1),(-1,0),(1,0)]
    reference_spaces = ((space0,),(space1,),(space2,))
    mesh(;
                node_coordinates,
                face_nodes,
                face_reference_id,
                reference_spaces,
                outward_normals
               )
end

function complexify(geom::UnitSimplex{3})
    Tv = real_type(options(geom))
    Ti = int_type(options(geom))
    Tr = reference_int_type(options(geom))
    geom0 = unit_simplex(Val(0);options=options(geom))
    geom1 = unit_simplex(Val(1);options=options(geom))
    geom2 = unit_simplex(Val(2);options=options(geom))
    geom3 = geom
    space0 = lagrange_space(geom0)
    space1 = lagrange_space(geom1)
    space2 = lagrange_space(geom2)
    space3 = lagrange_space(geom3)
    node_coordinates = SVector{3,Tv}[(0,0,0),(1,0,0),(0,1,0),(0,0,1)]
    face_nodes = [
                  Vector{Ti}[[1],[2],[3],[4]],
                  Vector{Ti}[[1,2],[1,3],[2,3],[1,4],[2,4],[3,4]],
                  Vector{Ti}[[1,2,3],[1,2,4],[1,3,4],[2,3,4]],
                  Vector{Ti}[[1,2,3,4]]
                 ]
    face_reference_id = [ones(Tr,4),ones(Tr,6),ones(Tr,4),ones(Tr,1)]
    n1 = sqrt(3)/3
    outwards_normals = SVector{3,Tv}[(0,0,-1),(0,-1,0),(-1,0,0),(n1,n1,n1)]
    reference_spaces = ((space0,),(space1,),(space2,),(space3,))
    mesh(;
                node_coordinates,
                face_nodes,
                face_reference_id,
                reference_spaces,
                outward_normals
               )
end

function complexify(geom::UnitNCube{3})
    Tv = real_type(options(geom))
    Ti = int_type(options(geom))
    Tr = reference_int_type(options(geom))
    geom0 = unit_n_cube(Val(0);options=options(geom))
    geom1 = unit_n_cube(Val(1);options=options(geom))
    geom2 = unit_n_cube(Val(2);options=options(geom))
    geom3 = geom
    space0 = lagrange_space(geom0)
    space1 = lagrange_space(geom1)
    space2 = lagrange_space(geom2)
    space3 = lagrange_space(geom3)
    node_coordinates = SVector{3,Tv}[(0,0,0),(1,0,0),(0,1,0),(1,1,0),(0,0,1),(1,0,1),(0,1,1),(1,1,1)]
    face_nodes = [
                  Vector{Ti}[[1],[2],[3],[4],[5],[6],[7],[8]],
                  Vector{Ti}[[1,2],[3,4],[1,3],[2,4],[5,6],[7,8],[5,7],[6,8],[1,5],[3,7],[2,6],[4,8]],
                  Vector{Ti}[[1,2,3,4],[5,6,7,8],[1,2,5,6],[3,4,7,8],[1,3,5,7],[2,4,6,8]],
                  Vector{Ti}[[1,2,3,4,5,6,7,8]],
                 ]
    face_reference_id = [ones(Tr,8),ones(Tr,12),ones(Tr,6),ones(Tr,1)]
    outwards_normals = SVector{3,Tv}[(0,0,-1),(0,0,1),(0,-1,0),(0,1,0),(-1,0,0),(1,0,0)]
    reference_spaces = ((space0,),(space1,),(space2,),(space3,))
    mesh(;
                node_coordinates,
                face_nodes,
                face_reference_id,
                reference_spaces,
                outward_normals
               )
end

function opposite_faces(geom::Union{UnitNCube{0},UnitSimplex{0}})
    Ti = int_type(options(geom))
    Vector{Ti}[[1]]
end

function opposite_faces(geom::Union{UnitNCube{1},UnitSimplex{1}})
    Ti = int_type(options(geom))
    Vector{Ti}[[2,1],[1]]
end

function opposite_faces(geom::UnitNCube{2})
    Ti = int_type(options(geom))
    Vector{Ti}[[4,3,2,1],[2,1,4,3],[1]]
end

function opposite_faces(geom::UnitNCube{3})
    Ti = int_type(options(geom))
    Vector{Ti}[[8,7,6,5,4,3,2,1],[6,5,8,7,2,1,4,3,12,11,10,9],[2,1,4,3,6,5],[1]]
end

function complexify(refface::LagrangeFaceSpace)
    geom = domain(refface)
    mesh = complexify(geom)
    D = num_dims(geom)
    order_inter = order(refface)
    round_coord(coord) = map(xi->round(Int,order_inter*xi),coord)
    node_coordinates = GT.node_coordinates(refface)
    node_coordinates_mapped = map(round_coord,node_coordinates)
    Ti = int_type(options(geom))
    dims = ntuple(d->d-1,Val(D+1))
    face_reference_id = GT.face_reference_id(mesh)
    reference_spaces = map(GT.reference_domains(mesh)) do rid_to_domain
        map(domain->lagrange_space(domain;order=order_inter),rid_to_domain)
    end
    face_nodes_tuple = map(dims) do d
        domain = GT.domain(mesh,d)
        rid_to_refqua = map(reference_domains(domain)) do refdom
            reffe = lagrange_space(refdom;order=order_inter)
            node_quadrature(reffe)
        end
        face_to_rid = face_reference_id[d+1]
        quadrature = GT.mesh_quadrature(;
            domain,
            reference_quadratures=rid_to_refqua,
            face_reference_id = face_to_rid
           )
        face_point_x = coordinate_accessor(quadrature)
        face_npoints = num_points_accessor(quadrature)
        nfaces = num_faces(mesh,d)
        dface_to_nodes = Vector{Vector{Ti}}(undef,nfaces)
        for face in 1:nfaces
            npoints = face_npoints(face)
            point_x = face_point_x(face)
            x_mapped = map(1:npoints) do point
                x = point_x(point)
                round_coord(x)
            end
            my_nodes = indexin(x_mapped,node_coordinates_mapped)
            dface_to_nodes[face] = my_nodes
        end
        dface_to_nodes
    end
    face_nodes = collect(face_nodes_tuple)
    outward_normals = GT.outward_normals(mesh)
    GT.mesh(;
        node_coordinates,
        face_nodes,
        face_reference_id,
        reference_spaces,
        outward_normals
       )
end


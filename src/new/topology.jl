
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

function complexify(mesh::AbstractMesh;glue=Val(false))
    Ti = int_type(options(mesh))
    T = JaggedArray{Ti,Ti}
    D = num_dims(mesh)
    oldface_to_newvertices = Vector{T}(undef,D+1)
    newvertex_to_oldfaces = Vector{T}(undef,D+1)
    newface_incidence = Matrix{T}(undef,D+1,D+1)
    nnewfaces = zeros(Int,D+1)
    newface_refid = Vector{Vector{Ti}}(undef,D+1)
    newreffaces = Vector{Any}(undef,D+1)
    newface_nodes = Vector{T}(undef,D+1)
    old_to_new = Vector{Vector{Ti}}(undef,D+1)
    node_to_newvertex, n_new_vertices = find_node_to_vertex(mesh) # Optimizable for linear meshes
    for d in 0:D
        oldface_to_newvertices[d+1] = fill_face_vertices(mesh,d,node_to_newvertex) # Optimizable for linear meshes
        newvertex_to_oldfaces[d+1] = generate_face_coboundary(oldface_to_newvertices[d+1],n_new_vertices) # Optimizable for linear meshes
    end
    newface_incidence[D+1,0+1] = oldface_to_newvertices[D+1]
    newface_incidence[0+1,D+1] = newvertex_to_oldfaces[D+1]
    nnewfaces[D+1] = length(oldface_to_newvertices[D+1])
    newface_refid[D+1] = face_reference_id(mesh,D)
    newreffaces[D+1] = reference_faces(mesh,D)
    newface_nodes[D+1] = face_nodes(mesh,D)
    old_to_new[D+1] = collect(Ti,1:length(newface_nodes[D+1]))
    # TODO optimize for d==0
    for d in (D-1):-1:0
        n = d+1
        new_nface_to_new_vertices = newface_incidence[n+1,0+1]
        new_vertex_to_new_nfaces = newface_incidence[0+1,n+1]
        old_dface_to_new_vertices = oldface_to_newvertices[d+1]
        new_vertex_to_old_dfaces = newvertex_to_oldfaces[d+1]
        new_nface_to_nrefid = newface_refid[n+1]
        old_dface_to_drefid = face_reference_id(mesh,d)
        drefid_to_ref_dface = reference_faces(mesh,d)
        old_dface_to_nodes = face_nodes(mesh,d)
        new_nface_to_nodes = newface_nodes[n+1]
        nrefid_to_ldface_to_lvertices = map(a->face_incidence(topology(boundary(geometry(a))),d,0),newreffaces[n+1])
        nrefid_to_ldface_to_lnodes = map(a->face_nodes(boundary(a),d),newreffaces[n+1])
        nrefid_to_ldface_to_drefrefid = map(a->face_reference_id(boundary(a),d),newreffaces[n+1])
        nrefid_to_drefrefid_to_ref_dface = map(a->reference_faces(boundary(a),d),newreffaces[n+1])
        new_nface_to_new_dfaces, n_new_dfaces, old_dface_to_new_dface = generate_face_boundary(
            new_nface_to_new_vertices,
            new_vertex_to_new_nfaces,
            old_dface_to_new_vertices,
            new_vertex_to_old_dfaces,
            new_nface_to_nrefid,
            nrefid_to_ldface_to_lvertices)
        new_dface_to_new_vertices = generate_face_vertices(
            n_new_dfaces,
            old_dface_to_new_dface,
            old_dface_to_new_vertices,
            new_nface_to_new_vertices,
            new_nface_to_new_dfaces,
            new_nface_to_nrefid,
            nrefid_to_ldface_to_lvertices)
        new_vertex_to_new_dfaces = generate_face_coboundary(new_dface_to_new_vertices,n_new_vertices)
        new_dface_to_nodes = generate_face_vertices(
            n_new_dfaces,
            old_dface_to_new_dface,
            old_dface_to_nodes,
            new_nface_to_nodes,
            new_nface_to_new_dfaces,
            new_nface_to_nrefid,
            nrefid_to_ldface_to_lnodes)
        new_dface_to_new_drefid, new_refid_to_ref_dface = generate_reference_faces(
            n_new_dfaces,
            old_dface_to_new_dface,
            old_dface_to_drefid,
            drefid_to_ref_dface,
            new_nface_to_new_dfaces,
            new_nface_to_nrefid,
            nrefid_to_ldface_to_drefrefid,
            nrefid_to_drefrefid_to_ref_dface)
        newface_incidence[n+1,d+1] = new_nface_to_new_dfaces
        newface_incidence[d+1,0+1] = new_dface_to_new_vertices
        newface_incidence[0+1,d+1] = new_vertex_to_new_dfaces
        newface_refid[d+1] = new_dface_to_new_drefid
        newreffaces[d+1] = new_refid_to_ref_dface
        nnewfaces[d+1] = n_new_dfaces
        newface_nodes[d+1] = new_dface_to_nodes
        old_to_new[d+1] = old_dface_to_new_dface
    end
    node_to_coords = node_coordinates(mesh)
    old_physical_faces = physical_faces(mesh)
    new_physical_faces = [ Dict{String,Vector{Int32}}() for d in 0:D] # TODO hardcoded
    for d in 0:D
        old_groups = old_physical_faces[d+1]
        for (group_name,old_group_faces) in old_groups
            new_group_faces = similar(old_group_faces)
            new_group_faces .= old_to_new[d+1][old_group_faces]
            new_physical_faces[d+1][group_name] = new_group_faces
        end
    end
    new_mesh = GT.mesh_from_arrays(
            node_to_coords,
            newface_nodes,
            newface_refid,
            Tuple(newreffaces);
            physical_faces = new_physical_faces,
            periodic_nodes = periodic_nodes(mesh),
            outwards_normals = outwards_normals(mesh)
           )
    if val_parameter(glue)
        new_mesh, old_to_new
    else
        new_mesh
    end
end


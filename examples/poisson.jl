import GalerkinToolkit as glk
import ForwardDiff
using StaticArrays
using LinearAlgebra
using SparseArrays
using WriteVTK

function add_default_params(params_in,params_default)
    UmD = setdiff(keys(params_in),keys(params_default))
    for k in UmD
        @warn "Parameter with key :$k is unused"
    end
    UiD = intersect(keys(params_default),keys(params_in))
    DmU = setdiff(keys(params_default),keys(params_in))
    a = [ k=>params_in[k] for k in UiD]
    b = [ k=>params_default[k] for k in DmU]
    Dict(vcat(a,b))
end

function default_params_poisson()
    params = Dict{Symbol,Any}()
    params[:mesh] = glk.cartesian_mesh((0,1,0,1),(10,10))
    params[:u] = (x) -> sum(x)
    params[:f] = (x) -> 0.0
    params[:g] = (x) -> 0.0
    params[:dirichlet_tags] = ["boundary"]
    params[:neumann_tags] = String[]
    params[:integration_degree] = 2
    params[:solve] = \
    params[:example_name] = "poisson"
    params
end

function poisson(params_in)

    # Process params
    params_default = default_params_poisson()
    params = add_default_params(params_in,params_default)

    # Get the parameters
    mesh = params[:mesh]
    u = params[:u]
    f = params[:f]
    g = params[:g]
    dirichlet_tags = params[:dirichlet_tags]
    degree = params[:integration_degree]
    solve = params[:solve]
    example_name = params[:example_name]
    neumann_tags = params[:neumann_tags]

    # Dirichlet values
    node_to_tag = zeros(glk.num_nodes(mesh))
    tag_to_name = dirichlet_tags
    glk.classify_mesh_nodes!(node_to_tag,mesh,tag_to_name)
    free_and_dirichlet_nodes = glk.partition_from_mask(i->i==0,node_to_tag)
    node_to_x = glk.node_coordinates(mesh)
    dirichlet_nodes = last(free_and_dirichlet_nodes)
    x_dirichlet = view(node_to_x,dirichlet_nodes)
    u_dirichlet = u.(x_dirichlet)

    # Reference cell
    d = glk.num_dims(mesh)
    ref_cells = glk.reference_faces(mesh,d)
    cell_to_rid = glk.face_reference_id(mesh,d)

    # Integration rule
    integration_rules = map(ref_cells) do ref_cell
        glk.quadrature(glk.geometry(ref_cell),degree)
    end
    w = map(glk.weights,integration_rules)
    q = map(glk.coordinates,integration_rules)

    # Shape functions
    ab = map(integration_rules,ref_cells) do integration_rule,ref_cell
        q = glk.coordinates(integration_rule)
        shape_functions = glk.shape_functions(ref_cell)
        m = length(q)
        n = glk.num_nodes(ref_cell)
        a1 = glk.tabulation_matrix!(shape_functions)(ForwardDiff.gradient,zeros(eltype(q),m,n),q)
        b1 = glk.tabulation_matrix!(shape_functions)(glk.value,zeros(m,n),q)
        (a1,b1)
    end
    ∇s = map(first,ab)
    s = map(last,ab)

    # Count coo entries
    n_coo = 0
    cell_to_nodes = glk.face_nodes(mesh,d)
    ncells = glk.num_faces(mesh,d)
    node_to_free_node = glk.permutation(free_and_dirichlet_nodes)
    nfree = length(first(free_and_dirichlet_nodes))
    for cell in 1:ncells
        nodes = cell_to_nodes[cell]
        for nj in nodes
            if node_to_free_node[nj] > nfree
                continue
            end
            for ni in nodes
                if node_to_free_node[ni] > nfree
                    continue
                end
                n_coo += 1
            end
        end
    end

    # Allocate coo values
    I_coo = zeros(Int32,n_coo)
    J_coo = zeros(Int32,n_coo)
    V_coo = zeros(Float64,n_coo)
    b = zeros(nfree)

    # Allocate auxiliary buffers
    ∇x = map(i->similar(i,size(i,2)),∇s)
    Aes = map(i->zeros(size(i,2),size(i,2)),∇s)
    fes = map(i->zeros(size(i,2)),∇s)
    ues = map(i->zeros(size(i,2)),∇s)

    Tx = SVector{d,Float64}
    TJ = typeof(zero(Tx)*zero(Tx)')

    # Fill coo values
    i_coo = 0
    for cell in 1:ncells
        rid = cell_to_rid[cell]
        nodes = cell_to_nodes[cell]
        Ae = Aes[rid]
        ue = ues[rid]
        fe = fes[rid]
        nl = length(nodes)
        ∇se = ∇s[rid]
        ∇xe = ∇x[rid]
        se = s[rid]
        we = w[rid]
        nq = length(we)
        fill!(Ae,zero(eltype(Ae)))
        fill!(fe,zero(eltype(Ae)))
        fill!(ue,zero(eltype(Ae)))
        for k in 1:nl
            nk = nodes[k]
            gk = node_to_free_node[nk]
            if  gk <= nfree
                continue
            end
            ue[k] = u_dirichlet[gk-nfree]
        end
        for iq in 1:nq
            Jt = zero(TJ) 
            xint = zero(Tx)
            for k in 1:nl
                x = node_to_x[nodes[k]]
                ∇sqx = ∇se[iq,k]
                sqx = se[iq,k]
                Jt += ∇sqx*x'
                xint += sqx*x
            end
            detJt = det(Jt)
            invJt = inv(Jt)
            dV = abs(detJt)*we[iq]
            for k in 1:nl
                ∇xe[k] = invJt*∇se[iq,k]
            end
            for j in 1:nl
                for i in 1:nl
                    Ae[i,j] += ∇xe[i]⋅∇xe[j]*dV
                end
            end
            fx = f(xint)
            for k in 1:nl
                fe[k] +=  fx*se[iq,k]*dV
            end
        end
        for i in 1:nl
            for j in 1:nl
                fe[i] -= Ae[i,j]*ue[j]
            end
        end
        for i in 1:nl
            ni = nodes[i]
            gi = node_to_free_node[ni]
            if gi > nfree
                continue
            end
            b[gi] += fe[i]
        end
        for j in 1:nl
            nj = nodes[j]
            gj = node_to_free_node[nj]
            if  gj > nfree
                continue
            end
            for i in 1:nl
                ni = nodes[i]
                gi = node_to_free_node[ni]
                if gi > nfree
                    continue
                end
                i_coo += 1
                I_coo[i_coo] = gi
                J_coo[i_coo] = gj
                V_coo[i_coo] = Ae[i,j]
            end
        end
    end

    # Neumann conditions
    if length(neumann_tags) != 0

        face_to_rid = glk.face_reference_id(mesh,d-1)
        face_to_nodes = glk.face_nodes(mesh,d-1)
        ref_faces = glk.reference_faces(mesh,d-1)
        n_ref_faces = length(ref_faces)

        # Integration rule
        integration_rules_f = map(ref_faces) do ref_face
            glk.quadrature(glk.geometry(ref_face),degree)
        end
        w_f = map(glk.weights,integration_rules_f)
        q_f = map(glk.coordinates,integration_rules_f)

        # Shape functions
        ab = map(integration_rules_f,ref_faces) do integration_rule,ref_cell
            q = glk.coordinates(integration_rule)
            shape_functions = glk.shape_functions(ref_cell)
            m = length(q)
            n = glk.num_nodes(ref_cell)
            a1 = glk.tabulation_matrix!(shape_functions)(ForwardDiff.gradient,zeros(eltype(q),m,n),q)
            b1 = glk.tabulation_matrix!(shape_functions)(glk.value,zeros(m,n),q)
            (a1,b1)
        end
        ∇s_f = map(first,ab)
        s_f = map(last,ab)

        tag_to_groups = glk.physical_groups(mesh,d-1)
        neum_face_to_face = reduce(union,map(tag->tag_to_groups[tag],neumann_tags))
        fes = map(i->zeros(size(i,2)),s_f)

        Ts = typeof(zero(SVector{d-1,Float64})*zero(Tx)')

        for face in neum_face_to_face
            rid = face_to_rid[face]
            nodes = face_to_nodes[face]
            ∇se = ∇s_f[rid]
            se = s_f[rid]
            we = w_f[rid]
            nq = length(we)
            fe = fes[rid]
            fill!(fe,zero(eltype(fe)))
            nl = length(nodes)
            for iq in 1:nq
                Jt = zero(Ts) 
                xint = zero(Tx)
                for k in 1:nl
                    x = node_to_x[nodes[k]]
                    ∇sqx = ∇se[iq,k]
                    sqx = se[iq,k]
                    Jt += ∇sqx*x'
                    xint += sqx*x
                end
                J = transpose(Jt)
                dS = sqrt(det(Jt*J))*we[iq]
                ∇ux = ForwardDiff.gradient(u,xint)
                gx = g(xint)
                for k in 1:nl
                    fe[k] +=  gx*se[iq,k]*dS
                end
            end
            for i in 1:nl
                ni = nodes[i]
                gi = node_to_free_node[ni]
                if gi > nfree
                    continue
                end
                b[gi] += fe[i]
            end
        end
    end

    A = sparse(I_coo,J_coo,V_coo,nfree,nfree)
    u_free = solve(A,b)
    nnodes = glk.num_nodes(mesh)
    uh = zeros(nnodes)
    uh[first(free_and_dirichlet_nodes)] = u_free
    uh[last(free_and_dirichlet_nodes)] = u_dirichlet

    eh1 = 0.0
    el2 = 0.0
    for cell in 1:ncells
        rid = cell_to_rid[cell]
        nodes = cell_to_nodes[cell]
        ue = ues[rid]
        nl = length(nodes)
        ∇se = ∇s[rid]
        ∇xe = ∇x[rid]
        se = s[rid]
        we = w[rid]
        nq = length(we)
        for k in 1:nl
            nk = nodes[k]
            ue[k] = uh[nk]
        end
        for iq in 1:nq
            Jt = zero(TJ) 
            xint = zero(Tx)
            for k in 1:nl
                x = node_to_x[nodes[k]]
                ∇sqx = ∇se[iq,k]
                sqx = se[iq,k]
                Jt += ∇sqx*x'
                xint += sqx*x
            end
            detJt = det(Jt)
            invJt = inv(Jt)
            dV = abs(detJt)*we[iq]
            for k in 1:nl
                ∇xe[k] = invJt*∇se[iq,k]
            end
            ux = u(xint)
            ∇ux = ForwardDiff.gradient(u,xint)
            ∇uhx = zero(∇ux)
            uhx = zero(ux)
            for k in 1:nl
                uek = ue[k]
                ∇uhx += uek*∇xe[k]
                uhx += uek*se[iq,k]
            end
            ∇ex = ∇ux - ∇uhx
            ex =  ux - uhx
            eh1 += (∇ex⋅∇ex + ex*ex)*dV
            el2 += ex*ex*dV
        end
    end

    results = Dict{Symbol,Any}()
    eh1 = sqrt(eh1)
    el2 = sqrt(el2)
    results[:eh1] = eh1
    results[:el2] = el2
    results[:nnodes] = nnodes
    results[:nfree] = nfree
    results[:ncells] = ncells

    vtk_grid(example_name,glk.vtk_args(mesh)...) do vtk
        glk.vtk_physical_groups!(vtk,mesh)
        vtk["tag"] = node_to_tag
        vtk["uh"] = uh
        vtk["dim"] = glk.face_dim(mesh)
    end
    results
end




function quantity(term)
    Quantity(term)
end

struct Quantity{A} <: AbstractQuantity
    term::A
end

term(a::Quantity) = a.term

function index(;
    domain = nothing,
    face=nothing,
    point=nothing,
    field=nothing,
    dof=nothing,
    face_around=nothing,
    face_around_dummy=nothing,
    prefix = gensym,
    dict=IdDict{Any,Symbol}(),
    )
    data = (;domain,face,point,field,dof,face_around,face_around_dummy,dict,prefix)
    Index(data)
end

struct Index{A} <: GT.AbstractType
    data::A
end

function generate_index(dom::AbstractDomain,form_arity=0;
    prefix=gensym,
    field=[prefix("field-in-axis-$i") for i in 1:form_arity],
    face_around = [prefix("face-around-in-axis-$a") for a in 1:form_arity],
    )
    d = num_dims(dom)
    D = num_dims(mesh(dom))
    face = [ prefix("$d-face-dummy") for d in 0:D ]
    point = prefix("point")
    dof = [prefix("dof-in-axis-$i") for i in 1:form_arity]
    face[d+1] = prefix("$d-face")
    face_around_dummy = [ prefix("$(d_around)face-around-$(d)face-dummy") for d in 0:D, d_around in 0:D ] 
    index(;domain=dom,face,point,field,dof,face_around,face_around_dummy,prefix)
end

domain(index::Index) = index.data.domain

function num_dims(index::Index)
    num_dims(index.data.domain)
end

function mesh(index::Index)
    index.data.domain.mesh
end

target_dim(index::Index) = num_dims(index.data.domain)

function form_arity(index::Index)
    length(face_around_index(index))
end

function index_storage(index)
    (;( key=>val for (val,key) in index.data.dict )...)
end

function unpack_index_storage(index,state)
    expr = Expr(:block)
    for k in Base.values(index.data.dict) |> collect |> sort
        push!(expr.args,:($k = $state.$k))
    end
    expr
end

function get_symbol!(index,val,name="";prefix=index.data.prefix)
    get!(index.data.dict,val,prefix(name))
end

function face_index(index,d)
    index.data.face[d+1]
end

max_dim(index::Index) = length(index.data.face)-1

function point_index(index)
    index.data.point
end

function face_around_term(index,d,D)
    domain = index.data.domain
    if GT.faces_around(domain) === nothing
        return nothing
    end
    if d == num_dims(index) && (d+1 == D)
        #NB this is not constant anymore
        #constant_term(first(GT.faces_around(domain)),index;compile_constant=true)
        faces_around = get_symbol!(index,GT.faces_around(domain),"faces_around")
        face_to_sface = inverse_faces(dom)
        face_to_sface_sym = get_symbol!(index,face_to_sface,"face_to_sface")
        expr = @term $faces_around[$face_to_sface_sym[$face]]
        p=one(eltype(faces_around))
        expr_term([d],expr,p,index)
    else
        nothing
    end
end

function face_around_index(index,a)
    index.data.face_around[a]
end

function face_around_index(index)
    index.data.face_around
end

function face_around_dummy_index(index,d,D)
    index.data.face_around_dummy[d+1,D+1]
end

function field_index(index,a)
    index.data.field[a]
end

function dof_index(index,a)
    index.data.dof[a]
end

function constant_quantity(v;kwargs...)
    quantity() do index
        constant_term(v,index;kwargs...)
        #expr = get_symbol!(index,v,"constant_quantity_value")
        #dim = ANY_DIM
        #(;expr,dim,prototype=v)
    end
end

function face_quantity(data,mesh::AbstractMesh,d;reference=Val(false))
    dim = d
    quantity() do index
        face = face_index(index,dim)
        p = zero(eltype(data))
        if val_parameter(reference)
            #reference_face_term(dim,rid_to_value,index)
            face_to_rid = get_symbol!(index,face_reference_id(mesh,dim),"face_to_rid")
            rid_to_value = get_symbol!(index,data,"rid_to_value")
            expr = @term $rid_to_value[$face_to_refid[$face]]
            expr_term([dim],expr,p,index)
        else
            #sface_to_value = data
            #face_to_sface = inverse_faces(dom)
            #physical_face_term(dim,sface_to_value,face_to_sface,index)
            data_sym = get_symbol!(index,data,"face_to_value_$(dim)d")
            expr = @term $data_sym[$face]
            expr_term([dim],expr,p,index)
        end
        #(;expr,dim,prototype=p)
    end
end

function face_quantity(data,dom::AbstractDomain;reference=Val(false))
    dim = num_dims(dom)
    quantity() do index
        face = face_index(index,dim)
        p = zero(eltype(data))
        if val_parameter(reference)
            #reference_face_term(dim,rid_to_value,index)
            face_to_rid = get_symbol!(index,face_reference_id(mesh(dom),dim),"face_to_rid")
            rid_to_value = get_symbol!(index,data,"rid_to_value")
            expr = @term $rid_to_value[$face_to_refid[$face]]
            expr_term([dim],expr,p,index)
        else
            #sface_to_value = data
            #face_to_sface = inverse_faces(dom)
            #physical_face_term(dim,sface_to_value,face_to_sface,index)
            face_to_sface = inverse_faces(dom)
            data_sym = get_symbol!(index,data,"face_to_value_$(dim)d")
            face_to_sface_sym = get_symbol!(index,face_to_sface,"face_to_sface_$(dim)d")
            expr = @term $data_sym[$face_to_sface_sym[$face]]
            expr_term([dim],expr,p,index)
        end
        #(;expr,dim,prototype=p)
    end
end

function term(a::AbstractQuantity,index)
    a.term(index)
end

function return_prototype(f,args...)
    f(args...)
end

function return_prototype(::typeof(getindex),a,b)
    zero(eltype(a))
end

function return_prototype(::typeof(getindex),a,b::Integer)
    if length(a) != 0
        first(a)
    else
        zero(eltype(a))
    end
end

function call(f,args...)
    f(args...)
end

function get_symbol!(index,val::typeof(call),name="";prefix=index.data.prefix)
    :(GalerkinToolkit.call)
end

function call(g,q::AbstractQuantity)
    quantity() do index
        unary_call_term(g,term(q,index))
        #g_sym = get_symbol!(index,g,"callee")
        #tq = term(q,index)
        #p = GT.return_prototype(g,prototype(tq))
        #expr = :($g_sym($(tq.expr)))
        #(;expr,dim=tq.dim,prototype=p)
    end
end

function call(g,a::AbstractQuantity,b::AbstractQuantity)
    quantity() do index
        ta = term(a,index)
        tb = term(b,index)
        t = binary_call_term(g,term(a,index),term(b,index))
        dims = union(free_dims(ta),free_dims(tb))
        if length(dims) <= 1
            return t
        end
        @assert length(dims) == 2
        dim = minimum(dims)
        dim2 = maximum(dims)
        cell_around = face_around_term(index,dim,dim2)
        if cell_around !== nothing
            boundary_term(dim,dim2,t,cell_around)
        else
            skeleton_term(dim,dim2,t)
        end
        #g_sym = get_symbol!(index,g,"callee")
        #ta = term(a,index)
        #tb = term(b,index)
        #p = return_prototype(g,prototype(ta),prototype(tb))
        #expr = :($g_sym($(ta.expr),$(tb.expr)))
        #dim = promote_dim(ta.dim,tb.dim)
        #p2 = p
        #if dim == nothing
        #    dim = min(ta.dim,tb.dim)
        #    dim2 = max(ta.dim,tb.dim)
        #    dummy = face_index(index,dim2)
        #    face = face_index(index,dim)
        #    #faces = actual_faces_index(index,dim,dim2)
        #    topo = topology(mesh(index))
        #    face_to_cells = get_symbol!(index,face_incidence(topo,dim,dim2),"face_to_cells")
        #    axis_to_face_around = face_around_index(index)
        #    #dom = target_domain(index)
        #    #gluable = faces_around(dom) !== nothing && num_dims(dom) == dim && dim2 == dim+1
        #    #if gluable
        #    #    face = face(index,dim)
        #    #    face_to_sface = get_symbol!(index,inverse_faces(dom))
        #    #    sface_to_cell_around = get_symbol!(index,faces_around(dom))
        #    #    topo = topology(mesh(dom))
        #    #    face_to_cells = get_symbol!(index,face_incidence(topo,dim,dim2))
        #    #    expr = @term begin
        #    #        sface = $face_to_sface[$face]
        #    #        cell_around = $face_to_cells[$sface]


        #    #        fun = $dummy -> $expr
        #    #        fun()
        #    #    end
        #    cells = @term $face_to_cells[$face]
        #    cell_around = face_around_term(index,dim,dim2)
        #    cell_around_sym = get_symbol!(index,cell_around,"cell_around")
        #    if isa(cell_around,Int)
        #        #expr = @term GalerkinToolkit.call($dummy -> $expr,$cells[$cell_around_sym])
        #        actual = :($cells[$cell_around_sym])
        #        expr = substitute(expr,dummy=>actual)
        #        p2 = p
        #    elseif isa(cell_around,AbstractArray)
        #        face_to_sface = get_symbol!(inverse_faces(domain(index)),"face_to_sface")
        #        actual = :($cells[$cell_around_sym[$face_to_sface[$face]]])
        #        expr = substitute(expr,dummy=>actual)
        #        #expr = @term GalerkinToolkit.call($dummy -> $expr,$cells[$cell_around_sym[$face_to_sface[$face]]])
        #        p2 = p
        #    elseif cell_around !== nothing
        #        error()
        #    elseif length(axis_to_face_around) == 0
        #        expr = @term begin
        #            fun = $dummy -> $expr
        #            map(fun,$cells)
        #        end
        #        p2 = [p,]
        #    else
        #        i = gensym("i-dummy")
        #        exprs = map(axis_to_face_around) do face_around
        #            :(LinearAlgebra.I[$i,$face_around])
        #        end
        #        if length(exprs) > 1
        #            deltas = Expr(:call,:tuple,exprs...)
        #            mask = :(prod($deltas))
        #        else
        #            mask = exprs[1]
        #        end
        #        expr = @term begin
        #            fun = ($i,$dummy) -> $expr*$mask
        #            map(fun,1:length($cells),$cells)
        #        end
        #        p2 = [p,]
        #    end
        #end
        #(;expr,dim,prototype=p2)
    end
end

function call(g,args::AbstractQuantity...)
    quantity() do index
        multivar_call_term(g,map(a->term(a,index),args))
        #g_sym = get_symbol!(index,v,gensym("callee"))
        #ts = map(a->term(a,index),args)
        #p = GT.return_prototype(g,(map(GT.prototype,ts)...))
        #dims = map(t->ts.dim,ts)
        #dim = promote_dim(dims...)
        #@assert dim !== nothing
        #exprs = map(t->ts.expr,ts)
        ## We use call since @rule is not able to
        ## identify function calls on variable slots
        #expr = @term call($g_sym,$(exprs...))
        #(;expr,dim,prototype=p)
    end
end

function (f::AbstractQuantity)(x::AbstractQuantity)
    call(call,f,x)
end

function physical_map(mesh::AbstractMesh,d)
    quantity() do index
        physical_map_term(d,index)
    end
end

function physical_map(mesh::AbstractMeshDomain{<:PMesh},d)
    quantity() do index
        dval = Val(val_parameter(d))
        map(x -> GT.physical_map_term(dval, index),partition(mesh))
    end
end

function inverse_physical_map(mesh::AbstractMesh,d)
    x0 = zero(SVector{val_parameter(d),Float64})
    x = constant_quantity(x0)
    phi = physical_map(mesh,d)
    call(inv_map,phi,x)
end

function point_quantity(data,dom::AbstractDomain;
    reference=Val(false),
    face_reference_id = face_reference_id(mesh(dom),num_dims(dom)),
    )
    dim = num_dims(dom)
    #p = zero(eltype(eltype(data)))
    quantity() do index
        if val_parameter(reference)
            reference_point_term(dim,data,face_reference_id,index)
        else
            physical_point_term(dim,data,inverse_faces(dom),index)
        end
        #face = face_index(index,dim)
        #point = point_index(index)
        #expr = if val_parameter(reference)
        #    face_to_rid = get_symbol!(index,face_reference_id(mesh(dom),dim),"face_to_rid")
        #    rid_to_points = get_symbol!(index,data,"rid_to_points")
        #    @term begin
        #        rid = $face_to_rid[$face]
        #        $rid_to_points[rid][$point]
        #    end
        #else
        #    face_to_sface = get_symbol!(index,inverse_faces(dom),"face_to_sface")
        #    sface_to_points = get_symbol!(index,data,"sface_to_points")
        #    @term begin
        #        sface = $face_to_sface[$face]
        #        $sface_to_points[sface][$point]
        #    end
        #end
        #(;expr,dim,prototype=p)
    end
end

function shape_function_quantity(data,dom::AbstractDomain;
    reference=Val(false),
    face_reference_id = face_reference_id(mesh(dom),num_dims(dom)),
    dof = gensym("dummy-dof"),
    )
    dim = num_dims(dom)
    #p = zero(eltype(eltype(data)))
    quantity() do index
        if val_parameter(reference)
            reference_shape_function_term(dim,data,face_reference_id,dof,index)
        else
            error("Not implemented, but possible to implement it.")
        end
    end
end

function form_argument(axis,field,data,dom::AbstractDomain;
        reference=Val(false),
        face_reference_id = face_reference_id(mesh(dom),num_dims(dom)),
    )

    dim = num_dims(dom)
    quantity() do index
        dof = dof_index(index,axis)
        s = if val_parameter(reference)
            reference_shape_function_term(dim,data,face_reference_id,dof,index)
        else
            error("Not implemented, but possible to implement it.")
        end
        evaluated = false
        form_argument_term(axis,field,s,evaluated)
    end
end

function reference_map(mesh::AbstractMesh,d,D)
    quantity() do index
        t = reference_map_term(d,D,index)
        #cell_around = face_around_term(index,d,D)
        #if cell_around !== nothing
        #    boundary_term(d,D,t,cell_around)
        #else
        #    skeleton_term(d,D,t)
        #end
    end
end

function face_node_coordinates(refDface,d)
    D = num_dims(refDface)
    if d == D
        [[node_coordinates(refDface)]]
    else
        boundary = refDface |> GT.domain |> GT.mesh
        lface_to_nodes = GT.face_nodes(boundary,d)
        node_to_coords = GT.node_coordinates(boundary)
        lface_to_lrefid = GT.face_reference_id(boundary,d)
        lrefid_to_lrefface = GT.reference_spaces(boundary,d)
        lrefid_to_perm_to_ids = map(GT.node_permutations,lrefid_to_lrefface)
        map(1:GT.num_faces(boundary,d)) do lface
            lrefid = lface_to_lrefid[lface]
            nodes = lface_to_nodes[lface]
            perm_to_ids = lrefid_to_perm_to_ids[lrefid]
            map(perm_to_ids) do ids
                coords = node_to_coords[nodes[ids]]
                coords
            end
        end
    end
end

function reference_map(refdface::AbstractFaceSpace,refDface::AbstractFaceSpace)
    d = num_dims(refdface)
    dof_to_f = shape_functions(refdface)
    boundary = refDface |> GT.domain |> GT.mesh
    lface_to_nodes = GT.face_nodes(boundary,d)
    node_to_coords = GT.node_coordinates(boundary)
    lface_to_lrefid = GT.face_reference_id(boundary,d)
    lrefid_to_lrefface = GT.reference_spaces(boundary,d)
    lrefid_to_perm_to_ids = map(GT.node_permutations,lrefid_to_lrefface)
    map(1:GT.num_faces(boundary,d)) do lface
        lrefid = lface_to_lrefid[lface]
        nodes = lface_to_nodes[lface]
        perm_to_ids = lrefid_to_perm_to_ids[lrefid]
        map(perm_to_ids) do ids
            dof_to_coeff = node_to_coords[nodes[ids]]
            ndofs = length(dof_to_coeff)
            x -> sum(dof->dof_to_coeff[dof]*dof_to_f[dof](x),1:ndofs)
        end
    end
end

#function face_local_map(refDface,refDface)
#end

function unit_normal(mesh::AbstractMesh,d)
    D = num_dims(mesh)
    @assert d == D-1
    phiDinv = inverse_physical_map(mesh,D)
    ctype_to_refface = GT.reference_spaces(mesh,D)
    Drid_to_lface_to_n_data= map(ctype_to_refface) do refface
        boundary = refface |> GT.domain |> GT.mesh
        boundary |> GT.normals # TODO also rename?
    end
    topo = topology(mesh)
    dface_to_Dfaces_data, dface_to_ldfaces_data = GT.face_incidence_ext(topo,d,D)
    quantity() do index
        dface = face_index(index,d)
        Dface = face_index(index,D)
        Dface_to_Drid = get_symbol!(index,face_reference_id(mesh,D),"Dface_to_Drid")
        dface_to_Dfaces = get_symbol!(index,dface_to_Dfaces_data,"dface_to_Dfaces")
        dface_to_ldfaces = get_symbol!(index,dface_to_ldfaces_data,"dface_to_ldfaces")
        Drid_to_lface_to_n = get_symbol!(index,Drid_to_lface_to_n_data,"Drid_to_lface_to_n")
        Dface_around = face_around_dummy_index(index,d,D)
        expr = @term begin
            Drid = $Dface_to_Drid[$Dface]
            Dfaces = $dface_to_Dfaces[$dface]
            ldface = $dface_to_ldfaces[$dface][$Dface_around]
            $Drid_to_lface_to_n[Drid][ldface]
        end
        p = Drid_to_lface_to_n_data |> eltype |> eltype |> zero
        n = expr_term([d,D],expr,p,index)
        n_ref = unit_normal_term(n)
        ϕinv_fun = term(phiDinv,index)
        n_phys = binary_call_term(∘,n_ref,ϕinv_fun)
        n_phys
    end
    #phiD = physical_map(mesh,D)
    #phiDinv = inverse_physical_map(mesh,D)
    #ctype_to_refface = GT.reference_spaces(mesh,D)
    #Drid_to_lface_to_n_data= map(ctype_to_refface) do refface
    #    boundary = refface |> GT.geometry |> GT.boundary
    #    boundary |> GT.normals # TODO also rename?
    #end
    #topo = topology(mesh)
    #dface_to_Dfaces_data, dface_to_ldfaces_data = GT.face_incidence_ext(topo,d,D)
    #quantity() do index
    #    dface = face_index(index,d)
    #    Dface = face_index(index,D)
    #    Dface_to_Drid = get_symbol!(index,face_reference_id(mesh,D),"Dface_to_Drid")
    #    dface_to_Dfaces = get_symbol!(index,dface_to_Dfaces_data,"dface_to_Dfaces")
    #    dface_to_ldfaces = get_symbol!(index,dface_to_ldfaces_data,"dface_to_ldfaces")
    #    Drid_to_lface_to_n = get_symbol!(index,Drid_to_lface_to_n_data,"Drid_to_lface_to_n")
    #    expr = @term begin
    #        Drid = $Dface_to_Drid[$Dface]
    #        Dfaces = $dface_to_Dfaces[$dface]
    #        Dface_around = GalerkinToolkit.find_face_adound($Dface,Dfaces)
    #        ldface = $dface_to_ldfaces[$dface][Dface_around]
    #        Drid_to_lface_to_n[Drif][ldface]
    #    end
    #    p = Drid_to_lface_to_n_data |> eltype |> eltype |> zero
    #    n = expr_term([d,D],expr,p,index)
    #    ϕ_fun = term(phiD,index)
    #    

    #    ϕinv_fun = term(phiDinv,index)
    #    n_ref = binary_call_term(map_unit_normal,ϕ_fun,n)
    #    n_phys = binary_call_term(∘,n_ref,ϕinv_fun)
    #    n_phys
    #end
end

function Base.:∘(a::AbstractQuantity,phi::AbstractQuantity)
    call(∘,a,phi)
end

function get_symbol!(index,val::typeof(Base.:∘),name="";prefix=index.data.prefix)
    :(Base.:∘)
end

# Operations

# Base

for op in (:+,:-,:sqrt,:abs,:abs2,:real,:imag,:conj,:transpose,:adjoint,:*,:/,:\,:^,:getindex)
    @eval begin
        function get_symbol!(index,val::typeof(Base.$op),name="";prefix=gensym)
            $( Expr(:quote,op) )
        end
    end
end

for op in (:+,:-,:sqrt,:abs,:abs2,:real,:imag,:conj,:transpose,:adjoint)
  @eval begin
      (Base.$op)(a::AbstractQuantity) = call(Base.$op,a)
  end
end

function Base.getindex(a::AbstractQuantity,b::AbstractQuantity)
    call(getindex,a,b)
end

function Base.getindex(a::AbstractQuantity,b::Number)
    Base.getindex(a,constant_quantity(b;compile_constant=true))
end

for op in (:+,:-,:*,:/,:\,:^)
  @eval begin
      (Base.$op)(a::AbstractQuantity,b::AbstractQuantity) = call(Base.$op,a,b)
      (Base.$op)(a::Number,b::AbstractQuantity) = call(Base.$op,GT.constant_quantity(a),b)
      (Base.$op)(a::AbstractQuantity,b::Number) = call(Base.$op,a,GT.constant_quantity(b))
  end
end

# LinearAlgebra

for op in (:inv,:det,:norm,:tr)
  @eval begin
      function get_symbol!(index,val::typeof(LinearAlgebra.$op),name="";prefix=gensym)
          $( Expr(:quote,op) )
      end
    (LinearAlgebra.$op)(a::AbstractQuantity) = call(LinearAlgebra.$op,a)
  end
end

for op in (:dot,:cross)
  @eval begin
      function get_symbol!(index,val::typeof(LinearAlgebra.$op),name="";prefix=gensym)
          $( Expr(:quote,op) )
      end
      (LinearAlgebra.$op)(a::AbstractQuantity,b::AbstractQuantity) = call(LinearAlgebra.$op,a,b)
      (LinearAlgebra.$op)(a::Number,b::AbstractQuantity) = call(LinearAlgebra.$op,GT.constant_quantity(a),b)
      (LinearAlgebra.$op)(a::AbstractQuantity,b::Number) = call(LinearAlgebra.$op,a,GT.constant_quantity(b))
  end
end

# ForwardDiff


#function call_function_symbol(g,args::AbstractQuantity...)
#    fs = map(GT.term,args)
#    domain = args |> first |> GT.domain
#    prototype = GT.return_prototype(g,(map(GT.prototype,args)...))
#    g_expr = nameof(g)
#    GT.quantity(prototype,domain) do index
#        f_exprs = map(f->f(index),fs)
#        :($g_expr($(f_exprs...)))
#    end
#end

#function call_function_symbol(g,args::AbstractQuantity{<:PMesh}...)
#    pargs = map(partition,args)
#
#    q = map(pargs...) do myargs...
#        call_function_symbol(g, myargs...)
#    end
#    domain = args |> first |> GT.domain
#    term = map(GT.term,q)
#    prototype = map(GT.prototype,q) |> PartitionedArrays.getany
#    GT.quantity(term,prototype,domain)
#end

for op in (:(ForwardDiff.gradient),:(ForwardDiff.jacobian),:(ForwardDiff.hessian))
  @eval begin
      function get_symbol!(index,val::typeof($op),name="";prefix=gensym)
          $( Expr(:quote,op) )
      end
      ($op)(a::AbstractQuantity,b::AbstractQuantity) = call($op,a,b)
  end
end

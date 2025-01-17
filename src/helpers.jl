
"""
    abstract type AbstractType end

Parent of all types defined in GalerkinToolkit.
"""
abstract type AbstractType end

function Base.show(io::IO,data::GT.AbstractType)
    print(io,"GalerkinToolkit.$(nameof(typeof(data)))(…)")
end

"""
    push(a,ai)

Like `push!`, but creates a new object to store the result. 
This function is used to push to immutable collections such as tuples.
"""
function push end

function push(a::AbstractVector,x)
    b = copy(a)
    push!(b,x)
    b
end

function push(a::Tuple,x)
    (a...,x)
end

"""
    val_parameter(a)

For `a::Val{A}` it returns `A`. Otherwise, it returns `a`.
"""
val_parameter(a) = a
val_parameter(::Val{a}) where a = a

"""
    options(;kwargs...) -> Options

Create an object representing the default options for the current simulation.
This object can be used as an optional argument in several object constructors in GalerkinToolkit,
such as the mesh constructors `cartesian_mesh` and `mesh_from_gmsh`.
In this case, the computations using the generated mesh, will use the given options by default.
"""
function options(;
    reference_int_type=Int16,
    int_type=Int32,
    global_int_type=Int,
    real_type=Float64,
    )
    contents = (;
                reference_int_type=Val(reference_int_type),
                int_type=Val(int_type),
                global_int_type=Val(global_int_type),
                real_type=Val(real_type),
               )
    Options(contents)
end

options(object::AbstractType) = object.options

"""
    struct Options{...} <: AbstractType

Type of the objects returned by function `options`.
All properties and type parameters are private.

# Basic queries

- [`reference_int_type`](@ref)
- [`int_type`](@ref)
- [`global_int_type`](@ref)
- [`real_type`](@ref)
"""
struct Options{A} <: AbstractType
    contents::A
end

"""
    reference_int_type(options::Options)

Return the type of the integers used to enumerate reference quantities.
"""
reference_int_type(options::Options) = val_parameter(options.contents.reference_int_type)

"""
    int_type(options::Options)

Return the default integer type used in the computation except for reference and global quantities.
"""
int_type(options::Options) = val_parameter(options.contents.int_type)

"""
    global_int_type(options::Options)

Return the type of the integers used to enumerate global quantities.
"""
global_int_type(options::Options) = val_parameter(options.contents.global_int_type)

"""
    real_type(options::Options)

Return the default real type used in the computation.
"""
real_type(options::Options) = val_parameter(options.contents.real_type)


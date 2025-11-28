
"""
    with_gmsh(f[;options])

A safe way of initialize and finalize the `gmsh` module. The given function
is called `f(gmsh)` on the `gmsh` module after is has been initialized. The module is finalized automatically when the function returns.

The optional keyword argument `options` is a vector for pairs `k=>v` containing gmesh options.
Each of these options are set with `gmsh.option.setNumber(k,v)` just after
gmsh has been initialized.

# Level

Beginner
"""
function with_gmsh end

"""
    mesh_from_msh(msh_file;kwargs...)

Create a mesh object from a `.msh` file found in path `msh_file`.

See also [`mesh_from_gmsh`](@ref) and [`with_gmsh`](@ref).

# Keyword arguments

- `complexify=true` [optional]: If `complexify==true`, the mesh will be completed with all low dimensional faces into a face complex.
- `renumber=true` [optional]: If `renumber==true`, then `gmsh.model.mesh.renumberNodes()` and `gmsh.model.mesh.renumberElements()` will be called.
-  Any other keyword argument will be passed to function [`with_gmsh`](@ref).

# Level

Beginner
"""
function mesh_from_msh end

"""
    mesh_from_gmsh(gmsh::Module;complexify=true)

Create a mesh objects from the current state of the `gmsh` module.
If `complexify==true`, the mesh will be completed with all low dimensional faces into a face complex.

See also [`mesh_from_msh`](@ref) and [`with_gmsh`](@ref).

# Level

Beginner
"""
function mesh_from_gmsh end

function default_gmsh_options end

function reference_face_from_gmsh_eltype end


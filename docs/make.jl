using GalerkinToolkit
using Documenter
using Literate

codefence = "```julia" => "```"
src_dir = joinpath(@__DIR__,"src") 
pdes_dir = joinpath(src_dir,"pdes_automatic") 
pdes = ["poisson","p_laplacian","elasticity","stokes"]
for pde in pdes
    file_jl = joinpath(pdes_dir,pde*".jl")
    Literate.markdown(file_jl,pdes_dir)#;codefence)
end
assembly_dir = joinpath(src_dir,"pdes_manual") 
assembly = ["poisson","p_laplacian"]
for pde in assembly
    file_jl = joinpath(assembly_dir,pde*".jl")
    Literate.markdown(file_jl,assembly_dir)#;codefence)
end
tooling_dir = joinpath(src_dir,"tooling") 
tooling = ["mesh_generation"]
for pde in tooling
    file_jl = joinpath(tooling_dir,pde*".jl")
    Literate.markdown(file_jl,tooling_dir)#;codefence)
end

DocMeta.setdocmeta!(GalerkinToolkit, :DocTestSetup, :(using GalerkinToolkit); recursive=true)

makedocs(;
    modules=[GalerkinToolkit],
    authors="Francesc Verdugo <f.verdugo.rojano@vu.nl> and contributors",
    sitename="GalerkinToolkit",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://GalerkinToolkit.github.io/GalerkinToolkit.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Users guide" => "users_guide.md",
        "Developers guide" => "developers_guide.md",
        "Examples" => [
                       "Introduction" => "examples.md",
                       "PDEs (automatic assembly)"=> map(pde->"pdes_automatic/$(pde).md",pdes),
                       "PDEs (manual assembly)" => map(pde->"pdes_manual/$(pde).md",assembly),
                       "Tooling" => map(pde->"tooling/$(pde).md",tooling),
                      ],
        "API reference" =>[
                       "Introduction" => "reference.md",
                       "Mesh" => "reference/mesh.md",
                       "Integration" => "reference/integration.md",
                       "Interpolation" => "reference/interpolation.md",
                      ],
        "refindex.md",
    ],
)

deploydocs(;
    repo="github.com/GalerkinToolkit/GalerkinToolkit.jl",
    devbranch="main",
)

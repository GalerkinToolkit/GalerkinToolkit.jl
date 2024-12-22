using GalerkinToolkit
using Documenter
using Literate

codefence = "```julia" => "```"
src_dir = joinpath(@__DIR__,"src") 
pdes_dir = joinpath(src_dir,"pdes") 
pdes = ["poisson"]
for pde in pdes
    file_jl = joinpath(pdes_dir,pde*".jl")
    Literate.markdown(file_jl,pdes_dir;codefence)
end
assembly_dir = joinpath(src_dir,"assembly") 
assembly = ["poisson"]
for pde in assembly
    file_jl = joinpath(assembly_dir,pde*".jl")
    Literate.markdown(file_jl,assembly_dir)#;codefence)
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
                       "PDEs" => map(pde->"pdes/$(pde).md",pdes),
                       "Manual assembly" => map(pde->"assembly/$(pde).md",assembly),
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

using GalerkinToolkit
using Documenter
using Literate

src_dir = joinpath(@__DIR__,"src") 
examples_dir = joinpath(src_dir,"pdes") 
#examples = [
#            "problem_types",
#            "methods",
#            "mesh_generation",
#            "interpolations",
#            "boundary_conditions",
#            "fields",
#            "solvers",
#            "posprocessing",
#            "visualization",
#           ]
examples = ["poisson"]
codefence = "```julia" => "```"
for example in examples
    file_jl = joinpath(examples_dir,example*".jl")
    Literate.markdown(file_jl,examples_dir;codefence)
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
                       "PDEs" => map(example->"pdes/$(example).md",examples),
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

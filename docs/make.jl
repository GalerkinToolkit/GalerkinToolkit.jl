using GalerkinToolkit
using Documenter
using Literate

src_dir = joinpath(@__DIR__,"src") 
examples_dir = joinpath(src_dir,"examples") 
examples = [
            "problem_types",
            "mesh_generation",
            "interpolations",
            "boundary_conditions",
            "fields",
            "solvers",
            "posprocessing",
            "visualization",
           ]
codefence = "```julia" => "```"
for example in examples
    file_jl = joinpath(examples_dir,example*".jl")
    Literate.markdown(file_jl,examples_dir)#;codefence)
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
        "Examples" => [ "Introduction" => "examples.md", map(example->"examples/$(example).md",examples)...],
        "Users manual" =>[
                       "Introduction" => "reference.md",
                       "Mesh" => "reference/mesh.md",
                       "Integration" => "reference/integration.md",
                       "Interpolation" => "reference/interpolation.md",
                      ],
        "Developers guide" => "developers_guide.md",
        "refindex.md",
    ],
)

deploydocs(;
    repo="github.com/GalerkinToolkit/GalerkinToolkit.jl",
    devbranch="main",
)

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
for example in examples
    file_jl = joinpath(examples_dir,example*".jl")
    Literate.markdown(file_jl,examples_dir)
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
        "Examples" => map(example->"examples/$(example).md",examples),
        "Manual" =>[
                       "Getting Started" => "manual/getting_started.md",
                      ],
        "Reference" =>[
                       "Mesh" => "reference/mesh.md",
                       "Integration" => "reference/integration.md",
                       "Interpolation" => "reference/interpolation.md",
                      ]
    ],
)

deploydocs(;
    repo="github.com/GalerkinToolkit/GalerkinToolkit.jl",
    devbranch="main",
)

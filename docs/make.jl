using GalerkinToolkit
using Documenter
using Literate

src_dir = joinpath(@__DIR__,"src") 
examples_jl  = joinpath(src_dir,"examples.jl")
Literate.markdown(examples_jl,src_dir)

DocMeta.setdocmeta!(GalerkinToolkit, :DocTestSetup, :(using GalerkinToolkit); recursive=true)

makedocs(;
    modules=[GalerkinToolkit],
    authors="Francesc Verdugo <f.verdugo.rojano@vu.nl> and contributors",
    sitename="GalerkinToolkit.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://GalerkinToolkit.github.io/GalerkinToolkit.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "examples.md",
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

using GalerkinToolkit
using Literate
using Documenter

DocMeta.setdocmeta!(GalerkinToolkit, :DocTestSetup, :(using GalerkinToolkit); recursive=true)

src_dir = joinpath(@__DIR__, "src")
tutorials_dir = joinpath(@__DIR__, "src/tutorials")

tutorials = ["poisson_tutorial.jl"]
for tutorial in tutorials
    tutorial_jl = joinpath(@__DIR__, tutorial)
    Literate.markdown(tutorial_jl, tutorials_dir)
end

makedocs(;
    modules=[GalerkinToolkit],
    authors="Francesc Verdugo <f.verdugo.rojano@vu.nl> and contributors",
    sitename="GalerkinToolkit.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://fverdugo.github.io/GalerkinToolkit.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Tutorials" =>[
            "Poisson equation" => "tutorials/poisson_tutorial.md",
        ],
        "Reference" =>[
                       "Geometry" => "reference/geometry.md",
                       "Integration" => "reference/integration.md",
                       "Interpolation" => "reference/interpolation.md",
                      ]
    ],
)

deploydocs(;
    repo="github.com/fverdugo/GalerkinToolkit.jl",
    devbranch="main",
)

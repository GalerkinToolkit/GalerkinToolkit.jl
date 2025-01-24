using GalerkinToolkit
using Documenter
using Literate

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
        "index.md",
        "tutorials.md",
        "examples.md",
        "manual.md",
        "reference.md",
        "refindex.md",
    ],
)

deploydocs(;
    repo="github.com/GalerkinToolkit/GalerkinToolkit.jl",
    devbranch="main",
)

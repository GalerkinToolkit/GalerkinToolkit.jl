using GalerkinToolkit
using Documenter

DocMeta.setdocmeta!(GalerkinToolkit, :DocTestSetup, :(using GalerkinToolkit); recursive=true)

makedocs(;
    modules=[GalerkinToolkit],
    authors="Francesc Verdugo <fverdugo@cimne.upc.edu> and contributors",
    sitename="GalerkinToolkit.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://fverdugo.github.io/GalerkinToolkit.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Reference" =>[
                       "Mesh interface" => "reference/mesh_interface.md"
                      ]
    ],
)

deploydocs(;
    repo="github.com/fverdugo/GalerkinToolkit.jl",
    devbranch="main",
)

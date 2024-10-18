using GalerkinToolkit
using Documenter

DocMeta.setdocmeta!(
    GalerkinToolkit,
    :DocTestSetup,
    :(using GalerkinToolkit);
    recursive = true,
)

makedocs(;
    modules = [GalerkinToolkit],
    authors = "Francesc Verdugo <fverdugo@cimne.upc.edu> and contributors",
    sitename = "GalerkinToolkit.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://GalerkinToolkit.github.io/GalerkinToolkit.jl",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Reference" => [
            "Mesh" => "reference/mesh.md",
            "Integration" => "reference/integration.md",
            "Interpolation" => "reference/interpolation.md",
        ],
    ],
)

deploydocs(; repo = "github.com/GalerkinToolkit/GalerkinToolkit.jl", devbranch = "main")

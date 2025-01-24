using GalerkinToolkit
using Documenter
using Literate

codefence = "```julia" => "```"
src_jl = joinpath(@__DIR__,"src","src_jl")
src_md = joinpath(@__DIR__,"src","src_md")
rm(src_md,force=true,recursive=true)
mkpath(src_md)

for file_jl in filter(f->f[end-2:end]==".jl",readdir(src_jl))
    Literate.markdown(joinpath(src_jl,file_jl),src_md)#;codefence)
end


DocMeta.setdocmeta!(GalerkinToolkit, :DocTestSetup, :(using GalerkinToolkit); recursive=true)

examples_pages = [
    "example_hello_world.md",
]
examples_pages = map(p->joinpath("src_md",p),examples_pages)
examples = ["examples.md",examples_pages...]

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
        "Examples" => examples,
        "manual.md",
        "reference.md",
        "refindex.md",
    ],
)

deploydocs(;
    repo="github.com/GalerkinToolkit/GalerkinToolkit.jl",
    devbranch="main",
)

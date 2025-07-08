module Make

using GalerkinToolkit
using Documenter
using Literate
using DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "src", "references.bib"))

const src_jl = joinpath(@__DIR__,"src","src_jl")
const src_md = joinpath(@__DIR__,"src","src_md")

function media()
    #rm(src_md,force=true,recursive=true)
    mkpath(src_md)

    for file_jl in filter(f->f[end-2:end]==".jl",readdir(src_jl))
        f = joinpath(src_jl,file_jl)
        println("Running file $f")
        m = Module()
        Base.include(m,f)
    end
end

function main(;debug=false)

    if debug
        for file_png in filter(f->f[end-3:end]==".png" || f[end-3:end]==".gif" ,readdir(src_jl))
            cp(joinpath(src_jl,file_png),joinpath(src_md,file_png),force=true)
        end
    else
        rm(src_md,force=true,recursive=true)
    end
    mkpath(src_md)

    codefence = "```julia" => "```"
    for file_jl in filter(f->f[end-2:end]==".jl",readdir(src_jl))
        f = joinpath(src_jl,file_jl)
        if debug
            Literate.markdown(f,src_md;codefence)
        else
            Literate.markdown(f,src_md)
        end
    end

    DocMeta.setdocmeta!(GalerkinToolkit, :DocTestSetup, :(using GalerkinToolkit); recursive=true)

    tutorials_pages = [
                       "tutorial_intro_to_fem.md",
                      ]
    tutorials_pages = map(p->joinpath("src_md",p),tutorials_pages)
    tutorials = ["tutorials.md",tutorials_pages...]

    examples_pages = [
                      "example_hello_world.md",
                      "example_poisson_equation.md",
                      "example_p_laplacian.md",
                      "example_stokes.md",
                      "example_transient_heat_eq.md",
                     ]
    examples_pages = map(p->joinpath("src_md",p),examples_pages)
    examples = ["examples.md",examples_pages...]

    manual_pages = [
                    joinpath("manual","introduction.md"),
                    joinpath("src_md","manual_geometry.md"),
                    joinpath("manual","for_developers.md")
                   ]
    manual = ["manual.md",manual_pages...]


    makedocs(;
             modules=[GalerkinToolkit],
             authors="Francesc Verdugo <f.verdugo.rojano@vu.nl> and contributors",
             sitename="GalerkinToolkit.jl",
             plugins=[bib],
             format=Documenter.HTML(;
                                    prettyurls=get(ENV, "CI", "false") == "true",
                                    canonical="https://GalerkinToolkit.github.io/GalerkinToolkit.jl",
                                    assets=["assets/citations.css","assets/favicon.ico"],
                                   ),
             pages=[
                    "index.md",
                    "Manual" => manual,
                    "Examples" => examples,
                    "Tutorials" => tutorials,
                    "reference.md",
                    "refindex.md",
                    "references.md",
                   ],
            )

    deploydocs(;
               repo="github.com/GalerkinToolkit/GalerkinToolkit.jl",
               devbranch="main",
              )

end

end # module

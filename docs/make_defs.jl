module Make

using GalerkinToolkit
using Documenter
using Literate

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
                      "example_hello_world_manual.md",
                      "example_poisson_equation.md",
                      "example_p_laplacian.md",
                      "example_p_laplacian_manual.md",
                      "example_stokes.md",
                     ]
    examples_pages = map(p->joinpath("src_md",p),examples_pages)
    examples = ["examples.md",examples_pages...]

    manual_pages = [
                    "getting_started.md",
                    "for_developers.md"
                   ]
    manual_pages = map(p->joinpath("manual",p),manual_pages)
    manual = ["manual.md",manual_pages...]


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
                    "index.md",
                    "Tutorials" => tutorials,
                    "Examples" => examples,
                    "Manual" => manual,
                    "reference.md",
                    "refindex.md",
                   ],
            )

    deploydocs(;
               repo="github.com/GalerkinToolkit/GalerkinToolkit.jl",
               devbranch="main",
              )

end

end # module

using Pkg

Pkg.activate(@__DIR__)
Pkg.develop(path=normpath(joinpath(@__DIR__,"..","..")))
Pkg.instantiate(;allow_autoprecomp=false)

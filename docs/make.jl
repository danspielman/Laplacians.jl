
using Documenter
using Laplacians


makedocs(modules=[Laplacians], doctest = false)

#=
deploydocs(
    deps = Deps.pip("pygments", "mkdocs", "mkdocs-material", "python-markdown-math"),
    repo   = "github.com/danspielman/Laplacians.jl.git"
#    julia  = "release"
)
=#

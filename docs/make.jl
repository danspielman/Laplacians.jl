
DOCUMENTER_DEBUG = true

using Documenter
using Laplacians


makedocs(modules=[Laplacians], doctest = false)

deploydocs(
           repo = "github.com/danspielman/Laplacians.jl.git",
           deps   = Deps.pip("mkdocs", "python-markdown-math"),
           julia  = "0.6.0"
)

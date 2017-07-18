
using Documenter
using Laplacians


makedocs(modules=[Laplacians], doctest = false, debug=true)

deploydocs(
           repo = "github.com/danspielman/Laplacians.jl.git",
           deps   = Deps.pip("mkdocs", "python-markdown-math"),
           julia  = "0.5.1"
)

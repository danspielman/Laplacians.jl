# copied this from LightGraphs

using Documenter
include("../src/Laplacians.jl")
using Laplacians

# index is equal to the README for the time being
cp(normpath(@__FILE__, "../../README.md"), normpath(@__FILE__, "../src/index.md"); remove_destination=true)

# same for contributing and license
#cp(normpath(@__FILE__, "../../CONTRIBUTING.md"), normpath(@__FILE__, "../src/contributing.md"); remove_destination=true)
cp(normpath(@__FILE__, "../../LICENSE.md"), normpath(@__FILE__, "../src/license.md"); remove_destination=true)

makedocs(modules=[Laplacians], doctest = false)

#=
deploydocs(
    deps = Deps.pip("pygments", "mkdocs", "mkdocs-material", "python-markdown-math"),
    repo   = "github.com/danspielman/Laplacians.jl.git"
#    julia  = "release"
)
=#

rm(normpath(@__FILE__, "../src/index.md"))
#rm(normpath(@__FILE__, "../src/contributing.md"))
rm(normpath(@__FILE__, "../src/license.md"))


DOCUMENTER_DEBUG = true

using Documenter
using Laplacians


makedocs(modules=[Laplacians],
    doctest = false,
    sitename = "Laplacians.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"),
    pages = [
        "About" => "index.md",
        "Manual" => Any[
            "Installation" => "Installation.md",
            "Examples" => "Examples.md",
            "Sparse matrices as graphs" => "CSCgraph.md",
            "Solving Linear Equations" => "usingSolvers.md",
            "Low Stretch Spanning Trees" => "LSST.md"
        ],
        "Developing" => "Developing.md",
        "API" => Any[
            "generators" => "graphGenerators.md",
            "operators"  => "operators.md",
            "graphUtils" => "graphUtils.md",
            "graphAlgs" => "graphAlgs.md",
            "IO" => "IO.md",
            "solvers" => "solvers.md",
            "sparsification" => "sparsification.md",
            "akpw" => "akpw.md",
            "treeAlgs" => "treeAlgs.md",
            "randTrees" => "randTrees.md",
            "localClustering" => "localClustering.md",
            "Private Functions" => "privateFuncs.md",
            "All of the above" => "indexOfAll.md"
        ]

    ]

    )

deploydocs(
           repo = "github.com/danspielman/Laplacians.jl.git",
#           deploy_config = Documenter.GitHubActions(),
)

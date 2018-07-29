using Test
using Laplacians
using LinearAlgebra
using Random
using SparseArrays
using LinearAlgebra
using RandomV06

L = Laplacians
include("testIncludes.jl")
include("test_graph_generators.jl")

include("testByFile.jl")

include("testByExport.jl")

include("testSolvers.jl")

include("testPCG.jl")




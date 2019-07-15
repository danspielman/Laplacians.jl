"""A package for graph computations related to graph Laplacians

Graphs are represented by sparse adjacency matrices, etc.
"""
module Laplacians


  function __init__()

    #=
    if !isdefined(Main, :LAPLACIANS_NOPLOT) & !haskey(ENV,"LAPLACIANS_NOPLOT")
        eval(Expr(:using, :PyPlot))
    end
    =#

    if isdefined(Main, :LAPLACIANS_AMG)
        eval(Expr(:using, :PyAMG))

    end
  end

  using RandomV06
  V06 = RandomV06.V06
  V07 = RandomV06.V07
  Vcur = RandomV06.Vcur

  using Arpack
  using SuiteSparse
  using DataStructures
  using SparseArrays
  using Random
  using LinearAlgebra
  using Statistics
  using Printf

  using DelimitedFiles

  using Plots

  include("IJV.jl")

  include("fastCSC.jl")
  export symPermuteCSC
  export symTransposeCSC
  export submatrixCSC

  include("graphUtils.jl")
  export deg
  export nbri
  export weighti
  export nbrs
  export wdeg
  export setValue
  export backIndices
  export flipIndex
  export findEntries
  export compConductance
  export getVolume
  export getObound

  include("graphGenerators.jl")
  export ring_graph
  export generalized_ring
  export rand_matching
  export rand_regular
  export grown_graph
  export grown_graph_d
  export pref_attach
  export hypercube
  export complete_binary_tree
  export complete_graph
  export empty_graph
  export path_graph

  #export wgrid2

  export grid2
  export grid2coords
  export grid3

  export rand_gen_ring
  export randperm
  export ErdosRenyi
  export ErdosRenyiCluster
  export ErdosRenyiClusterFix
  export pure_random_graph

  export chimera
  export rand_weight
  export wted_chimera, semiwted_chimera

  include("latinSquares.jl")
  export latin_square, latin_square_graph

  include("IO.jl")
  export writeIJV, read_graph

  include("graphOps.jl")

  export unweight, unweight!
  export mapweight
  export uniformWeight, uniformWeight!

  export edgeVertexMat, wtedEdgeVertexMat

  export power, thicken_once, thicken

  export product_graph
  export generalized_necklace
  export subsampleEdges

  export two_lift
  export join_graphs, join_graphs!, disjoin

  export plot_graph

  export shortIntGraph, floatGraph

  export lap
  export adj
  export spectral_coords
  export spectral_drawing

  export diagmat

  include("graphAlgs.jl")

  export components
  export biggestComp
  export vecToComps
  export isConnected

  export shortestPaths, shortestPathTree, pathFromParents
  export kruskal, prim

  include("treeAlgs.jl")

  export RootedTree
  export matToTree
  export matToTreeDepth
  export tarjanStretch
  export compDepth
  export comp_stretches
  export dfsOrder

  include("pcg.jl")

  export cg, cgSolver, cgLapSolver
  export pcg, pcgSolver, pcgLapSolver

  include("flow.jl")

  export maxflow


  include("akpw.jl")

  export akpw, akpwU


  include("localClustering.jl")

  export prn
  export apr
  export localImprove

  include("cutHeuristics.jl")

  export refineCut
  export dumbRefineCut

  include("randTrees.jl")
  export randishKruskal, randishPrim

  include("sampler.jl")
  include("fastSampler.jl")
  export FastSampler, sample, sampleMany
  export blockSample

  include("solverInterface.jl")
  export chol_sddm, chol_lap, lapWrapSDDM


  include("augTreeSolver.jl")

  export augmentTree, augTreePrecon, augTreeSddm
  export augTreeLapPrecon, augTreeLap, AugTreeParams, AugTreeParamsOld

  include("KMPSolver.jl")
  export KMPSDDMSolver
  export KMPLapSolver
  export KMPParams

  include("approxCholTypes.jl")
  include("approxChol.jl")
  export approxchol_lap, ApproxCholParams, approxchol_sddm

  include("fiedler.jl")
  export fiedler

  include("complexSolvers.jl")
  export SDDMSolvers
  export LapSolvers

  include("compare_solvers.jl")
  export SolverTest, speedTestLapSolvers

  include("conditionNumber.jl")
  export support, approxQual, conditionNumber

  include("sparsify.jl")
  export sparsify

  include("johnlind.jl")
  export johnlind

  include("toposort.jl")
  export toposort, dirEdgeVertexMat

  # include("isotonicIPM.jl")
  # export isotonicIPM, isotonicIPMrelEps

  include("harmonic.jl")
  export harmonic_interp

  include("from_cholmod.jl")
  export cholmod_perm, ask_cholmod

  include("deprecated.jl")

  include("graphGenGeom.jl")
  export plot_graph_weighted
  export ggrid2_ijv, ggrid2
  export ggrid2_checkered_ijv, ggrid2_checkered
  export ggrid2coords
  export getInterior2, getBoundary2
  export ggrid3_ijv, ggrid3
  export ggrid3_checkered_ijv, ggrid3_checkered
  export ggrid3coords
  export getInterior3, getBoundary3
  
end # module Laplacians.jl

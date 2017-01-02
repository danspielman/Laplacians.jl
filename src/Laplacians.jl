"""A package for graph computations related to graph Laplacians

Graphs are represented by sparse adjacency matrices, etc.
"""
module Laplacians

#=

Started by Dan Spielman.
Other contributors:


This represents graphs as sparse adjacency matrices.
The reasons are twofold:
  # 1. It is much faster than the Graphs.jl library, and
  # 2. It is more compatible with linear algebraic operations

This module begins by importing the packages it requires.
It then defines a few functions for dealing with the graphs:
deg(graph, node),  nbri(graph, node, i)  and weighti(graph, node i).
Using these is slower than what we actually do in the code.
But, you can start with this and then convert.

This then includes code from other files,
and exports the functions for which it seems appropriate.

=#


  function __init__()
    if !isdefined(Main, :LAPLACIANS_NOPLOT)
        eval(Expr(:using, :PyPlot))
    end

    if !isdefined(Main, :LAPLACIANS_NOAMG)
        eval(Expr(:using, :PyAMG))
    end
  end

  using DataStructures

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
  export readIJ
  export ringGraph
  export generalizedRing
  export randMatching
  export randRegular
  export grownGraph
  export grownGraphD
  export prefAttach
  export hyperCube
  export completeBinaryTree
  export completeGraph
  export pathGraph

  export wGrid2
  export wGrid3

  export grid2
  export grid2coords
  export grid3

  export randGenRing
  export randperm
  export ErdosRenyi
  export ErdosRenyiCluster
  export ErdosRenyiClusterFix
  export pureRandomGraph

  export chimera
  export randWeight
  export wtedChimera, semiWtedChimera

  include("IO.jl")
  export readIJ, readIJV, writeIJV

  include("graphOps.jl")

  export unweight, unweight!
  export mapweight
  export uniformWeight, uniformWeight!

  export edgeVertexMat

  export productGraph
  export generalizedNecklace
  export subsampleEdges

  export twoLift
  export joinGraphs, disjoin

  export plotGraph

  export shortIntGraph, floatGraph

  export lap
  export adj
  export spectralCoords
  export spectralDrawing

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
  export compStretches
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

  include("solverInterface.jl")
  export cholSDDM, cholLap, lapWrapSDDM


  include("augTreeSolver.jl")

  export augmentTree, augTreePrecon, augTreeSddm
  export augTreeLapPrecon, augTreeLap, AugTreeParams, AugTreeParamsOld

  include("externalSolvers.jl")
  export AMGSolver, AMGLapSolver

  include("KMPSolver.jl")
  export KMPSDDMSolver
  export KMPLapSolver
  export KMPParams

  include("complexSolvers.jl")
  export SDDMSolvers
  export LapSolvers

  include("johnlind.jl")
  export johnlind

  include("toposort.jl")
  export toposort, dirEdgeVertexMat

  # include("isotonicIPM.jl")
  # export isotonicIPM, isotonicIPMrelEps


end # module Laplacians.jl

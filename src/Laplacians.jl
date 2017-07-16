"""A package for graph computations related to graph Laplacians

Graphs are represented by sparse adjacency matrices, etc.
"""
module Laplacians


  function __init__()
    if !isdefined(Main, :LAPLACIANS_NOPLOT)
        eval(Expr(:using, :PyPlot))
    end

    if isdefined(Main, :LAPLACIANS_AMG)
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

  export edgeVertexMat, wtedEdgeVertexMat

  export power, thicken_once, thicken

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
  export blockSample

  include("solverInterface.jl")
  export cholSDDM, cholLap, lapWrapSDDM


  include("augTreeSolver.jl")

  export augmentTree, augTreePrecon, augTreeSddm
  export augTreeLapPrecon, augTreeLap, AugTreeParams, AugTreeParamsOld

  include("KMPSolver.jl")
  export KMPSDDMSolver
  export KMPLapSolver
  export KMPParams

  include("approxCholTypes.jl")
  include("approxChol.jl")
  export approxCholLap, ApproxCholParams

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


  include("from_cholmod.jl")
  export cholmod_perm, ask_cholmod

end # module Laplacians.jl

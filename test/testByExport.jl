# Test every function exported in Laplacians.jl by running it at least once
# Create this by copying over the export commands from Laplacians.jl

n = 101
a = wtedChimera(n,1)

# export symPermuteCSC

rp = randperm(n)
ap = symPermuteCSC(a,rp);
rpi = zeros(Int,n);
rpi[rp] = 1:n

@test sum(abs(a-symPermuteCSC(ap,rpi))) == 0

# export symTransposeCSC

a2 = triu(a) + (tril(a) .> 0)
@test sum(abs(symTransposeCSC(a2) - a2')) == 0

# export submatrixCSC

s = collect(1:10)
submatrixCSC(a,s)

  # export deg
  # export nbri
  # export weighti
  # export nbrs
  # export wdeg

x = 0
y = 0
z = 0
w = 0
for i in 1:n
    for j in 1:deg(a,i)
        x += weighti(a,i,j)
        y += a[i,nbri(a,i,j)]
    end
    for j in nbrs(a,i)
        z += a[i,j]
    end
    w += wdeg(a,i)
end
@test x == y == z == w == sum(a) 

# export setValue

setValue(a,1,1,0.0)

# export backIndices

b = backIndices(a)

# export flipIndex

b = flipIndex(a)

# export findEntries

u,v,w = findEntries(a)

# export compConductance

compConductance(a,collect(1:10))

# export getVolume

getVolume(a,collect(1:10))

# export getObound

getObound(a,collect(1:10))

  # export readIJ
  # export ringGraph
  # export generalizedRing
  # export randMatching
  # export randRegular
  # export grownGraph
  # export grownGraphD
  # export prefAttach
  # export hyperCube
  # export completeBinaryTree
  # export completeGraph
  # export pathGraph

  # export wGrid2
  # export wGrid3

  # export grid2
  # export grid2coords
  # export grid3

  # export randGenRing
  # export randperm
  # export ErdosRenyi
  # export ErdosRenyiCluster
  # export ErdosRenyiClusterFix
  # export pureRandomGraph

  # export chimera
  # export randWeight
  # export wtedChimera, semiWtedChimera

  # export readIJ, readIJV, writeIJV

  # export unweight, unweight!
  # export mapweight
  # export uniformWeight, uniformWeight!

  # export edgeVertexMat

  # export productGraph
  # export generalizedNecklace
  # export subsampleEdges

  # export twoLift
  # export joinGraphs, disjoin

  # export plotGraph

  # export shortIntGraph, floatGraph

  # export lap
  # export adj
  # export spectralCoords
  # export spectralDrawing

  # export toUnitVector

  # export diagmat

  # export components
  # export biggestComp
  # export vecToComps
  # export isConnected

  # export shortestPaths, shortestPathTree, pathFromParents
  # export kruskal, prim

  # export RootedTree
  # export matToTree
  # export matToTreeDepth
  # export tarjanStretch
  # export compDepth
  # export compStretches
  # export dfsOrder

  # export cg, cgSolver
  # export pcg, pcgSolver, pcgLapSolver

  # export maxflow

  # export akpw, akpwU

  # export prn
  # export apr
  # export localImprove

  # export refineCut
  # export dumb

  # export randishKruskal, randishPrim

  # export FastSampler, sample, sampleMany

  # export johnlind

  # export toposort, dirEdgeVertexMat


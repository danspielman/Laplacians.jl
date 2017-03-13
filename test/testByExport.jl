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

a2 = grownGraph(100,3)

# export grownGraphD

a2 = grownGraphD(100,3)

# export prefAttach

a2 = prefAttach(100,3,0.5)

# export hyperCube

a2 = hyperCube(3)

# export completeBinaryTree

a2 = completeBinaryTree(7)

# export completeGraph

a2 = completeGraph(7)

# export pathGraph

a2 = pathGraph(7)

# export wGrid2

a2 = wGrid2(3)

# export wGrid3

a3 = wGrid2(3)

# export grid2

a2 = grid2(3)

# export grid2coords

a2 = grid2coords(3)

# export grid3

a3 = grid3(3)

  # export randGenRing
  # export randperm

# export ErdosRenyi

a2 = ErdosRenyi(100,300)

# export ErdosRenyiCluster

a2 = ErdosRenyiCluster(100,4)

# export ErdosRenyiClusterFix

a2 = ErdosRenyiClusterFix(100,4)

# export pureRandomGraph

# export chimera
# export randWeight
# export wtedChimera, semiWtedChimera

for i in 1:5
    a2 = semiWtedChimera(10000,i)
end


# export readIJ, readIJV, writeIJV

n = 101
a = wtedChimera(n,1)
writeIJV("tmp.txt",a)
a2 = readIJV("tmp.txt")
@test sum(abs(a-a2)) == 0

# export unweight, unweight!

a2 = unweight(a2)
unweight!(a2)

  # export mapweight
  # export uniformWeight, uniformWeight!

a2 = uniformWeight(a2)
uniformWeight!(a2)

  # export edgeVertexMat

b = edgeVertexMat(a2)
@test sum(abs(b'*b - lap(unweight(a2)))) == 0

  # export productGraph
  # export generalizedNecklace
  # export subsampleEdges

subsampleEdges(a,0.5)

# export twoLift

twoLift(a)
twoLift(a,3)

  # export joinGraphs, disjoin

  # export plotGraph

  # export shortIntGraph, floatGraph

  # export lap
  # export adj
  # export spectralCoords
  # export spectralDrawing

spectralDrawing(a)

  # export toUnitVector

# export diagmat

diagmat(a)

  # export components
  # export biggestComp
  # export vecToComps
  # export isConnected

# export shortestPaths, shortestPathTree, pathFromParents

a = wtedChimera(102,1)
shortestPaths(a,1)
shortestPathTree(a,1)

Laplacians.intHeapSort(randn(10))

nh = Laplacians.intHeap(10)
for i in 1:10
    Laplacians.intHeapAdd!(nh, i, rand())
end
Laplacians.intHeapSort(nh)

  # export kruskal, prim


# export RootedTree
# export matToTree

# export matToTreeDepth

a = wtedChimera(101,1)
t = akpw(a)
tr = matToTree(t)
tr, d1 = matToTreeDepth(t);
d2 = Laplacians.treeDepthDFS(t)

  # export tarjanStretch
  # export compDepth
  # export compStretches
  # export dfsOrder

Laplacians.bfsOrder(t,1)

  # export cg, cgSolver
  # export pcg, pcgSolver, pcgLapSolver

  # export maxflow

for i in 1:10
    maxflow(wtedChimera(102,i),1,2)
end

  # export akpw, akpwU

akpw(wtedChimera(10000,1),ver=2)

  # export prn
  # export apr
  # export localImprove
  # export refineCut
  # export dumb

a = chimera(100, 3);
s = prn(a, [1,2,3], 0.2, 5);
conds = compConductance(a, s)
#println(conds, " ", length(s))
minEpsSigma = getVolume(a, s) / getVolume(a, setdiff(collect(1:max(a.n, a.m)), s));
cut, flow = localImprove(a, s, epsSigma = minEpsSigma);
condcut = compConductance(a, cut)
heur = refineCut(a, cut)
dumbRefineCut(a,collect(1:10))

#println(condcut, " ", length(cut))

  # export randishKruskal, randishPrim

  # export FastSampler, sample, sampleMany

  # export johnlind

johnlind(a)

  # export toposort, dirEdgeVertexMat


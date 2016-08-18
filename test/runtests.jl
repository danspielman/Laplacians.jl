using Base.Test
using Laplacians

# generate chimeric graphs, and test run many routines on them

function isTree(gr::SparseMatrixCSC)
    isConnected(gr) && (nnz(gr) == 2*(gr.n-1))
end



function testSolvers(a)

    n = a.n
    excess = zeros(n); excess[1] = excess[n] = 1e-5;
    la = lap(a);
    sdd = la + spdiagm(excess);
    b = rand(n); b = b - mean(b);

    for solver in SDDSolvers
        println(solver)
        @time f = solver(sdd, tol=1e-6, maxtime=0.001);
        @time x = f(b);
        println(norm(sdd * x - b) / norm(b))
        println()
    end

    for solver in LapSolvers
        println(solver)
        @time f = solver(a, tol=1e-6, maxtime=0.001);
        @time x = f(b);
        println(norm(la * x - b) / norm(b))
        println()
    end

end



n = 1000
tol = 1e-11

for i in 1:100
    println("wtedChimera($n, $i)")
    gr = wtedChimera(n,i)
    @test gr.n == n
    @test isConnected(gr)
    @test isTree(kruskal(gr))
    @test isTree(randishKruskal(gr))
    @test isTree(randishPrim(gr))
    @test abs(sum(kruskal(gr)) - sum(prim(gr))) < tol
    @test sum(kruskal(gr,kind=:min)) <= sum(kruskal(gr,kind=:max))

    testSolvers(gr)
end




# every single function/method should have at least one test!


# assume each vertex pair has at most ONE unique edge
function isSubsetOfEdges(subset::SparseMatrixCSC, superset::SparseMatrixCSC, marginOfError)
  nSubsetVertices = subset.n

  for vInd in 1:nSubsetVertices
    for eSubInd in subset.colptr[vInd]:(subset.colptr[vInd+1]-1)
      for eSupInd in superset.colptr[vInd]:(superset.colptr[vInd+1]-1)
        if superset.rowval[eSupInd] == subset.rowval[eSubInd]
          if abs(superset.nzval[eSupInd] - subset.nzval[eSubInd]) > marginOfError
            println("edge: v1: ", vInd, " v2: ", subset.rowval[eSubInd], " edgeWeight: ", subset.nzval[eSubInd])
            println("\t=/= edge: v1: ", vInd, " v2: ", superset.rowval[eSupInd], " edgeWeight: ", superset.nzval[eSupInd])
            return false
          end

          break
        end
      end
    end
  end #for

  return true

end #isSubsetOfEdges


# Test Function

# uses four different types of graphs: 
# dim will increase by "iterationsize" each time for each number of iterations
"""
testAKPW runs tests on akpw to make sure it runs correctly. It takes a few parameters:

startingDim is the dimension of the first set of graphs it will test (the number of nodes is dim squared)
iterationsize is the amount that dim will grow each iteration
numIterations is the number of iterations
kindVar, randomClustersVar, metisClusteringVar, exponentialX, and shuffleClustersVar each pertain to the equivalent options for akpw
marginOfError is the allowable margin of error for running edgeSubset tests (if the difference in edge weights between the tree and graph is above this margin, it will not pass the test)

In each iteration, 4 graphs (one grid, one product, one necklace and one chimera) each with
(startingDim+((numIterations-1)*iterationsize))^2 nodes are produced, and akpw is run on each. The resulting tree
is tested using two functions: isTree and isSubsetOfEdges. isTree tests to see if it is, in fact, a tree.
isSubsetOfEdges makes sure all of the edges in the tree are also edges in the original graph (this works with 
  akpw! as well, if you want to test that function instead, by using the marginOfError option).

for example:

  testAKPW(4, 2, 3) should 3 sets of four (grid, product, necklace and chimera) graphs: with 16, 36, and 64 nodes.
  None should produce any errors
"""
function testAKPW(startingDim, iterationSize, numIterations; kindVar=:max, randomClustersVar=false, metisClusteringVar=false, shuffleClustersVar=true, exponentialX=true, marginOfError=.0000001)

  if startingDim < 4
    println("please choose a higher startingDim (min: 4)")
    return
  end

  dim = startingDim

  for i in 1:numIterations

    a = grid2(dim)
    n = size(a)[1]
    (ai,aj,av) = findnz(triu(a))
    gridGraph = sparse(ai,aj,rand(size(av)),n,n)
    gridGraph = gridGraph + gridGraph';


    a = productGraph(generalizedRing(dim,[1 5]), grownGraphD(dim,3));
    n = size(a)[1]
    (ai,aj,av) = findnz(triu(a))
    productGraph1 = sparse(ai,aj,rand(size(av)),n,n)
    productGraph1 = productGraph1 + productGraph1'

    a = generalizedNecklace(generalizedRing(dim,[1 5]), grownGraphD(dim,3),5);
    n = size(a)[1]
    (ai,aj,av) = findnz(triu(a))
    necklaceGraph = sparse(ai,aj,rand(size(av)),n,n)
    necklaceGraph = necklaceGraph + necklaceGraph'

    chimGraph = wtedChimera(dim^2)

    rows, columns, edgeWeights = findnz(gridGraph)
    akpwTreeGrid = akpw(gridGraph, kind=kindVar, metisClustering = metisClusteringVar, randomClusters = randomClustersVar, shuffleClusters = shuffleClustersVar, exponentialX = exponentialX)
    oldGrid = sparse(rows, columns, edgeWeights)

    rows, columns, edgeWeights = findnz(productGraph1)
    akpwTreeProduct = akpw(productGraph1, kind=kindVar, metisClustering = metisClusteringVar, randomClusters = randomClustersVar, shuffleClusters = shuffleClustersVar, exponentialX = exponentialX)
    oldProductGraph = sparse(rows, columns, edgeWeights)

    rows, columns, edgeWeights = findnz(necklaceGraph)
    akpwTreeNecklace = akpw(necklaceGraph, kind=kindVar, metisClustering = metisClusteringVar, randomClusters = randomClustersVar, shuffleClusters = shuffleClustersVar, exponentialX = exponentialX)
    oldNecklace = sparse(rows, columns, edgeWeights)

    rows, columns, edgeWeights = findnz(chimGraph)
    akpwTreeChim = akpw(chimGraph, kind=kindVar, metisClustering = metisClusteringVar, randomClusters = randomClustersVar, shuffleClusters = shuffleClustersVar, exponentialX = exponentialX)
    oldChimera = sparse(rows, columns, edgeWeights)

    if !isTree(akpwTreeGrid)
      write(STDERR, "Grid is not a tree\n")
    end

    if !isSubsetOfEdges(akpwTreeGrid, oldGrid, marginOfError)
      write(STDERR, "akpwTreeGrid is not a subset of the edges in gridGraph\n")
    end


    if !isTree(akpwTreeProduct)
      write(STDERR, "productGraph is not a tree\n")
    end

    if !isSubsetOfEdges(akpwTreeProduct, oldProductGraph, marginOfError)
      write(STDERR, "akpwTreeProduct is not a subset of the edges in productGraph\n")
    end


    if !isTree(akpwTreeNecklace)
      write(STDERR, "necklace is not a tree\n")
    end

    if !isSubsetOfEdges(akpwTreeNecklace, oldNecklace, marginOfError)
      write(STDERR, "akpwTreeNecklace is not a subset of the edges in necklaceGraph\n")
    end


    if !isTree(akpwTreeChim)
      write(STDERR, "chimera is not a tree\n")
    end

    if !isSubsetOfEdges(akpwTreeChim, oldChimera, marginOfError)
      write(STDERR, "akpwTreeChim is not a subset of the edges in chimGraph\n")
    end


    dim += iterationSize

  end#for

  println("Tests PASS!")
  
end #testAKPW

include("testPCG.jl")

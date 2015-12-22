#=
Written by Jackson Thea for CPSC 490 Fall 2015

Provides
  akpw! and akpw (a wrapper for akpw!)

AKPW is an algorithm published by Alon, Karp, Peleg and West in 1995 to produce low stretch spanning trees. This
is my implementation of the algorithm. For more information on how to run AKPW, see the doc string above akpw!'s
function declaration.
=#



function normalizeEdgeWeights(mat::SparseMatrixCSC, kind, nEdges)

  if (kind == :max)
    for eInd in 1:nEdges
      mat.nzval[eInd] \= 1.0
    end #for
  end #if

  minEdgeWeight = minimum(mat.nzval)

  for eInd in 1:nEdges
    mat.nzval[eInd] /= minEdgeWeight
  end #for

  return minEdgeWeight
end #normalizeEdgeWeights


function divideEdgesIntoClasses(mat::SparseMatrixCSC, nEdges, rows, columns, edgeWeights)
  nVertices = mat.n
  x = exp(sqrt(log(nVertices) * log(log(nVertices))))
  y = x * log(nVertices) #placeholder for testing
  edgeClasses = zeros(Int64, nEdges)

  for eInd in 1:nEdges
    edgeClasses[eInd] = floor(log(y, edgeWeights[eInd])) + 1
  end #for

  return edgeClasses

end #divideEdgesIntoClasses



function unsortedArray(n)
  permutation = zeros(Int64, n)
  for i in 1:n
    permutation[i] = i
  end

  for i in 1:n
    j = rand(1:i)
    permutation[i] = permutation[j]
    permutation[j] = i
  end #for

  return permutation
end #unsortedArray


#could speed this up by creating linked list for clusterCount (cuz it'll usually be mostly zeros)
#this makes partitionList potentially unsorted. Will that mess shit up..?
# should be greatest total weight! not just the greatest number of neighbors
function reshuffleClusters(
  mat,
  partitionList,
  vertexToCluster, 
  starts, 
  vertexToClusterLocation, 
  finalRoundClusterVertices, 
  rows, 
  columns, 
  edgeWeights
)
  nVertices = mat.n
  visited = zeros(Bool, nVertices)
  nClusters = length(partitionList)

  for vInd in 1:nVertices
    validVertex = true
    if starts[vertexToCluster[vInd]] == vInd || !finalRoundClusterVertices[vInd]
      continue
    end #if

    clusterNeighborCount = Dict{Int64, Int64}()

    newCluster = vertexToCluster[vInd]
    newClusterCount = 0

    for eInd in mat.colptr[vInd]:(mat.colptr[vInd+1]-1) #eInd is the edge
      otherV = mat.rowval[eInd]

      # clusterNeighborCount[vertexToCluster[otherV]] += edgeWeights[eInd]
      if haskey(clusterNeighborCount, vertexToCluster[otherV])
        clusterNeighborCount[vertexToCluster[otherV]] += 1
      else
        clusterNeighborCount[vertexToCluster[otherV]] = 1
      end #if

      if clusterNeighborCount[vertexToCluster[otherV]] > newClusterCount
        newClusterCount = clusterNeighborCount[vertexToCluster[otherV]]
        newCluster = vertexToCluster[otherV]
      end

      # if any of neighbors
        # 1: in same cluster
        # 2: have no other neighbors in that cluster
        #   don’t shuffle vertex
      if vertexToCluster[otherV] == vertexToCluster[vInd]
        hasOtherNeighborsInOwnCluster = false
        for otherEInd in mat.colptr[otherV]:(mat.colptr[otherV+1]-1)
          otherOtherV = mat.rowval[otherEInd]
          if otherOtherV != vInd && vertexToCluster[otherOtherV] == vertexToCluster[vInd]
            hasOtherNeighborsInOwnCluster = true
          end #if
        end #for
        if !hasOtherNeighborsInOwnCluster
          validVertex = false
        end
      end #if

    end #for

    if !validVertex
      continue
    end

    if newCluster != vertexToCluster[vInd]
      
      #remove vInd from it's old clusterList
      partitionList[vertexToCluster[vInd]][vertexToClusterLocation[vInd]] = -1

      #push to new clusterList
      push!(partitionList[newCluster], vInd)

      # change in vertexToCluster
      vertexToCluster[vInd] = newCluster
    end #if
  end #for

  for i in 1:length(partitionList)
    sort!(partitionList[i])
    for j in 1:length(partitionList[i])
      if partitionList[i][j] != -1
        vertexToClusterLocation[partitionList[i][j]] = j
      end
    end #for
  end #for

  return partitionList, vertexToCluster, starts, vertexToClusterLocation
end #reshuffleClusters


# takes a matrix from which partition will be constructed (referred to as the parent matrix)
# a partition is a double array of nodes. each sub-array c (cluster) has n vertices (the n of the parent mat)
# and either a 1 or 0 indicating whether or not v is in cluster c
# this function also returns an array of edge indices (of the parent mat) with the edges in the tree added 



# only shuffle vertices at the last level of BFS cluster Tree
function partitionMatrix(
  mat::SparseMatrixCSC,
  bigIteration,
  edgeClasses,
  bigEdgeMapReversed,
  bigMatNVertices,
  rows,
  columns,
  edgeWeights,
  randomClusters,
  shuffleClusters,
  exponentialX
)

  nVertices = mat.n
  partitionList = Array{Int64, 1}[]
  nClusters = 0
  visited = zeros(Bool, nVertices)
  nAddedToClusters = 0
  vertexToCluster = zeros(Int64, nVertices)
  vertexToClusterLocation = zeros(Int64, nVertices) #returns the index where the vertex is located inside clusterList
  starts = Int64[]
  newVertices = zeros(Int64, nVertices)
  nNewVertices = 0
  newVerticesArray = zeros(Bool, nVertices)

  if randomClusters
    unsortedNodes = unsortedArray(nVertices)
  end

  finalRoundClusterVertices = zeros(Bool, nVertices)

  x::Float64

  if nVertices > 2
    if exponentialX
      x = exp(sqrt(log(nVertices) * log(log(nVertices))))
    else
      x = log(nVertices+1)/log(2)
    end
  else
    x = 5.0
  end

  lowestUnvisitedUnsortedIndex = 1

  while (nAddedToClusters != nVertices)
    nClusters += 1
    push!(partitionList, Int64[])

    # pick a random starting node between 1 and nVertices each time
    if randomClusters
      while visited[unsortedNodes[lowestUnvisitedUnsortedIndex]]
        lowestUnvisitedUnsortedIndex += 1
      end #while
      start = unsortedNodes[lowestUnvisitedUnsortedIndex]
    else
      # OR (if not random), just select next unvisited node
      while visited[lowestUnvisitedUnsortedIndex]
        lowestUnvisitedUnsortedIndex += 1
      end #while
      start = lowestUnvisitedUnsortedIndex
    end

    push!(starts, start)
    visited[start] = true

    nNodesInCluster = 1
    nAddedToClusters += 1

    push!(partitionList[nClusters], start)
    vertexToCluster[start] = nClusters
    vertexToClusterLocation[start] = nNodesInCluster

    lastnNodesInCluster = 0

    calculatedVolume = 0.0
    calculatedBoundary = 0.0

    for eInd in mat.colptr[start]:(mat.colptr[start+1]-1)
      calculatedVolume += 1.0
      calculatedBoundary += 1.0
    end #for


    #building a single cluster
    while !(nAddedToClusters == nVertices)

      if (calculatedVolume / calculatedBoundary >= x) || lastnNodesInCluster == nNodesInCluster
        break
      end #if

      lastnNodesInCluster = nNodesInCluster

      for i in 1:nNewVertices
        newVerticesArray[newVertices[i]] = 0
      end #for

      nNewVertices = 0
      for v in partitionList[nClusters] #for each vertex in the current cluster
        for eInd in mat.colptr[v]:(mat.colptr[v+1]-1) #eInd is the edge

          if edgeClasses[bigEdgeMapReversed[eInd]] <= bigIteration
            otherV = mat.rowval[eInd] #this gives us the other vertex
            #doing this step before for boundary... maybe can combine?
            if !visited[otherV]
              if !newVerticesArray[otherV]
                nNewVertices += 1
                newVertices[nNewVertices] = otherV
              end #if
              newVerticesArray[otherV] = true
            end #if
          end #if

        end #for
      end #for

      for nthNewVertex in 1:nNewVertices
        newV = newVertices[nthNewVertex]
        push!(partitionList[nClusters], newV)
        vertexToCluster[newV] = nClusters
        nAddedToClusters += 1
        nNodesInCluster += 1
        vertexToClusterLocation[newV] = nNodesInCluster
        visited[newV] = true

        #updating boundary and volume
        for eInd in mat.colptr[newV]:(mat.colptr[newV+1]-1)
          otherV = mat.rowval[eInd]
          if vertexToCluster[otherV] == vertexToCluster[newV]
            calculatedBoundary -= 1.0

            #LEAVE THIS IN IF Vol = sum of degrees of nodes. Take out of Vol = # edges
            calculatedVolume += 1.0 
          else
            calculatedVolume += 1.0
            calculatedBoundary += 1.0
          end #if
        end #for

      end #for

    end #while

    for nthNewVertex in 1:nNewVertices
      finalRoundClusterVertices[newVertices[nthNewVertex]] = true
    end #for 

  end #while

  if shuffleClusters
    return reshuffleClusters(mat, partitionList, vertexToCluster, starts, vertexToClusterLocation, finalRoundClusterVertices, rows, columns, edgeWeights)
  end

  return partitionList, vertexToCluster, starts, vertexToClusterLocation
end #partitionMatrix


# for a single cluster, updates vertexToCluster to include another cluster for each extra component
# inside the original cluster
function updateVertexToCluster{Tv,Ti}(
  mat::SparseMatrixCSC{Tv,Ti},
  nVertices,
  vertexToCluster,
  c,
  nClusters,
  partitionList
)

  nVerticesInCluster = length(partitionList[c]) 
  visited = zeros(Bool, nVertices)
  compList = Array{Int64, 1}[]
  order = Array(Int64, nVertices)

  numComps = 0

  for clusterInd in 1:nVerticesInCluster
    vInd = partitionList[c][clusterInd]
    if !visited[vInd]
      numComps += 1
      visited[vInd] = true
      if numComps > 1
        vertexToCluster[vInd] = nClusters + numComps-1
      end

      if mat.colptr[vInd+1] > mat.colptr[vInd]
        ptr = 1
        orderLen = 2
        order[ptr] = vInd

        while ptr < orderLen
          curNode = order[ptr]

          for ind in mat.colptr[curNode]:(mat.colptr[curNode+1]-1)
            nbr = mat.rowval[ind]
            if vertexToCluster[nbr] == c && !visited[nbr]
              visited[nbr] = true
              if numComps > 1
                vertexToCluster[nbr] = nClusters + numComps-1
              end
              order[orderLen] = nbr
              orderLen += 1
            end # if
          end # for
          ptr += 1
        end # while
      end # if
    end
  end #for

  if numComps > 1
    nClusters += numComps-1
  end

  return nClusters
  # return maximum(vertexToCluster)

end #clusterComponents



function partitionListFromMap(nVertices, nClusters, vertexToCluster)
  vertexToClusterLocation = zeros(Int64, nVertices)
  starts = zeros(Int64, nClusters)
  partitionList = Array{Int64, 1}[]

  for i in 1:nClusters
    push!(partitionList, Int64[])
  end

  for i in 1:nVertices
    thisCluster = vertexToCluster[i]
    push!(partitionList[thisCluster], i)
    nVerticesInCluster = length(partitionList[thisCluster])
    vertexToClusterLocation[i] = nVerticesInCluster
    if nVerticesInCluster == 1
      starts[thisCluster] = i
    end #if
  end

  return vertexToClusterLocation, starts, partitionList
end #partitionListFromMap


function metisPartition(mat, nClusters)
  nVertices = mat.n
  g = Graph(mat)

  if nVertices == 2 || nClusters >= nVertices || nClusters <= 2
    vertexToCluster = ones(Int64, nVertices)
    nClusters = 1
  else 
    vertexToCluster = partGraphRecursive(g, nClusters)[2]
    while minimum(vertexToCluster) > 1
      vertexToCluster = vertexToCluster.- 1
      nClusters -= 1
    end
  end #if

  vertexToClusterLocation, starts, partitionList = partitionListFromMap(nVertices, nClusters, vertexToCluster)
  for c in 1:nClusters
    nClusters = updateVertexToCluster(mat, nVertices, vertexToCluster, c, nClusters, partitionList) #what if cluster comps is a list of lists?
  end

  vertexToClusterLocation, starts, partitionList = partitionListFromMap(nVertices, nClusters, vertexToCluster)

  return partitionList, vertexToCluster, starts, vertexToClusterLocation, nClusters
end #metisPartition



#definitions for a potentially faster lightgraphs heap

# immutable DijkstraHeapEntry{T}
#   vertex::Int
#   dist::T
# end

# isless(e1::DijkstraHeapEntry, e2::DijkstraHeapEntry) = e1.dist < e2.dist


function shortestPathsForCluster(mat, clusterList, vertexToCluster, start, vertexToClusterLocation)
  nVerticesInMat = mat.n
  nVerticesInCluster = length(clusterList)

  visited = zeros(Bool, nVerticesInCluster) 

  nh = intHeap(nVerticesInCluster) 
  dists = nh.keys

  pArray = zeros(Int64, nVerticesInCluster)
  pArrayNonZeros = Int64[]

  # sizehint!(nh, nVerticesInCluster)

  intHeapAdd!(nh, vertexToClusterLocation[start], 0.0)

  pArray[vertexToClusterLocation[start]] = vertexToClusterLocation[start]

  while nh.nitems > 0
    vIndInCluster::Int64 = intHeapPop!(nh)
    visited[vIndInCluster] = true

    dv = dists[vIndInCluster]
    vInd = clusterList[vIndInCluster]

    if vInd == -1
      continue
    end

    for eInd in mat.colptr[vInd]:(mat.colptr[vInd+1]-1)
      otherVInd = mat.rowval[eInd]
      otherVIndInCluster = vertexToClusterLocation[otherVInd]

      if vertexToCluster[otherVInd] == vertexToCluster[vInd] && !visited[otherVIndInCluster]
        newdist = dv + mat.nzval[eInd]
        if newdist < dists[otherVIndInCluster]
          dists[otherVIndInCluster] = newdist
          intHeapAdd!(nh, otherVIndInCluster, newdist)
          pArray[otherVIndInCluster] = vIndInCluster
          push!(pArrayNonZeros, otherVIndInCluster)
        end # if
      end # if
    end # for

  end # while




  #trying the lightgraphs heap

  # dists = fill(Inf, nVerticesInCluster)
  # visited = zeros(Bool, nVerticesInCluster) 
  # pArray = zeros(Int64, nVerticesInCluster)
  # H = Vector{Float64}()
  # dists[vertexToClusterLocation[start]] = 0.0

  # sizehint!(H, nVerticesInCluster)


  # heappush!(H, DijkstraHeapEntry(vertexToClusterLocation[start], dists[vertexToClusterLocation[start]]))


  # while !isempty(H)
  #   hentry = heappop!(H)
  #   vIndInCluster = hentry.vertex
  #   vInd = clusterList[vIndInCluster]
  #   dv = dists[vIndInCluster]

  #   if vInd == -1
  #     continue
  #   end
  #   # vIndInCluster
  #   for eInd in mat.colptr[vInd]:(mat.colptr[vInd+1]-1)
  #     otherVInd = mat.rowval[eInd]
  #     otherVIndInCluster = vertexToClusterLocation[otherVInd]

  #     if !visited[otherVIndInCluster]
  #       dists[otherVIndInCluster] = dv
  #       pArray[otherVIndInCluster] = vIndInCluster
  #       visited[otherVIndInCluster] = true
  #       heappush!(H, DijkstraHeapEntry{T}(otherVIndInCluster, alt))
  #     else
  #       if dv < dists[otherVIndInCluster]
  #         dists[otherVIndInCluster] = dv
  #         pArray[otherVIndInCluster] = vIndInCluster
  #         heappush!(H, DijkstraHeapEntry{T}(otherVIndInCluster, dv))
  #       end
  #       if dv == dists[otherVIndInCluster]
  #           pathcounts[otherVIndInCluster] += pathcounts[vIndInCluster]
  #       end #if
  #     end #if/else
  #   end #for
  # end #while










  newTreeInds = Int64[]

  for vInd in pArrayNonZeros
    bigMatVInd = clusterList[vInd]
    for eInd in mat.colptr[bigMatVInd]:(mat.colptr[bigMatVInd+1]-1) #all edges connecting to vertex j
      otherV = mat.rowval[eInd]
      clusterIndOtherV = vertexToClusterLocation[otherV]
      if clusterIndOtherV == pArray[vInd] && vertexToCluster[otherV] == vertexToCluster[bigMatVInd]
        push!(newTreeInds, eInd)
      end #if
    end #for
  end #for 


  # helpful for debugging a bug that causes a non-connected tree

  # checking connectedness
  # for vInd in 1:length(visited)
  #   if !visited[vInd] && clusterList[vInd] != -1
  #     println("for cluster: ", vertexToCluster[clusterList[vInd]], " vertex: ", clusterList[vInd], " is not connected")
  #     println("vertex: ", clusterList[vInd], " is connected to: ")
  #     for eInd in mat.colptr[clusterList[vInd]]:(mat.colptr[clusterList[vInd]+1]-1)
  #       otherVInd = mat.rowval[eInd]
  #       println("\t otherV: ", otherVInd, " in cluster: ", vertexToCluster[otherVInd])
  #     end
  #   end #if
  # end #for

  return newTreeInds
end # shortestPaths


function collapsePartition(mat::SparseMatrixCSC, partitionList, map, nEdges, rows, columns, edgeWeights)

  nClusters = length(partitionList)
  newRows = Int64[]
  newColumns = Int64[]
  newEdgeWeights = Float64[]
  for i in 1:nEdges
    if map[rows[i]] != map[columns[i]]
      push!(newRows, map[rows[i]]) 
      push!(newColumns, map[columns[i]])
      push!(newEdgeWeights, edgeWeights[i])
    end #if
  end #for

  newMat = sparse(newRows, newColumns, newEdgeWeights, nClusters, nClusters, min)

  return newMat
end #collapsePartition


function getBigEdgeMapReversed(
  mat,
  newMat,
  bigMapD,
  nEdges,
  rows,
  columns,
  edgeWeights,
  newNEdges,
  newRows,
  newColumns,
  newEdgeWeights
)

  bigEdgeMapReversed = zeros(Int64, newNEdges)

  for eInd in 1:nEdges
    c = columns[eInd]
    r = rows[eInd]

    # this makes this function MUCH MUCH faster -> usually O(E) rather than O(E^2)... sort of
    if bigMapD[c] == bigMapD[r]
      continue
    end

    for e in newMat.colptr[bigMapD[c]]:(newMat.colptr[bigMapD[c]+1]-1) #eInd is the edge
      if newMat.rowval[e] == bigMapD[r] && (bigEdgeMapReversed[e] == 0 || edgeWeights[eInd] < edgeWeights[bigEdgeMapReversed[e]])
        bigEdgeMapReversed[e] = eInd
      end #if
    end #for
  end #for

  return bigEdgeMapReversed
end #getBigEdgeMapReversed

function denormalizeEdgeWeights(oldMinEdgeWeight, mat, kind, nEdges)
  for eInd in 1:nEdges
    mat.nzval[eInd] *= oldMinEdgeWeight
  end #for

  if (kind == :max)
    for eInd in 1:nEdges
      mat.nzval[eInd] \= 1.0
    end #for
  end #if
end #denormalizeEdgeWeights


function sparseMatrixFromTreeIndices(mat, treeInds)
  nVertices = mat.n

  rows, columns, edgeWeights = findnz(mat)

  nRows = Int64[]
  nColumns = Int64[]
  nEdgeWeights = Float64[]

  for i in 1:length(treeInds)
    if treeInds[i]
      push!(nRows, rows[i])
      push!(nColumns, columns[i])
      push!(nEdgeWeights, edgeWeights[i])
    end #if
  end #for

  tree = sparse(nRows,nColumns,nEdgeWeights,nVertices,nVertices, min)
  tree = tree + tree'
  return tree

end #sparseMatrixFromTreeIndices

"""
Constructs a low stretch tree using the Alon, Karp, Peleg, West algorithm. This version (akpw! instead of akpw)
modifies the graph slightly changing the edges weights, then changing them back, which may lead to floating
point imprecisions. akpw! is faster (about 10-20%), but akpw doesn't have float imprecisions.

The function has a few options:

kind: default is :max, which regards each edge weight as the inverse of its length (just like kruskal). If this is
  set to anything else (e.g. :min), it will regard edge weight as length

randomClusters: default is false. This means the partition function searches for the beginning of the next cluster
  in node order, rather than randomly choosing nodes. If this is set to false, it will
  randomly choose the next node. This slows down akpw, but may produce better stretch.

metisClustering: default is false. If this is set to false, the graph will be partitioned
  each time by metis, rather than by the akpw partitioning method.

shuffleClusters: default is true. This preserves the "reshuffleClusters" method after each each graph is
  partitioned into clusters. If set to false, the function will skip this step. May be faster
  but have worse stretch.

exponentialX: default is true, where the funciton exp(sqrt(log(nVertices) * log(log(nVertices)))) is used for X.
  If set fo false, the function log(nVertices+1)/log(2) will be used for X instead. 


EXAMPLE:

[2, 1]  =  0.631273
[3, 1]  =  0.40103
[1, 2]  =  0.631273
[4, 2]  =  0.147018
[1, 3]  =  0.40103
[4, 3]  =  0.772661
[2, 4]  =  0.147018
[3, 4]  =  0.772661

      |
      |
      V

[2, 1]  =  0.631273
[3, 1]  =  0.40103
[1, 2]  =  0.631273
[1, 3]  =  0.40103
[4, 3]  =  0.772661
[3, 4]  =  0.772661
"""
function akpw!(
  mat::SparseMatrixCSC;
  kind=:max,
  randomClusters= false,
  metisClustering= false,
  shuffleClusters= true,
  exponentialX= true
)

  # useful for debugging!
  # println("Starting up AKPW...")
  # info("Starting up AKPW...")

  nVertices = mat.n
  nEdges = nnz(mat)

  oldMinEdgeWeight = normalizeEdgeWeights(mat, kind, nEdges)
  rows, columns, edgeWeights = findnz(mat)

  newMat = mat
  newNEdges = nEdges
  newRows, newColumns, newEdgeWeights = rows, columns, edgeWeights

  bigMapD = zeros(Int64, nVertices) #assuming map does C->D, bigMapD does A->D.
  for i in 1:nVertices #the first iteration maps to iself
    bigMapD[i] = i
  end #for

  bigEdgeMapReversed = zeros(Int64, nEdges) #goes C->A with edges (picks min of the many possibilities)
  for i in 1:nEdges #the first iteration maps to iself
    bigEdgeMapReversed[i] = i
  end #for

  treeInds = zeros(Bool, nEdges)
  edgeClasses = divideEdgesIntoClasses(mat, nEdges, rows, columns, edgeWeights)

  bigIteration = 1
  while (true)
    if !metisClustering
      partitionList, vertexToCluster, starts, vertexToClusterLocation = partitionMatrix(newMat, bigIteration, edgeClasses, bigEdgeMapReversed, nVertices, newRows, newColumns, newEdgeWeights, randomClusters, shuffleClusters, exponentialX)
      nClusters = length(partitionList)
    else
      nClusters = Int64(floor((newMat.n)/900)+2)
      # nClusters = Int64(floor(log((newMat.n)+1)/log(2)))
      partitionList, vertexToCluster, starts, vertexToClusterLocation, nClusters = metisPartition(newMat, nClusters)
    end

    # useful for debugging!
    # println("nClusters ", nClusters)

    for i in 1:length(bigMapD)
      bigMapD[i] = vertexToCluster[bigMapD[i]]
    end #for


    #this, for each cluster, creates its own sparse matrix and runs shortest paths on that, then returns the 
    # original mapped edges corresponding to that shortest paths tree. Then it adds those indices to treeInds
    for cInd in 1:nClusters
      newTreeInds = shortestPathsForCluster(newMat, partitionList[cInd], vertexToCluster, starts[cInd], vertexToClusterLocation)

      for i in newTreeInds
        treeInds[bigEdgeMapReversed[i]] = true
      end #for
    end #for

    if (nClusters == 1)
      break
    end #if

    newMat = collapsePartition(newMat, partitionList, vertexToCluster, newNEdges, newRows, newColumns, newEdgeWeights)
    newNEdges = nnz(newMat)
    newRows, newColumns, newEdgeWeights = findnz(newMat)

    bigEdgeMapReversed = getBigEdgeMapReversed(mat, newMat, bigMapD, nEdges, rows, columns, edgeWeights, newNEdges, newRows, newColumns, newEdgeWeights)

    bigIteration += 1
  end #while

  denormalizeEdgeWeights(oldMinEdgeWeight, mat, kind, nEdges)

  finalTree = sparseMatrixFromTreeIndices(mat, treeInds)

  return finalTree
end #akpw!

"""
This is a wrapper for akpw!. It's slower, but won't modify the original graph. See akpw! documentation for more
details.
"""
function akpw(
  origMat::SparseMatrixCSC;
  kind=:max,
  randomClusters=false,
  metisClustering=false,
  shuffleClusters=true,
  exponentialX=true
)

  rows, columns, edgeWeights = findnz(origMat)
  mat = sparse(rows, columns, edgeWeights)

  return akpw!(mat, kind=kind, randomClusters=randomClusters, metisClustering=metisClustering, shuffleClusters=shuffleClusters, exponentialX=exponentialX)
end #akpw



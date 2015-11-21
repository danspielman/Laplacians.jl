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
    # println(j)
    permutation[i] = permutation[j]
    permutation[j] = i
  end #for

  return permutation
end #unsortedArray


#could speed this up by creating linked list for clusterCount (cuz it'll usually be mostly zeros)
#this makes partitionQueue potentially unsorted. Will that mess shit up..?
# what if this removes a vertex that's crucial for the connection of a cluster? Then it'll divide into 2 clusters.
#   but my code can't handle that... ------ suddenly they won't be connected, right?
# what if it was the starting vertex for a cluster..?
# should be greatest total weight! not just the greatest number of neighbors
function reshuffleClusters(mat, partitionQueue, vertexToCluster, starts, vertexToClusterLocation, finalRoundClusterVertices, rows, columns, edgeWeights)
  nVertices = mat.n
  visited = zeros(Bool, nVertices)
  nClusters = length(partitionQueue)
  clusterNeighborCount = zeros(Int64, nClusters)

  for vInd in 1:nVertices
    validVertex = true
    if starts[vertexToCluster[vInd]] == vInd || !finalRoundClusterVertices[vInd]
      continue
    end #if

    for i in 1:nClusters
      clusterNeighborCount[i] = 0
    end #for

    newCluster = vertexToCluster[vInd]
    newClusterCount = 0

    for eInd in mat.colptr[vInd]:(mat.colptr[vInd+1]-1) #eInd is the edge
      otherV = mat.rowval[eInd]

      # clusterNeighborCount[vertexToCluster[otherV]] += edgeWeights[eInd]
      clusterNeighborCount[vertexToCluster[otherV]] += 1

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
      
      #remove vInd from it's old clusterQueue
      partitionQueue[vertexToCluster[vInd]][vertexToClusterLocation[vInd]] = -1

      #push to new clusterQueue
      push!(partitionQueue[newCluster], vInd)

      # change in vertexToCluster
      vertexToCluster[vInd] = newCluster
    end #if
  end #for

  for i in 1:length(partitionQueue)
    sort!(partitionQueue[i])
    for j in 1:length(partitionQueue[i])
      if partitionQueue[i][j] != -1
        vertexToClusterLocation[partitionQueue[i][j]] = j
      end
    end #for
  end #for

  return partitionQueue, vertexToCluster, starts, vertexToClusterLocation
end #reshuffleClusters


# takes a matrix from which partition will be constructed (referred to as the parent matrix)
# a partition is a double array of nodes. each sub-array c (cluster) has n vertices (the n of the parent mat)
# and either a 1 or 0 indicating whether or not v is in cluster c
# this function also returns an array of edge indices (of the parent mat) with the edges in the tree added 



# only shuffle vertices at the last level of BFS cluster Tree
function partitionMatrix(mat::SparseMatrixCSC, bigIteration, edgeClasses, bigEdgeMapReversed, bigMatNVertices, rows, columns, edgeWeights, randomClusters, shuffleClusters)
  nVertices = mat.n
  partitionQueue = Array{Int64, 1}[]
  nClusters = 0
  visited = zeros(Bool, nVertices)
  nAddedToClusters = 0
  vertexToCluster = zeros(Int64, nVertices)
  vertexToClusterLocation = zeros(Int64, nVertices) #returns the index where the vertex is located inside clusterQueue
  starts = Int64[]
  newVertices = zeros(Int64, nVertices)
  nNewVertices = 0
  newVerticesArray = zeros(Bool, nVertices)

  if (randomClusters != :no)
    unsortedNodes = unsortedArray(nVertices)
  end

  finalRoundClusterVertices = zeros(Bool, nVertices)

  x::Float64

  if nVertices > 2
      # x = exp(sqrt(log(nVertices) * log(log(nVertices))))
      x = log(nVertices+1)/log(2)
      # x = log(nVertices)
      # x = 2.0*exp(sqrt(log(bigMatNVertices) * log(log(bigMatNVertices))))
  else
    x = 5.0
  end

  lowestUnvisitedUnsortedIndex = 1

  while (nAddedToClusters != nVertices)
    nClusters += 1
    push!(partitionQueue, Int64[])

    # pick a random starting node between 1 and nVertices each time
    if (randomClusters != :no)
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

    push!(partitionQueue[nClusters], start)
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
      for v in partitionQueue[nClusters] #for each vertex in the current cluster
        for eInd in mat.colptr[v]:(mat.colptr[v+1]-1) #eInd is the edge

          if edgeClasses[bigEdgeMapReversed[eInd]] <= bigIteration
            # println("HELLO")
            otherV = mat.rowval[eInd] #this gives us the other vertex
            #doing this step before for boundary... maybe can combine?
            if !visited[otherV]
              if !newVerticesArray[otherV]
                nNewVertices += 1
                newVertices[nNewVertices] = otherV
                # push!(newVertices, otherV)
              end #if
              newVerticesArray[otherV] = true
            end #if
          end #if

        end #for
      end #for

      for nthNewVertex in 1:nNewVertices
        newV = newVertices[nthNewVertex]
        push!(partitionQueue[nClusters], newV)
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

  if shuffleClusters == :yes
    return reshuffleClusters(mat, partitionQueue, vertexToCluster, starts, vertexToClusterLocation, finalRoundClusterVertices, rows, columns, edgeWeights)
  end

  return partitionQueue, vertexToCluster, starts, vertexToClusterLocation
end #partitionMatrix


function metisPartition(mat, nClusters)
  nVertices = mat.n
  g = Graph(mat)

  if nVertices == 2
    vertexToCluster = ones(Int64, nVertices)
    nClusters = 1
  else 
    vertexToCluster = partGraphRecursive(g, nClusters)[2]
    while minimum(vertexToCluster) > 1
      vertexToCluster = vertexToCluster.- 1
      nClusters -= 1
    end
  end #if

  vertexToClusterLocation = zeros(Int64, nVertices)
  starts = zeros(Int64, nClusters)
  partitionQueue = Array{Int64, 1}[]

  for i in 1:nClusters
    push!(partitionQueue, Int64[])
  end


  for i in 1:nVertices
    thisCluster = vertexToCluster[i]
    push!(partitionQueue[thisCluster], i)
    nVerticesInCluster = length(partitionQueue[thisCluster])
    vertexToClusterLocation[i] = nVerticesInCluster
    if nVerticesInCluster == 1
      starts[thisCluster] = i
    end #if
  end #for

  for i in 1:nClusters
    if starts[i] == 0
      println("************start[", i, "] = 0")
      println(partitionQueue[i])
    end
  end


  return partitionQueue, vertexToCluster, starts, vertexToClusterLocation, nClusters
end #metisPartition



function shortestPathsForCluster(mat, clusterQueue, vertexToCluster, start, vertexToClusterLocation)
  nVerticesInMat = mat.n
  nVerticesInCluster = length(clusterQueue)

  visited = zeros(Bool, nVerticesInCluster) 

  nh = intHeap(nVerticesInCluster) 
  dists = nh.keys

  pArray = zeros(Int64, nVerticesInCluster)
  pArrayNonZeros = Int64[]

  intHeapAdd!(nh, vertexToClusterLocation[start], 0.0)

  pArray[vertexToClusterLocation[start]] = vertexToClusterLocation[start]

  while nh.nitems > 0
    vIndInCluster::Int64 = intHeapPop!(nh)
    visited[vIndInCluster] = true

    dv = dists[vIndInCluster]
    vInd = clusterQueue[vIndInCluster]

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

  newTreeInds = Int64[]

  for vInd in pArrayNonZeros
    bigMatVInd = clusterQueue[vInd]
    for eInd in mat.colptr[bigMatVInd]:(mat.colptr[bigMatVInd+1]-1) #all edges connecting to vertex j
      otherV = mat.rowval[eInd]
      clusterIndOtherV = vertexToClusterLocation[otherV]
      if clusterIndOtherV == pArray[vInd] && vertexToCluster[otherV] == vertexToCluster[bigMatVInd]
        push!(newTreeInds, eInd)
      end #if
    end #for
  end #for 

  # checking connectedness
  for vInd in 1:length(visited)
    if !visited[vInd]
      if clusterQueue[vInd] == -1
        continue
      end
      println("for cluster: ", vertexToCluster[clusterQueue[vInd]], " vertex: ", clusterQueue[vInd], " is not connected")
      println("vertex: ", clusterQueue[vInd], " is connected to: ")
      for eInd in mat.colptr[clusterQueue[vInd]]:(mat.colptr[clusterQueue[vInd]+1]-1)
        otherVInd = mat.rowval[eInd]
        println("\t otherV: ", otherVInd, " in cluster: ", vertexToCluster[otherVInd])
      end
    end #if
  end #for

  return newTreeInds
end # shortestPaths


function collapsePartition(mat::SparseMatrixCSC, partitionQueue, map, nEdges, rows, columns, edgeWeights)

  nClusters = length(partitionQueue)
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


function getBigEdgeMapReversed(mat, newMat, bigMapD, nEdges, rows, columns, edgeWeights, newNEdges, newRows, newColumns, newEdgeWeights)

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


function akpw(mat::SparseMatrixCSC; kind=:max, randomClusters=:no, metisClustering=:no, shuffleClusters=:yes)
  println("Starting up AKPW...")

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
    if (metisClustering == :no)
      partitionQueue, vertexToCluster, starts, vertexToClusterLocation = partitionMatrix(newMat, bigIteration, edgeClasses, bigEdgeMapReversed, nVertices, newRows, newColumns, newEdgeWeights, randomClusters, shuffleClusters)
      nClusters = length(partitionQueue)
    else
      nClusters = Int64(floor(log((newMat.n)+1)/log(2)))
      partitionQueue, vertexToCluster, starts, vertexToClusterLocation, nClusters = metisPartition(newMat, nClusters)
    end

    println("nClusters ", nClusters)

    for i = 1:length(bigMapD)
      bigMapD[i] = vertexToCluster[bigMapD[i]]
    end #for


    #this, for each cluster, creates its own sparse matrix and runs shortest paths on that, then returns the 
    # original mapped edges corresponding to that shortest paths tree. Then it adds those indices to treeInds
    for cInd in 1:nClusters
      newTreeInds = shortestPathsForCluster(newMat, partitionQueue[cInd], vertexToCluster, starts[cInd], vertexToClusterLocation)

      for i in newTreeInds
        treeInds[bigEdgeMapReversed[i]] = true
      end #for
    end #for

    if (nClusters == 1)
      break
    end #if

    newMat = collapsePartition(newMat, partitionQueue, vertexToCluster, newNEdges, newRows, newColumns, newEdgeWeights)
    newNEdges = nnz(newMat)
    newRows, newColumns, newEdgeWeights = findnz(newMat)

    bigEdgeMapReversed = getBigEdgeMapReversed(mat, newMat, bigMapD, nEdges, rows, columns, edgeWeights, newNEdges, newRows, newColumns, newEdgeWeights)

    bigIteration += 1
  end #while

  denormalizeEdgeWeights(oldMinEdgeWeight, mat, kind, nEdges)

  finalTree = sparseMatrixFromTreeIndices(mat, treeInds)
  return finalTree
end #akpw



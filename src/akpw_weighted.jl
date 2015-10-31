#returns volume = sum of all the degrees of nodes in cluster
function getVolume(mat::SparseMatrixCSC, cluster)
  count = 0
  for v in cluster
    for ind in mat.colptr[v]:(mat.colptr[v+1]-1)
      count += 1
    end #for
  end #for
  return count
end #getVolume

#returns boundary = # edges leaving cluster
function getBoundary(mat::SparseMatrixCSC, cluster, vertexToCluster)
  count = 0
  for v in cluster
    for ind in mat.colptr[v]:(mat.colptr[v+1]-1)
      if vertexToCluster[v] != vertexToCluster[mat.rowval[ind]]
        count += 1
      end #if
    end #for
  end #for
  return count
end #getBoundary


function normalizeEdgeWeights(mat::SparseMatrixCSC, kind)
  rows, columns, edgeWeights = findnz(mat)
  if (kind == :max)
    edgeWeights = edgeWeights.\ 1
  end #if

  minEdgeWeight = minimum(edgeWeights)
  edgeWeights = edgeWeights./ minEdgeWeight

  newMat = sparse(rows, columns, edgeWeights)
  return minEdgeWeight, newMat
end #normalizeEdgeWeights


function divideEdgesIntoClasses(mat::SparseMatrixCSC)
  rows, columns, edgeWeights = findnz(mat)
  nVertices = mat.n
  nEdges = nnz(mat)
  x = exp(sqrt(log(nVertices) * log(log(nVertices))))
  y = x * log(nVertices) #placeholder for testing
  edgeClasses = zeros(Int64, nEdges)

  for eInd in 1:nEdges
    edgeClasses[eInd] = floor(log(y, edgeWeights[eInd])) + 1
  end #for
  

  return edgeClasses

end #divideEdgesIntoClasses


# takes a matrix from which partition will be constructed (referred to as the parent matrix)
# a partition is a double array of nodes. each sub-array c (cluster) has n vertices (the n of the parent mat)
# and either a 1 or 0 indicating whether or not v is in cluster c
# this function also returns an array of edge indices (of the parent mat) with the edges in the tree added 
function partitionMatrix(mat::SparseMatrixCSC, bigIteration, edgeClasses, bigEdgeMapReversed, bigMatNVertices)
  nVertices = mat.n
  partition = Array{Bool,1}[] #is this right..? it's an array of an array of bools...
  partitionQueue = Array{Int64, 1}[]
  nClusters = 0
  visited = zeros(Bool, nVertices)
  nAddedToClusters = 0
  vertexToCluster = zeros(Int64, nVertices)
  starts = Int64[]

  println(nVertices)
  if nVertices > 2
      # x = exp(sqrt(log(nNodesInCluster) * log(log(nNodesInCluster))))
      x = log(nVertices+1)/log(2)
      # x = log(nVertices)
      # x = 2.0*exp(sqrt(log(nNodesInCluster) * log(log(nNodesInCluster))))
      # x = 5.0
  else
    x = 5
  end
  println(x)

  # build all the clusters
  while (nAddedToClusters != nVertices)
    nClusters += 1
    push!(partitionQueue, Int64[])

    #select arbitrary unvisited node
    start = 1
    for node in 1:nVertices
      if !visited[node]
        start = node
        push!(starts, start)
        break
      end #if
    end #for
    visited[start] = true

    push!(partitionQueue[nClusters], start)
    vertexToCluster[start] = nClusters

    nNodesInCluster = 1
    nAddedToClusters += 1

    lastnNodesInCluster = 0

    #building a single cluster
    while !(nAddedToClusters == nVertices)

      volume = getVolume(mat, partitionQueue[nClusters])
      boundary = getBoundary(mat, partitionQueue[nClusters], vertexToCluster)

      if (convert(Float64, volume) / convert(Float64, boundary) >= x) || lastnNodesInCluster == nNodesInCluster
        break
      end #if

      # else 
      lastnNodesInCluster = nNodesInCluster

      newVertices = Int64[]
      newVerticesArray = zeros(Bool, nVertices)
      for v in partitionQueue[nClusters] #for each vertex in the current cluster
        for eInd in mat.colptr[v]:(mat.colptr[v+1]-1) #eInd is the edge

          if edgeClasses[bigEdgeMapReversed[eInd]] <= bigIteration
            otherV = mat.rowval[eInd] #this gives us the other vertex
            #doing this step before for boundary... maybe can combine?
            if !visited[otherV]
              if !newVerticesArray[otherV]
                push!(newVertices, otherV)
              end #if
              newVerticesArray[otherV] = true
            end #if
          end #if

        end #for
      end #for

      for v in newVertices
        push!(partitionQueue[nClusters], v)
        vertexToCluster[v] = nClusters
        nAddedToClusters += 1
        nNodesInCluster += 1
        visited[v] = true
      end #for

    end #while
  end #while

  return partitionQueue, vertexToCluster, starts#, treeInds

end #partitionMatrix


function shortestPathsForCluster(mat, clusterQueue, vertexToCluster, start)
  rows = Int64[]
  columns = Int64[]
  edgeWeights = Float64[]
  miniMap = Int64[]


  for vInd in clusterQueue
    for eInd in mat.colptr[vInd]:(mat.colptr[vInd+1]-1) #eInd is the edge
      otherV = mat.rowval[eInd] #this gives us the other vertex
      if vertexToCluster[otherV] == vertexToCluster[vInd]
        push!(rows, vInd)
        push!(columns, otherV)
        push!(edgeWeights, mat.nzval[eInd])
        push!(miniMap, eInd)
      end #if
    end #for
  end #for

  for eInd in 1:length(edgeWeights)
    edgeWeights[eInd] \= 1.0
  end #for

  if length(rows) < 1
    return Int64[]
  end #if

  clusterMat = sparse(rows, columns, edgeWeights)
  dists, pArray, pArrayNonZeros = shortestPaths(clusterMat, start)

  newTreeInds = Int64[]

  for vInd in pArrayNonZeros
    for eInd in mat.colptr[vInd]:(mat.colptr[vInd+1]-1) #all edges connecting to vertex j
      otherV = mat.rowval[eInd]
      if otherV == pArray[vInd]
        push!(newTreeInds, eInd)
      end #if
    end #for
  end #for 

  return newTreeInds
end #shortestPathsForCluster



function collapsePartition(mat::SparseMatrixCSC, partitionQueue, map)

  nClusters = length(partitionQueue)
  nEdges = nnz(mat)
  newRows = Int64[]
  newColumns = Int64[]
  newEdgeWeights = Float64[]
  rows, columns, edgeWeights = findnz(mat)
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


function getBigEdgeMapReversed(mat, newMat, bigMapD)
  nEdges = nnz(mat)
  rows, columns, edgeWeights = findnz(mat)

  newNEdges = nnz(newMat)
  newRows, newColumns, newEdgeWeights = findnz(newMat)

  bigEdgeMapReversed = zeros(Int64, newNEdges)

  for eInd in 1:nEdges
    c = columns[eInd]
    r = rows[eInd]
    for e in newMat.colptr[bigMapD[c]]:(newMat.colptr[bigMapD[c]+1]-1) #eInd is the edge
      otherV = newMat.rowval[e] #this gives us the other vertex
      if otherV == bigMapD[r] && (bigEdgeMapReversed[e] == 0 || edgeWeights[eInd] < edgeWeights[bigEdgeMapReversed[e]])
        bigEdgeMapReversed[e] = eInd
      end #if
    end #for
  end #for

  return bigEdgeMapReversed
end #getBigEdgeMapReversed






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


function deNormalizeEdges(mat, oldMinEdgeWeight, kind)
  rows, columns, edgeWeights = findnz(mat)
  edgeWeights = edgeWeights.* oldMinEdgeWeight
  if (kind == :max)
    edgeWeights = edgeWeights.\ 1
  end #if 

  newMat = sparse(rows, columns, edgeWeights)
  return newMat
end #deNormalizeEdges




function akpw(mat::SparseMatrixCSC; kind=:min)
  nVertices = mat.n
  nEdges = nnz(mat)
  println(mat)

  oldMinEdgeWeight, mat = normalizeEdgeWeights(mat, kind)
  newMat = mat

  bigMapD = zeros(Int64, nVertices) #assuming map does C->D, bigMapD does A->D.
  for i in 1:nVertices #the first iteration maps to iself
    bigMapD[i] = i
  end #for

  bigEdgeMapReversed = zeros(Int64, nEdges) #goes C->A with edges (picks min of the many possibilities)
  for i in 1:nEdges #the first iteration maps to iself
    bigEdgeMapReversed[i] = i
  end #for

  treeInds = zeros(Bool, nEdges)
  edgeClasses = divideEdgesIntoClasses(mat)
  println("classes ", edgeClasses)
  println(newMat)



  bigIteration = 1
  while (true)
    partitionQueue, vertexToCluster, starts = partitionMatrix(newMat, bigIteration, edgeClasses, bigEdgeMapReversed, nVertices)

    println("nClusters ", length(partitionQueue))

    for i = 1:length(bigMapD)
      bigMapD[i] = vertexToCluster[bigMapD[i]]
    end #for

    println("made map for ", bigIteration)

    # hangs here
    time = 0
    #this, for each cluster, creates its own sparse matrix and runs shortest paths on that, then returns the 
    # original mapped edges corresponding to that shortest paths tree. Then it adds those indices to treeInds
    for cInd in 1:length(partitionQueue)

      # tic()
      newTreeInds = shortestPathsForCluster(newMat, partitionQueue[cInd], vertexToCluster, starts[cInd])
      # time += toq()

      for i in newTreeInds
        treeInds[bigEdgeMapReversed[i]] = true
      end #for

    end #for
    println(time)

    if (length(partitionQueue) == 1)
      break
    end #if

    println("locked in paths for ", bigIteration)

    newMat = collapsePartition(newMat, partitionQueue, vertexToCluster)

    println("collapsed ", bigIteration)

    bigEdgeMapReversed = getBigEdgeMapReversed(mat, newMat, bigMapD)

    println("big edge map for  ", bigIteration)

    bigIteration += 1
  end #while


  finalTree = sparseMatrixFromTreeIndices(mat, treeInds)

  return deNormalizeEdges(finalTree, oldMinEdgeWeight, kind)
end #akpw

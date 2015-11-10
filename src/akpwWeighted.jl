#returns volume = sum of all the degrees of nodes in cluster
function getVolume(mat::SparseMatrixCSC, cluster)
  count = 0.0
  for v in cluster
    for ind in mat.colptr[v]:(mat.colptr[v+1]-1)
      count += 1.0
    end #for
  end #for
  return count
end #getVolume

#returns boundary = # edges leaving cluster
function getBoundary(mat::SparseMatrixCSC, cluster, vertexToCluster)
  count = 0.0
  for v in cluster
    for ind in mat.colptr[v]:(mat.colptr[v+1]-1)
      if vertexToCluster[v] != vertexToCluster[mat.rowval[ind]]
        count += 1.0
      end #if
    end #for
  end #for
  return count
end #getBoundary


function normalizeEdgeWeights(mat::SparseMatrixCSC, kind)

  # println(mat)
  rows, columns, edgeWeights = findnz(mat)
  if (kind == :max)
    edgeWeights = edgeWeights.\ 1.0
  end #if

  minEdgeWeight = minimum(edgeWeights)
  edgeWeights = edgeWeights./ minEdgeWeight

  newMat = sparse(rows, columns, edgeWeights)

  # println(newMat)

  # nVertices = mat.n

  # for vInd in 1:nVertices
  #   for eInd in mat.colptr[vInd]:(mat.colptr[vInd+1]-1)
  #     if (kind == :max)
  #       mat.rowval[eInd] \= 1.0
  #     end #if

  #     mat.rowval[eInd] /= minEdgeWeight
  #   end #for
  # end #for

  # newMat = sparse(rows, columns, edgeWeights)

  # println(mat)



  # for vInd in 1:nVertices
  #   for eInd in mat.colptr[vInd]:(mat.colptr[vInd+1]-1)
  #     if mat[eInd] != newMat[eInd]
  #       println("FALSE")
  #       break
  #     end #if
  #   end #for
  # end #for

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





#could speed this up by creating linked list for clusterCount (cuz it'll usually be mostly zeros)
#this makes partitionQueue potentially unsorted. Will that mess shit up..?
# what if this removes a vertex that's crucial for the connection of a cluster? Then it'll divide into 2 clusters.
#   but my code can't handle that... ------ suddenly they won't be connected, right?
# what if it was the starting vertex for a cluster..?
# should be greatest total weight! not just the greatest number of neighbors
function reshuffleClusters(mat, partitionQueue, vertexToCluster, starts, vertexToClusterLocation, finalRoundClusterVertices)
  nVertices = mat.n
  visited = zeros(Bool, nVertices)
  nClusters = length(partitionQueue)
  clusterNeighborCount = zeros(Float64, nClusters)
  rows, columns, edgeWeights = findnz(mat)

  # println(mat)

  # println("partitionQueue: ", partitionQueue)
  # println("vertexToCluster: ", vertexToCluster)
  # println("starts: ", starts)
  # println("vertexToClusterLocation: ", vertexToClusterLocation)


  for vInd in 1:nVertices
    # println("vertex: ", vInd)
    if starts[vertexToCluster[vInd]] == vInd || !finalRoundClusterVertices[vInd]
      # println("is a start/not final round: ", vInd)
      continue
    end #if

    for i in 1:nClusters
      clusterNeighborCount[i] = 0
    end #for

    for eInd in mat.colptr[vInd]:(mat.colptr[vInd+1]-1) #eInd is the edge
      otherV = mat.rowval[eInd]
      clusterNeighborCount[vertexToCluster[otherV]] += edgeWeights[eInd]
    end #for

    # println("vertex: ", vInd)
    # println("cluster count: ",clusterNeighborCount)

    newCluster = indmax(clusterNeighborCount)

    if newCluster != vertexToCluster[vInd]
      

      # println("SWITCHING Vertex: ", vInd, " from cluster ", vertexToCluster[vInd], " to ", newCluster)
      #remove from old clusterQueue
      oldClusterSubQueue = Int64[]
      for i in partitionQueue[vertexToCluster[vInd]]
        if i != vInd
          push!(oldClusterSubQueue, i)
        end #if
      end #for

      # println(partitionQueue)
      # println(oldClusterSubQueue)
      # println(vertexToCluster[vInd])

      partitionQueue[vertexToCluster[vInd]] = oldClusterSubQueue

      # println("HELLO")


      #push to new clusterQueue
      push!(partitionQueue[newCluster], vInd)

      # change in vertexToClusterLocation
      # vertexToClusterLocation[vInd] = length(partitionQueue[vertexToCluster[vInd]])

      # if was a start, change that
      # starts[vertexToCluster[vInd]] = clusterQueue[vertexToCluster[vInd]][2]..?

      # change in vertexToCluster
      vertexToCluster[vInd] = newCluster
    end #if
  end #for

  for i in 1:length(partitionQueue)
    sort!(partitionQueue[i])
    for j in 1:length(partitionQueue[i])
      vertexToClusterLocation[partitionQueue[i][j]] = j
    end #for
  end #for

  # partitionQueue

  # println("partitionQueue: ", partitionQueue)
  # println("vertexToCluster: ", vertexToCluster)
  # println("starts: ", starts)
  # println("vertexToClusterLocation: ", vertexToClusterLocation)


  return partitionQueue, vertexToCluster, starts, vertexToClusterLocation
end #reshuffleClusters


# takes a matrix from which partition will be constructed (referred to as the parent matrix)
# a partition is a double array of nodes. each sub-array c (cluster) has n vertices (the n of the parent mat)
# and either a 1 or 0 indicating whether or not v is in cluster c
# this function also returns an array of edge indices (of the parent mat) with the edges in the tree added 



# only shuffle vertices at the last level of BFS cluster Tree
function partitionMatrix(mat::SparseMatrixCSC, bigIteration, edgeClasses, bigEdgeMapReversed, bigMatNVertices)
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

  finalRoundClusterVertices = zeros(Bool, nVertices)

  x::Float64

  if nVertices > 2
      # x = exp(sqrt(log(nVertices) * log(log(nVertices))))
      x = log(nVertices+1)/log(2)
      # x = log(nVertices)
      # x = 2.0*exp(sqrt(log(nNodesInCluster) * log(log(nNodesInCluster))))
      # x = 5.0
  else
    x = 5.0
  end

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

    nNodesInCluster = 1
    nAddedToClusters += 1

    push!(partitionQueue[nClusters], start)
    vertexToCluster[start] = nClusters
    vertexToClusterLocation[start] = nNodesInCluster

    lastnNodesInCluster = 0

    #building a single cluster
    while !(nAddedToClusters == nVertices)

      volume = getVolume(mat, partitionQueue[nClusters])
      boundary = getBoundary(mat, partitionQueue[nClusters], vertexToCluster)

      if (volume / boundary >= x) || lastnNodesInCluster == nNodesInCluster
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

      # println("adding: ", newVertices," (", nNewVertices, ") to cluster: ", nClusters)

      for nthNewVertex in 1:nNewVertices
        push!(partitionQueue[nClusters], newVertices[nthNewVertex])
        vertexToCluster[newVertices[nthNewVertex]] = nClusters
        nAddedToClusters += 1
        nNodesInCluster += 1
        vertexToClusterLocation[newVertices[nthNewVertex]] = nNodesInCluster
        visited[newVertices[nthNewVertex]] = true
      end #for

      # for i in 1:nNewVertices
      #   newVerticesArray[newVertices[i]] = 0
      # end #for

    end #while

    for nthNewVertex in 1:nNewVertices
      finalRoundClusterVertices[newVertices[nthNewVertex]] = true
    end #for 

  end #while

  # println("finalRoundClusterVertices: ", finalRoundClusterVertices)

  # return reshuffleClusters(mat, partitionQueue, vertexToCluster, starts, vertexToClusterLocation, finalRoundClusterVertices)
  return partitionQueue, vertexToCluster, starts, vertexToClusterLocation
end #partitionMatrix



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
      if clusterIndOtherV == pArray[vInd]
        push!(newTreeInds, eInd)
      end #if
    end #for
  end #for 

  return newTreeInds
end # shortestPaths


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


function akpw(mat::SparseMatrixCSC; kind=:min)
  println("Starting up AKPW...")

  nVertices = mat.n
  nEdges = nnz(mat)

  unNormalizedMat = mat

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

  bigIteration = 1
  while (true)
    partitionQueue, vertexToCluster, starts, vertexToClusterLocation = partitionMatrix(newMat, bigIteration, edgeClasses, bigEdgeMapReversed, nVertices)

    println("nClusters ", length(partitionQueue))

    for i = 1:length(bigMapD)
      bigMapD[i] = vertexToCluster[bigMapD[i]]
    end #for

    
    # time = 0


    #this, for each cluster, creates its own sparse matrix and runs shortest paths on that, then returns the 
    # original mapped edges corresponding to that shortest paths tree. Then it adds those indices to treeInds
    for cInd in 1:length(partitionQueue)

      # tic()
      newTreeInds = shortestPathsForCluster(newMat, partitionQueue[cInd], vertexToCluster, starts[cInd], vertexToClusterLocation)
      # time += toq()

      for i in newTreeInds
        treeInds[bigEdgeMapReversed[i]] = true
      end #for

    end #for
    # println(time)

    if (length(partitionQueue) == 1)
      break
    end #if

    newMat = collapsePartition(newMat, partitionQueue, vertexToCluster)

    bigEdgeMapReversed = getBigEdgeMapReversed(mat, newMat, bigMapD)

    bigIteration += 1
  end #while


  finalTree = sparseMatrixFromTreeIndices(unNormalizedMat, treeInds)

  return finalTree
end #akpw










#to replace array of linkedLists (pushing to arrays)
# clusters = zeros(n, 1)
# clusters [1][2]...[n1][n1+1]...[n1+n2]... -> if n1 nodes in cluster 1, n2 nodes in cluster 2
  # saving time

#run against compute stretches code!










# type intHeap{Tkey,Tind}
#   keys::Array{Tkey,1}
#   heap::Array{Tind,1}
#   index::Array{Tind,1}
#   nitems::Tind
# end #intHeap

# intHeap(n::Int64) = intHeap(Inf*ones(Float64,n),-ones(Int64,n),zeros(Int64,n),0)
# intHeap(n::Int32) = intHeap(Inf*ones(Float32,n),-ones(Int32,n),zeros(Int32,n),0)

# function intHeapAdd!{Tkey,Tind}(nh::intHeap, node::Tind, key::Tkey)
#   if nh.index[node] > 0 # if already in the heap
#     if key < nh.keys[node]
#       intHeapSet!(nh, node, key)
#     end

#   else # if it really is new

#     nhp = nh.nitems+1

#     nh.keys[node] = key
#     nh.heap[nhp] = node
#     nh.index[node] = nhp
#     nh.nitems = nhp

#     intHeapUp!(nh, node)

#   end
# end # intHeapAdd!

# function intHeapDown!{Tind}(nh::intHeap, node::Tind)
#   pos = nh.index[node]
#   key = nh.keys[node]
#   leftPos = pos*2
#   moved = true
#   @inbounds while (leftPos <= nh.nitems) && moved
#     moved = false
#     rightPos = pos*2+1

#     if rightPos > nh.nitems
#       childPos = leftPos
#       childNode = nh.heap[childPos]
#       childKey = nh.keys[childNode]
#     else
#       leftNode = nh.heap[leftPos]
#       leftKey = nh.keys[leftNode]
#       rightNode = nh.heap[rightPos]
#       rightKey = nh.keys[rightNode]

#       if leftKey < rightKey
#         childPos = leftPos
#         childNode = leftNode
#         childKey = leftKey
#       else
#         childPos = rightPos
#         childNode = rightNode
#         childKey = rightKey
#       end
#     end

#     if childKey < key
#       nh.heap[childPos] = node
#       nh.heap[pos] = childNode
#       nh.index[node] = childPos
#       nh.index[childNode] = pos

#       pos = childPos
#       leftPos = pos*2
#       moved = true
#     end

#   end #while
# end # intHeapDown!

# function intHeapPop!(nh::intHeap)
#   minNode = nh.heap[1]

#   nh.index[minNode] = 0

#   @inbounds if (nh.nitems > 1)
#     node = nh.heap[nh.nitems]
#     nh.heap[1] = node
#     nh.index[node] = 1
#     intHeapDown!(nh, node)
#   end
#   nh.nitems = nh.nitems - 1

#   return minNode
# end # intHeapPop!

# function intHeapUp!{Tind}(nh::intHeap, node::Tind)
#   pos = nh.index[node]
#   moved = true

#   @inbounds while (pos > 1) && moved
#     key = nh.keys[node]

#     parentPos = div(pos,2)
#     parentNode = nh.heap[parentPos]
#     parentKey = nh.keys[parentNode]

#     moved = false

#     if (parentKey > key)
#       nh.heap[parentPos] = node
#       nh.heap[pos] = parentNode
#       nh.index[node] = parentPos
#       nh.index[parentNode] = pos
#       pos = parentPos
#       moved = true
#     end
#   end

# end # intHeapUp!

# function intHeapSort(x::Array{Float64,1})
#   n = length(x)
#   nh = intHeap(n)

#   @inbounds for i in 1:n
#     intHeapAdd!(nh, i, x[i])
#   end

#   out = zeros(Float64,n)
#   @inbounds for i in 1:n
#     out[i] = nh.keys[intHeapPop!(nh)]
#   end

#   return out

# end # intHeapSort


# function intHeapSort(nh::intHeap)
#   n = length(nh.keys)

#   out = zeros(Float64,n)
#   for i in 1:n
#     out[i] = nh.keys[intHeapPop!(nh)]
#   end

#   return out

# end # intHeapSort

# function intHeapSet!{Tkey,Tind}(nh::intHeap, node::Tind, key::Tkey)
#   oldKey = nh.keys[node]
#   nh.keys[node] = key

#   if (key < oldKey)
#     intHeapUp!(nh,node)
#   else
#     intHeapDown!(nh,node)
#   end
# end # intHeapSet!

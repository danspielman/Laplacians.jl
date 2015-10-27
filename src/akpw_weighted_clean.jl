# include("julia/yinsGraph/graphAlgs.jl")
# include("julia/yinsGraph/graphGenerators.jl")
# include("julia/yinsGraph/treeAlgs.jl")

using yinsGraph








#returns volume = sum of all the degrees of nodes in cluster
function getVolume(mat::SparseMatrixCSC, cluster)
  count = 0
  for v in 1:length(cluster)
    if cluster[v] == 1
      for ind in mat.colptr[v]:(mat.colptr[v+1]-1)
        count += 1
      end #for
    end #if
  end #for
  return count
end #getVolume

#returns boundary = # edges leaving cluster
function getBoundary(mat::SparseMatrixCSC, cluster)
  count = 0
  for v in 1:length(cluster)
    if cluster[v] == 1
      for ind in mat.colptr[v]:(mat.colptr[v+1]-1)
        if cluster[mat.rowval[ind]] == 0
          count += 1
        end #if
      end #for
    end #if
  end #for
  return count
end #getBoundary


function normalizeEdgeWeights(mat::SparseMatrixCSC)
  rows, columns, edgeWeights = findnz(mat)
  minEdgeWeight = minimum(edgeWeights)
  newEdgeWeights = zeros(Float64, length(edgeWeights))
  newEdgeWeights = edgeWeights./ minEdgeWeight
  # for i in 1:length(edgeWeights)
  #   newEdgeWeights[i] = edgeWeights[i] / minEdgeWeight
  # end #for

  newMat = sparse(rows, columns, newEdgeWeights)
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
function partitionMatrix(mat::SparseMatrixCSC, bigIteration, edgeClasses, bigEdgeMapReversed)
  nVertices = mat.n
  partition = Array{Bool,1}[] #is this right..? it's an array of an array of bools...
  nClusters = 0
  visited = zeros(Bool, nVertices)
  nAddedToClusters = 0
  vertexToCluster = zeros(Int64, nVertices)
  starts = Int64[]

  # build all the clusters
  while (nAddedToClusters != nVertices)
    nClusters += 1
    push!(partition, zeros(Bool, nVertices))

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

    partition[nClusters][start] = true
    vertexToCluster[start] = nClusters

    nNodesInCluster = 1
    nAddedToClusters += 1

    lastnNodesInCluster = 0

    #building a single cluster
    while !(nAddedToClusters == nVertices)
      volume = getVolume(mat, partition[nClusters])
      boundary = getBoundary(mat, partition[nClusters])

      if nNodesInCluster > 2
        x = exp(sqrt(log(nNodesInCluster) * log(log(nNodesInCluster))))
      else
        x = 10
      end

      if ((convert(Float64, volume)-convert(Float64, boundary)) / convert(Float64, boundary) >= x) #|| lastnNodesInCluster == nNodesInCluster
        break
      # end #if

      else 
        newVertices = zeros(Bool, nVertices)
        for v in 1:length(partition[nClusters]) #for each vertex in the current cluster
          if partition[nClusters][v] == true
            for eInd in mat.colptr[v]:(mat.colptr[v+1]-1) #eInd is the edge

              if edgeClasses[bigEdgeMapReversed[eInd]] <= bigIteration
                otherV = mat.rowval[eInd] #this gives us the other vertex
                #doing this step before for boundary... maybe can combine?
                if !visited[otherV]
                  newVertices[otherV] = true
                end #if
              end #if

            end #for
          end #if
        end #for

        for v in 1:length(newVertices)
          if newVertices[v] == 1
            partition[nClusters][v] = true
            vertexToCluster[v] = nClusters
            nAddedToClusters += 1
            nNodesInCluster += 1
            visited[v] = true
          end #if
        end #for

      end #if/else

      stop = true
      for i in newVertices
        if i
          stop = false
          break
        end
      end
      if stop
        break
      end

    end #while
  end #while

  return partition, vertexToCluster, starts#, treeInds

end #partitionMatrix





function shortestPathsForCluster(mat, cluster, start)
  nVertices = mat.n
  rows = Int64[]
  columns = Int64[]
  edgeWeights = Float64[]

  for vInd in 1:length(cluster)
    if cluster[vInd]
      for eInd in mat.colptr[vInd]:(mat.colptr[vInd+1]-1) #eInd is the edge
        otherV = mat.rowval[eInd] #this gives us the other vertex
        if cluster[otherV]
          push!(rows, vInd)
          push!(columns, otherV)
          push!(edgeWeights, mat.nzval[eInd])
        end #if
      end #for
    end #for
  end #for

  for eInd in 1:length(edgeWeights)
    edgeWeights[eInd] \= 1.0
  end #for

  clusterMat = sparse(rows, columns, edgeWeights, nVertices, nVertices, min)
  spReturn = shortestPaths(clusterMat, start)
  pArray = spReturn[2]

  return pArray
end #shortestPathsForCluster



function collapsePartition(mat::SparseMatrixCSC, partition, map)
  nVertices = mat.n
  nClusters = length(partition)
  nEdges = nnz(mat)
  iVector = Int64[]
  jVector = Int64[]
  vVector = Float64[]


  for c in 1:nClusters
    for v in 1:nVertices #all vertices in cluster i
      if partition[c][v] == true
        for eInd in mat.colptr[v]:(mat.colptr[v+1]-1) #all edges connecting to vertex j
          otherV = mat.rowval[eInd]

          #if the two ends of the edge are in the same cluster don't add edge
          if partition[c][otherV] == 1
            continue

          #otherwise... add edge to matrix
          else
            push!(iVector, c) 
            push!(jVector, map[otherV])
            push!(vVector, mat.nzval[eInd])
          end #if
        end
      end #if
    end #for
  end #for


  newMat = sparse(iVector, jVector, vVector, nClusters, nClusters, min)

  return newMat
end #collapsePartition


function getBigEdgeMapReversed(mat, newMat, bigMapD)
  nEdges = nnz(mat)
  rows, columns, edgeWeights = findnz(mat)

  newNEdges = nnz(newMat)
  newWows, newColumns, newEdgeWeights = findnz(newMat)

  bigEdgeMapReversed = zeros(Int64, newNEdges)

    for newEInd in 1:newNEdges
      minEdgeWeight = Inf
      minEdgeIndex = -1

      for eInd in 1:nEdges
        c = columns[eInd]
        r = rows[eInd]
        if newWows[newEInd] == bigMapD[r] && newColumns[newEInd] == bigMapD[c] && edgeWeights[eInd] < minEdgeWeight
          minEdgeWeight = edgeWeights[eInd]
          minEdgeIndex = eInd
        end #for
      end #for columns

      bigEdgeMapReversed[newEInd] = minEdgeIndex
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


function deNormalizeEdges(mat, oldMinEdgeWeight)
  rows, columns, edgeWeights = findnz(mat)
  newEdgeWeights = zeros(Float64, length(edgeWeights))
  newEdgeWeights = edgeWeights.* oldMinEdgeWeight
  # for i in 1:length(edgeWeights)
  #   newEdgeWeights[i] = edgeWeights[i] / minEdgeWeight
  # end #for

  newMat = sparse(rows, columns, newEdgeWeights)
  return newMat
end #deNormalizeEdges




function akpw(mat::SparseMatrixCSC)
  nVertices = mat.n
  nEdges = nnz(mat)

  oldMinEdgeWeight, mat = normalizeEdgeWeights(mat)
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

    partition, vertexToCluster, starts = partitionMatrix(newMat, bigIteration, edgeClasses, bigEdgeMapReversed)

    for i = 1:length(bigMapD)
      bigMapD[i] = vertexToCluster[bigMapD[i]]
    end #for

    #this, for each cluster, creates its own sparse matrix and runs shortest paths on that, then returns the 
    # original mapped edges corresponding to that shortest paths tree. Then it adds those indices to treeInds
    for cInd in 1:length(partition)

      pArray = shortestPathsForCluster(newMat, partition[cInd], starts[cInd])

      for vInd in 1:length(pArray)
        if pArray[vInd] == 0 #|| vInd == pArray[vInd]
          continue
        end #if
        for eInd in newMat.colptr[vInd]:(newMat.colptr[vInd+1]-1) #all edges connecting to vertex j
          otherV = newMat.rowval[eInd]
          if otherV == pArray[vInd]
            treeInds[bigEdgeMapReversed[eInd]] = true
          end #if
        end #for
      end #for 

    end #for

    if (length(partition) == 1)
      break
    end #if

    newMat = collapsePartition(newMat, partition, vertexToCluster)

    bigEdgeMapReversed = getBigEdgeMapReversed(mat, newMat, bigMapD)

    bigIteration += 1
  end #while


  finalTree = sparseMatrixFromTreeIndices(mat, treeInds)
  # return finalTree

  return deNormalizeEdges(finalTree, oldMinEdgeWeight)
end #akpw

# a = grid2(3)
# n = size(a)[1]
# (ai,aj,av) = findnz(triu(a))
# edgeWeights = [5.135, 26.1885, 514.0135, 221.832, 1029.5675, 413.3675, 155.077, 83.7005, 110.4025, 16.432, 87.8085, 529.932]
# ar = sparse(ai,aj,edgeWeights,n,n)
# ar = ar + ar';
# stretchTree = akpw(ar)
# println(stretchTree)

# a = grid2(100)
# n = size(a)[1]
# (ai,aj,av) = findnz(triu(a))
# ar = sparse(ai,aj,edgeWeights,n,n)
# ar = ar + ar';
# @time stretchTree = akpw(a)






# Optimizations:

# 1. store the clusters in a better way. This is efficient for some things, but not, say, collapsePartition.
#   -for collapsePartition, instead of going through all of the vertices for each cluster, we should just go through
#     all of the edges. This is inefficient now because 
#   -Am i already doing this with the CD map?
#   -should partition construct the map? (then map won't take O(nVert*nClusters))

# 2. instead of partition tree inds being a nVertices length matrix, it should just be the length of number of
# edges we will add. then can have the index of the original 

# eliminate bigMapD
# 4. make edgeMapReversed inside collapse function?

#add reverse v->Cluster map






#the shortest path tree will not discriminate between edge classes because it would be silly to do so
  # (if there is an edge I haven't considered yet, but it makes a better shortest path tree, I should use it)






# KNOWN ISSUES
# the x(n) condition should be boundary/internal edges, not boundary/volume -- check with Spielman




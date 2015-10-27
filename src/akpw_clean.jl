include("julia/yinsGraph/graphAlgs.jl")
include("julia/yinsGraph/graphGenerators.jl")
# include("julia/yinsGraph/treeAlgs.jl")




#partition into clusters of radii O(x(n)lnn) where 1/x(n) edges border the cluster
  # means radius <= ln(n)/#edges on border
#construct shortest path spanning tree in each cluster
  #what vertex is the base for the shortest spanning tree??
#add edges of spanning tree to output tree T
#recurse after collapsing each cluster into one vertex




# QUESTIONS:

#store cluster as array..?
#yes

# store new combined graph as a new graph though.
#better to just copy onto a new sparseMatrix..? Then we have to keep track of all of those matrices..?

#show proposal

#what's the start for the shortest path tree??







#AFTER TALKING WITH SPIELMAN:
# imagine unweighted graph!!!!!
  # start with some node
  # breathfirst search on neighbors of nodes (n, n2, n3 each is group of neighbors)
  # go until 
  # number of nodes in cluster / number of neighbors of cluster is small (<= ln(n)/r).
      # this ratio is a function of R somehow
      #their x(n) is 
      #grow ball until nodes in cluster / number of neighbors <= ln(n)/r


      # |
      # V
  #start with some arbitrary node
  #volume = sum of all the degrees of nodes in cluster
  #boundary = # edges leaving cluster
  #grow balls (using basically breadth-first search -- gives us shortest path) until "boundary" of the ball/volume <= 1/x(n).
    # use their function of x(n)? make it parameter of code to be fixed later..?
  #don't really need radius function

  # in unweighted, two big issued:
    # 1. what is X(n)?
    # 2. how to we collapse..?
    # don't do data structures and stuff yet..?


    #to compute 






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


# takes a matrix from which partition will be constructed (referred to as the parent matrix)
# a partition is a double array of nodes. each sub-array c (cluster) has n vertices (the n of the parent mat)
# and either a 1 or 0 indicating whether or not v is in cluster c
# this function also returns an array of edge indices (of the parent mat) with the edges in the tree added 
function partitionMatrix(mat::SparseMatrixCSC)
  nVertices = mat.n
  partition = Array{Float64,1}[]
  nClusters = 0
  visited = zeros(Bool, nVertices)
  nAddedToClusters = 0
  treeInds = zeros(Int64, nVertices)


  # build all the clusters
  while (nAddedToClusters != nVertices)
    nClusters += 1
    push!(partition, zeros(Bool, nVertices))

    #select arbitrary unvisited node
    start = 1
    for node in 1:nVertices
      if !visited[node]
        start = node
        break
      end #if
    end #for
    visited[start] = true

    partition[nClusters][start] = true

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

      if (convert(Float64, volume) / convert(Float64, boundary) >= x) || lastnNodesInCluster == nNodesInCluster
        break
      else 
        newVertices = zeros(Bool, nVertices)
        for v in 1:length(partition[nClusters])
          if partition[nClusters][v] == true
            for eInd in mat.colptr[v]:(mat.colptr[v+1]-1) #eInd is the edge
              otherV = mat.rowval[eInd] #this gives us the other vertex

              #doing this step before for boundary... maybe can combine?
              if !visited[otherV]
                newVertices[otherV] = true

                treeInds[mat.rowval[eInd]] = eInd #..?
              end #if
            end #for
          end #if
        end #for

        for v in 1:length(newVertices)
          if newVertices[v] == 1
            partition[nClusters][v] = true
            nAddedToClusters += 1
            nNodesInCluster += 1
            visited[v] = true
          end #if
        end #for
      end #if/else

      lastnNodesInCluster = nNodesInCluster
    end #while
  end #while

  return partition, treeInds

end #partitionMatrix



function collapsePartition(mat::SparseMatrixCSC, partition, map)
  nVertices = mat.n
  nClusters = length(partition)
  iVector = Int64[]
  jVector = Int64[]
  vVector = Int64[]
  for c in 1:nClusters
    for v in 1:nVertices #all vertices in cluster i
      if partition[c][v] == 1
        for eInd in mat.colptr[v]:(mat.colptr[v+1]-1) #all edges connecting to vertex j
          if v == 49
            # println("BOOOOP: ",eInd)
          end
          otherV = mat.rowval[eInd]

          #if the two ends of the edge are in the same cluster don't add edge
          if partition[c][otherV] == 1
            continue

          #otherwise... add edge to matrix
          else
            push!(iVector, c) 
            push!(jVector, map[otherV])
            push!(vVector, 1)

            # iVector[]
            push!(iVector, map[otherV]) 
            push!(jVector, c)
            push!(vVector, 1)
          end #if
        end
      end #if
    end #for
  end #for

  newMat = sparse(iVector, jVector, vVector, nClusters, nClusters)

  return newMat
end #collapsePartition


#start with some arbitrary node
  #volume = sum of all the degrees of nodes in cluster
  #boundary = # edges leaving cluster
  #grow balls (using basically breadth-first search -- gives us shortest path) until "boundary" of the ball/volume <= 1/x(n).
    # use their function of x(n)? make it parameter of code to be fixed later..?
    # i'm converting it to volume/boundary >= x(n)
  #don't really need radius function

  # in unweighted, two big issued:
    # 1. what is X(n)?
    # 2. how to we collapse..?
    # don't do data structures and stuff yet..?

  # To collapse, will build a while new matrix, and keep a mapping to the original matrix
    # mapping: n-degree array. each value is the cluster index (v's new compressed value)

function akpw(mat::SparseMatrixCSC)
  nVertices = mat.n
  nEdges = nnz(mat)
  treeInds = zeros(Int64, nEdges)

  newMat = mat
  bigMapC = zeros(Int64, nVertices) #assuming map does C->D, bigMapC does A->C
  bigMapD = zeros(Int64, nVertices) #assuming map does C->D, bigMapD does A->D. don't think I need this
  map = zeros(Int64, nVertices)

  first = true

  while (true)
    partition, partitionTreeInds = partitionMatrix(newMat)

    if !first
      for i = 1:length(bigMapD)
        bigMapC[i] = bigMapD[i]
      end #for
    end #if

    #making the new map
    for c = 1:length(partition)
      for v = 1:length(partition[c])
        # println("checking: ", c, ", ", v)
        if partition[c][v] == 1
          # println("printing: ", v)
          map[v] = c
        end #if
      end #for
    end #for

    if first
      for i = 1:length(bigMapD)
        bigMapD[i] = map[i]
        bigMapC[i] = map[i]
      end #for

    else
      for i = 1:length(bigMapD)
        bigMapD[i] = map[bigMapD[i]]
      end #for
    end #if


    #CAN WE SPEED THIS UP BY USING partition INSTEAD OF iterating through the maps..?
    # append every edge in partitionTree to treeInds
    for e in 1:length(partitionTreeInds)
      if partitionTreeInds[e] == 0
        continue
      end

      if first
        treeInds[partitionTreeInds[e]] = 1
      else

        a = findn(newMat)
        newr = a[1][partitionTreeInds[e]]
        newc = a[2][partitionTreeInds[e]]


        for i = 1:length(bigMapC)
          found = false

          #this finds an edge in big mat corresponding to e in the partitionTree
          if bigMapC[i] == newr
            for eInd in mat.colptr[i]:(mat.colptr[i+1]-1)
              edgeV = mat.rowval[eInd]
              if bigMapC[edgeV] == newc
                found = true
                # println("treeInd: ", eInd)
                treeInds[eInd] = 1
                break
              end #if
            end #for
          end #if

          if found
            break
          end #if
        end #for
      end #if/else


    end #for
    

    if (length(partition) == 1)
      break
    end #if

    newMat = collapsePartition(newMat, partition, map)

    first = false

  end #while

  t2 = zeros(Int64, nVertices - 1)
  t2Counter = 1

  for i in 1:length(treeInds)
    if treeInds[i] == 1
      t2[t2Counter] = i
      t2Counter += 1
    end
  end

  (ai,aj,av) = findnz(mat);
  tr2 = sparse(ai[t2],aj[t2],av[t2],nVertices,nVertices)
  tr2 = tr2 + tr2';

  return tr2

end #akpw


a = grid2(125)
n = size(a)[1]
(ai,aj,av) = findnz(triu(a))

@time stretchTree = akpw(a)
# primTree = prim(a)
# println("done")
# println(stretchTree)
# println(compStretches(stretchTree, a))
# println(compStretches(stretchTree, a))




# BUGS:

# why do clusters have weird edge weights
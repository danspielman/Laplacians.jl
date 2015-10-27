

#using PyCall
#VERSION

#include("julia/yinsGraph/graphAlgs.jl")
#include("julia/yinsGraph/graphGenerators.jl")

#x = linspace(0,2*pi,1000); y = sin(3*x + 4*cos(2*x));


#a = grid2(1000)
#n = size(a)[1]
#(ai,aj,av) = findnz(triu(a))
#ar = sparse(ai,aj,rand(size(av)),n,n)
#ar = ar + ar'


import yinsGraph.intHeap
import yinsGraph: intHeapAdd!, intHeapPop!, intHeapSet!




function prim(mat::SparseMatrixCSC)

#  tic()

  nVertices = mat.n
  nh = intHeap(nVertices)
  #   visited = empty array
  visited = zeros(Bool, nVertices)
  associatedEdges = zeros(Int64, nVertices) #this is the edge to be associated with each vertex eventually. 0 if unused, -1 if visited

  #   tree = [], gonna by filled with edge-indices
  treeInds = zeros(Int64, nVertices)

#   visit some arbitrary node
  visited[1] = true
  associatedEdges[1] = -1
#   add all the vertices (except first one) to the heap, their values all infinite
  for vInd in 2:nVertices
    intHeapAdd!(nh, vInd, Inf)
  end #for


#  toc()
#  tic()

#   set all vertices to heap that border the first vertex
  for eInd in mat.colptr[1]:(mat.colptr[2]-1)
    #for each edge connected to first vertex, reset its other vertex to the value of the edge between them
    intHeapSet!(nh, mat.rowval[eInd], mat.nzval[eInd])
    associatedEdges[mat.rowval[eInd]] = eInd #this won't work if there are multiple edges between two vertices... is that possible?
  end #for

#   UNTIL visited IS FULL (has as many nodes as mat = mat.n) -- lets assume mat is connected
#   aka add an edge n-1 times
  for x in 1:nVertices-1
#     sort the vertices by edges connected to T
#     pick smallest one such the the other node isn't also in T
    vInd = intHeapPop!(nh)
    while visited[vInd]
      vInd = intHeapPop!(nh) #pop returns the heap value, not the vertex index..?
    end #while
    visited[vInd] = true
    treeInds[vInd] = associatedEdges[vInd]
    associatedEdges[vInd] = -1

    for eInd in mat.colptr[vInd]:(mat.colptr[vInd+1]-1) #we've just added a vertex, now going through vertices next to this one
      edgeWeight = mat.nzval[eInd] #the edge value for each neighboring vertex
      otherVert = mat.rowval[eInd]
      previousEdgeIndex = associatedEdges[otherVert]
      if previousEdgeIndex == -1
        continue
      end #if
#       print(previousEdgeIndex)
      if (previousEdgeIndex == 0 || edgeWeight < mat.nzval[previousEdgeIndex])
        #if the edge weight is less than the vertex's previously associated edge, or edge hasn't been touching a tree-v before

        intHeapSet!(nh, otherVert, edgeWeight)
        associatedEdges[otherVert] = eInd #then reset it's associated edge
      end #if
    end # for
  end #for

#  toc()
#  tic()
  
#   finalTree =

  t2 = treeInds[2:nVertices];
    
  (ai,aj,av) = findnz(mat);
  tr2 = sparse(ai[t2],aj[t2],av[t2],nVertices,nVertices)
  tr2 = tr2 + tr2';
    
#  toc()

  return tr2


end #prim


function prim2{Tval,Tind}(mat::SparseMatrixCSC{Tval,Tind})

  tic()

  nVertices = mat.n
  nh = intHeap(nVertices)
  #   visited = empty array
  visited = zeros(Bool, nVertices)
  associatedEdges = zeros(Tind, nVertices) #this is the edge to be associated with each vertex eventually. 0 if unused, -1 if visited

  #   tree = [], gonna by filled with edge-indices
  treeInds = zeros(Tind, nVertices)

#   visit some arbitrary node
  visited[1] = true
  associatedEdges[1] = -1
#   add all the vertices (except first one) to the heap, their values all infinite
  for vInd in 2:nVertices
    intHeapAdd!(nh, vInd, Inf)
  end #for


  toc()
  tic()

#   set all vertices to heap that border the first vertex
  for eInd in mat.colptr[1]:(mat.colptr[2]-1)
    #for each edge connected to first vertex, reset its other vertex to the value of the edge between them
    intHeapSet!(nh, mat.rowval[eInd], mat.nzval[eInd])
    associatedEdges[mat.rowval[eInd]] = eInd #this won't work if there are multiple edges between two vertices... is that possible?
  end #for

#   UNTIL visited IS FULL (has as many nodes as mat = mat.n) -- lets assume mat is connected
#   aka add an edge n-1 times
  for x in 1:nVertices-1
#     sort the vertices by edges connected to T
#     pick smallest one such the the other node isn't also in T
    vInd = intHeapPop!(nh)
    while visited[vInd]
      vInd = intHeapPop!(nh) #pop returns the heap value, not the vertex index..?
    end #while
    visited[vInd] = true
    treeInds[vInd] = associatedEdges[vInd]
    associatedEdges[vInd] = -1

    for eInd in mat.colptr[vInd]:(mat.colptr[vInd+1]-1) #we've just added a vertex, now going through vertices next to this one
      edgeWeight = mat.nzval[eInd] #the edge value for each neighboring vertex
      otherVert = mat.rowval[eInd]
      previousEdgeIndex = associatedEdges[otherVert]
      if previousEdgeIndex == -1
        continue
      end #if
#       print(previousEdgeIndex)
      if (previousEdgeIndex == 0 || edgeWeight < mat.nzval[previousEdgeIndex])
        #if the edge weight is less than the vertex's previously associated edge, or edge hasn't been touching a tree-v before

        intHeapSet!(nh, otherVert, edgeWeight)
        associatedEdges[otherVert] = eInd #then reset it's associated edge
      end #if
    end # for
  end #for

  toc()
  tic()
  
#   finalTree =

  t2 = treeInds[2:end];
    
  (ai,aj,av) = findnz(a);
  tr2 = sparse(ai[t2],aj[t2],av[t2],nVertices,nVertices)
  tr2 = tr2 + tr2';
    
  toc()

  return tr2


end #prim

#@time tree = prim(ar)


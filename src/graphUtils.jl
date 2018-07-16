# using these is a little slower than the way we do it in the code, but ok
# can say:
# for i in 1:deg(mat, v)
#   nbr = nbri(mat, v, i)
#   wt = weighti(mat, v, i)
#
# but, it is faster to write:
#  for ind in colptr[curNode]:(colptr[curNode+1]-1)
#    nbr = rowval[ind]

deg(mat::SparseMatrixCSC{Tv,Ti}, v::Ti) where {Tv,Ti} = mat.colptr[v+1]-mat.colptr[v]
nbri(mat::SparseMatrixCSC{Tv,Ti}, v::Ti, i::Ti) where {Tv,Ti} = mat.rowval[mat.colptr[v]+i-1]
weighti(mat::SparseMatrixCSC{Tv,Ti}, v::Ti, i::Ti) where {Tv,Ti} = mat.nzval[mat.colptr[v]+i-1]
nbrs(mat::SparseMatrixCSC{Tv,Ti}, v::Ti) where {Tv,Ti} = mat.rowval[mat.colptr[v]:(mat.colptr[v+1]-1)]


" Finds the weighted degree of a vertex in the graph "
function wdeg(mat::SparseMatrixCSC{Tv,Ti}, v::Ti) where {Tv,Ti}
  sum = 0
  for i in 1:deg(mat,v)
    sum = sum + weighti(mat,v,i)
  end
  return sum
end


" Sets the value of a certain edge in a sparse graph; value can be 0 without the edges dissapearing "
function setValue(mat::SparseMatrixCSC{Tv,Ti}, v::Ti, i::Ti, a::Tv) where {Tv,Ti}
  mat.nzval[mat.colptr[v]+i-1] = a
end # setValue


""" For a symmetric matrix, this gives the correspondance between pairs of entries in an ijv.
So, ai[ind] = aj[flip[ind]].  For example, 

~~~
(ai,aj,av) = findnz(a);
fl = flipIndex(a)
ind = 10
@show backind = fl[10]
@show [ai[ind], aj[ind], av[ind]]
@show [ai[backind], aj[backind], av[backind]];

backind = fl[10] = 4
[ai[ind],aj[ind],av[ind]] = [2.0,4.0,0.7]
[ai[backind],aj[backind],av[backind]] = [4.0,2.0,0.7]
~~~
"""
function flipIndex(a::SparseMatrixCSC{Tval,Tind}) where {Tval,Tind}

  b = SparseMatrixCSC(a.m, a.n, copy(a.colptr), copy(a.rowval), collect(UnitRange{Tind}(1,nnz(a))) );
  bakMat = sparse(b');
  return bakMat.nzval

end


" Computes the back indices in a graph in O(M+N). works if for every edge (u,v), (v,u) is also in the graph "
function backIndices(G::SparseMatrixCSC{Tv,Ti}) where {Tv,Ti}
  n = max(G.n, G.m)

  # initialize values
  u,v,w = findEntries(G)

  backInfo = sparse(u, v, [(0,0) for i in 1:length(u)])
  backIndices = sparse(u, v, [-1 for i in 1:length(u)])
  index = zeros(Ti, n)

  # compute for each vertex the incoming edges
  for u in 1:n
    for i in 1:deg(G, u)
      v = nbri(G, u, i)
      index[v] = index[v] + 1
      setValue(backInfo, v, index[v], (u, i))
    end
  end

  # create the back edges
  incoming = zeros(Ti, n)
  for u in 1:n
    for i in 1:deg(G, u)
      incoming[weighti(backInfo, u, i)[1]] = weighti(backInfo, u, i)[2]
    end

    for i in 1:deg(G, u)
      v = nbri(G, u, i)
      setValue(backIndices, v, incoming[v], i)
    end
  end

  return backIndices

end # backIndices

" Same as the above, but now the graph is in adjacency list form "
function backIndices(G::Array{Array{Tuple{Tv1,Tv2},1},1}) where {Tv1,Tv2}
  n = length(G)

  # initialize values
  backInfo = [Tuple{Int64,Int64}[] for i in 1:n]
  backIndices = [(Int64)[] for i in 1:n]
  for u in 1:length(G)
    listU = G[u]
    backInfo[u] = [(0,0) for i in 1:length(listU)]
    backIndices[u] = [-1 for i in 1:length(listU)]
  end
  index = zeros(Int64, n)

  # compute for each vertex the incoming edges
  for u in 1:n
    for i in 1:length(G[u])
      v,w = G[u][i]
      index[v] = index[v] + 1
      backInfo[v][index[v]] = (u, i)
    end
  end

  # create the back edges
  incoming = zeros(Int64, n)
  for u in 1:n
    for i in 1:length(G[u])
      incoming[backInfo[u][i][1]] = backInfo[u][i][2]
    end

    for i in 1:length(G[u])
      v,w = G[u][i]
      backIndices[v][incoming[v]] = i
    end
  end

  return backIndices

end # backIndices


" Similar to findnz, but also returns 0 entries that have an edge in the sparse matrix "
function findEntries(G::SparseMatrixCSC{Tv,Ti}) where {Tv,Ti}

  n = max(G.n, G.m)

  u = copy(G.rowval)
  v = zeros(Ti, length(u))
  for i in 1:n
    for j in (G.colptr[i]):(G.colptr[i + 1] - 1)
      v[j] = i
    end
  end
  w = copy(G.nzval)

  return u,v,w

end # findEntries


""" 
  Returns the quality of the cut for a given graph and a given cut set s.
  the result will be |outgoing edges| / min(|vertices in set|, |N - vertices in set|)
"""
function compConductance(G::SparseMatrixCSC{Tv,Ti}, s::Array{Int64,1}) where {Tv, Ti}

  n = max(G.n, G.m)

  vols = getVolume(G, s)
  obound = getObound(G, s)

  return obound / vols

end # compConductance


" Computes the volume of subset s in an unweighted graph G "
function getVolume(G::SparseMatrixCSC{Tv,Ti}, s::Array{Int64,1}) where {Tv,Ti}

  vol = 0.0
  for i in 1:length(s)
    deggsi = 0
    for j in 1:deg(G, s[i])
      deggsi = deggsi + weighti(G,s[i],j)
    end
    # vol = vol + deg(G, s[i])
    vol = vol + deggsi
  end

  return vol

end # getVolume


" Computes the number of edges leaving s "
function getObound(G::SparseMatrixCSC{Tv,Ti}, s::Array{Int64,1}) where {Tv,Ti}

  sets = IntSet(s)

  obound = 0
  for u in s
    for i in 1:deg(G, u)
      v = nbri(G, u, i)
      if !(v in sets)
        obound = obound + weighti(G, u, i)
      end
    end
  end

  return obound

end #getObound
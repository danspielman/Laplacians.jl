# using these is a little slower than the way we do it in the code, but ok
# can say:
# for i in 1:deg(mat, v)
#   nbr = nbri(mat, v, i)
#   wt = weighti(mat, v, i)
#
# but, it is faster to write:
#  for ind in colptr[curNode]:(colptr[curNode+1]-1)
#    nbr = rowval[ind]

deg{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti}, v::Ti) = mat.colptr[v+1]-mat.colptr[v]
nbri{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti}, v::Ti, i::Ti) = mat.rowval[mat.colptr[v]+i-1]
weighti{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti}, v::Ti, i::Ti) = mat.nzval[mat.colptr[v]+i-1]
nbrs{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti}, v::Ti) = mat.rowval[mat.colptr[v]:(mat.colptr[v+1]-1)]

function setValue{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti}, v::Ti, i::Ti, a::Tv)
  mat.nzval[mat.colptr[v]+i-1] = a
end

" computes the back indices in a graph in O(M+N) "
function backIndices{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti})
  n = max(G.n, G.m)

  # initialize values
  u,v,w = findnz(G)
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

end

# This is code for generating random trees of various sorts.
# The hope is that one might have low stretch,
# although none is guaranteed to.

using DataStructures

# sample edges with probability proportional to their weight

function randishKruskal{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti})
    n = size(mat)[1]
    (ai,aj,av) = findnz(triu(mat))
    m = length(av)
    
    s = Sampler(av)

    comps = IntDisjointSets(n)

    treeinds = zeros(Int64,n-1)
    numintree = 0
    for j in 1:m
        i = pop!(s)
        if !DataStructures.in_same_set(comps,ai[i],aj[i])
            numintree = numintree+1
            treeinds[numintree] = i
            DataStructures.union!(comps,ai[i],aj[i])
        end
    end

    tree = sparse(ai[treeinds],aj[treeinds],av[treeinds],n,n)
    tree = tree + tree'

    return tree

end

function randishPrim{Tval,Tind}(mat::SparseMatrixCSC{Tval,Tind})

  n = mat.n
  m = nnz(mat)

  (ai, aj, av) = findnz(mat)
  flipInd = flipIndex(mat)

    
  visited = zeros(Bool,n)

  s = Sampler(m)

  treeEdges = zeros(Tind,n-1)
  numEdges = 0    

  # should randomize the initial choice
  v = one(Tind)

  # always add edge by highest vertex number first
  for ind in mat.colptr[v]:(mat.colptr[v+1]-1)
      wt = mat.nzval[ind]
      push!(s, ind, wt)
  end

  visited[v] = true

  while s.nitems > 0

    edge = pop!(s)

    v = aj[edge]
    if visited[v]
        v = ai[edge]
    end

    if !visited[v]

        numEdges += 1
        treeEdges[numEdges] = edge

        
        for ind in mat.colptr[v]:(mat.colptr[v+1]-1)
            u = mat.rowval[ind]

            wt = mat.nzval[ind]
            useind = v < u ? ind : flipInd[ind]
            if !visited[u]
                push!(s, useind, wt)

            end
        end

        visited[v] = true
    end

  end

  treeEdges = treeEdges[1:numEdges]

  tr = sparse(ai[treeEdges],aj[treeEdges],av[treeEdges],n,n)
  tr = tr + tr';

  return tr

end # randishPrim


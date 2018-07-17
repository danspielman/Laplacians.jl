
# Rooted tree
#   parent of root is self
#   the kids of node i are children[kidsPtr[i] : kidsPtr[i]-numKids[i]-1 ]
#  the weigts of edges from parents to kids are indexed similarly.
#  that is, weights holds the weight of an edge to the parent.
#
#  We require that the children appears in a dfs order
#

struct RootedTree{Tval,Tind}
  root::Tind
  parent::Array{Tind,1}
  children::Array{Tind,1}
  weights::Array{Tval,1}
  numKids::Array{Tind,1}
  kidsPtr::Array{Tind,1}
end # RootedTree

# the following two routines will return the children of a node,
# and the weights on the edges to them.
# but, using these can be slow as it seems to cause memory allocation

kids(tr::RootedTree, i) = tr.children[tr.kidsPtr[i]:(tr.kidsPtr[i]+tr.numKids[i] - 1)]
weights(tr::RootedTree, i) = tr.weights[tr.kidsPtr[i]:(tr.kidsPtr[i]+tr.numKids[i] - 1)]



# matToTree takes a weighted tree in mat format, and produces
# a RootedTree from it

matToTree(mat::SparseMatrixCSC{Tv,Ti}) where {Tv,Ti} = matToTree(mat, 1::Ti)

function matToTree(mat::SparseMatrixCSC{Tv,Ti}, root::Ti) where {Tv,Ti}

  n = mat.n
  visited = zeros(Bool,n)

  children = zeros(Ti,n)
  weights = zeros(Tv,n)
  parent = zeros(Ti,n)
  numKids = zeros(Ti,n)
  kidsPtr = zeros(Ti,n)
    
  children[1] = root
  visited[root] = true
  parent[root] = root
  kidsPtr[root] = 2
    
  degRoot = deg(mat,root)
  numKids[root] = degRoot

  children[1 .+ collect(1:degRoot)] = mat.rowval[mat.colptr[root]:(mat.colptr[root+1]-1)]
  weights[1 .+ collect(1:degRoot)] = mat.nzval[mat.colptr[root]:(mat.colptr[root+1]-1)]
   
  numIn = 1+degRoot

  parent[mat.rowval[mat.colptr[root]:(mat.colptr[root+1]-1)]] .= root
  visited[mat.rowval[mat.colptr[root]:(mat.colptr[root+1]-1)]] .= true

  ptr = 2
    
  while ptr <= numIn
    v = children[ptr]
    numKids[v] = deg(mat,v)-1

    if numKids[v] > 0
        kidsPtr[v] = numIn+1
    end

    for ind in mat.colptr[v]:(mat.colptr[v+1]-1)
      nbr = mat.rowval[ind]
      if !visited[nbr]
          numIn = numIn+1
          visited[nbr] = true
          children[numIn] = nbr
          weights[numIn] = mat.nzval[ind]

          parent[nbr] = v
      end # if
    end # for
 
    ptr = ptr+1

  end # while

  if numIn < n
      error("graph is not connected")
  end
    
  return RootedTree(root, parent, children, weights, numKids, kidsPtr)

end

# create a rooted tree structure, and compute the depth as we go
# this is faster than doing it independently

matToTreeDepth(mat::SparseMatrixCSC{Tv,Ti}) where {Tv,Ti} = matToTreeDepth(mat, 1::Ti)

function matToTreeDepth(mat::SparseMatrixCSC{Tv,Ti}, root::Ti) where {Tv,Ti}

  n = mat.n
  visited = zeros(Bool,n)

  children = zeros(Ti,n)
  weights = zeros(Tv,n)
  parent = zeros(Ti,n)
  numKids = zeros(Ti,n)
  kidsPtr = zeros(Ti,n)
    
  depth = zeros(Tv,n)

  children[1] = root
  visited[root] = true
  parent[root] = root
  kidsPtr[root] = 2
    
  degRoot = deg(mat,root)
  numKids[root] = degRoot

  kids = mat.rowval[mat.colptr[root]:(mat.colptr[root+1]-1)]
  wts =  mat.nzval[mat.colptr[root]:(mat.colptr[root+1]-1)]    

  children[1 .+ collect(1:degRoot)] .= kids
  weights[1 .+ collect(1:degRoot)] .= wts
   
  numIn = 1+degRoot

  parent[kids] .= root
  visited[kids] .= true
  depth[kids] = 1 ./ wts

  ptr = 2
    
  while ptr <= numIn
    v = children[ptr]
    numKids[v] = deg(mat,v)-1

    if numKids[v] > 0
        kidsPtr[v] = numIn+1
    end

    for ind in mat.colptr[v]:(mat.colptr[v+1]-1)
      nbr = mat.rowval[ind]
      if !visited[nbr]
          numIn = numIn+1
          visited[nbr] = true
          children[numIn] = nbr
          weights[numIn] = mat.nzval[ind]

          parent[nbr] = v
          depth[nbr] = depth[v] + 1/mat.nzval[ind]
          
      end # if
    end # for
 
    ptr = ptr+1

  end # while

  if numIn < n
      error("graph is not connected")
  end
    
  return RootedTree(root, parent, children, weights, numKids, kidsPtr), depth

end


import Base.pop!

# put the nodes of a tree in a dfs order,
# so the last node we see is a leaf,
# the first is the root,
# and each subtree is completed before a parallel subtree
#
# we need this for the non-recursive tarjanLCA
function dfsOrder(tr::RootedTree)
    ord = Int64[]
    stk = Int64[]
    push!(stk,tr.root)
    while ~isempty(stk)
        u = pop!(stk)
        push!(ord,u)
        for  vi in tr.kidsPtr[u]:(tr.numKids[u] + tr.kidsPtr[u] - 1)
            v = tr.children[vi]
            push!(stk,v)
        end
    end
    return ord

end


# put the nodes of a tree in a dfs order,
# so the last node we see is a leaf,
# the first is the root,
# and each subtree is completed before a parallel subtree
#
# this matrix should be the adj matrix of a tree
function dfsOrder(t::SparseMatrixCSC{Tv,Ti}; start::Ti = 1) where {Tv,Ti}

    n = size(t,1)
    seen = zeros(Bool,n)
    ord = Int64[]
    stk = Int64[]
    push!(stk,start)
    while ~isempty(stk)
        u = pop!(stk)
        push!(ord,u)
        seen[u] = true
        for vi in t.colptr[u]:(t.colptr[u+1]-1)
            v = t.rowval[vi]
            if !seen[v]
                push!(stk,v)
            end
        end
    end
    return ord

end


# put the nodes of a tree in a dfs order,
# so the first node we see is a leaf,
# the last is the root,
# and each subtree is completed before a parallel subtree
#
# we need this for the non-recursive tarjanLCA
#=
function dfsOrder(t::RootedTree)

    n = size(t.children,1)
    
    sz = ones(Int64,n);
    for vi in n:-1:2
        v = t.children[vi]
        par = t.parent[v]
        sz[par] += sz[v]
    end

    numLeft = zeros(Int64,n)
    for ui in 1:n
        u = t.children[ui]
        cnt = numLeft[u]
        for vi in t.kidsPtr[u]:(t.numKids[u] + t.kidsPtr[u] - 1)
            v = t.children[vi]
            numLeft[v] = cnt
            cnt += sz[v]
        end
    end
    
    ordInd = numLeft + sz

    ord = zeros(Int64,n)
    for i in 1:n
        ord[ordInd[i]] = i
    end

    return ord
    
end
=#

# this is intended to be used with a tree
function bfsOrder(mat::SparseMatrixCSC{Tv,Ti}, start::Ti) where {Tv,Ti}
  n = mat.n
  visited = zeros(Bool,n)

  order = zeros(Ti,n)
  order[1] = start
  ptr::Ti = 1
  numInOrder::Ti = 1
  visited[start] = true
    
  while ptr <= numInOrder
    v = order[ptr]

    for ind in mat.colptr[v]:(mat.colptr[v+1]-1)
      nbr = mat.rowval[ind]
      if !visited[nbr]
          numInOrder = numInOrder+1
          order[numInOrder] = nbr
          visited[nbr] = true
      end # if
    end # for
 
    ptr = ptr+1

  end # while

  if numInOrder < n
      error("graph is not connected")
  end
    
  return order

end # bfsOrder




#-----------------------------
#
# code for computing stretch



# tarjan stretch uses tarjans offline lca algorithm to compute stretches

# it seems that this cannot be made much faster:
# almost all of the time is in the disjoint sets code.
# and, when I just run that without the rest, it is most of the time.
# this will save time by computing stretch directly,
# by using the depth
function tarjanStretch(t::RootedTree{Tv,Ti}, mat::SparseMatrixCSC{Tv,Ti}, depth::Array{Tv,1}) where {Tv,Ti}
    n = length(t.parent)
    su = IntDisjointSets(n)

    ancestor = collect(1:n)

    answer = zeros(Tv,nnz(mat))

    seen = zeros(Bool, n)

    tarjanStretchSub(t.root, t, mat, ancestor, answer, seen, su, depth)

    stretches = copy(mat)
    for i in 1:length(stretches.nzval)
        stretches.nzval[i] = answer[i]
    end
    stretches = stretches + stretches'
    
    return stretches

end # tarjanStretch

function tarjanStretchSub(u::Ti, t::RootedTree{Tv,Ti}, mat::SparseMatrixCSC{Tv,Ti}, ancestor::Array{Ti,1},
                          answer::Array{Tv,1}, seen::Array{Bool,1}, su::IntDisjointSets, depth::Array{Tv,1}) where {Tv,Ti}    

    ord = dfsOrder(t)

    # traverse nodes from leaves back to root
    for vi in size(ord,1):-1:1
        v = ord[vi]

        par = t.parent[v]

        # just for debugging
        if seen[par]
            error("saw parent!")
        end

        seen[v] = true

        # ancestor[DataStructures.find_root(su, v)] = par
        
        for ind in mat.colptr[v]:(mat.colptr[v+1]-1)
            w = mat.rowval[ind]
            if seen[w]
                answer[ind] = mat.nzval[ind]*(depth[v] + depth[w] - 2*depth[ancestor[DataStructures.find_root(su,w)]])
                # println(u, " ", v, " : ", answer[ind])
            end # can fill u-v query
        end # over queries

        DataStructures.union!(su, par, v)
        
        ancestor[DataStructures.find_root(su, par)] = par
        

     end # for vi, v
        
end # TarjanStretchSub



"""Compute the stretched of every edge in `mat` with respect to the tree `tree`.
Returns the answer as a sparse matrix with the same nonzero structure as `mat`.
Assumes that `mat` is symmetric.
`tree` should be the adjacency matrix of a spanning tree."""
function compStretches(tree::SparseMatrixCSC{Tv,Ti}, mat::SparseMatrixCSC{Tv,Ti}) where {Tv,Ti}

    t, depth = matToTreeDepth(tree)

    stretches = tarjanStretch(t,mat,depth)
    return stretches
    
end # compStretches


"""Compute the vector of depths in a tree that is in DFS order,
*with the root at the first position, and the leaves at the end*
"""
function treeDepthDFS(tree::SparseMatrixCSC{Tv,Ti}) where {Tv,Ti}
    n = tree.n

    depth = zeros(Tv,n)
    root = 1
    kids = tree.rowval[tree.colptr[root]:(tree.colptr[root+1]-1)]
    wts =  tree.nzval[tree.colptr[root]:(tree.colptr[root+1]-1)]    
    depth[kids] = 1 ./wts

    for v in 2:n

        for ind in (tree.colptr[v]+1):(tree.colptr[v+1]-1)
            kid = tree.rowval[ind]
            depth[kid] = depth[v] + 1/tree.nzval[ind]
        end
    end

    return depth
end







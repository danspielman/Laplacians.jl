
# Rooted tree
#   parent of root is self
#   the kids of node i are children[kidsPtr[i] : kidsPtr[i]-numKids[i]-1 ]
#  the weigts of edges from parents to kids are indexed similarly.
#  that is, weights holds the weight of an edge to the parent.
#
#  We require that the children appears in a dfs order
#

type RootedTree{Tval,Tind}
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

matToTree{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti}) = matToTree(mat, 1::Ti)

function matToTree{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti}, root::Ti)

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

  children[1 + collect(1:degRoot)] = mat.rowval[mat.colptr[root]:(mat.colptr[root+1]-1)]
  weights[1 + collect(1:degRoot)] = mat.nzval[mat.colptr[root]:(mat.colptr[root+1]-1)]
   
  numIn = 1+degRoot

  parent[mat.rowval[mat.colptr[root]:(mat.colptr[root+1]-1)]] = root
  visited[mat.rowval[mat.colptr[root]:(mat.colptr[root+1]-1)]] = true

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

matToTreeDepth{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti}) = matToTreeDepth(mat, 1::Ti)

function matToTreeDepth{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti}, root::Ti)

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

  children[1 + collect(1:degRoot)] = kids
  weights[1 + collect(1:degRoot)] = wts
   
  numIn = 1+degRoot

  parent[kids] = root
  visited[kids] = true
  depth[kids] = 1./wts

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
    ord = Array{Int64}(0)
    stk = Array{Int64}(0)
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
function dfsOrder{Tv,Ti}(t::SparseMatrixCSC{Tv,Ti}; start::Ti = 1)

    n = size(t,1)
    seen = zeros(Bool,n)
    ord = Array{Int64}(0)
    stk = Array{Int64}(0)
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
function bfsOrder{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti}, start::Ti)
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



# compute the depth of every node in a rooted tree,
# taking the weight as the reciprocal of distance
function compDepth{Tv,Ti}(t::RootedTree{Tv,Ti})
    n = length(t.parent)
    depth = zeros(Tv,n)

    depth[t.root] = 0

    for v in t.children
        for i in t.kidsPtr[v]:(t.numKids[v] + t.kidsPtr[v] - 1)
            depth[t.children[i]] = depth[v] + 1/t.weights[i]
        end
    end

    return depth
end # compDepth


# tarjan stretch uses tarjans offline lca algorithm to compute stretches

# it seems that this cannot be made much faster:
# almost all of the time is in the disjoint sets code.
# and, when I just run that without the rest, it is most of the time.
# this will save time by computing stretch directly,
# by using the depth
function tarjanStretch{Tv,Ti}(t::RootedTree{Tv,Ti}, mat::SparseMatrixCSC{Tv,Ti}, depth::Array{Tv,1})
    n = length(t.parent)
    su = IntDisjointSets(n)

    ancestor = collect(1:n)

    answer = zeros(Tv,nnz(mat))

    seen = zeros(Bool, n)

    tarjanStretchSub(t.root, t, mat, ancestor, answer, seen, su, depth)

    stretches = copy(mat)
    stretches.nzval = answer
    stretches = stretches + stretches'
    
    return stretches

end # tarjanStretch

function tarjanStretchSub{Tv,Ti}(u::Ti, t::RootedTree{Tv,Ti}, mat::SparseMatrixCSC{Tv,Ti}, ancestor::Array{Ti,1},
                                 answer::Array{Tv,1}, seen::Array{Bool,1}, su::IntDisjointSets, depth::Array{Tv,1})    

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



function compStretches{Tv,Ti}(t::RootedTree{Tv,Ti}, mat::SparseMatrixCSC{Tv,Ti})
    n = length(t.order)
    
    depth = compDepth(t)

    stretches = tarjanStretch(t,mat,depth)
    return stretches

end # compStretches

"""Compute the stretched of every edge in `mat` with respect to the tree `tree`.
Returns the answer as a sparse matrix with the same nonzero structure as `mat`.
Assumes that `mat` is symmetric.
`tree` should be the adjacency matrix of a spanning tree."""
function compStretches{Tv,Ti}(tree::SparseMatrixCSC{Tv,Ti}, mat::SparseMatrixCSC{Tv,Ti})

    t, depth = matToTreeDepth(tree)

    stretches = tarjanStretch(t,mat,depth)
    return stretches
    
end # compStretches


"""Compute the vector of depths in a tree that is in DFS order,
*with the root at the first position, and the leaves at the end*
"""
function treeDepthDFS{Tv,Ti}(tree::SparseMatrixCSC{Tv,Ti})
    n = tree.n

    depth = zeros(Tv,n)
    root = 1
    kids = tree.rowval[tree.colptr[root]:(tree.colptr[root+1]-1)]
    wts =  tree.nzval[tree.colptr[root]:(tree.colptr[root+1]-1)]    
    depth[kids] = 1./wts

    for v in 2:n

        for ind in (tree.colptr[v]+1):(tree.colptr[v+1]-1)
            kid = tree.rowval[ind]
            depth[kid] = depth[v] + 1/tree.nzval[ind]
        end
    end

    return depth
end


"""Compute the stretched of every edge in `mat` with respect to the tree `tree`.
Returns the answer as a sparse matrix with the same nonzero structure as `mat`.
Assumes that `mat` is symmetric.
`tree` should be the adjacency matrix of a spanning tree, 
*ordered by DFS so that every parent comes before its children in the order*"""
function compStretchesDFS{Tv,Ti}(tree::SparseMatrixCSC{Tv,Ti}, mat::SparseMatrixCSC{Tv,Ti})


    n = size(tree,1)
    su = IntDisjointSets(n)

    ancestor = collect(1:n)
    answer = zeros(Tv,nnz(mat))
    seen = zeros(Bool, n)

    depth = treeDepthDFS(tree)
    
    # traverse nodes from leaves back to root
    for v in n:-1:1

        if v > 1
            par = tree.rowval[tree.colptr[v]]
        else
            par = 1
        end

        # just for debugging
        if seen[par]
            error("saw parent!")
        end

        seen[v] = true

        for ind in mat.colptr[v]:(mat.colptr[v+1]-1)
            w = mat.rowval[ind]
            if seen[w]
                answer[ind] = mat.nzval[ind]*(depth[v] + depth[w] - 2*depth[ancestor[DataStructures.find_root(su,w)]])
            end # can fill u-v query
        end # over queries

        DataStructures.union!(su, par, v)
        
        ancestor[DataStructures.find_root(su, par)] = par
        

    end # for v

    stretches = copy(mat)
    stretches.nzval = answer
    stretches = stretches + stretches'

    return stretches
    
end # compStretchesDFS




#-----------------------------------
#
# extra stuff we don't need

# this is tarjans offline lca algorithm
# we don't use it.
#
# it seems that this cannot be made much faster:
# almost all of the time is in the disjoint sets code.
# and, when I just run that without the rest, it is most of the time.
function tarjanOLCA{Tv,Ti}(t::RootedTree{Tv,Ti}, mat::SparseMatrixCSC{Tv,Ti})
    n = length(t.order)
#    su = SetUnion(n)
    su = IntDisjointSets(n)

    ancestor = zeros(Ti,n)

    answer = zeros(Ti,nnz(mat))

    seen = zeros(Bool, n)

    tarjanSub(t.root, t, mat, ancestor, answer, seen, su)

    ansmat = copy(mat)
    ansmat.nzval = answer
    ansmat = ansmat + ansmat'
    
    return ansmat

end # tarjanOLCA

function tarjanSub{Tv,Ti}(u::Ti, t::RootedTree{Tv,Ti}, mat::SparseMatrixCSC{Tv,Ti}, ancestor::Array{Ti,1}, answer::Array{Ti,1}, seen::Array{Bool,1}, su::IntDisjointSets)    

        ancestor[u] = u

        for v in t.children[u]
            tarjanSub(v, t, mat, ancestor, answer, seen, su)
            DataStructures.union!(su, u, v)
            ancestor[DataStructures.find_root(su, u)] = u
#            ancestor[DataStructures.find!(su, u)] = u
        end

        seen[u] = true

        for ind in mat.colptr[u]:(mat.colptr[u+1]-1)
            v = mat.rowval[ind]
            if seen[v]
#                answer[ind] = ancestor[find!(su,v)]
                answer[ind] = ancestor[DataStructures.find_root(su,v)]
                # println(u, " ", v, " : ", answer[ind])
            end # can fill u-v query
        end # over queries

end # TarjanOLCA



# union-find code, modified from https://en.wikipedia.org/wiki/Tarjan%27s_off-line_lowest_common_ancestors_algorithm

#=
type SetUnion{Ti}
  parent::Array{Ti,1}
  rank::Array{Ti,1}
end # setUnion

SetUnion(n::Int32) = SetUnion([1:n], zeros(Int32,n))
SetUnion(n::Int64) = SetUnion([1:n], zeros(Int64,n))

# recursive path compression
function find!{Ti}(su::SetUnion{Ti}, x::Ti)
     if su.parent[x] == x
        return x
     else
        su.parent[x] = find!(su,su.parent[x])
        return su.parent[x]
     end
end


function union!{Ti}(su::SetUnion{Ti}, x::Ti, y::Ti)
     xroot = find!(su, x)
     yroot = find!(su, y)
     if su.rank[xroot] > su.rank[yroot]
         su.parent[yroot] = xroot
     elseif su.rank[xroot] < su.rank[yroot]
         su.parent[xroot] = yroot
     elseif xroot != yroot
         su.parent[yroot] = xroot
         su.rank[xroot] = su.rank[xroot] + 1
     end
end
  
=#

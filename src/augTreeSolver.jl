#=

Code for solving Laplacian and Diagnally Dominant Systems
by augmented spanning tree preconditioners.

Started by Dan Spielman

=#




#=========================================
  Augmented Spanning Tree Preconditioner
=========================================#



"""
    B = augmentTree{Tv,Ti}(tree, A, k)


Takes as input a tree and an adjacency matrix of a graph.
It then computes the stretch of every edge of the graph wrt
the tree.  It then adds back the k edges of highest stretch,
and k edges sampled according to stretch
"""
function augmentTree{Tv,Ti}(tree::SparseMatrixCSC{Tv,Ti}, A::SparseMatrixCSC{Tv,Ti}, k::Ti)

    st = compStretches(tree, A)

    # just to be safe, remove the tree from this
    #=
    (ti,tj) = findnz(tree)
    for i in 1:length(ti)   # most of the time
        st[ti[i],tj[i]] = 0
        st[tj[i],ti[i]] = 0
    end
    =#
    
    (ai,aj,av) = findnz(triu(st))

    ord = sortperm(av, rev=true)

    edgeinds = zeros(Bool,length(av))
    for i in 1:k
        edgeinds[ord[i]] = true
    end

    s = sum(av[ord[(k+1):end]])
    probs = av * k / s
    probs[ord[1:k]] = 0
    edgeinds[rand(length(av)) .< probs] = true

    augi = ai[edgeinds]
    augj = aj[edgeinds]
    augm = length(augi)
    augv = zeros(Tv, augm)
    for i in 1:augm,
        augv = A[augi[i],augj[i]]
    end

    n = size(A)[1]
    aug = sparse(augi, augj, augv, n, n)
    aug = aug + aug'

    return tree + aug
    
end


    
"""
    pre = augTreePrecon{Tv,Ti}(ddmat::SparseMatrixCSC{Tv,Ti}; treeAlg=akpw)

This is an augmented spanning tree preconditioner for diagonally dominant
linear systems.  It takes as optional input a tree growing algorithm.
It adds back 2sqrt(n) edges via augmentTree: the sqrt(n) of highest stretch
and another sqrt(n) sampled according to stretch.
For most purposes, one should directly call `augTreeSolver`.
"""
function augTreePrecon{Tv,Ti}(ddmat::SparseMatrixCSC{Tv,Ti}; treeAlg=akpw)

  adjmat = -triu(ddmat,1)
  adjmat = adjmat + adjmat'

  tree = treeAlg(adjmat)

  n = size(ddmat)[1]

  augtree = augmentTree(tree,adjmat,convert(Int,round(sqrt(n))))

  Dx = spdiagm(ddmat*ones(n))

  augDD = Dx + spdiagm(augtree*ones(n)) - augtree

  F = cholfact(augDD)

  return x -> (F\x)

end

"""
    solver = augTreeSolver(sddm; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[], treeAlg=akpw)

An "augmented spanning tree" solver for positive definite diagonally dominant matrices.
It works by adding edges to a low stretch spanning tree.  It calls `augTreePrecon` to form
the preconditioner.




"""
function augTreeSolver{Tv,Ti}(sddm::SparseMatrixCSC{Tv,Ti}; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[], treeAlg=akpw)

    F = augTreePrecon(sddm, treeAlg=treeAlg)
    tol_=tol
    maxits_=maxits
    maxtime_=maxtime
    verbose_=verbose
    pcgIts_=pcgIts
    
    f(b;tol=tol_,maxits=maxits_, maxtime=maxtime_, verbose=verbose_,pcgIts=pcgIts_) = pcg(sddm, b, F, tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose, pcgIts=pcgIts)
  
    return f

end

"""This is an augmented spanning tree preconditioner for Laplacians.
It takes as optional input a tree growing algorithm.
It adds back 2sqrt(n) edges via `augmentTree`: the sqrt(n) of highest stretch
and another sqrt(n) sampled according to stretch.
For most purposes, one should directly call `augTreeLapSolver`."""
function augTreeLapPrecon{Tv,Ti}(ddmat::SparseMatrixCSC{Tv,Ti}; treeAlg=akpw)

  adjmat = -triu(ddmat,1)
  adjmat = adjmat + adjmat'

  tree = treeAlg(adjmat)

  n = size(ddmat)[1]

  augtree = augmentTree(tree,adjmat,convert(Int,round(sqrt(n))))

  #Dx = spdiagm(ddmat*ones(n))

  #augDD = Dx + spdiagm(augtree*ones(n)) - augtree

  F = cholLap(augtree)

  return F

end

"""
    solver = augTreeLapSolver(A; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[], treeAlg=akpw)

An "augmented spanning tree" solver for Laplacian matrices.
It works by adding edges to a low stretch spanning tree.  It calls `augTreeLapPrecon` to form the preconditioner. In line with other solver, it takes as input the adjacency matrix of the system.
"""
function augTreeLapSolver{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti}; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[], treeAlg=akpw)

    return lapWrapComponents(augTreeLapSolver1, a, verbose=verbose, tol=tol, maxits=maxits, maxtime=maxtime, pcgIts=pcgIts, treeAlg=treeAlg)


end
     
function augTreeLapSolver1{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti}; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[], treeAlg=akpw)

  la = lap(a)

  F = augTreeLapPrecon(la, treeAlg=treeAlg)

  tol_ =tol
  maxits_ =maxits
  maxtime_ =maxtime
  verbose_ =verbose
  pcgIts_ =pcgIts


  f(b;tol=tol_,maxits=maxits_, maxtime=maxtime_, verbose=verbose_, pcgIts=pcgIts_) = pcg(la, b-mean(b), F, tol=tol, maxits=maxits, maxtime=maxtime, pcgIts=pcgIts, verbose=verbose)
    
  return f

end
     

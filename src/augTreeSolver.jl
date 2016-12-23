#=

Code for solving Laplacian and Diagnally Dominant Systems
by augmented spanning tree preconditioners.

Started by Dan Spielman

=#




#=========================================
  Augmented Spanning Tree Preconditioner
=========================================#



"""Takes as input a tree and an adjacency matrix of a graph.
It then computes the stretch of every edge of the graph wrt
the tree.  It then adds back the k edges of highest stretch,
and k edges sampled according to stretch"""
function augmentTree{Tv,Ti}(tree::SparseMatrixCSC{Tv,Ti}, mat::SparseMatrixCSC{Tv,Ti}, k::Ti)

    st = compStretches(tree, mat)

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
        augv = mat[augi[i],augj[i]]
    end

    n = size(mat)[1]
    aug = sparse(augi, augj, augv, n, n)
    aug = aug + aug'

    return tree + aug
    
end


    
"""This is an augmented spanning tree preconditioner for diagonally dominant
linear systems.  It takes as optional input a tree growing algorithm.
It adds back 2sqrt(n) edges via augmentTree: the sqrt(n) of highest stretch
and another sqrt(n) sampled according to stretch.
For most purposes, one should directly call `augTreeSolver`."""
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
An "augmented spanning tree" solver for positive definite diagonally dominant matrices.
It works by adding edges to a low stretch spanning tree.  It calls `augTreePrecon` to form
the preconditioner.

~~~julia
 augTreeSolver{Tv,Ti}(ddmat::SparseMatrixCSC{Tv,Ti}; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, treeAlg=akpw)
~~~
"""
function augTreeSolver{Tv,Ti}(ddmat::SparseMatrixCSC{Tv,Ti}; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, treeAlg=akpw)

    F = augTreePrecon(ddmat, treeAlg=treeAlg)
    tol_=tol
    maxits_=maxits
    maxtime_=maxtime
    verbose_=verbose
    
    f(b;tol=tol_,maxits=maxits_, maxtime=maxtime_, verbose=verbose_) = pcg(ddmat, b, F, tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose)
  
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
An "augmented spanning tree" solver for Laplacian matrices.
It works by adding edges to a low stretch spanning tree.  It calls `augTreeLapPrecon` to form
the preconditioner. In line with other solver, it takes as input the adjacency matrix of the system.

~~~julia
 augTreeLapSolver{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti}; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, treeAlg=akpw)
~~~
"""
function augTreeLapSolver{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti}; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, treeAlg=akpw)

  la = lap(a)

  F = augTreeLapPrecon(la, treeAlg=treeAlg)

  tol_ =tol
  maxits_ =maxits
  maxtime_ =maxtime
  verbose_ =verbose

  f(b;tol=tol_,maxits=maxits_, maxtime=maxtime_, verbose=verbose_) = pcg(la, b, F, tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose)
    
  return f

end
     
"""
A wrapper for the PyAMG solver.

~~~julia
 amgSolver{Tv,Ti}(ddmat::SparseMatrixCSC{Tv,Ti}; tol::Float64=1e-6, maxits=Inf, maxtime=Inf, verbose=false)
~~~
"""
function AMGSolver{Tv,Ti}(ddmat::SparseMatrixCSC{Tv,Ti}; tol::Float64=1e-6, maxits=Inf, maxtime=Inf, verbose=false)

  amg = PyAMG.RugeStubenSolver(ddmat);
  M = PyAMG.aspreconditioner(amg);
  function F(b)
    return M \ b;
  end

    tol_=tol
    maxits_=maxits
    maxtime_=maxtime
    verbose_=verbose

    
  f(b;tol=tol_,maxits=maxits_, maxtime=maxtime_, verbose=verbose_) = pcg(ddmat, b, F, tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose)

  return f
  
end


"""
A wrapper for the PyAMG solver. In line with our other solvers, takes in an adjacency matrix.

~~~julia
 amgSolver{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti}; tol::Float64=1e-6, maxits=Inf, maxtime=Inf, verbose=false)
~~~
"""
function AMGLapSolver{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti}; tol::Float64=1e-6, maxits=Inf, maxtime=Inf, verbose=false)

  la = lap(a)

  amg = PyAMG.RugeStubenSolver(la);
  M = PyAMG.aspreconditioner(amg);
  function F(b)
    return M \ b;
  end

    tol_=tol
    maxits_=maxits
    maxtime_=maxtime
    verbose_=verbose

    
  f(b;tol=tol_,maxits=maxits_, maxtime=maxtime_, verbose=verbose_) = pcg(la, b, F, tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose)

  return f
  
end

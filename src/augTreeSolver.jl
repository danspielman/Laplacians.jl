#=

Code for solving Laplacian and Diagnally Dominant Systems
by augmented spanning tree preconditioners.

Started by Dan Spielman

=#




#=========================================
  Augmented Spanning Tree Preconditioner
=========================================#


type AugTreeParams
    treeAlg::Function
    opt::Bool
    nnzL_fac::Float64
    flops::Float64
    verbose::Bool
end

AugTreeParams() = AugTreeParams(akpw, true, 4.0, 200.0, false)
AugTreeParamsOld() = AugTreeParams(akpw, false, 0, 0, false)


"""
    B = augmentTree{Tv,Ti}(tree, A, k)


Takes as input a tree and an adjacency matrix of a graph.
It then computes the stretch of every edge of the graph wrt
the tree.  It then adds back the k edges of highest stretch,
and k edges sampled according to stretch.

This is the old alg.  We now recommend using augmentTreeOpt.
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
    B = augmentTreeOpt{Tv,Ti}(tree, A, params)


Takes as input a tree and an adjacency matrix of a graph.
It then computes the stretch of every edge of the graph wrt
the tree.  It uses cholmod to decide how many edge to add back,
shooting for nnzL_fac times n entries in the factored augmented tree,
with a number of flops to factor equal to nnz(a)*flops_fac.
The edges to add back are then choen at random.
"""
function augmentTreeOpt{Tv,Ti}(tree::SparseMatrixCSC{Tv,Ti}, A::SparseMatrixCSC{Tv,Ti}; params=AugTreeParams())

    Aminus = A - tree
    
    ai,aj,av = findnz(triu(Aminus))

    n = A.n
    m = length(ai)

    st = compStretches(tree, Aminus)
    _,_,sv = findnz(triu(st))

    
    r = -log.(rand(m)) ./ sv
    ord = sortperm(r)

    nnzLTooBig(nnzL) = (nnzL-2*(n-1)) > n*params.nnzL_fac
    nnzLTooSmall(nnzL) = (nnzL-2*(n-1)) < n*params.nnzL_fac/2
    flopsTooBig(flops) = flops > (n+m)*params.flops
    flopsTooSmall(flops) = flops < (n+m)*params.flops/4
    
    k = Int(round(2*sqrt(n)))

    first = true
    direction = 0
    done = false

    while ~done
    
        edgeinds = ord[1:min(k,m)]
        augi = ai[edgeinds]
        augj = aj[edgeinds]
        augv = av[edgeinds]
        
        n = size(A,1)
        aug = sparse(augi, augj, augv, n, n)
        aug = aug + aug'

        augTree = tree+aug

        nnzL, flops = ask_cholmod(lap(augTree))

        if first
            first = false
            if ~(nnzLTooBig(nnzL) || nnzLTooSmall(nnzL) ||
                 flopsTooBig(flops) || flopsTooSmall(flops))
                done = true
            elseif nnzLTooBig(nnzL) || flopsTooBig(flops)
                k = div(k,2)
                direction = -1
            else
                k = k * 2
                direction = 1
                if k >= m
                    done = true
                end
                
            end
        else
            if direction == -1
                if nnzLTooBig(nnzL) || flopsTooBig(flops)
                    k = div(k,2)
                else
                    done = true
                end
                
            else # direction == 1

                if nnzLTooSmall(nnzL) && flopsTooSmall(flops)
                    k = k * 2
                    if k >= m
                        done = true
                    end
                else
                    done = true
                end
            end
        end

        if done
            if params.verbose
                println("Increased k by $(k/round(2*sqrt(n))) to $(k).")
            end

            return augTree
        end
    end
end


    
"""
    pre = augTreePrecon{Tv,Ti}(ddmat::SparseMatrixCSC{Tv,Ti}; params=AugTreeParams())

This is an augmented spanning tree preconditioner for diagonally dominant
linear systems.  It takes as optional input a tree growing algorithm.
It adds back 2sqrt(n) edges via augmentTree: the sqrt(n) of highest stretch
and another sqrt(n) sampled according to stretch.
For most purposes, one should directly call `augTreeSolver`.
"""
function augTreePrecon{Tv,Ti}(ddmat::SparseMatrixCSC{Tv,Ti};  params=AugTreeParams())

  adjmat = -triu(ddmat,1)
  adjmat = adjmat + adjmat'

  tree = params.treeAlg(adjmat)

  n = size(ddmat)[1]


  if params.opt  
      augtree = augmentTreeOpt(tree,adjmat,params=params)
  else
      augtree = augmentTree(tree,adjmat,convert(Int,round(sqrt(n))))
  end
    

  Dx = spdiagm(ddmat*ones(n))

  augDD = Dx + spdiagm(augtree*ones(n)) - augtree

  F = cholfact(augDD)

  return x -> (F\x)

end

"""
    solver = augTreeSddm(sddm; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[],  params=AugTreeParams())

An "augmented spanning tree" solver for positive definite diagonally dominant matrices.  It works by adding edges to a low stretch spanning tree.  It calls `augTreePrecon` to form the preconditioner.  `params` has entries

* `params.treeAlg` default to `akpw`
* `params.opt` if true, it interacts with cholmod to choose a good number of edges to add back.  If false, it adds back 2*sqrt(n).
"""
function augTreeSddm{Tv,Ti}(sddm::SparseMatrixCSC{Tv,Ti}; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[],  params=AugTreeParams())

    t1 = time()
    
    F = augTreePrecon(sddm; params=params)
    tol_=tol
    maxits_=maxits
    maxtime_=maxtime
    verbose_=verbose
    pcgIts_=pcgIts
    
    f(b;tol=tol_,maxits=maxits_, maxtime=maxtime_, verbose=verbose_,pcgIts=pcgIts_) = pcg(sddm, b, F, tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose, pcgIts=pcgIts)

    if verbose
        println("Solver build time: ", round((time() - t1),3), " seconds.")
    end
  
    return f

end

"""
    pre = augTreeLapPrecon{Tv,Ti}(A; params=AugTreeParams())

This is an augmented spanning tree preconditioner for Laplacians.
It takes as optional input a tree growing algorithm.
It adds back 2sqrt(n) edges via `augmentTree`: the sqrt(n) of highest stretch
and another sqrt(n) sampled according to stretch.
For most purposes, one should directly call `augTreeLapSolver`."""
function augTreeLapPrecon{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti}; params=AugTreeParams())

    tree = params.treeAlg(a)

    if params.opt  
        augtree = augmentTreeOpt(tree,a, params=params)
    else
        augtree = augmentTree(tree,a,convert(Int,round(sqrt(a.n))))
    end

    F = (h -> h)

    try
        F = cholLap(augtree)
    catch
        println("cholfact failed in augTreeLapPrecon.  Going to backup routine")
        F = augTreeFactor(augtree, tree)
    end

    return F

end

"""
    solver = augTreeLap(A; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[], params=AugTreeParams())

An "augmented spanning tree" solver for Laplacians.  It works by adding edges to a low stretch spanning tree.  It calls `augTreePrecon` to form the preconditioner.  `params` has entries

* `params.treeAlg` default to `akpw`
* `params.opt` if true, it interacts with cholmod to choose a good number of edges to add back.  If false, it adds back 2*sqrt(n).
"""
function augTreeLap{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti}; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[], params=AugTreeParams())

    return lapWrapComponents(augTreeLap1, a, verbose=verbose, tol=tol, maxits=maxits, maxtime=maxtime, pcgIts=pcgIts, params=params)


end
     
function augTreeLap1{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti}; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[], params=AugTreeParams())

  F = augTreeLapPrecon(a, params=params)

  tol_ =tol
  maxits_ =maxits
  maxtime_ =maxtime
  verbose_ =verbose
  pcgIts_ =pcgIts

  la = lap(a)

  f(b;tol=tol_,maxits=maxits_, maxtime=maxtime_, verbose=verbose_, pcgIts=pcgIts_) = pcg(la, b-mean(b), F, tol=tol, maxits=maxits, maxtime=maxtime, pcgIts=pcgIts, verbose=verbose)
    
  return f

end
     

function augTreeFactor(a, tree; verbose=false,
                      tol::Real=1e-6, maxits::Integer=1000, maxtime=Inf, pcgIts=Int[])
    n = size(a,1);

    ord::Array{Int64,1} = Laplacians.dfsOrder(tree)

    # these lines could be MUCH faster
    aord = symPermuteCSC(a,ord)
    tord = symPermuteCSC(tree,ord)
    
    laord = lap(aord)

    Aminus = aord - tord
    
    (aminusi, aminusj, aminusv) = findnz(Aminus)

    marked = ones(Int64,n)
    marked[aminusi] = 0
    marked[aminusj] = 0

    if (nnz(a) == 2*(n - 1))
        if verbose
            println("The graph is a tree.")
        end

        # mark something so this code stops somewhere
        i = 1
        j = nbri(aord, i, 1)
        marked[i] = 0
        marked[j] = 0
    end



    elims1, elims2, ind::Array{Int64,1}, subtree = elimDeg12(tord, marked)
    
    map = zeros(Int64,n)

    n1 = length(ind)

    map[ind] = collect(1:n1)
    mapi = map[aminusi]
    mapj = map[aminusj]

    rest = sparse(mapi, mapj, aminusv, n1,n1)
    a1 = rest + subtree
    la1 = lap(a1)

    fsub = cholLap(a1)


    f1 = function(b::Array{Float64,1})

        subMean!(b) # b = b - mean(b)
        
        y = forwardSolve(b, elims1, elims2)
        ys = y[ind]

        xs = fsub(ys)
        
        x = zeros(Float64,n)
        x[ind] = xs

        backSolve(x, y, elims1, elims2)
        subMean!(x) # x = x - mean(x)

        return x
    end


    f = function(b)

        bord = b[ord] 

        xord = f1(bord)
        
        #xord = pcg(laord, bord, f1, tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose, pcgIts=pcgIts)

        x = zeros(Float64,n)
        x[ord] = xord - mean(xord)

        return x
    end

    return f

end


#=
  A hybrid that reduces an augmented spanning tree,
  and then uses a variant of the sampling solver.

  Uses a lot of code from KMPSolver

=#

type AugTreeHybridParams
    n0::Int64      # the number of edges at which to go direct
    frac::Float64  # add frac*n edges to the tree
end

AugTreeHybridParams() = AugTreeHybridParams(100, 0.1)

function AugTreeHybridLap(a; verbose=false,
                      tol::Real=1e-6, maxits::Integer=1000, maxtime=Inf, pcgIts=Int[], params::AugTreeHybridParams=AugTreeHybridParams())


        return lapWrapComponents(AugTreeHybridLap1, a, verbose=verbose, tol=tol, maxits=maxits, maxtime=maxtime, pcgIts=pcgIts, params=params)

end

function AugTreeHybridLap1(a; verbose=false,
                      tol::Real=1e-6, maxits::Integer=1000, maxtime=Inf, pcgIts=Int[], params::AugTreeHybridParams=AugTreeHybridParams())

    if (a.n <= params.n0)
        if verbose
            println("The graph is small.  Solve directly")
        end
        
        return cholLap(a)
    end

    if (nnz(a) == 2*(a.n - 1))
        if verbose
            println("The graph is a tree.  Solve directly")
        end
        
        return cholLap(a)
    end

    tree = akpw(a)
    if verbose
        println("akpw stretch : ", sum(compStretches(tree,a))/nnz(a))
    end

    n = size(a,1);

    ord::Array{Int64,1} = Laplacians.dfsOrder(tree)

    # these lines could be MUCH faster
    aord = symPermuteCSC(a,ord)
    tord = symPermuteCSC(tree,ord)
    
    la = lap(aord)


end



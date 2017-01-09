#=
  A hybrid that reduces an augmented spanning tree,
  and then uses a variant of the sampling solver.

  Uses a lot of code from KMPSolver

=#

#include("types.jl")

AugTreeHybridParams() = AugTreeHybridParams(100, 0.05, 10, 0.5, true)

import Laplacians.IJVS
import Laplacians.elimDeg12
import Laplacians.subMean!
import Laplacians.forwardSolve
import Laplacians.backSolve


function augTreeHybridLap(a; verbose=false,
                      tol::Real=1e-6, maxits::Integer=1000, maxtime=Inf, pcgIts=Int[], params::AugTreeHybridParams=AugTreeHybridParams())

        return lapWrapComponents(augTreeHybridLap1, a, verbose=verbose, tol=tol, maxits=maxits, maxtime=maxtime, pcgIts=pcgIts, params=params)

end

function augTreeHybridLap1(a; verbose=false,
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

    n = size(a,1);

    ord::Array{Int64,1} = Laplacians.dfsOrder(tree)

    # these lines could be MUCH faster
    aord = symPermuteCSC(a,ord)
    tord = symPermuteCSC(tree,ord)
    
    laord = lap(aord)

    # Note: this really needs an ijv, not an ijvs
    ijvs1 = augmentingEdges(aord, tord, params)  

    marked = ones(Int64,n)
    marked[ijvs1.i] = 0
    marked[ijvs1.j] = 0

    elims1, elims2, ind::Array{Int64,1}, subtree = elimDeg12(tord, marked)

    map = zeros(Int64,n)

    n1 = length(ind)
    
    map[ind] = collect(1:n1)
    ijvs1.i = map[ijvs1.i]
    ijvs1.j = map[ijvs1.j]

    rest = sparse(ijvs1.i,ijvs1.j,ijvs1.v,n1,n1)
    a1 = rest + rest' + subtree
    la1 = lap(a1)

    tol_=tol
    maxits_=maxits
    maxtime_=maxtime
    verbose_=verbose
    pcgIts_=pcgIts


    fsub = simpleSamplerLap1(a1,tol=params.subTol,maxits=params.subIts)


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


    f = function(b; tol=tol_, maxits=maxits_, maxtime=maxtime_, verbose=verbose_, pcgIts=pcgIts_)

        bord = b[ord] - mean(b)
        
        xord = pcg(laord, bord, f1, tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose, pcgIts=pcgIts)

        x = zeros(Float64,n)
        x[ord] = xord
        subMean!(x) # x = x - mean(x)
        return x
    end


        
    return f

end


function augmentingEdges(a, tree, params)

    Aminus = a - tree
    
    ai,aj,av = findnz(triu(Aminus))

    m = length(ai)

    st = compStretches(tree, Aminus)
    _,_,sv = findnz(triu(st))

    k = min(m, round(Int,params.frac*a.n))

    sampleProbs = computePforX(sv,k)
    edgeinds = find(rand(size(sampleProbs)) .< sampleProbs)
    augi = ai[edgeinds]
    augj = aj[edgeinds]
    augv = av[edgeinds] 

    if params.verbose
        println("Expected $(sum(sampleProbs)) edges.  Got $(length(edgeinds)) edges.")
    end

    ijvs = IJVS(augi, augj, augv, augv)

    return ijvs
end


"""
    p = computePforX(x,k::Int)

For sampling k items with probabilities proportional to x.
We choose r according to rand(), and then keep those for which r <= p.
For example

```julia
x = collect(1:10).^2
p = computePforX(x,5)
s = find(rand(size(p)) .< p)
```   
""" 
function computePforX{Tv}(x::Array{Tv},k::Int)

    if (k == 0) || (length(x) == 0)
        return Tv[]
    end
    
    
    xs = sort(x, rev=true)

    csx = sum(xs)
    i = 0
    while xs[i+1]*(k-i) > csx
        csx -= xs[i+1]
        i = i + 1
    end    
    fac = (k-i)/csx 

    p = map(z->min(1,fac*z),x)
     
    return p
end



#=
  TEST FAILING: wtedChimera(1000, 5)
   from runtests.
  So, this is removed from the list of solvers for now.

    A mix between KMP, augTree and samplingSolver.
    At the top level, we perform a sampling scheme similar to what we see in augTreeSolver. After 
    sampling and eliminating degree 1 and degree 2 vertices, we perform a low accuracy solve using
    the samplingSolver, without blowing up the low stretch tree inside. 

Calls a bunch of routines from KMPSolver.jl
=#

global HYBRID_SAVEMATS=false

global HYBRID_MATS=[]
global HYBRID_FS=[]

"""Parameters for the hybrid solver"""
type hybridParams
    frac::Float64   # fraction to decrease at each level
    iters::Int64    # iters of PCG to apply between levels
    treeAlg         # :akpw or :rand
    n0::Int64       # the number of edges at which to go direct

    ssParams::samplingParams
end

defaultHybridParams = hybridParams(1/200, 15, :akpw, 600, 
                        samplingParams(0.5,0.2,1.0,1000,20,false,false,1e-3))



"""
A mix between KMP, augTree and samplingSolver.
At the top level, we perform a sampling scheme similar to what we see in augTreeSolver. After 
sampling and eliminating degree 1 and degree 2 vertices, we perform a low accuracy solve using
the samplingSolver, without blowing up the low stretch tree inside. 

Solves linear equations in symmetric, diagonally dominant matrices with non-positive off-diagonals.

~~~julia
hybridSDDSolver(mat; verbose=false, tol::Real=1e-2, maxits::Integer=1000, maxtime=Inf, params::hybridParams=defaultHybridParams)
~~~
"""
function hybridSDDSolver(mat; verbose=false,
                      tol::Real=1e-2, maxits::Integer=1000, maxtime=Inf, params::hybridParams=defaultHybridParams)

    n = size(mat,1)
    s = mat*ones(n)

    dmat = diag(mat)
    s = sparse(max(s,0) .* (s .> (dmat*1e-12)))

    if (s == 0)
        error("Matrix was not diagonally dominant.")
    end
    
    # Force symmetric and diagonal zero
    a = triu(abs(mat),1)
    a = a + a'
    
    a1 = [sparse([0 s']); [s a]]
    
    f1 = hybridLapSolver(a1, verbose=verbose, tol=tol, maxits=maxits, maxtime=maxtime, params=params)

    f = function(b::Array{Float64,1})

        b1 = [-sum(b);b]
        x1 = f1(b1)
        x = x1[2:end] - x1[1]
        
        return x
        
    end
    
    return f
end


"""
A mix between KMP, augTree and samplingSolver.
At the top level, we perform a sampling scheme similar to what we see in augTreeSolver. After 
sampling and eliminating degree 1 and degree 2 vertices, we perform a low accuracy solve using
the samplingSolver, without blowing up the low stretch tree inside. 

Solves linear equations in the Laplacian of graph with adjacency matrix `a`.


~~~julia
hybridLapSolver(a; verbose=false, tol::Real=1e-2, maxits::Integer=1000, maxtime=Inf, params::hybridParams=defaultHybridParams)
~~~
"""
function hybridLapSolver(a; verbose=false,
                      tol::Real=1e-2, maxits::Integer=1000, maxtime=Inf, params::hybridParams=defaultHybridParams)

    if (minimum(a) < 0)
        error("The adjacency matrix cannot have negative entries.")
    end

    co = components(a)
    if maximum(co) == 1
        return hybridLapSolver1(a, verbose=verbose, tol=tol, maxits=maxits, maxtime=maxtime, params=params)
    else
        comps = vecToComps(co)

        if verbose
            println("The graph has $(length(comps)) components.")
        end

        solvers = []
        for i in 1:length(comps)
            ind = comps[i]
            
            asub = a[ind,ind]

            if (length(ind) == 1)
                error("Node $ind has no edges.")
            end
            
            subSolver = hybridLapSolver1(asub, verbose=verbose, tol=tol, maxits=maxits, maxtime=maxtime, params=params)

            push!(solvers, subSolver)
        end

        return Laplacians.blockSolver(comps,solvers)
    end
    
end


# hybridLapSolver drops right in to this after doing some checks and splitting on components
function hybridLapSolver1(a; verbose=false,
                      tol::Real=1e-2, maxits::Integer=1000, maxtime=Inf, params::hybridParams=defaultHybridParams)

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

    if params.treeAlg == :rand
        tree = randishPrim(a)
    else
        tree = akpw(a)
    end

    n = size(a,1);

    # if for some reason the graph is a tree, this will fail
    # because the stretches will sum to zero
    # so, default to a direct method

    ord::Array{Int64,1} = Laplacians.dfsOrder(tree)

    # these lines could be MUCH faster
    aord = symPermuteCSC(a,ord)
    tord = symPermuteCSC(tree,ord)
    
    la = lap(aord)

    if HYBRID_SAVEMATS
        HYBRID_MATS = []
        push!(HYBRID_MATS,la)
    end

    fsub = hybridLapPrecon(aord, tord, params, verbose=verbose)

    f = function(b::Array{Float64,1})
        bord = b[ord]
        
        xord = pcg(la, bord, fsub, tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose)

        x = zeros(Float64,n)
        x[ord] = xord
        subMean!(x) # x = x - mean(x)
        return x
    end

    if HYBRID_SAVEMATS
        HYBRID_FS = []
        push!(HYBRID_FS,f)
    end

    return f

end


function hybridLapPrecon(a, tree, params; verbose=false)
    n = size(a,1);

    rest = a-tree;
    st = compStretches(tree,rest);
    
    (ai,aj,av) = findnz(triu(rest))
    (si,sj,sv) = findnz(triu(st))

    ijvs = IJVS(ai,aj,av,sv)

    f = hybridLapPreconSub(tree, ijvs, 0, params, verbose=verbose)
    
    return f
end


function hybridLapPreconSub(tree, ijvs::IJVS, level::Int64, params::hybridParams; verbose=false)

    # problem: are forming la before sampling.  should be other way around, at least for top level!
    # that is, we are constructing Heavy, and I don't want to!

    m = size(ijvs.i,1)
    n = size(tree,1)

    if verbose
        println("level ", level, ". Dimension ", n, " off-tree edges : ", m)
    end

    # if is nothing in ijvs
    if m == 0
        return cholLap(tree)
    end

    ijvs1 = stretchSample(ijvs,params.frac)
    
    if level > 0

        rest = sparse(ijvs1.i,ijvs1.j,ijvs1.v,n,n)
        adjMat = rest + rest' + tree
        la = lap(rest + rest' + tree)

        F = samplingLapSolver(adjMat, tol=1e-1, params=params.ssParams)

        f = function(b::Array{Float64,1})
            auxb = copy(b);
            subMean!(auxb)

        	x = F(auxb)
            subMean!(x)

        	return x
        end
    else

        marked = ones(Int64,n)
        marked[ijvs1.i] = 0
        marked[ijvs1.j] = 0

        elims1, elims2, ind::Array{Int64,1}, subtree = elimDeg12(tree, marked)

        map = zeros(Int64,n)

        n1 = length(ind)
        
        map[ind] = collect(1:n1)
        ijvs1.i = map[ijvs1.i]
        ijvs1.j = map[ijvs1.j]

        rest = sparse(ijvs1.i,ijvs1.j,ijvs1.v,n1,n1)
        la1 = lap(rest + rest' + subtree)

        fsub = hybridLapPreconSub(subtree, ijvs1, level+1, params, verbose=verbose)

        f = function(b::Array{Float64,1})
            subMean!(b) # b = b - mean(b)

            y = forwardSolve(b, elims1, elims2)
            ys = y[ind]

            xs = pcg(la1, ys, fsub, tol=0, maxits=params.iters)
            
            x = zeros(Float64,n)
            x[ind] = xs

            backSolve(x, y, elims1, elims2)
            subMean!(x) # x = x - mean(x)

            return x
        end
    end

    if HYBRID_SAVEMATS
        push!(HYBRID_FS,f)
    end

    return f
    
end


#=
    Sample k = m * frac edges in the following way:
        1. take the highest stretch 1/8k edges
        2. sample 7/8k edges proportional to stretch
=#
function stretchSample(ijvs::IJVS,frac::Float64)

    m = size(ijvs.i,1)
    k = ceil(Int64, m * frac)

    take1 = ceil(Int64, 3/8 * k)
    take2 = ceil(Int64, 5/8 * k)

    # make sure we sample less than m edges
    take2 = min(take2, m - take1)

    ord = sortperm(ijvs.s, rev=true)
    edgeinds = zeros(Bool,length(ijvs.i))

    # sample take2 edges form take1+1 to m proportional to stretch
    s = sum(ijvs.s[ord[(take1+1):end]])
    probs = ijvs.s * take2 / s

    # take the first k edges (the ones of highest stretch)
    probs[ord[1:take1]] = 1

    edgeinds[rand(m) .< probs] = true

    sampi = ijvs.i[edgeinds]
    sampj = ijvs.j[edgeinds]
    sampv = ijvs.v[edgeinds]
    samps = ijvs.s[edgeinds]

    return IJVS(sampi,sampj,sampv, samps)

end

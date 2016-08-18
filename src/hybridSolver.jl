#=
    A mix between KMP, augTree and samplingSolver.
    At the top level, we perform a sampling scheme similar to what we see in augTreeSolver. After 
    sampling and eliminating degree 1 and degree 2 vertices, we perform a low accuracy solve using
    the samplingSolver, without blowing up the low stretch tree inside. 
=#

global HYBRID_SAVEMATS=false

global HYBRID_MATS=[]
global HYBRID_FS=[]

type IJVS
    i::Array{Int64,1}
    j::Array{Int64,1}
    v::Array{Float64,1}  # wt of edge
    s::Array{Float64,1}  # stretch of edge
end


type hybridParams
    frac::Float64   # fraction to decrease at each level
    iters::Int64    # iters of PCG to apply between levels
    treeAlg         # :akpw or :rand
    n0::Int64       # the number of edges at which to go direct

    ssParams::samplingParams
end

defaultHybridParams = hybridParams(1/200, 15, :akpw, 600, 
                        samplingParams(0.5,0.2,1.0,1000,20,false,false,false,1e-3))

# this is just for Laplacians, not general SDD
immutable elimLeafNode
    nodeid::Int64
    parent::Int64
    wtDeg::Float64
end

# this is just for Laplacians, not general SDD
immutable elimDeg2Node
    nodeid::Int64
    nbr1::Int64
    nbr2::Int64
    wt1::Float64
    wt2::Float64
end


# The tree must be in DFS order
# marked is 1 if flagged for possible elimination,
# and set to 2 if we do eliminate it
# is just for Laplacians matrices, not general SDD    
function elimDeg12(t, marked)

    # make sure root is not marked
    marked[1] = 0

    n = size(t,1)

    deg = Array{Int64}(n)
    for v in 1:n
        deg[v] = t.colptr[v+1] - t.colptr[v]
    end

    elims1 = Array{elimLeafNode}(0)

    for v in n:-1:2

        if (deg[v] == 1 && marked[v] == 1)
            parent = t.rowval[t.colptr[v]];
            wt = t.nzval[t.colptr[v]];
            push!(elims1,elimLeafNode(v,parent,wt))

            deg[parent] = deg[parent] - 1
            marked[v] = 2
            deg[v] = 0
        end
    end

    elims2 = Array{elimDeg2Node}(0)

    subt = triu(t)

    for v in n:-1:2

        if (deg[v] == 2 && marked[v] == 1)

            parent = t.rowval[t.colptr[v]];

            # to ident the child, enumerate to find one uneliminated, which we check by marked
            kidind = t.colptr[v]+1
            kid = t.rowval[kidind]
            while deg[kid] == 0
                kidind = kidind+1
                kid = t.rowval[kidind]

                if kidind >= t.colptr[v+1]
                    error("went of the end of the kid list without finding node to elim from")
                end
            end


            wt1 = t.nzval[t.colptr[v]];
            wt2 = t.nzval[kidind]

            push!(elims2,elimDeg2Node(v,parent,kid,wt1,wt2))
            marked[v] = 2

            newwt = 1/(1/wt1 + 1/wt2)

            # now that we've found the kid, go up the chain until done
            while (deg[parent] == 2 && marked[parent] == 1)
                v = parent
                parent = t.rowval[t.colptr[v]];
                wt1 = t.nzval[t.colptr[v]];
                wt2 = newwt

                push!(elims2,elimDeg2Node(v,parent,kid,wt1,wt2))
                marked[v] = 2
                
                newwt = 1/(1/wt1 + 1/wt2)
            end

            # now, hack the tree to adjust parent and wt of node kid
            subt.rowval[subt.colptr[kid]] = parent
            subt.nzval[subt.colptr[kid]] = newwt

        end
    end

    subt = subt + subt'

    ind = find(marked.<2)
    subt = subt[ind,ind]
    
    return elims1, elims2, ind, subt
end




function forwardSolve(b, elims1, elims2)

    y = copy(b)

    for i in 1:length(elims1)
        y[elims1[i].parent] += y[elims1[i].nodeid]
    end

    for i in 1:length(elims2)
        wtsum = elims2[i].wt1 + elims2[i].wt2
        y[elims2[i].nbr1] += y[elims2[i].nodeid]*elims2[i].wt1 / wtsum
        y[elims2[i].nbr2] += y[elims2[i].nodeid]*elims2[i].wt2 / wtsum
    end

    return y
    
end

function backSolve(x, y, elims1, elims2)
    
    for i in length(elims2):-1:1
        node = elims2[i].nodeid
        wtsum = elims2[i].wt1 + elims2[i].wt2

        x[node] = (elims2[i].wt1*x[elims2[i].nbr1] + elims2[i].wt2*x[elims2[i].nbr2] + y[node])/wtsum
    end

    
    for i in length(elims1):-1:1
        node = elims1[i].nodeid
        x[node] = x[elims1[i].parent] + y[node]/elims1[i].wtDeg
    end

end


# subtract off the mean from a vector, in place
function subMean!(x::Array{Float64,1})
    n = size(x,1)
    mn = mean(x)
    for i in 1:n,
        x[i] = x[i] - mn
    end
end



"""Solves linear equations in symmetric, diagonally dominant matrices with non-positive off-diagonals."""
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


"""Solves linear equations in the Laplacian of graph with adjacency matrix `a`."""
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
        
        return lapWrapSolver(cholfact, lap(a))
    end

    if (nnz(a) == 2*(a.n - 1))
        if verbose
            println("The graph is a tree.  Solve directly")
        end
        
        return lapWrapSolver(cholfact, lap(a))
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
        la = lap(tree)
        return lapWrapSolver(cholfact, la)
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

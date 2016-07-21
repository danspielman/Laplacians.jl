#=

An implementation of the Laplacians and SDD solvers of Koutis, Miller and Peng
=#

include("fastCSC.jl")
include("samplingSolver.jl")

global KMP_SAVEMATS=false

global KMP_MATS=[]
global KMP_FS=[]

type IJVS
    i::Array{Int64,1}
    j::Array{Int64,1}
    v::Array{Float64,1}  # wt of edge
    s::Array{Float64,1}  # stretch of edge
end


type KMPparams
    frac::Float64  # fraction to decrease at each level
    iters::Int64   # iters of PCG to apply between levels
    treeScale::Float64 # scale tree by treeScale*log_2 (n) * aveStretch
    n0::Int64 # the number of edges at which to go direct
    treeAlg # :akpw or :rand
end

defaultKMPparams = KMPparams(1/36, 6, 0.125, 600, :akpw)

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
function KMPSDDSolver(mat; verbose=false,
                      tol::Real=1e-2, maxits::Integer=1000,  params::KMPparams=defaultKMPparams)

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
    
    f1 = KMPLapSolver(a1, verbose=verbose, tol=tol, maxits=maxits, params=params)

    f = function(b::Array{Float64,1})

        b1 = [-sum(b);b]
        x1 = f1(b1)
        x = x1[2:end] - x1[1]
        
        return x
        
    end
    
    return f
end


"""Solves linear equations in the Laplacian of graph with adjacency matrix `a`."""
function KMPLapSolver(a; verbose=false,
                      tol::Real=1e-2, maxits::Integer=1000,  params::KMPparams=defaultKMPparams)

    if (minimum(a) < 0)
        error("The adjacency matrix cannot have negative entries.")
    end

    co = components(a)
    if maximum(co) == 1
        return KMPLapSolver1(a, verbose=verbose, tol=tol, maxits=maxits,  params=params)
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
            
            subSolver = KMPLapSolver1(asub, verbose=verbose, tol=tol, maxits=maxits,  params=params)

            push!(solvers, subSolver)
        end

        return Laplacians.blockSolver(comps,solvers)
    end
    
end



# KMPLapSolver drops right in to this after doing some checks and splitting on components
function KMPLapSolver1(a; verbose=false,
                      tol::Real=1e-2, maxits::Integer=1000,  params::KMPparams=defaultKMPparams)

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

    if KMP_SAVEMATS
        KMP_MATS = []
        push!(KMP_MATS,la)
    end

    fsub = KMPLapPrecon(aord, tord, params, verbose=verbose)

#    f = function(b::Array{Float64,1}; tol::Real=tol, maxits::Integer=maxits)
    f = function(b::Array{Float64,1})
        bord = b[ord]
        
        xord = pcg(la, bord, fsub, tol=tol, maxits=maxits, verbose=verbose)

        x = zeros(Float64,n)
        x[ord] = xord
        subMean!(x) # x = x - mean(x)
        return x
    end

    if KMP_SAVEMATS
        KMP_FS = []
        push!(KMP_FS,f)
    end

    return f

end


function KMPLapPrecon(a, t, params; verbose=false)
    n = size(a,1);

    rest = a-t;

    st = compStretches(t,rest);
    aveStretch = sum(st)/nnz(rest)

    if params.treeScale > 0
    
        targetStretch = 1/(params.treeScale*log(n)/log(2))

        fac = aveStretch/targetStretch
        tree = fac*t;

    else

        targetStretch = 1.0

        fac = aveStretch/targetStretch
        tree = t;

    end

    if verbose
        println("aveStretch : ", aveStretch, " fac : ", fac)
    end
    
    (ai,aj,av) = findnz(triu(rest))
    (si,sj,sv) = findnz(triu(st))
    sv = sv ./ fac

    ijvs = IJVS(ai,aj,av,sv)

    f = KMPLapPreconSub(tree, ijvs, targetStretch, 0, params, verbose=verbose)
    
    return f
end


function KMPLapPreconSub(tree, ijvs::IJVS, targetStretch::Float64, level::Int64, params::KMPparams; verbose=false)

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

    ijvs1 = stretchSample(ijvs,targetStretch,params.frac)
    
    if (length(ijvs1.i) <= params.n0 || level > 0)

        rest = sparse(ijvs1.i,ijvs1.j,ijvs1.v,n,n)
        adj = rest + rest' + tree
        la = lap(rest + rest' + tree)

        F = samplingSolver(adj, tol=1e-1, eps=0.5, sampConst=0.02, beta=1e3, verbose=verbose)
        # F = amgSolver(la)

        f = function(b::Array{Float64,1})
        	subMean!(b)
        	x = F(b)
        	subMean!(x)

        	return x
        end

        # println(adj.n, " ", length(f(rand(adj.n))))

        # f = lapWrapSolver(cholfact, la)
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

        fsub = KMPLapPreconSub(subtree, ijvs1, targetStretch, level+1, params, verbose=verbose)

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

    if KMP_SAVEMATS
        push!(KMP_FS,f)
    end

    return f
    
end


#=
goal is to downsample edges by frac.
However, those whose stretches are larger than stretchTarget
might just get their weight reduced instead

This will mostly be the "uniform sampling" envisioned by KMP.
=#
function stretchSample(ijvs::IJVS,stretchTarget::Float64,frac::Float64)

    m = size(ijvs.i,1)
    k = ceil(Int64, m * frac / 2)

    ord = sortperm(ijvs.s, rev=true)
    edgeinds = zeros(Bool,length(ijvs.i))

    # sample k edges form k+1 to m proportional to stretch
    s = sum(ijvs.s[ord[(k+1):end]])
    probs = ijvs.s * k / s

    # take the first k edges (the ones of highest stretch)
    probs[ord[1:k]] = 1

    edgeinds[rand(m) .< probs] = true

    sampi = ijvs.i[edgeinds]
    sampj = ijvs.j[edgeinds]
    sampv = ijvs.v[edgeinds]
    samps = ijvs.s[edgeinds]

    return IJVS(sampi,sampj,sampv, samps)

end

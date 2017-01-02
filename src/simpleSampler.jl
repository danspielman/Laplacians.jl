#=
    An implementation of the linear system solver in https://arxiv.org/abs/1605.02353 by Rasmus Kyng and Sushant Sachdeva.

    This is an optimistic heuristic implementation by DAS.

    In addition to the setup in the paper, we also use a low stretch tree to approximate effective 
    resistances on edges. To perform well cache wise, we implement a cache friendly list of
    linked lists - found in revampedLinkedListFloatStorage.jl 
    
    The parameters used in the solver are described below.


=#

"""
Parameters for the sampling solver.
"""
type SimpleSamplerParams{Tv,Ti}
    startingSize::Ti    # the initial size of the linked list storage structure
    blockSize::Ti       # the size of each consecutive block of memory assigned to a certain element

    # debug parameters, we can get rid of them later
    verboseSS::Bool
    returnCN::Bool
    CNTol::Tv
    perm::Symbol    # options :tree, :amd
    fixA::Bool
end

SimpleSamplerParams() = SimpleSamplerParams(1000, 20,
                                       false,false,1e-3,:tree,false)
#=
""" 
    solver = samplingSDDMSolver(sddm)

An implementation of the linear system solver in https://arxiv.org/abs/1605.02353 by Rasmus Kyng and Sushant Sachdeva. In addition to the setup in the paper, we also use a low stretch tree to approximate effective resistances on edges. To perform well cache wise, we implement a cache friendly list of linked lists - found in revampedLinkedListFloatStorage.jl 
    """

function simpleSamplerSDDM{Tv,Ti}(sddm::SparseMatrixCSC{Tv,Ti}; tol::Tv=1e-6, maxits=1000, maxtime=Inf, verbose=false, pcgIts=Int[], params::samplingParams{Tv,Ti}=defaultSamplingParams)

    # srand(1234)

    adjMat,diag = adj(sddm)

    a = extendMatrix(adjMat,diag)
    n = a.n

    F,gOp,_,_,ord,cn,cntime = buildSolver(a, params=params)

    if params.verboseSS
        println()
        println("Eps error: ", checkError(gOp))
        println("Condition number: ", cn)
        println()
    end 

    tol_=tol
    maxits_=maxits
    maxtime_=maxtime
    verbose_=verbose
    pcgIts_=pcgIts

    la = lap(symPermuteCSC(a, ord))
    f = function(b; tol=tol_, maxits=maxits_, maxtime=maxtime_, verbose=verbose_, pcgIts=pcgIts_)
        #= 
            We need to add an extra entry to b to make it match the size of a. The extra vertex in a is
            vertex n, thus, we will add the new entry in b on position n as well.
        =#
        auxb = copy(b)
        if norm(diag) != 0
            push!(auxb, -sum(auxb))
        end

        ret = pcg(la, auxb[ord], F, tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose, pcgIts=pcgIts_)
        ret = ret[invperm(ord)]

        # We want to discard the nth element of ret (which corresponds to the first element in the permutation)
        ret = ret - ret[n]

        if norm(diag) != 0
            pop!(ret)
        end

        return ret
    end
    
    return f

end
=#


""" 
    solver = samplingLapSolver(A)

An implementation of the linear system solver in https://arxiv.org/abs/1605.02353 by Rasmus Kyng and Sushant Sachdeva. In addition to the setup in the paper, we also use a low stretch tree to approximate effective resistances on edges. To perform well cache wise, we implement a cache friendly list of linked lists - found in revampedLinkedListFloatStorage.jl 
"""
function simpleSamplerLap{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti}; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[], params::SimpleSamplerParams{Tv,Ti}=SimpleSamplerParams())

    return lapWrapComponents(simpleSamplerLap1, a, verbose=verbose, tol=tol, maxits=maxits, maxtime=maxtime, pcgIts=pcgIts, params=params)


end


function simpleSamplerLap1{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti}; tol::Tv=1e-6, maxits=1000, maxtime=Inf, verbose=false, pcgIts=Int[], params::SimpleSamplerParams{Tv,Ti}=defaultsimpleSamplerParams)

    # srand(1234)

	n = a.n;

    F,gOp,_,_,ord,cn,cntime = buildSolver(a, params=params)

    if params.verboseSS
        println()
        println("Eps error: ", checkError(gOp))
        println("Condition number: ", cn)
        println()
    end 

    la = lap(symPermuteCSC(a, ord))

    tol_=tol
    maxits_=maxits
    maxtime_=maxtime
    verbose_=verbose
    pcgIts_=pcgIts


    f = function(b; tol=tol_, maxits=maxits_, maxtime=maxtime_, verbose=verbose_, pcgIts=pcgIts_)

        ret = pcg(la, b[ord] - mean(b), F, tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose, pcgIts=pcgIts)
        return ret[invperm(ord)]
    end
    
    return f

end

# Add a new vertex to a with weights to the other vertices corresponding to diagonal surplus weight
function extendMatrix{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti}, diag::Array{Tv,1})

    if norm(diag) == 0
        return a
    end
    
    n = a.n
    u,v,w = findnz(a)
    for i in 1:n
        if diag[i] > 0
            push!(u, i)
            push!(v, n + 1)
            push!(w, diag[i])

            push!(u, n + 1)
            push!(v, i)
            push!(w, diag[i])
        end
    end
    
    return sparse(u,v,w)

end

function buildSolver{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti};
                            params::SimpleSamplerParams{Tv,Ti}=SimpleSamplerParams())

    # compute rho
    n = a.n;

    if params.verboseSS
	    tic()
    end
    #=
        Note: The n'th vertex will correspond to the vertex used to solve Laplacian systems
        with extra weight on the diagonal. Thus, we want to eliminate it last.
    =#

    tree = akpw(a);

    if params.perm == :tree
        ord = reverse!(dfsOrder(tree, start = n));
    else
        ord = cholmod_perm(lap(a))
    end
    
    a = symPermuteCSC(a, ord)
    tree = symPermuteCSC(tree, ord)

    if params.verboseSS
	    println("Nonzeros in a = ", nnz(a))
    end

    if params.verboseSS
	    print("Time to build the tree: ")
	    toc()
    end

    if params.verboseSS
	    tic()
    end

    # Get u and d such that ut d u = -a (doesn't affect solver)

    Ut,d = simpleSampleTree(a, tree, params.startingSize, params.blockSize, params.verboseSS)



    U = Ut'

    if params.verboseSS
	print("Time to factorize: ")
	toc()
    end



    if params.verboseSS
	    println("nnz in U matrix: ", length(U.data.nzval))
    end

    # Create the solver function
    f = function(b::Array{Float64,1})
        # center
        res = copy(b)
        res = res - sum(res) / n

        # forward solve
        res = Ut \ res

        # diag inverse
        for i in 1:(n - 1)
            res[i] = res[i] / d[i]
        end

        # backward solve
        res = U \ res

        # center
        res = res - sum(res) / n
        
        return res
    end

    # Create the error check function
    la = lap(a)   
    g = function(b)
        res = copy(b)   
        res[n] = 0
            
        # diag sqrt inverse 
        for i in 1:(n - 1)
            res[i] = res[i] * d[i]^(-1/2)
        end

        # backward solve #TODO?
        res = U \ res

        # apply lapl
        res = la * res

        # forward solve #TODO?
        res = Ut \ res

        # diag sqrt inverse
        for i in 1:(n - 1)
            res[i] = res[i] * d[i]^(-1/2)
        end

        # subtract identity, except we haven't zeroed out last coord
        res = res - b 

        #zero out last coord
        res[n] = 0 #TODO?
            
        return res
    end
    gOp = SqLinOp(true,1.0,n,g)

    if params.returnCN
	    tic()
        cn = condNumber(lap(a2), U, d, tol=params.CNTol)
        print("computing the condition number takes: ")
	    cntime = toc()
    end

    if params.returnCN
        return f,gOp,U,d,ord,cn,cntime
    else
        return f,gOp,U,d,ord,(0.0,0.0),0.0
    end
end

function simpleSampleTree{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti}, tree::SparseMatrixCSC{Tv,Ti}, startingSize::Ti, blockSize::Ti, verbose::Bool=false)

    n = a.n

    # some extra memory to be used later in the algorithm. this can be later pulled out of this function
    # into an external recipient, to be used on subsequent runs of the solver
    auxVal = zeros(Tv, n)                       # used to sum weights from multiedges
    auxMult = zeros(Tv, n)                      # used to count the number of multiedges

    wNeigh = zeros(Tv, n)
    multNeigh = zeros(Tv, n)
    indNeigh = zeros(Ti, n)
    
    ut = Array{Tuple{Tv,Ti},1}[[] for i in 1:n]  # the lower triangular u matrix part of u d u'
    d = zeros(Tv, n)                            # the d matrix part of u d u'

    # neigh[i] = the list of neighbors for vertex i with their corresponding weights
    # note neigh[i] only stores neighbors j such that j > i
    # neigh[i][1] is weight, [2] is number of multi-edges, [3] is neighboring vertex

    neigh = llsInit(a, startingSize = startingSize, blockSize = blockSize)

    # gather the info in a and put it into neigh and w
    for i in 1:length(a.colptr) - 1
        for j in a.colptr[i]:a.colptr[i + 1] - 1
            if a.rowval[j] > i
                llsAdd(neigh, i, (a.nzval[j], 1.0, a.rowval[j]))
            end
        end
    end

#    tr = matToTree(tree,n)


    if params.fixA
        if params.verboseSS
            println("fix A")
        end
        
        check = a
    else
        if params.verboseSS
            println("fix tree")
        end

        check = tree
        
    end
    

    # Now, for every i, we will compute the i'th column in U
    for i in 1:(n-1)
        # We will get rid of duplicate edges
        # wSum - sum of weights of edges
        # multSum - sum of number of edges (including multiedges)
        # numPurged - the size in use of wNeigh, multNeigh and indNeigh
        # wNeigh - list of weights correspongind to each neighbors
        # multNeigh - list of number of multiedges to each neighbor
        # indNeigh - the indices of the neighboring vertices
        wSum, multSum, numPurged = llsPurge(neigh, i, auxVal, auxMult, wNeigh, multNeigh, indNeigh, rho = 0.0) # also a lot of time

        # need to divide weights by the diagonal entry
        for j in 1:numPurged
            push!(ut[i], (-wNeigh[j] / wSum, indNeigh[j]))
        end
        push!(ut[i], (1, i)) #diag term

        d[i] = wSum

        deg = numPurged

        nbrs = indNeigh[1:deg]
        
        # handle all tree edges
        for j in 1:(deg-1)
            for k in (j+1):deg
                posj = indNeigh[j]
                posk = indNeigh[k]

                if check[posj,posk] > 0

#                if tr.parent[posj] == posk || tr.parent[posk] == posj

                    # swap so posj is smaller
                    if posk < posj  
                        j, k = k, j
                        posj, posk = posk, posj
                    end

                    wj = wNeigh[j]                
                    wk = wNeigh[k]

                    sampScaling = wSum
                    
                    llsAdd(neigh, posj, (wj * wk / sampScaling, 1.0, posk))
                end
            end
        end

        if deg > 0
            wSamp = FastSampler(wNeigh[1:deg])
        
            jSamples = sampleMany(wSamp, deg)
            kSamples = randperm(deg)
        end

        # now propagate the clique to the neighbors of i
        for l in 1:deg
            
            j = jSamples[l]
            k = kSamples[l]

            posj = indNeigh[j]
            posk = indNeigh[k]

            if (j != k) && (check[posj,posk] == 0)

                # swap so posj is smaller
                if posk < posj  
                    j, k = k, j
                    posj, posk = posk, posj
                end

                wj = wNeigh[j]                
                wk = wNeigh[k]

                sampScaling = wj + wk 
                
                llsAdd(neigh, posj, (wj * wk / sampScaling, 1.0, posk))
            end
        end  

    end

    # add the last diagonal term
    push!(ut[n], (1, n))
    d[n] = 0

    if verbose
	    println()
	    println("The actual size is ", neigh.size * neigh.blockSize)
	    println()
    end

    return constructLowerTriangularMat(ut), d
end


# u is an array of arrays of tuples. to be useful, we need to convert it to a lowerTriangular matrix
# looks like 10% of the construction time, and could be sped up a lot
function constructLowerTriangularMat{Tv,Ti}(u::Array{Array{Tuple{Tv,Ti},1},1})
    n = length(u)

    nnz = 0
    for i in 1:n
        nnz = nnz + length(u[i])
    end

    colptr = Array{Ti,1}(n + 1)
    rowval = Array{Ti,1}(nnz)
    nzval = Array{Tv,1}(nnz)

    colptr[1] = 1
    for i in 1:n
        colptr[i + 1] = colptr[i] + length(u[i])
    end
    index = copy(colptr)

    # We know that in u the values aren't necessarily ordered by row. So, we do a count sort-like algorith to keep linear time.
    helper = Array{Tuple{Tv,Ti},1}[[] for i in 1:n]
    for i in 1:n
        for j in 1:length(u[i])
            row = u[i][j][2]
            col = i
            val = u[i][j][1]
            push!(helper[row], (val, col))
        end
    end

    for i in 1:n
        for j in 1:length(helper[i])
            row = i
            col = helper[i][j][2]
            val = helper[i][j][1]

            rowval[index[col]] = row
            nzval[index[col]] = val
            index[col] = index[col] + 1
        end
    end

    return LowerTriangular(SparseMatrixCSC(n, n, colptr, rowval, nzval))
end

function checkError{Tv,Ti}(gOp::SqLinOp{Tv,Ti}; tol::Float64 = 0.0)
    return eigs(gOp;nev=1,which=:LM,tol=tol)[1][1]
end

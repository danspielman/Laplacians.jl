# Rasmus' linear system sovler. 

#=
    Some historically good parameters:
    sampConst in [0.001, 0.03]
    beta in [30,200]
=#

using Laplacians

include("fastSampler.jl")
include("johnlind.jl")
include("revampedLinkedListFloatStorage.jl")
include("sqLinOpWrapper.jl")
include("fastCSC.jl")
include("condNumber.jl")

function samplingSolver{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti}, diag::Array{Tv,1};
                                tol::Tv=1e-6, maxits=100, maxtime=Inf, verbose::Bool = false,
                                eps::Tv = 0.5, sampConst::Tv = 0.02, beta::Tv = 1000.0,
                                JLEps::Tv = 0.5, JLSolver=(la -> augTreeSolver(la,tol=1e-2,maxits=1000,maxtime=10)),
                                startingSize::Ti = 1000, blockSize::Ti = 20,
                                returnCN::Bool = false)

    a = extendedLaplacian(a, diag)
    n = a.n

    F,gOp,_,_,ord,cn,cntime = buildSolver(a, eps = eps, sampConst = sampConst, beta = beta, 
        JLEps = JLEps, JLSolver = JLSolver,
        startingSize = startingSize, blockSize = blockSize, 
        returnCN = returnCN, verbose = verbose)

    if verbose
        println()
        println("Eps error: ", checkError(gOp))
        println("Condition number: ", cn)
        println()
    end 

    invperm = zeros(Int64, n)
    for i in 1:n
    	invperm[ord[i]] = i
    end

    la = lap(sympermute(a, ord))
    function f(b)
        #= 
            We need to add an extra entry to b to make it match the size of a. The extra vertex in a is
            vertex n, thus, we will add the new entry in b on position n as well.
        =#
        auxb = copy(b)
        push!(auxb, -sum(auxb))

        ret = pcg(la, auxb[ord], F, tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose)
        ret = ret[invperm]

        # We want to discard the nth element of ret (which corresponds to the first element in the permutation)
        ret = ret - ret[n]

        return ret[1:(n-1)]
    end
    
    return f

end

# Add a new vertex to a with weights to the other vertices corresponding to diagonal surplus weight
function extendedLaplacian{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti}, diag::Array{Tv,1})

    n = a.n
    u,v,w = findnz(a)
    for i in 1:n
        push!(u, i)
        push!(v, n + 1)
        push!(w, diag[i])

        push!(u, n + 1)
        push!(v, i)
        push!(w, diag[i])
    end
    
    return sparse(u,v,w)

end

function buildSolver{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti};
                            eps::Tv = 0.5, sampConst::Tv = 0.02, beta::Tv = 1000.0,
                            JLEps::Tv = 0.5, JLSolver=(la -> augTreeSolver(la,tol=1e-2,maxits=1000,maxtime=10)),
                            startingSize::Ti = 1000, blockSize::Ti = 20,
                            returnCN::Bool = false, verbose::Bool = false)

    # compute rho
    n = a.n;
    rho::Tv = sampConst * log(n) ^ 2 / eps ^ 2
    if verbose
	    println("rho = ", rho)
    end

    if verbose
	    tic()
    end
    #=
        Note: The n'th vertex will correspond to the vertex used to solve Laplacian systems
        with extra weight on the diagonal. Thus, we want to eliminate it last.
    =#

    # augment the tree
    tree = akpw(a);
    augTree = augmentTree(tree,a,convert(Int,round(sqrt(n))))

    ord = reverse!(dfsOrder(tree, start = n));

    a = sympermute(a, ord)
    augTree = sympermute(augTree, ord)

    # Blow up the augTree
    a2 = copy(a)

    a = a + (beta - 1) * augTree
    xhat = johnlind(a, eps=JLEps, solver=JLSolver, retXhat=true)

    # stretch = rho * compStretches(beta * tree, a)
    # stretch.nzval = min(rho, stretch.nzval)

  	if verbose
	    print("Time to build the augTree and run JL : ")
	    toc()
    end

    # Get u and d such that ut d u = -a (doesn't affect solver)
    Ut,d = samplingLDL(a, xhat, rho, startingSize, blockSize, verbose)
    U = Ut'

    if verbose
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
    g = function(b::Array{Float64,1})
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

    if returnCN
	    tic()
        cn = condNumber(lap(a2), U, d)
        print("computing the condition number ", cn, " takes: ")
	    cntime = toc()
    end

    if verbose
        print("Preconditioner has ", nnz(U.data), " nonzeros")
    end

    if returnCN
        return f,gOp,U,d,ord,cn,cntime
    else
        return f,gOp,U,d,ord,(0.0, 0.0),0.0
    end
end

# a is an adjacency matrix
function samplingLDL{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti}, xhat::Array{Tv,2}, rho::Tv,
    startingSize::Ti, blockSize::Ti, verbose::Bool=false)

    n = a.n

    # later will have to do a permutation here, for now consider the matrix is already permuted

    # some extra memory to be used later in the algorithm. this can be later pulled out of this function
    # into an external recipient, to be used on subsequent runs of the solver
    auxVal = zeros(Tv, n)                       # used to sum weights from multiedges
    auxMult = zeros(Tv, n)                      # used to count the number of multiedges

    wNeigh = zeros(Tv, n)
    multNeigh = zeros(Tv, n)
    indNeigh = zeros(Ti, n)
    
    u = Array{Tuple{Tv,Ti},1}[[] for i in 1:n]  # the lower triangular u matrix part of u d u'
    d = zeros(Tv, n)                            # the d matrix part of u d u'

    # neigh[i] = the list of neighbors for vertex i with their corresponding weights
    # note neigh[i] only stores neighbors j such that j > i
    # neigh[i][1] is weight, [2] is number of multi-edges, [3] is neighboring vertex

    neigh = llsInit(a, startingSize = startingSize, blockSize = blockSize)

    totalLev = 0

    # gather the info in a and put it into neigh and w
    for i in 1:length(a.colptr) - 1
        # push!(neigh, Tuple{Tv,Ti,Ti}[])

        for j in a.colptr[i]:a.colptr[i + 1] - 1
            if a.rowval[j] > i
                p = a.rowval[j]
                q = i
                lev = a.nzval[j] * norm(xhat[p,:] - xhat[q,:])^2

                totalLev += lev

                llsAdd(neigh, i, (a.nzval[j], min(rho, rho * lev), a.rowval[j]))
            end
        end
    end

    if verbose
        println("Sum of leverage scores: ", totalLev)
    end

    # Now, for every i, we will compute the i'th column in U
    for i in 1:(n - 1)
        # We will get rid of duplicate edges
        # wSum - sum of weights of edges
        # multSum - sum of number of edges (including multiedges)
        # numPurged - the size in use of wNeigh, multNeigh and indNeigh
        # wNeigh - list of weights correspongind to each neighbors
        # multNeigh - list of number of multiedges to each neighbor
        # indNeigh - the indices of the neighboring vertices
        wSum, multSum, numPurged = llsPurge(neigh, i, auxVal, auxMult, wNeigh, multNeigh, indNeigh, rho = rho)

        # need to divide weights by the diagonal entry
        for j in 1:numPurged
            push!(u[i], (-wNeigh[j] / wSum, indNeigh[j]))
        end
        push!(u[i], (1, i)) #diag term
        d[i] = wSum

        multSum = ceil(Int64, multSum)
        wSamp = FastSampler(wNeigh[1:numPurged])
        multSamp = FastSampler(multNeigh[1:numPurged])
        
        jSamples = sampleMany(wSamp, multSum)
        kSamples = sampleMany(multSamp, multSum)
        
        # now propagate the clique to the neighbors of i
        for l in 1:multSum
            # newSeed = rand(UInt32)
            # srand(newSeed)
            
            j = jSamples[l]
            k = kSamples[l]

            if j != k
                posj = indNeigh[j]
                posk = indNeigh[k]

                # swap so posj is smaller
                if posk < posj  
                    j, k = k, j
                    posj, posk = posk, posj
                end

                wj = wNeigh[j]                
                wk = wNeigh[k]

                sampScaling = wj * multNeigh[k] + wk * multNeigh[j]
                
                llsAdd(neigh, posj, (wj * wk / sampScaling, 1.0, posk))
            end
        end  
    end

    # add the last diagonal term
    push!(u[n], (1, n))
    d[n] = 0

    if verbose
	    println()
	    println("The size of the linked list data structure is ", neigh.size * neigh.blockSize)
	    println()
    end

    return constructLowerTriangularMat(u), d
end

# u is an array of arrays of tuples. to be useful, we need to convert it to a lowerTriangular matrix
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
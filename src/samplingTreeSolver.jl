# Rasmus' linear system sovler. 

using Laplacians

include("fastSampler.jl")
include("linkedListStorage.jl")
include("sqLinOpWrapper.jl")

function samplingSolver{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti}; tol::Float64=1e-6, maxits::Int64=100, 
                                eps::Float64 = 0.5, sampConst::Float64 = 10.0, k::Float64 = 5, beta::Float64 = 5)

    a2 = copy(a)

    F,_,_,_,ord = buildSolver(a, eps = eps, sampConst = sampConst, k = k, beta = beta)

    invperm = collect(1:n)
    sort!(invperm, by=x->ord[x])

    la = lap(a[ord,ord])
    function f(b)
        ret = pcg(la, b[ord], F, tol=tol, maxits=maxits, verbose=true)
        return ret[invperm]
    end
    
    return f

end

function checkError{Tv,Ti}(gOp::SqLinOp{Tv,Ti})
    return eigs(gOp;nev=1,which=:LM)[1][1]
end

function buildSolver{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti}; 
                            eps::Float64 = 0.5, sampConst::Float64 = 10.0, k::Float64 = 5, beta::Float64 = 6)

    n = a.n;
    rho = ceil(Ti, sampConst * log(n) ^ 2 / eps ^ 2)

    tic()

    tree = randishKruskal(a);
    ord = reverse!(dfsOrder(tree));

    a = a[ord, ord];
    tree = tree[ord, ord];
    rootedTree = matToTree(tree);

    # blow up the tree by beta
    a = a + beta * tree

    depth = compDepth(rootedTree);
    stretch = tarjanStretch(rootedTree, a, depth);

    # I think Dan's code also multiplies stretch by edge weight
    stretch = ceil(Int64, stretch / k);

    # set the tree strech to rho
    #### NOT EFFICIENT
    for u in 1:n
        for i in 1:deg(tree,u)
            v = nbri(tree,u,i)

            stretch[u,v] = rho
        end
    end

    print("Time to build the tree and compute the stretch: ")
    toc()

    # Get u and d such that u d u' = -a (doesn't affect solver)
    U,d = samplingLDL(a, stretch, rho)
    Ut = U'

    # Create the solver function
    f = function(b::Array{Float64,1})
        # center
        res = copy(b)
        res = res - sum(res) / n

        # forward solve
        res = U \ res

        # diag inverse
        for i in 1:(n - 1)
            res[i] = res[i] / d[i]
        end

        # backward solve
        res = Ut \ res

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
        res = Ut \ res

        # apply lapl
        res = la * res

        # forward solve #TODO?
        res = U \ res

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

    return f,gOp,U,d,ord
end

# a is an adjacency matrix
function samplingLDL{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti}, stretch::SparseMatrixCSC{Ti,Ti}, rho::Ti)
    n = a.n

    # later will have to do a permutation here, for now consider the matrix is already permuted

    # some extra memory to be used later in the algorithm. this can be later pulled out of this function
    # into an external recipient, to be used on subsequent runs of the solver
    auxVal = zeros(Tv, n)                       # used to sum weights from multiedges
    auxMult = zeros(Ti, n)                      # used to counte the number of multiedges

    wNeigh = zeros(Tv, n)
    multNeigh = zeros(Ti, n)
    indNeigh = zeros(Ti, n)
    
    u = Array{Tuple{Tv,Ti},1}[[] for i in 1:n]  # the lower triangular u matrix part of u d u'
    d = zeros(Tv, n)                            # the d matrix part of u d u'

    # neigh[i] = the list of neighbors for vertex i with their corresponding weights
    # note neigh[i] only stores neighbors j such that j > i
    # neigh[i][1] is weight, [2] is number of multi-edges, [3] is neighboring vertex

    # TODO: change after this works for the current behavior
    neigh = llsInit(a, rho)

    # gather the info in a and put it into neigh and w
    for i in 1:length(a.colptr) - 1
        # push!(neigh, Tuple{Tv,Ti,Ti}[])

        for j in a.colptr[i]:a.colptr[i + 1] - 1
            if a.rowval[j] > i
                # set the value min between rho and the stretch
                llsAdd(neigh, i, (a.nzval[j], min(rho, stretch.nzval[j]), a.rowval[j]))
            end
        end
    end

    # Now, for every i, we will compute the i'th column in U
    for i in 1:(n - 1)
        # We will get rid of duplicate edges
        # wSum -  sum of weights of edges
        # wNeigh - list of weights correspongind to each neighbors
        # multSum - sum of number of edges (including multiedges)
        # multNeigh - list of number of multiedges to each neighbor
        # indNeigh - the indices of the neighboring vertices

        wSum, multSum, numPurged = llsPurge(neigh, i, auxVal, auxMult, wNeigh, multNeigh, indNeigh)
        
        # println(i, " ", numPurged)
        # println(wNeigh')
        # println(multNeigh')
        # println(indNeigh')
        # println()

        # need to divide weights by the diagonal entry
        for j in 1:numPurged
            push!(u[i], (-wNeigh[j] / wSum, indNeigh[j]))
        end
        push!(u[i], (1, i)) #diag term
        d[i] = wSum

        # wSamp = sampler(wNeigh)
        # multSamp = sampler(convert(Array{Tv,1}, multNeigh))
        
        # jSamples = sampleMany(wSamp, multSum)
        # kSamples = sampleMany(multSamp, multSum)

        wSamp = sampler(wNeigh[1:numPurged])
        jSamples = sampleMany(wSamp, multSum)

        kSamples = Ti[]
        ind = 0
        for j in 1:numPurged
            for k in 1:multNeigh[j]
                push!(kSamples, j)
            end
        end
        
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
                
                # push!(neigh[posj], (wj * wk / sampScaling, 1, posk))
                llsAdd(neigh, posj, (wj * wk / sampScaling, 1, posk))
            end
        end  
    end

    # add the last diagonal term
    push!(u[n], (1, n))
    d[n] = 0

    return constructLowerTriangularMat(u), d
end

# u is an array of arrays of tuples. to be useful, we need to convert it to a lowerTriangular matrix
function constructLowerTriangularMat{Tv,Ti}(u::Array{Array{Tuple{Tv,Ti},1},1})
    n = length(u)

    nnz::Ti = 0
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
    index::Array{Ti,1} = copy(colptr)

    # We know that in u the values aren't necessarily ordered by row. So, we do a count sort-like algorith to keep linear time.
    helper = Array{Tuple{Tv,Ti},1}[[] for i in 1:n]
    for i in 1:n
        for j in 1:length(u[i])
            row::Ti = u[i][j][2]
            col::Ti = i
            val::Tv = u[i][j][1]
            push!(helper[row], (val, col))
        end
    end

    for i in 1:n
        for j in 1:length(helper[i])
            row::Ti = i
            col::Ti = helper[i][j][2]
            val::Tv = helper[i][j][1]

            rowval[index[col]] = row
            nzval[index[col]] = val
            index[col] = index[col] + 1
        end
    end

    return LowerTriangular(SparseMatrixCSC(n, n, colptr, rowval, nzval))
end
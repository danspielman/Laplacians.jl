# Rasmus' linear system sovler. 

using Laplacians

include("fastSampler.jl")
include("sqLinOpWrapper.jl")

function samplingSolver{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti}; tol::Float64=1e-6, maxits::Int64=100, eps::Float64 = 0.5, sampConst::Float64 = 10.0)

    F = buildSolver(a, eps = eps, sampConst = sampConst)

    la = lap(a)
    f(b) = pcg(la, b, F, tol=tol, maxits=maxits, verbose=true)

    return f

end

function checkError{Tv,Ti}(gOp::SqLinOp{Tv,Ti})
    return eigs(gOp;nev=1,which=:LM)[1][1]
end

function buildSolver{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti}; eps::Float64 = 0.5, sampConst::Float64 = 10.0)

    # Get u and d such that u d u' = -a (doesn't affect solver)
    u,d = samplingLDL(a, eps, sampConst)

    # Initialize the data structures
    n = length(d)
    m = n

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

    # Get U as sparse matrix
    U = LowerTriangular(SparseMatrixCSC(n, n, colptr, rowval, nzval))
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

    return f,gOp,U,d
end

# a is an adjacency matrix
function samplingLDL{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti}, eps::Float64, sampConst::Float64)
    n = a.n

    rho = ceil(Ti, sampConst * log(n) ^ 2 / eps ^ 2)

    # later will have to do a permutation here, for now consider the matrix is already permuted

    # some extra memory to be used later in the algorithm. this can be later pulled out of this function
    # into an external recipient, to be used on subsequent runs of the solver
    auxVal = zeros(Tv, n)                       # used to sum weights from multiedges
    auxMult = zeros(Ti, n)                      # used to counte the number of multiedges
    
    u = Array{Tuple{Tv,Ti},1}[[] for i in 1:n]  # the lower triangular u matrix part of u d u'
    d = zeros(Tv, n)                            # the d matrix part of u d u'

    # neigh[i] = the list of neighbors for vertex i with their corresponding weights
    # note neigh[i] only stores neighbors j such that j > i
    # neigh[i][1] is weight, [2] is number of multi-edges, [3] is neighboring vertex
    neigh = Array{Tuple{Tv,Ti,Ti},1}[]

    # gather the info in a and put it into neigh and w
    for i in 1:length(a.colptr) - 1
        push!(neigh, Tuple{Tv,Ti,Ti}[])

        for j in a.colptr[i]:a.colptr[i + 1] - 1
            if a.rowval[j] > i
                push!(neigh[i], (a.nzval[j], rho, a.rowval[j]))
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
        wSum, wNeigh, multSum, multNeigh, indNeigh = purge(i, neigh[i], auxVal, auxMult)
        
        # need to divide weights by the diagonal entry
        for j in 1:length(indNeigh)
            push!(u[i], (-wNeigh[j] / wSum, indNeigh[j]))
        end
        push!(u[i], (1, i)) #diag term
        d[i] = wSum

        # wSamp = sampler(wNeigh)
        # multSamp = sampler(convert(Array{Tv,1}, multNeigh))
        
        # jSamples = sampleMany(wSamp, multSum)
        # kSamples = sampleMany(multSamp, multSum)

        wSamp = sampler(wNeigh)
        jSamples = sampleMany(wSamp, multSum)

        kSamples = Ti[]
        ind = 0
        for j in 1:length(indNeigh)
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
                
                push!(neigh[posj], (wj * wk / sampScaling, 1, posk))
            end
        end  
    end

    # add the last diagonal term
    push!(u[n], (1, n))
    d[n] = 0

    return u, d
    # return u, d, sparse(tomatrix(u, d))
end

# see description in samplingSolver
function purge{Tv,Ti}(col, v::Array{Tuple{Tv,Ti,Ti},1}, auxVal::Array{Tv,1}, auxMult::Array{Ti,1})

    multSum::Ti = 0
    diag::Tv = 0
    for i in 1:length(v)
        neigh = v[i][3]
        w = v[i][1]
        e = v[i][2]

        auxVal[neigh] += w
        diag += w
        auxMult[neigh] += e
        multSum += e
    end

    res = Tv[]
    mult = Ti[]
    ind = Ti[]
    
    for i in 1:length(v)
        neigh = v[i][3]

        if auxVal[neigh] != 0
            @assert(col < neigh, "col < neigh in purge")
            if neigh == col
                @assert(false, "col = neigh in purge")
            else
                push!(res, auxVal[neigh])
                push!(mult, auxMult[neigh])
                push!(ind, neigh)
                auxVal[neigh] = 0
                auxMult[neigh] = 0
            end
        end
    end

    return diag, res, multSum, mult, ind
end

# this is a debug function, it returns u' * d * u
# function tomatrix{Tv,Ti}(u::Array{Array{Tuple{Tv,Ti},1}}, d::Array{Tv,1})
#     n = length(d)
#     matu = zeros(n, n)
#     matd = zeros(n, n)

#     for i in 1:n
#       matd[i, i] = d[i]
#     end

#     for i in 1:n
#       for j in 1:length(u[i])
#           pos = u[i][j][2]
#           val = u[i][j][1]
#           matu[pos, i] = val
#       end
#     end

#     return matu * matd * matu'

# end

# Rasmus' linear system sovler. 

using Laplacians

include("fastSampler.jl")

# a is a laplacian, b is the target answer vector. We return x
function samplingSolver{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti}, b::Array{Tv,1})
    n = a.n
    #a = lap(a)

    eps::Tv = 1/2
    sampConst::Tv = 10
    # Theorem ? 9*4
    # TODO
    rho = ceil(Ti,sampConst*log(n)^2/eps^2)

    # later will have to do a permutation here, for now consider the matrix is already permuted

    # some extra memory to be used later in the algorithm. this can be later pulled out of this function
    # into an external recipient, to be used on subsequent runs of the solver
    auxVal = zeros(Tv, n) 						# used to sum weights from multiedges
    auxMult = zeros(Ti, n) 						# used to counte the number of multiedges
    
    u = Array{Tuple{Tv,Ti},1}[[] for i in 1:n]	# the lower triangular u matrix part of u d u'
    d = zeros(Tv, n)							# the d matrix part of u d u'

    # neigh[i] = the list of neighbors for vertex i with their corresponding weights
    # note neigh[i] only stores neighbors j such that j > i
    # neigh[i][1] is weight, [2] is number of multi-edges, [3] is neighboring vertex
    neigh = Array{Tuple{Tv,Ti,Ti},1}[]

    # gather the info in a and put it into neigh and w
    for i in 1:length(a.colptr) - 1
		push!(neigh, [])

		for j in a.colptr[i]:a.colptr[i + 1] - 1
		    if a.rowval[j] > i
				push!(neigh[i], (a.nzval[j], rho, a.rowval[j]))
		    end
		end
    end

    # Now, for every i, we will compute the i'th column in U
    for i in 1:(n-1)
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

        wSamp = sampler(wNeigh)
        multSamp = sampler(convert(Array{Tv,1}, multNeigh))

		# now propagate the clique to the neighbors of i
        for l in 1:multSum
            j = sample(wSamp)
            k = sample(multSamp)
            if j != k
                if indNeigh[k] < indNeigh[j]  #swap so posj is smaller
                    j,k = k,j
                end
                
                posj = indNeigh[j]
                wj = wNeigh[j]
                
                posk = indNeigh[k]
                wk = wNeigh[k]

                assert(posj < posk) #remove eventually

                sampScaling = wj*multNeigh[k] + wk*multNeigh[j]
                
                push!(neigh[posj], (wj * wk/sampScaling,1,posk))
            end
        end
    end

    push!(u[n], (1, n)) #diag term
    d[n] = 0

    return u, d, sparse(tomatrix(u, d))
end

# see description in samplingSolver
function purge{Tv,Ti}(col, v::Array{Tuple{Tv,Ti,Ti},1}, auxVal::Array{Tv,1}, auxMult::Array{Ti,1})

    for i in 1:length(v)
    	neigh = v[i][3]
		auxVal[neigh] = 0
        auxMult[neigh] = 0
    end
    # RAS: Don't we maintain as invariant that this is zeroed out at the end?

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
	        assert(col < neigh)
		    if neigh == col
				assert(false)
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
function tomatrix{Tv,Ti}(u::Array{Array{Tuple{Tv,Ti},1}}, d::Array{Tv,1})
    n = length(d)
    matu = zeros(n, n)
    matd = zeros(n, n)

    for i in 1:n
		matd[i, i] = d[i]
    end

    for i in 1:n
		for j in 1:length(u[i])
		    pos = u[i][j][2]
		    val = u[i][j][1]
		    matu[i, pos] = val
		end
    end

    return matu' * matd * matu

end

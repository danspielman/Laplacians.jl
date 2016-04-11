# Rasmus' linear system sovler. 

using Laplacians

include("sampler.jl")

# a is an adjacency matrix here, we can rewrite this to have a laplacian fed in
# right now the solver jsut performs and ldl decomposition of the matrix a
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
    auxVal = zeros(Tv, n) 				# used to get rid of multiedges
    auxMult = zeros(Ti, n) 				# used to get rid of multiedges
    
    u = Array{Tuple{Tv,Ti},1}[[] for i in 1:n]	# the lower triangular u matrix part of u d u'
    d = zeros(Tv, n)				# the d matrix part of u d u'

    # neigh[i] = the list of neighbors for vertex i with their corresponding weights
    # note, neigh[i] only stores neighbors j such that j >= i (>= because we also want the diagonal)
    # neigh[i][1] is weight, [2] is number of multi-edges, [3] is neigh vertex
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

    # now, for every i, retain the ith column times
    for i in 1:(n-1)
	# get rid of duplicate vertices, and return the list of neighbors and the diagonal entry
	# updtNeigh, dEntry, multNeigh, multSum = purge(i, neigh[i], auxVal, auxMult)
        
        wSum, wNeigh, multSum, multNeigh, indNeigh = purge(i, neigh[i], auxVal, auxMult)
        # println(i,"      ",indNeigh)
        # println(i," ofo: ",find(indNeigh .<= i))
        # println(i," ok:  ",isSortedAdjList(neigh))
        # if wSum <= 0 # debugging
        #     @printf "index = %f \n" i
        # end
        
	# need to divide weights by the diagonal entry
	for j in 1:length(indNeigh)
	    push!(u[i], (-wNeigh[j] / wSum, indNeigh[j]))
	end
        push!(u[i], (1, i)) #diag term
	d[i] = wSum

        wSamp = Sampler(wNeigh)
        multSamp = Sampler(convert(Array{Tv,1},multNeigh))

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

# gets rid of duplicate entries in v. returns v with no duplicate entries and the diagonal entry for the current column
function purge{Tv,Ti}(col, v::Array{Tuple{Tv,Ti,Ti},1}, auxVal::Array{Tv,1}, auxMult::Array{Ti,1})

    for i in 1:length(v)
	auxVal[v[i][3]] = 0
        auxMult[v[i][3]] = 0
    end
    # RAS: Don't we maintain as invariant that this is zeroed out at the end?

    multSum::Ti = 0
    diag::Tv = 0
    for i in 1:length(v)
	auxVal[v[i][3]] += v[i][1]
        diag += v[i][1]
        auxMult[v[i][3]] += v[i][2]
        multSum += v[i][2]
    end

    res = Tv[]
    mult = Ti[]  
    ind = Ti[]
    
    for i in 1:length(v)
	if auxVal[v[i][3]] != 0
            assert(col < v[i][3])
	    if v[i][3] == col
		assert(false)
	    else
                push!(res, auxVal[v[i][3]])
                push!(mult, auxMult[v[i][3]])
                push!(ind, v[i][3])
		auxVal[v[i][3]] = 0
                auxMult[v[i][3]] = 0
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


function isSortedAdjList{Tv,Ti}(neigh::Array{Array{Tuple{Tv,Ti,Ti},1},1})
    n = length(neigh)
    isSorted = 0
    for i in 1:n
	for j in 1:length(neigh[i])
	    pos = neigh[i][j][3]
            posOK = (pos > i)
	    if !posOK
                isSorted +=1
            end
	end
    end
    return isSorted
end


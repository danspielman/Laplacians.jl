# Rasmus & Sushant's linear system sovler. 

using Laplacians

# a is an adjacency matrix here, we can rewrite this to have a laplacian fed in
# right now the solver jsut performs and ldl decomposition of the matrix a
function basicSolver{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti}, b::Array{Tv,1})
    n = a.n

    eps::Tv = 1/2
    sampConst::Tv = 1
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
    # note, neigh[i] only stores neighbors j such that j > i
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
        
        indNeigh, wSum, wNeigh, multSum, multNeigh = purge(i, neigh[i], auxVal, auxMult)
        # println(i,"  ",indNeigh)
        # println(i,"  ",find(indNeigh .<= i))
	# in the future, from updtNeigh we will do a sampling - for now, we just use it as is

	# need to divide weights by the diagonal entry
	for j in 1:length(indNeigh)
	    push!(u[i], (-wNeigh[j] / wSum, indNeigh[j]))
	end
        push!(u[i], (1, i)) #diag term
	d[i] = wSum
        
	for j in 1:length(indNeigh)
	    for k in 1:length(indNeigh)
		posj = indNeigh[j]
		wj = wNeigh[j] / wSum

		posk = indNeigh[k]
		wk = wNeigh[k] / wSum

                if posj < posk
		    push!(neigh[posj], (wSum * wj * wk,1,posk))
		end
	    end
	end
    end

    push!(u[n], (1, n)) #diag term
    d[n] = 0

    # important - here u has 0 on the diagonal. tomatrix adds 1ns for testing
    return u, d, sparse(tomatrix(u, d))

end

# gets rid of duplicate entries in v. returns v with no duplicate entries and the diagonal entry for the current column
function purge{Tv,Ti}(col, v::Array{Tuple{Tv,Ti,Ti},1}, auxVal::Array{Tv,1}, auxMult::Array{Ti,1})

    for i in 1:length(v)
	auxVal[v[i][3]] = 0
        auxMult[v[i][3]] = 0
    end
    # TODO: Don't we maintain as invariant that this is zeroed out at the end?

    multSum::Ti = 0
    diag::Tv = 0
    for i in 1:length(v)
	auxVal[v[i][3]] += v[i][1]
        diag += v[i][1]
        # TODO wrong below
        auxMult[v[i][3]] += 1
        multSum += 1
    end

    res = Tv[]
    ind = Ti[]
    mult = Tv[]   

    for i in 1:length(v)
	if auxVal[v[i][3]] != 0
            # assert(v[i][3] != col) #There should not be a diagonal entry
            push!(res,auxVal[v[i][3]])
            push!(ind,v[i][3])
            push!(mult,auxMult[v[i][3]])
	    auxVal[v[i][3]] = 0
            auxMult[v[i][3]] = 0
	end
    end

    return ind, diag, res, multSum, mult
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

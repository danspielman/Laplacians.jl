#=
Compute the condition number given the initial graph and the preconditioner.
Uses cholfact, so, if the initial graph is not sdd or if inverse has a lot of nonzeros this fails. 
=#
function condNumber{Tv,Ti}(la::SparseMatrixCSC{Tv,Ti}, prec::SparseMatrixCSC{Tv,Ti}; tol = 1e-6)

	cf = cholfact(prec)
	return computeCN(la, LowerTriangular(sparse(cf[:L]))', ones(Tv, la.n), tol = tol)

end

# M = Ut * d * U
function condNumber{Tv,Ti}(la::SparseMatrixCSC{Tv,Ti}, U::UpperTriangular{Tv,SparseMatrixCSC{Tv,Ti}}, d::Array{Tv,1}; tol = 1e-6)

	return computeCN(la, U, d, tol = tol)

end

# M = Ut * U
function condNumber{Tv,Ti}(la::SparseMatrixCSC{Tv,Ti}, U::UpperTriangular{Tv,SparseMatrixCSC{Tv,Ti}}; tol = 1e-6)

	return computeCN(la, U, ones(a.n), tol = tol)

end

function computeCN{Tv,Ti}(la::SparseMatrixCSC{Tv,Ti}, U::UpperTriangular{Tv,SparseMatrixCSC{Tv,Ti}}, d::Array{Tv,1}; tol = 1e-6)

	n = la.n

	# get the max eigenvalue of d^-1/2 * L^-1 * a * U^-1 * d^-1/2
	f = function(b::Array{Float64,1})
		res = copy(b)

		# diag sqrt inverse 
		res = res .* d.^-0.5
		res[n] = 0

		res = U \ res
		res = la * res
		res = U' \ res

		# diag sqrt inverse 
        res = res .* d.^-0.5
        res[n] = 0

		return res
	end
	fOp = SqLinOp(true,1.0,n,f)
	lmax = eigs(fOp;nev=1,which=:LM,tol=0.0)[1][1]

	# get the max eigenvalue of 1 - M / lambdaMax
	g = function(b::Array{Float64,1})
		res2 = copy(b)
		res = copy(b) / lmax

		# diag sqrt inverse
		res = res .* d.^-0.5
		res[n] = 0

		res = U \ res
		res = la * res
		res = U' \ res

		# diag sqrt inverse
		res = res .* d.^-0.5
		res[n] = res2[n] = 0

		return res2 - res
	end
	gOp = SqLinOp(true,1.0,n,g)

	R = checkError(gOp, tol = tol)

    Kmin = 1 / (1 - R)
    Kmax = 1 / (1 - R - tol)

    return (Kmin, Kmax)

end
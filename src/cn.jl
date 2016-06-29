#=
Compute the condition number given the initial graph and the preconditioner.
Uses cholfact, so, if the initial graph is not sdd or if inverse has a lot of nonzeros this fails. 
=#
function cn{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti}, prec::SparseMatrixCSC{Tv,Ti})

	cf = cholfact(prec)

end

function cn{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti}, L, U, d)

end

function cn{Tv,Ti}(a::SparseMatrixCSC, L, U)
end

function computeCN{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti}, L, U; tol = 2e-6)

	n = a.n

	# get the max eigenvalue of L^-1 * a * U^-1
	f = function(b::Array{Float64,1})
		res = copy(b)
		res[n] = 0

		res = U \ res
		res = a * res
		res = L \ res

		res[n] = 0
		return res
	end
	fOp = SqLinOp(true,1.0,n,f)
	lmax = eigs(fOp;nev=1,which=:LM,tol=0.0)[1][1]

	# get the max eigenvalue of 1 - M / lambdaMax
	g = function(b::Array{Float64,1})
		res2 = copy(b)
		res = copy(b) / lmax

		res[n] = 0

		res = U \ res
		res = a * res
		res = L \ res

		res[n] = res2[n] = 0

		return res2 - res
	end
	gOp = SqLinOp(true,1.0,n,g)

	R = checkError(hOp, tol = tol)

    Kmin = 1 / (1 - R)
    Kmax = 1 / (1 - R - tol)

    return (Kmin, Kmax)

end
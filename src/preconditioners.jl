
# applies the jacobi preconditioner of the laplacian (diagonal entries)
function jacobi{Tv,Ti}(la::SparseMatrixCSC{Tv,Ti})
	
	d = diag(la)

	f = function(b::Array{Float64,1})
		res = copy(b) - mean(b)

		res = res ./ d

		return res - mean(res)
	end

	return f

end

# applies the gauss siedel preconditioner. precond = (D + L) * D^-1 * (D + U)
function gs{Tv,Ti}(la::SparseMatrixCSC{Tv,Ti})

	m1 = LowerTriangular(la)
	m2 = 1 ./ diag(la)
	m3 = UpperTriangular(la)

	f = function(b::Array{Float64,1})
		res = copy(b) - mean(b)

		res = m3 \ res
		res = res ./ m2
		res = m1 \ res

		return res - mean(res)
	end

	return f

end
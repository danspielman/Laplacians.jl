# get the max eigenvalue using the power method
function powerIteration{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti}; tol=1e-3, maxit = Inf, verbose=false)

	n = a.n

	notConverged = true
	b = rand(n)

	iter = 0
	while iter < maxit
		iter = iter + 1

		newb = a * b / norm(a * b)
		b = newb

		lam = norm(a * b) / norm(b)
		if norm(a * b - lam * b) < tol
			notConverged = false

			if verbose
				println("Converged after ", iter, " iterations with error ", norm(a * b - lam * b))
			end

			return norm(a * b) / norm(b)
		end

	end

	if verbose
		lam = norm(a * b) / norm(b)
		println("Stopped after ", iter, " iterations with error ", norm(a * b - lam * b))
	end

	return b

end
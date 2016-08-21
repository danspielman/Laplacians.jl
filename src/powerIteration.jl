# get the max eigenvalue using the power method
function powerIteration{Tv,Ti}(op::SqLinOp{Tv,Ti}; tol=1e-3, maxit = Inf, verbose=false)

	n = op.n

	b = rand(n)
	prevLam = 0
	lam = 0
	err = 0

	iter = 0
	while iter < maxit
		iter = iter + 1

		b = op.multFn(b) / norm(op.multFn(b))

		prevLam = lam
		lam = norm(op.multFn(b)) / norm(b)
		err = abs(prevLam - lam)

		if err < tol && iter > 1
			if verbose
				println("Converged after ", iter, " iterations with error ", err)
			end

			return lam
		end

	end

	if verbose
		println("Failed to converge after ", iter, " iterations with error ", err)
	end

	return b

end
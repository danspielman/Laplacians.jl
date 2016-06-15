# Implements the Johnson-Lindenstauss resistance upperbounding
# TODO: not optimized for speed - let's see how resistance estimates improve the number of nonzeros at the end

function johnlind{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti}, eps::Float64)

	n = a.n
	# m = ceil(Int64, length(a.nzval) / 2)
	m = length(a.nzval)
	dhat = ceil(Int64, 24 * log(m) / eps^2)

	println("dhat = ", dhat)

	P = ones(dhat, m)
	for i in 1:dhat
		for j in 1:m
			if rand() < 1/2
				P[i,j] = -P[i,j]
			end
		end
	end
	# we want ||P * x|| = ||x||
	P = P / sqrt(dhat)

	# compute W, B
	W = speye(m,m) * 0
	B = speye(m,n) * 0

	pos = 1
	ind = 0
	for i in 1:length(a.nzval)
		while a.colptr[pos + 1] <= i
			pos = pos + 1
		end

		p = a.rowval[i]
		q = pos
		w = a.nzval[i]

		if p > q
			aux = p
			p = q
			q = aux
		end

		ind = ind + 1

		W[ind,ind] = w
		B[ind,p] = 1
		B[ind,q] = -1
	end

	# be L+ be = be L+ L L+ be = be L+ B' W B L+ be

	# check if B * W * B' - L = 0
	println("the diff is: ", maximum(B' * W * B / 2 - lap(a)))

	# P * W^(1/2) * B. Solve for each line. dims are dhat x n
	bs = P * sqrt(W) * B / sqrt(2)

	la = lap(a)
	f = lapWrapSolver(augTreeSolver,la,tol=1e-6,maxits=1000)

	xs = zeros(dhat, n)
	for i in 1:dhat
		b = reshape(bs[i,:], n)
		b = b - mean(b)
		xs[i,:] = f(b)
	end

	# let's compute xhat[i] = xs * unitVec(i)
	xhat = zeros(n, dhat)
	for i in 1:n
		unitVec = zeros(n)
		unitVec[i] = 1

		xhat[i,:] = xs * unitVec
	end

	println(size(xhat))

	# compute the effective resistance
	reff = copy(a)
	pos = 1
	for i in 1:m
		while a.colptr[pos + 1] <= i
			pos = pos + 1
		end

		p = a.rowval[i]
		q = pos

		reff.nzval[i] = norm(xhat[p,:] - xhat[q,:])^2
	end

	return reff

end
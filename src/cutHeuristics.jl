
"""
	Modifies a cluster by adding/removing vertices based on Deg_external - Deg_Internal
	Each vertex can be added in/removed only once
	Uses O(M + N) memory
"""
function refineCut{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Array{Int64,1})

	n = max(G.n, G.m)

	# inset - keep track of vertices in the low conductance set
	inset = zeros(Bool, n)
	for u in s
		inset[u] = true
	end

	# pq_in - costs are maintained to add elements to the set
	# pq_out - costs are maintained to remove elements from the set
	pq_in = Collections.PriorityQueue{Int64,Int64,Base.Order.ForwardOrdering}()
	pq_out = Collections.PriorityQueue{Int64,Int64,Base.Order.ReverseOrdering{Base.Order.ForwardOrdering}}(Base.Order.Reverse)

	for u in 1:n
		cost = 0

		for i in 1:deg(G,u)
			v = nbri(G,u,i)
			if !inset[v]
				cost = cost + weighti(G,u,i)
			else
				cost = cost - weighti(G,u,i)
			end
		end

		pq_in[u] = cost
		pq_out[u] = cost
	end

	improve = true
	tries = 0
	while improve && tries < n
		improve = false

		# try to add a new vertex into the set
		u,contrib = Base.Collections.peek(pq_in)
		if !inset[u] && contrib < 0
			inset[u] = true
			improve = true

			for i in 1:deg(G,u)
				v = nbri(G,u,i)
				pq_in[v] = pq_in[v] - weighti(G,u,i)
				pq_out[v] = pq_out[v] - weighti(G,u,i)
			end

			Collections.dequeue!(pq_in)
			pq_out[u] = typemax(Int64)
			Collections.dequeue!(pq_out)
		end

		# try to remove a vertex from the set
		u,contrib = Base.Collections.peek(pq_out)
		if inset[u] && contrib > 0
			inset[u] = false
			improve = true

			for i in 1:deg(G,u)
				v = nbri(G,u,i)
				pq_in[v] = pq_in[v] + weighti(G,u,i)
				pq_out[v] = pq_out[v] + weighti(G,u,i)
			end

			Collections.dequeue!(pq_out)
			pq_in[u] = typemin(Int64)
			Collections.dequeue!(pq_in)
		end
	end

	news = (Int64)[]
	for i in 1:n
		if inset[i]
			push!(news, i)
		end
	end

	return news
end

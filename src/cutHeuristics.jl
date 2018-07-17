
"""
	Modify a cluster by adding or removing vertices by picking at each step 
	the vertex that has the maximum value of (Deg_external - Deg_Internal).
	Each vertex can be added in/removed only once.
"""
function refineCut(G::SparseMatrixCSC{Tv,Ti}, s::Array{Int64,1}) where {Tv,Ti}

	n = max(G.n, G.m)

	# inset - keep track of vertices in the low conductance set
	inset = zeros(Bool, n)
	for u in s
		inset[u] = true
	end

	# pq_in - costs are maintained to add elements to the set
	# pq_out - costs are maintained to remove elements from the set
	pq_in = DataStructures.PriorityQueue(Base.Order.Reverse)
	pq_out = DataStructures.PriorityQueue(Base.Order.Forward)

	for u in 1:n
		cost = 0

		for i in 1:deg(G,u)
			v = nbri(G,u,i)
			if inset[v]
				cost = cost + weighti(G,u,i)
			else
				cost = cost - weighti(G,u,i)
			end
		end

		if inset[u] == false
			pq_in[u] = cost
		end
		if inset[u] == true
			pq_out[u] = cost
		end
	end

	improve = true
	tries = 0
	while improve && tries < n
		improve = false

		# try to add a new vertex into the set
		if !isempty(pq_in)
			u,contrib = DataStructures.peek(pq_in)

			if contrib >= 0
				improve = true
				DataStructures.dequeue!(pq_in)
				inset[u] = true

				for i in 1:deg(G,u)
					v = nbri(G,u,i)

					if v in keys(pq_in)
						pq_in[v] = pq_in[v] + weighti(G,u,i)
					end
					if v in keys(pq_out)
						pq_out[v] = pq_out[v] + weighti(G,u,i)
					end
				end
			end
		end

		# try to remove a vertex from the set
		if !isempty(pq_out)
			u,contrib = DataStructures.peek(pq_out)
			if contrib < 0
				improve = true
				DataStructures.dequeue!(pq_out)
				inset[u] = false

				for i in 1:deg(G,u)
					v = nbri(G,u,i)
					if v in keys(pq_in)
						pq_in[v] = pq_in[v] - weighti(G,u,i)
					end
					if v in keys(pq_out)
						pq_out[v] = pq_out[v] - weighti(G,u,i)
					end
				end
			end
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

"""
	Modify a cluster by passing through all the vertices exactly once and 
	adding/removing them based on the value of (Deg_external - Deg_Internal).
"""
function dumbRefineCut(G::SparseMatrixCSC{Tv,Ti}, s::Array{Int64,1}) where {Tv,Ti}

	n = max(G.n, G.m)

	news = IntSet(s)

	for v in 1:n
		nrbsA = 0
		nrbsnotA = 0
		for i in 1:deg(G, v)
			u = nbri(G, v, i)
			if u in news
				nrbsA = nrbsA + weighti(G, v, i)
			else
				nrbsnotA = nrbsnotA + weighti(G, v, i)
			end
		end

		if nrbsA >= nrbsnotA
			push!(news, v)
		else
			if v in news
				pop!(news, v)
			end
		end
	end

	return collect(news)

end

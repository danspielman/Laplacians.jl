
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
	pq_in = Collections.PriorityQueue{Int64,Int64,Base.Order.ReverseOrdering{Base.Order.ForwardOrdering}}(Base.Order.Reverse) 
	pq_out = Collections.PriorityQueue{Int64,Int64,Base.Order.ForwardOrdering}()

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
			u,contrib = Base.Collections.peek(pq_in)

			if contrib >= 0
				improve = true
				Collections.dequeue!(pq_in)
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
			u,contrib = Base.Collections.peek(pq_out)
			if contrib < 0
				improve = true
				Collections.dequeue!(pq_out)
				inset[u] = false

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
	end

	news = (Int64)[]
	for i in 1:n
		if inset[i]
			push!(news, i)
		end
	end

	return news
end

function dumb{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Array{Int64,1})

	n = max(G.n, G.m)

	news = copy(s)

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
			for i in 1:length(news)
				if news[i] == v
					news[i] = news[length(news)]
					pop!(news)
					break
				end
			end
		end
	end

	return news

end

" the LocalImprove function, from the Orrechia-Zhu paper "
function localImprove{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, A::Array{Int64,1}, epsSigma::Float64, eps::Float64, err=1e-10) 
  #=
    Notes: err < 1e-13 breaks the code. Precision is 1e-16, so 3 bits of precision are lost somewhere...

  =#

  # show warning for epsSigma
  minEpsSigma = getVolume(G, A) / getVolume(G, setdiff(collect(1:max(G.n, G.m)), A))
  if epsSigma < minEpsSigma
      print_with_color(:red, "eps-sigma should be greater than ")
      println(minEpsSigma)
  end

  alphaMin = 0.0
  alphaMax = 1.0
  while alphaMax - err > alphaMin
    alpha = (alphaMin + alphaMax) / 2

    if abs(localFlow(G, A, alpha, epsSigma)[2] - getVolume(G, A)) < err
      alphaMin = alpha
    else
      alphaMax = alpha
    end
  end

  return localFlow(G, A, alphaMax, epsSigma)

end # localImprove

" the LocalFlow function, from the Orecchia-Zhu paper "
function localFlow{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, A::Array{Int64,1}, alpha::Float64, epsSigma::Float64)

  # compute the number of vertices
  n = max(G.n, G.m)

  # compute the number of Iterations for which we'll do flow
  # as this number can become really big, we'll upper bound it by 2^32
  volA = getVolume(G, A)
  I = 5 / alpha * log(3 * volA / alpha)

  B = (Ti)[]
  totalflow = 0
  GPrime = createGPrime(G, A, alpha, epsSigma)
  s = max(G.n, G.m) + 1
  t = max(G.n, G.m) + 2

  #=
    "considered" will keep track of the vertices interesting for the flow graph.
    Initially: {A} + {Neighbors(A)}
  =#
  considered = zeros(Bool, n + 2)
  considered[s] = 1
  considered[t] = 1
  for i in 1:length(A)
    considered[A[i]] = 1
  end
  for i in 1:length(A)
    u = A[i]
    for j in 1:deg(GPrime, u)
      v = nbri(GPrime, u, j)
      considered[v] = 1
    end
  end

  # do the iterative part of the algorithm: lines 4-9
  backInd = backIndices(GPrime)
  for iter in 1:floor(Int64, I)
    saturated, flowInc = localBlockFlow(GPrime, backInd, s, t, considered)
    if flowInc == 0
      break
    end

    # update the total flow
    totalflow = totalflow + flowInc

    # update the vertices consiedered for the flow
    for u in saturated
      for i in 1:deg(GPrime, u)
        v = nbri(GPrime, u, i)
        considered[v] = 1
      end
    end
  end

  return getCutSet(GPrime, s, t), totalflow

end # localFlow

" computes the volume of subset A in G "
function getVolume{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, A::Array{Int64,1})

  vol = 0
  for i in 1:length(A)
    vol = vol + deg(G, A[i])
  end

  return vol

end # getVolume

" creates the GPrime graph, following Definition 3.1 "
function createGPrime{Tv,Ti}(G::SparseMatrixCSC{Tv, Ti}, A::Array{Int64, 1}, alpha::Float64, epsSigma::Float64)

  u,v,w = findEntries(G)

  # set costs of edges to 0
  for i in 1:length(w)
    w[i] = 1 / alpha
  end

  # we'll consider the source to be n+1, and the sink to be n+2
  n = max(G.n, G.m)
  s = n + 1
  t = n + 2

  # add edges from s to A with weight deg()
  inA = zeros(Int64, n)
  for i in 1:length(A)
    inA[A[i]] = 1

    push!(u, s)
    push!(v, A[i])
    push!(w, deg(G, A[i]))

    push!(u, A[i])
    push!(v, s)
    push!(w, deg(G, A[i]))
  end

  # add edges from V-A to t with weight epsSigma * deg()
  for i in 1:n
    if inA[i] == 0
      push!(u, t)
      push!(v, i)
      push!(w, epsSigma * deg(G, i))

      push!(u, i)
      push!(v, t)
      push!(w, epsSigma * deg(G, i))
    end
  end

  #=
    Now that we've created the graph, make edges from s to the graph and from
    the graph to t directed. We're doing this now, because we need edges of
    weight 0 for the flow, and if we use sparse() with weight 0 the edges will
    get ignored.
  =#
  GPrime = sparse(u, v, w)

  # make edges from s directed
  for i in 1:length(A)
    u = A[i]
    for j in 1:deg(GPrime, u)
      v = nbri(GPrime, u, j)
      if v == s
        setValue(GPrime, u, j, convert(Tv, 0))
      end
    end
  end

  # make edges to t directed
  for i in 1:deg(GPrime, t)
    setValue(GPrime, t, i, convert(Tv, 0))
  end

  return GPrime

end # createGPrime

" computes block flow between s and t only on nodes in the set <considered> "
function localBlockFlow{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, backInd::SparseMatrixCSC{Ti,Ti}, s::Int64, t::Int64, considered::Array{Bool,1})

  # Note: this function is overall very similar to the flow function 

  n = max(G.n, G.m)
  inQ = zeros(Bool, n)
  prev = [(0,0) for i in 1:n]
  Q = zeros(Ti, n)

  # do one iteration of block-flow
  left = 0
  right = 1
  Q[right] = s
  inQ[s] = 1

  while left < right
    left = left + 1
    u = Q[left]
    if (u == t)
      continue
    end

    for i in 1:deg(G, u)
      v = nbri(G, u, i)

      if considered[v] == true && inQ[v] == false && weighti(G, u, i) > 0
        right = right + 1
        Q[right] = v
        inQ[v] = true
        prev[v] = (u, i)
      end
    end
  end

  if inQ[t] == false
    return [], 0
  else
    totalflow = 0
    saturated = []

    for i in 1:deg(G, t)
      # ignore the current vertex if it's not in considered
      if considered[nbri(G, t, i)] == false || inQ[nbri(G, t, i)] == false
        continue
      end

      # see if there is any flow going from u to t
      prev[t] = (nbri(G, t, i), weighti(backInd, t, i))
      u = t

      currentflow = typemax(Tv)
      while u != s
        prevU = prev[u][1]
        indexPrevU = prev[u][2]
        currentflow = min(currentflow, weighti(G, prevU, indexPrevU))
        u = prevU
      end

      if currentflow != 0
        # do flow and update edges
        u = t
        while u != s
          prevU = prev[u][1]
          indexPrevU = prev[u][2]
          indexU = weighti(backInd, prevU, indexPrevU)

          capPrevUU = weighti(G, prevU, indexPrevU)
          capUPrevU = weighti(G, u, indexU)

          setValue(G, prevU, indexPrevU, capPrevUU - currentflow)
          setValue(G, u, indexU, capUPrevU + currentflow)

          u = prevU
        end
      end

      # check if the vertex to the sink is saturated
      if weighti(G, nbri(G, t, i), weighti(backInd, t, i)) == 0
        push!(saturated, nbri(G, t, i))
      end

      totalflow = totalflow + currentflow
    end

    return saturated, totalflow
  end

end # localBlockFlow

" get the min cut from the source - return all vertices in the cut besides the source "
function getCutSet{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Int64, t::Int64)

  n = max(G.n, G.m)
  Q = zeros(Int64, n)
  inQ = zeros(Bool, n)

  left = 0
  right = 1
  Q[right] = s
  inQ[s] = true

  while left < right
    left = left + 1
    u = Q[left]

    for i in 1:deg(G, u)
      v = nbri(G, u, i)

      if inQ[v] == 0 && weighti(G, u, i) > 0
        right = right + 1
        Q[right] = v
        inQ[v] = 1
      end
    end
  end

  # consider only vertices in the initial graph into the getCutSet
  cut = (Int64)[]
  for i in 1:right
    if Q[i] != s && Q[i] != t
      push!(cut, Q[i])
    end
  end

  return cut

end # getCutSet

""" 
  Returns the quality of the cut for a given *unweighted* graph and a given cut set.
  Result will be |outgoing edges| / min(|vertices in set|, |N - vertices in set|)
"""
function compConductance{Tv, Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Array{Int64,1})

	n = max(G.n, G.m)

	volIn = getVolume(G, s)
	volOut = getVolume(G, setdiff(collect(1:n), s))

	ins = zeros(Bool, n)
	for u in s
		ins[u] = true
	end

	obound = 0
	for u in s
		for i in 1:deg(G, u)
			v = nbri(G, u, i)
			if !ins[v]
				obound = obound + 1
			end
		end
	end

	return obound / min(volIn, volOut)

end # compConductance



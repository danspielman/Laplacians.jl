""" 
  localImprove{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, A::Array{Int64,1}; epsSigma=-1.0, err=1e-10, maxSize = max(G.n, G.m)

The LocalImprove function, from the Orrechia-Zhu paper. Given a graph and an initial set, finds a set of smaller conductance
based on the starting set using a localized version of max-flow. 

G is the given graph, A is the initial set 
epsSigma is a measure of the quality of the returning set (the smaller the better). It's defaulted to volume(A) / volume(V\A)
err is the numerical error considered throughout the algorithm. It's defaulted to 1e-10
maxSize is the maximum allowed size for the flow graph at any iteration of the algorithm. It's defaulted to |V|
"""
function localImprove{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, A::Array{Int64,1}; epsSigma=-1.0, err=1e-10, maxSize = max(G.n, G.m)) 
  #=
    Notes: err smaller than 1e-13 might give weird behavior
  =#

  n = max(G.n, G.m)

  # set epsSigma if it was not set
  minEpsSigma = getVolume(G, A) / getVolume(G, setdiff(collect(1:n), A))
  if epsSigma == -1
    epsSigma = minEpsSigma
  end

  # show warnings for epsSigma
  if epsSigma < minEpsSigma
    print_with_color(:red, "eps-sigma should be greater than ")
    println(minEpsSigma)
  end

  # show warnings for maxSize
  nbrs = Set{Int64}()
  for i in 1:length(A)
    push!(nbrs, A[i])
    for j in 1:deg(G,A[i])
      push!(nbrs, nbri(G,A[i],j))
    end
  end
  if length(nbrs) > maxSize
    print_with_color(:red, "maxSize should be at least ", length(nbrs), " .....finishing execution")
    return [],0
  end

  # the actual algorithm
  alphaMin = 0.0
  alphaMax = 1.0
  while alphaMax - err > alphaMin
    alpha = (alphaMin + alphaMax) / 2

    if abs(localFlow(G, A, alpha, epsSigma, maxSize)[2] - getVolume(G, A)) < err
      alphaMin = alpha
    else
      alphaMax = alpha
    end
  end

  return localFlow(G, A, alphaMax, epsSigma, maxSize)

end # localImprove


" The LocalFlow function, from the Orecchia-Zhu paper "
function localFlow{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, A::Array{Int64,1}, alpha::Float64, epsSigma::Float64, maxSize = max(G.n, G.m))

  # compute the number of vertices
  n = max(G.n, G.m)

  # compute the number of Iterations for which we'll do flow
  volA = getVolume(G, A)
  I = 5 / alpha * log(3 * volA / alpha)

  #=
    Initially: {A} + {Neighbors(A)}
    newID will represent the index in GPrime of a certain vertex in G
    oldID will represent the index in G from a vertex in GPrime
  =#

  newID = Dict{Int64,Int64}()
  oldID = Dict{Int64,Int64}()
  s = maxSize + 1
  t = maxSize + 2

  # initialize GPrime with the input set, A. The current index for new vertices will be "current"
  GPrime = initGPrime(G, A, newID, oldID, alpha, maxSize)
  current = length(A)

  # "considered" will keep track of the vertices interesting for the flow graph.
  considered = Set{Int64}()
  for i in 1:length(A)
    push!(considered, A[i])
  end
  for i in 1:length(A)
    u = A[i]
    for j in 1:deg(G,u)
      v = nbri(G,u,j)

      if v in considered == false
        push!(considered, v)

        # add this new vertex to GPrime
        current = current + 1

        addToGPrime(G, GPrime, newID, oldID, v, current, alpha, epsSigma, maxSize)
      end
    end
  end

  # do the iterative part of the algorithm: lines 4-9
  totalflow = 0
  for iter in 1:floor(Int64, I)
    saturated, flowInc = localBlockFlow(GPrime, s, t)
    saturated = [oldID[u] for u in saturated]
    if flowInc == 0
      break
    end

    # update the total flow
    totalflow = totalflow + flowInc

    # update the vertices considered for the flow
    newbatch = Set{Int64}()
    for u in saturated
      for i in 1:deg(G, u)
        v = nbri(G, u, i)
        
        if v in considered == false
          push!(considered, v)
          push!(newbatch, v)
        end
      end
    end

    # see if the size of the graph exceeds the maximum size imposed by the user
    if (length(considered) > maxSize)
      return [oldID[u] for u in getCutSet(GPrime, s, t)], getVolume(G, A)
    end

    # update GPrime
    for u in newbatch
      current = current + 1
      addToGPrime(G, GPrime, newID, oldID, u, current, alpha, epsSigma, maxSize)
    end

  end

  return [oldID[u] for u in getCutSet(GPrime, s, t)], totalflow

end # localFlow


" Initialize GPrime with the set A and edges of type s->u"
function initGPrime{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, A::Array{Int64,1}, newID::Dict{Int64,Int64}, oldID::Dict{Int64,Int64}, alpha::Float64, maxSize::Int64)
  GPrime = [Tuple{Int64,Float64}[] for i in 1:(maxSize + 2)]
  s = maxSize + 1

  for i in 1:length(A)
    newID[A[i]] = i
    oldID[i] = A[i]

    # add a directed edges from s to i
    push!(GPrime[s], (newID[A[i]], wdeg(G, A[i])))
    push!(GPrime[newID[A[i]]], (s, 0.0))
  end

  # add edges between the vertices in A
  for u in A
    for i in 1:deg(G,u)
      v = nbri(G,u,i)
      if haskey(newID, v)
        push!(GPrime[newID[u]], (newID[v], weighti(G,u,i) / alpha))
      end
    end
  end

  return GPrime
end # initGPrime


" Add a new vertex to GPrime "
function addToGPrime{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, GPrime::Array{Array{Tuple{Int64,Float64},1},1}, 
                            newID::Dict{Int64,Int64}, oldID::Dict{Int64,Int64},
                            u::Int64, newU::Int64, alpha::Float64, epsSigma::Float64, maxSize::Int64)
  t = maxSize + 2

  # 1. store indices for the new vertex
  newID[u] = newU
  oldID[newU] = u

  # 2. add an edge from newU to t (directed edge)
  push!(GPrime[newU], (t, epsSigma * wdeg(G, u)))
  push!(GPrime[t], (newU, 0))

  # 3. add weights between newU and his neighbors already evaluated
  for i in 1:deg(G, u)
    v = nbri(G, u, i)
    if haskey(newID, v)
      newV = newID[v]
      push!(GPrime[newU], (newV, weighti(G, u, i) / alpha))
      push!(GPrime[newV], (newU, weighti(G, u, i) / alpha))
    end
  end
end


" Compute block flow between s and t"
function localBlockFlow(G::Array{Array{Tuple{Int64,Float64},1},1}, s::Int64, t::Int64)

  # backInd = backIndices(G)

  n = length(G)
  inQ = zeros(Bool, n)
  prev = [(0,0) for i in 1:n]
  Q = zeros(Int64, n)

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

    for i in 1:length(G[u])
      v,w = G[u][i]

      if inQ[v] == false && w > 0
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
    for i in 1:length(G[t])
      u,w = G[t][i]
      if inQ[u] == false
        continue
      end

      # see the index of t in u's list
      indexT = 1
      for j in 1:length(G[u])
        if G[u][j][1] == t
          indexT = j
          break
        end
      end

      # indexT = backInd[t][i]
      prev[t] = (u, indexT)
      u = t

      # compute the augmentation flow
      currentflow = typemax(Float64)
      while u != s
        prevU = prev[u][1]
        indexPrevU = prev[u][2]
        currentflow = min(currentflow, G[prevU][indexPrevU][2])
        u = prevU
      end

      if currentflow != 0
        # update edges with augmentation flow
        u = t
        while u != s
          prevU = prev[u][1]
          indexPrevU = prev[u][2]
          # indexU = backInd[prevU][indexPrevU]

          indexU = 1
          for j in 1:length(u)
            if G[u][j][1] == prevU
              indexU = j
              break
            end
          end

          # between prevU and u
          G[prevU][indexPrevU] = (G[prevU][indexPrevU][1], G[prevU][indexPrevU][2] - currentflow)

          # between u and prevU
          G[u][indexU] = (G[u][indexU][1], G[u][indexU][2] + currentflow)

          u = prevU
        end
      end

      # check if the vertex to the sink is saturated
      u = G[t][i][1]
      if G[u][indexT][2] == 0
        push!(saturated, u)
      end

      totalflow = totalflow + currentflow
    end

    return saturated, totalflow
  end

end # localBlockFlow


" Get the min cut from the source - return all vertices in the cut besides the source "
function getCutSet(G::Array{Array{Tuple{Int64,Float64},1},1}, s::Int64, t::Int64)

  n = length(G)
  Q = zeros(Int64, n)
  inQ = zeros(Bool, n)

  left = 0
  right = 1
  Q[right] = s
  inQ[s] = true

  while left < right
    left = left + 1
    u = Q[left]

    for i in 1:length(G[u])
      v,w = G[u][i]

      if inQ[v] == 0 && w > 0
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
  prn{Tv, Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Array{Int64,1}, phi::Float64, b::Int64)

The PageRank-Nibble cutting algorithm from the Anderson/Chung/Lang paper\n
s is a set of starting vertices, phi is a constant in (0, 1], and b is an integer in [1, [log m]]

phi is a bound on the quality of the conductance of the cut - the smaller the phi, the higher the quality. 
b is used to handle precision throughout the algorithm - the higher the b, the greater the precision.
"""
function prn{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Array{Int64,1}, phi::Float64, b::Int64)

  m = div(nnz(G), 2)

  # show warnings
  if (!(0 < phi && phi <= 1))
    print_with_color(:red, "phi should be in (0, 1]\n")
  end
  if (!(1 <= b && b <= round(Int64, log(m))))
    print_with_color(:red, string("b should be in [1, ", string(round(Int64, log(m))), "]\n"))
  end

  alpha = phi * phi / (225 * log(100 * sqrt(m)))
  eps = 1 / (2^b * 48 * log(m))

  p = apr(G, s, alpha, eps)

  # iterate through support set and find a cut-set S that fits our constraints
  S = Set(Int64[])
  volS = 0
  volSc = 2 * m
  obound = 0

  while !isempty(p)
    u = Collections.dequeue!(p)
    push!(S, u)

    # update volumes and number of connecting edges
    volS = volS + deg(G, u)
    volSc = volSc - deg(G, u)
    for j in 1:deg(G, u)
      v = nbri(G, u, j)

      if v in S
        obound = obound - 1
      else
        obound = obound + 1
      end
    end

    # next are ignoring conditions

    # ignore from conductance
    conductance = obound / (min(volS, volSc))
    if (conductance >= phi)
      continue
    end

    # ignore from set size
    if !(2^b < volS && volS < 2.0/3.0 * (2 * m))
      continue
    end

    return collect(S)
  end

  return Int64[]

end # prn_local

""" 
Computes an approximate page rank vector from a starting set s, an alpha and an epsilon
The algorithm follows the Anderson,Chung,Lang paper and Dan Spielman's lecture notes
"""
function apr{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Array{Int64,1}, alpha::Float64, eps::Float64)

  # p is ordered by p[i] / deg(i) (ordering required for prn)
  # r is oredered by r[i] - eps * deg(G,u). r is initially mimicking the unit vector
  p = Collections.PriorityQueue{Int64,Float64,Base.Order.ReverseOrdering{Base.Order.ForwardOrdering}}(Base.Order.Reverse) 
  pelems = Set(s)

  r = Collections.PriorityQueue{Int64,Float64,Base.Order.ReverseOrdering{Base.Order.ForwardOrdering}}(Base.Order.Reverse) 
  relems = Set(s)

  for u in s
    p[u] = 0
    r[u] = 1 / length(s)^2 - eps * deg(G,u)
  end

  # we are moving mass from a node u only if it has more mass than eps * deg(G,u)
  # making the inequality strictly greater than 0 deals with u being an isolated node
  while Base.Collections.peek(r)[2] > 0

    u,ru = Base.Collections.peek(r)
    ru = ru + eps * deg(G,u)

    # check if u is in the priority queue for p
    if u in pelems == false
      p[u] = 0
      push!(pelems, u)
    end

    # update p[u] & r[u]
    p[u] = p[u] + alpha * ru / deg(G,u)
    r[u] = -eps * deg(G,u) # this means it's set to 0

    for i in 1:deg(G,u)
      v = nbri(G,u,i)

      # check if v is in the priority queue for r
      if v in relems == false
        r[v] = -eps * deg(G,v) # this means it's set to 0
        push!(relems, v)
      end

      # update u's neighbors in r
      r[v] = r[v] + (1 - alpha) * ru / deg(G, u)
    end

  end

  return p

end # apr_local
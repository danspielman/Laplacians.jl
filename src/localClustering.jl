"""
  prn_local{Tv, Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Array{Int64,1}, phi::Float64, b::Int64)

  the PageRank-Nibble cutting algorithm from the Anderson/Chung/Lang paper\n
  s is a set of starting vertices, phi is a constant in (0, 1], and b is an integer in [1, [log m]]

  phi is a bound on the quality of the conductance of the cut - the smaller the phi, the higher the quality
  b is used to handle precision throughout the algorithm - the higher the b, the smaller the eps
"""
function prn_local{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Array{Int64,1}, phi::Float64, b::Int64)

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

  p = apr_local(G, s, alpha, eps)

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
  computes an approximate page rank vector from a starting set s, an alpha and an epsilon
  algorithm follows the Anderson,Chung,Lang paper and Dan Spielman's notes
"""
function apr_local{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Array{Int64,1}, alpha::Float64, eps::Float64)

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
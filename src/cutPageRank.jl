"""
  prn{Tv, Ti}(G::SparseMatrixCSC{Tv,Ti}, v::Array{Int64,1}, phi::Float64, b::Int64)

  the PageRank-Nibble cutting algorithm from the Anderson/Chung/Lang paper\n
  v is a set of starting vertices, phi is a constant in (0, 1], and b is an integer in [1, [log m]]

  phi is a bound on the quality of the conductance of the cut - the smaller the p, the higher the quality
  b is used to handle precision throughout the algorithm - the higher the b, the smaller the eps
"""
function prn{Tv, Ti}(G::SparseMatrixCSC{Tv,Ti}, v::Array{Int64,1}, phi::Float64, b::Int64)

  g = findnz(G)
  m = length(g[1]) / 2

  # show warnings
  if (!(0 < phi && phi <= 1))
    print_with_color(:red, "phi should be in (0, 1]\n")
  end
  if (!(1 <= b && b <= round(Int64, log(m))))
    print_with_color(:red, string("b should be in [1, ", string(round(Int64, log(m))), "]\n"))
  end

  alpha = phi * phi / (225 * log(100 * sqrt(m)))
  eps = 1 / (2^b * 48 * log(m))

  p = apr(G, toUnitVector(v, G.n), alpha, eps)

  # order p by p[i] / deg[i]
  supp = Tuple{Float64, Int64}[]
  for i in 1:G.n
    if (isequal(0, p[i]) == false)
      push!(supp, (p[i] / deg(G, i), i))
    end
  end
  supp = sort(supp)

  # iterate through support set and find a cut-set S that fits our constraints
  S = Int64[]
  inS = zeros(G.n)
  volS = 0
  volSc = 2 * m
  obound = 0

  for j in length(supp):-1:1
    u = supp[j][2]
    push!(S, supp[j][2])
    inS[u] = 1

    # update volumes and number of connecting edges
    volS = volS + deg(G, u)
    volSc = volSc - deg(G, u)
    for i in 1:deg(G, u)
      v = nbri(G, u, i)
      if inS[v] == 0
        obound = obound + 1
      else
        obound = obound - 1
      end
    end

    conductance = obound / (min(volS, volSc))
    if (conductance >= phi)
      continue
    end

    # next up - the constraint on the volume of S
    if !(2^b < volS && volS < 2.0/3.0 * (2 * m))
      continue
    end

    return S
  end

  return Int64[]

end # prnibble

" computes an approximate page rank vector from a starting vector s, an alpha and an epsilon "
function apr{Tv, Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Array{Float64,1}, alpha::Float64, eps::Float64)

  p = zeros(length(s))
  r = copy(s)

  Q = Queue(Int64)
  inQ = zeros(G.n)
  for u in 1:G.n
    if r[u] >= eps * deg(G, u)
      enqueue!(Q, u)
    end
  end

  while isempty(Q) == false
    u = dequeue!(Q)
    inQ[u] = 0

    pushpr(G, u, alpha, p, r)
    for i in 1:deg(G, u)
      v = nbri(G, u, i)

      if inQ[v] == 0 && r[v] >= eps * deg(G, v)
        enqueue!(Q, v)
        inQ[v] = 1
      end
    end

    # check if u is still greater than eps * deg(G, u)
    if r[u] >= eps * deg(G, u)
      enqueue!(Q, u)
      inQ[u] = 1
    end
  end

  return p

end # apr

" implements the push operation for the aproximate page rank vector "
function pushpr{Tv, Ti}(G::SparseMatrixCSC{Tv,Ti}, u::Int64, alpha::Float64, p::Array{Float64,1}, r::Array{Float64,1})

  pu = p[u]
  ru = r[u]

  p[u] = pu + alpha * ru
  r[u] = 0

  for i in 1:deg(G, u)
    v = nbri(G, u, i)
    r[v] = r[v] + (1 - alpha) * ru / deg(G, u)
  end

end # pushpr

" computes a page rank vector satisfying p = a/n * 1 + (1 - a) * W * p "
function pr{Tv, Ti}(G::SparseMatrixCSC{Tv,Ti}, alpha::Float64)
  # compute the walk matrix
  u,v,w = findnz(G)
  for i in 1:length(u)
    if (u[i] == v[i])
      w[i] = 0
    else
      w[i] = w[i] / deg(G, v[i])
    end
  end
  W = sparse(u, v, w)

  p = \(W - eye(G.n) / (1 - alpha), -alpha / (1 - alpha) / G.n * ones(G.n))

end # pr

" computes the personal page rank vector from a starting vector s and an alpha; operates with lazy walk matrix "
function ppr{Tv, Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Array{Float64,1}, alpha::Float64)

  # LW = 1/2 * (I + A * D^-1)
  LW = 0.5 * (speye(G.n) + G * inv(full(diagmat(G))))

  # p = inv(I - (1 - alpha) * W) * alpha * s
  p = \(((speye(G.n) - (1 - alpha) * LW) / alpha), s)

  return p

end # ppr

function ppr{Tv, Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Array{Float64,1}, alpha::Float64, niter::Int64)

  # LW = 1/2 * (I + A * D^-1)
  LW = 0.5 * (speye(G.n) + G * inv(full(diagmat(G))))

  # s is the vector of dried paint, r is the vector of wet paint
  r = copy(s)
  s = zeros(G.n)

  for iter in 1:niter
    news = s + alpha * r
    newr = (1 - alpha) * LW * r

    s = copy(news)
    r = copy(newr)
  end

  return s

end # ppr


# lex:
#
# contains algorithms in the lex-minimizer paper
# includes iterative lex algorithm
#
# Started by xiao.shi@yale.edu on Sep 22, 2015
# Contributers:

import Laplacians.intHeap
import Laplacians.intHeapAdd!, Laplacians.intHeapPop!

LEX_DEBUG = false # debug flag for the lex module
LEX_EPS = 1e-12 # default error tolerance

function setLexDebugFlag(f::Bool)
  global LEX_DEBUG = f
end

#=
Simulate IterLex on Uniformly Weighted Graphs

numIter: number of iterations;

A: n by n adjacency matrix of the graph;

isTerm: boolean array[n] of whether each vertex is a terminal;

initVal: array[n] of initial voltage assignments;

=#
function simIterLexUnwtd{Tv, Ti}(numIter::Int64,
                                 A::SparseMatrixCSC{Tv, Ti},
                                 isTerm::Array{Bool, 1},
                                 initVal::Array{Float64, 1}, )
  n = A.n
  val = copy(initVal)
  nextVal = zeros(Float64, n)

  for t = 1:numIter
    # if the bits representation of val and nextVal are the same for
    # ever vertex, then there is no point in keeping iterating.
    progress = false
    for u = 1:n
      if (!isTerm[u])
        nbrs = A.rowval[A.colptr[u]:(A.colptr[u + 1] - 1)]
        maxNeighbor = maximum(val[nbrs])
        minNeighbor = minimum(val[nbrs])
        nextVal[u] = minNeighbor + (maxNeighbor - minNeighbor) / 2.0
        if (bits(val[u]) != bits(nextVal[u]))
          progress = true
        end
      else
        nextVal[u] = val[u]
      end
    end

    if (!progress)
      @printf("INFO: simIterLexUnwtd: terminating early after %d iterations,
              as numerical error prevents further progress.\n", t)
      return val
    end

    tmp = val
    val = nextVal
    nextVal = tmp
  end
  return val
end

function simIterLex{Tv<:Float64, Ti}(numIter::Int64,
                                     A::SparseMatrixCSC{Tv, Ti},
                                     isTerm::Array{Bool, 1},
                                     initVal::Array{Float64, 1}, )
  n = A.n
  val = copy(initVal)
  nextVal = zeros(Float64, n)

  for t = 1:numIter
    # if the bits representation of val and nextVal are the same for
    # ever vertex, then there is no point in keeping iterating.
    progress = false
    for u = 1:n
      if (!isTerm[u])
        nbrs = A.rowval[A.colptr[u]:(A.colptr[u + 1] - 1)]
        w = A.nzval[A.colptr[u]:(A.colptr[u + 1] - 1)] # weights
        deg = length(nbrs)
        maxGrad = -Inf
        maxI = 0
        maxJ = 0
        for i in 1:deg
          for j in 1:deg
            grad = (val[nbrs[i]] - val[nbrs[j]]) / (1/w[i] + 1/w[j])
            if (grad > maxGrad)
              maxGrad = grad
              maxI = i
              maxJ = j
            end
          end
        end

        wi = w[maxI]
        vi = val[nbrs[maxI]]
        wj = w[maxJ]
        vj = val[nbrs[maxJ]]
        nextVal[u] = (wi * vi + wj * vj) / (wi + wj)
        if (bits(val[u]) != bits(nextVal[u]))
          progress = true
        end
      else
        nextVal[u] = val[u]
      end
    end

    if (!progress)
      @printf("INFO: simIterLexUnwtd: terminating early after %d iterations,
              as numerical error prevents further progress.\n", t)
      return val
    end

    tmp = val
    val = nextVal
    nextVal = tmp
  end
  return val
end

#=
check the correctness of lex assignment for uniformly weighted graphs

lex: the assignment to check;

eps: absolute error tolerance;

fatal: if true, throws error and halt; if false, return false;

=#
function checkLexUnwtd{Tv, Ti}(A::SparseMatrixCSC{Tv, Ti},
                               isTerm::Array{Bool, 1},
                               initVal::Array{Float64, 1},
                               lex::Array{Float64, 1};
                               eps::Float64 = LEX_EPS,
                               fatal::Bool = true)
  n = A.n
  for i in 1:n
    if (isTerm[i])
      expected = initVal[i]
    else
      nbrs = A.rowval[A.colptr[i]:A.colptr[i+1]-1]
      maxval = maximum(lex[nbrs])
      minval = minimum(lex[nbrs])
      expected = minval + (maxval - minval) / 2.0
    end
    if (abs(lex[i] - expected) > eps)
      if (fatal)
        error(@sprintf("The average of the max and min neighbors is %f
                       at vertex %d, but got %f\n", expected, i, lex[i]))
      end
      return false
    end
  end
  return true
end

#=
check the correctness of lex assignment for uniformly weighted graphs

lex: the assignment to check;

eps: absolute error tolerance;

fatal: if true, throws error and halt; if false, return false;

=#
function checkLex{Tv, Ti}(A::SparseMatrixCSC{Tv, Ti},
                          isTerm::Array{Bool, 1},
                          initVal::Array{Float64, 1},
                          lex::Array{Float64, 1};
                          eps::Float64 = LEX_EPS,
                          fatal::Bool = true)
  n = A.n

  for i in 1:n
    correct = true
    if (isTerm[i])
      correct = (lex[i] == initVal[i])
    else
      nbrs = A.rowval[A.colptr[i]:A.colptr[i+1]-1]
      w = A.nzval[A.colptr[i]:A.colptr[i+1]-1] # weights
      deg = length(nbrs)
      grads = zeros(deg)
      for j in 1:deg
        grads[j] = (lex[nbrs[j]] - lex[i]) * w[j]
      end

      maxgrad = maximum(grads)
      mingrad = minimum(grads)

      if (maxgrad + mingrad > eps)
        correct = false
      end
    end
    if (!correct)
      if (fatal)
        error(@sprintf("The average of the max and min neighbors is %f
                       at vertex %d, but got %f\n", expected, i, lex[i]))
      end
      return false
    end
  end
  return true
end

#==========
 LEX TESTS
===========#

function simIterLexTest()
  simIterLexTestPath(4)
  simIterLexTestPath(40)
end

function simIterLexTestPath(n::Int64)
  if (n < 3)
    error("n must be at least 3.")
  end

  numIter = 1000 * n * n

  # path graph with n verticies
  Pn = pathGraph(n)
  # terminal assignments
  isTerm = zeros(Bool, n)
  isTerm[1] = true
  isTerm[n] = true

  initVal = zeros(Float64, n)
  initVal[n] = n - 1

  val = simIterLexUnwtd(numIter, Pn, isTerm, initVal)
  # correct = checkLexUnwtd(Pn, isTerm, initVal, val, fatal=false)

  correct = true
  for i in 1:n
    expected = (i - 1)
    if (abs(val[i] - expected) > LEX_EPS)
      println(val[i], expected)
      println(val)
      error("simIterLex is not correct.")
    end
  end
end

# returns the maximum gradient of an edge in an assignment;
# for ease of verifying CompInfMin
function MaxEdgeGrad(A, val)
  maxGrad = -Inf
  maxX = 0
  maxY = 0
  for x in 1:A.n
    for i in A.colptr[x]:A.colptr[x+1]-1
      y = A.rowval[i]
      grad = abs(val[x]-val[y]) / A[x,y]
      if (grad > maxGrad)
        maxGrad = grad
        maxX = x
        maxY = y
      end
    end
  end
  return maxX, maxY, maxGrad
end


#===============
 Lex Algorithms
===============#

# Shortest Paths _without_ going past terminals.
# Treat weights as reciprocals of distances.
function termFreeShortestPaths{Tv,Ti}(mat::SparseMatrixCSC{Tv, Ti},
                                      start::Ti,
                                      isTerm = zeros(Bool, mat.n), )
  n = mat.n
  visited = zeros(Bool,n)

  nh = intHeap(n)
  dists = nh.keys

  pArray = [1:n...]

  intHeapAdd!(nh, start, 0.0)
  pArray[start] = start

  while nh.nitems > 0
    v::Ti = intHeapPop!(nh)
    visited[v] = true

    dv = dists[v]
    if (isTerm[v] && v != start)
      continue
    end
    for ind in mat.colptr[v]:(mat.colptr[v+1]-1)
      nbr = mat.rowval[ind]
      if !visited[nbr]
        newdist = dv + 1.0 / mat.nzval[ind]
        if newdist < dists[nbr]
          dists[nbr] = newdist
          intHeapAdd!(nh,nbr,newdist)
          pArray[nbr] = v
        end # if
      end # if
    end # for

  end # while
  return copy(dists), pArray
end # termFreeShortestPaths

function ModDijkstra{Tv, Ti}(G::SparseMatrixCSC{Tv, Ti},
                             isTerm::Array{Bool, 1},
                             initVal::Array{Float64, 1},
                             alpha::Float64, )
  n = G.n
  finished = zeros(Bool, n)

  nh = intHeap(n)
  key = nh.keys
  val = copy(initVal)

  parent = [1:n...]
  for i = 1:n
    if (isTerm[i])
      key[i] = initVal[i]
    end
    intHeapAdd!(nh, i, key[i]);
  end

  while (nh.nitems > 0)
    u = intHeapPop!(nh)
    val[u] = key[u]
    finished[u] = true

    for ind in G.colptr[u]:(G.colptr[u+1]-1)
      nbr = G.rowval[ind]
      if (G.nzval[ind] > 0 && !finished[nbr] && !isTerm[nbr])
        newkey = val[u] + alpha / G.nzval[ind]
        if (key[nbr] > newkey)
          intHeapAdd!(nh, nbr, newkey)
          parent[nbr] = u
        end
      end
    end
  end

  return val, parent
end

function CompVLow{Tv, Ti}(G::SparseMatrixCSC{Tv, Ti},
                          isTerm::Array{Bool, 1},
                          initVal::Array{Float64, 1},
                          alpha::Float64)
  return ModDijkstra(G, isTerm, initVal, alpha)
end

function CompVHigh{Tv, Ti}(G::SparseMatrixCSC{Tv, Ti},
                           isTerm::Array{Bool, 1},
                           initVal::Array{Float64, 1},
                           alpha::Float64)
  newInitVal = -copy(initVal)
  temp, parent = ModDijkstra(G, isTerm, newInitVal, alpha)
  return -temp, parent
end

# returns a vertex set (of non terminals)
function CompHighPressGraph(G, isTerm, initVal, alpha)
  vLow, LParent = CompVLow(G, isTerm, initVal, alpha)
  vHigh, HParent = CompVHigh(G, isTerm, initVal, alpha)

  find(vHigh .> vLow * (1 + LEX_EPS))
end


# visited and val will change!
function fixPath!(A, isTerm, val, path)
  l = length(path)
  if (l <= 2)
    return
  end
  if (!isTerm[path[1]])
    error(@sprintf("Value of vertex %d is not set!", path[1]))
  end
  if (!isTerm[path[l]])
    error(@sprintf("Value of vertex %d is not set!", path[l]))
  end

  dist = Array(Float64, l-1)
  cdist = Array(Float64, l-1)
  for i in 1:l-1
    if (A[path[i], path[i+1]] == 0)
      error(@sprintf("No edge between vertices %d and %d!", path[i], path[i+1]))
    end
    dist[i] = A[path[i], path[i+1]]
  end
  cdist = cumsum!(cdist, dist)
  sum = cdist[l-1]

  v1 = val[path[1]]
  v2 = val[path[l]]
  for i in 1:l-2
    d1 = cdist[i]
    d2 = sum - cdist[i]
    val[path[i + 1]] = (d1 * v2 + d2 * v1) / (d1 + d2)
    isTerm[path[i + 1]] = true
  end
end

# check the correctness of lex assignment
function CheckLex(A, isTerm, initVal, lex, eps = 1e-14, fatal = true)
  n = A.n
  correct = true
  for i in 1:n
    if (isTerm[i])
      expected = initVal[i]
    else
      nbrs = A.rowval[A.colptr[i]:A.colptr[i+1]-1]
      maxval = maximum(lex[nbrs])
      minval = minimum(lex[nbrs])
      expected = minval + (maxval - minval) / 2.0
    end
    if (abs(lex[i]-expected) > eps)
      if (fatal)
        @printf("Expect value %f at vertex %d, but got %f\n", expected, i, lex[i])
      end
      return false
    end
  end
  return true
end

# StarSteepestPath: all terminal nodes in T connects to a single non-terminal node u.
# Assignments are v; distances from u to terminals are d.
#
# Naive version: enumerate pairs of vertices t1, t2 in T
# that maximizes (v(t1)-v(t2))/(d(t1)-d(t2)). O(|T|^2)
function NaiveStarSteepestPath(T, v, d)
  max = -Inf
  res1 = 0
  res2 = 0
  for i in 1:length(T)
    for j in (i+1):length(T)
      if (d[i] == Inf || d[j] == Inf)
        continue
      end
      grad = (v[i] - v[j]) / (d[i] + d[j])
      if (grad > max)
        max = grad
        res1 = i
        res2 = j
      end
      if (-grad > max)
        max = -grad
        res1 = j
        res2 = i
      end
    end
  end
  return T[res1], T[res2]
end

# grad error tolerance
gradTol = 1e-14

# Improved version. Expected O(|T|). (Alg 10 in Lex Paper)
function StarSteepestPath(T, v, d)
  if (LEX_DEBUG)
    println("[Star] ", T', " ", v', " ", d')
  end
  numT = length(T)
  if (numT != length(v) || numT != length(d))
    error("The dimensions of T, v, and d do not match!")
  end

  if (numT < 1)
    error("No terminal in the star graph, cannot process!")
  end
  if (numT == 1)
    return (T[1], T[1])
  end
  if (numT == 2)
    return (v[1] > v[2])? (T[1],T[2]) : (T[2], T[1])
  end

  i = rand(1:length(T))
  alpha = -Inf
  i2 = 0
  for j in 1:length(T)
    grad = abs(v[i] - v[j]) / (d[i] + d[j])
    if (grad > alpha)
      alpha = grad
      i2 = j
    end
  end

  vPlusArr = copy(v)
  vMinusArr = copy(v)
  for j in 1:length(T)
    vPlusArr[j] = v[j] + alpha * d[j]
    vMinusArr[j] = v[j] - alpha * d[j]
  end

  vLow = minimum(vPlusArr)
  vHigh = maximum(vMinusArr)

  T2 = falses(length(T))
  for j in 1:length(T)
    if (vMinusArr[j] > vLow * (1 + gradTol)
        || vPlusArr[j] < vHigh * (1 - gradTol))
      T2[j] = true
    end
  end
  if (findfirst(T2, true) == 0 || vLow >= vHigh)
    return (v[i] > v[i2])? (T[i],T[i2]) : (T[i2], T[i])
  else
    if (findfirst(T2, false) == 0)
      println(T, v, d)
      error("Infinite recursion!")
    end
    return StarSteepestPath(T[T2], v[T2], d[T2])
  end
end

# VertexSteepestPath. (Alg 9 in Lex Paper)
# does NOT work for directed graphs
function VertexSteepestPath(A, cpnts, isTerm, val, x)
  n = A.n
  d, parents = termFreeShortestPaths(A, x, isTerm)
  if (isTerm[x])
    maxGrad = -Inf
    y = x

    # at this point the graph may not be connected
    for i in 1:n
      if (isTerm[i] && cpnts[i] == cpnts[x])
        grad = abs(val[x]-val[i])/(d[i])
        if (grad > maxGrad)
          maxGrad = grad
          y = i
        end
      end
    end

    path = pathFromParents(parents, y)

    if (LEX_DEBUG)
      println(x, " ", y, " ", path)
    end

    if (val[x] > val[y])
      return reverse(path)
    else
      return path
    end
  else
    # filter out terminals in other components
    filter = copy(isTerm)
    for i in find(isTerm)
      if (d[i] == Inf)
        filter[i] = false
      end
    end
    if (LEX_DEBUG)
      println(x, "; filter = ", find(filter))
    end

    t1, t2 = StarSteepestPath(find(filter), val[filter], d[filter])
    p1 = pathFromParents(parents, t1)
    p2 = pathFromParents(parents, t2)
    if (LEX_DEBUG)
      println(t1, " ", t2, " ", p1', " ",  p2', " ", val[isTerm], " ", d[isTerm])
    end
    return vcat(p1[1:end-1], reverse(p2))
  end
end

# given a path from terminal to terminal, compute its gradient
function gradOfPath(A, isTerm, val, path)
  l = length(path)
  if (l < 2)
    return 0.0
  end
  @assert isTerm[path[1]] @sprintf("Node %d is not a terminal.", path[1])
  @assert isTerm[path[l]] @sprintf("Node %d is not a terminal.", path[l])

  dist = 0.0
  for i in 1:l-1
    # verify that it is actually a path
    @assert A[path[i], path[i + 1]] != 0.0
    @sprintf("Vertices %d and %d are not connected.", path[i], path[i + 1])

    dist += A[path[i], path[i + 1]]
  end

  dVolt = val[path[1]] - val[path[l]]
  return dVolt/dist
end # gradOfPath

function SteepestPath(A, cpnts, isTerm, val, vert=collect(1:n))

  n = A.n

  # do the following in each component, recurse at the same time

  # randomly pick an edge (x1, x2), and a vertex x3
  ind = rand(1:A.colptr[end]-1)
  x1 = findfirst(x -> x > ind, A.colptr) - 1
  x2 = A.rowval[ind]
  x3 = rand(vert)

  maxPath = VertexSteepestPath(A, cpnts, isTerm, val, x1)
  if (LEX_DEBUG)
    println(isTerm)
  end
  if (LEX_DEBUG)
    print("[SteepestPath] x1: ")
    print(x1)
    println(maxPath')
  end
  maxGrad = gradOfPath(A, isTerm, val, maxPath)

  p2 = VertexSteepestPath(A, cpnts, isTerm, val, x2)
  if (LEX_DEBUG)
    print("[SteepestPath] x2: ")
    print(x2)
    println(p2')

  end
  g2 = gradOfPath(A, isTerm, val, p2)
  if (g2 > maxGrad)
    maxPath = p2
    maxGrad = g2
  end

  p3 = VertexSteepestPath(A, cpnts, isTerm, val, x3)
  if (LEX_DEBUG)
    print("[SteepestPath] x3: ")
    print(x3)
    println(p3')
  end
  g3 = gradOfPath(A, isTerm, val, p3)

  if (g3 > maxGrad)
    maxPath = p3
    maxGrad = g3
  end

  hPresSet = CompHighPressGraph(A, isTerm, val, maxGrad)

  if (LEX_DEBUG)
    print("[SteepestPath] maxPath = ")
    println(maxPath')
    print("[SteepestPath] high pressure graph then becomes ")
    println(hPresSet')
  end

  if (length(hPresSet) == 0)
    return maxPath
  end

  # construct G' to recurse
  # TODO: make this cleaner and much faster, preferably without copying
  # by using a filter or something.
  t = copy(isTerm)
  t[hPresSet] = true
  newA = copy(A)
  for i = 1:n
    for j = 1:n
      if (!t[i] || !t[j])
        newA[i,j] = 0.0
      end
    end
  end

  newCpnts = components(newA)
  newVertices = []
  append!(newVertices, find(isTerm))
  append!(newVertices, hPresSet)
  return SteepestPath(newA, newCpnts, isTerm, val, newVertices)
end # SteepestPath

# remove (x,y) in E(G) such that both x and y are terminals;
# returns max_{x,y} grad[x,y]
function removeTerm2TermEdges!{Tv, Ti}(G::SparseMatrixCSC{Tv, Ti},
                                       isTerm::Array{Bool, 1},
                                       val::Array{Float64, 1}, )
  n = G.n
  terms = find(isTerm)
  alpha = -Inf

  for u in terms
    for v in terms
      if (G[u,v] > 0)
        grad = abs(val[u] - val[v]) / G[u, v]
        alpha = max(alpha, grad)
        # remove terminal to terminal edges
        G[u,v] = 0.0
        G[v,u] = 0.0
      end
    end
  end

  return alpha
end


# G may change (terminal to terminal edges will be removed)
# returns a new assignment
function CompInfMin{Tv, Ti}(G::SparseMatrixCSC{Tv, Ti},
                            isTerm::Array{Bool, 1},
                            val::Array{Float64, 1}, )
  alpha = removeTerm2TermEdges!(G, isTerm, val)

  cpnts = components(G)

  p = SteepestPath(G, cpnts, isTerm, val)
  grad = gradOfPath(G, isTerm, val, p)
  alpha = max(alpha, grad)

  vl, vlp = CompVLow(G, isTerm, val, alpha)
  vh, vhp = CompVHigh(G, isTerm, val, alpha)

  newVal = copy(val)
  for i in 1:G.n
    if (!isTerm[i])
      newVal[i] = vl[i] + (vh[i] - vl[i]) / 2.0
    end
  end

  return newVal
end # CompInfMin


function CompLexMin(A, isTerm, initVal)
  visited = copy(isTerm)

  B = sparse(copy(A))
  val = copy(initVal)

  # we keep going until every vertex has an assignment
  while (findfirst(visited, false) != 0)
    # remove terminal to terminal edges
    removeTerm2TermEdges!(B, visited, val)

    cpnts = components(B)
    p = SteepestPath(B, cpnts, visited, val)

    # visited will change here:
    fixPath!(B, visited, val, p)

    if (LEX_DEBUG)
      print("[CompLexMin] found steepest path: ")
      println(p')
      println(val')
    end
  end
  return val
end

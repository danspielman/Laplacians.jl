# lex:
#
# contains algorithms in the lex-minimizer paper
# includes iterative lex algorithm
#
# Started by xiao.shi@yale.edu on Sep 22, 2015
# Contributers:

include("graphGenerators.jl")
include("graphAlg.jl")

# takes an array of terminal values (null means non-terminal nodes)
# and returns an array of (node, value) pairs
function termValsToPairs(T::Array)
  n = size(T, 1)
  termPairs = Array[]
  for i = 1:n
    if (T[i] != null)
      push!(termPairs, [i, T[i]])
    end
  end
  return termPairs
end

#=
numIter: number of iterations

A: n by n adjacency matrix of the graph

initVal: 1 by n matrix of initial voltage assignments

isTerm: 1 by n boolean matrix of whether each vertex is a terminal
(value cannot change)

=#
function sim(numIter, A, initVal, isTerm)
  # A should be square
  (n, m) = size(A)
  val = copy(initVal)
  maxV = maximum(initVal)
  minV = minimum(initVal)

  for i = 1:numIter
    nextVal = copy(val)
    for u = 1:n
      if (!isTerm[u])
        maxNeighbor = minV
        minNeighbor = maxV
        for v = 1:n
          if (A[u,v] != 0)
            maxNeighbor = max(maxNeighbor, val[v])
            minNeighbor = min(minNeighbor, val[v])
          end
        end
        nextVal[u] = minNeighbor + (maxNeighbor - minNeighbor) / 2.0
      end
    end
    val = nextVal
  end
  return val
end

function testIterLexOnPath(n, numIter)
  IterLexPathCheckN(n)
j
  # path graph with n verticies
  A = grid2(1, n)
  # terminal assignments
  isTerm = falses(n)
  isTerm[1] = true
  isTerm[n] = true

  termVal = zeros(1,n)
  termVal[n] = 1

  # initial assignments
  initVal = copy(termVal)

  return sim(numIter, A, initVal, isTerm)
end

#=
  Experiments on path graphs

  For Pn, we consider the path with n vertices. Vertex 1 and n are
  terminals: T = {1,n}. The initial values for vertex n is 1, and 0
  otherwise.

  The final solution (lex-minimizer): voltage at vertex i is
  (i-1)/(n-1)

  We then consider the difference between current voltage assignment
  and the minimizer. Call it y_t at iteration t. We get the following
  relation:

  y_t = (1/2)A_{P_{n-2}} * y_{t-1}
=#

function IterLexPathCheckN(n::Int64)
  if (n < 3)
    error("n is too small, must be at least 3.")
  end
end

# returns the matrix (1/2)A_{P_{n-2}}
function IterLexPathOperator(n::Int64)
  IterLexPathCheckN(n)
  H = sparse(diagm(vec(ones(n-3,1)),1))
  return (H' + H)/2
end

# returns a n*1 vector v*(i) = (i-1)/(n-1)
function IterLexPathLexMinimizer(n::Int64)
  IterLexPathCheckN(n)
  return linspace(0, n-1, n) / (n-1)
end

# returns y0, a (n-2)*1 vector whose ith component is -i/(n-1)
function IterLexPathInitVal(n::Int64)
  IterLexPathCheckN(n)
  return linspace(-1, -n+2, n-2) / (n-1)
end

function IterLexPathSim(n::Int64, numIter::Int64)
  IterLexPathCheckN(n)
  assert(numIter >= 0)
  W = IterLexPathOperator(n)
  y = IterLexPathInitVal(n)

  if (numIter == 0)
    return y
  end
  return W^numIter * y
end

function IterLexPathVal(n::Int64, numIter::Int64)
  IterLexPathCheckN(n)
  assert(numIter >= 0)
  return [0; IterLexPathSim(n, numIter); 0] + IterLexPathLexMinimizer(n)
end


#= Lex Algorithms =#


function ModDijkstra(G, isTerm, initVal, alpha)
  n = G.n
  finished = zeros(Bool, n)

  nh = intHeap(n)
  key = nh.keys
  val = initVal

  parent = zeros(n)
  for i = 1:n
    if (isTerm[i])
      key[i] = initVal[i]*1.0
    end
    intHeapAdd!(nh, i, key[i]);
  end

  while (nh.nitems > 0)
    u = intHeapPop!(nh)
    val[u] = key[u]
    finished[u] = true
    for y = 1:n
      if (G[u,y]>0 && !finished[y] && !isTerm[y])
        newkey = val[u] + alpha*G[u,y]
        if (key[y]>newkey)
          intHeapAdd!(nh, y, newkey)
          parent[y] = u
        end
      end
    end
  end

  return val, parent
end

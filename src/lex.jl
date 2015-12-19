# lex:
#
# contains algorithms in the lex-minimizer paper
# includes iterative lex algorithm
#
# Started by xiao.shi@yale.edu on Sep 22, 2015
# Contributers:

LEX_DEBUG = true # debug flag for the lex module
LEX_EPS = 1e-12 # default error tolerance

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
      @printf("INFO: simIterLexUnwtd: terminating early after %d iterations, as numerical error prevents further progress.\n", t)
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
    if (abs(lex[i] - expected) > eps)
      if (fatal)
        error(@sprintf("The average of the max and min neighbors is %f at vertex %d, but got %f\n", expected, i, lex[i]))
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

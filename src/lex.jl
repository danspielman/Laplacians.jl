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

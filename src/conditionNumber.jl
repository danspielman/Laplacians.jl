#=

Code for checking how well one graph approximates another,
and for checking the quality of a preconditioner.

=#

"""
    sup12, sup21 = support(graph1, graph2; tol=1e-5)

Computes the support of graph1 wrt graph2, and the other way around.
It is randomized, so you might want to run it again if you don't trust the answers.
"""
function support(a1,a2; tol=1e-5)
    la1 = lap(a1)
    la2 = lap(a2)
    f1 = approxchol_lap(a1,tol=tol)
    f2 = approxchol_lap(a2,tol=tol)
    op12 = SqLinOp(false,1.0,size(a1,1),b->la2*f1(b))
    op21 = SqLinOp(false,1.0,size(a2,1),b->la1*f2(b))

    if isConnected(a1)
      sup12 = abs(eigs(op12;nev=1,which=:LM,tol=tol)[1][1])
    else
      sup12 = Inf
    end
    if isConnected(a2)
      sup21 = abs(eigs(op21;nev=1,which=:LM,tol=tol)[1][1])
    else
      sup21 = Inf
    end


    return sup12, sup21
end

"""
    eps = approxQual(graph1, graph2; tol=1e-5, verbose=false)

Computes the eps for which graph1 and graph2 are eps approximations of each other.
That is, L1 <= (1+eps) L2, and vice versa.

It is randomized, so you might want to run it again if you don't trust the answers.
"""
function approxQual(a1,a2; verbose=false, tol=1e-5)
    sup12, sup21 = support(a1, a2, tol=tol)

    if verbose
        println("support12: ", sup12, ", support21: ", sup21)
    end

    return max(sup12-1, sup21-1)
end

"""
    kapps = conditionNumber(graph1, graph2; tol=1e-5, verbose=false)

Computes the relative condition number of graph1 and graph2.

It is randomized, so you might want to run it again if you don't trust the answers.
"""
function conditionNumber(a1::SparseMatrixCSC, a2::SparseMatrixCSC; tol=1e-5, verbose=false)
    sup12, sup21 = support(a1, a2, tol=tol)

    if verbose
        println("support12: ", sup12, ", support21: ", sup21)
    end

    return sup12*sup21
end

"""
    kappa = conditionNumber(graph, precon; tol=1e-5, verbose=false)

Computes the relative condition number of graph and a preconditioning function.

It is randomized, so you might want to run it again if you don't trust the answers.
"""
function conditionNumber(a::SparseMatrixCSC, f::Function; tol=1e-5, verbose=false)

  la = lap(a)

  op = SqLinOp(false,1.0,size(a,1),b->la*f(b))

  upper = abs(eigs(op;nev=1,which=:LM,tol=tol)[1][1])

  op = SqLinOp(false,1.0,size(a,1),b->(upper*b - la*f(b)))

  lower = upper - abs(eigs(op;nev=2,which=:LM,tol=tol)[1][2])

  if verbose
      println("lower: ", lower, ", upper: ", upper)
  end
  return upper/lower
end

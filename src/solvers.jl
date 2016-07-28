#=

Code for solving Laplacian and Diagnally Dominant Systems
This is not all of it, just the short and sweet routines.

Started by Dan Spielman
Contributors: ???

  lapWrapSolver: takes a solver for DD systems, and uses it to solve a lap system in la
  lapWrapSolver(solver, la::AbstractArray)
  lapWrapSolver(solver)
  lapWrapSolver(solver, la::AbstractArray, b) = lapWrapSolver(solver,la)(b)    

  For example, to make a Cholesky-based solver for Laplacians, we created
  lapChol = lapWrapSolver(cholfact)


  augmentTree : takes a spanning tree, a graph, and adds in 2k edges of high stretch
    takes both as sparse matrices


  augTreeSolver : a solver that shouldn't suck

=#

include("sampler.jl")
include("randTrees.jl")



#==================================
  Wrapping solver for Laplacians
==================================#

#=
function lapWrapSolver(solver, la::AbstractMatrix)
    N = size(la)[1]
    lasub = la[1:(N-1),1:(N-1)]
    subSolver = solver(lasub);
    
    f = function(b)
        b = b - mean(b)
        bs = b[1:(N-1)]
        if isa(subSolver,Factorization)
            xs = subSolver \ bs
        else
            xs = subSolver(bs)
        end
        x = [xs;0]
        x = x - mean(x)
        return x
    end
    
    return f
end
=#


"""Apply the ith solver on the ith component"""
function blockSolver(comps, solvers)

    f = function(b)
        x = zeros(size(b))
        for i in 1:length(comps)
            ind = comps[i]
            bi = b[ind]
            if isa(solvers[i],Factorization)
                x[ind] = solvers[i] \ bi
            else
                x[ind] = (solvers[i])(bi)
            end
        end
        return x
    end
        
end



"""Takes a solver for solving nonsingular sdd systems,
and returns a solver for solving Laplacian systems.
The optional args tol and maxits are not necessarily taken by
all solvers.  But, if they are, one can pass them here"""
function lapWrapSolver(solver, la::AbstractArray; tol::Real=0.0, maxits::Integer=0)
    N = size(la)[1]

    lasub = la[1:(N-1),1:(N-1)]

    # need to be sure that removing does not disconnect.
    # or, if it does, solve on the components!

    co = components(lasub)
    if maximum(co) == 1

        if tol > 0
            if maxits > 0
                subSolver = solver(lasub, tol=tol, maxits=maxits);
            else
                subSolver = solver(lasub, tol=tol);
            end
        else
            if maxits > 0
                subSolver = solver(lasub, maxits=maxits);
            else
                subSolver = solver(lasub);
            end
        end

    else

        comps = vecToComps(co)

        solvers = []
        for i in 1:length(comps)
            ind = comps[i]
            
            lasubsub = lasub[ind,ind]

            if (length(ind) == 1)
                ssubSolver = x->(x/lasubsub[1,1])
            
            elseif (length(ind) < 50)
                ssubSolver = cholfact(lasubsub)

            else

                if tol > 0
                    if maxits > 0
                        ssubSolver = solver(lasubsub, tol=tol, maxits=maxits);
                    else
                        ssubSolver = solver(lasubsub, tol=tol);
                    end
                else
                    if maxits > 0
                        ssubSolver = solver(lasubsub, maxits=maxits);
                    else
                        ssubSolver = solver(lasubsub);
                    end
                end

            end
            push!(solvers, ssubSolver)
        end


        subSolver = blockSolver(comps,solvers)

    end
        

    f = function(b)
        b = b - mean(b)
        bs = b[1:(N-1)]
        if isa(subSolver,Factorization)
            xs = subSolver \ bs
        else
            xs = subSolver(bs)
        end
        x = [xs;0]
        x = x - mean(x)
        return x
    end
    
    return f
end


#=
function lapWrapSolver(solver)
    f = function(la::AbstractMatrix)
        return lapWrapSolver(solver, la)
    end
    return f
end
=#

function lapWrapSolver(solver; tol::Real=0.0, maxits::Integer=0)
    f = function(la::AbstractArray)
        return lapWrapSolver(solver, la, tol=tol, maxits=maxits)
    end
    return f
end


"""lapChol uses cholesky factorization to solver systems in Laplacians"""    
lapChol = lapWrapSolver(cholfact)

#lapWrapSolver(solver, la::AbstractMatrix, b) = lapWrapSolver(solver,la)(b)    

lapWrapSolver(solver, la::AbstractArray, b; tol::Real=0.0, maxits::Integer=0) = lapWrapSolver(solver,la, tol=tol, maxits=maxits)(b)    
    




#=========================================
  Augmented Spanning Tree Preconditioner
=========================================#



"""Takes as input a tree and an adjacency matrix of a graph.
It then computes the stretch of every edge of the graph wrt
the tree.  It then adds back the k edges of highest stretch,
and k edges sampled according to stretch"""
function augmentTree{Tv,Ti}(tree::SparseMatrixCSC{Tv,Ti}, mat::SparseMatrixCSC{Tv,Ti}, k::Ti)

    st = compStretches(tree, mat)

    # just to be safe, remove the tree from this
    #=
    (ti,tj) = findnz(tree)
    for i in 1:length(ti)   # most of the time
        st[ti[i],tj[i]] = 0
        st[tj[i],ti[i]] = 0
    end
    =#
    
    (ai,aj,av) = findnz(triu(st))

    ord = sortperm(av, rev=true)

    edgeinds = zeros(Bool,length(av))
    for i in 1:k
        edgeinds[ord[i]] = true
    end

    s = sum(av[ord[(k+1):end]])
    probs = av * k / s
    probs[ord[1:k]] = 0
    edgeinds[rand(length(av)) .< probs] = true

    augi = ai[edgeinds]
    augj = aj[edgeinds]
    augm = length(augi)
    augv = zeros(Tv, augm)
    for i in 1:augm,
        augv = mat[augi[i],augj[i]]
    end

    n = size(mat)[1]
    aug = sparse(augi, augj, augv, n, n)
    aug = aug + aug'

    return tree + aug
    
end


    
"""This is an augmented spanning tree preconditioner for diagonally dominant
linear systems.  It takes as optional input a tree growing algorithm.
It adds back 2sqrt(n) edges via augmentTree: the sqrt(n) of highest stretch
and another sqrt(n) sampled according to stretch.
For most purposes, one should directly call `augTreeSolver`."""
function augTreePrecon{Tv,Ti}(ddmat::SparseMatrixCSC{Tv,Ti}; treeAlg=akpw)

  adjmat = -triu(ddmat,1)
  adjmat = adjmat + adjmat'

  tree = treeAlg(adjmat)

  n = size(ddmat)[1]

  augtree = augmentTree(tree,adjmat,convert(Int,round(sqrt(n))))

  Dx = spdiagm(ddmat*ones(n))

  augDD = Dx + spdiagm(augtree*ones(n)) - augtree

  F = cholfact(augDD)

  return x -> (F\x)

end

"""
An "augmented spanning tree" solver for positive definite diagonally dominant matrices.
It works by adding edges to a low stretch spanning tree.  It calls `augTreePrecon` to form
the preconditioner.

~~~julia
 augTreeSolver{Tv,Ti}(ddmat::SparseMatrixCSC{Tv,Ti}; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, treeAlg=akpw)
~~~
"""
function augTreeSolver{Tv,Ti}(ddmat::SparseMatrixCSC{Tv,Ti}; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, treeAlg=akpw)

  F = augTreePrecon(ddmat, treeAlg=treeAlg)

  f(b;maxits=maxits, maxtime=maxtime, verbose=verbose) = pcg(ddmat, b, F, tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose)
  
    
  return f

end

"""This is an augmented spanning tree preconditioner for Laplacians.
It takes as optional input a tree growing algorithm.
It adds back 2sqrt(n) edges via `augmentTree`: the sqrt(n) of highest stretch
and another sqrt(n) sampled according to stretch.
For most purposes, one should directly call `augTreeLapSolver`."""
function augTreeLapPrecon{Tv,Ti}(ddmat::SparseMatrixCSC{Tv,Ti}; treeAlg=akpw)

  adjmat = -triu(ddmat,1)
  adjmat = adjmat + adjmat'

  tree = treeAlg(adjmat)

  n = size(ddmat)[1]

  augtree = augmentTree(tree,adjmat,convert(Int,round(sqrt(n))))

  Dx = spdiagm(ddmat*ones(n))

  augDD = Dx + spdiagm(augtree*ones(n)) - augtree

  F = lapWrapSolver(cholfact,augDD)

  return F

end

"""
An "augmented spanning tree" solver for Laplacian matrices.
It works by adding edges to a low stretch spanning tree.  It calls `augTreeLapPrecon` to form
the preconditioner.

~~~julia
 augTreeLapSolver{Tv,Ti}(ddmat::SparseMatrixCSC{Tv,Ti}; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, treeAlg=akpw)
~~~
"""
function augTreeLapSolver{Tv,Ti}(ddmat::SparseMatrixCSC{Tv,Ti}; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, treeAlg=akpw)

  F = augTreeLapPrecon(ddmat, treeAlg=treeAlg)

  f(b;maxits=maxits, maxtime=maxtime, verbose=verbose) = pcg(ddmat, b, F, tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose)
    
  return f

end
     
"""
A wrapper for the PyAMG solver.

~~~julia
 amgSolver{Tv,Ti}(ddmat::SparseMatrixCSC{Tv,Ti}; tol::Float64=1e-6, maxits=Inf, maxtime=Inf, verbose=false)
~~~
"""
function amgSolver{Tv,Ti}(ddmat::SparseMatrixCSC{Tv,Ti}; tol::Float64=1e-6, maxits=Inf, maxtime=Inf, verbose=false)

  amg = PyAMG.RugeStubenSolver(ddmat);
  M = PyAMG.aspreconditioner(amg);
  function F(b)
    return M \ b;
  end

  f(b;maxits=maxits, maxtime=maxtime, verbose=verbose) = pcg(ddmat, b, F, tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose)

  return f
  
end

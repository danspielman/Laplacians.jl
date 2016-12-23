"""
A wrapper for the PyAMG solver.

~~~julia
 amgSolver{Tv,Ti}(ddmat::SparseMatrixCSC{Tv,Ti}; tol::Float64=1e-6, maxits=Inf, maxtime=Inf, pcgIts=Int[], verbose=false)
~~~
"""
function AMGSolver{Tv,Ti}(ddmat::SparseMatrixCSC{Tv,Ti}; tol::Float64=1e-6, maxits=Inf, maxtime=Inf, pcgIts=Int[], verbose=false)

  amg = PyAMG.RugeStubenSolver(ddmat);
  M = PyAMG.aspreconditioner(amg);
  function F(b)
    return M \ b;
  end

    tol_=tol
    maxits_=maxits
    maxtime_=maxtime
    verbose_=verbose
    pcgIts_=pcgIts

    
  f(b;tol=tol_,maxits=maxits_, maxtime=maxtime_, pcgIts=pcgIts_, verbose=verbose_) = pcg(ddmat, b, F, tol=tol, maxits=maxits, maxtime=maxtime, pcgIts=pcgIts, verbose=verbose)

  return f
  
end


"""
A wrapper for the PyAMG solver. In line with our other solvers, takes in an adjacency matrix.

~~~julia
 amgSolver{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti}; tol::Float64=1e-6, maxits=Inf, maxtime=Inf, pcgIts=Int[], verbose=false)
~~~
"""
function AMGLapSolver{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti}; tol::Float64=1e-6, maxits=Inf, maxtime=Inf, pcgIts=Int[], verbose=false)

  la = lap(a)

  amg = PyAMG.RugeStubenSolver(la);
  M = PyAMG.aspreconditioner(amg);
  function F(b)
    return M \ b;
  end

    tol_=tol
    maxits_=maxits
    maxtime_=maxtime
    verbose_=verbose
    pcgIts_=pcgIts

    
  f(b;tol=tol_,maxits=maxits_, maxtime=maxtime_, pcgIts=pcgIts_, verbose=verbose_) = pcg(la, b, F, tol=tol, maxits=maxits, maxtime=maxtime, pcgIts=pcgIts, verbose=verbose)

  return f
  
end

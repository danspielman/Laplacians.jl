#=

This file contains code for calling Matlab solvers from inside Julia.
To the extent possible, it uses the same interfaces as the solvers in Laplacians.

They were mainly written to allow speed testing of linear solvers.

=#

using MATLAB
test_var_ = 1
@mput test_var_
mat"""
disp('Connected to Matlab!')
"""


"""
    f = matlab_ichol(sdd)

Computes a zero-fill incomplete Cholesky factorization of sdd, and returns a function that solves systems of linear equations in it.  Uses the symrcm ordering of the matrix.  Requires that MATLAB.jl be included.
"""
function matlab_ichol(sdd)

    if ~isdefined(:MATLAB)
        error("Type using MATLAB before calling this routine.")
    end

    n = sdd.n
    
    @mput sdd
    
    mat"""
    p = symrcm(sdd);
    laperm = sdd(p,p);
    L = ichol(laperm)
    """

    @mget p
    @mget L

    p = convert(Array{Int,1},vec(p))
    pi = zeros(Int,n)
    pi[p] = collect(1:n)
    
    Ltri = LowerTriangular(L)
    Utri = Ltri'
    
    # Create the solver function
    f = function(b::Array{Float64,1})

        y = Utri \ (Ltri \ b[p])
        return y[pi]

    end

    return f
end

"""
    f = matlab_ichol_sddm(sddm; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false)

Computes a zero-fill incomplete Cholesky factorization of sddm, and returns a function that solves systems of linear equations in it.  Uses the symrcm ordering of the matrix.  Requires that MATLAB.jl be included.
"""    
function matlab_ichol_sddm(sddm::SparseMatrixCSC; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[])

    F = matlab_ichol(sddm)
    tol_=tol
    maxits_=maxits
    maxtime_=maxtime
    verbose_=verbose
    pcgIts_=pcgIts
    
    f(b;tol=tol_,maxits=maxits_, maxtime=maxtime_, verbose=verbose_,pcgIts=pcgIts_) = pcg(sddm, b, F, tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose, pcgIts=pcgIts)
  
    return f

end


"""
    f = matlab_ichol_lap(A; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false)

Wrapper for matlab_ichol_sddm.
"""    
matlab_ichol_lap = lapWrapSDDM(matlab_ichol_sddm)
    

"""
    x = matlabCmgSolver(mat, b; tol::Real=1e-6, maxits=10000)

This runs Koutis's CMG solver.  You must have installed the solver, and it must be on Matlab's default path.  This routine does not implement all of our preferred interface.  Use the same solver for sddm and Laplacian matrices.
"""
function matlabCmgSolver(mat::SparseMatrixCSC, b; tol::Real=1e-6, maxits=10000)

    maxits = min(10000,maxits)
    
    @mput mat
    @mput b
    @mput tol
    @mput maxits

    mat"""
    pfun = cmg_sdd(mat);

    [x,flag] = pcg(mat, b, tol, maxits, pfun);
    """

    @mget x
    
    return x
end

"""
    solver = matlabCmgSolver(mat; tol::Real=1e-6, maxits=10000)

This runs Koutis's CMG solver.  You must have installed the solver, and it must be on Matlab's default path.  This routine does not implement all of our preferred interface.  Use the same solver for sddm and Laplacian matrices.

Note that this does not build the solver.  Rather, it just makes sure that the call to `solver` will.
"""
function matlabCmgSolver(mat::SparseMatrixCSC; tol::Real=1e-6, maxits=1000, params...)

    tol_=tol
    maxits_=maxits
    
    f(b; tol=tol_, maxits=maxits_, params...) = matlabCmgSolver(mat, b; tol=tol, maxits=maxits)

    return f
end


"""
    solver = matlabCmgLap(a; tol::Real=1e-6, maxits=10000)

This runs Koutis's CMG solver.  You must have installed the solver, and it must be on Matlab's default path.  This routine does not implement all of our preferred interface.  

The input `a` should be the adjacency matrix.  This will find connected components, and then call the sddm version of the solver.

Note that this does not build the solver.  Rather, it just makes sure that the call to `solver` will.
"""
matlabCmgLap = lapWrapSDDM(matlabCmgSolver)

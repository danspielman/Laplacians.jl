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
    f = matlab_CMG(sdd)

This uses Matlab to call Koutis's CMG code to construct a preconditioner.
The function f that is returned applied the preconditioner.
"""
function matlab_CMG(sdd)

    if ~isdefined(:MATLAB)
        error("Type using MATLAB before calling this routine.")
    end

    n = sdd.n
    
    @mput sdd
    
    mat"""
    cmgOpts.display = 0;
    pfun = cmg_sdd(sdd,cmgOpts);
    """

    @mget pfun
    
    # Create the solver function
    f = function(b::Array{Float64,1})

        y = Utri \ (Ltri \ b[p])
        return y[pi]

    end

    return f
end


#=
"""
    x, err, cmg_build, cmg_solve = matlabCMGLapSolver(a, b; tol, maxits)

Calls an incomplete Cholesky factorization solver from Matlab.
Before using this, you need to include setUpMatlab.jl.
This requires Yiannis Koutis's CMG solver.
"""
function matlabCMGLapSolver(a, b; tol::Real=1e-6, maxits=10000)
    
    la = lap(a)
    opts = Opts(tol,maxits)

    @mput opts
    @mput la
    @mput b

    mat"""
    tic; f = lapCmgSolver(la,[],opts); cmg_build_time = toc

    tic; x = f(b); cmg_solve_time = toc;
    err = norm(la*x-b)/norm(b);
    """

    @mget x
    @mget err
    @mget cmg_build_time
    @mget cmg_solve_time
    
    return x, err, cmg_build_time, cmg_solve_time
end
=#

#=
"""
    x, err, cmg_build, cmg_solve = matlabCMGDDSolver(a, b; tol, maxits)

Calls an incomplete Cholesky factorization solver from Matlab.
Before using this, you need to include setUpMatlab.jl.
This requires Yiannis Koutis's CMG solver.
"""
function matlabCMGSDDSolver(sdd, b; tol::Real=1e-6, maxits=10000)
    
    opts = Opts(tol,maxits)

    @mput opts
    @mput sdd
    @mput b

    mat"""
    tic; f = cmgSolver(sdd,[],opts); cmg_build_time = toc

    tic; x = f(b); cmg_solve_time = toc;
    err = norm(sdd*x-b)/norm(b);
    """

    @mget x
    @mget err
    @mget cmg_build_time
    @mget cmg_solve_time
    
    return x, err, cmg_build_time, cmg_solve_time
end

=#

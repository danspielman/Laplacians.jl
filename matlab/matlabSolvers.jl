#=

This file contains code for calling Matlab solvers from inside Julia.
To the extent possible, it uses the same interfaces as the solvers in Laplacians.

They were mainly written to allow speed testing of linear solvers.

=#

immutable Opts
    tol::Float64
    maxit::Int
end
opts = Opts(1e-3,10000)


"""
    x, err, icc_build, icc_solve = matlabIccLapSolver(a, b; tol, maxits)

Calls an incomplete Cholesky factorization solver from Matlab.
Before using this, you need to include setUpMatlab.jl
"""
function matlabIccLapSolver(a, b; tol::Real=1e-6, maxits=10000)
    
    la = lap(a)
    opts = Opts(tol,maxits)

    @mput opts
    @mput la
    @mput b

    mat"""
    tic; f = iccSolver(la,[],opts); icc_build_time = toc

    tic; x = f(b); icc_solve_time = toc;
    err = norm(la*x-b)/norm(b);
    """

    @mget x
    @mget err
    @mget icc_build_time
    @mget icc_solve_time
    
    return x, err, icc_build_time, icc_solve_time
end

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
    tic; f = cmgSolver(la,[],opts); cmg_build_time = toc

    tic; x = f(b); cmg_solve_time = toc;
    err = norm(la*x-b)/norm(b);
    """

    @mget x
    @mget err
    @mget cmg_build_time
    @mget cmg_solve_time
    
    return x, err, cmg_build_time, cmg_solve_time
end

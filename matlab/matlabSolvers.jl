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
    tic; f = lapIccSolver(la,[],opts); icc_build_time = toc

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

"""
    x, err, icc_build, icc_solve = matlabIccSDDSolver(a, b; tol, maxits)

Calls an incomplete Cholesky factorization solver from Matlab.
Before using this, you need to include setUpMatlab.jl
"""
function matlabIccSDDSolver(sdd, b; tol::Real=1e-6, maxits=10000)
    
    opts = Opts(tol,maxits)

    @mput opts
    @mput sdd
    @mput b

    mat"""
    tic; f = iccSolver(sdd,[],opts); icc_build_time = toc

    tic; x = f(b); icc_solve_time = toc;
    err = norm(sdd*x-b)/norm(b);
    """

    @mget x
    @mget err
    @mget icc_build_time
    @mget icc_solve_time
    
    return x, err, icc_build_time, icc_solve_time
end

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

function matlab_ichol(sdd)

    n = sdd.n
    
    @mput sdd
    
    mat"""
    p = symrcm(sdd);
    laperm = sdd(p,p);
    L = ichol(laperm)
    """

    @mget p
    @mget L

    @show n
    @show length(p)
    
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

function matlab_ichol_sdd(sdd::SparseMatrixCSC; tol::Real=1e-6, maxits=Inf, maxtime=Inf)

    F = matlab_ichol(sdd)
    tol_=tol
    maxits_=maxits
    maxtime_=maxtime
    
    f(b;tol=tol_,maxits=maxits_, maxtime=maxtime_) = pcg(sdd, b, F, tol=tol, maxits=maxits, maxtime=maxtime, verbose=true)
  
    return f

end


function matlab_ichol_lap(la::SparseMatrixCSC; tol::Real=1e-6, maxits=Inf)

    tol_=tol
    maxits_=maxits
    
    return lapWrapSolver(matlab_ichol_sdd, la, tol=tol_, maxits=maxits_)
  
end
    

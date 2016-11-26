
using DataFrames

include("$(Pkg.dir("Laplacians"))/matlab/setUpMatlab.jl")
include("$(Pkg.dir("Laplacians"))/matlab/matlabSolvers.jl")

global speedDF = DataFrame()


"""

    function speedTestLapSolver{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti}, b::Array{Tv,1}; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false)

Runs many Laplacians solvers.  Puts the build and solve time results into solverDF.  Returns the answer given by KMP.

"""
function speedTestLapSolver{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti}, b::Array{Tv,1}; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false)


    la = lap(a)
    
    # augTree
    
    tic()
    f = augTreeLapSolver(a, tol=tol, maxtime=maxtime)
    aug_build = toq()

    tic()
    x = f(b)
    aug_solve = toq()

    aug_err = norm(la * x - b) / norm(b)

    # KMP
    
    tic()
    f = KMPLapSolver(a, tol=tol, maxtime=maxtime)
    kmp_build = toq()

    tic()
    xhat = f(b)
    kmp_solve = toq()

    kmp_err = norm(la * xhat - b) / norm(b)

    # AMG
    
    tic()
    f = AMGLapSolver(a, tol=tol, maxtime=maxtime)
    amg_build = toq()

    tic()
    xhat = f(b)
    amg_solve = toq()

    amg_err = norm(la * xhat - b) / norm(b)

    # Matlab

    _, cmg_err, cmg_build, cmg_solve = matlabCMGLapSolver(a, b; tol=tol)
    _, icc_err, icc_build, icc_solve = matlabIccLapSolver(a, b; tol=tol)

    nv = size(a,1)
    ne = nnz(a)
    hash_a = hash(a)
    hash_b = hash(b)
    
    row = DataFrame(nv = nv, ne = ne, hash_a = hash_a, hash_b = hash_b,
            kmp_build = kmp_build, kmp_solve = kmp_solve, kmp_err = kmp_err,
            aug_build = aug_build, aug_solve = aug_solve, aug_err = aug_err,
            amg_build = amg_build, amg_solve = amg_solve, amg_err = amg_err,
            icc_build = icc_build, icc_solve = icc_solve, icc_err = icc_err,
            cmg_build = cmg_build, cmg_solve = cmg_solve, cmg_err = cmg_err)


    global speedDF
    
    if size(speedDF,1) == 0
        speedDF = row
    else
        append!(speedDF,row)
    end
    
    return x

end



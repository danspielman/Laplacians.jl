
using DataFrames

include("$(Pkg.dir("Laplacians"))/matlab/setUpMatlab.jl")
include("$(Pkg.dir("Laplacians"))/matlab/matlabSolvers.jl")
include("$(Pkg.dir("Laplacians"))/src/isotonicIPM.jl")

global speedDF = DataFrame()


"""

    function speedTestLapSolver{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti}, b::Array{Tv,1}; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false)

Runs many Laplacians solvers.  Puts the build and solve time results into solverDF.  Returns the answer given by KMP.

"""
function speedTestLapSolver{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti}, b::Array{Tv,1}; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false)


    b = b - mean(b)
    
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


"""

    function speedTestSDDSolver{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti}, b::Array{Tv,1}; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false)

Runs many SDD solvers.  Puts the build and solve time results into solverDF.  Returns the answer given by KMP.

"""
function speedTestSDDSolver{Tv,Ti}(sdd::SparseMatrixCSC{Tv,Ti}, b::Array{Tv,1}; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false)



    
    # augTree
    
    tic()
    f = augTreeSolver(sdd, tol=tol, maxtime=maxtime)
    aug_build = toq()

    tic()
    x = f(b)
    aug_solve = toq()

    aug_err = norm(sdd * x - b) / norm(b)

    # KMP
    
    tic()
    f = KMPSDDSolver(sdd, tol=tol, maxtime=maxtime)
    kmp_build = toq()

    tic()
    xhat = f(b)
    kmp_solve = toq()

    kmp_err = norm(sdd * xhat - b) / norm(b)

    # AMG
    
    tic()
    f = AMGSolver(sdd, tol=tol, maxtime=maxtime)
    amg_build = toq()

    tic()
    xhat = f(b)
    amg_solve = toq()

    amg_err = norm(sdd * xhat - b) / norm(b)

    # Matlab

    _, cmg_err, cmg_build, cmg_solve = matlabCMGSDDSolver(sdd, b; tol=tol)
    _, icc_err, icc_build, icc_solve = matlabIccSDDSolver(sdd, b; tol=tol)

    nv = size(sdd,1)
    ne = nnz(sdd)-nv
    hash_a = hash(sdd)
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


function speedTestSDDSolver{Tv,Ti}(sdd::SparseMatrixCSC{Tv,Ti}; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false)
    f(b) = speedTestSDDSolver(sdd, b, tol=tol, maxits=maxits, maxtime=maxtime, verbose=verbose)
    return f
end

"""
    plotTests(DF::DataFrame)

Positive means first is better
"""
function plotTests(DF::DataFrame)
    
    solver_names = ["kmp","aug","amg", "icc", "cmg"]
    solver_syms = [:kmp_solve, :aug_solve, :amg_solve, :icc_solve, :cmg_solve]

    k = 0
    ns = length(solver_names)
    z = zeros(length(DF[solver_syms[1]]))
    
    for i in 1:ns
        for j in 1:ns
            k = k + 1
            if i != j
                subplot(ns,ns,k)
                
                plot(sort(log(DF[solver_syms[j]] ./ DF[solver_syms[i]])))
                plot(z, color=:red)
                ax = gca()
                ax[:xaxis][:set_visible](false)
                title("$(solver_names[i]) vs $(solver_names[j])")
                
            end
        end
    end
    
end


function testInApps(a)

    n = size(a,1)
    v = randn(n)
    v = v / norm(v)
    v = v - mean(v)
    prev = zeros(n)
    it = 0
    while (norm(v-prev) > 0.05)
        prev = copy(v)
        v = speedTestLapSolver(a,v,tol=1e-6)
        v = v - mean(v)
        v = v / norm(v)
        it += 1
        println(it, " : ", norm(v - prev))
    end

    
    p = sortperm(v);
    ap = triu(a[p,p])
    vp = v[p];
    iso = isotonicIPM(ap,vp,eps=.01,solver=(H -> speedTestSDDSolver(H,tol=1e-1,maxits=100000)))

end

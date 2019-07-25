#==========================================================
Code for comparing solvers, with time limits

This will run tests of solvers inside a thread so that it can impose a
time limit.

For this reason, julia should either be started with workers "julia -p 2",
or workers should be added before including Laplacians, like
addprocs()
@everywhere using Laplacians
@everywere include(this_file)

Routines for plotting results and producing reports are in
compare_solvers_reports.jl

The timeout features won't work with the matlab solvers.


===========================================================#

using Statistics

"""
    function speedTestLapTL{Tv,Ti}(solvers, dic, a::SparseMatrixCSC{Tv,Ti}, b::Array{Tv,1}; tol::Real=1e-2, maxits=Inf, maxtime=Inf, verbose=false)

Runs many Laplacians solvers.  Puts the build and solve time results into a dictionary dic.  It would be easiest to look at it via DataFrame(dic).  Returns the answer from the last solver.  `solvers` should be an array of `SolverTest`.

Uses threads to time limit all solvers after the first to at most 10x the time.
"""
function speedTestLapTL(solvers, dic, a::SparseMatrixCSC{Tv,Ti}, b::Array{Tv,1}; tol::Real = 1e-2, maxits = 10000, maxtime = 1000, verbose = false, testName = "") where {Tv,Ti}

    b = b - mean(b) * ones(size(b))

    la = lap(a)

    it = Int[1]

    # init the dict, if it is the first time
    initDictCol!(dic, "nv", Int)
    initDictCol!(dic, "ne", Int)
    initDictCol!(dic, "hash_a", UInt64)
    initDictCol!(dic, "testName", String)

    solvecol(name) = "$(name)_solve"
    buildcol(name) = "$(name)_build"
    totcol(name) = "$(name)_tot"
    itscol(name) = "$(name)_its"
    errcol(name) = "$(name)_err"

    dic["names"] = [t.name for t in solvers]

    for solver in solvers
        name = solver.name
        initDictCol!(dic, solvecol(name), Float64)
        initDictCol!(dic, buildcol(name), Float64)
        initDictCol!(dic, totcol(name), Float64)
        initDictCol!(dic, itscol(name), Float64)
        initDictCol!(dic, errcol(name), Float64)
    end

    nv = size(a, 1)
    ne = nnz(a)
    hash_a = hash(a)

    push!(dic["nv"], nv)
    push!(dic["ne"], ne)
    push!(dic["hash_a"], hash_a)
    push!(dic["testName"], testName)

    x = []

    i = 1

    solverTest = solvers[i]

    if verbose
        println()
        println(solverTest.name)
    end


    t0 = time()
    ret = testSolver(solverTest.solver, a, b, tol, maxits, verbose)
    tl = 5 + 10 * (time() - t0)
    if verbose
        println("time limit: ", tl)
    end

    pushSpeedResult!(dic, solverTest.name, ret)

    for i in 2:length(solvers)
        solverTest = solvers[i]

        if verbose
            println()
            println(solverTest.name)
        end

        #time_limit = 10.0

        ret = testSolverTimed(tl, solverTest.solver, a, b, tol, maxits, verbose)

        pushSpeedResult!(dic, solverTest.name, ret)

    end

    return x

end


function testSolverTimed(timelimit, solver, a, b, tol, maxits, verbose)

    if length(workers()) < 1
        error("Run addprocs() before setting up this code.")
    end
    proc = workers()[1]
    remotecall_fetch(gc, proc)


    c = Channel(1)
    @async put!(c, remotecall_fetch(()->(testSolver(solver, a, b, tol, maxits, verbose)), proc))
    t0 = time()
    while !isready(c) && (time() - t0) < timelimit
        sleep(0.1)
    end

    ret = (Inf, Inf, Inf, Inf)

    if isready(c)
        ret = fetch(c)

    else
        println("Interrupt process at time", time() - t0)
        interrupt(proc)
        sleep(1)
        ret = (Inf, Inf, Inf, Inf)

    end

    return ret

end

import Laplacians.initDictCol!
import Laplacians.testSolver
import Laplacians.testSolverSddm
import Laplacians.pushSpeedResult!

"""
Runs many Laplacians solvers.  Puts the build and solve time results into a dictionary dic.  It would be easiest to look at it via DataFrame(dic).  Returns the answer from the last solver.  `solvers` should be an array of `SolverTest`.

Also compares them against the solvers we have in matlab, with a time limit of 10x the first solver here.
"""
function testVMatlabLap(solvers, dic::Dict, a::SparseMatrixCSC{Tv,Ti}, b::Array{Tv,1};
  tol::Real = 1e-8, maxits = 1000, maxtime = 1000, verbose = false, testName = "",
                        test_icc = false, test_cmg = false, test_lamg = false, test_muelubelos = false, test_hypre = false, tl_fac = 10, tl = 0) where {Tv,Ti}

    b = b .- mean(b)

    #la = lap(a) #TODO RAT we do this later, so don't need here

    it = Int[1]

    # init the dict, if it is the first time
    initDictCol!(dic, "nv", Int)
    initDictCol!(dic, "ne", Int)
    initDictCol!(dic, "hash_a", UInt64)
    initDictCol!(dic, "testName", String)

    solvecol(name) = "$(name)_solve"
    buildcol(name) = "$(name)_build"
    totcol(name) = "$(name)_tot"
    itscol(name) = "$(name)_its"
    errcol(name) = "$(name)_err"

    dic["names"] = String[]
    for t in solvers
        push!(dic["names"], t.name)
    end

    if test_hypre
        push!(dic["names"], "hypre")
    end
    if test_cmg
        push!(dic["names"], "cmg")
    end
    if test_icc
        push!(dic["names"], "icc")
    end
    if test_lamg
        push!(dic["names"], "lamg")
    end
    if test_muelubelos
        push!(dic["names"], "muelubelos")
    end

    for name in dic["names"]
        initDictCol!(dic, solvecol(name), Float64)
        initDictCol!(dic, buildcol(name), Float64)
        initDictCol!(dic, totcol(name), Float64)
        initDictCol!(dic, itscol(name), Float64)
        initDictCol!(dic, errcol(name), Float64)
    end

    nv = size(a, 1)
    ne = nnz(a)
    hash_a = hash(a)

    push!(dic["nv"], nv)
    push!(dic["ne"], ne)
    push!(dic["hash_a"], hash_a)
    push!(dic["testName"], testName)

    x = []

    #tl = 0

    for i in 1:length(solvers)
        solverTest = solvers[i]

        if verbose
            println()
            println(solverTest.name)
        end

        ret = testSolver(solverTest.solver, a, b, tol, maxits, verbose)

        if i == 1 && tl == 0
            x = ret[5]
            tl = round(Int, 30 + tl_fac * (ret[1] + ret[2]))
        end


        pushSpeedResult!(dic, solverTest.name, ret)

    end

    println()

    if tl == 0
        error("tl is zero")
    end

    if test_hypre || test_cmg || test_icc || test_lamg || test_muelubelos
        la = lap(a)


        if test_hypre
            if verbose
                println("--------------")
                println("hypre")
            end
        
            ret = timeLimitHypre(tl, la, b; verbose = true, num_procs = 2)
            pushSpeedResult!(dic, "hypre", ret)
        end
 
        if test_cmg
            if verbose
                println("--------------")
                println("cmg")
            end

            ret = timeLimitCmg(tl, la, b, verbose = true);
            pushSpeedResult!(dic, "cmg", ret)
        end

        if test_icc
            if verbose
                println("--------------")
                println("icc")
            end
            ret = timeLimitIcc(tl, la, b, verbose = true);
            pushSpeedResult!(dic, "icc", ret)
        end

        if test_lamg
            if verbose
                println("--------------")
                println("lamg")
            end
            ret = timeLimitLamg(tl, la, b, verbose = true);
            pushSpeedResult!(dic, "lamg", ret)
        end

        if test_muelubelos
            if verbose
                println("--------------")
                println("muelubelos")
            end
            ret = timeLimitMueluBelos(tl, la, b, verbose = true);
            pushSpeedResult!(dic, "muelubelos", ret)
        end
    end
    

    return x

end



function testVMatlabLap(solvers::Array, dic::Dict, maker::Function; testName = "")
    println(testName)
    a = maker()
    n = size(a, 1)
    b = randn(n);
    b = b - mean(b) * ones(size(b));
    testVMatlabLap(solvers, dic, a, b, testName = testName, verbose = true)
end


function testVMatlabSddm(solvers, dic::Dict, sddmat::SparseMatrixCSC{Tv,Ti}, b::Array{Tv,1};
  tol::Real = 1e-8, maxits = 1000, maxtime = 1000, verbose = false, testName = "",
  test_icc = true, test_cmg = true, test_lamg = true, test_hypre = false, tl_fac = 10) where {Tv,Ti}


    it = Int[1]

    # init the dict, if it is the first time
    initDictCol!(dic, "nv", Int)
    initDictCol!(dic, "ne", Int)
    initDictCol!(dic, "hash_a", UInt64)
    initDictCol!(dic, "testName", String)

    solvecol(name) = "$(name)_solve"
    buildcol(name) = "$(name)_build"
    totcol(name) = "$(name)_tot"
    itscol(name) = "$(name)_its"
    errcol(name) = "$(name)_err"

    dic["names"] = Array{String}(0)
    for t in solvers
        push!(dic["names"], t.name)
    end

    if test_hypre
        push!(dic["names"], "hypre")
    end
    if test_cmg
        push!(dic["names"], "cmg")
    end
    if test_icc
        push!(dic["names"], "icc")
    end
    if test_lamg
        push!(dic["names"], "lamg")
    end

    for name in dic["names"]
        initDictCol!(dic, solvecol(name), Float64)
        initDictCol!(dic, buildcol(name), Float64)
        initDictCol!(dic, totcol(name), Float64)
        initDictCol!(dic, itscol(name), Float64)
        initDictCol!(dic, errcol(name), Float64)
    end

    nv = size(sddmat, 1)
    ne = nnz(sddmat)
    hash_a = hash(sddmat)

    push!(dic["nv"], nv)
    push!(dic["ne"], ne)
    push!(dic["hash_a"], hash_a)
    push!(dic["testName"], testName)

    x = []

    tl = 0

    for i in 1:length(solvers)
        solverTest = solvers[i]

        if verbose
            println()
            println(solverTest.name)
        end

        ret = testSolver(solverTest.solver, sddmat, b, tol, maxits, verbose)

        if i == 1
            x = ret[5]
            tl = round(Int, 30 + tl_fac * (ret[1] + ret[2]))
        end


        pushSpeedResult!(dic, solverTest.name, ret)

    end

    if tl == 0
        error("tl is zero")
    end

    if test_hypre
        if verbose
            println("--------------")
            println("hypre")
        end
    
        ret = timeLimitHypre(tl, sddmat, b; verbose = true, num_procs = 2)
        pushSpeedResult!(dic, "hypre", ret)
    end

    if test_cmg
        if verbose
            println("cmg")
        end

        ret = timeLimitCmg(tl, sddmat, b, verbose = true);
        pushSpeedResult!(dic, "cmg", ret)
    end

    if test_icc
        if verbose
            println("icc")
        end
        ret = timeLimitIcc(tl, sddmat, b, verbose = true);
        pushSpeedResult!(dic, "icc", ret)
    end

    if test_lamg
        if verbose
            println("lamg")
        end
        ret = timeLimitLamgSddm(tl, sddmat, b, verbose = true);
        pushSpeedResult!(dic, "lamg", ret)
    end

    return x

end

function testSddm(solvers, dic::Dict, sddmmat::SparseMatrixCSC{Tv,Ti}, b::Array{Tv,1};
  tol::Real = 1e-8, maxits = 1000, maxtime = 1000, verbose = true, testName = "",
   test_hypre = true, test_icc = false, test_cmg = false, test_lamg = false, tl_fac = 10) where {Tv,Ti}


    it = Int[1]

    # init the dict, if it is the first time
    initDictCol!(dic, "nv", Int)
    initDictCol!(dic, "ne", Int)
    initDictCol!(dic, "hash_a", UInt64)
    initDictCol!(dic, "testName", String)

    solvecol(name) = "$(name)_solve"
    buildcol(name) = "$(name)_build"
    totcol(name) = "$(name)_tot"
    itscol(name) = "$(name)_its"
    errcol(name) = "$(name)_err"

    dic["names"] = String[]
    for t in solvers
        push!(dic["names"], t.name)
    end


    if test_hypre
        push!(dic["names"], "hypre")
    end
    if test_cmg
        push!(dic["names"], "cmg")
    end
    if test_icc
        push!(dic["names"], "icc")
    end
    if test_lamg
        push!(dic["names"], "lamg")
    end

    for name in dic["names"]
        initDictCol!(dic, solvecol(name), Float64)
        initDictCol!(dic, buildcol(name), Float64)
        initDictCol!(dic, totcol(name), Float64)
        initDictCol!(dic, itscol(name), Float64)
        initDictCol!(dic, errcol(name), Float64)
    end

    nv = size(sddmmat, 1)
    ne = nnz(sddmmat)
    hash_a = hash(sddmmat)

    push!(dic["nv"], nv)
    push!(dic["ne"], ne)
    push!(dic["hash_a"], hash_a)
    push!(dic["testName"], testName)

    x = []

    tl = 0

    for i in 1:length(solvers)
        solverTest = solvers[i]

        if verbose
            println("--------------")
            println(solverTest.name)
        end

        ret = testSolverSddm(solverTest.solver, sddmmat, b, tol, maxits, verbose)

        if i == 1
            x = ret[5]
            tl = round(Int, 30 + tl_fac * (ret[1] + ret[2]))
        end


        pushSpeedResult!(dic, solverTest.name, ret)

    end

    if tl == 0
        error("tl is zero")
    end

    if test_hypre
        if verbose
            println("--------------")
            println("hypre")
        end
      
        ret = timeLimitHypre(tl, sddmmat, b; verbose = false, num_procs = 2)
        pushSpeedResult!(dic, "hypre", ret)
    end

    if test_cmg
        if verbose
            println("--------------")
            println("cmg")
        end

        ret = timeLimitCmg(tl, sddmmat, b, verbose = true);
        pushSpeedResult!(dic, "cmg", ret)
    end

    if test_icc
        if verbose
            println("--------------")
            println("icc")
        end
        ret = timeLimitIcc(tl, sddmmat, b, verbose = true);
        pushSpeedResult!(dic, "icc", ret)
    end

    if test_lamg
        if verbose
            println("--------------")
            println("lamg")
        end
        ret = timeLimitLamgSddm(tl, sddmmat, b, verbose = true);
        pushSpeedResult!(dic, "lamg", ret)
    end

    return x

end

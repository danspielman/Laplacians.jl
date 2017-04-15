#==========================================================
Code for comparing solvers

This will run tests of solvers inside a thread so that it can impose a
time limit.

For this reason, julia should either be started with workers "julia -p 2",
or workers should be added before including Laplacians, like
addprocs()
@everywhere using Laplacians
@everywere inlclude(this_file)

Routines for plotting results and producing reports are in
compare_solvers_reports.jl

===========================================================#


"""
    SolverTest(solver, name)

Encloses a solver with its name, so that we can compare it in tests
"""
type SolverTest
    solver::Function
    name::String
end


"""
    initDictCol!(dic, name, typ)

For a dictionary in which each key indexes an array.
If dic does not contain an entry of `name`, create with set to `Array(typ,0)`.
"""
function initDictCol!(dic, name, typ)
    if ~haskey(dic,name)
        dic[name] = Array(typ,0)
    end
end



"""
    function speedTestLapSolvers{Tv,Ti}(solvers, dic, a::SparseMatrixCSC{Tv,Ti}, b::Array{Tv,1}; tol::Real=1e-2, maxits=Inf, maxtime=Inf, verbose=false)

Runs many Laplacians solvers.  Puts the build and solve time results into a dictionary dic.  It would be easiest to look at it via DataFrame(dic).  Returns the answer from the last solver.  `solvers` should be an array of `SolverTest`.
"""
function speedTestLapSolvers{Tv,Ti}(solvers, dic, a::SparseMatrixCSC{Tv,Ti}, b::Array{Tv,1}; tol::Real=1e-2, maxits=10000, maxtime=1000, verbose=false)

    b = b - mean(b)

    la = lap(a)

    it = Int[1]

    # init the dict, if it is the first time
    initDictCol!(dic, "nv", Int)
    initDictCol!(dic, "ne", Int)
    initDictCol!(dic, "hash_a", UInt64)

    solvecol(name) = "$(name)_solve"
    buildcol(name) = "$(name)_build"
    totcol(name) = "$(name)_tot"
    itscol(name) = "$(name)_its"
    errcol(name) = "$(name)_err"

    dic["names"] = [t.name for t in tests]

    for solver in solvers
        name = solver.name
        initDictCol!(dic, solvecol(name), Float64)
        initDictCol!(dic, buildcol(name), Float64)
        initDictCol!(dic, totcol(name), Float64)
        initDictCol!(dic, itscol(name), Float64)
        initDictCol!(dic, errcol(name), Float64)
    end

    nv = size(a,1)
    ne = nnz(a)
    hash_a = hash(a)

    push!(dic["nv"],nv)
    push!(dic["ne"],ne)
    push!(dic["hash_a"],hash_a)

    x = []

    for i in 1:length(solvers)
        solverTest = solvers[i]

        if verbose
            println()
            println(solverTest.name)
        end

        time_limit = 10.0
        if i == 1
          t0 = time()
          ret = testSolver(solverTest.solver, a, b, tol, maxits, verbose)
          time_limit = 10*(time()-t0)
        else
          ret = testSolverTimed(time_limit, solverTest.solver, a, b, tol, maxits, verbose)
        end

        name = solverTest.name
        push!(dic[solvecol(name)],ret[1])
        push!(dic[buildcol(name)],ret[2])
        push!(dic[totcol(name)],ret[1]+ret[2])
        push!(dic[itscol(name)],ret[3])
        push!(dic[errcol(name)],ret[4])

    end

    return x

end


function testSolver(solver, a, b, tol, maxits, verbose)

  try

    gc()
    tic()
    f = solver(a, tol=tol, maxits=maxits, verbose=verbose)
    build_time = toq()

    it = [0]
    gc()
    tic()
    x = f(b, pcgIts = it, tol=tol, maxits=maxits, verbose=verbose)
    solve_time = toq()

    err = norm(lap(a) * x - b) / norm(b)

    ret = (solve_time, build_time, it[1], err)
    if verbose
      println(ret)
    end
    return ret
  catch
    println("Solver Error.")
    return (Inf, Inf, Inf, Inf)
  end

end

function testSolverTimed(timelimit, solver, a, b, tol, maxits, verbose)

  #@everywhere ts(a,b,c,d,e,f) = testSolver(a,b,c,d,e,f)

  if length(workers()) < 1
    error("Run addprocs() before setting up this code.")
  end
  proc = workers()[1]

  c = Channel(1)
  @async put!(c, remotecall_fetch(()->(testSolver(solver, a, b, tol, maxits, verbose)),proc))
  t0 = time()
  while !isready(c) && (time() - t0) < timelimit
    sleep(0.1)
  end

  ret = (Inf,Inf,Inf,Inf)

  if isready(c)
    ret = fetch(c)

  else
    println("Interrupt process!")
    interrupt(proc)
    sleep(1)
    ret = (Inf,Inf,Inf,Inf)

  end

  return ret

end

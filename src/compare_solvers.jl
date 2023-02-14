#==========================================================
Code for comparing solvers

For comparing possibly flaky code, we need an ability to stop the code
if it runs for too long.
Code that does that is in the compare directory.

It uses threads for Julia code,
and has special code for the matlab solvers

===========================================================#


"""
    SolverTest(solver, name)

Encloses a solver with its name, so that we can compare it in tests
"""
struct SolverTest
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
        dic[name] = typ[]
    end
end



"""
`ret` is the answer returned by a speed test.
This pushed it into the dictionary on which we are storing the tests.
"""
function pushSpeedResult!(dic, name, ret)
    push!(dic["$(name)_solve"],ret[1])
    push!(dic["$(name)_build"],ret[2])
    push!(dic["$(name)_tot"],ret[1]+ret[2])
    push!(dic["$(name)_its"],ret[3])
    push!(dic["$(name)_err"],ret[4])
end


"""
    function speedTestLapSolvers{Tv,Ti}(solvers, dic, a::SparseMatrixCSC{Tv,Ti}, b::Array{Tv,1}; tol::Real=1e-2, maxits=Inf, maxtime=Inf, verbose=false)

Runs many Laplacians solvers.  Puts the build and solve time results into a dictionary dic.  It would be easiest to look at it via DataFrame(dic).  Returns the answer from the last solver.  `solvers` should be an array of `SolverTest`.
"""
function speedTestLapSolvers(solvers, dic, a::SparseMatrixCSC{Tv,Ti}, b::Array{Tv,1}; tol::Real=1e-2, maxits=10000, maxtime=1000, verbose=false, testName="") where {Tv,Ti}

    b = b .- mean(b)

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

    nv = size(a,1)
    ne = nnz(a)
    hash_a = hash(a)

    push!(dic["nv"],nv)
    push!(dic["ne"],ne)
    push!(dic["hash_a"],hash_a)
    push!(dic["testName"],testName)
    
    x = []

    for i in 1:length(solvers)
        solverTest = solvers[i]

        if verbose
            println()
            println(solverTest.name)
        end

        ret = testSolver(solverTest.solver, a, b, tol, maxits, verbose)

        if i == 1
            x = ret[5]
        end
        
        
        pushSpeedResult!(dic, solverTest.name, ret)
    end

    return x

end


function testSolver(solver, a, b, tol, maxits, verbose)

  try

    GC.gc()
    t0 = time()
    f = solver(a, tol=tol, maxits=maxits, verbose=verbose)
    build_time = time() - t0

    it = [0]
    GC.gc()
    
    t0 = time()
    x = f(b, pcgIts = it, tol=tol, maxits=maxits, verbose=verbose)
    solve_time = time() - t0

    err = norm(lap(a) * x .- b) / norm(b)

    ret = (solve_time, build_time, it[1], err, x)
    if verbose
      println("Solve time, build time, iter, err:", (solve_time, build_time, it[1], err))
    end
    return ret
  catch
    println("Solver Error.")
    return (Inf, Inf, Inf, Inf, Inf)
  end

end


function testSolverSddm(solver, M, b, tol, maxits, verbose)

    try
  
      GC.gc()
      t0 = time()
      f = solver(M, tol=tol, maxits=maxits, verbose=verbose)
      build_time = time() - t0
  
      it = [0]
      GC.gc()
      
      t0 = time()
      x = f(b, pcgIts = it, tol=tol, maxits=maxits, verbose=verbose)
      solve_time = time() - t0
  
      err = norm(M * x .- b) / norm(b)
  
      ret = (solve_time, build_time, it[1], err, x)
      if verbose
        println("Solve time, build time, iter, err:",(solve_time, build_time, it[1], err))
      end
      return ret
    catch
      println("Solver Error.")
      return (Inf, Inf, Inf, Inf, Inf)
    end
  
  end
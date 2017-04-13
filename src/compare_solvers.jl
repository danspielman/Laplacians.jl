#==========================================================
Code for comparing solvers
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
function speedTestLapSolvers{Tv,Ti}(solvers, dic, a::SparseMatrixCSC{Tv,Ti}, b::Array{Tv,1}; tol::Real=1e-2, maxits=Inf, maxtime=Inf, verbose=false)

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
    
    for solverTest in solvers
        tic()
        f = solverTest.solver(a, tol=tol, maxtime=maxtime, maxits=maxits, verbose=verbose)
        build_time = toq()

        tic()
        x = f(b, pcgIts = it)
        solve_time = toq()

        err = norm(la * x - b) / norm(b)

        name = solverTest.name
        push!(dic[solvecol(name)],solve_time)
        push!(dic[buildcol(name)],build_time)
        push!(dic[totcol(name)],solve_time + build_time)
        push!(dic[itscol(name)],it[1])
        push!(dic[errcol(name)],err)

    end

    return x

end


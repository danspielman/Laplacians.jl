

function compareSolversOut(n, nruns; maxtime=Inf)

    t1 = time()

    # first, force a compile
    tab = compareSolvers(1000,2)

    tab = compareSolvers(n, nruns, maxtime=maxtime)
    fn = string("compSolvers", n , "x", nruns, ".csv")
    writecsv(fn,tab)

    t2 = time() - t1
    println("Took ", t2, " seconds.")
end


function compareSolvers(n, nruns; maxtime=Inf)

    akpwBuild = Array(Float64,nruns)
    akpwSolve = Array(Float64,nruns)
    primBuild = Array(Float64,nruns)
    primSolve = Array(Float64,nruns)
    augBuild = Array(Float64,nruns)
    augSolve = Array(Float64,nruns)

    nedges = Array(Int64,nruns)
    nlist = Array(Int64,nruns)
    ilist = Array(Int64,nruns)

    for i in 1:nruns
        a = wted_chimera(n,i)
        la = lap(a)
        
        b = randn(n)
        b = b - mean(b)
        
        nlist[i] = n
        ilist[i] = i
        nedges[i] = nnz(a)
        
        tic()
        f = treeSolver(a,akpw)
        t = toq()
        akpwBuild[i] = t
        
        tic()
        x = f(b, maxtime=maxtime)
        t = toq()
        akpwSolve[i] = t
        
        tic()
        f = treeSolver(a,prim)
        t = toq()
        primBuild[i] = t
        
        tic()
        x = f(b, maxtime=maxtime)
        t = toq()
        primSolve[i] = t
        
        tic()
        f = augTreeLapSolver(la)
        t = toq()
        augBuild[i] = t
        
        tic()
        x = f(b, maxtime=maxtime)
        t = toq()
        augSolve[i] = t
        
        print(".")
    end

    tab = ["nlist"; nlist]
    tab = [tab ["ilist"; ilist]]
    tab = [tab ["nedges"; nedges]]
    tab = [tab ["akpwBuild"; akpwBuild]]
    tab = [tab ["akpwSolve"; akpwSolve]]
    tab = [tab ["primBuild"; primBuild]]
    tab = [tab ["primSolve"; primSolve]]
    tab = [tab ["augBuild"; augBuild]]
    tab = [tab ["augSolve"; augSolve]]

    #=
    tab = makeColumn(:nlist)
    tab = [tab makeColumn(:ilist)]
    tab = [tab makeColumn(:nedges)]
    tab = [tab makeColumn(:akpwBuild)]    
    tab = [tab makeColumn(:akpwSolve)]  
    tab = [tab makeColumn(:primBuild)]    
    tab = [tab makeColumn(:primSolve)]  
    tab = [tab makeColumn(:augBuild)]    
    tab = [tab makeColumn(:augSolve)]  
    =#

    return tab
end


function treeSolver(a,treeAlg;maxtime=Inf)
    t = treeAlg(a)
    la = lap(a)
    lt = lap(t)
    f = pcgLapSolver(la,lt;maxtime=maxtime)
end

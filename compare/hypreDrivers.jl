#=
  TRILINOS_HOME needs to be set to where Trilinos lives for using MATLAB,
  for example you could put the following in .bash_profile
  export TRILINOS_HOME="/Users/janedoe/Trilinos_collection/Trilinos_gcc"

  WARNING: currently needs the above line AND NO / (slash) at end, which is not common?

  and need gtimeout, which can get from brew install coreutils
=#

using MatrixMarket
using CSV

include("./MatrixMarketVectors.jl")

function timeLimitHypre(limit, M, b; tol::Real=1e-8, maxits=1000, verbose=false, num_procs=2)
    
    solverHome = ENV["HYPRE_HOME"]
    scriptpath = "$(solverHome)/src/test/ij_print"
    
    matname = "tmpFromJulia_mat"
    vecname = "tmpFromJulia_vec"

    hypreExportMatrixVector(matname,M,vecname,b,num_procs=num_procs)

    tmpOutFileName = "tmpToJulia.csv"

    bt = Inf
    st = Inf
    err = Inf
    iter = Inf

    #cmd = `gtimeout $(limit) $(scriptpath) --verbose=false --muelu-xml=$(scripxmlsettings) --tol=1e-6 --max-iters=$(maxits) --filepath=$(matpath) --rhsfile=$(vecpath) --outputfile=$(tmpOutFileName)`

    cmd = `gtimeout $(limit)  mpirun -np 2 $(scriptpath) -solver 1 -fromfile $(matname) -rhsparcsrfile $(vecname) -print`
    
    try
        run(cmd)
        results = CSV.read(tmpOutFileName)
        
        bt = results[1,:TimeSetup]
        st = results[1,:TimeSolve]
        iter = 0 #we're not recording this atm
        err = results[1,:relresidual]
    catch e
        errtrace = backtrace()
        msg = sprint(showerror, e, errtrace)
        println(msg)
        println("Muelu Belos script died")
    end
        
    if verbose
      println("Build Time: ", bt)
      println("Solve Time: ", st)
      #println("Iterations: ", iter)
      println("error: ", err)
      #d1 = DateTime(start,"d-u-yyyy H:M:S")
      #println("Time to load and start Matlab: $(d1-DateTime(t0))/1000)")
    end

    return (st, bt, iter, err)
end

#=
  TRILINOS_HOME needs to be set to where Trilinos lives for using MATLAB,
  for example you could put the following in .bash_profile
  export TRILINOS_HOME="/Users/janedoe/Trilinos_collection/Trilinos_gcc"

  and need gtimeout, which can get from brew install coreutils
=#

using MatrixMarket
using CSV

include("./MatrixMarketVectors.jl")

function timeLimitMueluBelos(limit, la, b; tol::Real=1e-8, maxits=1000, verbose=false)
    
    trilinosHome = ENV["TRILINOS_HOME"]
    scriptpath = "$trilinosHome/script_belos-muelu_timed/muelu_belos.exe"
    scripxmlsettings = "$trilinosHome/script_belos-muelu_timed/SA_3LVL.xml"
    
    tmpOutFileName = "tmpToJulia.csv"
    matpath = "tmpFromJulia_mat.mtx"
    vecpath = "tmpFromJulia_vec.mtx"

    mmwrite(matpath,la)
    mmwritevec(vecpath,b)

    bt = Inf
    st = Inf
    err = Inf
    iter = Inf

    cmd = `gtimeout $(limit) $(scriptpath) --verbose=false --muelu-xml=$(scripxmlsettings) --tol=1e-6 --max-iters=$(maxits) --filepath=$(matpath) --rhsfile=$(vecpath) --outputfile=$(tmpOutFileName)`

    try
        run(cmd)
        results = CSV.read(tmpOutFileName)
        
        bt = results[1,:TimeSetup]
        st = results[1,:TimeSolve]
        iter = 0 #we're not recording this atm
        err = results[1,:relresidual]
    catch
        println("Muelu Belos Script Died")
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

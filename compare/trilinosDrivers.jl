#=
  TRILINOS_HOME needs to be set to where Trilinos lives for using MATLAB,
  for example you could put the following in .bash_profile
  export TRILINOS_HOME="/Users/me/Trilinos_collection/Trilinos_gcc"

  and need gtimeout, which can get from brew install coreutils
=#

function timeLimitMueluBelos(limit, la, b; tol::Real=1e-8, maxits=1000, verbose=false)

    mf = MatFile("fromJulia.mat","w")
    put_variable(mf, "la", la)
    put_variable(mf, "b", b)
    put_variable(mf, "tol", tol)
    put_variable(mf, "maxits", maxits)
    close(mf)
    
    fn = "$(Pkg.dir("Laplacians"))/matlab/timeLamg.m"
    mat = ENV["MATLAB_HOME"]
    matlab = "$(mat)/bin/matlab"
    cmd = `gtimeout $(limit) $(matlab) -nojvm \< $(fn)`

    t0 = now()

    bt = Inf
    st = Inf
    err = Inf
    iter = Inf
    start = -Inf

    try
        run(cmd)
    catch
        println("Matlab Died")
    end

    try
        mf = MatFile("toJulia.mat")

        bt = get_variable(mf,:bt)
        st = get_variable(mf,:st)
        err = get_variable(mf,:relres)
        iter = get_variable(mf,:iter)
        start = get_variable(mf,:startTime)
        close(mf)
    catch
        println("No .mat file generated")
    end
    


    if verbose
      println("Build Time: ", bt)
      println("Solve Time: ", st)
      println("Iterations: ", iter)
      println("error: ", err)
      d1 = DateTime(start,"d-u-yyyy H:M:S")
      println("Time to load and start Matlab: $(d1-DateTime(t0))/1000)")
    end


    return (st, bt, iter, err)

end

#=
  MATLAB_HOME needs to be set to where Matlab lives for using MATLAB,
  for example you could put the following in .bash_profile
  export MATLAB_HOME="/Applications/MATLAB_R2017a.app"

  and need gtimeout, which can get from brew install coreutils
=#

using MATLAB
using Laplacians

function hypreExportVector(output_filename, b; tol::Real=1e-8, maxits=1000)

    mf = MatFile("julia2matlab2hypre_vector.mat","w")
    #put_variable(mf, "M", M)
    put_variable(mf, "b", b)
    put_variable(mf, "tol", tol)
    put_variable(mf, "maxits", maxits)
    put_variable(mf, "num_procs", 2)
    put_variable(mf, "output_filename", "$(output_filename)")
    close(mf)

    lapdir = dirname(pathof(Laplacians))
    
    # fn = "$(hypretestdir)/matlab/matlab2hypreParVectorsScript.m"
    fn = "$(lapdir)/../hypre/matlab2hypreParVectorsScript.m"
    mat = ENV["MATLAB_HOME"]
    matlab = "$(mat)/bin/matlab"
    cmd = `$(matlab) -nojvm \< $(fn)`
    try
        run(cmd)
    catch e
        errtrace = backtrace()
        msg = sprint(showerror, e, errtrace)
        println(msg)
        println("Matlab Died")
    end
end
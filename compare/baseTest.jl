# this is a script that will test the solvers on some basic graphs.
# It requires MATLAB, CMG and LAMG.

using Laplacians
using MATLAB

include("/Users/spielman/Laplacians/compare/matlabSafe.jl")
include("/Users/spielman/Laplacians/compare/compare_solvers_TL.jl")

ac_deg = function(a; verbose=false, args...)
    approxCholLap(a; params=ApproxCholParams(:deg), verbose=verbose, args...)
end
ac_wdeg = function(a; verbose=false, args...)
    approxCholLap(a; params=ApproxCholParams(:wdeg), verbose=verbose, args...)
end


test1 = SolverTest(ac_deg, "ac")
test2 = SolverTest(ac_wdeg, "ac_w")
test3 = SolverTest(cholLap, "chol")
test4 = SolverTest(cgLapSolver, "cg")
tests = [test1 test2 test3 test4]


# warm up the code
a = chimera(10000,1)
b = randn(10000)
b = b - mean(b)
for solver in tests
    f = solver.solver(a,verbose=true)
    x = f(b)
end

dic = Dict()

using JLD

hostname = split(gethostname(),'.')[1]
dic["host"]  = hostname

fn = "$(hostname)_base.jld"

function vary_tests(tests, dic, a, grName)
  b = randn(a.n)
  b = b - mean(b)

  testVMatlabLap(tests, dic, a, b, verbose=true, tol=1e-8, testName=grName)
  for i in 1:5
    srand(i)
    tn = string("wt",i," ",grName)
    aw = Laplacians.randWeightSub(a)
    testVMatlabLap(tests, dic, aw, b, verbose=true, tol=1e-8, testName=tn)
  end


  at = thicken(a)

  testVMatlabLap(tests, dic, at, b, verbose=true, tol=1e-8,
    testName=string("thick ",grName))

  for i in 1:5
      srand(i)
      tn = string("wt",i," thick ",grName)
      aw = Laplacians.randWeightSub(at)
      testVMatlabLap(tests, dic, aw, b, verbose=true, tol=1e-8, testName=tn)
  end
end

vary_tests(tests,dic,grid2(1000), "grid2(1000)")
save(fn,"dic",dic)

vary_tests(tests,dic,grid3(100), "grid3(100)")
save(fn,"dic",dic)

vary_tests(tests,dic,randRegular(1000000,3), "randRegular(1000000,3)")
save(fn,"dic",dic)

vary_tests(tests,dic,grownGraph(1000000,2), "grownGraph(1000000,2)")
save(fn,"dic",dic)

vary_tests(tests,dic,hyperCube(20), "hyperCube(20)")
save(fn,"dic",dic)

vary_tests(tests,dic, randGenRing(1000000,6), "randGenRing(1000000,6)")
save(fn,"dic",dic)


vary_tests(tests,dic, ErdosRenyiClusterFix(1000000,2), "ErdosRenyiClusterFix(1000000,2)")
save(fn,"dic",dic)

exit()

#=

This is a script that we use to run speed tests on two-dim grids.
Run like julia st_grid2.jl 1
Which will run for approximately 1 hour, increasing the sizes of the graphs.

=#

hours = parse(ARGS[1])

fn = "st_grid2_h$(hours).jld"
println(fn)

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


test_ac = SolverTest(ac_deg, "ac")
test_acw = SolverTest(ac_wdeg, "ac_fast")
test_amg = SolverTest(AMGLapSolver, "pyamg")
test_chol = SolverTest(cholLap, "chol")
    
tests = [test_ac test_acw test_amg test_chol]

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

n = 10000
t0 = time()

while time() - t0 < 60*60*hours

    println("-----------")
    println(n)

    nuse = round(Int,sqrt(n))
    a = grid2(nuse)

    tn = "grid2($(nuse))"
    b = randn(size(a,1))
    b = b - mean(b)
    b = b / norm(b)
    x = testVMatlabLap(tests, dic, a, b, verbose=true, tol=1e-8, testName=tn, test_icc=false)

    save(fn,"dic",dic)

    n = 2*n

end

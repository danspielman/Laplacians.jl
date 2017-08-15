#=

This is a script that we use to run speed tests on two-dim grids.
Run like julia st_grid3.jl 1
Which will run for approximately 1 hour, increasing the sizes of the graphs.

=#

hours = parse(ARGS[1])

fn = "st_grid3_h$(hours).jld"
println(fn)

using Laplacians
using MATLAB

lapdir = Pkg.dir("Laplacians")

include("$(lapdir)/compare/matlabSafe.jl")
include("$(lapdir)/compare/compare_solvers_TL.jl")


ac_deg = function(a; verbose=false, args...)
    approxCholLap(a; params=ApproxCholParams(:deg), verbose=verbose, args...)
end
ac_wdeg = function(a; verbose=false, args...)
    approxCholLap(a; params=ApproxCholParams(:wdeg), verbose=verbose, args...)
end


test_ac = SolverTest(ac_deg, "ac")
test_acw = SolverTest(ac_wdeg, "ac_fast")
test_amg = SolverTest(AMGLapSolver, "pyamg")
    
tests = [test_ac test_acw test_amg]

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

    n = 2*n
    println("-----------")
    println(n)

    nuse = round(Int,n^(1/3))
    a = grid3(nuse)

    tn = "grid3($(nuse))"
    b = randn(size(a,1))
    b = b - mean(b)
    b = b / norm(b)
    x = testVMatlabLap(tests, dic, a, b, verbose=true, tol=1e-8, testName=tn, test_icc=true)

    save(fn,"dic",dic)

end

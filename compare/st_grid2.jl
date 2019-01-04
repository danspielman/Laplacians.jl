#=

This is a script that we use to run speed tests on two-dim grids.
Run like julia st_grid2.jl 1
Which will run for approximately 1 hour, increasing the sizes of the graphs.

=#

hours = Base.parse(Float64,ARGS[1])

fn = "st_grid2_h$(hours).jld"
println(fn)

using Laplacians
using MATLAB
using SparseArrays
using Statistics
using LinearAlgebra

lapdir = dirname(pathof(Laplacians))

include("$(lapdir)/../compare/matlabSafe.jl")
include("$(lapdir)/../compare/trilinosDrivers.jl")
include("$(lapdir)/../compare/compare_solvers_TL.jl")

ac_deg = function(a; verbose=false, args...)
    approxchol_lap(a; params=ApproxCholParams(:deg), verbose=verbose, args...)
end
ac_wdeg = function(a; verbose=false, args...)
    approxchol_lap(a; params=ApproxCholParams(:wdeg), verbose=verbose, args...)
end


test_ac = SolverTest(ac_deg, "ac")
test_acw = SolverTest(ac_wdeg, "ac_fast")
#test_amg = SolverTest(AMGLapSolver, "pyamg")
test_chol = SolverTest(chol_lap, "chol")
    
#tests = [test_ac test_acw test_amg test_chol]
tests = [test_ac test_acw test_chol]
    
# warm up the code
a = chimera(10000,1)
b = randn(10000)
b = b - mean(b)*ones(size(b))
for solver in tests
    f = solver.solver(a,verbose=true)
    x = f(b)
end

dic = Dict()

using JLD2

n = 10000
t0 = time()

while time() - t0 < 60*60*hours

    println("-----------")
    println(n)

    nuse = round(Int,sqrt(n))
    a = grid2(nuse)

    tn = "grid2($(nuse))"
    b = randn(size(a,1))
    b = b - mean(b)*ones(size(b))
    b = b / norm(b)
    x = testVMatlabLap(tests, dic, a, b, verbose=true, tol=1e-8, testName=tn, test_icc=true, test_cmg=false, test_lamg=false, test_muelubelos=false )

    @save fn dic

    global n = 2*n

end
